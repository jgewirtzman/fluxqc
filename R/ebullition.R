#' Detect ebullition (bubble) events as abrupt upward jumps in a closure trace
#'
#' Within one chamber placement over water or saturated soil, diffusive flux
#' produces a smooth rise while bubbles produce abrupt step-ups. This function
#' trims the first and last `trim.secs`, flags every point-to-point increase
#' `>= jump.thresh` as a jump, builds a "de-ebulliated" trace by subtracting
#' the cumulative jump magnitude from all later points, and marks a
#' `buffer.secs` neighbourhood around each jump during which the headspace is
#' still equilibrating.
#'
#' Limitations (documented, not solved): the threshold is a single-sample
#' criterion, so slow bubble trains that raise the concentration over several
#' samples are missed, and analyzer artefacts (pressure transients, breaths on
#' an inlet, valve switches) look exactly like bubbles. In the source campaign
#' every automatically detected event was reviewed by eye and only confirmed
#' traces were kept; treat the output as candidates. `jump.thresh` is in the
#' units of `gastype` and, with `scale.by.dt = TRUE`, is scaled up in
#' proportion to the logging interval beyond 5 s because longer intervals
#' accumulate more diffusive change per step.
#'
#' @param dataframe data.frame with a time column and the gas column for one
#'   placement (a single `UniqueID`).
#' @param gastype gas column, e.g. `"CH4dry_ppb"`.
#' @param time.col POSIXct time column.
#' @param jump.thresh minimum upward point-to-point change that counts as a
#'   jump, in the units of `gastype` (100 ppb = 0.1 ppm, the source default).
#' @param trim.secs seconds discarded at each end of the trace.
#' @param buffer.secs half-width of the equilibration window excluded around
#'   each jump when fitting the diffusive slope.
#' @param scale.by.dt scale `jump.thresh` by `max(1, dt / 5)` where `dt` is the
#'   median logging interval.
#' @return list with `trace` (the trimmed data with `Etime_eb`, `dC`,
#'   `is_jump`, `cum_jump`, `C_deebull`, `near_jump`), `jumps` (rows of
#'   `trace` where `is_jump`), and `summary` (one-row data.frame: `n_points`,
#'   `duration_sec`, `total_rise`, `n_jumps`, `total_ebullitive`,
#'   `max_jump`, `ebullitive_fraction`, `jump_thresh_used`, `dt`).
#' @seealso [split_ebullition()]
#' @examples
#' set.seed(5)
#' t0 <- as.POSIXct("2022-10-23 11:00:00", tz = "UTC")
#' ch4 <- 2000 + 0.5 * (0:599) + rnorm(600, 0, 1)
#' ch4[300:600] <- ch4[300:600] + 400       # one bubble at 300 s
#' d <- data.frame(POSIX.time = t0 + 0:599, CH4dry_ppb = ch4)
#' eb <- detect_ebullition(d, "CH4dry_ppb", jump.thresh = 100)
#' eb$summary
#' @export
detect_ebullition <- function(dataframe, gastype = "CH4dry_ppb", time.col = "POSIX.time",
                              jump.thresh = 100, trim.secs = 20, buffer.secs = 15,
                              scale.by.dt = TRUE) {
  for (col in c(gastype, time.col)) if (!col %in% names(dataframe)) stop("column '", col, "' not found")
  d <- dataframe[order(dataframe[[time.col]]), , drop = FALSE]
  d <- d[!is.na(d[[gastype]]) & !is.na(d[[time.col]]), , drop = FALSE]
  el <- as.numeric(difftime(d[[time.col]], min(d[[time.col]]), units = "secs"))
  keep <- el >= trim.secs & el <= max(el) - trim.secs
  d <- d[keep, , drop = FALSE]
  empty_summary <- function(jt, dt) data.frame(n_points = nrow(d), duration_sec = 0, total_rise = NA_real_,
                                               n_jumps = 0L, total_ebullitive = 0, max_jump = NA_real_,
                                               ebullitive_fraction = NA_real_, jump_thresh_used = jt, dt = dt)
  if (nrow(d) < 3) return(list(trace = d, jumps = d[0, ], summary = empty_summary(jump.thresh, NA_real_)))
  d$Etime_eb <- as.numeric(difftime(d[[time.col]], min(d[[time.col]]), units = "secs"))
  dt <- stats::median(diff(d$Etime_eb))
  jt <- if (scale.by.dt) jump.thresh * max(1, dt / 5) else jump.thresh
  x <- d[[gastype]]
  d$dC <- c(NA, diff(x))
  d$is_jump <- !is.na(d$dC) & d$dC >= jt
  d$cum_jump <- cumsum(ifelse(d$is_jump, d$dC, 0))
  d$C_deebull <- x - d$cum_jump
  jt_times <- d$Etime_eb[d$is_jump]
  d$near_jump <- if (length(jt_times)) vapply(d$Etime_eb, function(t) any(abs(t - jt_times) <= buffer.secs), logical(1)) else FALSE
  jumps <- d[d$is_jump, , drop = FALSE]
  total_rise <- x[length(x)] - x[1]
  total_eb <- sum(jumps$dC)
  summary <- data.frame(
    n_points = nrow(d), duration_sec = max(d$Etime_eb), total_rise = total_rise,
    n_jumps = nrow(jumps), total_ebullitive = total_eb,
    max_jump = if (nrow(jumps)) max(jumps$dC) else NA_real_,
    ebullitive_fraction = if (total_rise > 0) total_eb / total_rise else NA_real_,
    jump_thresh_used = jt, dt = dt)
  rownames(d) <- NULL
  list(trace = d, jumps = jumps, summary = summary)
}

#' Partition a closure into diffusive and ebullitive flux
#'
#' Takes the output of [detect_ebullition()] and computes (1) the diffusive
#' flux from a linear fit to the de-ebulliated trace, excluding the
#' equilibration buffer around each jump, and (2) the ebullitive flux as the
#' summed jump magnitude divided by the trace duration. Both are converted to
#' area fluxes with [flux_term()], so the units follow the gas column
#' (nmol m-2 s-1 for a ppb column).
#'
#' This is a deliberately simple, transparent partition. It uses a linear
#' (not Hutchinson-Mosier) diffusive fit; if you need goFlux's model
#' selection, build a goFlux data frame from `eb$trace` with `C_deebull` as
#' the gas column, `flag = !near_jump`, and run [goFlux::goFlux()] on it. The
#' ebullitive rate assumes the bubbles observed during the closure are
#' representative of the sampling period, which for episodic ebullition is
#' rarely true for a single 5-10 min placement.
#'
#' @param eb output of [detect_ebullition()].
#' @inheritParams flux_term
#' @param Vtot,Area,Tcham,Pcham chamber volume (L), area (cm2), temperature
#'   (C) and pressure (kPa), as in a goFlux auxfile.
#' @return one-row data.frame: `diffusive_slope` (gas units s-1),
#'   `diffusive_r2`, `diffusive_flux`, `ebullitive_rate` (gas units s-1),
#'   `ebullitive_flux`, `total_flux`, `ebullitive_fraction`, `n_jumps`,
#'   `n_points_clean`, `duration_sec`, `flux.term`.
#' @examples
#' set.seed(5)
#' t0 <- as.POSIXct("2022-10-23 11:00:00", tz = "UTC")
#' ch4 <- 2000 + 0.5 * (0:599) + rnorm(600, 0, 1)
#' ch4[300:600] <- ch4[300:600] + 400
#' d <- data.frame(POSIX.time = t0 + 0:599, CH4dry_ppb = ch4)
#' eb <- detect_ebullition(d, "CH4dry_ppb", jump.thresh = 100)
#' split_ebullition(eb, Vtot = 4.318, Area = 324.3, Tcham = 28, Pcham = 101.325)
#' @export
split_ebullition <- function(eb, Vtot, Area, Tcham, Pcham, H2O_mol = 0) {
  tr <- eb$trace
  ft <- flux_term(Vtot, Pcham, Area, Tcham, H2O_mol)
  clean <- tr[!tr$near_jump, , drop = FALSE]
  if (nrow(clean) < 3) {
    slope <- NA_real_; r2 <- NA_real_
  } else {
    m <- stats::lm(C_deebull ~ Etime_eb, data = clean)
    slope <- unname(stats::coef(m)[2]); r2 <- summary(m)$r.squared
  }
  dur <- eb$summary$duration_sec
  eb_rate <- if (eb$summary$n_jumps > 0 && dur > 0) eb$summary$total_ebullitive / dur else 0
  diff_flux <- slope * ft; eb_flux <- eb_rate * ft
  total <- diff_flux + eb_flux
  data.frame(diffusive_slope = slope, diffusive_r2 = r2, diffusive_flux = diff_flux,
             ebullitive_rate = eb_rate, ebullitive_flux = eb_flux, total_flux = total,
             ebullitive_fraction = if (isTRUE(total > 0)) eb_flux / total else NA_real_,
             n_jumps = eb$summary$n_jumps, n_points_clean = nrow(clean),
             duration_sec = dur, flux.term = ft)
}
