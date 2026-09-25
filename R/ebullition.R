# =============================================================================
# Ebullition add-on. Nothing in this file is called by the core pipeline
# (goFlux -> best.flux -> flag_detection -> qc_screens) unless the user opts in
# with process_fluxes(ebullition = "flag" | "split") or calls these functions
# directly. Users of upland (unsaturated) surfaces never need it.
# =============================================================================

#' Detect ebullition (bubble) events within chamber closures
#'
#' Ebullition add-on (opt-in). Over water, saturated soil or flooded
#' vegetation, diffusive flux
#' produces a smooth concentration rise while bubbles produce abrupt step-ups.
#' This function works on a click.peak2-style data frame (several closures,
#' one per `UniqueID`), uses the flagged rows of each closure, and returns the
#' events, the annotated traces and a per-closure summary. Pass the result to
#' [split_ebullition()] for a de-ebulliated goFlux input (diffusive term) and
#' the ebullitive flux, or let [process_fluxes()] do that with
#' `ebullition = "flag"` or `"split"`.
#'
#' @details
#' **Methods.**
#' \describe{
#'   \item{`"iqr"` (default)}{the variable ebullition filter of Hoffmann et al.
#'     (2017). Their R script (Zenodo record 53168, `R-Script-CH4_version_1.0.R`,
#'     lines 365-366 and 457-458, read for this implementation) computes the
#'     point-to-point changes \eqn{\Delta C} of a closure and the fences
#'     \eqn{Q_1 - J \cdot IQR} and \eqn{Q_3 + J \cdot IQR} with `Jumper` J =
#'     0.25, and accepts a moving-window linear fit only when every change in
#'     the window lies inside the fences; the flux is then taken from accepted
#'     windows. fluxqc applies the same fences per closure and calls a change
#'     above the upper fence an ebullition event (with `direction = "both"`
#'     changes below the lower fence are also excluded from the diffusive fit,
#'     but never counted as ebullition mass). Because fluxqc sums the events
#'     as bubble mass rather than merely selecting windows, it needs a fence
#'     that Gaussian noise alone does not cross: the conventional `k = 1.5`
#'     fence lies about 3 SD of the differences out and, on a 600-point 1 Hz
#'     closure, produces roughly one false event per closure from noise
#'     alone, so the default is the "far out" fence `k = 3` (about 4.7 SD,
#'     ~1e-6 per point). Set `k = 0.25` to reproduce Hoffmann's fences
#'     exactly, `k = 1.5` for the conventional fence.}
#'   \item{`"threshold"`}{any upward change of at least `jump.thresh` (gas
#'     units; scaled by `max(1, dt/5)` when `scale.by.dt = TRUE` because
#'     longer logging intervals accumulate more diffusive change per step);
#'     the fixed criterion used in the BlueFlux processing.}
#'   \item{`"runvar"`}{running-variance method after Sø et al. (2024, JGR
#'     Biogeosciences, \doi{10.1029/2024JG008035}): the variance of the
#'     point-to-point changes in a moving window of `runvar_window` samples is
#'     compared with `runvar_cutoff` (gas units squared); contiguous stretches
#'     above the cutoff are events, and each event's mass is the concentration
#'     rise across the stretch. Implemented from the published description
#'     (event = period where running variance exceeds a user cutoff, about
#'     0.5); the authors' code was not available to verify the exact window
#'     definition, so treat the defaults as a starting point.}
#' }
#'
#' **How other packages handle ebullition** (checked September 2026 from the
#' sources or manuals named; "not verified" where a claim could not be
#' checked):
#' \itemize{
#'   \item goFlux 0.2.0: none. The source tree contains no ebullition or
#'     piecewise handling; one LM and one HM fit per closure (verified from the
#'     source snapshot).
#'   \item Hoffmann et al. (2017): the IQR filter described above (verified
#'     from the Zenodo script).
#'   \item `fluxCH4` (github.com/mthayne527/fluxCH4, standalone scripts, not
#'     a package): fits a GAM per `id` (diffusive if R2 > 0.99), otherwise
#'     trains gradient-boosted models with `dismo::gbm.step`, smooths the
#'     residuals, applies a dynamic threshold (SD x multiplier) and a
#'     concentration-difference quantile test to flag points, groups flagged
#'     points into sequences with a 5-point buffer and classifies sequences of
#'     at least `min_ebullition_sequence` as events (verified from its
#'     README).
#'   \item CRAN `flux` (Jurasinski et al.; archived from CRAN 2026-04-15):
#'     `flux.odae()` fits all sub-models from n down to `min.allowed` points
#'     and keeps the best by nRMSE, i.e. it detects and *excludes* outlying
#'     concentration values to find stable linear conditions; it quantifies no
#'     ebullition (verified from the `flux.Rd` manual page).
#'   \item `gasfluxes` 0.7-2: manual offers linear, robust linear (MASS),
#'     HMR and NDFE fits; no ebullition or piecewise handling mentioned
#'     (verified from the CRAN manual).
#'   \item `HMR` 1.0.5: a prefiltering test classifies a series as 'Noise';
#'     no ebullition handling (verified from the CRAN manual).
#'   \item `fluxfinder` 1.2.2 (Wilson et al. 2024, \doi{10.1029/2024JG008208}):
#'     `ffi_fit_models()`, `ffi_hm()`, `ffi_qaqc()`; no ebullition or
#'     piecewise handling mentioned (verified from the CRAN manual).
#'   \item `RespChamberProc` (GitHub): uses `changepoint`/`segmented` to detect
#'     the initial lag (breakpoint) of a closure, not ebullition (verified from
#'     `fluxEstimates.R`).
#'   \item `FluxCalR`: no ebullition mention in its README (verified); the
#'     full source was not inspected (not verified further).
#' }
#'
#' **Limitations.** All three criteria are local tests: slow bubble trains
#' that raise the concentration over many samples are missed by `"iqr"` and
#' `"threshold"`, and analyzer artefacts (pressure transients, breaths on an
#' inlet, valve switches) look exactly like bubbles. Review `events` against
#' the traces before reporting.
#'
#' @param dataframe click.peak2-style data.frame (`UniqueID`, `POSIX.time`,
#'   `flag`, `Etime`, the gas column and chamber columns).
#' @param gastype gas column, e.g. `"CH4dry_ppb"`.
#' @param method `"iqr"`, `"threshold"` or `"runvar"` (see Details).
#' @param by grouping column (default `"UniqueID"`).
#' @param k fence multiplier for `"iqr"` (default 3; 1.5 is the conventional
#'   outlier fence; Hoffmann et al. used 0.25 for window selection).
#' @param jump.thresh minimum upward change for `"threshold"`, gas units.
#' @param runvar_window,runvar_cutoff window (samples) and variance cutoff
#'   (gas units squared) for `"runvar"`.
#' @param direction `"up"` (default) flags upward changes only; `"both"` also
#'   removes downward outliers from the diffusive fit.
#' @param buffer.secs half-width of the equilibration window marked
#'   `near_jump` around each event (excluded from the diffusive fit).
#' @param merge.consecutive treat a run of consecutive flagged samples as one
#'   event (its mass is the summed change, its time the first sample); a
#'   bubble burst often spans 2-10 samples at 1 Hz. `FALSE` reports every
#'   sample separately, as the BlueFlux campaign script did.
#' @param flag.only use only rows with `flag == 1` (default).
#' @param scale.by.dt scale `jump.thresh` by the logging interval.
#' @param time.col POSIXct time column.
#' @return list of class `fluxqc_ebullition` with `events` (one row per
#'   event: `UniqueID`, `POSIX.time`, `Etime`, `dC` (jump size, gas units),
#'   `threshold`), `trace` (rows used, with `dC`, `is_jump`, `is_outlier`,
#'   `cum_jump`, `C_deebull`, `near_jump`) and `summary` (per closure:
#'   `n_points`, `duration_sec`, `total_rise`, `n_events`, `total_ebullitive`,
#'   `max_jump`, `ebullitive_fraction`, `threshold_used`, `dt`, `method`).
#' @references
#' Hoffmann, M., Schulz-Hanke, M., Garcia Alba, J., Jurisch, N., Hagemann, U.,
#' Sachs, T., Sommer, M., & Augustin, J. (2017). A simple calculation algorithm
#' to separate high-resolution CH4 flux measurements into ebullition- and
#' diffusion-derived components. Atmospheric Measurement Techniques, 10,
#' 109-118. \doi{10.5194/amt-10-109-2017}. Script: Zenodo record 53168.
#'
#' So, J. S., et al. (2024). Journal of Geophysical Research: Biogeosciences.
#' \doi{10.1029/2024JG008035}
#' @family ebullition (add-on)
#' @seealso [split_ebullition()], [merge_ebullition()], [process_fluxes()]
#' @examples
#' d <- example_closure()
#' d$CH4dry_ppb[d$t_s >= 300] <- d$CH4dry_ppb[d$t_s >= 300] + 60   # inject one bubble
#' eb <- detect_ebullition(d, "CH4dry_ppb")
#' eb$summary; eb$events
#' @export
detect_ebullition <- function(dataframe, gastype = "CH4dry_ppb",
                              method = c("iqr", "threshold", "runvar"), by = "UniqueID",
                              k = 3, jump.thresh = 100, runvar_window = 5, runvar_cutoff = 0.5,
                              direction = c("up", "both"), buffer.secs = 15,
                              merge.consecutive = TRUE,
                              flag.only = TRUE, scale.by.dt = TRUE, time.col = "POSIX.time") {
  method <- match.arg(method); direction <- match.arg(direction)
  for (col in c(gastype, time.col, by)) if (!col %in% names(dataframe)) stop("column '", col, "' not found")
  d <- dataframe
  if (flag.only) {
    if (!"flag" %in% names(d)) stop("'dataframe' has no 'flag' column; set flag.only = FALSE")
    d <- d[!is.na(d$flag) & d$flag == 1, , drop = FALSE]
  }
  d <- d[!is.na(d[[gastype]]) & !is.na(d[[time.col]]), , drop = FALSE]
  d <- d[order(d[[by]], d[[time.col]]), , drop = FALSE]
  sp <- split(d, d[[by]])
  traces <- list(); summ <- list()
  for (uid in names(sp)) {
    x <- sp[[uid]]
    et <- as.numeric(difftime(x[[time.col]], min(x[[time.col]]), units = "secs"))
    n <- nrow(x)
    dt <- if (n > 1) stats::median(diff(et)) else NA_real_
    y <- x[[gastype]]
    dC <- c(NA, diff(y))
    up <- rep(FALSE, n); dn <- rep(FALSE, n); mass <- rep(0, n); thr <- NA_real_
    if (n >= 4) {
      if (method == "iqr") {
        q <- stats::quantile(dC, c(0.25, 0.75), na.rm = TRUE); iqr <- unname(q[2] - q[1])
        thr <- unname(q[2]) + k * iqr
        up <- !is.na(dC) & dC > thr
        dn <- direction == "both" & !is.na(dC) & dC < unname(q[1]) - k * iqr
        mass[up] <- dC[up]
      } else if (method == "threshold") {
        thr <- if (scale.by.dt && !is.na(dt)) jump.thresh * max(1, dt / 5) else jump.thresh
        up <- !is.na(dC) & dC >= thr
        dn <- direction == "both" & !is.na(dC) & dC <= -thr
        mass[up] <- dC[up]
      } else {                                   # runvar
        thr <- runvar_cutoff
        w <- max(2, min(runvar_window, n - 1))
        dd <- dC[-1]
        rv <- rep(NA_real_, n)
        for (i in seq_len(length(dd) - w + 1)) rv[i + w] <- stats::var(dd[i:(i + w - 1)])
        hot <- !is.na(rv) & rv > thr
        r <- rle(hot)
        ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
        for (j in which(r$values)) {
          a <- max(starts[j] - w, 1); b <- ends[j]        # window covers the w points before the flag
          rise <- y[b] - y[a]
          if (rise > 0) { up[a] <- TRUE; mass[a] <- rise }
          dn[a:b] <- TRUE                                 # exclude the whole stretch from the diffusive fit
        }
        x$runvar <- rv
      }
    }
    # a bubble burst spans several consecutive samples: report it as ONE event
    # whose mass is the summed change, located at its first sample
    if (merge.consecutive && any(up) && method != "runvar") {
      r <- rle(up); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
      for (j in which(r$values & r$lengths > 1)) {
        mass[starts[j]] <- sum(mass[starts[j]:ends[j]]); mass[(starts[j] + 1):ends[j]] <- 0
        up[(starts[j] + 1):ends[j]] <- FALSE
        dn[(starts[j] + 1):ends[j]] <- TRUE            # still excluded from the diffusive fit
      }
    }
    x$Etime_eb <- et; x$dC <- dC
    x$is_jump <- up
    x$is_outlier <- up | dn
    x$cum_jump <- cumsum(mass)
    x$C_deebull <- y - x$cum_jump
    jt <- et[x$is_outlier]
    x$near_jump <- if (length(jt)) vapply(et, function(t) any(abs(t - jt) <= buffer.secs), logical(1)) else rep(FALSE, n)
    x$eb_threshold <- thr
    x$jump_mass <- mass
    traces[[uid]] <- x
    tot_rise <- if (n) y[n] - y[1] else NA_real_
    tot_eb <- sum(mass)
    summ[[uid]] <- data.frame(
      UniqueID = uid, n_points = n, duration_sec = if (n) max(et) else 0, total_rise = tot_rise,
      n_events = sum(up), total_ebullitive = tot_eb,
      max_jump = if (any(up)) max(mass[up]) else NA_real_,
      ebullitive_fraction = if (isTRUE(tot_rise > 0)) tot_eb / tot_rise else NA_real_,
      threshold_used = thr, dt = dt, method = method, stringsAsFactors = FALSE)
  }
  trace <- do.call(rbind, traces); rownames(trace) <- NULL
  ev <- trace[trace$is_jump, c(by, time.col, "Etime_eb", "jump_mass", "eb_threshold"), drop = FALSE]
  names(ev) <- c(by, time.col, "Etime", "dC", "threshold"); rownames(ev) <- NULL
  summary <- do.call(rbind, summ); rownames(summary) <- NULL
  structure(list(events = ev, trace = trace, summary = summary, gastype = gastype, by = by,
                 method = method),
            class = "fluxqc_ebullition")
}

#' Split closures into a de-ebulliated goFlux input and an ebullitive flux
#'
#' Ebullition add-on (opt-in). Takes the output of [detect_ebullition()] and
#' returns (a) a data frame in goFlux shape in which the gas column holds the
#' *de-ebulliated* trace (the cumulative event mass subtracted from every
#' point after each event) and the equilibration buffer around each event has
#' `flag = 0`, so that `goFlux::goFlux()` and `goFlux::best.flux()` run on it
#' unchanged give the diffusive term with goFlux's own model selection; and
#' (b) the per-closure ebullitive flux, the summed event mass divided by the
#' closure duration and converted with [flux_term()] from the chamber
#' columns, in the same units as the diffusive flux (nmol m-2 s-1 for a ppb
#' column).
#'
#' The ebullitive rate assumes the bubbles caught during the closure are
#' representative of the sampling period, which for episodic ebullition is
#' rarely true for a single 5-10 min placement.
#'
#' @param eb output of [detect_ebullition()].
#' @param H2O_col water column used for the flux term (mean of the flagged
#'   rows, ppm); 0 if absent.
#' @return list with `diffusive` (goFlux-ready data.frame; the raw gas is kept
#'   as `<gastype>_raw`) and `ebullitive` (data.frame: `UniqueID`, `n_events`,
#'   `duration_sec`, `total_ebullitive`, `ebullitive_rate`, `flux.term`,
#'   `flux_ebullitive`).
#' @family ebullition (add-on)
#' @seealso [detect_ebullition()], [merge_ebullition()]
#' @examples
#' d <- example_closure()
#' d$CH4dry_ppb[d$t_s >= 300] <- d$CH4dry_ppb[d$t_s >= 300] + 60
#' eb <- detect_ebullition(d, "CH4dry_ppb")
#' sp <- split_ebullition(eb)
#' sp$ebullitive
#' diff_best <- goFlux::best.flux(goFlux::goFlux(sp$diffusive, "CH4dry_ppb"))
#' out <- merge_ebullition(diff_best, sp$ebullitive)
#' out[, c("best.flux", "flux_ebullitive", "flux_total", "n_events")]
#' @export
split_ebullition <- function(eb, H2O_col = "H2O_ppm") {
  if (!inherits(eb, "fluxqc_ebullition")) stop("'eb' must come from detect_ebullition()")
  tr <- eb$trace; g <- eb$gastype; by <- eb$by
  diffusive <- tr
  diffusive[[paste0(g, "_raw")]] <- tr[[g]]
  diffusive[[g]] <- tr$C_deebull
  diffusive$flag <- ifelse(tr$near_jump, 0, 1)
  diffusive <- diffusive[, setdiff(names(diffusive), c("Etime_eb", "dC", "is_jump", "is_outlier",
                                                       "cum_jump", "C_deebull", "near_jump",
                                                       "eb_threshold", "jump_mass", "runvar"))]
  sp <- split(tr, tr[[by]])
  ft <- vapply(sp, function(x) {
    need <- c("Vtot", "Pcham", "Area", "Tcham")
    if (!all(need %in% names(x))) return(NA_real_)
    h2o <- if (H2O_col %in% names(x)) mean(x[[H2O_col]], na.rm = TRUE) / 1e6 else 0
    if (is.na(h2o)) h2o <- 0
    flux_term(x$Vtot[1], x$Pcham[1], x$Area[1], x$Tcham[1], h2o)
  }, numeric(1))
  if (any(is.na(ft))) warning("chamber columns (Vtot, Pcham, Area, Tcham) missing for some closures; flux_ebullitive is NA")
  s <- eb$summary
  s$ebullitive_rate <- ifelse(s$n_events > 0 & s$duration_sec > 0, s$total_ebullitive / s$duration_sec, 0)
  s$flux.term <- ft[s$UniqueID]
  s$flux_ebullitive <- s$ebullitive_rate * s$flux.term
  ebul <- s[, c("UniqueID", "n_events", "duration_sec", "total_ebullitive", "ebullitive_rate",
                "flux.term", "flux_ebullitive")]
  rownames(ebul) <- NULL
  list(diffusive = diffusive, ebullitive = ebul)
}

#' Append the ebullitive component to a best.flux table
#'
#' Ebullition add-on (opt-in).
#'
#' @param best.flux data.frame from [goFlux::best.flux()] run on the
#'   `diffusive` element of [split_ebullition()].
#' @param ebullitive the `ebullitive` element of [split_ebullition()].
#' @return `best.flux` with `flux_ebullitive`, `flux_total`
#'   (`best.flux + flux_ebullitive`), `n_events` and `ebullitive_fraction`
#'   appended; closures absent from `ebullitive` get 0 events.
#' @family ebullition (add-on)
#' @examples
#' # see split_ebullition()
#' @export
merge_ebullition <- function(best.flux, ebullitive) {
  i <- match(best.flux$UniqueID, ebullitive$UniqueID)
  best.flux$flux_ebullitive <- ifelse(is.na(i), 0, ebullitive$flux_ebullitive[i])
  best.flux$n_events <- ifelse(is.na(i), 0L, as.integer(ebullitive$n_events[i]))
  best.flux$flux_total <- best.flux$best.flux + best.flux$flux_ebullitive
  best.flux$ebullitive_fraction <- ifelse(best.flux$flux_total > 0,
                                          best.flux$flux_ebullitive / best.flux$flux_total, NA_real_)
  best.flux
}

# Internal dispatcher used by process_fluxes() ONLY when ebullition != "none".
# Returns list(fit_data, events, summary, ebullitive) where fit_data is the
# concentration data to hand to goFlux (de-ebulliated for "split").
apply_ebullition <- function(d, gastype, mode, method, ids, args, H2O_col) {
  eb <- do.call(detect_ebullition, c(list(dataframe = d[d$UniqueID %in% ids, ], gastype = gastype,
                                          method = method), args))
  out <- list(fit_data = d, events = eb$events, summary = eb$summary, ebullitive = NULL)
  if (mode == "split") {
    sp <- split_ebullition(eb, H2O_col = H2O_col)
    keep <- intersect(names(d), names(sp$diffusive))
    fit <- rbind(d[!d$UniqueID %in% ids, keep, drop = FALSE], sp$diffusive[, keep, drop = FALSE])
    fit <- fit[order(match(fit$UniqueID, unique(d$UniqueID)), fit$POSIX.time), ]
    rownames(fit) <- NULL
    out$fit_data <- fit; out$ebullitive <- sp$ebullitive
  }
  out
}

# Internal: append the add-on's columns to the flux table
append_ebullition_cols <- function(best, ebres, ids, mode, method) {
  best$ebullition_handling <- ifelse(best$UniqueID %in% ids, mode, "none")
  i <- match(best$UniqueID, ebres$summary$UniqueID)
  best$n_ebullition_events <- ifelse(is.na(i), 0L, as.integer(ebres$summary$n_events[i]))
  best$ebullition_flag <- best$n_ebullition_events > 0
  if (mode == "split") {
    j <- match(best$UniqueID, ebres$ebullitive$UniqueID)
    best$flux_diffusive <- best$best.flux
    best$flux_ebullitive <- ifelse(is.na(j), 0, ebres$ebullitive$flux_ebullitive[j])
    best$flux_total <- best$flux_diffusive + best$flux_ebullitive
    best$ebullition_method <- ifelse(best$UniqueID %in% ids, method, NA_character_)
  }
  best
}
