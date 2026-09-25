#' Minimum detectable flux by several published methods (vectorised core)
#'
#' All methods share the chamber form \eqn{\mathrm{MDF} = (\Delta C_{min} / t)
#' \cdot \mathrm{flux.term}}, where `flux.term` = (V/A)(P/RT) (see
#' [flux_term()]) and `t` is the closure length in seconds. They differ in
#' what counts as the minimum detectable concentration change
#' \eqn{\Delta C_{min}}:
#'
#' \describe{
#'   \item{`goflux`}{\eqn{p / t \cdot \mathrm{flux.term}} with `p` = the
#'     instrument precision given to goFlux at import (`prec`). This is
#'     goFlux's built-in [goFlux::MDF()], whose documentation cites
#'     Christiansen et al. (2015) for it.}
#'   \item{`wass{90,95,99}`}{\eqn{z \, \sigma_{global} / t \cdot \mathrm{flux.term}},
#'     the form of Wassmann et al. (2018) -- a *global* (instrument-level)
#'     empirical precision times a multiplier -- generalised to a stated
#'     two-sided confidence: `wass90` uses `qnorm(0.95)` = 1.645, `wass95`
#'     `qnorm(0.975)` = 1.96, `wass99` `qnorm(0.995)`. Wassmann et al.
#'     themselves used the SD of replicate ambient-air analyses with a fixed
#'     k = 3 ("99 percent", citing the MIRSA guidelines); the Monte-Carlo
#'     limits of Parkin, Venterea & Hargreaves (2012) at alpha = 0.05 give the
#'     same dependence on precision and closure time. `wass95` with a
#'     record-wide MAD sigma and `t` in seconds is the method of the companion
#'     data-processing manuscript and of the guidelines paper.}
#'   \item{`chr{90,95,99}`}{\eqn{\sigma_i \cdot 3 \cdot t_{crit}(df) / t \cdot
#'     \mathrm{flux.term}} with the *per-measurement* precision
#'     \eqn{\sigma_i}, `df = nb.obs - 2` and a two-sided t quantile: a
#'     conservative per-closure variant kept for comparison (the column prefix
#'     `chr` is historical). It is NOT the Christiansen et al. (2015) MDF:
#'     their section 2.5 defines MDF as analytic precision divided by
#'     enclosure time (k = 1), and the "SD x 3 x t99" there is a gas-chromatograph
#'     method quantification limit from 15 replicate standards (Corley 2003).
#'     Nickerson (2016) likewise uses k = 1 with an analytical accuracy.}
#'   \item{`2sigma`}{\eqn{2 \sigma_i / t \cdot \mathrm{flux.term}}, the
#'     illustration used in Table 2.5 of the guidelines paper.}
#' }
#'
#' All these single-sample criteria are conservative: they ignore the averaging
#' benefit of fitting a slope to `t` points. The regression standard error of
#' the slope (`LM.SE` in goFlux) is typically several times smaller.
#' [flag_detection()] is the user-facing wrapper for a `best.flux()` table.
#'
#' @param flux.term numeric vector; goFlux `flux.term` column (mol m-2 per
#'   unit mole fraction).
#' @param t numeric vector; closure length in seconds (goFlux uses
#'   `max(Etime) + 1`, which equals `nb.obs` at 1 Hz).
#' @param prec instrument precision (goFlux `prec`), same units as the gas.
#' @param sigma per-measurement empirical precision (from [precision_allan()]
#'   or [precision_by_id()]), same units as the gas.
#' @param sigma.global global (instrument-level) empirical precision, e.g.
#'   [precision_mad()] over a whole record.
#' @param nb.obs number of observations in the fit (for the t quantile);
#'   defaults to `t`.
#' @param conf two-sided confidence levels for the Wassmann and
#'   Christiansen-style variants.
#' @param chr.mult multiplier in the Christiansen-style variant (3).
#' @return data.frame with one row per input element and one column per
#'   method: `MDF_goflux`, `MDF_wass90`, ..., `MDF_chr99`, `MDF_2sigma`.
#'   Methods whose inputs were not supplied are `NA`. Fluxes and MDFs are in
#'   the units implied by the gas column (umol m-2 s-1 for ppm, nmol m-2 s-1
#'   for ppb).
#' @references
#' Christiansen, J. R., Outhwaite, J., & Smukler, S. M. (2015). Comparison of
#' CO2, CH4 and N2O soil-atmosphere exchange measured in static chambers with
#' cavity ring-down spectroscopy and gas chromatography. Agricultural and
#' Forest Meteorology, 211-212, 48-57. \doi{10.1016/j.agrformet.2015.06.004}
#'
#' Parkin, T. B., Venterea, R. T., & Hargreaves, S. K. (2012). Calculating the
#' detection limits of chamber-based soil greenhouse gas flux measurements.
#' Journal of Environmental Quality, 41, 705-715. \doi{10.2134/jeq2011.0394}
#'
#' Wassmann, R., Alberto, M. C., Tirol-Padre, A., Hoang, N. T., Romasanta, R.,
#' Centeno, C. A., & Sander, B. O. (2018). Increasing sensitivity of methane
#' emission measurements in rice through deployment of 'closed chambers' at
#' nighttime. PLoS ONE, 13(2), e0191352. \doi{10.1371/journal.pone.0191352}
#' @seealso [flag_detection()]
#' @examples
#' ft <- flux_term(V_L = 7.08, P_kPa = 101.3, A_cm2 = 4254, T_C = 26)
#' mdf_multi(flux.term = ft, t = 540, prec = 0.9, sigma = 1.1, sigma.global = 1.2)
#' @export
mdf_multi <- function(flux.term, t, prec = NULL, sigma = NULL, sigma.global = NULL,
                      nb.obs = t, conf = c(0.90, 0.95, 0.99), chr.mult = 3) {
  n <- max(length(flux.term), length(t))
  rec <- function(v) if (is.null(v)) rep(NA_real_, n) else rep_len(as.numeric(v), n)
  flux.term <- rec(flux.term); t <- rec(t); prec <- rec(prec); sigma <- rec(sigma)
  sigma.global <- rec(sigma.global); nb.obs <- rec(nb.obs)
  out <- data.frame(MDF_goflux = prec / t * flux.term)
  lab <- function(c) sprintf("%02d", round(100 * c))
  for (c in conf) {
    z <- stats::qnorm(1 - (1 - c) / 2)
    out[[paste0("MDF_wass", lab(c))]] <- z * sigma.global / t * flux.term
  }
  df <- pmax(nb.obs - 2, 1)
  for (c in conf) {
    tc <- stats::qt(1 - (1 - c) / 2, df = df)
    out[[paste0("MDF_chr", lab(c))]] <- chr.mult * sigma * tc / t * flux.term
  }
  out$MDF_2sigma <- 2 * sigma / t * flux.term
  out[] <- lapply(out, abs)
  out
}

#' Append an empirical precision, MDF and detection flag to a best.flux table
#'
#' Adds exactly four columns to a [goFlux::best.flux()] (or
#' [goFlux::goFlux()]) table: `sigma_emp` (the precision by the chosen
#' estimator), `MDF_emp` (the minimum detectable flux by the chosen method),
#' `below_MDF_emp` (logical) and `MDF_emp_method` (a string such as
#' `"wassmann | mad | 0.95"`), plus `detected_emp` and `det_class_emp`
#' (`"emission"`, `"uptake"`, `"below detection"`). goFlux's own columns
#' (`MDF`, `prec`, `best.flux`, `quality.check`, ...) are never modified or
#' duplicated; fluxqc only appends. Nothing is deleted: rows are flagged.
#'
#' The default precision and MDF form are those of the companion
#' data-processing manuscript (`mdf_FINAL_precision_and_detection.R`):
#' \eqn{\sigma} = MAD(\eqn{\Delta x})/\eqn{\sqrt 2} over the *entire* 1 Hz record
#' of each field period, and \eqn{\mathrm{MDF} = z \, \sigma / t \cdot
#' \mathrm{flux.term}}. The stated confidence is 95 percent two-sided,
#' `z = qnorm(0.975)` = 1.96, i.e. the 2-sigma criterion (the companion
#' manuscript used `qnorm(0.95)` = 1.645, its "90% confidence"; that is
#' `conf = 0.90`). Hence `precision = "mad"`, `mdf = "wassmann"`, `conf = 0.95`.
#' This is applied *after* `best.flux()`; the `prec` argument given to goFlux's
#' `import.*()` functions feeds only goFlux's own `MDF` column.
#'
#' @section Precision estimators (`precision`):
#' \describe{
#'   \item{`"mad"` (default)}{\eqn{\sigma = MAD(\Delta x)/\sqrt{2}}
#'     ([precision_mad()]) over *all* rows of `traces` (shoulders and dead
#'     bands included), computed once per `group` (or once for the whole
#'     table when `group = NULL`). At 1 Hz every second-to-second difference
#'     is instrument-dominated, a constant flux shifts all differences equally
#'     and cancels, and MAD ignores the isolated large differences from
#'     moves, breaths and chamber changes. Every closure in a group gets the
#'     same `sigma_emp`.}
#'   \item{`"allan"`}{per closure, \eqn{\sigma = SD(\Delta x)/\sqrt{2}} on the
#'     flagged rows ([precision_allan()]).}
#'   \item{`"rolling"`}{per closure, quietest-window scan of the flagged rows
#'     ([precision_rolling()]).}
#'   \item{`"datasheet"`}{the instrument precision goFlux carries in its
#'     `prec` column.}
#'   \item{`"custom"`}{a per-closure value supplied through `sigma` (e.g. a
#'     bench sigma from a closed-loop test, Box 2 of the guidelines paper).}
#' }
#'
#' @section MDF methods (`mdf`):
#' \describe{
#'   \item{`"wassmann"` (default)}{Wassmann et al. (2018):
#'     \eqn{z \, \sigma_{global} / t \cdot \mathrm{flux.term}} with the two-sided
#'     normal quantile `z = qnorm(1 - (1 - conf)/2)` and
#'     \eqn{\sigma_{global}} = `sigma.global` if given, else the median of
#'     the chosen estimator across the table (for `"mad"` that is the
#'     record-wide value itself). The record is checked for a change of
#'     logging interval ([precision_mad_runs()]); if one is found, each closure
#'     gets the sigma of the interval it was logged at (from
#'     [closure_seconds()]) and a warning is given.}
#'   \item{`"christiansen"`}{\eqn{\sigma_i \cdot 3 \cdot t_{crit}(df) / t
#'     \cdot \mathrm{flux.term}}, per closure, `df = nb.obs - 2`, two-sided t
#'     quantile: a conservative per-closure variant kept for comparison. The
#'     name is historical and misleading -- see [mdf_multi()]: Christiansen et
#'     al. (2015) define MDF = precision / t (k = 1); their 3 x SD x t99 is a
#'     GC quantification limit, not a per-closure flux threshold.}
#'   \item{`"two_sigma"`}{\eqn{2 \sigma_i / t \cdot \mathrm{flux.term}}, the
#'     illustration of Table 2.5 of the guidelines paper.}
#'   \item{`"goflux"`}{\eqn{\sigma / t \cdot \mathrm{flux.term}}; with
#'     `precision = "datasheet"` this equals goFlux's own `MDF` column
#'     (`p / t * flux.term`, goFlux citing Christiansen et al. 2015).}
#' }
#' `conf` applies to `wassmann` and `christiansen` only and is two-sided
#' (0.95 -> `qnorm(0.975)` = 1.96, the 2-sigma criterion; 0.90 -> `qnorm(0.95)`).
#'
#' @section Extras:
#' `extra = list(precision = ..., mdf = ..., conf = ...)` (or `extra = "all"`
#' for every estimator, method and 0.90/0.95/0.99) additionally appends
#' `sigma_<precision>` for each extra estimator and
#' `MDF_<mdf>_<precision>[_<conf>]` (two-digit percentage) for each
#' combination, e.g. `MDF_christiansen_allan_95`, for cross-study comparison.
#' These are comparison columns only; the reported values stay in the `*_emp`
#' columns.
#'
#' @param flux.result data.frame from [goFlux::best.flux()] or
#'   [goFlux::goFlux()]; needs `UniqueID`, `flux.term`, `nb.obs` and the flux
#'   column.
#' @param traces the concentration data the fluxes were computed from
#'   (output of [goFlux::click.peak2()], [click_peak2_stacked()],
#'   [auto_id_rise()] or [windows_from_table()]); for `"mad"` the whole record
#'   is used, for `allan`/`rolling` its flagged rows.
#' @param gastype gas column of `traces`; inferred when `traces` holds exactly
#'   one recognised gas column.
#' @param precision one estimator (see Details).
#' @param mdf one MDF method (see Details).
#' @param conf two-sided confidence level for `wassmann` / `christiansen`.
#' @param group name of a column in `flux.result` defining the field periods
#'   over which `"mad"` is computed (campaign, day, instrument); `NULL` = one
#'   record. Trace rows are assigned to groups through their `UniqueID`.
#' @param extra `NULL` (default), `"all"`, or a list with any of `precision`,
#'   `mdf`, `conf` (character/numeric vectors) for comparison columns.
#' @param sigma optional custom per-closure precision: a data.frame with
#'   `UniqueID` and `sigma`, or a numeric vector aligned with `flux.result`
#'   (used by `precision = "custom"`).
#' @param sigma.global optional global precision for `wassmann` (a number, or
#'   a named vector/list keyed by estimator name).
#' @param t closure length in seconds. Default `NULL`: taken from `traces`
#'   through [closure_seconds()] (span of the flagged rows plus one logging
#'   interval, which equals goFlux's `max(Etime) + 1` at 1 Hz); without
#'   `traces` it falls back to `nb.obs`, which assumes 1 Hz logging and is
#'   wrong by the logging interval for records logged every 5 or 10 s
#'   (a warning says so).
#' @param flux.col flux column to test.
#' @param rolling.window window (samples) for `"rolling"`.
#' @param chr.mult multiplier in the Christiansen-style method (3).
#' @return `flux.result` with `sigma_emp`, `MDF_emp`, `below_MDF_emp`,
#'   `MDF_emp_method`, `detected_emp`, `det_class_emp` appended (and the
#'   `extra` columns if requested). The global precisions used are stored in
#'   the attribute `"sigma_global"`.
#' @references
#' Christiansen, J. R., Outhwaite, J., & Smukler, S. M. (2015). Agricultural and
#' Forest Meteorology, 211-212, 48-57. \doi{10.1016/j.agrformet.2015.06.004}
#'
#' Wassmann, R., et al. (2018). PLoS ONE, 13(2), e0191352.
#' \doi{10.1371/journal.pone.0191352}
#' @examples
#' d <- example_closure()
#' best <- goFlux::best.flux(goFlux::goFlux(d, "CH4dry_ppb"))
#' out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb")
#' out[, c("best.flux", "MDF", "sigma_emp", "MDF_emp", "below_MDF_emp",
#'         "MDF_emp_method", "det_class_emp")]
#' # comparison columns
#' cmp <- flag_detection(best, traces = d, gastype = "CH4dry_ppb",
#'                       extra = list(precision = c("allan", "datasheet"),
#'                                    mdf = c("goflux", "christiansen", "two_sigma")))
#' grep("^MDF_", names(cmp), value = TRUE)
#' @export
flag_detection <- function(flux.result, traces = NULL, gastype = NULL,
                           precision = "mad", mdf = "wassmann", conf = 0.95, group = NULL,
                           extra = NULL, sigma = NULL, sigma.global = NULL, t = NULL,
                           flux.col = "best.flux", rolling.window = 30, chr.mult = 3) {
  fr <- as.data.frame(flux.result)
  for (col in c("UniqueID", "flux.term", "nb.obs", flux.col))
    if (!col %in% names(fr)) stop("'flux.result' has no column '", col, "'")
  prec_ok <- c("mad", "allan", "rolling", "datasheet", "custom")
  mdf_ok <- c("wassmann", "christiansen", "two_sigma", "goflux")
  precision <- match.arg(precision, prec_ok); mdf <- match.arg(mdf, mdf_ok)
  if (length(conf) != 1) stop("'conf' must be a single number; use 'extra' for several")
  ex <- NULL
  if (identical(extra, "all")) ex <- list(precision = setdiff(prec_ok, if (is.null(sigma)) "custom" else character(0)),
                                          mdf = mdf_ok, conf = c(0.90, 0.95, 0.99))
  else if (is.list(extra)) {
    ex <- list(precision = if (is.null(extra$precision)) precision else match.arg(extra$precision, prec_ok, several.ok = TRUE),
               mdf = if (is.null(extra$mdf)) mdf else match.arg(extra$mdf, mdf_ok, several.ok = TRUE),
               conf = if (is.null(extra$conf)) conf else extra$conf)
  } else if (!is.null(extra)) stop("'extra' must be NULL, \"all\" or a list")
  all_prec <- unique(c(precision, ex$precision))
  needs_traces <- intersect(all_prec, c("mad", "allan", "rolling"))
  if (is.null(traces) && length(needs_traces))
    stop("precision estimator(s) ", paste(needs_traces, collapse = ", "), " need 'traces'")
  if ("custom" %in% all_prec && is.null(sigma)) stop("precision = \"custom\" needs 'sigma'")
  if (!is.null(traces)) {
    if (is.null(gastype)) {
      g <- intersect(goflux_gas_columns(), names(traces))
      if (length(g) != 1) stop("give 'gastype'; traces hold ", length(g), " recognised gas columns")
      gastype <- g
    }
    if (!gastype %in% names(traces)) stop("column '", gastype, "' not in 'traces'")
  }
  if (!is.null(group) && !group %in% names(fr)) stop("'group' column '", group, "' not in 'flux.result'")
  if (is.null(t)) {
    if (!is.null(traces) && all(c("UniqueID", "POSIX.time", "flag") %in% names(traces))) {
      cs <- closure_seconds(traces)
      t <- cs$t[match(fr$UniqueID, cs$UniqueID)]
      miss <- is.na(t)
      if (any(miss)) t[miss] <- fr$nb.obs[miss]
    } else {
      t <- fr$nb.obs
      warning("'t' taken as nb.obs, which assumes 1 Hz logging; pass 't' or 'traces' ",
              "for records logged at other intervals")
    }
  }
  n <- nrow(fr)

  sig <- list()
  for (p in all_prec) {
    sig[[p]] <- switch(p,
      mad = {
        # whole record per group (the FINAL-script method), not per closure;
        # split by logging interval if the record changed interval mid-way
        mad_for <- function(rows) {
          r <- precision_mad_runs(traces[[gastype]][rows], traces$POSIX.time[rows])
          if (nrow(r) == 1 || all(is.na(r$dt_s))) return(rep(r$sigma[1], n))
          cs <- closure_seconds(traces[rows, ])
          dtc <- cs$dt[match(fr$UniqueID, cs$UniqueID)]
          r$sigma[vapply(dtc, function(v) if (is.na(v)) 1L else which.min(abs(r$dt_s - v)), integer(1))]
        }
        if (is.null(group)) mad_for(rep(TRUE, nrow(traces)))
        else {
          gv <- as.character(fr[[group]])
          tg <- gv[match(traces$UniqueID, fr$UniqueID)]
          out <- rep(NA_real_, n)
          for (g in unique(gv)) out[gv == g] <- mad_for(!is.na(tg) & tg == g)[gv == g]
          out
        }
      },
      datasheet = {
        if (!"prec" %in% names(fr)) stop("'datasheet' precision needs a 'prec' column")
        fr$prec
      },
      allan = {
        s <- precision_by_id(traces, gastype, method = "allan")
        s$sigma[match(fr$UniqueID, s$UniqueID)]
      },
      rolling = {
        tr <- traces[!is.na(traces$flag) & traces$flag == 1, ]
        sp <- split(tr[[gastype]], tr$UniqueID)
        r <- vapply(sp, function(v) suppressWarnings(precision_rolling(v, rolling.window)$sd), numeric(1))
        r[match(fr$UniqueID, names(sp))]
      },
      custom = if (is.data.frame(sigma)) sigma$sigma[match(fr$UniqueID, sigma$UniqueID)]
               else rep_len(as.numeric(sigma), n))
  }
  glob <- list()
  for (p in names(sig)) {
    glob[[p]] <- if (is.null(sigma.global)) stats::median(sig[[p]], na.rm = TRUE)
                 else if (length(sigma.global) == 1 && is.null(names(sigma.global))) as.numeric(sigma.global)
                 else if (p %in% names(sigma.global)) as.numeric(sigma.global[[p]])
                 else stats::median(sig[[p]], na.rm = TRUE)
  }
  df <- pmax(fr$nb.obs - 2, 1)
  q2 <- function(c) 1 - (1 - c) / 2                          # two-sided level -> upper quantile
  one_mdf <- function(m, p, c) abs(switch(m,
    goflux = sig[[p]] / t * fr$flux.term,
    two_sigma = 2 * sig[[p]] / t * fr$flux.term,
    wassmann = if (is.null(group) || p != "mad") stats::qnorm(q2(c)) * glob[[p]] / t * fr$flux.term
               else stats::qnorm(q2(c)) * sig[[p]] / t * fr$flux.term,   # per-group record-wide sigma
    christiansen = chr.mult * sig[[p]] * stats::qt(q2(c), df = df) / t * fr$flux.term))
  flux <- fr[[flux.col]]

  m <- one_mdf(mdf, precision, conf)
  fr$sigma_emp <- sig[[precision]]
  fr$MDF_emp <- m
  fr$below_MDF_emp <- ifelse(is.na(m), NA, abs(flux) < m)
  fr$MDF_emp_method <- if (mdf %in% c("christiansen", "wassmann")) paste(mdf, precision, conf, sep = " | ")
                       else paste(mdf, precision, sep = " | ")
  fr$detected_emp <- ifelse(is.na(m), NA, abs(flux) > m)
  fr$det_class_emp <- ifelse(is.na(m), NA_character_,
                             ifelse(flux > m, "emission", ifelse(flux < -m, "uptake", "below detection")))

  if (!is.null(ex)) {
    lab <- function(c) sprintf("%02d", round(100 * c))
    for (p in ex$precision) fr[[paste0("sigma_", p)]] <- sig[[p]]
    for (mm in ex$mdf) for (p in ex$precision) {
      if (mm %in% c("goflux", "two_sigma")) fr[[paste("MDF", mm, p, sep = "_")]] <- one_mdf(mm, p, NA)
      else for (c in ex$conf) fr[[paste("MDF", mm, p, lab(c), sep = "_")]] <- one_mdf(mm, p, c)
    }
  }
  attr(fr, "sigma_global") <- unlist(glob)
  fr
}

#' Closure length in seconds from flagged concentration data
#'
#' goFlux's `nb.obs` counts samples, not seconds; the two differ by the logging
#' interval for analyzers set to report every 5 or 10 s. This returns, per
#' `UniqueID`, the span of the flagged rows plus one logging interval (the
#' median spacing of `POSIX.time` inside the window), which equals goFlux's
#' `max(Etime) + 1` for 1 Hz data.
#'
#' @param traces data.frame with `UniqueID`, `POSIX.time` and `flag`.
#' @return data.frame with `UniqueID`, `t` (seconds), `dt` (logging interval,
#'   seconds) and `n` (flagged samples).
#' @examples
#' d <- example_closure()
#' closure_seconds(d)
#' @export
closure_seconds <- function(traces) {
  for (col in c("UniqueID", "POSIX.time", "flag"))
    if (!col %in% names(traces)) stop("'traces' has no column '", col, "'")
  tr <- traces[!is.na(traces$flag) & traces$flag == 1, ]
  sp <- split(tr$POSIX.time, as.character(tr$UniqueID))
  one <- function(p) {
    p <- sort(as.numeric(p)); n <- length(p)
    if (n < 2) return(c(NA_real_, NA_real_, n))
    dt <- stats::median(diff(p))
    c(diff(range(p)) + dt, dt, n)
  }
  m <- t(vapply(sp, one, numeric(3)))
  data.frame(UniqueID = names(sp), t = m[, 1], dt = m[, 2], n = as.integer(m[, 3]),
             stringsAsFactors = FALSE, row.names = NULL)
}
