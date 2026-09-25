#' Minimum detectable flux by several published methods (vectorised core)
#'
#' All methods share the chamber form \eqn{\mathrm{MDF} = (\Delta C_{min} / t)
#' \cdot \mathrm{flux.term}}, where `flux.term` = (V/A)(P/RT) (see
#' [flux_term()]) and `t` is the closure length in seconds. They differ in
#' what counts as the minimum detectable concentration change
#' \eqn{\Delta C_{min}}:
#'
#' \describe{
#'   \item{`goflux`}{\eqn{\Delta C_{min} = p}, the manufacturer precision
#'     (goFlux's [goFlux::MDF()] with `prec`). No confidence level.}
#'   \item{`wass{90,95,99}`}{Wassmann et al. (2018): \eqn{z_{\alpha} \sigma_{global}},
#'     a *global* empirical precision for the instrument times a normal
#'     quantile (one-sided 90/95/99 percent).}
#'   \item{`chr{90,95,99}`}{Christiansen et al. (2015) style:
#'     \eqn{3 \sigma_i t_{\alpha, n-2}}, the *per-measurement* empirical
#'     precision times 3 times a Student-t quantile with `nb.obs - 2` degrees
#'     of freedom.}
#'   \item{`2sigma`}{the simple criterion used in the companion guidelines
#'     paper (Table 2.5): \eqn{\mathrm{MDF} = (2\sigma_i / t)(V/A)(P/RT)}.}
#' }
#'
#' All these single-sample criteria are conservative: they ignore the averaging
#' benefit of fitting a slope to `t` points. The regression standard error of
#' the slope (`LM.SE` in goFlux) is typically several times smaller.
#' [flag_detection()] is the user-facing wrapper that computes chosen
#' combinations of precision estimator and MDF method for a `best.flux()`
#' table.
#'
#' @param flux.term numeric vector; goFlux `flux.term` column (mol m-2 per
#'   unit mole fraction).
#' @param t numeric vector; closure length in seconds (goFlux uses
#'   `max(Etime) + 1`, which equals `nb.obs` at 1 Hz).
#' @param prec manufacturer precision (goFlux `prec`), same units as the gas.
#' @param sigma per-measurement empirical precision (from [precision_allan()]
#'   or [precision_by_id()]), same units as the gas.
#' @param sigma.global global (instrument-level) empirical precision, e.g. the
#'   campaign median of `sigma` or the output of [precision_mad()] on a whole
#'   record.
#' @param nb.obs number of observations in the fit (for the t quantile);
#'   defaults to `t`.
#' @param conf confidence levels for the Wassmann and Christiansen variants.
#' @param chr.mult multiplier in the Christiansen variant (3).
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
    z <- stats::qnorm(c)
    out[[paste0("MDF_wass", lab(c))]] <- z * sigma.global / t * flux.term
  }
  df <- pmax(nb.obs - 2, 1)
  for (c in conf) {
    tc <- stats::qt(c, df = df)
    out[[paste0("MDF_chr", lab(c))]] <- chr.mult * sigma * tc / t * flux.term
  }
  out$MDF_2sigma <- 2 * sigma / t * flux.term
  out[] <- lapply(out, abs)
  out
}

#' Append precision estimates, MDFs and detection flags to a best.flux table
#'
#' Computes every requested combination of a *precision estimator* and an
#' *MDF method* (at every requested confidence level where the method has one)
#' for each row of a [goFlux::best.flux()] or [goFlux::goFlux()] result,
#' appends them as columns, and then sets the single reported `MDF`,
#' `below_MDF` and `MDF_method` from `primary`. Nothing is deleted: rows are
#' flagged only. This is applied *after* `best.flux()`; the `prec` argument
#' given to goFlux's `import.*()` functions only feeds goFlux's own `MDF`
#' column (equal to the `goflux` / `datasheet` combination here).
#'
#' @section Precision estimators (`precision`):
#' \describe{
#'   \item{`"allan"`}{per closure, \eqn{\sigma = SD(\Delta x)/\sqrt{2}} on the
#'     flagged rows of `traces` ([precision_allan()]).}
#'   \item{`"mad"`}{per closure, \eqn{\sigma = MAD(\Delta x)/\sqrt{2}}
#'     ([precision_mad()]); robust to isolated jumps.}
#'   \item{`"rolling"`}{per closure, quietest-window scan of the flagged rows
#'     ([precision_rolling()]).}
#'   \item{`"datasheet"`}{the manufacturer precision goFlux already carries in
#'     its `prec` column (set with the `prec` argument of `import.*()`).}
#' }
#' A custom per-closure estimate (e.g. a bench sigma from a closed-loop test,
#' Box 2 of the guidelines paper) can be supplied through `sigma` and is
#' labelled `"custom"`.
#'
#' @section MDF methods (`mdf`):
#' \describe{
#'   \item{`"goflux"`}{\eqn{\sigma / t \cdot \mathrm{flux.term}}; with
#'     `precision = "datasheet"` this is goFlux's own `MDF` column.}
#'   \item{`"wassmann"`}{Wassmann et al. (2018):
#'     \eqn{z_{conf} \, \sigma_{global} / t \cdot \mathrm{flux.term}} where
#'     \eqn{\sigma_{global}} is `sigma.global` if given, else the median of the
#'     per-closure estimate across the table (for `"datasheet"`, the `prec`
#'     value itself).}
#'   \item{`"christiansen"`}{Christiansen et al. (2015):
#'     \eqn{3 \, \sigma_i \, t_{conf, n-2} / t \cdot \mathrm{flux.term}} with the
#'     per-closure \eqn{\sigma_i}.}
#'   \item{`"two_sigma"`}{\eqn{2 \sigma_i / t \cdot \mathrm{flux.term}}
#'     (guidelines paper, Table 2.5).}
#' }
#'
#' @section Column naming:
#' `sigma_<precision>` for each estimator; `MDF_<mdf>_<precision>` for
#' methods without a confidence level (`goflux`, `two_sigma`) and
#' `MDF_<mdf>_<precision>_<conf>` with `conf` as a two-digit percentage for
#' `wassmann` and `christiansen`, e.g. `MDF_christiansen_allan_95`. Each MDF
#' column has a matching `below_...` logical. goFlux's original `MDF` column
#' is kept as `MDF_goflux_orig`; the new `MDF`, `below_MDF`, `MDF_method`,
#' `detected` and `det_class` (`"emission"`, `"uptake"`, `"below detection"`)
#' come from `primary`.
#'
#' @param flux.result data.frame from [goFlux::best.flux()] or
#'   [goFlux::goFlux()]; needs `UniqueID`, `flux.term`, `nb.obs` and the flux
#'   column.
#' @param traces the flagged concentration data the fluxes were computed from
#'   (output of [goFlux::click.peak2()], [click_peak2_stacked()],
#'   [auto_id_rise()] or [windows_from_table()]); required for the
#'   `allan`, `mad` and `rolling` estimators.
#' @param gastype gas column of `traces`; inferred when `traces` holds exactly
#'   one recognised gas column.
#' @param precision character vector of estimators to compute (see Details).
#'   The default computes all four; when `traces` is missing and `precision`
#'   was not given explicitly, only `"datasheet"` is used.
#' @param mdf character vector of MDF methods to compute (see Details).
#' @param conf confidence levels for `wassmann` and `christiansen`.
#' @param primary list with elements `precision`, `mdf` and (where relevant)
#'   `conf` naming the one combination reported in `MDF` / `below_MDF`.
#' @param sigma optional custom per-closure precision: a data.frame with
#'   `UniqueID` and `sigma`, or a numeric vector aligned with `flux.result`.
#' @param sigma.global optional global precision for the Wassmann method: a
#'   single number, or a named list/vector keyed by estimator name.
#' @param t closure length in seconds; default `nb.obs` (1 Hz).
#' @param flux.col flux column to test.
#' @param rolling.window window (samples) for `"rolling"`.
#' @param chr.mult multiplier in the Christiansen method (3).
#' @return `flux.result` with the added columns; the global precisions used
#'   are stored in the attribute `"sigma_global"`.
#' @references
#' Christiansen, J. R., Outhwaite, J., & Smukler, S. M. (2015). Agricultural and
#' Forest Meteorology, 211-212, 48-57. \doi{10.1016/j.agrformet.2015.06.004}
#'
#' Wassmann, R., et al. (2018). PLoS ONE, 13(2), e0191352.
#' \doi{10.1371/journal.pone.0191352}
#' @examples
#' d <- example_closure()
#' res  <- goFlux::goFlux(d, "CH4dry_ppb")
#' best <- goFlux::best.flux(res)
#' out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb",
#'                       precision = c("allan", "datasheet"),
#'                       mdf = c("goflux", "christiansen", "two_sigma"), conf = 0.95,
#'                       primary = list(precision = "allan", mdf = "christiansen", conf = 0.95))
#' out[, c("sigma_allan", "sigma_datasheet", "MDF_goflux_datasheet",
#'         "MDF_christiansen_allan_95", "MDF_two_sigma_allan", "MDF", "MDF_method", "det_class")]
#' @export
flag_detection <- function(flux.result, traces = NULL, gastype = NULL,
                           precision = c("allan", "mad", "rolling", "datasheet"),
                           mdf = c("goflux", "wassmann", "christiansen", "two_sigma"),
                           conf = c(0.90, 0.95, 0.99),
                           primary = list(precision = "allan", mdf = "christiansen", conf = 0.95),
                           sigma = NULL, sigma.global = NULL, t = NULL,
                           flux.col = "best.flux", rolling.window = 30, chr.mult = 3) {
  fr <- as.data.frame(flux.result)
  for (col in c("UniqueID", "flux.term", "nb.obs", flux.col))
    if (!col %in% names(fr)) stop("'flux.result' has no column '", col, "'")
  explicit <- !missing(precision)
  precision <- match.arg(precision, several.ok = TRUE)
  mdf <- match.arg(mdf, several.ok = TRUE)
  n <- nrow(fr)
  needs_traces <- intersect(precision, c("allan", "mad", "rolling"))
  if (is.null(traces) && length(needs_traces)) {
    if (explicit) stop("precision estimator(s) ", paste(needs_traces, collapse = ", "),
                       " need 'traces'")
    precision <- setdiff(precision, needs_traces)
    if (!length(precision)) precision <- "datasheet"
    message("no 'traces' given: using precision = \"", paste(precision, collapse = "\", \""), "\"")
  }
  if (!is.null(traces)) {
    if (is.null(gastype)) {
      g <- intersect(goflux_gas_columns(), names(traces))
      if (length(g) != 1) stop("give 'gastype'; traces hold ", length(g), " recognised gas columns")
      gastype <- g
    }
    if (!gastype %in% names(traces)) stop("column '", gastype, "' not in 'traces'")
  }
  if (is.null(t)) t <- fr$nb.obs

  # ---- per-closure sigma for every estimator -------------------------------
  sig <- list()
  for (p in precision) {
    sig[[p]] <- switch(p,
      datasheet = {
        if (!"prec" %in% names(fr)) stop("'datasheet' precision needs a 'prec' column")
        fr$prec
      },
      allan = , mad = {
        s <- precision_by_id(traces, gastype, method = p)
        s$sigma[match(fr$UniqueID, s$UniqueID)]
      },
      rolling = {
        tr <- traces[!is.na(traces$flag) & traces$flag == 1, ]
        sp <- split(tr[[gastype]], tr$UniqueID)
        r <- vapply(sp, function(v) suppressWarnings(precision_rolling(v, rolling.window)$sd), numeric(1))
        r[match(fr$UniqueID, names(sp))]
      })
  }
  if (!is.null(sigma)) {
    sig$custom <- if (is.data.frame(sigma)) sigma$sigma[match(fr$UniqueID, sigma$UniqueID)]
                  else rep_len(as.numeric(sigma), n)
  }
  for (p in names(sig)) fr[[paste0("sigma_", p)]] <- sig[[p]]

  # ---- global sigma per estimator (Wassmann) --------------------------------
  glob <- list()
  for (p in names(sig)) {
    glob[[p]] <- if (is.null(sigma.global)) stats::median(sig[[p]], na.rm = TRUE)
                 else if (length(sigma.global) == 1 && is.null(names(sigma.global))) as.numeric(sigma.global)
                 else if (p %in% names(sigma.global)) as.numeric(sigma.global[[p]])
                 else stats::median(sig[[p]], na.rm = TRUE)
  }

  # ---- MDF columns ----------------------------------------------------------
  if ("MDF" %in% names(fr) && !"MDF_goflux_orig" %in% names(fr)) fr$MDF_goflux_orig <- fr$MDF
  flux <- fr[[flux.col]]
  lab <- function(c) sprintf("%02d", round(100 * c))
  df <- pmax(fr$nb.obs - 2, 1)
  add <- function(name, value) {
    fr[[name]] <<- abs(value)
    fr[[sub("^MDF_", "below_", name)]] <<- ifelse(is.na(value), NA, abs(flux) < abs(value))
  }
  for (m in mdf) for (p in names(sig)) {
    s <- sig[[p]]
    if (m %in% c("goflux", "two_sigma")) {
      add(paste("MDF", m, p, sep = "_"), if (m == "goflux") s / t * fr$flux.term else 2 * s / t * fr$flux.term)
    } else for (c in conf) {
      v <- if (m == "wassmann") stats::qnorm(c) * glob[[p]] / t * fr$flux.term
           else chr.mult * s * stats::qt(c, df = df) / t * fr$flux.term
      add(paste("MDF", m, p, lab(c), sep = "_"), v)
    }
  }

  # ---- primary --------------------------------------------------------------
  if (is.null(primary$precision) || is.null(primary$mdf)) stop("'primary' needs 'precision' and 'mdf'")
  pcol <- if (primary$mdf %in% c("goflux", "two_sigma")) paste("MDF", primary$mdf, primary$precision, sep = "_")
          else {
            if (is.null(primary$conf)) stop("'primary$conf' is needed for ", primary$mdf)
            paste("MDF", primary$mdf, primary$precision, lab(primary$conf), sep = "_")
          }
  if (!pcol %in% names(fr))
    stop("primary combination '", pcol, "' was not computed; check 'precision', 'mdf' and 'conf'")
  m <- fr[[pcol]]
  fr$MDF <- m
  fr$below_MDF <- ifelse(is.na(m), NA, abs(flux) < m)
  fr$MDF_method <- sub("^MDF_", "", pcol)
  fr$detected <- ifelse(is.na(m), NA, abs(flux) > m)
  fr$det_class <- ifelse(is.na(m), NA_character_,
                         ifelse(flux > m, "emission", ifelse(flux < -m, "uptake", "below detection")))
  attr(fr, "sigma_global") <- unlist(glob)
  fr
}
