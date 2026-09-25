#' Minimum detectable flux by several published methods
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
#' @seealso [flag_detection()] to append these and detection flags to a
#'   `best.flux` output.
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

#' Append multi-method MDFs and below-detection flags to a best.flux output
#'
#' Computes [mdf_multi()] for every row of a [goFlux::best.flux()] (or
#' [goFlux::goFlux()]) result and appends the MDF columns, a logical
#' `below_<method>` flag per method, and a summary classification based on
#' one `primary` method: `detected` (logical) and `det_class`
#' (`"emission"`, `"uptake"` or `"below detection"`).
#'
#' @param flux.result data.frame from [goFlux::best.flux()] or
#'   [goFlux::goFlux()]; must contain `UniqueID`, `flux.term`, `nb.obs` and the
#'   flux column.
#' @param sigma per-measurement precision: either a data.frame with
#'   `UniqueID` and `sigma` (output of [precision_by_id()]), a numeric vector
#'   aligned with the rows of `flux.result`, or the name of a column in
#'   `flux.result`. `NULL` disables the per-measurement methods.
#' @param sigma.global scalar global precision for the Wassmann variants;
#'   default is the median of `sigma` when that is supplied.
#' @param t closure length in seconds; default `nb.obs` (1 Hz).
#' @param prec manufacturer precision; default the `prec` column of
#'   `flux.result`.
#' @param flux.col flux column to test (`"best.flux"`, or `"LM.flux"` /
#'   `"HM.flux"` on a raw `goFlux()` output).
#' @param primary method used for `detected` / `det_class`: one of
#'   `"2sigma"`, `"goflux"`, `"wass90"`, `"wass95"`, `"wass99"`, `"chr90"`,
#'   `"chr95"`, `"chr99"`. Default `"2sigma"` when `sigma` is supplied,
#'   otherwise `"goflux"`.
#' @inheritParams mdf_multi
#' @return `flux.result` with the added columns.
#' @examples
#' d <- example_closure()
#' res <- goFlux::goFlux(d, "CH4dry_ppb")
#' best <- goFlux::best.flux(res)
#' sig <- precision_by_id(d, "CH4dry_ppb")
#' out <- flag_detection(best, sigma = sig)
#' out[, c("UniqueID", "best.flux", "MDF_goflux", "MDF_2sigma", "MDF_chr95",
#'         "below_2sigma", "det_class")]
#' @export
flag_detection <- function(flux.result, sigma = NULL, sigma.global = NULL, t = NULL,
                           prec = NULL, flux.col = "best.flux", primary = NULL,
                           conf = c(0.90, 0.95, 0.99), chr.mult = 3) {
  fr <- flux.result
  for (col in c("UniqueID", "flux.term", "nb.obs", flux.col))
    if (!col %in% names(fr)) stop("'flux.result' has no column '", col, "'")
  n <- nrow(fr)
  # per-measurement sigma
  sig <- NULL
  if (!is.null(sigma)) {
    if (is.data.frame(sigma)) {
      if (!all(c("UniqueID", "sigma") %in% names(sigma)))
        stop("'sigma' data.frame needs columns 'UniqueID' and 'sigma'")
      sig <- sigma$sigma[match(fr$UniqueID, sigma$UniqueID)]
    } else if (is.character(sigma) && length(sigma) == 1 && sigma %in% names(fr)) {
      sig <- fr[[sigma]]
    } else if (is.numeric(sigma)) {
      sig <- rep_len(sigma, n)
    } else stop("'sigma' must be a data.frame, a numeric vector or a column name")
  }
  if (is.null(sigma.global) && !is.null(sig)) sigma.global <- stats::median(sig, na.rm = TRUE)
  if (is.null(t)) t <- fr$nb.obs
  if (is.null(prec)) prec <- if ("prec" %in% names(fr)) fr$prec else NULL
  m <- mdf_multi(flux.term = fr$flux.term, t = t, prec = prec, sigma = sig,
                 sigma.global = sigma.global, nb.obs = fr$nb.obs, conf = conf,
                 chr.mult = chr.mult)
  if (!is.null(sig)) fr$sigma <- sig
  fr <- cbind(fr, m)
  flux <- fr[[flux.col]]
  for (col in names(m)) {
    fr[[sub("^MDF_", "below_", col)]] <- ifelse(is.na(m[[col]]), NA, abs(flux) < m[[col]])
  }
  if (is.null(primary)) primary <- if (!is.null(sig)) "2sigma" else "goflux"
  pcol <- paste0("MDF_", primary)
  if (!pcol %in% names(m)) stop("unknown 'primary' method: ", primary)
  mdf <- fr[[pcol]]
  fr$detected <- ifelse(is.na(mdf), NA, abs(flux) > mdf)
  fr$det_class <- ifelse(is.na(mdf), NA_character_,
                         ifelse(flux > mdf, "emission",
                                ifelse(flux < -mdf, "uptake", "below detection")))
  fr$det_method <- primary
  fr
}
