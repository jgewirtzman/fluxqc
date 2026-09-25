#' Empirical instrument precision from first differences (Allan deviation at tau = 1 sample)
#'
#' During a chamber closure the concentration trends monotonically; the trend
#' *is* the flux signal, so the raw standard deviation of the record is useless
#' as a noise estimate. First-differencing removes any smooth trend:
#' \eqn{\delta_i = x_{i+1} - x_i} is (nearly) constant for a linear rise, while
#' white instrument noise \eqn{\epsilon} does not cancel,
#' \eqn{\mathrm{Var}(\epsilon_{i+1} - \epsilon_i) = 2\sigma^2}. Hence
#'
#' \deqn{\sigma = \mathrm{SD}(\Delta x) / \sqrt{2}}
#'
#' recovers the 1-sigma instrument noise at the logging interval regardless of
#' the trend shape. This is the Allan deviation at tau = one sample (1 s at 1 Hz).
#' It is the quantity goFlux's `prec` argument should carry when datasheet
#' precision is not trusted.
#'
#' @param x numeric concentration series (one measurement window, one gas),
#'   in the units the precision should be reported in.
#' @return numeric scalar; `NA` if fewer than three non-missing values.
#' @seealso [precision_mad()] for a robust variant, [precision_rolling()] for
#'   a whole-record scan, [precision_by_id()] to apply per `UniqueID`.
#' @examples
#' set.seed(2)
#' x <- 2000 + 0.05 * (1:600) + rnorm(600, 0, 1.2)   # rising CH4 with 1.2 ppb noise
#' precision_allan(x)   # ~1.2
#' sd(x)                # dominated by the trend, ~9
#' @export
precision_allan <- function(x) {
  x <- x[!is.na(x)]
  d <- diff(x)
  if (length(d) < 2) return(NA_real_)
  stats::sd(d) / sqrt(2)
}

#' Robust precision from the MAD of first differences
#'
#' Same idea as [precision_allan()] but with the median absolute deviation in
#' place of the standard deviation:
#'
#' \deqn{\sigma = \mathrm{MAD}(\Delta x) / \sqrt{2}}
#'
#' where `MAD` is scaled by 1.4826 so that it estimates the standard deviation
#' for Gaussian noise. Because it is built on medians it is unaffected until
#' more than half the differences are contaminated, so it can be applied to an
#' *entire* field record (moves, breaths, chamber swaps and all) rather than to
#' hand-selected quiet periods: a constant rate of change shifts every
#' difference equally and cancels, and isolated large jumps are ignored. In the
#' source campaigns the whole-record MAD estimate agreed to within 15 percent
#' with the in-closure Allan estimate from disjoint data.
#'
#' @inheritParams precision_allan
#' @param constant scale factor passed to [stats::mad()] (1.4826 = Gaussian
#'   consistency).
#' @return numeric scalar; `NA` if fewer than three non-missing values.
#' @examples
#' set.seed(3)
#' x <- 2000 + rnorm(3000, 0, 1.2)
#' x[sample(3000, 200)] <- x[sample(3000, 200)] + 50   # 7% gross outliers
#' precision_allan(x)   # inflated
#' precision_mad(x)     # ~1.2
#' @export
precision_mad <- function(x, constant = 1.4826) {
  x <- x[!is.na(x)]
  d <- diff(x)
  if (length(d) < 2) return(NA_real_)
  stats::mad(d, constant = constant) / sqrt(2)
}

#' Precision from the quietest windows of a continuous record
#'
#' Scans a continuous concentration record with a rolling window, fits a line
#' in each window, scores every window by normalised `|slope| + residual SD`,
#' and reports the median residual SD of the lowest-scoring `quantile.cutoff`
#' fraction of windows (flat *and* quiet). This gives a global per-instrument
#' precision from field data without a reference gas.
#'
#' Caveats (learned the hard way, see the source repositories): a 30-sample
#' window sits past the Allan minimum of these analyzers (a few seconds), so
#' it can capture real concentration change, and selecting an extreme quantile
#' biases the estimate low by roughly 20 percent even on pure white noise.
#' [precision_mad()] on the whole record is usually a better default; use this
#' function for comparison.
#'
#' @param x numeric concentration series (continuous record).
#' @param window window length in samples (equals seconds at 1 Hz).
#' @param quantile.cutoff fraction of windows retained (lowest scores).
#' @return list with `sd` (the precision estimate), `n_windows` and
#'   `n_selected`.
#' @examples
#' set.seed(4)
#' x <- c(2000 + rnorm(600, 0, 1), 2000 + 0.5 * (1:300) + rnorm(300, 0, 1))
#' precision_rolling(x, window = 30)
#' @export
precision_rolling <- function(x, window = 30, quantile.cutoff = 0.05) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < window + 1) {
    warning("series too short for the rolling window (", n, " points)")
    return(list(sd = NA_real_, n_windows = 0L, n_selected = 0L))
  }
  n_windows <- n - window + 1L
  tw <- seq_len(window) - 1
  X <- cbind(1, tw)
  slopes <- numeric(n_windows); resid_sd <- numeric(n_windows)
  for (i in seq_len(n_windows)) {
    xw <- x[i:(i + window - 1L)]
    if (anyNA(xw)) { slopes[i] <- NA; resid_sd[i] <- NA; next }
    fit <- stats::.lm.fit(X, xw)
    slopes[i] <- abs(fit$coefficients[2])
    resid_sd[i] <- stats::sd(fit$residuals)
  }
  norm01 <- function(v) { r <- range(v, na.rm = TRUE); if (diff(r) == 0) v * 0 else (v - r[1]) / diff(r) }
  score <- norm01(slopes) + norm01(resid_sd)
  thr <- stats::quantile(score, probs = quantile.cutoff, na.rm = TRUE)
  sel <- which(score <= thr)
  list(sd = stats::median(resid_sd[sel]), n_windows = n_windows, n_selected = length(sel))
}

#' Per-measurement precision for every UniqueID in a manID data frame
#'
#' Applies [precision_allan()] or [precision_mad()] to the flagged rows of
#' each measurement in a click.peak2-style data frame.
#'
#' @param dataframe data.frame with `UniqueID`, `flag` and the gas column
#'   (output of [goFlux::click.peak2()], [click_peak2_stacked()],
#'   [auto_id_rise()] or [windows_from_table()]).
#' @param gastype gas column, e.g. `"CH4dry_ppb"`.
#' @param method `"allan"` (default) or `"mad"`.
#' @param flag.only use only rows with `flag == 1` (default); set `FALSE` to
#'   use the whole window including shoulders.
#' @return data.frame with `UniqueID`, `sigma` (same units as `gastype`),
#'   `n` (observations used) and `method`.
#' @examples
#' d <- example_closure()
#' precision_by_id(d, "CH4dry_ppb")
#' precision_by_id(d, "CO2dry_ppm", method = "mad")
#' @export
precision_by_id <- function(dataframe, gastype, method = c("allan", "mad"),
                            flag.only = TRUE) {
  method <- match.arg(method)
  if (!gastype %in% names(dataframe)) stop("column '", gastype, "' not found")
  d <- dataframe
  if (flag.only) {
    if (!"flag" %in% names(d)) stop("'dataframe' has no 'flag' column; set flag.only = FALSE")
    d <- d[!is.na(d$flag) & d$flag == 1, ]
  }
  f <- if (method == "allan") precision_allan else precision_mad
  sp <- split(d[[gastype]], d$UniqueID)
  out <- data.frame(UniqueID = names(sp),
                    sigma = vapply(sp, f, numeric(1)),
                    n = vapply(sp, function(v) sum(!is.na(v)), integer(1)),
                    method = method, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Robust precision per run of constant logging interval
#'
#' [precision_mad()] over a whole record assumes one logging interval. An
#' analyzer whose interval was changed mid-record (1 s to 10 s, say) would mix
#' two noise scales in one MAD, and the result would describe neither. This
#' splits the record wherever the spacing of `posix` changes by more than
#' `tol` (a fraction), computes [precision_mad()] on each run, merges runs
#' with the same interval (gaps between files do not count as a change), and
#' returns one row per interval. A warning is given when more than one
#' interval is present. Intervals are rounded to the nearest 0.5 s.
#'
#' @param x numeric concentration series.
#' @param posix POSIXct timestamps of `x`.
#' @param constant scale factor passed to [stats::mad()].
#' @param tol fractional change in spacing that starts a new run.
#' @return data.frame with `dt_s` (logging interval, s), `n` (samples) and
#'   `sigma` (the MAD precision at that interval).
#' @examples
#' t0 <- as.POSIXct("2024-06-01", tz = "UTC")
#' p <- c(t0 + 0:599, t0 + 600 + 10 * (0:299))
#' x <- 2000 + c(rnorm(600, 0, 1), rnorm(300, 0, 0.4))
#' precision_mad_runs(x, p)
#' @export
precision_mad_runs <- function(x, posix, constant = 1.4826, tol = 0.2) {
  ok <- is.finite(x) & !is.na(posix); x <- x[ok]; posix <- posix[ok]
  d <- diff(as.numeric(posix)); dx <- diff(x)
  valid <- is.finite(d) & d > 0                     # duplicates / backwards jumps are not intervals
  if (sum(valid) < 3) return(data.frame(dt_s = NA_real_, n = length(x), sigma = NA_real_))
  dv <- d[valid]
  # a run breaks where the spacing changes by more than tol (fraction)
  brk <- c(TRUE, abs(diff(dv)) > tol * pmax(dv[-1], 1e-9))
  run <- cumsum(brk)
  keep <- run %in% which(tabulate(run) >= 3)        # a lone gap between files is not a run
  idx <- which(valid)[keep]
  if (length(idx) < 3) return(data.frame(dt_s = NA_real_, n = length(x), sigma = NA_real_))
  dt_r <- round(d[idx] * 2) / 2                     # nearest 0.5 s
  # one pooled MAD per interval (runs with the same interval are pooled, not averaged)
  out <- do.call(rbind, lapply(split(idx, dt_r), function(i)
    if (length(i) < 3) NULL else                      # an interval seen < 3 times is a gap, not a run
    data.frame(dt_s = round(d[i[1]] * 2) / 2, n = length(i) + 1,
               sigma = stats::mad(dx[i], constant = constant) / sqrt(2))))
  if (is.null(out)) return(data.frame(dt_s = NA_real_, n = length(x), sigma = NA_real_))
  rownames(out) <- NULL
  if (nrow(out) > 1)
    warning("record has ", nrow(out), " logging intervals (", paste(out$dt_s, collapse = ", "),
            " s); precision is reported per interval")
  out
}

#' Datasheet precision rescaled to a logging interval
#'
#' Manufacturers quote precision at 1 s. A record logged every `dt_s` seconds
#' averages that noise down, so the fair reference for an empirical precision
#' measured on such a record is `spec_1s / sqrt(dt_s)` (white noise). Compare
#' [precision_mad()] of a 10 s record with this, not with the 1 s figure;
#' otherwise a unit can appear to beat its datasheet when it is several times
#' worse. Not valid once `dt_s` exceeds the analyzer's Allan minimum.
#'
#' @param spec_1s datasheet precision at 1 s.
#' @param dt_s logging interval in seconds.
#' @return numeric, same units as `spec_1s`.
#' @examples
#' spec_at_interval(0.9, c(1, 5, 10))
#' @export
spec_at_interval <- function(spec_1s, dt_s) spec_1s / sqrt(dt_s)
