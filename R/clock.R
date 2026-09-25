#' Estimate the constant offset between notebook closure times and analyzer time
#'
#' Field notebooks record closure start times on a watch or phone; the
#' analyzer stamps data with its own clock, which may be minutes off (time
#' zone slips, drift, a reset). This function finds the constant offset
#' `delta` such that `notebook + delta` lines up with the concentration-rise
#' onsets in the analyzer record. For every candidate offset in `search`
#' (seconds, stepped by `step`) it computes, for each notebook start, an
#' onset score at `start + delta`: the rise from the 7 s around that instant
#' to the mean of the `window` seconds *after* it, minus the absolute change
#' from the `window` seconds *before* it (large only where a rise begins from
#' a flat baseline; near zero in the middle of a steady rise; strongly
#' negative at the drop when the chamber is lifted), and scores the offset by
#' the median across closures (a cross-correlation of the notebook start
#' impulses with the rise onsets). A sealed chamber always makes CO2 start
#' rising at closure, so the score peaks at the true offset. The score curve
#' is returned so a flat or multi-modal curve can be recognised. This is a
#' non-interactive replacement for the slider-based Shiny alignment app used
#' in the BlueFlux rescue (`offset_finder_app.R`), which relied on visual
#' alignment. It works with any goFlux importer's output.
#'
#' Tools for a rarer problem, analyzer clock *resets* that interleave two
#' time streams in one raw LGR file, are shipped as unexported scripts in
#' `system.file("rescue", package = "stemflux")`.
#'
#' @param traces data.frame with `POSIX.time` and the `gas` column (a whole
#'   imported record, e.g. from [goFlux::import.UGGA()]).
#' @param notebook_times POSIXct vector of notebook closure start times, or a
#'   data.frame with a `start.time` column (an auxfile).
#' @param gas concentration column used for the alignment (CO2 recommended).
#' @param search numeric length 2; offset range searched, seconds (analyzer
#'   time minus notebook time).
#' @param step search step, seconds.
#' @param window length of the before/after averaging windows, seconds.
#' @param plot draw the score curve and the aligned traces.
#' @param auxfile optional auxfile whose `start.time` column is shifted by the
#'   best offset (original kept as `start.time_notebook`).
#' @return list with `offset` (seconds to *add* to notebook times), `score`
#'   (onset score at the best offset, gas units), `scores` (data.frame of
#'   `offset`, `score` for the whole curve), `n_closures` used, and, if
#'   `auxfile` was given, `auxfile` with corrected `start.time`.
#' @examples
#' set.seed(7)
#' t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")
#' tt <- t0 + 0:3599
#' co2 <- 410 + rnorm(3600, 0, 0.5)
#' true_starts <- t0 + c(300, 1200, 2100, 3000)          # analyzer time
#' for (s in as.numeric(difftime(true_starts, t0, units = "secs")))
#'   co2[(s + 1):(s + 300)] <- co2[(s + 1):(s + 300)] + 0.3 * (1:300)
#' traces <- data.frame(POSIX.time = tt, CO2dry_ppm = co2)
#' notebook <- true_starts - 137                          # watch was 137 s behind
#' fc <- find_clock_offset(traces, notebook, search = c(-600, 600), plot = FALSE)
#' fc$offset
#' @export
find_clock_offset <- function(traces, notebook_times, gas = "CO2dry_ppm",
                              search = c(-1800, 1800), step = 1, window = 60,
                              plot = TRUE, auxfile = NULL) {
  if (!all(c("POSIX.time", gas) %in% names(traces))) stop("'traces' needs 'POSIX.time' and '", gas, "'")
  nb <- if (is.data.frame(notebook_times)) notebook_times$start.time else notebook_times
  if (!inherits(nb, "POSIXct")) stop("'notebook_times' must be POSIXct (or a data.frame with a POSIXct 'start.time')")
  nb <- nb[!is.na(nb)]
  if (!length(nb)) stop("no notebook times")
  tt <- as.numeric(traces$POSIX.time); y <- as.numeric(traces[[gas]])
  ok <- !is.na(tt) & !is.na(y); tt <- tt[ok]; y <- y[ok]
  o <- order(tt); tt <- tt[o]; y <- y[o]
  t0 <- floor(min(tt)); t1 <- ceiling(max(tt))
  grid <- seq(t0, t1, by = 1)
  yg <- stats::approx(tt, y, xout = grid, rule = 1, ties = mean)$y
  has <- !is.na(yg)
  cs <- cumsum(ifelse(has, yg, 0)); cn <- cumsum(has)
  N <- length(grid)
  wmean <- function(a, b) {                                # mean of yg over grid indices a..b
    a <- pmax(a, 1); b <- pmin(b, N)
    bad <- a > b
    a2 <- pmin(a, N); b2 <- pmax(b, 1)
    s <- cs[b2] - ifelse(a2 > 1, cs[pmax(a2 - 1, 1)], 0); k <- cn[b2] - ifelse(a2 > 1, cn[pmax(a2 - 1, 1)], 0)
    r <- ifelse(k > 0, s / k, NA); r[bad] <- NA; r
  }
  offsets <- seq(search[1], search[2], by = step)
  nbn <- as.numeric(nb)
  score <- vapply(offsets, function(d) {
    idx <- round(nbn + d - t0) + 1
    # onset score: rise after the candidate start, minus any change before it;
    # large only where a rise begins from a flat baseline, ~0 mid-rise,
    # strongly negative at the drop when the chamber is lifted
    ctr <- wmean(idx - 3, idx + 3)
    onset <- (wmean(idx, idx + window) - ctr) - abs(ctr - wmean(idx - window, idx))
    if (all(is.na(onset))) NA_real_ else stats::median(onset, na.rm = TRUE)
  }, numeric(1))
  if (all(is.na(score))) stop("notebook times never fall inside the record for any offset in 'search'")
  best <- offsets[which.max(score)]
  out <- list(offset = best, score = max(score, na.rm = TRUE),
              scores = data.frame(offset = offsets, score = score), n_closures = length(nb))
  if (!is.null(auxfile)) {
    if (!"start.time" %in% names(auxfile)) stop("'auxfile' needs a 'start.time' column")
    auxfile$start.time_notebook <- auxfile$start.time
    auxfile$start.time <- auxfile$start.time + best
    out$auxfile <- auxfile
  }
  if (plot) {
    op <- par(mfrow = c(2, 1), mar = c(4, 4.5, 2.5, 1)); on.exit(par(op))
    plot(offsets, score, type = "l", xlab = "offset added to notebook time (s)",
         ylab = paste0("median onset score (", gas, ")"),
         main = sprintf("best offset = %g s (n = %d closures)", best, length(nb)))
    abline(v = best, col = "red"); abline(h = 0, col = "grey60", lty = 2)
    show <- nb[seq_len(min(6, length(nb)))]
    rel <- seq(-2 * window, 4 * window)
    ymat <- sapply(show, function(s) yg[pmin(pmax(round(as.numeric(s) + best - t0) + 1 + rel, 1), N)])
    graphics::matplot(rel, ymat, type = "l", lty = 1, xlab = "seconds from corrected notebook start",
                      ylab = gas, main = "traces aligned at notebook start + offset")
    abline(v = 0, col = "red")
  }
  out
}
