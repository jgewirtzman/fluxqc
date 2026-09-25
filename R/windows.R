#' Flag an observation window the way click.peak2 does
#'
#' Adds the five columns that [goFlux::click.peak2()] adds to an observation
#' window (`flag`, `Etime`, `start.time_corr`, `end.time_corr`,
#' `obs.length_corr`) given explicit start and end times. This is the common
#' back end of [click_peak2_stacked()], [auto_id_rise()] and
#' [windows_from_table()].
#'
#' @param d data.frame with a `POSIX.time` column (one element of the list
#'   returned by [goFlux::obs.win()]).
#' @param start,end POSIXct; the corrected window limits (inclusive).
#' @return `d` with the columns `flag` (1 inside the window, 0 outside),
#'   `Etime` (seconds since `start`), `start.time_corr`, `end.time_corr` and
#'   `obs.length_corr` (seconds).
#' @examples
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' d <- data.frame(POSIX.time = t0 + 0:299, CO2dry_ppm = 400 + 0.2 * (0:299))
#' w <- flag_window(d, t0 + 30, t0 + 200)
#' table(w$flag); range(w$Etime[w$flag == 1])
#' @export
flag_window <- function(d, start, end) {
  if (!"POSIX.time" %in% names(d)) stop("'d' must contain a 'POSIX.time' column")
  if (length(start) != 1 || length(end) != 1) stop("'start' and 'end' must be length 1")
  if (end <= start) stop("'end' must be after 'start'")
  inside <- d$POSIX.time >= start & d$POSIX.time <= end
  d$flag <- as.numeric(inside)
  d$Etime <- as.numeric(d$POSIX.time - start, units = "secs")
  d$start.time_corr <- start
  d$end.time_corr <- end
  d$obs.length_corr <- as.numeric(end - start, units = "secs")
  d
}

#' Build a manID object from a table of explicit window times
#'
#' When closure windows are known from notebooks, screenshots or a previous
#' processing round, this builds the same object [goFlux::click.peak2()] would
#' produce, without any clicking. Windows are in instrument time.
#'
#' @param ow.list list of data.frames from [goFlux::obs.win()].
#' @param windows data.frame with columns `UniqueID`, `start` and `end`.
#'   `start`/`end` may be POSIXct or character strings parseable by
#'   [as.POSIXct()] in time zone `tz`.
#' @param tz time zone used to parse character `start`/`end`.
#' @param warn.length windows with fewer flagged observations trigger a warning.
#' @return data.frame; all windows in `ow.list` whose `UniqueID` appears in
#'   `windows`, bound together with the click.peak2 columns added. Windows
#'   without a matching row in `windows` are dropped with a message.
#' @seealso [flag_window()], [auto_id_rise()], [click_peak2_stacked()]
#' @examples
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' ow <- list(data.frame(UniqueID = "A", POSIX.time = t0 + 0:400,
#'                       CO2dry_ppm = 400 + 0.3 * (0:400)))
#' win <- data.frame(UniqueID = "A", start = "2024-06-01 10:01:00",
#'                   end = "2024-06-01 10:05:00")
#' manID <- windows_from_table(ow, win)
#' sum(manID$flag); unique(manID$obs.length_corr)
#' @export
windows_from_table <- function(ow.list, windows, tz = "UTC", warn.length = 60) {
  req <- c("UniqueID", "start", "end")
  if (!all(req %in% names(windows)))
    stop("'windows' must have columns: ", paste(req, collapse = ", "))
  if (!inherits(windows$start, "POSIXct")) windows$start <- as.POSIXct(windows$start, tz = tz)
  if (!inherits(windows$end, "POSIXct"))   windows$end   <- as.POSIXct(windows$end, tz = tz)
  if (any(is.na(windows$start) | is.na(windows$end)))
    stop("could not parse some 'start'/'end' values")
  if (any(duplicated(windows$UniqueID))) stop("'windows$UniqueID' must be unique")
  if (is.data.frame(ow.list)) ow.list <- split(ow.list, ow.list$UniqueID)

  out <- vector("list", length(ow.list))
  for (k in seq_along(ow.list)) {
    d <- ow.list[[k]]
    uid <- unique(d$UniqueID)
    w <- windows[windows$UniqueID == uid, ]
    if (nrow(w) == 0) { message("No window for ", uid, "; dropped"); next }
    d2 <- flag_window(d, w$start[1], w$end[1])
    n <- sum(d2$flag)
    if (n < warn.length)
      warning(uid, ": only ", n, " observations in window (< warn.length = ",
              warn.length, ")", call. = FALSE)
    out[[k]] <- d2
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.frame())
  res <- do.call(rbind, lapply(out, as.data.frame))
  rownames(res) <- NULL
  res
}

#' Find the longest sustained CO2 rise in a trace
#'
#' Scans a concentration series for the longest stretch during which the
#' chamber was evidently closed: CO2 rises by more than `rise.ppm` over
#' `rise.secs` and then keeps climbing (never dropping more than `drop.ppm`
#' below its running maximum) with no gaps longer than `gap.secs`.
#'
#' @param t POSIXct (or numeric seconds) time stamps.
#' @param conc numeric concentration (ppm for CO2).
#' @param rise.ppm minimum increase over `rise.secs` to start a candidate rise.
#' @param rise.secs look-ahead interval for the start test, seconds.
#' @param drop.ppm tolerated drop below the running maximum before the rise ends.
#' @param min.dur minimum duration of an accepted rise, seconds.
#' @param gap.secs a gap between consecutive samples longer than this ends the rise.
#' @param conc.range values outside this range are discarded before scanning.
#' @param min.n minimum number of valid samples required to attempt a scan.
#' @return `NULL` if no rise is found; otherwise a list with `start`, `end`
#'   (same class as `t`), `dur` (seconds) and `dconc` (rise in concentration
#'   units).
#' @examples
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' tt <- t0 + 0:900
#' co2 <- c(rep(410, 200), 410 + 0.3 * (0:400), rep(530, 300)) + rnorm(901, 0, 0.3)
#' find_rise(tt, co2)
#' @export
find_rise <- function(t, conc, rise.ppm = 6, rise.secs = 60, drop.ppm = 8,
                      min.dur = 90, gap.secs = 5, conc.range = c(300, 20000),
                      min.n = 120) {
  ok <- !is.na(conc) & !is.na(t) & conc > conc.range[1] & conc < conc.range[2]
  t <- t[ok]; conc <- conc[ok]
  n <- length(t)
  if (n < min.n) return(NULL)
  ts <- as.numeric(t)
  best <- NULL
  i <- 1L
  while (i <= n) {
    k <- findInterval(ts[i] + rise.secs, ts)   # last index with time <= t_i + rise.secs
    if (k <= i) k <- min(i + 1L, n)
    if (k >= n && ts[k] - ts[i] < rise.secs) break
    if (conc[k] - conc[i] > rise.ppm) {
      j <- i; peak <- conc[i]
      while (j < n && conc[j + 1] >= peak - drop.ppm &&
             (ts[j + 1] - ts[j]) < gap.secs) {
        j <- j + 1L; peak <- max(peak, conc[j])
      }
      dur <- ts[j] - ts[i]
      if (dur >= min.dur && (is.null(best) || dur > best$dur))
        best <- list(start = t[i], end = t[j], dur = dur, dconc = conc[j] - conc[i])
      i <- j + 10L
    } else i <- i + 5L
  }
  best
}

#' Non-interactive window selection from the sustained CO2 rise
#'
#' A stand-in for [goFlux::click.peak2()] that can run under `Rscript`. For
#' every observation window it (1) finds the longest sustained rise of
#' `detect.gas` with [find_rise()] (the chamber is closed, so CO2 accumulates
#' from respiration), (2) trims the first `trim.secs` (mixing) and caps the
#' window at the notebook end time plus `slack.secs` or at `max.secs`,
#' whichever is shorter, and (3) falls back to notebook start + `trim.secs`
#' .. notebook end (or `fallback.secs`) when no rise is found. It writes one
#' diagnostic PNG per measurement, with one stacked panel per gas in `gases`,
#' so the automatic windows can be reviewed and redone interactively where
#' needed.
#'
#' @inheritParams windows_from_table
#' @param detect.gas column used to find the rise (default `"CO2dry_ppm"`).
#' @param gases gas columns drawn in the diagnostic plot, one stacked panel
#'   each; `"auto"` (default) uses every recognised goFlux gas column present.
#' @param end.times optional notebook end times: a data.frame with columns
#'   `UniqueID` and `end.time` (POSIXct or character in `tz`), or a named
#'   POSIXct vector. Missing entries are allowed.
#' @param trim.secs seconds removed from the start of the detected rise.
#' @param max.secs maximum window length in seconds.
#' @param min.secs windows shorter than this are extended to `min.secs + 30`
#'   and tagged `"short"` in the log.
#' @param slack.secs seconds allowed past the notebook end time.
#' @param fallback.secs window length used when neither a rise nor a notebook
#'   end time is available.
#' @param plot.dir directory for diagnostic PNGs; `NULL` disables plotting.
#' @param ... further arguments passed to [find_rise()].
#' @return data.frame in click.peak2 format (all windows bound together), with
#'   an attribute `"log"` holding a data.frame (`UniqueID`, `method`,
#'   `auto_start`, `auto_end`, `n_flag`, `dCO2`) describing what was chosen.
#' @seealso [find_rise()], [flag_window()], [click_peak2_stacked()]
#' @examples
#' t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
#' tt <- t0 + 0:900
#' set.seed(1)
#' ow <- list(data.frame(UniqueID = "tree1", POSIX.time = tt,
#'   CO2dry_ppm = c(rep(410, 200), 410 + 0.3 * (0:400), rep(530, 300)) + rnorm(901, 0, 0.3),
#'   CH4dry_ppb = 2000 + 0.01 * (0:900) + rnorm(901, 0, 0.5),
#'   start.time = t0 + 180, obs.length = 600))
#' manID <- auto_id_rise(ow, plot.dir = NULL)
#' attr(manID, "log")
#' @export
auto_id_rise <- function(ow.list, detect.gas = "CO2dry_ppm", gases = "auto",
                         end.times = NULL, trim.secs = 30, max.secs = 600,
                         min.secs = 60, slack.secs = 60, fallback.secs = 300,
                         plot.dir = NULL, tz = "UTC", warn.length = 60, ...) {
  if (is.data.frame(ow.list)) ow.list <- split(ow.list, ow.list$UniqueID)
  en <- NULL
  if (!is.null(end.times)) {
    if (is.data.frame(end.times)) {
      e <- end.times$end.time
      if (!inherits(e, "POSIXct")) e <- as.POSIXct(e, tz = tz)
      en <- e; names(en) <- end.times$UniqueID
    } else {
      en <- end.times
      if (!inherits(en, "POSIXct")) en <- as.POSIXct(en, tz = tz)
    }
  }
  if (!is.null(plot.dir)) dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

  out <- vector("list", length(ow.list)); log <- vector("list", length(ow.list))
  for (k in seq_along(ow.list)) {
    d <- ow.list[[k]]; uid <- unique(d$UniqueID)[1]
    st <- unique(d$start.time)[1]
    if (is.null(st) || is.na(st)) stop(uid, ": no 'start.time' column in window")
    if (!detect.gas %in% names(d)) stop("'detect.gas' column '", detect.gas, "' not found")
    en_uid <- if (!is.null(en) && uid %in% names(en) && !is.na(en[[uid]])) en[[uid]] else NA
    r <- find_rise(d$POSIX.time, d[[detect.gas]], ...)
    method <- "auto_rise"
    if (is.null(r)) {
      s0 <- st + trim.secs
      e0 <- if (!is.na(en_uid)) en_uid else st + fallback.secs
      method <- "notebook_fallback"
    } else {
      s0 <- r$start + trim.secs
      e0 <- r$end
      if (!is.na(en_uid)) e0 <- min(e0, en_uid + slack.secs)
      e0 <- min(e0, s0 + max.secs)
    }
    if (as.numeric(e0 - s0, units = "secs") < min.secs) {
      e0 <- s0 + min.secs + 30; method <- paste(method, "short")
    }
    d2 <- flag_window(d, s0, e0)
    n_flag <- sum(d2$flag)
    if (n_flag < warn.length)
      warning(uid, ": only ", n_flag, " observations in window (< warn.length = ",
              warn.length, ")", call. = FALSE)
    out[[k]] <- d2
    log[[k]] <- data.frame(UniqueID = uid, method = method,
                           auto_start = format(s0, "%Y-%m-%d %H:%M:%S"),
                           auto_end = format(e0, "%Y-%m-%d %H:%M:%S"),
                           n_flag = n_flag,
                           dCO2 = if (is.null(r)) NA_real_ else round(r$dconc, 1),
                           stringsAsFactors = FALSE)
    if (!is.null(plot.dir)) {
      fl <- d2$flag == 1
      g <- resolve_gases(gases, d)
      png(file.path(plot.dir, paste0(sprintf("%03d", k), "_", uid, ".png")), 900, 350 * length(g))
      op <- par(mfrow = c(length(g), 1), mar = c(3, 4.5, 2.5, 1))
      for (i in seq_along(g)) {
        plot(d$POSIX.time, d[[g[i]]], col = ifelse(fl, "red", "grey40"), pch = 16,
             cex = 0.5, ylab = g[i], xlab = "",
             main = if (i == 1) paste(uid, "-", method) else "")
        abline(v = st, col = "blue")
        if (!is.na(en_uid)) abline(v = en_uid, col = "blue", lty = 2)
      }
      par(op); dev.off()
    }
  }
  res <- do.call(rbind, lapply(out, as.data.frame)); rownames(res) <- NULL
  attr(res, "log") <- do.call(rbind, log)
  res
}
