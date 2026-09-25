#' Manual window selection with stacked CO2 and CH4 panels
#'
#' A drop-in variant of [goFlux::click.peak2()] that shows two stacked panels,
#' `gastype` (top, clicked) and `gastype2` (bottom, shown), sharing one time
#' axis, so each window can be judged with both gas traces in view. You click
#' the start and the end of the measurement on the top panel; the same window
#' is applied to both gases and goFlux is then run per gas on the shared
#' `flag`. The output has exactly the structure of `click.peak2()`.
#'
#' Like `click.peak2()`, this uses [graphics::identify()] in a pop-up device
#' opened with [grDevices::dev.new()], so it must be run in an interactive
#' session with a graphics device that supports `identify()` (not the RStudio
#' plot pane). Under `Rscript` use [auto_id_rise()] or [windows_from_table()].
#'
#' @param ow.list list of data.frames from [goFlux::obs.win()].
#' @param gastype gas shown in the top panel and used for clicking.
#' @param gastype2 gas shown in the bottom panel.
#' @param seq integer vector; which elements of `ow.list` to process
#'   (default all).
#' @param ref.secs optional typical closure length; a dashed guide line is
#'   drawn at `start.time + ref.secs`.
#' @param sleep seconds the validation plot stays open (0 to skip it).
#' @param plot.lim,plot.lim2 y-axis clipping limits for the top and bottom panel.
#' @param warn.length windows with fewer flagged observations trigger a warning.
#' @param save.plots optional path of a PDF that receives all validation panels.
#' @param width,height size of the pop-up window in inches.
#' @param abline draw the recorded `start.time` as a blue line.
#' @param abline_corr kept for signature compatibility with `click.peak2()`.
#' @return data.frame; all processed windows bound together, with the columns
#'   `flag`, `Etime`, `start.time_corr`, `end.time_corr` and
#'   `obs.length_corr` added, exactly as `click.peak2()` returns. Windows where
#'   identification failed are dropped with a warning.
#' @seealso [goFlux::click.peak2()], [auto_id_rise()], [windows_from_table()]
#' @examples
#' \dontrun{
#' ow <- goFlux::obs.win(inputfile = imp, auxfile = aux, obs.length = 600, shoulder = 120)
#' manID <- click_peak2_stacked(ow, ref.secs = 600, save.plots = "windows.pdf")
#' CO2_flux <- goFlux::goFlux(manID, "CO2dry_ppm")
#' CH4_flux <- goFlux::goFlux(manID, "CH4dry_ppb")
#' }
#' @export
click_peak2_stacked <- function(ow.list, gastype = "CO2dry_ppm", gastype2 = "CH4dry_ppb",
                                seq = NULL, ref.secs = NULL, sleep = 3,
                                plot.lim = c(300, 5000), plot.lim2 = c(1500, 100000),
                                warn.length = 60, save.plots = NULL,
                                width = 7, height = 6, abline = TRUE, abline_corr = TRUE) {
  if (is.data.frame(ow.list)) ow.list <- split(ow.list, ow.list$UniqueID)
  if (is.null(seq)) seq <- seq_along(ow.list)
  ow.list.name <- deparse(substitute(ow.list))
  sleeploop <- function(x) { p <- proc.time(); Sys.sleep(x); proc.time() - p }
  default.par <- par(no.readonly = TRUE)
  on.exit(par(default.par)); on.exit(Sys.unsetenv("TZ"), add = TRUE)
  ow.corr.ls <- vector("list", length(ow.list))

  for (ow in seq) {
    d0 <- ow.list[[ow]]
    g1 <- d0[[gastype]]; g2 <- d0[[gastype2]]
    time.meas <- d0$POSIX.time
    uid <- unique(d0$UniqueID)[1]
    start.time <- unique(d0$start.time)[1]
    yl1 <- stacked_ylim(g1, plot.lim); yl2 <- stacked_ylim(g2, plot.lim2)
    Sys.setenv(TZ = attr(time.meas, "tzone"))

    rownum <- NULL; identify.error <- NULL
    tryCatch({
      dev.new(noRStudioGD = TRUE, width = width, height = height)
      scr <- draw_stacked(time.meas, g1, g2, yl1, yl2, gastype, gastype2,
                          paste0(uid, "  -  click START then END here (", gastype, ")"),
                          v_blue = if (abline) start.time else NULL,
                          v_guide = if (!is.null(ref.secs)) start.time + ref.secs else NULL,
                          time.axis = TRUE)
      screen(scr[1], new = FALSE)
      rownum <- identify(time.meas, g1, pos = FALSE, n = 2, plot = TRUE,
                         atpen = FALSE, offset = 0.5, tolerance = 0.25)
    }, error = function(e) identify.error <<- e)
    try(close.screen(all.screens = TRUE), silent = TRUE)
    if (!is.null(dev.list())) dev.off()

    if (!is.null(identify.error)) {
      warning(ow.list.name, "[[", ow, "]] ", uid, ": ", conditionMessage(identify.error),
              call. = FALSE); next
    }
    if (length(rownum) < 2) {
      warning(ow.list.name, "[[", ow, "]] ", uid, ": fewer than two points clicked; skipped",
              call. = FALSE); next
    }
    rownum <- sort(rownum)
    d <- flag_window(d0, time.meas[rownum[1]], time.meas[rownum[2]])
    ow.corr.ls[[ow]] <- d
    n_flag <- sum(d$flag == 1)
    if (n_flag < warn.length)
      warning(ow.list.name, "[[", ow, "]] ", uid, ": only ", n_flag,
              " obs (< warn.length=", warn.length, ")", call. = FALSE)
    else message(ow.list.name, "[[", ow, "]] ", uid, ": good window (", n_flag, " obs)")

    if (!is.null(sleep) && sleep > 0) {
      dev.new(noRStudioGD = TRUE, width = width, height = height)
      draw_stacked(d$Etime, g1, g2, yl1, yl2, gastype, gastype2,
                   paste0(uid, "   (red = selected window)"),
                   v_red = c(0, max(d$Etime[d$flag == 1])), col = d$flag + 1)
      sleeploop(sleep)
      try(close.screen(all.screens = TRUE), silent = TRUE)
      if (!is.null(dev.list())) dev.off()
    }
  }

  done <- ow.corr.ls[!vapply(ow.corr.ls, is.null, logical(1))]
  if (!is.null(save.plots) && length(done) > 0) {
    pdf(file = save.plots, width = width, height = height)
    for (d in done) {
      draw_stacked(d$Etime, d[[gastype]], d[[gastype2]],
                   stacked_ylim(d[[gastype]], plot.lim), stacked_ylim(d[[gastype2]], plot.lim2),
                   gastype, gastype2,
                   paste0(unique(d$UniqueID)[1], "   (red = selected window)"),
                   v_red = c(0, max(d$Etime[d$flag == 1])), col = d$flag + 1)
      close.screen(all.screens = TRUE)
    }
    dev.off()
  }
  if (length(done) == 0) return(data.frame())
  res <- do.call(rbind, lapply(done, as.data.frame)); rownames(res) <- NULL
  res
}

# y-limits like click.peak2: data range + 5% pad, clipped to lim
stacked_ylim <- function(v, lim) {
  ymax <- max(v, na.rm = TRUE); ymin <- min(v, na.rm = TRUE)
  yd <- ymax - ymin; lo <- ymin - yd * 0.05; hi <- ymax + yd * 0.05
  c(max(lo, lim[1]), min(hi, lim[2]))
}

# Draw the two stacked panels on the current device using split.screen();
# returns the screen numbers so the caller can re-activate the top screen.
draw_stacked <- function(x, y1, y2, yl1, yl2, lab1, lab2, top.title,
                         v_blue = NULL, v_red = NULL, v_guide = NULL,
                         time.axis = FALSE, col = "grey30") {
  vlines <- function() {
    if (!is.null(v_blue))  abline(v = v_blue,  col = "blue")
    if (!is.null(v_guide)) abline(v = v_guide, col = "grey55", lty = 2)
    if (!is.null(v_red))   abline(v = v_red,   col = "red")
  }
  scr <- split.screen(c(2, 1))
  screen(scr[1]); par(mar = c(2.2, 4.6, 3.0, 1))
  plot(x, y1, col = col, xlab = "", ylab = lab1, xaxt = "n", ylim = yl1)
  vlines()
  if (time.axis) axis.POSIXct(1, at = pretty(x, 12), format = "%H:%M:%S") else axis(1)
  title(main = top.title, font.main = 2)
  screen(scr[2]); par(mar = c(3.2, 4.6, 1.6, 1))
  plot(x, y2, col = col, xlab = if (time.axis) "Time" else "Etime (s)",
       ylab = lab2, xaxt = "n", ylim = yl2)
  vlines()
  if (time.axis) axis.POSIXct(1, at = pretty(x, 12), format = "%H:%M:%S") else axis(1)
  invisible(scr)
}
