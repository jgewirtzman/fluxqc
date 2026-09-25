#' Manual window selection with any number of stacked gas panels
#'
#' A drop-in variant of [goFlux::click.peak2()] that shows one panel per gas
#' column in `gases`, stacked on a shared time axis, so each window can be
#' judged with every trace in view (CO2, CH4, N2O, H2O, ... in any order). You
#' click the start and the end of the measurement on the *top* panel (the
#' first gas listed); the same window is applied to all gases and goFlux is
#' then run per gas on the shared `flag`. The output has exactly the structure
#' of `click.peak2()`.
#'
#' Like `click.peak2()`, this uses [graphics::identify()] in a pop-up device
#' opened with [grDevices::dev.new()], so it must be run in an interactive
#' session with a graphics device that supports `identify()` (not the RStudio
#' plot pane). Under `Rscript` use [auto_id_rise()] or [windows_from_table()].
#'
#' @param ow.list list of data.frames from [goFlux::obs.win()].
#' @param gases character vector of gas columns to stack, top to bottom; the
#'   first is the one clicked. `"auto"` uses every recognised goFlux gas
#'   column present (`CO2dry_ppm`, `CH4dry_ppb`, `N2Odry_ppb`, `NH3dry_ppb`,
#'   `COdry_ppb`, `H2O_ppm`, ...), in that order.
#' @param seq integer vector; which elements of `ow.list` to process
#'   (default all).
#' @param ref.secs optional typical closure length; a dashed guide line is
#'   drawn at `start.time + ref.secs`.
#' @param sleep seconds the validation plot stays open (0 to skip it).
#' @param plot.lim y-axis clipping limits: `NULL` (data range + 5 percent), a
#'   numeric `c(lo, hi)` applied to every panel, or a named list of such
#'   vectors keyed by gas column (unnamed gases are unclipped).
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
#' manID <- click_peak2_stacked(ow, gases = c("CO2dry_ppm", "CH4dry_ppb", "N2Odry_ppb"),
#'                              ref.secs = 600, save.plots = "windows.pdf")
#' manID <- click_peak2_stacked(ow, gases = "auto")
#' }
#' @export
click_peak2_stacked <- function(ow.list, gases = c("CO2dry_ppm", "CH4dry_ppb"),
                                seq = NULL, ref.secs = NULL, sleep = 3, plot.lim = NULL,
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
    g <- resolve_gases(gases, d0)
    ys <- lapply(g, function(x) d0[[x]])
    time.meas <- d0$POSIX.time
    uid <- unique(d0$UniqueID)[1]
    start.time <- unique(d0$start.time)[1]
    yl <- stacked_ylims(ys, g, plot.lim)
    Sys.setenv(TZ = attr(time.meas, "tzone"))

    rownum <- NULL; identify.error <- NULL
    tryCatch({
      dev.new(noRStudioGD = TRUE, width = width, height = height)
      scr <- draw_stacked(time.meas, ys, g, yl,
                          paste0(uid, "  -  click START then END here (", g[1], ")"),
                          v_blue = if (abline) start.time else NULL,
                          v_guide = if (!is.null(ref.secs)) start.time + ref.secs else NULL,
                          time.axis = TRUE)
      screen(scr[1], new = FALSE)
      rownum <- identify(time.meas, ys[[1]], pos = FALSE, n = 2, plot = TRUE,
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
      draw_stacked(d$Etime, ys, g, yl, paste0(uid, "   (red = selected window)"),
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
      g <- resolve_gases(gases, d); ys <- lapply(g, function(x) d[[x]])
      draw_stacked(d$Etime, ys, g, stacked_ylims(ys, g, plot.lim),
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

# gas columns goFlux knows about, in display order
goflux_gas_columns <- function() {
  c("CO2dry_ppm", "CH4dry_ppb", "N2Odry_ppb", "NH3dry_ppb", "COdry_ppb",
    "CO2_ppm", "CH4_ppb", "N2O_ppb", "H2O_ppm")
}

# Resolve gases = "auto" or a character vector against the columns of d
resolve_gases <- function(gases, d) {
  if (identical(gases, "auto")) {
    g <- intersect(goflux_gas_columns(), names(d))
    if (!length(g)) stop("no recognised gas column found in the data; give 'gases' explicitly")
    return(g)
  }
  miss <- setdiff(gases, names(d))
  if (length(miss)) stop("gas column(s) not found: ", paste(miss, collapse = ", "))
  gases
}

# y-limits like click.peak2: data range + 5% pad, clipped to lim (NULL = no clip)
stacked_ylim <- function(v, lim = NULL) {
  ymax <- max(v, na.rm = TRUE); ymin <- min(v, na.rm = TRUE)
  yd <- ymax - ymin; if (yd == 0) yd <- 1
  lo <- ymin - yd * 0.05; hi <- ymax + yd * 0.05
  if (is.null(lim)) c(lo, hi) else c(max(lo, lim[1]), min(hi, lim[2]))
}

stacked_ylims <- function(ys, gases, plot.lim) {
  lapply(seq_along(ys), function(i) {
    lim <- if (is.null(plot.lim)) NULL
           else if (is.list(plot.lim)) plot.lim[[gases[i]]]
           else plot.lim
    stacked_ylim(ys[[i]], lim)
  })
}

# Draw one panel per gas on the current device using split.screen();
# returns the screen numbers so the caller can re-activate the top screen.
draw_stacked <- function(x, ys, labs, ylims, top.title, v_blue = NULL, v_red = NULL,
                         v_guide = NULL, time.axis = FALSE, col = "grey30") {
  n <- length(ys)
  vlines <- function() {
    if (!is.null(v_blue))  abline(v = v_blue,  col = "blue")
    if (!is.null(v_guide)) abline(v = v_guide, col = "grey55", lty = 2)
    if (!is.null(v_red))   abline(v = v_red,   col = "red")
  }
  scr <- split.screen(c(n, 1))
  for (i in seq_len(n)) {
    screen(scr[i])
    par(mar = c(if (i == n) 3.2 else 2.2, 4.6, if (i == 1) 3.0 else 1.2, 1))
    plot(x, ys[[i]], col = col, xlab = "", ylab = labs[i], xaxt = "n", ylim = ylims[[i]])
    vlines()
    if (time.axis) axis.POSIXct(1, at = pretty(x, 12), format = "%H:%M:%S") else axis(1)
    if (i == 1) title(main = top.title, font.main = 2)
    if (i == n) title(xlab = if (time.axis) "Time" else "Etime (s)")
  }
  invisible(scr)
}
