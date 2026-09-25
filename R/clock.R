#' Detect analyzer clock resets from the raw write order of an LGR/UGGA file
#'
#' Some ABB/LGR (UGGA, GLA131) files contain two interleaved time streams: the
#' analyzer clock reset mid-file (typically after a power loss or GPS sync) and
#' rows written later carry *earlier* time stamps. goFlux's import sorts by
#' time, so the two streams become interleaved and the true write order is
#' lost. This function reads the raw text file in write order and labels every
#' data line as `"ahead"` (at or advancing the running maximum time stamp) or
#' `"behind"` (more than `tolerance` seconds behind the running maximum, i.e.
#' from the reset stream). It works whether the streams are sequential blocks
#' or interspersed row by row.
#'
#' @param filepath path to a raw LGR `.txt` data file (`_f` file).
#' @param tolerance seconds behind the running maximum that mark a row as
#'   belonging to the reset stream.
#' @param date.pattern regular expression identifying data lines (default:
#'   lines starting with `M/D/YYYY H:MM:SS`).
#' @param time.format format string for [as.POSIXct()] applied to the first
#'   comma-separated field.
#' @return list with `has_reset` (logical), `streams` (character vector, one
#'   per data line), `stream_labels`, `header_lines`, `data_lines` and
#'   `timestamps` (POSIXct). Use [write_stream_files()] to split the file.
#' @seealso [write_stream_files()], [separate_window_streams()],
#'   [obs_win_deinterleave()], [click_flux_2d()]
#' @examples
#' f <- tempfile(fileext = ".txt")
#' hdr <- c("SN: 123", "SysTime, [CH4]d_ppm, [CO2]d_ppm")
#' t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")
#' good <- sprintf("%s, 2.0, 400", format(t0 + 0:9, "%m/%d/%Y %H:%M:%S"))
#' reset <- sprintf("%s, 2.1, 401", format(t0 - 600 + 0:9, "%m/%d/%Y %H:%M:%S"))
#' writeLines(c(hdr, good, reset), f)
#' cs <- detect_clock_streams(f)
#' cs$has_reset; table(cs$streams)
#' @export
detect_clock_streams <- function(filepath, tolerance = 30,
                                 date.pattern = "^\\s*\\d{1,2}/\\d{1,2}/\\d{4}\\s+\\d{1,2}:\\d{2}:\\d{2}",
                                 time.format = "%m/%d/%Y %H:%M:%OS") {
  raw_lines <- readLines(filepath, warn = FALSE)
  is_data <- grepl(date.pattern, raw_lines)
  empty <- function(header, data) list(has_reset = FALSE, streams = rep("ahead", length(data)),
                                       stream_labels = "ahead", header_lines = header,
                                       data_lines = data, timestamps = as.POSIXct(character(0), tz = "UTC"))
  if (sum(is_data) < 2) return(empty(raw_lines[!is_data], raw_lines[is_data]))
  first_data <- min(which(is_data))
  header_lines <- raw_lines[seq_len(first_data - 1)]
  data_lines <- raw_lines[is_data]
  ts <- suppressWarnings(as.POSIXct(trimws(sub(",.*", "", data_lines)),
                                    format = time.format, tz = "UTC"))
  valid <- !is.na(ts)
  if (sum(valid) < 2) return(empty(header_lines, data_lines))
  streams <- rep("ahead", length(data_lines))
  vts <- as.numeric(ts[valid])
  behind <- vts < (cummax(vts) - tolerance)
  streams[which(valid)[behind]] <- "behind"
  has_reset <- any(behind)
  list(has_reset = has_reset, streams = streams,
       stream_labels = if (has_reset) c("ahead", "behind") else "ahead",
       header_lines = header_lines, data_lines = data_lines, timestamps = ts)
}

#' Write the "ahead" and "behind" streams of a clock-reset file to separate files
#'
#' Companion to [detect_clock_streams()]: writes one file per stream (header
#' lines plus that stream's data lines) into `dir`, so each can be imported
#' independently with [goFlux::import.UGGA()] and tagged.
#'
#' @param streams output of [detect_clock_streams()].
#' @param filename base name used for the output files.
#' @param dir output directory (created if needed).
#' @return named character vector of the files written (names are the stream
#'   labels).
#' @examples
#' f <- tempfile(fileext = ".txt")
#' t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")
#' writeLines(c("SysTime, [CH4]d_ppm",
#'              sprintf("%s, 2.0", format(t0 + 0:9, "%m/%d/%Y %H:%M:%S")),
#'              sprintf("%s, 2.1", format(t0 - 600 + 0:9, "%m/%d/%Y %H:%M:%S"))), f)
#' out <- write_stream_files(detect_clock_streams(f), "micro_f0001.txt", tempdir())
#' out
#' @export
write_stream_files <- function(streams, filename, dir = tempfile("streams_")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out <- character(0)
  for (label in streams$stream_labels) {
    lines <- c(streams$header_lines, streams$data_lines[streams$streams == label])
    fn <- file.path(dir, sub("\\.txt$", paste0("_", label, ".txt"), filename))
    writeLines(lines, fn)
    out[label] <- fn
  }
  out
}

#' Separate two interleaved concentration streams within one observation window
#'
#' Simple state machine for a short (about 10 min) window in which two data
#' streams alternate: a point-to-point jump larger than `threshold` switches
#' the stream label. Within a single window the two streams are well
#' separated, so this works where a global separation would fail at
#' convergence points or chamber events.
#'
#' @param window_data data.frame for one window (must contain `gas_col`).
#' @param gas_col concentration column used to detect jumps; `NULL` picks the
#'   first column matching `CO2.*dry`.
#' @param threshold jump size that switches streams; default 3 times the
#'   median absolute point-to-point difference. The default assumes most
#'   consecutive rows belong to the same stream; when the two streams
#'   alternate on every row the median jump is itself large and you must
#'   supply a `threshold` (e.g. a value between the within-stream and the
#'   between-stream jump sizes).
#' @return `window_data` with a `stream_id` column (`"stream_a"` /
#'   `"stream_b"`) and an attribute `"roughness"` giving the mean absolute
#'   difference of each stream (the smoother stream is usually the real one).
#' @examples
#' t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")
#' a <- 400 + 0.2 * (0:199); b <- 900 + 0.2 * (0:199)
#' d <- data.frame(POSIX.time = t0 + rep(0:199, each = 2),
#'                 CO2dry_ppm = as.vector(rbind(a, b)))
#' s <- separate_window_streams(d)
#' table(s$stream_id); attr(s, "roughness")
#' @export
separate_window_streams <- function(window_data, gas_col = NULL, threshold = NULL) {
  gas_col <- pick_gas_col(window_data, gas_col)
  conc <- window_data[[gas_col]]
  n <- length(conc)
  if (n < 4) { window_data$stream_id <- rep("stream_a", n); return(window_data) }
  d <- abs(diff(conc))
  if (is.null(threshold)) threshold <- stats::median(d, na.rm = TRUE) * 3
  lab <- character(n); lab[1] <- "stream_a"
  for (i in 2:n) {
    lab[i] <- if (!is.na(d[i - 1]) && d[i - 1] > threshold)
      if (lab[i - 1] == "stream_a") "stream_b" else "stream_a" else lab[i - 1]
  }
  window_data$stream_id <- lab
  rough <- vapply(c("stream_a", "stream_b"), function(s) {
    v <- conc[lab == s]; if (length(v) < 2) NA_real_ else mean(abs(diff(v)), na.rm = TRUE)
  }, numeric(1))
  attr(window_data, "roughness") <- rough
  window_data
}

pick_gas_col <- function(d, gas_col) {
  if (!is.null(gas_col) && gas_col %in% names(d)) return(gas_col)
  cand <- names(d)[grep("CO2.*dry|CO2d", names(d), ignore.case = TRUE)]
  if (length(cand) == 0) stop("cannot find a gas column; available: ", paste(names(d), collapse = ", "))
  cand[1]
}

#' Build observation windows from interleaved data, de-interleaving per window
#'
#' An alternative to [goFlux::obs.win()] for records that contain two
#' interleaved clock streams (see [detect_clock_streams()]). For each auxfile
#' row it slices the raw data to `[start - shoulder, start + obs_length +
#' shoulder]`, separates the streams within that window with
#' [separate_window_streams()], keeps one of them, and adds `Etime` and the
#' chamber columns from the auxfile so the result can go straight to
#' [goFlux::click.peak2()], [click_peak2_stacked()] or [windows_from_table()].
#'
#' @param raw_data imported data with `POSIX.time` and gas columns; tagging
#'   columns (`source_file`, `serial_number`, `clock_stream`, `row_order`,
#'   `stream_id`) are stripped.
#' @param auxfile data.frame with `UniqueID`, `start.time` and chamber columns
#'   (`Area`, `offset`, `Vcham`, `Vtube`, `Vinst`, `Vtot`, `Tcham`, `Pcham`;
#'   only those present are copied).
#' @param gas_col column used to separate streams.
#' @param obs_length,shoulder window length and shoulder in seconds.
#' @param threshold passed to [separate_window_streams()].
#' @param keep `"smoother"` (default: the stream with the lower roughness),
#'   `"stream_a"` or `"stream_b"`.
#' @param tz time zone used if `start.time` is character.
#' @return named list of data.frames (one per `UniqueID`), like `obs.win()`.
#' @examples
#' t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")
#' a <- 400 + 0.2 * (0:599); b <- 900 + 5 * sin(0:599)
#' raw <- data.frame(POSIX.time = t0 + rep(0:599, each = 2),
#'                   CO2dry_ppm = as.vector(rbind(a, b)))
#' aux <- data.frame(UniqueID = "w1", start.time = t0 + 120, Area = 324, Vtot = 4.3,
#'                   Tcham = 25, Pcham = 101.3)
#' ow <- obs_win_deinterleave(raw, aux, obs_length = 300, shoulder = 60)
#' nrow(ow$w1); range(ow$w1$CO2dry_ppm)
#' @export
obs_win_deinterleave <- function(raw_data, auxfile, gas_col = "CO2dry_ppm",
                                 obs_length = 300, shoulder = 300, threshold = NULL,
                                 keep = c("smoother", "stream_a", "stream_b"), tz = "UTC") {
  keep <- match.arg(keep)
  gas_col <- pick_gas_col(raw_data, gas_col)
  extra <- intersect(names(raw_data), c("source_file", "serial_number", "clock_stream",
                                        "row_order", "stream_id"))
  cham_cols <- intersect(names(auxfile), c("Area", "offset", "Vcham", "Vtube", "Vinst",
                                           "Vtot", "Tcham", "Pcham"))
  out <- list()
  for (i in seq_len(nrow(auxfile))) {
    row <- auxfile[i, ]
    uid <- as.character(row$UniqueID)
    start_t <- row$start.time
    if (!inherits(start_t, "POSIXct")) start_t <- as.POSIXct(start_t, tz = tz)
    off <- if ("offset" %in% names(row) && !is.na(row$offset)) row$offset else 0
    adj <- start_t + off
    w <- raw_data[raw_data$POSIX.time >= adj - shoulder &
                  raw_data$POSIX.time <= adj + obs_length + shoulder, , drop = FALSE]
    if (nrow(w) < 5) { message(uid, ": only ", nrow(w), " rows in window; skipped"); next }
    w <- separate_window_streams(w, gas_col = gas_col, threshold = threshold)
    if (length(unique(w$stream_id)) > 1) {
      rough <- attr(w, "roughness")
      chosen <- switch(keep,
                       smoother = if (isTRUE(rough["stream_a"] <= rough["stream_b"])) "stream_a" else "stream_b",
                       keep)
      w <- w[w$stream_id == chosen, , drop = FALSE]
    }
    w <- w[, setdiff(names(w), c(extra, "stream_id")), drop = FALSE]
    w$UniqueID <- uid
    w$start.time <- adj
    w$obs.length <- obs_length
    w$Etime <- as.numeric(difftime(w$POSIX.time, adj, units = "secs"))
    for (cc in cham_cols) w[[cc]] <- row[[cc]]
    rownames(w) <- NULL
    out[[uid]] <- w
  }
  out
}

#' Interactive 2-D (time x concentration) window selection with immediate flux calculation
#'
#' Replaces the `obs.win() -> click.peak2() -> goFlux()` sequence with one
#' interactive step: for each auxfile row the full window is plotted, you
#' click two opposite corners of a box around the good data (selecting both
#' the time range *and* the concentration band, which naturally excludes an
#' interleaved second stream or an outlier band), and [goFlux::goFlux()] is
#' run on the selected points at once. Blue lines show the expected start and
#' end. Right-click / Escape skips a measurement.
#'
#' Requires an interactive graphics device that supports
#' [graphics::locator()]; the function opens one with [grDevices::dev.new()]
#' (on macOS this is a Quartz window, which unlike the RStudio pane supports
#' `locator()`).
#'
#' @inheritParams obs_win_deinterleave
#' @param gas_col gas column to plot and to pass to `goFlux()`.
#' @param H2O_col water column passed to `goFlux()`; if absent from the data a
#'   zero column is added (no water correction).
#' @param plot.lim optional y limits.
#' @param save.plots optional directory receiving one PNG per selection.
#' @param ... further arguments to [goFlux::goFlux()] (e.g. `warn.length`).
#' @return data.frame of `goFlux()` results (pass to [goFlux::best.flux()]),
#'   or `NULL` if nothing was processed. The selected points are returned in
#'   the attribute `"selected"` (a click.peak2-style data frame with
#'   `flag = 1`).
#' @examples
#' \dontrun{
#' res <- click_flux_2d(raw, aux, gas_col = "CO2dry_ppm", obs_length = 300,
#'                      shoulder = 300, plot.lim = c(200, 1200))
#' best <- goFlux::best.flux(res)
#' }
#' @export
click_flux_2d <- function(raw_data, auxfile, gas_col = "CO2dry_ppm", H2O_col = "H2O_ppm",
                          obs_length = 300, shoulder = 300, plot.lim = NULL,
                          save.plots = NULL, tz = "UTC", ...) {
  gas_col <- pick_gas_col(raw_data, gas_col)
  extra <- intersect(names(raw_data), c("source_file", "serial_number", "clock_stream",
                                        "row_order", "stream_id"))
  cham_cols <- intersect(names(auxfile), c("Area", "offset", "Vcham", "Vtube", "Vinst",
                                           "Vtot", "Tcham", "Pcham"))
  if (!is.null(save.plots)) dir.create(save.plots, showWarnings = FALSE, recursive = TRUE)
  dev.new(noRStudioGD = TRUE, width = 10, height = 5)
  idev <- dev.cur()
  on.exit(try(dev.off(idev), silent = TRUE))
  message("Click BOTTOM-LEFT then TOP-RIGHT of the good data; right-click/Esc to skip.")
  results <- list(); selected_all <- list()
  draw <- function(w, conc, ylim, uid, i, sel = NULL, box = NULL) {
    plot(w$POSIX.time, conc, main = paste0("[", i, "/", nrow(auxfile), "]  ", uid),
         xlab = "Time", ylab = gas_col, pch = 16, cex = 0.5, ylim = ylim,
         col = if (is.null(sel)) "black" else "gray70")
    if (!is.null(sel)) points(sel$POSIX.time, sel[[gas_col]], pch = 16, cex = 0.5, col = "dodgerblue")
    if (!is.null(box)) rect(box[1], box[3], box[2], box[4], border = "red", lwd = 2)
  }
  for (i in seq_len(nrow(auxfile))) {
    row <- auxfile[i, ]; uid <- as.character(row$UniqueID)
    start_t <- row$start.time
    if (!inherits(start_t, "POSIXct")) start_t <- as.POSIXct(start_t, tz = tz)
    off <- if ("offset" %in% names(row) && !is.na(row$offset)) row$offset else 0
    adj <- start_t + off
    w <- raw_data[raw_data$POSIX.time >= adj - shoulder &
                  raw_data$POSIX.time <= adj + obs_length + shoulder, , drop = FALSE]
    w <- w[, setdiff(names(w), extra), drop = FALSE]
    if (nrow(w) < 5) { message(uid, ": only ", nrow(w), " rows in window; skipped"); next }
    if (!(idev %in% dev.list())) { dev.new(noRStudioGD = TRUE, width = 10, height = 5); idev <- dev.cur() }
    dev.set(idev)
    conc <- w[[gas_col]]
    ylim <- if (!is.null(plot.lim)) plot.lim else range(conc, na.rm = TRUE) * c(0.98, 1.02)
    draw(w, conc, ylim, uid, i)
    abline(v = adj, col = "blue", lwd = 1.5); abline(v = adj + obs_length, col = "blue", lwd = 1.5, lty = 2)
    pts <- tryCatch(locator(2), error = function(e) NULL)
    if (is.null(pts) || length(pts$x) < 2) { message(uid, ": skipped"); next }
    x1 <- as.POSIXct(min(pts$x), origin = "1970-01-01", tz = tz)
    x2 <- as.POSIXct(max(pts$x), origin = "1970-01-01", tz = tz)
    y1 <- min(pts$y); y2 <- max(pts$y)
    sel <- w[w$POSIX.time >= x1 & w$POSIX.time <= x2 & conc >= y1 & conc <= y2, , drop = FALSE]
    box <- c(as.numeric(x1), as.numeric(x2), y1, y2)
    draw(w, conc, ylim, uid, i, sel, box)
    if (!is.null(save.plots)) {
      tryCatch({
        png(file.path(save.plots, paste0(uid, "_selected.png")), width = 1000, height = 500, res = 120)
        draw(w, conc, ylim, uid, i, sel, box); dev.off()
      }, error = function(e) { try(dev.off(), silent = TRUE); message("could not save plot: ", conditionMessage(e)) })
      dev.set(idev)
    }
    if (nrow(sel) < 5) { message(uid, ": too few points selected; skipped"); next }
    sel$UniqueID <- uid; sel$flag <- 1
    sel$Etime <- as.numeric(difftime(sel$POSIX.time, adj, units = "secs"))
    sel$start.time <- adj; sel$obs.length <- obs_length
    sel$start.time_corr <- x1; sel$end.time_corr <- x2
    sel$obs.length_corr <- as.numeric(difftime(x2, x1, units = "secs"))
    for (cc in cham_cols) sel[[cc]] <- row[[cc]]
    if (!H2O_col %in% names(sel)) sel[[H2O_col]] <- 0
    selected_all[[uid]] <- sel
    r <- tryCatch(goFlux::goFlux(sel, gas_col, H2O_col = H2O_col, ...),
                  error = function(e) { message(uid, ": goFlux error: ", conditionMessage(e)); NULL })
    if (!is.null(r)) results[[uid]] <- r
  }
  if (length(results) == 0) { message("No measurements processed."); return(NULL) }
  out <- do.call(rbind, results); rownames(out) <- NULL
  attr(out, "selected") <- do.call(rbind, lapply(selected_all, as.data.frame))
  out
}
