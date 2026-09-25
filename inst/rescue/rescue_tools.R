# =============================================================================
# rescue_tools.R -- one-off rescue tools for LGR clock resets and interleaved
# streams (BlueFlux ground campaign, 2022-2023). NOT part of the stemflux
# namespace; see README.md in this folder. Source with:
#   source(system.file("rescue", "rescue_tools.R", package = "stemflux"))
# Requires stemflux (for detect_clock_streams()) and goFlux.
# =============================================================================

# Write the "ahead" and "behind" streams of a clock-reset file to separate files
# so each can be imported independently with goFlux::import.UGGA().
#   streams: output of stemflux::detect_clock_streams()
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

pick_gas_col <- function(d, gas_col) {
  if (!is.null(gas_col) && gas_col %in% names(d)) return(gas_col)
  cand <- names(d)[grep("CO2.*dry|CO2d", names(d), ignore.case = TRUE)]
  if (length(cand) == 0) stop("cannot find a gas column; available: ", paste(names(d), collapse = ", "))
  cand[1]
}

# Label two interleaved streams within ONE short window: a point-to-point jump
# larger than `threshold` switches the stream label. The default threshold
# (3 x median jump) assumes most consecutive rows belong to the same stream;
# when the streams alternate on every row, supply `threshold` explicitly.
# Returns window_data with a stream_id column and attr "roughness".
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

# obs.win() replacement for interleaved records: slice each auxfile row to
# [start - shoulder, start + obs_length + shoulder], separate the streams within
# that window, keep one (the smoother by default), add Etime and chamber columns.
# Returns a named list of data frames (one per UniqueID), like obs.win().
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

# Interactive: for each auxfile row plot the window, click two opposite corners
# of a box around the good points (time range AND concentration band), run
# goFlux::goFlux() on the selection at once. Needs a device supporting locator().
# Returns a data frame of goFlux results (pass to goFlux::best.flux()), with the
# selected points in attr "selected".
click_flux_2d <- function(raw_data, auxfile, gas_col = "CO2dry_ppm", H2O_col = "H2O_ppm",
                          obs_length = 300, shoulder = 300, plot.lim = NULL,
                          save.plots = NULL, tz = "UTC", ...) {
  gas_col <- pick_gas_col(raw_data, gas_col)
  extra <- intersect(names(raw_data), c("source_file", "serial_number", "clock_stream",
                                        "row_order", "stream_id"))
  cham_cols <- intersect(names(auxfile), c("Area", "offset", "Vcham", "Vtube", "Vinst",
                                           "Vtot", "Tcham", "Pcham"))
  if (!is.null(save.plots)) dir.create(save.plots, showWarnings = FALSE, recursive = TRUE)
  grDevices::dev.new(noRStudioGD = TRUE, width = 10, height = 5)
  idev <- grDevices::dev.cur()
  on.exit(try(grDevices::dev.off(idev), silent = TRUE))
  message("Click BOTTOM-LEFT then TOP-RIGHT of the good data; right-click/Esc to skip.")
  results <- list(); selected_all <- list()
  draw <- function(w, conc, ylim, uid, i, sel = NULL, box = NULL) {
    plot(w$POSIX.time, conc, main = paste0("[", i, "/", nrow(auxfile), "]  ", uid),
         xlab = "Time", ylab = gas_col, pch = 16, cex = 0.5, ylim = ylim,
         col = if (is.null(sel)) "black" else "gray70")
    if (!is.null(sel)) graphics::points(sel$POSIX.time, sel[[gas_col]], pch = 16, cex = 0.5, col = "dodgerblue")
    if (!is.null(box)) graphics::rect(box[1], box[3], box[2], box[4], border = "red", lwd = 2)
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
    if (!(idev %in% grDevices::dev.list())) { grDevices::dev.new(noRStudioGD = TRUE, width = 10, height = 5); idev <- grDevices::dev.cur() }
    grDevices::dev.set(idev)
    conc <- w[[gas_col]]
    ylim <- if (!is.null(plot.lim)) plot.lim else range(conc, na.rm = TRUE) * c(0.98, 1.02)
    draw(w, conc, ylim, uid, i)
    graphics::abline(v = adj, col = "blue", lwd = 1.5); graphics::abline(v = adj + obs_length, col = "blue", lwd = 1.5, lty = 2)
    pts <- tryCatch(graphics::locator(2), error = function(e) NULL)
    if (is.null(pts) || length(pts$x) < 2) { message(uid, ": skipped"); next }
    x1 <- as.POSIXct(min(pts$x), origin = "1970-01-01", tz = tz)
    x2 <- as.POSIXct(max(pts$x), origin = "1970-01-01", tz = tz)
    y1 <- min(pts$y); y2 <- max(pts$y)
    sel <- w[w$POSIX.time >= x1 & w$POSIX.time <= x2 & conc >= y1 & conc <= y2, , drop = FALSE]
    box <- c(as.numeric(x1), as.numeric(x2), y1, y2)
    draw(w, conc, ylim, uid, i, sel, box)
    if (!is.null(save.plots)) {
      tryCatch({
        grDevices::png(file.path(save.plots, paste0(uid, "_selected.png")), width = 1000, height = 500, res = 120)
        draw(w, conc, ylim, uid, i, sel, box); grDevices::dev.off()
      }, error = function(e) { try(grDevices::dev.off(), silent = TRUE); message("could not save plot: ", conditionMessage(e)) })
      grDevices::dev.set(idev)
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
