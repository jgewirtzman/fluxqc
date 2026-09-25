#' Run the core chamber-flux pipeline for one gas in one call
#'
#' The path is `goFlux::goFlux()` -> `goFlux::best.flux()` ->
#' [flag_detection()] -> [qc_screens()], with every option recorded in a
#' `settings` element for the methods section. goFlux's columns are returned
#' exactly as `best.flux()` produced them; fluxqc only appends columns
#' (`sigma_emp`, `MDF_emp`, `below_MDF_emp`, `MDF_emp_method`,
#' `detected_emp`, `det_class_emp`, `qc_*`). Ebullition is not handled here:
#' closures over water or saturated surfaces go to [goFlux::goAquaFlux()]
#' directly.
#'
#' @param manID_data flagged concentration data for all closures (output of
#'   [goFlux::click.peak2()], [click_peak2_stacked()], [auto_id_rise()] or
#'   [windows_from_table()]).
#' @param aux optional auxfile (one row per `UniqueID`); used for the
#'   `group` column of the QC screens, which is merged into the flux table
#'   by `UniqueID`.
#' @param gastype gas column to process (one gas per call).
#' @param precision,mdf,conf,extra passed to [flag_detection()]; the defaults are the
#'   companion manuscript's method (record-wide MAD precision per `group`, Wassmann
#'   precision and the Wassmann form at 95 percent two-sided, z = qnorm(0.975) = 1.96,
#'   i.e. the 2-sigma criterion).
#' @param qc list of screens passed to [qc_screens()] (`screens` argument), or
#'   `FALSE` to skip QC. Each screen is a toggle; set `co2_tracer = FALSE` for
#'   surfaces without respiratory CO2 efflux (open water, lit foliage, dead
#'   wood).
#' @param group column name (in `aux` or the flux table) defining the field
#'   periods for the record-wide `"mad"` precision in [flag_detection()] and the
#'   `c0` median in [qc_screens()].
#' @param co2 optional CO2 flux table for the `co2_tracer` screen.
#' @param H2O_col water column passed to `goFlux()`; a zero column is added if
#'   absent.
#' @param best.flux_args list of arguments to [goFlux::best.flux()].
#' @param ... further arguments to [goFlux::goFlux()] (e.g. `warn.length`,
#'   `k.min`).
#' @return list of class `fluxqc_result` with `fluxes` (one row per
#'   `UniqueID`), `traces` (the concentration data fitted), `gastype` and
#'   `settings` (every option used).
#' @seealso [write_outputs()]
#' @examples
#' d <- example_closure()
#' res <- process_fluxes(d, gastype = "CH4dry_ppb", qc = list(c0 = TRUE, convex = TRUE))
#' res
#' res$fluxes[, c("UniqueID", "best.flux", "MDF", "MDF_emp", "MDF_emp_method",
#'                "det_class_emp", "qc_any")]
#' @export
process_fluxes <- function(manID_data, aux = NULL, gastype = "CH4dry_ppb",
                           precision = "mad", mdf = "wassmann", conf = 0.95, extra = NULL,
                           qc = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE,
                                     min_window = TRUE, ambient_start = TRUE, noisy = TRUE),
                           group = NULL, co2 = NULL,
                           H2O_col = "H2O_ppm", best.flux_args = list(), ...) {
  d <- as.data.frame(manID_data)
  if (!gastype %in% names(d)) stop("column '", gastype, "' not in 'manID_data'")
  if (!"flag" %in% names(d)) stop("'manID_data' needs a 'flag' column (run window selection first)")
  if (!H2O_col %in% names(d)) d[[H2O_col]] <- 0
  ids <- unique(as.character(d$UniqueID))

  gf <- goFlux::goFlux(d, gastype, H2O_col = H2O_col, ...)
  best <- as.data.frame(do.call(goFlux::best.flux, c(list(flux.result = gf), best.flux_args)))
  if (!is.null(aux)) {
    extra_cols <- setdiff(intersect(group, names(aux)), names(best))
    i <- match(best$UniqueID, aux$UniqueID)
    for (cc in extra_cols) best[[cc]] <- aux[[cc]][i]
  }
  best <- flag_detection(best, traces = d, gastype = gastype, precision = precision,
                         mdf = mdf, conf = conf, group = group, extra = extra)
  if (!isFALSE(qc)) {
    best <- qc_screens(best, traces = d, screens = qc, group = group, gastype = gastype,
                       co2 = co2, verbose = FALSE)
  }

  settings <- list(
    fluxqc_version = as.character(utils::packageVersion("fluxqc")),
    goFlux_version = as.character(utils::packageVersion("goFlux")),
    gastype = gastype, H2O_col = H2O_col,
    precision = precision, mdf = mdf, conf = conf,
    extra = if (is.null(extra)) NULL else if (identical(extra, "all")) "all" else extra,
    sigma_global = attr(best, "sigma_global"),
    qc = if (isFALSE(qc)) NULL else attr(best, "qc_screens"),
    group = group,
    goFlux_args = list(...), best.flux_args = best.flux_args,
    n_closures = length(ids), date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  structure(list(fluxes = best, traces = d, gastype = gastype, settings = settings),
            class = "fluxqc_result")
}

#' @export
print.fluxqc_result <- function(x, ...) {
  cat("fluxqc result:", nrow(x$fluxes), "closures,", x$gastype, "\n")
  cat("  MDF_emp:", x$fluxes$MDF_emp_method[1], "| detected:", sum(x$fluxes$detected_emp, na.rm = TRUE),
      "| qc_any:", if ("qc_any" %in% names(x$fluxes)) sum(x$fluxes$qc_any, na.rm = TRUE) else "not run", "\n")
  invisible(x)
}

#' Write the outputs of process_fluxes() to a folder
#'
#' Writes `fluxes.csv` (one row per closure), `settings.json` (every option
#' used, for the methods section) and, if `ggplot2` is installed,
#' `flux_plots.pdf` with one page per closure from [goFlux::flux.plot()] via
#' [goFlux::flux2pdf()].
#'
#' @param result output of [process_fluxes()].
#' @param dir output folder (created if needed).
#' @param plots write the per-closure plot PDF.
#' @param ... further arguments to [goFlux::flux.plot()].
#' @return the paths written, invisibly (named character vector).
#' @examples
#' d <- example_closure()
#' res <- process_fluxes(d, gastype = "CH4dry_ppb", qc = list(c0 = TRUE))
#' out <- write_outputs(res, tempfile("fluxqc_"), plots = FALSE)
#' basename(out)
#' @export
write_outputs <- function(result, dir, plots = TRUE, ...) {
  if (!inherits(result, "fluxqc_result")) stop("'result' must come from process_fluxes()")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fmt <- function(df) { df[] <- lapply(df, function(v) if (inherits(v, "POSIXct")) format(v, "%Y-%m-%d %H:%M:%S") else v); df }
  paths <- character(0)
  f <- file.path(dir, "fluxes.csv"); utils::write.csv(fmt(result$fluxes), f, row.names = FALSE); paths["fluxes"] <- f
  f <- file.path(dir, "settings.json"); writeLines(to_json(result$settings), f); paths["settings"] <- f
  if (plots && requireNamespace("ggplot2", quietly = TRUE)) {
    pl <- tryCatch(goFlux::flux.plot(result$fluxes, result$traces, result$gastype, ...),
                   error = function(e) { message("flux.plot failed: ", conditionMessage(e)); NULL })
    if (!is.null(pl)) {
      f <- file.path(dir, "flux_plots.pdf")
      invisible(utils::capture.output(goFlux::flux2pdf(pl, outfile = f)))
      paths["plots"] <- f
    }
  }
  invisible(paths)
}

# Minimal JSON writer (no jsonlite dependency) for the settings list
to_json <- function(x, indent = 0) {
  pad <- strrep("  ", indent)
  esc <- function(s) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", s)), "\"")
  atom <- function(v) {
    if (length(v) == 0) return("[]")
    f <- function(e) if (is.na(e)) "null" else if (is.logical(e)) tolower(as.character(e))
                     else if (is.numeric(e)) format(e, digits = 10) else esc(as.character(e))
    if (length(v) == 1 && is.null(names(v))) f(v)
    else if (!is.null(names(v))) paste0("{", paste(paste0(esc(names(v)), ": ", vapply(v, f, character(1))), collapse = ", "), "}")
    else paste0("[", paste(vapply(v, f, character(1)), collapse = ", "), "]")
  }
  if (is.null(x)) return("null")
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    if (length(x) == 0) return("{}")
    if (is.null(names(x))) return(paste0("[", paste(vapply(x, to_json, character(1), indent = indent + 1), collapse = ", "), "]"))
    items <- vapply(names(x), function(nm) paste0(pad, "  ", esc(nm), ": ", to_json(x[[nm]], indent + 1)), character(1))
    return(paste0("{\n", paste(items, collapse = ",\n"), "\n", pad, "}"))
  }
  atom(x)
}
