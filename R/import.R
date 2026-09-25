#' Import a date subset of a huge LI-7810 .data file
#'
#' LI-COR LI-7810 analyzers write one continuous tab-delimited `.data` file
#' that can grow to hundreds of megabytes over a campaign; reading all of it
#' into R with [goFlux::import.LI7810()] is slow and memory-hungry when only a
#' few field days are needed.
#'
#' This function is a thin shim over the proposed upstream change (see
#' `inst/upstream/patches/import.LI7810-dates.patch`): if the installed
#' goFlux `import.LI7810()` has a `dates` argument, the call is forwarded to
#' it. Otherwise the file is subset here (header lines plus the `DATA` rows
#' whose `DATE` column is in `dates`; `awk` when available, an R stream
#' otherwise), written to `outfile`, and `import.LI7810()` is run on the
#' subset.
#'
#' @param inputfile path to the raw `.data` file.
#' @param dates character or Date vector of dates to keep, formatted as they
#'   appear in the file's `DATE` column (`YYYY-MM-DD` by default).
#' @param outfile path of the subset file (fallback path only); default a
#'   temporary file.
#' @param header.lines number of header lines copied verbatim (the LI-7810
#'   format has 7: five metadata lines, `DATAH` and `DATAU`).
#' @param date.col 1-based index of the `DATE` column in `DATA` rows (7 in
#'   the standard format).
#' @param import if `TRUE` (default) the data frame is returned; if `FALSE`
#'   the path of the subset file (fallback path only).
#' @param use.awk use `awk` if found on the PATH (fallback path only).
#' @param chunk.size lines per chunk for the R fallback.
#' @param ... arguments passed to [goFlux::import.LI7810()] (`date.format`,
#'   `timezone`, `prec`, `keep_all`, ...).
#' @return the imported data frame (with attribute `"subset.file"` on the
#'   fallback path), or the subset file path when `import = FALSE`.
#' @examples
#' f <- system.file("extdata", "LI7810", "LI7810.data", package = "goFlux")
#' if (nzchar(f)) {
#'   imp <- import_li7810_subset(f, dates = "2022-12-05", timezone = "UTC")
#'   nrow(imp)
#' }
#' @export
import_li7810_subset <- function(inputfile, dates, outfile = tempfile(fileext = ".data"),
                                 header.lines = 7, date.col = 7, import = TRUE,
                                 use.awk = TRUE, chunk.size = 200000, ...) {
  if (!file.exists(inputfile)) stop("file not found: ", inputfile)
  dates <- as.character(dates)
  # goFlux's importers inspect match.call()[[1]] and fail when called as
  # goFlux::import.LI7810(); bind locally so the call head is a plain symbol.
  importer <- get("import.LI7810", envir = asNamespace("goFlux"))
  if ("dates" %in% names(formals(importer)) && import) {
    return(importer(inputfile = inputfile, dates = dates, ...))
  }
  awk <- if (use.awk) Sys.which("awk") else ""
  if (nzchar(awk)) {
    prog <- sprintf('NR<=%d {print; next} $1=="DATA" && $%d ~ /^(%s)$/ {print}',
                    header.lines, date.col, paste(dates, collapse = "|"))
    status <- system2(awk, c("-F", shQuote("\t"), shQuote(prog), shQuote(inputfile)),
                      stdout = outfile)
    if (status != 0) stop("awk subsetting failed (status ", status, ")")
  } else {
    con <- file(inputfile, "r")
    out <- file(outfile, "w")
    tryCatch({
      hdr <- readLines(con, n = header.lines); writeLines(hdr, out)
      repeat {
        ch <- readLines(con, n = chunk.size)
        if (length(ch) == 0) break
        f <- strsplit(ch, "\t", fixed = TRUE)
        keep <- vapply(f, function(v) length(v) >= date.col && v[1] == "DATA" && v[date.col] %in% dates, logical(1))
        if (any(keep)) writeLines(ch[keep], out)
      }
    }, finally = { close(con); close(out) })
  }
  n <- length(readLines(outfile)) - header.lines
  if (n <= 0) warning("no DATA rows matched dates: ", paste(dates, collapse = ", "))
  message("Subset written: ", n, " data rows -> ", outfile)
  if (!import) return(outfile)
  res <- importer(inputfile = outfile, ...)
  attr(res, "subset.file") <- outfile
  res
}

#' Import a folder of raw analyzer files without goFlux's import2RData wrapper
#'
#' [goFlux::import2RData()] loops over files with a `pbapply` progress bar and
#' a `withCallingHandlers()` that calls `invokeRestart("muffleError")`, a
#' restart that does not exist, so any error inside a single file import is
#' turned into `no restart 'muffleError' found` and the whole run stops. It
#' also writes its `.RData` output into `getwd()/RData` and, for some
#' instruments, does not search subfolders. This function instead:
#'
#' * finds the files recursively (with an optional file-name `pattern`, e.g.
#'   `"_f\\\\d+\\\\.txt$"` for LGR flux files) and skips empty files;
#' * calls the per-file goFlux importer (`import.UGGA()`, `import.LI7810()`,
#'   ...) directly with `save = FALSE`, so nothing is written to disk;
#' * catches errors per file, reports them, and continues;
#' * returns one merged data frame (with a `source_file` column) or a list.
#'
#' @param path folder containing the raw files.
#' @param instrument goFlux instrument name, i.e. the suffix of an exported
#'   `goFlux::import.<instrument>()` function (`"UGGA"`, `"LI7810"`,
#'   `"G2508"`, ...).
#' @param pattern optional regular expression for file names.
#' @param recursive search subfolders.
#' @param min.size files smaller than this many bytes are skipped.
#' @param merge return one data frame (`TRUE`) or a named list (`FALSE`).
#' @param ... arguments passed to the importer (`date.format`, `timezone`,
#'   `prec`, `keep_all`, ...); only those the importer accepts are forwarded.
#' @return data.frame (with `source_file`) or list of data frames. Files that
#'   failed are listed in the attribute `"failed"`.
#' @examples
#' d <- system.file("extdata", "LI7810", package = "goFlux")
#' if (nzchar(d)) {
#'   imp <- import2RData_flat(d, "LI7810", pattern = "\\.data$", timezone = "UTC")
#'   nrow(imp); unique(imp$source_file)
#' }
#' @export
import2RData_flat <- function(path, instrument, pattern = NULL, recursive = TRUE,
                              min.size = 1, merge = TRUE, ...) {
  fun_name <- paste0("import.", instrument)
  if (!exists(fun_name, envir = asNamespace("goFlux")))
    stop("goFlux has no importer '", fun_name, "'")
  fun <- get(fun_name, envir = asNamespace("goFlux"))
  files <- list.files(path, pattern = pattern, recursive = recursive, full.names = TRUE)
  files <- files[!dir.exists(files)]
  files <- files[!grepl("^\\._", basename(files))]
  files <- files[file.size(files) >= min.size]
  if (length(files) == 0) stop("no files found in ", path)
  dots <- list(...)
  ok <- intersect(names(dots), names(formals(fun)))
  args <- dots[ok]
  if ("save" %in% names(formals(fun))) args$save <- FALSE
  out <- list(); failed <- character(0)
  for (f in files) {
    # call by name (not by function object): goFlux importers inspect
    # match.call()[[1]] and need a plain symbol there
    r <- tryCatch(do.call(fun_name, c(list(inputfile = f), args), envir = asNamespace("goFlux")),
                  error = function(e) { message(basename(f), ": ", conditionMessage(e)); NULL })
    if (is.null(r) || nrow(r) == 0) { failed <- c(failed, f); next }
    r$source_file <- basename(f)
    out[[basename(f)]] <- as.data.frame(r)
  }
  message(length(out), " file(s) imported, ", length(failed), " failed")
  if (!merge) { attr(out, "failed") <- failed; return(out) }
  if (length(out) == 0) stop("no file could be imported")
  res <- do.call(rbind, out); rownames(res) <- NULL
  res <- res[order(res$POSIX.time), ]
  rownames(res) <- NULL
  attr(res, "failed") <- failed
  res
}
