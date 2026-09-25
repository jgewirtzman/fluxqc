#' Chamber-contamination screen on the starting concentration (C0)
#'
#' Flags measurements whose chamber starting concentration exceeds
#' `k` times the median starting concentration of their campaign (or group).
#' A chamber that starts far above ambient was not sampling ambient air at
#' closure: it was contaminated or carried over from a previous
#' high-concentration sample. This is a *physical* diagnostic of what went
#' wrong, not a statistical threshold on the flux value, and so is a legitimate
#' hard removal. Statistical filtering of plausible flux values is not done
#' here because it biases campaign means.
#'
#' The default `k = 1.5` was chosen where the source campaigns showed a wide
#' gap in the C0 distribution between contaminated chambers (5-25 times the
#' median) and everything else (below 1.3 times); any multiplier across that
#' gap rejected the same measurements. Check the sensitivity table printed
#' with `verbose = TRUE` for your own data.
#'
#' Rows are retained and flagged, never deleted: filter on `qc_c0_pass`
#' downstream. Rows without a `C0` value get `qc_c0_pass = NA` (unscreened,
#' not verified clean).
#'
#' @param flux.result data.frame from [goFlux::goFlux()] or
#'   [goFlux::best.flux()] with a `C0` column (goFlux's fitted starting
#'   concentration).
#' @param k multiplier of the group median above which a measurement fails.
#' @param C0.col name of the starting-concentration column.
#' @param group optional name of a column defining campaigns/groups; the
#'   median is computed within each group. `NULL` uses all rows.
#' @param verbose print a summary and a threshold-sensitivity line.
#' @return `flux.result` with columns `qc_c0_ratio` (C0 / group median),
#'   `qc_c0_pass` (logical, `NA` if C0 missing) and `qc_c0_reason`
#'   (character, `NA` when passing). The thresholds used are stored in the
#'   attribute `"qc_c0_threshold"`.
#' @examples
#' res <- data.frame(UniqueID = paste0("m", 1:6), C0 = c(2100, 2050, 1980, 2300, 11700, 2200),
#'                   best.flux = c(0.1, -0.2, 0.05, 0.3, 31.6, 0.2))
#' out <- qc_c0_screen(res)
#' out[, c("UniqueID", "C0", "qc_c0_ratio", "qc_c0_pass")]
#' @export
qc_c0_screen <- function(flux.result, k = 1.5, C0.col = "C0", group = NULL,
                         verbose = TRUE) {
  if (!C0.col %in% names(flux.result)) stop("column '", C0.col, "' not found")
  c0 <- flux.result[[C0.col]]
  g <- if (is.null(group)) rep("all", nrow(flux.result)) else as.character(flux.result[[group]])
  med <- tapply(c0, g, stats::median, na.rm = TRUE)
  ratio <- as.numeric(c0) / as.numeric(med[g])
  pass <- ifelse(is.finite(ratio), ratio <= k, NA)
  reason <- ifelse(!is.na(pass) & !pass,
                   sprintf("C0 %.0f = %.1fx group median; chamber not at ambient at closure",
                           c0, ratio), NA_character_)
  flux.result$qc_c0_ratio <- as.numeric(ratio)
  flux.result$qc_c0_pass <- pass
  flux.result$qc_c0_reason <- reason
  thr <- data.frame(group = names(med), median_C0 = as.numeric(med),
                    threshold = as.numeric(med) * k, stringsAsFactors = FALSE)
  attr(flux.result, "qc_c0_threshold") <- thr
  if (verbose) {
    for (i in seq_len(nrow(thr))) {
      sel <- g == thr$group[i]
      message(sprintf("C0 screen [%s]: n = %d, median C0 = %.1f, threshold = %.1f (%.1fx), rejected = %d, unscreened = %d",
                      thr$group[i], sum(sel), thr$median_C0[i], thr$threshold[i], k,
                      sum(!pass[sel], na.rm = TRUE), sum(is.na(pass[sel]))))
      sens <- vapply(c(1.2, 1.5, 2, 5, 10), function(m) sum(ratio[sel] > m, na.rm = TRUE), numeric(1))
      message("  threshold sensitivity: ",
              paste(sprintf("%.1fx->%d", c(1.2, 1.5, 2, 5, 10), sens), collapse = "  "))
    }
  }
  flux.result
}
