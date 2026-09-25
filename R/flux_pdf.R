#' Save goFlux flux plots to a multi-page PDF (works with ggplot2 >= 3.5)
#'
#' A replacement for [goFlux::flux2pdf()]. goFlux v0.2.0 validated the plot
#' list with `sapply(plot.list, class)[1, ]`, which breaks once ggplot objects
#' carry more than two classes (ggplot2 >= 3.5 / the S7-based 4.0 series) and
#' in the released version the check is commented out. `flux_pdf()` checks
#' each element with `inherits(x, "ggplot")`, groups the plots of one
#' `UniqueID` (e.g. CO2 and CH4) on one page, and needs no `rlist` or
#' `pbapply`.
#'
#' @param plot.list list of ggplot objects from [goFlux::flux.plot()]; plots
#'   from several gases can be concatenated with `c()`.
#' @param outfile path of the PDF to write. Default: `<name of plot.list>.pdf`
#'   in the current working directory.
#' @param width,height page size in inches.
#' @param group.by how to group plots on pages: `"UniqueID"` (default, read
#'   from each plot's `plot_env`, as goFlux does) or `"none"` (one plot per
#'   page).
#' @return the output path, invisibly.
#' @examples
#' \dontrun{
#' CO2_plots <- goFlux::flux.plot(CO2_best, manID, "CO2dry_ppm")
#' CH4_plots <- goFlux::flux.plot(CH4_best, manID, "CH4dry_ppb")
#' flux_pdf(c(CO2_plots, CH4_plots), "flux_results.pdf")
#' }
#' @export
flux_pdf <- function(plot.list, outfile = NULL, width = 11.6, height = 8.2,
                     group.by = c("UniqueID", "none")) {
  group.by <- match.arg(group.by)
  if (missing(plot.list)) stop("'plot.list' is required")
  if (!is.list(plot.list)) stop("'plot.list' must be a list")
  is_gg <- vapply(plot.list, function(x) inherits(x, "ggplot"), logical(1))
  if (!all(is_gg)) stop("all elements of 'plot.list' must be ggplot objects")
  if (is.null(outfile)) outfile <- file.path(getwd(), paste0(deparse(substitute(plot.list)), ".pdf"))
  if (!is.character(outfile)) stop("'outfile' must be a character string")
  if (!is.numeric(width) || !is.numeric(height)) stop("'width' and 'height' must be numeric")

  uids <- vapply(seq_along(plot.list), function(i) plot_uid(plot.list[[i]], i), character(1))
  groups <- if (group.by == "UniqueID") split(plot.list, factor(uids, levels = unique(uids)))
            else stats::setNames(lapply(plot.list, list), uids)
  if (length(groups) > 10000) stop("the outfile would contain more than 10,000 pages")

  pages <- lapply(seq_along(groups), function(k) {
    p <- groups[[k]]
    title <- grid::textGrob(paste("Unique ID:", names(groups)[k]), gp = grid::gpar(fontsize = 16))
    footnote <- grid::textGrob(paste("page", k, "of", length(groups)),
                               gp = grid::gpar(fontface = 3, fontsize = 12), hjust = 1)
    n.plot <- length(p)
    nrow <- if (n.plot <= 2) 1 else if (n.plot <= 4) 2 else 4
    ncol <- if (n.plot <= 1) 1 else 2
    gridExtra::marrangeGrob(grobs = p, ncol = ncol, nrow = nrow, top = title, bottom = footnote)
  })
  pdf(file = outfile, width = width, height = height)
  on.exit(dev.off())
  for (pg in pages) print(pg)
  invisible(outfile)
}

# Extract the UniqueID that goFlux::flux.plot() leaves in the plot environment;
# fall back to the plot index.
plot_uid <- function(p, i) {
  uid <- tryCatch({
    env <- if (isS4(p) || inherits(p, "S7_object")) attr(p, "plot_env") else p$plot_env
    if (is.null(env)) env <- tryCatch(p$plot_env, error = function(e) NULL)
    if (is.environment(env) && exists("UniqueID", envir = env, inherits = FALSE))
      as.character(get("UniqueID", envir = env))[1] else NA_character_
  }, error = function(e) NA_character_)
  if (is.na(uid) || !nzchar(uid)) uid <- paste0("plot_", i)
  uid
}
