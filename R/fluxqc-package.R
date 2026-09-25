#' fluxqc: detection limits, quality screens and ebullition on top of goFlux
#'
#' `fluxqc` works on any closure processed with the \pkg{goFlux} package (see
#' [goFlux::goFlux()]), whatever the surface enclosed (soil, water, tree
#' stems, branches, leaves, ...), diffusive or ebullitive, and slots into a
#' standard goFlux pipeline
#' (`import -> obs.win -> click.peak2 -> goFlux -> best.flux -> flux.plot`).
#'
#' All functions keep goFlux's column conventions (`UniqueID`, `POSIX.time`,
#' `start.time`, `obs.length`, `Vtot`, `Area`, `Pcham`, `Tcham`, `flag`,
#' `Etime`, `start.time_corr`, `end.time_corr`, `obs.length_corr`,
#' `flux.term`, `nb.obs`, `MDF`, `prec`, `best.flux`, ...) and never modify a
#' goFlux column: fluxqc only appends.
#'
#' @section Core path:
#' `goFlux::goFlux()` -> `goFlux::best.flux()` -> [flag_detection()] ->
#' [qc_screens()], wrapped by [process_fluxes()] and written out by
#' [write_outputs()]. Nothing else runs unless asked for.
#'
#' @section Function map:
#' \describe{
#'   \item{Window selection}{[click_peak2_stacked()], [auto_id_rise()],
#'     [windows_from_table()], [flag_window()], [find_rise()]}
#'   \item{Empirical precision}{[precision_allan()], [precision_mad()],
#'     [precision_rolling()], [precision_by_id()]}
#'   \item{Detection limits}{[flag_detection()], [mdf_multi()], [flux_term()]}
#'   \item{Quality screens}{[qc_screens()], [qc_c0_screen()]}
#'   \item{Pipeline and outputs}{[process_fluxes()], [write_outputs()], [flux_pdf()]}
#'   \item{Import helpers}{[import_li7810_subset()], [import2RData_flat()]}
#'   \item{Utility (not on the default path)}{[find_clock_offset()]}
#'   \item{Ebullition add-on (opt-in)}{[detect_ebullition()],
#'     [split_ebullition()], [merge_ebullition()]; `vignette("ebullition")`}
#'   \item{Example data}{[fluxqc_examples()], [example_closure()]}
#'   \item{Unexported one-off rescue scripts}{`system.file("rescue", package = "fluxqc")`
#'     (LGR clock-reset detection and interleaved-stream tools)}
#' }
#'
#' @references
#' Rheault, K., Christiansen, J. R., & Larsen, K. S. (2024). goFlux: A user-friendly
#' way to calculate GHG fluxes yourself, regardless of user experience. Journal
#' of Open Source Software, 9(96), 6393. \doi{10.21105/joss.06393}
#'
#' Christiansen, J. R., Outhwaite, J., & Smukler, S. M. (2015). Comparison of
#' CO2, CH4 and N2O soil-atmosphere exchange measured in static chambers with
#' cavity ring-down spectroscopy and gas chromatography. Agricultural and
#' Forest Meteorology, 211-212, 48-57. \doi{10.1016/j.agrformet.2015.06.004}
#'
#' Wassmann, R., et al. (2018). Increasing sensitivity of methane emission
#' measurements in rice through deployment of 'closed chambers' at nighttime.
#' PLoS ONE, 13(2), e0191352. \doi{10.1371/journal.pone.0191352}
#'
#' Hoffmann, M., et al. (2017). A simple calculation algorithm to separate
#' high-resolution CH4 flux measurements into ebullition- and diffusion-derived
#' components. Atmospheric Measurement Techniques, 10, 109-118.
#' \doi{10.5194/amt-10-109-2017}
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom goFlux goFlux best.flux
#' @importFrom stats sd mad median lm coef quantile qt qnorm na.omit
#' @importFrom graphics plot abline axis axis.POSIXct par title identify
#'   split.screen screen close.screen legend points rect locator matplot
#' @importFrom grDevices dev.new dev.off dev.list dev.cur dev.set pdf png
#' @importFrom utils read.csv write.csv head tail
NULL
