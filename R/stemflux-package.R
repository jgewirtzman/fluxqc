#' stemflux: extensions to goFlux for stem, soil and water chamber fluxes
#'
#' `stemflux` collects tools written while processing chamber measurements of
#' CO2 and CH4 with the \pkg{goFlux} package (see [goFlux::goFlux()]) and
#' packages them so they slot into a standard goFlux pipeline
#' (`import -> obs.win -> click.peak2 -> goFlux -> best.flux -> flux.plot`).
#'
#' All functions keep goFlux's column conventions (`UniqueID`, `POSIX.time`,
#' `start.time`, `obs.length`, `Vtot`, `Area`, `Pcham`, `Tcham`, `flag`,
#' `Etime`, `start.time_corr`, `end.time_corr`, `obs.length_corr`,
#' `flux.term`, `nb.obs`, `MDF`, `prec`, `best.flux`, ...).
#'
#' @section Function map:
#' \describe{
#'   \item{Window selection}{[click_peak2_stacked()], [auto_id_rise()],
#'     [windows_from_table()], [flag_window()], [find_rise()]}
#'   \item{Empirical precision}{[precision_allan()], [precision_mad()],
#'     [precision_rolling()], [precision_by_id()]}
#'   \item{Detection limits}{[mdf_multi()], [flag_detection()], [flux_term()]}
#'   \item{Quality control}{[qc_c0_screen()], [flux_pdf()]}
#'   \item{Analyzer clock problems}{[detect_clock_streams()],
#'     [separate_window_streams()], [obs_win_deinterleave()], [click_flux_2d()]}
#'   \item{Ebullition}{[detect_ebullition()], [split_ebullition()]}
#'   \item{Import helpers}{[import_li7810_subset()], [import2RData_flat()]}
#'   \item{Example data}{[example_closure()]}
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
#' Forest Meteorology, 211-212, 48-57.
#'
#' Wassmann, R., et al. (2018). Introducing a new tool for greenhouse gas
#' calculation tailored for cropland: rationale, operational framework and
#' potential application. Carbon Management, 9(1), 79-88.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom goFlux goFlux best.flux
#' @importFrom stats sd mad median lm coef quantile qt qnorm na.omit
#' @importFrom graphics plot abline axis axis.POSIXct par title identify
#'   split.screen screen close.screen legend points rect locator
#' @importFrom grDevices dev.new dev.off dev.list dev.cur dev.set pdf png
#' @importFrom utils read.csv write.csv head tail
NULL
