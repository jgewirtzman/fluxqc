#' One real chamber closure, ready for goFlux
#'
#' Loads the example closure shipped in `inst/extdata/example_closure.csv` and
#' returns it as a goFlux-ready data frame (the same structure as the output of
#' [goFlux::click.peak2()] or [windows_from_table()]).
#'
#' The closure is `2021-05-18_5086_I_7.0`: a 600 s semi-rigid sleeve chamber
#' measurement on a tree stem at Yale-Myers Forest (Connecticut, USA) recorded
#' at 1 Hz with an ABB/LGR micro-portable UGGA (GLA131). It is a clean,
#' high-signal closure in both CO2 and CH4. Chamber geometry comes from the
#' campaign auxfile: `Area` 4254.121 cm2, `Vtot` 7.082698 L, `Tcham` 26.02 C,
#' `Pcham` 101.325 kPa. The raw record holds dry mole fractions only, so
#' `H2O_ppm` is set to 0 (no water-vapour correction of the flux term).
#'
#' @param dead.band seconds discarded at the start of the closure (the
#'   mixing "dead band"). Rows with `t_s >= dead.band` get `flag = 1`.
#' @param prec manufacturer precision `c(CO2 ppm, CH4 ppb, H2O ppm)` written to
#'   the `*_prec` columns, as in the source repository's `import.UGGA()` call.
#' @return data.frame with columns `UniqueID`, `POSIX.time`, `t_s`,
#'   `CO2dry_ppm`, `CH4dry_ppb`, `H2O_ppm`, `start.time`, `obs.length`, `Area`,
#'   `Vtot`, `Tcham`, `Pcham`, `CO2_prec`, `CH4_prec`, `H2O_prec`, `flag`,
#'   `Etime`, `start.time_corr`, `end.time_corr`, `obs.length_corr`.
#' @examples
#' d <- example_closure()
#' nrow(d); sum(d$flag)
#' plot(d$t_s, d$CH4dry_ppb, col = d$flag + 1, pch = 16, cex = 0.4)
#' @export
example_closure <- function(dead.band = 60, prec = c(0.35, 0.9, 200)) {
  f <- system.file("extdata", "example_closure.csv", package = "stemflux")
  if (f == "") stop("example_closure.csv not found in the installed package")
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$POSIX.time <- as.POSIXct(d$POSIX.time, tz = "UTC",
                             format = "%Y-%m-%d %H:%M:%OS")
  start <- as.POSIXct("2021-05-18 13:36:00", tz = "UTC")
  d$H2O_ppm    <- 0
  d$start.time <- start
  d$obs.length <- 600
  d$Area  <- 4254.121
  d$Vtot  <- 7.082698
  d$Tcham <- 26.02
  d$Pcham <- 101.325
  d$CO2_prec <- prec[1]
  d$CH4_prec <- prec[2]
  d$H2O_prec <- prec[3]
  flag_window(d, start + dead.band, start + max(d$t_s))
}
