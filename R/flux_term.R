#' Chamber flux term (V/A)(P/RT), as used by goFlux
#'
#' Converts a concentration slope (mole fraction per second) into a flux per
#' unit area. This is the same formula as goFlux's internal (unexported)
#' `flux.term()`:
#'
#' \deqn{\mathrm{flux.term} = \frac{V \, P \, (1 - x_{H_2O})}{R \, A \, T}}
#'
#' with `V` in litres, `P` in kPa, `A` in cm2, `T` in degrees Celsius and
#' `R = 8.314` J mol-1 K-1. The litre/kPa/cm2 unit choices cancel so the result
#' is in mol m-2 per unit mole fraction. Multiplying a slope in ppm s-1
#' (umol mol-1 s-1) therefore gives umol m-2 s-1, and a slope in ppb s-1
#' (nmol mol-1 s-1) gives nmol m-2 s-1.
#'
#' @param V_L total system volume in litres (goFlux `Vtot`).
#' @param P_kPa chamber pressure in kPa (goFlux `Pcham`).
#' @param A_cm2 chamber footprint area in cm2 (goFlux `Area`).
#' @param T_C chamber air temperature in degrees Celsius (goFlux `Tcham`).
#' @param H2O_mol water vapour mole fraction (mol mol-1) used to correct to a
#'   dry-air basis. Default 0 (no correction).
#' @return numeric; mol m-2 per unit mole fraction.
#' @examples
#' # 1 L chamber on a 100 cm2 footprint (V/A = 0.1 m), 20 C, 101.3 kPa
#' ft <- flux_term(1, 101.3, 100, 20)
#' ft
#' # a CH4 slope of 0.05 ppb/s becomes
#' 0.05 * ft   # nmol m-2 s-1
#' @export
flux_term <- function(V_L, P_kPa, A_cm2, T_C, H2O_mol = 0) {
  (V_L * P_kPa * (1 - H2O_mol)) / (8.314 * (A_cm2 / 10000) * (T_C + 273.15))
}
