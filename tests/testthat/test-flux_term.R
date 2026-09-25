test_that("flux_term matches goFlux's internal formula", {
  ft <- flux_term(V_L = 7.082698, P_kPa = 101.325, A_cm2 = 4254.121, T_C = 26.02)
  expect_equal(ft, goFlux:::flux.term(7.082698, 101.325, 4254.121, 26.02, 0))
  # V/A = 0.1 m at STP-ish: (0.1 m) * P/RT ~ 4.15 mol m-2
  expect_equal(flux_term(1, 101.3, 100, 20), 0.1 * 101300 / (8.314 * 293.15), tolerance = 1e-6)
})

test_that("water correction reduces the flux term", {
  expect_lt(flux_term(1, 101.3, 100, 20, H2O_mol = 0.02), flux_term(1, 101.3, 100, 20))
})
