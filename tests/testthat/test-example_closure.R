test_that("example closure loads with goFlux columns", {
  d <- example_closure()
  expect_equal(nrow(d), 603)
  need <- c("UniqueID", "POSIX.time", "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm", "start.time",
            "obs.length", "Area", "Vtot", "Tcham", "Pcham", "flag", "Etime",
            "start.time_corr", "end.time_corr", "obs.length_corr", "CH4_prec", "CO2_prec")
  expect_true(all(need %in% names(d)))
  expect_s3_class(d$POSIX.time, "POSIXct")
  expect_equal(sum(d$flag), sum(d$t_s >= 60))
  # sampling is ~0.994 s, so the first flagged sample sits just past the dead band
  expect_gte(d$Etime[d$flag == 1][1], 0)
  expect_lt(d$Etime[d$flag == 1][1], 1)
})

test_that("example closure runs through goFlux and best.flux", {
  d <- example_closure()
  res <- suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))
  best <- suppressWarnings(goFlux::best.flux(res))
  expect_equal(nrow(best), 1)
  expect_gt(best$LM.r2, 0.99)
  expect_gt(best$best.flux, 0.2)
  expect_equal(best$flux.term, flux_term(7.082698, 101.325, 4254.121, 26.02))
})
