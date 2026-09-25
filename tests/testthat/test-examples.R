test_that("bundled examples load with real geometry and goFlux columns", {
  ex <- stemflux_examples()
  expect_named(ex, c("data", "aux"))
  expect_equal(ex$aux$example, c("a", "b", "c", "d", "e", "f"))
  expect_equal(ex$aux$kind, c(rep("real", 5), "synthetic"))
  expect_equal(sort(unique(ex$data$UniqueID)), sort(ex$aux$UniqueID))
  need <- c("UniqueID", "POSIX.time", "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm", "start.time", "obs.length",
            "Area", "Vtot", "Tcham", "Pcham", "flag", "Etime", "start.time_corr", "end.time_corr",
            "obs.length_corr", "CH4_prec", "CO2_prec")
  expect_true(all(need %in% names(ex$data)))
  expect_equal(ex$aux$Area[1], 4254.121, tolerance = 1e-6); expect_equal(ex$aux$Vtot[1], 7.082698, tolerance = 1e-6)
  expect_equal(ex$aux$Area[5], 324.3); expect_equal(ex$aux$Vtot[5], 4.318)
  expect_equal(as.vector(table(ex$data$UniqueID)[ex$aux$UniqueID]), c(1207L, 1208L, 1073L, 1207L, 278L, 603L))
  # shipped flags: the clicked windows for a-d, whole trace for e, 60 s dead band for f
  fl <- tapply(ex$data$flag, ex$data$UniqueID, sum)[ex$aux$UniqueID]
  expect_equal(as.vector(fl), c(770, 673, 351, 587, 278, 542))
  a <- ex$data[ex$data$UniqueID == "2021-05-18_5086_I_7.0", ]
  expect_equal(as.numeric(difftime(a$start.time_corr[1], a$start.time[1], units = "secs")), 91)
  expect_equal(unique(a$obs.length_corr), 765)
  expect_equal(min(a$Etime[a$flag == 1]), 0)
  # a numeric dead.band re-flags from start.time instead
  re <- stemflux_examples("a", dead.band = 60)$data
  expect_equal(sum(re$flag), sum(re$t_s >= 60))
  one <- stemflux_examples("c")
  expect_equal(unique(one$data$UniqueID), "2021-02-25_1885_I_5_1")
  expect_equal(nrow(stemflux_examples("2020-12-09_5091_I_5")$aux), 1)
  expect_error(stemflux_examples("z"), "no example")
  dg <- read.csv(system.file("extdata", "stemflux_examples_diagnostics.csv", package = "stemflux"))
  expect_equal(dg$UniqueID, ex$aux$UniqueID[1:4])
  expect_equal(dg$class_FINAL, c("emission", "below detection", "below detection", "uptake"))
})

test_that("example_closure() is the alias for (a)", {
  d <- example_closure()
  expect_equal(unique(d$UniqueID), "2021-05-18_5086_I_7.0")
  expect_equal(nrow(d), 1207)
  expect_equal(sum(d$flag), 770)
  expect_equal(d$Tcham[1], 26.02)
})

test_that("the examples separate as intended through process_fluxes (Figure 4 logic)", {
  ex <- stemflux_examples(c("a", "b", "c", "d", "f"))
  res <- suppressWarnings(process_fluxes(ex$data, gastype = "CH4dry_ppb",
                                         qc = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE, min_window = list(secs = 120))))
  fx <- res$fluxes; rownames(fx) <- fx$UniqueID
  a <- fx["2021-05-18_5086_I_7.0", ]; b <- fx["2020-09-23_5092_I_6", ]; cc <- fx["2021-02-25_1885_I_5_1", ]
  d <- fx["2020-12-09_5091_I_5", ]; f <- fx["SYN_uptake_from_a", ]
  expect_equal(a$det_class_emp, "emission"); expect_false(a$qc_co2_tracer)
  expect_equal(a$best.flux, 0.254, tolerance = 0.02)               # matches flux_FINAL
  expect_equal(b$det_class_emp, "below detection"); expect_false(b$qc_co2_tracer)
  expect_equal(cc$det_class_emp, "below detection"); expect_true(cc$qc_co2_tracer)   # failed closure: no CO2 rise
  expect_true(d$qc_co2_tracer)                                     # ambiguous: CO2 flat
  expect_equal(f$det_class_emp, "uptake"); expect_false(f$qc_co2_tracer)
  # goFlux's own columns intact for all
  expect_true(all(c("best.flux", "MDF", "quality.check") %in% names(fx)))
})

test_that("example (e) carries real ebullition events", {
  e <- stemflux_examples("e")$data
  eb <- detect_ebullition(e, "CH4dry_ppb")                       # consecutive samples merged: one burst
  expect_equal(eb$summary$n_events, 1)
  expect_gt(eb$summary$total_ebullitive, 500)
  expect_equal(eb$summary$n_points, 278)
  raw <- detect_ebullition(e, "CH4dry_ppb", merge.consecutive = FALSE)
  expect_gte(raw$summary$n_events, 3)                             # the campaign script counted 5 samples
  expect_equal(raw$summary$total_ebullitive, eb$summary$total_ebullitive, tolerance = 1e-6)
  camp <- detect_ebullition(e, "CH4dry_ppb", method = "threshold", jump.thresh = 100, merge.consecutive = FALSE)
  expect_equal(camp$summary$n_events, 5)                          # reproduces the campaign's count
  expect_equal(camp$summary$total_ebullitive, 676, tolerance = 0.01)
})
