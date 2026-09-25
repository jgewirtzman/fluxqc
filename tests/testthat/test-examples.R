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
  expect_equal(as.vector(table(ex$data$UniqueID)[ex$aux$UniqueID]), c(603L, 603L, 469L, 603L, 278L, 603L))
  # dead band: 60 s for a-d,f; 0 for e (already trimmed by the source script)
  fl <- tapply(ex$data$flag, ex$data$UniqueID, sum)[ex$aux$UniqueID]
  expect_equal(unname(fl[5]), 278)
  expect_lt(unname(fl[1]), 603)
  one <- stemflux_examples("c")
  expect_equal(unique(one$data$UniqueID), "2021-02-25_1641_WD_5")
  expect_equal(nrow(stemflux_examples("2020-09-21_3589_U_6")$aux), 1)
  expect_error(stemflux_examples("z"), "no example")
})

test_that("example_closure() is the alias for (a)", {
  d <- example_closure()
  expect_equal(unique(d$UniqueID), "2021-05-18_5086_I_7.0")
  expect_equal(nrow(d), 603)
  expect_equal(sum(d$flag), sum(d$t_s >= 60))
  expect_equal(d$Tcham[1], 26.02)
})

test_that("the six examples separate as intended through process_fluxes (Figure 4 logic)", {
  ex <- stemflux_examples(c("a", "b", "c", "d", "f"))
  res <- suppressWarnings(process_fluxes(ex$data, gastype = "CH4dry_ppb",
                                         qc = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE, min_window = list(secs = 120))))
  fx <- res$fluxes; rownames(fx) <- fx$UniqueID
  a <- fx["2021-05-18_5086_I_7.0", ]; b <- fx["2020-09-23_5092_I_6", ]; cc <- fx["2021-02-25_1641_WD_5", ]
  d <- fx["2020-09-21_3589_U_6", ]; f <- fx["SYN_uptake_from_a", ]
  expect_equal(a$det_class_emp, "emission"); expect_false(a$qc_co2_tracer)
  expect_equal(b$det_class_emp, "below detection"); expect_false(b$qc_co2_tracer)
  expect_true(cc$qc_co2_tracer)                                    # failed closure: no CO2 rise
  expect_true(d$detected_emp); expect_true(d$qc_co2_tracer)        # ambiguous: CH4 rises, CO2 flat
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
