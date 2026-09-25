keys <- c("emission_stem_semirigid", "emission_upland_stem", "below_detection_stem", "failed_closure_stem",
          "ambiguous_stem", "uptake_stem", "uptake_soil", "high_flux_wetland_stem", "ebullition_floating",
          "emission_li7810_stem")

test_that("bundled examples load, anonymized, with real geometry and goFlux columns", {
  ex <- fluxqc_examples()
  expect_named(ex, c("data", "aux"))
  expect_equal(ex$aux$UniqueID, keys)
  expect_true(all(ex$aux$kind == "real"))
  expect_equal(sort(unique(ex$data$UniqueID)), sort(keys))
  need <- c("UniqueID", "POSIX.time", "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm", "start.time", "obs.length",
            "Area", "Vtot", "Tcham", "Pcham", "flag", "Etime", "start.time_corr", "end.time_corr",
            "obs.length_corr", "CH4_prec", "CO2_prec")
  expect_true(all(need %in% names(ex$data)))
  # anonymized: synthetic calendar, no identifiers in column names
  expect_true(all(format(ex$data$POSIX.time, "%Y") == "2000"))
  expect_false(any(c("plot", "tree_id", "species", "DATE_TIME", "File name") %in% names(ex$data)))
  # real geometry
  expect_equal(ex$aux$Area[1], 4254.121, tolerance = 1e-6); expect_equal(ex$aux$Vtot[1], 7.082698, tolerance = 1e-6)
  expect_equal(ex$aux$Area[ex$aux$UniqueID == "ebullition_floating"], 324.3)
  expect_equal(ex$aux$Area[ex$aux$UniqueID == "high_flux_wetland_stem"], 107.68)
  # shipped (clicked) flags and segment sizes match the diagnostics table
  dg <- read.csv(system.file("extdata", "fluxqc_examples_diagnostics.csv", package = "fluxqc"), stringsAsFactors = FALSE)
  expect_equal(dg$key, keys)
  fl <- tapply(ex$data$flag, ex$data$UniqueID, sum)[keys]; nn <- table(ex$data$UniqueID)[keys]
  expect_equal(as.vector(fl), dg$n_flagged); expect_equal(as.vector(nn), dg$n_segment)
  a <- ex$data[ex$data$UniqueID == "emission_stem_semirigid", ]
  expect_equal(as.numeric(difftime(a$start.time_corr[1], a$start.time[1], units = "secs")), 91)
  expect_equal(unique(a$obs.length_corr), 765)
  expect_equal(min(a$Etime[a$flag == 1]), 0)
  # two analyzer families, sampling intervals preserved
  expect_true(any(grepl("LI-7810", dg$analyzer)) && any(grepl("UGGA", dg$analyzer)))
  expect_equal(dg$sampling_interval_s[dg$key == "emission_upland_stem"], 5)
  expect_equal(dg$sampling_interval_s[dg$key == "high_flux_wetland_stem"], 1)
  # continuity of every clicked window (recomputed, not just read from the table)
  cont <- t(sapply(split(ex$data, ex$data$UniqueID)[keys], function(x) { w <- x[x$flag == 1, ]; w <- w[order(w$POSIX.time), ]
    dt <- as.numeric(diff(w$POSIX.time), units = "secs"); dc <- abs(diff(w$CH4dry_ppb))
    c(gap = max(dt) / median(dt), step = max(dc) / mad(dc)) }))
  expect_true(all(cont[, "gap"] <= 2))
  expect_true(all(cont[keys != "ebullition_floating", "step"] <= 6))   # Gaussian extreme over 200-800 differences
  expect_gt(cont["ebullition_floating", "step"], 20)                    # the burst
  expect_equal(dg$max_gap_s, unname(cont[, "gap"] * dg$sampling_interval_s), tolerance = 0.05)
  # a numeric dead.band re-flags from start.time instead
  re <- fluxqc_examples("emission_stem_semirigid", dead.band = 60)$data
  expect_equal(sum(re$flag), sum(re$t_s >= 60))
  expect_equal(nrow(fluxqc_examples(c("uptake_soil", "uptake_stem"))$aux), 2)
  expect_error(fluxqc_examples("z"), "no example")
})

test_that("example_closure() is the alias for emission_stem_semirigid", {
  d <- example_closure()
  expect_equal(unique(d$UniqueID), "emission_stem_semirigid")
  expect_equal(nrow(d), 1207)
  expect_equal(sum(d$flag), 770)
  expect_equal(d$Tcham[1], 26.02)
})

test_that("the teaching cases separate through process_fluxes (Figure 4 logic)", {
  ex <- fluxqc_examples(c("emission_stem_semirigid", "below_detection_stem", "failed_closure_stem",
                            "ambiguous_stem", "uptake_stem", "uptake_soil", "emission_upland_stem"))
  res <- suppressWarnings(process_fluxes(ex$data, gastype = "CH4dry_ppb",
                                         qc = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE, min_window = list(secs = 120), noisy = TRUE)))
  fx <- res$fluxes; rownames(fx) <- fx$UniqueID
  expect_equal(fx["emission_stem_semirigid", "det_class_emp"], "emission"); expect_false(fx["emission_stem_semirigid", "qc_co2_tracer"])
  expect_equal(fx["emission_stem_semirigid", "best.flux"], 0.254, tolerance = 0.02)
  expect_equal(fx["emission_upland_stem", "det_class_emp"], "emission"); expect_false(fx["emission_upland_stem", "qc_co2_tracer"])
  expect_equal(fx["below_detection_stem", "det_class_emp"], "below detection"); expect_false(fx["below_detection_stem", "qc_co2_tracer"])
  expect_equal(fx["failed_closure_stem", "det_class_emp"], "below detection"); expect_true(fx["failed_closure_stem", "qc_co2_tracer"])
  expect_true(fx["ambiguous_stem", "qc_co2_tracer"])
  expect_equal(fx["uptake_stem", "det_class_emp"], "uptake"); expect_false(fx["uptake_stem", "qc_co2_tracer"])
  expect_equal(fx["uptake_soil", "det_class_emp"], "uptake"); expect_lt(fx["uptake_soil", "best.flux"], -1)
  expect_true(all(c("best.flux", "MDF", "quality.check") %in% names(fx)))
  # noisy screen (default in process_fluxes): which of the teaching cases are noisier than the set sigma
  expect_true("qc_noisy" %in% names(fx))
  expect_true(all(!is.na(fx$qc_noisy)))
})

test_that("the two other analyzer formats run through goFlux", {
  ex <- fluxqc_examples(c("high_flux_wetland_stem", "emission_li7810_stem"))
  res <- suppressWarnings(process_fluxes(ex$data, gastype = "CH4dry_ppb", qc = FALSE))
  fx <- res$fluxes; rownames(fx) <- fx$UniqueID
  expect_gt(fx["high_flux_wetland_stem", "best.flux"], 15)
  expect_equal(fx["emission_li7810_stem", "best.flux"], 0.653, tolerance = 0.05)
})

test_that("ebullition_floating carries one real bubble burst (goFlux::goAquaFlux)", {
  e <- fluxqc_examples("ebullition_floating")$data
  invisible(capture.output(r <- suppressWarnings(goFlux::goAquaFlux(e, "CH4dry_ppb"))))
  expect_equal(nrow(r$bubbles), 1L)
  expect_equal(r$bubbles$magnitude, 1089, tolerance = 0.02)   # goFlux's settled step, ppb
  expect_gt(r$flux_summary$flux_ebullition, 15)
  expect_true(r$bubbles$start < 30)                           # the burst begins ~42 s; window opens earlier
})
