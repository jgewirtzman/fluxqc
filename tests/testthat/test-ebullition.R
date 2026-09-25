ex <- fluxqc_examples(c("ebullition_floating", "emission_stem_semirigid"))

test_that("detect_ebullition wraps goAquaFlux and returns goFlux's bubbles", {
  invisible(capture.output(eb <- suppressWarnings(detect_ebullition(ex$data, "CH4dry_ppb"))))
  expect_s3_class(eb, "fluxqc_ebullition")
  expect_named(eb, c("events", "summary", "diffusive", "gastype"))
  s <- eb$summary; rownames(s) <- s$UniqueID
  expect_equal(s["ebullition_floating", "n_events"], 1L)
  expect_equal(s["ebullition_floating", "total_ebullitive"], 1089, tolerance = 0.02)   # goFlux's settled step (ppb)
  expect_gt(s["ebullition_floating", "flux_ebullitive"], 15)
  expect_true(is.na(s["ebullition_floating", "flux_diffusive"]))                      # bubble at 23 s: too few pre-bubble obs
  expect_equal(s["emission_stem_semirigid", "n_events"], 0L)
  expect_equal(s["emission_stem_semirigid", "flux_ebullitive"], 0)
  expect_equal(s["emission_stem_semirigid", "flux_diffusive"], 0.254, tolerance = 0.01)
  expect_equal(s["emission_stem_semirigid", "flux_total"], s["emission_stem_semirigid", "flux_diffusive"])
  expect_equal(nrow(eb$events), 1)
  expect_equal(eb$events$UniqueID, "ebullition_floating")
  expect_true(all(c("start", "end", "magnitude", "SE") %in% names(eb$events)))
  expect_true(eb$events$start < 30 && eb$events$end > 40)
  expect_true("best.flux" %in% names(eb$diffusive))
  expect_error(detect_ebullition(ex$data[, setdiff(names(ex$data), "flag")], "CH4dry_ppb"), "flag")
})

test_that("goAquaFlux arguments pass through", {
  invisible(capture.output(eb <- suppressWarnings(detect_ebullition(ex$data, "CH4dry_ppb",
                                                                     use_bubble_detection = FALSE))))
  expect_equal(eb$summary$n_events, c(0L, 0L))
  expect_equal(nrow(eb$events), 0)
})
