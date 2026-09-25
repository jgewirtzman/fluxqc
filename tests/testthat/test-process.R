d <- example_closure()

test_that("process_fluxes core path: no ebullition columns, events NULL", {
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", precision = c("allan", "datasheet"),
                                         mdf = c("goflux", "christiansen"), conf = 0.95,
                                         qc = list(c0 = TRUE, convex = TRUE, min_window = list(secs = 120))))
  expect_s3_class(res, "stemflux_result")
  expect_named(res, c("fluxes", "events", "traces", "gastype", "settings"))
  expect_null(res$events)
  fx <- res$fluxes
  expect_equal(nrow(fx), 1)
  expect_true(all(c("best.flux", "MDF", "MDF_method", "det_class", "qc_c0", "qc_convex", "qc_any") %in% names(fx)))
  expect_false(any(grepl("ebullition|flux_total|flux_ebullitive", names(fx))))
  expect_equal(fx$MDF_method, "christiansen_allan_95")
  expect_equal(res$settings$ebullition, "none")
  expect_true(is.na(res$settings$ebullition_method))
  expect_equal(res$settings$precision, c("allan", "datasheet"))
  expect_output(print(res), "stemflux result")
  expect_false(grepl("ebullition", paste(capture.output(print(res)), collapse = " ")))
  # qc = FALSE skips screens
  r2 <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                                        primary = list(precision = "allan", mdf = "two_sigma"), qc = FALSE))
  expect_false("qc_any" %in% names(r2$fluxes))
  expect_null(r2$settings$qc)
})

test_that("process_fluxes ebullition add-on: flag and split, per-closure control", {
  d2 <- d; d2$UniqueID <- "bubbly"
  d2$CH4dry_ppb[d2$t_s >= 300] <- d2$CH4dry_ppb[d2$t_s >= 300] + 60
  both <- rbind(d, d2)
  fl <- suppressWarnings(process_fluxes(both, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                                        primary = list(precision = "allan", mdf = "two_sigma"), qc = FALSE,
                                        ebullition = "flag"))
  expect_true(all(c("n_ebullition_events", "ebullition_flag", "ebullition_handling") %in% names(fl$fluxes)))
  expect_false("flux_total" %in% names(fl$fluxes))
  expect_equal(fl$fluxes$n_ebullition_events[fl$fluxes$UniqueID == "bubbly"], 1)
  expect_equal(nrow(fl$events), 1)
  # split only for the flagged closure, via aux column
  aux <- data.frame(UniqueID = c("2021-05-18_5086_I_7.0", "bubbly"), ebullition = c(FALSE, TRUE))
  sp <- suppressWarnings(process_fluxes(both, aux = aux, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                                        primary = list(precision = "allan", mdf = "two_sigma"), qc = FALSE,
                                        ebullition = "split"))
  fx <- sp$fluxes
  expect_equal(fx$ebullition_handling[order(fx$UniqueID)], c("none", "split"))
  expect_true(all(c("flux_diffusive", "flux_ebullitive", "flux_total", "ebullition_method") %in% names(fx)))
  b <- fx[fx$UniqueID == "bubbly", ]; o <- fx[fx$UniqueID != "bubbly", ]
  expect_gt(b$flux_ebullitive, 0); expect_equal(o$flux_ebullitive, 0)
  expect_equal(fx$flux_total, fx$flux_diffusive + fx$flux_ebullitive)
  expect_equal(b$ebullition_method, "iqr"); expect_true(is.na(o$ebullition_method))
  # diffusive term on the de-ebulliated trace is close to the clean closure's flux
  expect_equal(b$flux_diffusive, o$flux_diffusive, tolerance = 0.1)
  expect_equal(sp$settings$ebullition_ids, "bubbly")
  expect_true(all(sp$traces$CH4dry_ppb[sp$traces$UniqueID == "bubbly"] < 2400))
  # ids vector instead of aux
  sp2 <- suppressWarnings(process_fluxes(both, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                                         primary = list(precision = "allan", mdf = "two_sigma"), qc = FALSE,
                                         ebullition = "split", ebullition_ids = "bubbly", ebullition_method = "threshold",
                                         ebullition_args = list(jump.thresh = 30)))
  expect_equal(sp2$fluxes$ebullition_method[sp2$fluxes$UniqueID == "bubbly"], "threshold")
})

test_that("write_outputs writes the documented files", {
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", precision = "allan", mdf = "christiansen",
                                         conf = 0.95, qc = list(c0 = TRUE)))
  dir <- tempfile("out_")
  p <- write_outputs(res, dir, plots = requireNamespace("ggplot2", quietly = TRUE))
  expect_true(all(c("fluxes", "settings") %in% names(p)))
  expect_false("events" %in% names(p))
  expect_true(file.exists(file.path(dir, "fluxes.csv")))
  js <- readLines(file.path(dir, "settings.json"))
  expect_true(any(grepl("\"ebullition\": \"none\"", js)))
  expect_true(any(grepl("christiansen", js)))
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_true(file.exists(file.path(dir, "flux_plots.pdf")))
  fx <- read.csv(file.path(dir, "fluxes.csv"))
  expect_equal(nrow(fx), 1)
  # with the add-on, the events file appears
  d2 <- d; d2$CH4dry_ppb[d2$t_s >= 300] <- d2$CH4dry_ppb[d2$t_s >= 300] + 60
  r2 <- suppressWarnings(process_fluxes(d2, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                                        primary = list(precision = "allan", mdf = "two_sigma"), qc = FALSE,
                                        ebullition = "split"))
  p2 <- write_outputs(r2, tempfile("out2_"), plots = FALSE)
  expect_true(file.exists(p2["events"]))
  expect_equal(nrow(read.csv(p2["events"])), 1)
  expect_error(write_outputs(list(), dir), "process_fluxes")
})
