d <- example_closure()
best_ref <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))))

test_that("process_fluxes core path: appends only, goFlux columns untouched, events NULL", {
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb",
                                         qc = list(c0 = TRUE, convex = TRUE, min_window = list(secs = 120))))
  expect_s3_class(res, "fluxqc_result")
  expect_named(res, c("fluxes", "events", "traces", "gastype", "settings"))
  expect_null(res$events)
  fx <- res$fluxes
  expect_equal(nrow(fx), 1)
  # every goFlux column present and identical
  expect_length(setdiff(names(best_ref), names(fx)), 0)
  for (nm in names(best_ref)) expect_identical(fx[[nm]], as.data.frame(best_ref)[[nm]], info = nm)
  added <- setdiff(names(fx), names(best_ref))
  expect_true(all(grepl("_emp$|_emp_method$|^qc_", added)))
  expect_true(all(c("sigma_emp", "MDF_emp", "below_MDF_emp", "MDF_emp_method", "detected_emp",
                    "det_class_emp", "qc_c0", "qc_convex", "qc_any") %in% names(fx)))
  expect_false(any(grepl("ebullition|flux_total|flux_ebullitive|MDF_goflux_orig", names(fx))))
  expect_equal(fx$MDF_emp_method, "wassmann | mad | 0.95")
  expect_equal(res$settings$ebullition, "none")
  expect_null(res$settings$extra)
  expect_output(print(res), "fluxqc result")
  expect_false(grepl("ebullition", paste(capture.output(print(res)), collapse = " ")))
  # extras through the wrapper
  r2 <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", qc = FALSE,
                                        extra = list(precision = "datasheet", mdf = "goflux")))
  expect_true("MDF_goflux_datasheet" %in% names(r2$fluxes))
  expect_false("qc_any" %in% names(r2$fluxes))
  expect_null(r2$settings$qc)
})

test_that("goFlux::flux.plot runs unchanged on process_fluxes()$fluxes", {
  skip_if_not_installed("ggplot2")
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", qc = list(c0 = TRUE)))
  pl <- suppressWarnings(goFlux::flux.plot(res$fluxes, res$traces, "CH4dry_ppb"))
  expect_length(pl, 1)
  expect_true(inherits(pl[[1]], "ggplot"))
  ref <- suppressWarnings(goFlux::flux.plot(best_ref, d, "CH4dry_ppb"))
  expect_equal(fluxqc:::plot_uid(pl[[1]], 1), fluxqc:::plot_uid(ref[[1]], 1))
})

test_that("process_fluxes ebullition add-on: flag and split, per-closure control", {
  d2 <- d; d2$UniqueID <- "bubbly"
  d2$CH4dry_ppb[d2$t_s >= 300] <- d2$CH4dry_ppb[d2$t_s >= 300] + 60
  both <- rbind(d, d2)
  fl <- suppressWarnings(process_fluxes(both, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "flag"))
  expect_true(all(c("n_ebullition_events", "ebullition_flag", "ebullition_handling") %in% names(fl$fluxes)))
  expect_false("flux_total" %in% names(fl$fluxes))
  expect_equal(fl$fluxes$n_ebullition_events[fl$fluxes$UniqueID == "bubbly"], 1)
  expect_equal(nrow(fl$events), 1)
  aux <- data.frame(UniqueID = c("emission_stem_semirigid", "bubbly"), ebullition = c(FALSE, TRUE))
  sp <- suppressWarnings(process_fluxes(both, aux = aux, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "split"))
  fx <- sp$fluxes
  expect_equal(fx$ebullition_handling[match(c("emission_stem_semirigid", "bubbly"), fx$UniqueID)], c("none", "split"))
  expect_true(all(c("flux_diffusive", "flux_ebullitive", "flux_total", "ebullition_method") %in% names(fx)))
  b <- fx[fx$UniqueID == "bubbly", ]; o <- fx[fx$UniqueID != "bubbly", ]
  expect_gt(b$flux_ebullitive, 0); expect_equal(o$flux_ebullitive, 0)
  expect_equal(fx$flux_total, fx$flux_diffusive + fx$flux_ebullitive)
  expect_equal(b$flux_diffusive, b$best.flux)              # best.flux itself untouched
  expect_equal(b$ebullition_method, "iqr"); expect_true(is.na(o$ebullition_method))
  expect_equal(b$flux_diffusive, o$flux_diffusive, tolerance = 0.1)
  expect_equal(sp$settings$ebullition_ids, "bubbly")
  sp2 <- suppressWarnings(process_fluxes(both, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "split",
                                         ebullition_ids = "bubbly", ebullition_method = "threshold",
                                         ebullition_args = list(jump.thresh = 30)))
  expect_equal(sp2$fluxes$ebullition_method[sp2$fluxes$UniqueID == "bubbly"], "threshold")
})

test_that("write_outputs writes the documented files", {
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", qc = list(c0 = TRUE)))
  dir <- tempfile("out_")
  p <- write_outputs(res, dir, plots = requireNamespace("ggplot2", quietly = TRUE))
  expect_true(all(c("fluxes", "settings") %in% names(p)))
  expect_false("events" %in% names(p))
  expect_true(file.exists(file.path(dir, "fluxes.csv")))
  js <- readLines(file.path(dir, "settings.json"))
  expect_true(any(grepl("\"ebullition\": \"none\"", js)))
  expect_true(any(grepl("\"mdf\": \"wassmann\"", js)))
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_true(file.exists(file.path(dir, "flux_plots.pdf")))
  expect_equal(nrow(read.csv(file.path(dir, "fluxes.csv"))), 1)
  d2 <- d; d2$CH4dry_ppb[d2$t_s >= 300] <- d2$CH4dry_ppb[d2$t_s >= 300] + 60
  r2 <- suppressWarnings(process_fluxes(d2, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "split"))
  p2 <- write_outputs(r2, tempfile("out2_"), plots = FALSE)
  expect_true(file.exists(p2["events"]))
  expect_equal(nrow(read.csv(p2["events"])), 1)
  expect_error(write_outputs(list(), dir), "process_fluxes")
})
