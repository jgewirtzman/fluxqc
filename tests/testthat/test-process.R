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

test_that("process_fluxes ebullition add-on dispatches to goAquaFlux (flag and split, per-closure control)", {
  ex <- fluxqc_examples(c("emission_stem_semirigid", "ebullition_floating")); both <- ex$data
  invisible(capture.output(fl <- suppressWarnings(process_fluxes(both, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "flag"))))
  expect_true(all(c("n_ebullition_events", "ebullition_flag", "ebullition_handling") %in% names(fl$fluxes)))
  expect_false("flux_total" %in% names(fl$fluxes))
  expect_equal(fl$fluxes$n_ebullition_events[fl$fluxes$UniqueID == "ebullition_floating"], 1L)
  expect_equal(nrow(fl$events), 1)
  expect_true(all(c("start", "end", "magnitude") %in% names(fl$events)))
  # goFlux columns untouched: best.flux is still the whole-trace fit
  ref <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(both, "CH4dry_ppb"))))
  expect_equal(fl$fluxes$best.flux[match(ref$UniqueID, fl$fluxes$UniqueID)], ref$best.flux)
  aux <- data.frame(UniqueID = c("emission_stem_semirigid", "ebullition_floating"), ebullition = c(FALSE, TRUE))
  invisible(capture.output(sp <- suppressWarnings(process_fluxes(both, aux = aux, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "split"))))
  fx <- sp$fluxes
  expect_equal(fx$ebullition_handling[match(c("emission_stem_semirigid", "ebullition_floating"), fx$UniqueID)], c("none", "split"))
  expect_true(all(c("flux_diffusive", "flux_ebullitive", "flux_total", "SE_ebullitive", "n_obs_diffusive",
                    "first_bubble_time", "ebullition_method") %in% names(fx)))
  b <- fx[fx$UniqueID == "ebullition_floating", ]; o <- fx[fx$UniqueID != "ebullition_floating", ]
  expect_gt(b$flux_ebullitive, 15); expect_equal(o$flux_ebullitive, 0)
  expect_true(is.na(b$flux_diffusive))                    # early burst: goAquaFlux has too few pre-bubble points
  expect_true(is.na(o$flux_diffusive))                    # not in ebullition_ids: no partition columns
  expect_equal(b$ebullition_method, "goAquaFlux"); expect_true(is.na(o$ebullition_method))
  expect_equal(sp$settings$ebullition_ids, "ebullition_floating")
  invisible(capture.output(sp2 <- suppressWarnings(process_fluxes(both, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "split",
                                         ebullition_ids = "emission_stem_semirigid"))))
  s2 <- sp2$fluxes[sp2$fluxes$UniqueID == "emission_stem_semirigid", ]
  expect_equal(s2$n_ebullition_events, 0L); expect_equal(s2$flux_total, s2$flux_diffusive)
  expect_equal(s2$flux_diffusive, s2$best.flux, tolerance = 0.01)
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
  d2 <- fluxqc_examples("ebullition_floating")$data
  invisible(capture.output(r2 <- suppressWarnings(process_fluxes(d2, gastype = "CH4dry_ppb", qc = FALSE, ebullition = "split"))))
  p2 <- write_outputs(r2, tempfile("out2_"), plots = FALSE)
  expect_true(file.exists(p2["events"]))
  expect_equal(nrow(read.csv(p2["events"])), 1)
  expect_error(write_outputs(list(), dir), "process_fluxes")
})
