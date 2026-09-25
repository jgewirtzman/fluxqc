d <- example_closure()
best_ref <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))))

test_that("process_fluxes core path: appends only, goFlux columns untouched, events NULL", {
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb",
                                         qc = list(c0 = TRUE, convex = TRUE, min_window = list(secs = 120))))
  expect_s3_class(res, "fluxqc_result")
  expect_named(res, c("fluxes", "traces", "gastype", "settings"))
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
  expect_equal(ggplot2::ggplot_build(pl[[1]])$data, ggplot2::ggplot_build(ref[[1]])$data)
})

test_that("write_outputs writes the documented files", {
  res <- suppressWarnings(process_fluxes(d, gastype = "CH4dry_ppb", qc = list(c0 = TRUE)))
  dir <- tempfile("out_")
  p <- write_outputs(res, dir, plots = requireNamespace("ggplot2", quietly = TRUE))
  expect_true(all(c("fluxes", "settings") %in% names(p)))
  expect_true(file.exists(file.path(dir, "fluxes.csv")))
  js <- readLines(file.path(dir, "settings.json"))
  expect_true(any(grepl("\"mdf\": \"wassmann\"", js)))
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_true(file.exists(file.path(dir, "flux_plots.pdf")))
  expect_equal(nrow(read.csv(file.path(dir, "fluxes.csv"))), 1)
  expect_error(write_outputs(list(), dir), "process_fluxes")
})
