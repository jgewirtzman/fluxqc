d <- example_closure()
res <- suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))
best <- suppressWarnings(goFlux::best.flux(res))

test_that("mdf_multi reproduces each formula", {
  ft <- 0.6782307; t <- 540; prec <- 0.9; sig <- 2.2; sg <- 1.2
  m <- mdf_multi(flux.term = ft, t = t, prec = prec, sigma = sig, sigma.global = sg)
  expect_equal(m$MDF_goflux, prec / t * ft)
  expect_equal(m$MDF_wass95, qnorm(0.95) * sg / t * ft)
  expect_equal(m$MDF_chr99, 3 * sig * qt(0.99, df = t - 2) / t * ft)
  expect_equal(m$MDF_2sigma, 2 * sig / t * ft)
})

test_that("flag_detection default adds exactly the *_emp columns and touches nothing else", {
  out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb")
  added <- setdiff(names(out), names(best))
  expect_setequal(added, c("sigma_emp", "MDF_emp", "below_MDF_emp", "MDF_emp_method",
                           "detected_emp", "det_class_emp"))
  expect_length(setdiff(names(best), names(out)), 0)
  for (nm in names(best)) expect_identical(out[[nm]], best[[nm]], info = nm)
  ft <- best$flux.term; t <- best$nb.obs
  expect_equal(out$sigma_emp, precision_allan(d$CH4dry_ppb[d$flag == 1]))
  expect_equal(out$MDF_emp, 3 * out$sigma_emp * qt(0.95, t - 2) / t * ft)
  expect_equal(out$MDF_emp_method, "christiansen | allan | 0.95")
  expect_false(out$below_MDF_emp); expect_true(out$detected_emp); expect_equal(out$det_class_emp, "emission")
  expect_equal(unname(attr(out, "sigma_global")["allan"]), out$sigma_emp)
})

test_that("other estimators and methods", {
  ft <- best$flux.term; t <- best$nb.obs
  o <- flag_detection(best, precision = "datasheet", mdf = "goflux")
  expect_equal(o$MDF_emp, 0.9 / t * ft)
  expect_equal(o$MDF_emp, best$MDF, tolerance = 0.01)      # same formula as goFlux's own MDF
  expect_equal(o$MDF_emp_method, "goflux | datasheet")
  o2 <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "mad", mdf = "two_sigma")
  expect_equal(o2$MDF_emp, 2 * precision_mad(d$CH4dry_ppb[d$flag == 1]) / t * ft)
  o3 <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", mdf = "wassmann", conf = 0.9, sigma.global = 1.2)
  expect_equal(o3$MDF_emp, qnorm(0.9) * 1.2 / t * ft)
  expect_equal(o3$MDF_emp_method, "wassmann | allan | 0.9")
  cu <- flag_detection(best, precision = "custom", mdf = "goflux", sigma = data.frame(UniqueID = best$UniqueID, sigma = 5))
  expect_equal(cu$MDF_emp, 5 / t * ft)
  o4 <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "rolling")
  expect_true(is.finite(o4$sigma_emp))
  expect_error(flag_detection(best), "need 'traces'")
  expect_error(flag_detection(best, precision = "custom", mdf = "goflux"), "needs 'sigma'")
  expect_error(flag_detection(best, traces = d, gastype = "CH4dry_ppb", conf = c(0.9, 0.95)), "single")
  fake <- best; fake$best.flux <- -1e-4
  expect_equal(flag_detection(fake, traces = d, gastype = "CH4dry_ppb")$det_class_emp, "below detection")
  fake$best.flux <- -1
  expect_equal(flag_detection(fake, traces = d, gastype = "CH4dry_ppb")$det_class_emp, "uptake")
})

test_that("extra appends comparison columns with the documented names", {
  out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb",
                        extra = list(precision = c("allan", "mad", "datasheet"),
                                     mdf = c("goflux", "wassmann", "christiansen", "two_sigma"),
                                     conf = c(0.90, 0.95)))
  mdf_cols <- grep("^MDF_[a-z]", names(out), value = TRUE)
  mdf_cols <- setdiff(mdf_cols, c("MDF_emp", "MDF_emp_method"))
  # 3 precisions x (goflux + two_sigma + 2 x wassmann + 2 x christiansen) = 18
  expect_length(mdf_cols, 18)
  expect_true(all(c("MDF_goflux_datasheet", "MDF_wassmann_allan_90", "MDF_christiansen_mad_95",
                    "MDF_two_sigma_allan", "sigma_mad", "sigma_datasheet") %in% names(out)))
  ft <- best$flux.term; t <- best$nb.obs
  expect_equal(out$MDF_goflux_datasheet, 0.9 / t * ft)
  expect_equal(out$MDF_christiansen_allan_95, out$MDF_emp)
  expect_equal(out$MDF_wassmann_datasheet_90, qnorm(0.90) * 0.9 / t * ft)
  expect_false(any(grepl("^below_MDF_[a-z]+_", names(out))))
  all <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", extra = "all")
  expect_true(all(c("MDF_wassmann_rolling_99", "MDF_two_sigma_datasheet") %in% names(all)))
  expect_false("sigma_custom" %in% names(all))
  expect_error(flag_detection(best, traces = d, gastype = "CH4dry_ppb", extra = 3), "extra")
})
