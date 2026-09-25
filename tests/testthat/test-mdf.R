d <- example_closure()
res <- suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))
best <- suppressWarnings(goFlux::best.flux(res))
ft <- best$flux.term; t <- best$nb.obs

test_that("mdf_multi reproduces each formula with two-sided quantiles", {
  m <- mdf_multi(flux.term = ft, t = 540, prec = 0.9, sigma = 2.2, sigma.global = 1.2)
  expect_equal(m$MDF_goflux, 0.9 / 540 * ft)
  expect_equal(m$MDF_wass90, qnorm(0.95) * 1.2 / 540 * ft)      # "90%" = qnorm(0.95), as in the FINAL script
  expect_equal(m$MDF_wass95, qnorm(0.975) * 1.2 / 540 * ft)
  expect_equal(m$MDF_chr99, 3 * 2.2 * qt(0.995, df = 538) / 540 * ft)
  expect_equal(m$MDF_2sigma, 2 * 2.2 / 540 * ft)
})

test_that("default = record-wide MAD, Wassmann form at 95% two-sided (z = 1.96, 2 sigma)", {
  out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb")
  added <- setdiff(names(out), names(best))
  expect_setequal(added, c("sigma_emp", "MDF_emp", "below_MDF_emp", "MDF_emp_method",
                           "detected_emp", "det_class_emp"))
  for (nm in names(best)) expect_identical(out[[nm]], best[[nm]], info = nm)
  expect_equal(out$sigma_emp, precision_mad(d$CH4dry_ppb))          # ALL rows, not just flagged
  expect_equal(out$MDF_emp, qnorm(0.975) * out$sigma_emp / t * ft)
  expect_equal(out$MDF_emp_method, "wassmann | mad | 0.95")
  expect_false(out$below_MDF_emp); expect_equal(out$det_class_emp, "emission")
})

test_that("mad is computed per group over the whole record", {
  d2 <- d; d2$UniqueID <- "second"; d2$CH4dry_ppb <- d2$CH4dry_ppb + rnorm(nrow(d2), 0, 5)
  both <- rbind(d, d2)
  b2 <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(both, "CH4dry_ppb"))))
  b2$camp <- c("x", "y")[match(b2$UniqueID, c("2021-05-18_5086_I_7.0", "second"))]
  g <- flag_detection(b2, traces = both, gastype = "CH4dry_ppb", group = "camp")
  expect_equal(g$sigma_emp[g$UniqueID == "2021-05-18_5086_I_7.0"], precision_mad(d$CH4dry_ppb))
  expect_equal(g$sigma_emp[g$UniqueID == "second"], precision_mad(d2$CH4dry_ppb))
  expect_gt(g$sigma_emp[g$UniqueID == "second"], g$sigma_emp[g$UniqueID != "second"])
  u <- flag_detection(b2, traces = both, gastype = "CH4dry_ppb")     # no group: one record
  expect_equal(u$sigma_emp[1], u$sigma_emp[2])
  expect_equal(u$sigma_emp[1], precision_mad(both$CH4dry_ppb))
  expect_error(flag_detection(b2, traces = both, gastype = "CH4dry_ppb", group = "nope"), "group")
})

test_that("other estimators and methods", {
  o <- flag_detection(best, precision = "datasheet", mdf = "goflux")
  expect_equal(o$MDF_emp, 0.9 / t * ft)
  expect_equal(o$MDF_emp, best$MDF, tolerance = 0.01)
  expect_equal(o$MDF_emp_method, "goflux | datasheet")
  o2 <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "allan", mdf = "christiansen", conf = 0.95)
  sa <- precision_allan(d$CH4dry_ppb[d$flag == 1])
  expect_equal(o2$sigma_emp, sa)
  expect_equal(o2$MDF_emp, 3 * sa * qt(0.975, t - 2) / t * ft)
  expect_equal(o2$MDF_emp_method, "christiansen | allan | 0.95")
  o3 <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "allan", mdf = "wassmann",
                       conf = 0.9, sigma.global = 1.2)
  expect_equal(o3$MDF_emp, qnorm(0.95) * 1.2 / t * ft)
  cu <- flag_detection(best, precision = "custom", mdf = "goflux", sigma = data.frame(UniqueID = best$UniqueID, sigma = 5))
  expect_equal(cu$MDF_emp, 5 / t * ft)
  o4 <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "rolling", mdf = "two_sigma")
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
  mdf_cols <- setdiff(grep("^MDF_[a-z]", names(out), value = TRUE), c("MDF_emp", "MDF_emp_method"))
  expect_length(mdf_cols, 18)
  expect_true(all(c("MDF_goflux_datasheet", "MDF_wassmann_allan_90", "MDF_christiansen_mad_95",
                    "MDF_two_sigma_allan", "sigma_mad", "sigma_datasheet") %in% names(out)))
  expect_equal(out$MDF_goflux_datasheet, 0.9 / t * ft)
  expect_equal(out$MDF_wassmann_mad_95, out$MDF_emp)
  expect_equal(out$MDF_wassmann_datasheet_90, qnorm(0.95) * 0.9 / t * ft)
  expect_equal(out$MDF_christiansen_allan_95, 3 * out$sigma_allan * qt(0.975, t - 2) / t * ft)
  all <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", extra = "all")
  expect_true(all(c("MDF_wassmann_rolling_99", "MDF_two_sigma_datasheet") %in% names(all)))
  expect_error(flag_detection(best, traces = d, gastype = "CH4dry_ppb", extra = 3), "extra")
})
