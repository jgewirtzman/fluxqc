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
  m2 <- mdf_multi(flux.term = c(ft, ft), t = c(300, 600), prec = 0.9)
  expect_true(all(is.na(m2$MDF_2sigma)))
  expect_equal(m2$MDF_goflux[1], 2 * m2$MDF_goflux[2])
})

test_that("flag_detection computes every requested combination with documented names", {
  out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb",
                        precision = c("allan", "mad", "rolling", "datasheet"),
                        mdf = c("goflux", "wassmann", "christiansen", "two_sigma"),
                        conf = c(0.90, 0.95), primary = list(precision = "allan", mdf = "christiansen", conf = 0.95))
  expect_true(all(c("sigma_allan", "sigma_mad", "sigma_rolling", "sigma_datasheet") %in% names(out)))
  expect_equal(out$sigma_datasheet, 0.9)
  expect_equal(out$sigma_allan, precision_allan(d$CH4dry_ppb[d$flag == 1]))
  expect_equal(out$sigma_mad, precision_mad(d$CH4dry_ppb[d$flag == 1]))
  mdf_cols <- grep("^MDF_", names(out), value = TRUE)
  # 4 precisions x (goflux + two_sigma + 2 conf x wassmann + 2 conf x christiansen) = 24, + orig
  expect_equal(length(setdiff(mdf_cols, c("MDF_goflux_orig", "MDF_method"))), 24)
  expect_true(all(c("MDF_goflux_datasheet", "MDF_wassmann_allan_90", "MDF_christiansen_allan_95",
                    "MDF_two_sigma_mad", "below_christiansen_allan_95") %in% names(out)))
  # formulas
  ft <- best$flux.term; t <- best$nb.obs
  expect_equal(out$MDF_goflux_datasheet, 0.9 / t * ft)
  expect_equal(out$MDF_two_sigma_allan, 2 * out$sigma_allan / t * ft)
  expect_equal(out$MDF_christiansen_allan_95, 3 * out$sigma_allan * qt(0.95, t - 2) / t * ft)
  expect_equal(out$MDF_wassmann_datasheet_90, qnorm(0.90) * 0.9 / t * ft)  # global = median = prec
  # primary
  expect_equal(out$MDF, out$MDF_christiansen_allan_95)
  expect_equal(out$MDF_method, "christiansen_allan_95")
  expect_false(out$below_MDF); expect_true(out$detected); expect_equal(out$det_class, "emission")
  expect_equal(out$MDF_goflux_orig, best$MDF)
  expect_named(attr(out, "sigma_global"), c("allan", "mad", "rolling", "datasheet"))
})

test_that("flag_detection primary selection, custom sigma and sigma.global", {
  out <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "allan",
                        mdf = "two_sigma", primary = list(precision = "allan", mdf = "two_sigma"))
  expect_equal(out$MDF_method, "two_sigma_allan")
  expect_error(flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                              primary = list(precision = "datasheet", mdf = "goflux")), "not computed")
  cu <- flag_detection(best, precision = "datasheet", mdf = "goflux", sigma = data.frame(UniqueID = best$UniqueID, sigma = 5),
                       primary = list(precision = "custom", mdf = "goflux"))
  expect_equal(cu$MDF, 5 / best$nb.obs * best$flux.term)
  w <- flag_detection(best, traces = d, gastype = "CH4dry_ppb", precision = "allan", mdf = "wassmann",
                      conf = 0.95, sigma.global = 1.2, primary = list(precision = "allan", mdf = "wassmann", conf = 0.95))
  expect_equal(w$MDF, qnorm(0.95) * 1.2 / best$nb.obs * best$flux.term)
})

test_that("flag_detection without traces falls back to datasheet and errors when explicit", {
  expect_message(o <- flag_detection(best, primary = list(precision = "datasheet", mdf = "goflux")), "datasheet")
  expect_equal(o$MDF_method, "goflux_datasheet")
  expect_equal(o$MDF, best$MDF, tolerance = 0.01)
  expect_error(flag_detection(best, precision = "allan"), "need 'traces'")
  fake <- best; fake$best.flux <- -1e-4
  expect_equal(flag_detection(fake, traces = d, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                              primary = list(precision = "allan", mdf = "two_sigma"))$det_class, "below detection")
  fake$best.flux <- -1
  expect_equal(flag_detection(fake, traces = d, gastype = "CH4dry_ppb", precision = "allan", mdf = "two_sigma",
                              primary = list(precision = "allan", mdf = "two_sigma"))$det_class, "uptake")
})
