test_that("mdf_multi reproduces each formula", {
  ft <- 0.6782307; t <- 540; prec <- 0.9; sig <- 2.2; sg <- 1.2
  m <- mdf_multi(flux.term = ft, t = t, prec = prec, sigma = sig, sigma.global = sg)
  expect_equal(m$MDF_goflux, prec / t * ft)
  expect_equal(m$MDF_wass95, qnorm(0.95) * sg / t * ft)
  expect_equal(m$MDF_chr99, 3 * sig * qt(0.99, df = t - 2) / t * ft)
  expect_equal(m$MDF_2sigma, 2 * sig / t * ft)
  expect_true(all(m$MDF_wass90 < m$MDF_wass95, m$MDF_wass95 < m$MDF_wass99))
  # missing inputs give NA, vectorised over rows
  m2 <- mdf_multi(flux.term = c(ft, ft), t = c(300, 600), prec = 0.9)
  expect_equal(nrow(m2), 2)
  expect_true(all(is.na(m2$MDF_2sigma)))
  expect_equal(m2$MDF_goflux[1], 2 * m2$MDF_goflux[2])
})

test_that("mdf_multi agrees with goFlux's MDF at t = nb.obs", {
  d <- example_closure()
  res <- suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))
  m <- mdf_multi(res$flux.term, t = max(d$Etime[d$flag == 1]) + 1, prec = res$prec)
  expect_equal(m$MDF_goflux, res$MDF, tolerance = 1e-6)
})

test_that("flag_detection appends MDFs and flags", {
  d <- example_closure()
  res <- suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))
  best <- suppressWarnings(goFlux::best.flux(res))
  sig <- precision_by_id(d, "CH4dry_ppb")
  out <- flag_detection(best, sigma = sig)
  expect_true(all(c("MDF_goflux", "MDF_wass95", "MDF_chr95", "MDF_2sigma", "below_2sigma",
                    "detected", "det_class", "sigma") %in% names(out)))
  expect_equal(out$sigma, sig$sigma)
  expect_true(out$detected)
  expect_equal(out$det_class, "emission")
  expect_equal(out$det_method, "2sigma")
  expect_false(out$below_chr99)
  # a tiny flux is below detection / uptake
  fake <- best; fake$best.flux <- -1e-4
  o2 <- flag_detection(fake, sigma = sig$sigma)
  expect_equal(o2$det_class, "below detection")
  fake$best.flux <- -1
  expect_equal(flag_detection(fake, sigma = sig$sigma)$det_class, "uptake")
  # no sigma: goflux method only
  o3 <- flag_detection(best)
  expect_equal(o3$det_method, "goflux")
  expect_true(is.na(o3$MDF_2sigma))
  expect_error(flag_detection(best, sigma = sig, primary = "nope"), "unknown")
  expect_error(flag_detection(best[, setdiff(names(best), "flux.term")]), "flux.term")
})
