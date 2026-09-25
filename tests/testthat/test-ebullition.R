make_trace <- function(jump_at = 300, jump = 400, n = 600) {
  set.seed(5)
  t0 <- as.POSIXct("2022-10-23 11:00:00", tz = "UTC")
  ch4 <- 2000 + 0.5 * (0:(n - 1)) + rnorm(n, 0, 1)
  if (!is.null(jump_at)) ch4[jump_at:n] <- ch4[jump_at:n] + jump
  data.frame(POSIX.time = t0 + 0:(n - 1), CH4dry_ppb = ch4)
}

test_that("detect_ebullition finds a single jump and de-ebulliates", {
  d <- make_trace()
  eb <- detect_ebullition(d, "CH4dry_ppb", jump.thresh = 100)
  expect_equal(eb$summary$n_jumps, 1)
  expect_equal(eb$summary$total_ebullitive, 400, tolerance = 0.02)
  expect_equal(eb$summary$jump_thresh_used, 100)   # dt = 1 s, no scaling
  expect_equal(nrow(eb$trace), 560)   # elapsed 20..579 of 0..599
  # de-ebulliated trace has no step
  expect_lt(max(abs(diff(eb$trace$C_deebull))), 10)
  expect_equal(sum(eb$trace$near_jump), 31)
  expect_true(eb$summary$ebullitive_fraction > 0.5)
  # threshold scales with logging interval
  d5 <- d[seq(1, 600, by = 5), ]
  expect_equal(detect_ebullition(d5, "CH4dry_ppb", jump.thresh = 100)$summary$jump_thresh_used, 100)
  d10 <- d[seq(1, 600, by = 10), ]
  expect_equal(detect_ebullition(d10, "CH4dry_ppb", jump.thresh = 100)$summary$jump_thresh_used, 200)
})

test_that("no jumps in a diffusive trace", {
  eb <- detect_ebullition(make_trace(NULL), "CH4dry_ppb", jump.thresh = 100)
  expect_equal(eb$summary$n_jumps, 0)
  expect_equal(nrow(eb$jumps), 0)
  expect_false(any(eb$trace$near_jump))
})

test_that("split_ebullition partitions diffusive and ebullitive flux", {
  d <- make_trace()
  eb <- detect_ebullition(d, "CH4dry_ppb", jump.thresh = 100)
  p <- split_ebullition(eb, Vtot = 4.318, Area = 324.3, Tcham = 28, Pcham = 101.325)
  ft <- flux_term(4.318, 101.325, 324.3, 28)
  expect_equal(p$flux.term, ft)
  expect_equal(p$diffusive_slope, 0.5, tolerance = 0.02)
  expect_equal(p$diffusive_flux, p$diffusive_slope * ft)
  expect_equal(p$ebullitive_rate, eb$summary$total_ebullitive / eb$summary$duration_sec)
  expect_equal(p$total_flux, p$diffusive_flux + p$ebullitive_flux)
  expect_gt(p$ebullitive_fraction, 0.5)
  p0 <- split_ebullition(detect_ebullition(make_trace(NULL), "CH4dry_ppb", jump.thresh = 100),
                         4.318, 324.3, 28, 101.325)
  expect_equal(p0$ebullitive_flux, 0)
  expect_equal(p0$n_jumps, 0)
})
