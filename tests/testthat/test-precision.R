test_that("precision_allan recovers white noise under a linear trend", {
  set.seed(10)
  x <- 2000 + 0.5 * (1:5000) + rnorm(5000, 0, 1.5)
  expect_equal(precision_allan(x), 1.5, tolerance = 0.05)
  expect_gt(sd(x), 100)
  expect_true(is.na(precision_allan(c(1, 2))))
  expect_equal(precision_allan(c(NA, 1, 2, 3, 4)), 0)
})

test_that("precision_mad is robust to gross outliers", {
  set.seed(11)
  x <- 2000 + rnorm(5000, 0, 1.2)
  x[sample(5000, 300)] <- x[sample(5000, 300)] + 60
  expect_gt(precision_allan(x), 5)
  # 6% contaminated points -> 12% contaminated differences; MAD is inflated only mildly
  expect_lt(abs(precision_mad(x) - 1.2), 0.25)
  # a constant drift cancels exactly
  expect_equal(precision_mad(x + 0.13 * seq_along(x)), precision_mad(x), tolerance = 1e-6)
})

test_that("precision_rolling finds the quiet part of a record", {
  set.seed(12)
  x <- c(2000 + rnorm(600, 0, 1), 2000 + 0.5 * (1:300) + rnorm(300, 0, 1))
  r <- precision_rolling(x, window = 30)
  expect_named(r, c("sd", "n_windows", "n_selected"))
  expect_equal(r$n_windows, 871)
  expect_lt(r$sd, 1.2)
  expect_warning(precision_rolling(1:10, window = 30), "too short")
})

test_that("precision_by_id summarises per UniqueID", {
  d <- example_closure()
  d2 <- d; d2$UniqueID <- "second"
  s <- precision_by_id(rbind(d, d2), "CH4dry_ppb")
  expect_equal(nrow(s), 2)
  expect_equal(s$n, c(542, 542))
  expect_equal(s$sigma[1], precision_allan(d$CH4dry_ppb[d$flag == 1]))
  sm <- precision_by_id(d, "CH4dry_ppb", method = "mad", flag.only = FALSE)
  expect_equal(sm$n, 603)
  expect_equal(sm$method, "mad")
  expect_error(precision_by_id(d, "nope"), "not found")
})
