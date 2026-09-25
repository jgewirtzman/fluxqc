make_trace <- function(uid = "w1", jump_at = 300, jump = 400, n = 600, seed = 5) {
  set.seed(seed)
  t0 <- as.POSIXct("2022-10-23 11:00:00", tz = "UTC")
  ch4 <- 2000 + 0.5 * (0:(n - 1)) + rnorm(n, 0, 1)
  if (!is.null(jump_at)) ch4[jump_at:n] <- ch4[jump_at:n] + jump
  data.frame(UniqueID = uid, POSIX.time = t0 + 0:(n - 1), CH4dry_ppb = ch4, H2O_ppm = 0, CH4_prec = 0.9,
             flag = 1, Etime = 0:(n - 1), Vtot = 4.318, Area = 324.3, Tcham = 28, Pcham = 101.325)
}

test_that("iqr method (Hoffmann fences) finds one injected jump per closure", {
  d <- rbind(make_trace("w1"), make_trace("w2", jump_at = 200, jump = 250, seed = 6),
             make_trace("w3", jump_at = NULL, seed = 8))
  eb <- detect_ebullition(d, "CH4dry_ppb", method = "iqr")
  expect_equal(eb$summary$method, rep("iqr", 3))
  # k = 1.5 (conventional fence) lets Gaussian noise through; k = 3 (default) does not
  expect_gt(sum(detect_ebullition(d, "CH4dry_ppb", method = "iqr", k = 1.5)$summary$n_events), 2)
  expect_s3_class(eb, "fluxqc_ebullition")
  expect_equal(eb$summary$n_events, c(1, 1, 0))
  expect_equal(eb$summary$total_ebullitive[1:2], c(400, 250), tolerance = 0.02)
  expect_equal(nrow(eb$events), 2)
  expect_equal(eb$events$UniqueID, c("w1", "w2"))
  expect_equal(eb$events$Etime[1], 299)
  # de-ebulliated trace has no step; buffer marked
  tr1 <- eb$trace[eb$trace$UniqueID == "w1", ]
  expect_lt(max(abs(diff(tr1$C_deebull))), 10)
  expect_equal(sum(tr1$near_jump), 31)
  # Hoffmann's own fence multiplier is allowed
  eb2 <- detect_ebullition(d, "CH4dry_ppb", method = "iqr", k = 0.25)
  expect_gte(eb2$summary$n_events[1], 1)
})

test_that("threshold and runvar methods", {
  d <- make_trace()
  th <- detect_ebullition(d, "CH4dry_ppb", method = "threshold", jump.thresh = 100)
  expect_equal(th$summary$n_events, 1)
  expect_equal(th$summary$threshold_used, 100)
  d10 <- d[seq(1, 600, by = 10), ]
  expect_equal(detect_ebullition(d10, "CH4dry_ppb", method = "threshold", jump.thresh = 100)$summary$threshold_used, 200)
  rv <- detect_ebullition(d, "CH4dry_ppb", method = "runvar", runvar_window = 5, runvar_cutoff = 50)
  expect_equal(rv$summary$n_events, 1)
  expect_equal(rv$summary$total_ebullitive, 400, tolerance = 0.05)
  expect_true("runvar" %in% names(rv$trace))
  none <- detect_ebullition(make_trace(jump_at = NULL), "CH4dry_ppb", method = "runvar", runvar_cutoff = 50)
  expect_equal(none$summary$n_events, 0)
  expect_equal(nrow(none$events), 0)
})

test_that("split_ebullition yields goFlux-ready data and merge appends totals", {
  d <- rbind(make_trace("w1"), make_trace("w3", jump_at = NULL, seed = 8))
  eb <- detect_ebullition(d, "CH4dry_ppb")
  sp <- split_ebullition(eb)
  expect_true(all(c("CH4dry_ppb", "CH4dry_ppb_raw", "flag", "Etime", "Vtot") %in% names(sp$diffusive)))
  expect_false(any(c("C_deebull", "is_jump") %in% names(sp$diffusive)))
  expect_equal(sum(sp$diffusive$flag == 0), 31)
  ft <- flux_term(4.318, 101.325, 324.3, 28)
  expect_equal(sp$ebullitive$flux.term, c(ft, ft))
  expect_equal(sp$ebullitive$flux_ebullitive[2], 0)
  expect_equal(sp$ebullitive$flux_ebullitive[1], sp$ebullitive$total_ebullitive[1] / 599 * ft)
  best <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(sp$diffusive, "CH4dry_ppb"))))
  expect_equal(best$LM.slope, c(0.5, 0.5), tolerance = 0.03)
  m <- merge_ebullition(best, sp$ebullitive)
  expect_equal(m$flux_total, m$best.flux + m$flux_ebullitive)
  expect_equal(m$n_events, c(1L, 0L))
  expect_gt(m$ebullitive_fraction[1], 0.5)
  expect_error(split_ebullition(list()), "detect_ebullition")
  # missing chamber columns -> NA with warning
  expect_warning(sp2 <- split_ebullition(detect_ebullition(d[, setdiff(names(d), "Vtot")], "CH4dry_ppb")), "chamber columns")
  expect_true(all(is.na(sp2$ebullitive$flux_ebullitive)))
})
