d <- example_closure()
best <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))))
co2  <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(d, "CO2dry_ppm"))))

test_that("qc_screens runs every screen on the example closure and never drops rows", {
  out <- qc_screens(best, traces = d, gastype = "CH4dry_ppb", co2 = co2,
                    screens = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE,
                                   min_window = list(secs = 120), ambient_start = TRUE,
                                   leak = list(rate = 1e-6)), verbose = FALSE)
  expect_equal(nrow(out), 1)
  expect_true(all(c("qc_c0", "qc_co2_tracer", "qc_convex", "qc_min_window", "qc_ambient_start",
                    "qc_leak", "qc_any", "qc_note") %in% names(out)))
  expect_false(out$qc_c0); expect_false(out$qc_co2_tracer); expect_false(out$qc_min_window)
  expect_equal(out$qc_co2_slope, co2$LM.slope)
  expect_true(is.logical(out$qc_convex))
  expect_false(is.na(out$qc_ambient_start))     # the obs.win context supplies pre-closure rows
  expect_false(out$qc_leak)                     # 1e-6 x range x flux.term << goFlux MDF
  expect_equal(names(attr(out, "qc_screens")), c("c0", "co2_tracer", "convex", "min_window", "ambient_start", "leak"))
})

test_that("individual screens fire as designed", {
  # c0: five closures, one contaminated; grouped
  fx <- data.frame(UniqueID = paste0("m", 1:6), C0 = c(2100, 2050, 1980, 2300, 11700, NA),
                   camp = c("a", "a", "a", "a", "a", "b"), nb.obs = c(300, 300, 300, 30, 300, 300),
                   HM.k = c(0.01, 0.01, -0.02, 0.01, 0.01, NA), LM.slope = c(0.5, -0.1, 0.5, 0.5, 0.5, 0.5))
  out <- qc_screens(fx, screens = list(c0 = list(k = 1.5), min_window = list(secs = 60), convex = TRUE),
                    group = "camp", verbose = FALSE)
  expect_equal(out$qc_c0, c(FALSE, FALSE, FALSE, FALSE, TRUE, NA))
  expect_equal(out$qc_min_window, c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(out$qc_convex, c(FALSE, FALSE, TRUE, FALSE, FALSE, NA))   # HM.k < 0 without traces
  expect_equal(out$qc_any, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_match(out$qc_note[5], "C0 >")
  expect_match(out$qc_note[3], "convex")
  expect_true(is.na(out$qc_note[1]))
  # co2_tracer with tissue exemption, fluxes is the CO2 table itself
  out2 <- qc_screens(fx, gastype = "CO2dry_ppm", tissue = c("live", "live", "live", "dead", "live", "live"),
                     screens = list(co2_tracer = list(min_slope = 0)), verbose = FALSE)
  expect_equal(out2$qc_co2_tracer, c(FALSE, TRUE, FALSE, NA, FALSE, FALSE))   # row 4 is dead tissue
  out3 <- qc_screens(fx, gastype = "CO2dry_ppm", tissue = c("live", "dead", "live", "live", "live", "live"),
                     screens = list(co2_tracer = TRUE), verbose = FALSE)
  expect_true(is.na(out3$qc_co2_tracer[2]))
  expect_error(qc_screens(fx, screens = list(bogus = TRUE)), "unknown screen")
})

test_that("convex, ambient_start and leak screens use traces", {
  t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")
  set.seed(3)
  mk <- function(uid, y) data.frame(UniqueID = uid, POSIX.time = t0 + (-60):299, CH4dry_ppb = y,
                                    flag = c(rep(0, 60), rep(1, 300)), Etime = (-60):299, obs.length_corr = 299,
                                    start.time = t0)
  lin  <- mk("lin",  c(rep(2000, 60), 2000 + 0.5 * (0:299)) + rnorm(360, 0, 1))
  conv <- mk("conv", c(rep(2000, 60), 2000 + 0.002 * (0:299)^2) + rnorm(360, 0, 1))
  late <- mk("late", c(rep(2000, 60), 2100 + 0.5 * (0:299)) + rnorm(360, 0, 1))   # already enriched at start
  tr <- rbind(lin, conv, late)
  fx <- data.frame(UniqueID = c("lin", "conv", "late"), flux.term = 0.68, MDF = 999, MDF_emp = 0.02,
                   C0 = 2000, Ct = 2150)   # leak compares with MDF_emp, not goFlux's MDF
  out <- qc_screens(fx, traces = tr, gastype = "CH4dry_ppb",
                    screens = list(convex = TRUE, ambient_start = list(sigma = 3, secs = 10), leak = list(rate = 1e-4)),
                    verbose = FALSE)
  expect_equal(out$qc_convex, c(FALSE, TRUE, FALSE))
  expect_equal(out$qc_ambient_start, c(FALSE, FALSE, TRUE))
  # leak: rate * range * flux.term vs MDF
  rng <- diff(range(lin$CH4dry_ppb[lin$flag == 1]))
  expect_equal(out$qc_leak_apparent_flux[1], 1e-4 * rng * 0.68)
  expect_equal(out$qc_leak, out$qc_leak_apparent_flux >= 0.02)
})

test_that("qc_c0_screen wrapper keeps the 0.1.0 columns", {
  res <- data.frame(UniqueID = paste0("m", 1:6), C0 = c(2100, 2050, 1980, 2300, 11700, NA))
  out <- qc_c0_screen(res, verbose = FALSE)
  expect_equal(out$qc_c0_pass, c(TRUE, TRUE, TRUE, TRUE, FALSE, NA))
  expect_match(out$qc_c0_reason[5], "not at ambient")
  expect_equal(out$qc_c0_ratio[5], 11700 / median(c(2100, 2050, 1980, 2300, 11700)))
  expect_error(qc_c0_screen(res, C0.col = "nope"), "not found")
})

test_that("noisy screen: per-closure sigma vs the campaign sigma behind MDF_emp", {
  d2 <- d; d2$UniqueID <- "noisy_one"; set.seed(9); d2$CH4dry_ppb <- d2$CH4dry_ppb + rnorm(nrow(d2), 0, 6)
  both <- rbind(d, d2)
  b <- suppressWarnings(goFlux::best.flux(suppressWarnings(goFlux::goFlux(both, "CH4dry_ppb"))))
  b <- flag_detection(b, traces = both, gastype = "CH4dry_ppb")           # sigma_emp = campaign MAD (one group)
  out <- qc_screens(b, traces = both, gastype = "CH4dry_ppb", screens = list(noisy = list(k = 1.5)), verbose = FALSE)
  expect_equal(out$qc_noisy_sigma, c(precision_mad(d$CH4dry_ppb[d$flag == 1]), precision_mad(d2$CH4dry_ppb[d2$flag == 1])))
  expect_equal(out$qc_noisy_ratio, out$qc_noisy_sigma / out$sigma_emp)
  expect_equal(out$qc_noisy, c(FALSE, TRUE))
  expect_match(out$qc_note[2], "closure sigma = .* ppb = .*x campaign sigma \\(disturbed closure\\?\\)")
  expect_true(is.na(out$qc_note[1]))
  # without sigma_emp the campaign sigma is recomputed per group from the whole record
  b2 <- b[, setdiff(names(b), "sigma_emp")]; b2$camp <- c("x", "y")
  o2 <- qc_screens(b2, traces = both, gastype = "CH4dry_ppb", group = "camp", screens = list(noisy = TRUE), verbose = FALSE)
  expect_equal(o2$qc_noisy_ratio[2], precision_mad(d2$CH4dry_ppb[d2$flag == 1]) / precision_mad(d2$CH4dry_ppb))
  expect_false(o2$qc_noisy[2])                                             # its own group: not noisier than itself
  # no traces -> NA
  expect_true(is.na(qc_screens(b, screens = list(noisy = TRUE), verbose = FALSE)$qc_noisy[1]))
})
