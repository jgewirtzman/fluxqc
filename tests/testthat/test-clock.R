t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")

write_reset_file <- function() {
  f <- tempfile(fileext = ".txt")
  hdr <- c("SN: 123", "SysTime, [CH4]d_ppm, [CO2]d_ppm")
  good <- sprintf("%s, 2.0, 400", format(t0 + 0:19, "%m/%d/%Y %H:%M:%S"))
  reset <- sprintf("%s, 2.1, 401", format(t0 - 600 + 0:9, "%m/%d/%Y %H:%M:%S"))
  # interspersed row by row
  writeLines(c(hdr, as.vector(rbind(good[1:10], reset)), good[11:20]), f)
  f
}

test_that("detect_clock_streams labels the reset stream", {
  f <- write_reset_file()
  cs <- detect_clock_streams(f)
  expect_true(cs$has_reset)
  expect_equal(as.vector(table(cs$streams)[c("ahead", "behind")]), c(20L, 10L))
  expect_length(cs$header_lines, 2)
  expect_equal(cs$stream_labels, c("ahead", "behind"))
  # no reset
  f2 <- tempfile(fileext = ".txt")
  writeLines(c("SysTime, x", sprintf("%s, 1", format(t0 + 0:9, "%m/%d/%Y %H:%M:%S"))), f2)
  expect_false(detect_clock_streams(f2)$has_reset)
  out <- write_stream_files(cs, "micro_f0001.txt", tempfile())
  expect_equal(names(out), c("ahead", "behind"))
  expect_equal(length(readLines(out["behind"])), 12)
})

test_that("separate_window_streams splits interleaved streams", {
  a <- 400 + 0.2 * (0:199); b <- 900 + 0.2 * (0:199)
  d <- data.frame(POSIX.time = t0 + rep(0:199, each = 2), CO2dry_ppm = as.vector(rbind(a, b)))
  s <- separate_window_streams(d, threshold = 100)
  expect_equal(as.vector(table(s$stream_id)), c(200L, 200L))
  # with every row alternating, the default threshold (3 x median jump) never switches
  expect_equal(unique(separate_window_streams(d)$stream_id), "stream_a")
  expect_equal(s$CO2dry_ppm[s$stream_id == "stream_a"], a)
  expect_true(all(is.finite(attr(s, "roughness"))))
  tiny <- separate_window_streams(d[1:3, ])
  expect_equal(unique(tiny$stream_id), "stream_a")
})

test_that("obs_win_deinterleave keeps the smoother stream and adds chamber columns", {
  a <- 400 + 0.2 * (0:599); b <- 900 + 5 * sin(0:599)
  raw <- data.frame(POSIX.time = t0 + rep(0:599, each = 2), CO2dry_ppm = as.vector(rbind(a, b)),
                    source_file = "x.txt")
  aux <- data.frame(UniqueID = "w1", start.time = t0 + 120, Area = 324, Vtot = 4.3,
                    Tcham = 25, Pcham = 101.3)
  # fully row-by-row interleaving: the default 3 x median threshold is too large, so give one
  ow <- obs_win_deinterleave(raw, aux, obs_length = 300, shoulder = 60, threshold = 100)
  expect_named(ow, "w1")
  w <- ow$w1
  expect_equal(nrow(w), 421)
  expect_true(all(w$CO2dry_ppm < 600))
  expect_false("source_file" %in% names(w))
  expect_true(all(c("Etime", "Area", "Vtot", "Tcham", "Pcham", "start.time", "obs.length") %in% names(w)))
  expect_equal(w$Etime[1], -60)
  ob <- obs_win_deinterleave(raw, aux, obs_length = 300, shoulder = 60, keep = "stream_b", threshold = 100)
  expect_true(all(ob$w1$CO2dry_ppm > 800))
  expect_message(obs_win_deinterleave(raw, data.frame(UniqueID = "z", start.time = t0 + 5000),
                                      obs_length = 300, shoulder = 60), "skipped")
})
