t0 <- as.POSIXct("2024-06-01 10:00:00", tz = "UTC")

make_window <- function(uid = "A", n = 900, start = 180) {
  set.seed(1)
  data.frame(UniqueID = uid, POSIX.time = t0 + 0:(n - 1),
             CO2dry_ppm = c(rep(410, 200), 410 + 0.3 * (0:400), rep(530, n - 601)) + rnorm(n, 0, 0.3),
             CH4dry_ppb = 2000 + 0.01 * (0:(n - 1)) + rnorm(n, 0, 0.5),
             start.time = t0 + start, obs.length = 600)
}

test_that("flag_window adds the click.peak2 columns", {
  d <- data.frame(POSIX.time = t0 + 0:299, CO2dry_ppm = 400 + 0.2 * (0:299))
  w <- flag_window(d, t0 + 30, t0 + 200)
  expect_equal(sum(w$flag), 171)
  expect_equal(w$Etime[w$flag == 1][1], 0)
  expect_equal(unique(w$obs.length_corr), 170)
  expect_equal(unique(w$start.time_corr), t0 + 30)
  expect_error(flag_window(d, t0 + 200, t0 + 30), "after")
})

test_that("windows_from_table builds a manID and warns on short windows", {
  ow <- list(make_window("A"), make_window("B"))
  win <- data.frame(UniqueID = c("A", "B"),
                    start = c("2024-06-01 10:01:00", "2024-06-01 10:02:00"),
                    end = c("2024-06-01 10:05:00", "2024-06-01 10:02:30"))
  expect_warning(m <- windows_from_table(ow, win), "only 31")
  expect_equal(sort(unique(m$UniqueID)), c("A", "B"))
  expect_equal(sum(m$flag[m$UniqueID == "A"]), 241)
  expect_true(all(c("flag", "Etime", "start.time_corr", "end.time_corr", "obs.length_corr") %in% names(m)))
  # missing UniqueID is dropped with a message
  expect_message(m2 <- windows_from_table(ow, win[1, ]), "No window for B")
  expect_equal(unique(m2$UniqueID), "A")
})

test_that("find_rise locates the sustained CO2 rise", {
  d <- make_window()
  r <- find_rise(d$POSIX.time, d$CO2dry_ppm)
  expect_false(is.null(r))
  expect_lt(abs(as.numeric(r$start - (t0 + 200), units = "secs")), 40)
  expect_gt(r$dur, 300)
  expect_gt(r$dconc, 100)
  expect_null(find_rise(d$POSIX.time, rep(410, nrow(d)) + rnorm(nrow(d), 0, 0.3)))
})

test_that("auto_id_rise returns click.peak2 structure and a log", {
  ow <- list(make_window("A"), make_window("B"))
  ow[[2]]$CO2dry_ppm <- 410 + rnorm(900, 0, 0.3)   # no rise -> fallback
  ends <- data.frame(UniqueID = "B", end.time = "2024-06-01 10:07:00")
  dir <- tempfile("autoid_");
  m <- auto_id_rise(ow, end.times = ends, plot.dir = dir)
  lg <- attr(m, "log")
  expect_equal(lg$method, c("auto_rise", "notebook_fallback"))
  expect_true(all(lg$n_flag > 60))
  expect_equal(length(list.files(dir, pattern = "png$")), 2)
  expect_true(all(c("flag", "Etime", "start.time_corr", "end.time_corr", "obs.length_corr") %in% names(m)))
  # fallback window = start + 30 .. notebook end
  expect_equal(unique(m$start.time_corr[m$UniqueID == "B"]), t0 + 180 + 30)
  expect_equal(unique(m$end.time_corr[m$UniqueID == "B"]), t0 + 420)
  # rise window is capped at max.secs
  expect_lte(unique(m$obs.length_corr[m$UniqueID == "A"]), 600)
})
