t0 <- as.POSIXct("2023-03-15 10:00:00", tz = "UTC")

test_that("find_clock_offset recovers a known offset and corrects an auxfile", {
  set.seed(7)
  tt <- t0 + 0:3599
  co2 <- 410 + rnorm(3600, 0, 0.5)
  true_starts <- t0 + c(300, 1200, 2100, 3000)
  for (s in as.numeric(difftime(true_starts, t0, units = "secs")))
    co2[(s + 1):(s + 300)] <- co2[(s + 1):(s + 300)] + 0.3 * (1:300)
  traces <- data.frame(POSIX.time = tt, CO2dry_ppm = co2)
  notebook <- true_starts - 137
  pdf(NULL); on.exit(dev.off())
  aux <- data.frame(UniqueID = paste0("c", 1:4), start.time = notebook)
  fc <- find_clock_offset(traces, aux, search = c(-600, 600), step = 1, auxfile = aux, plot = TRUE)
  expect_lte(abs(fc$offset - 137), 3)
  expect_equal(nrow(fc$scores), 1201)
  expect_equal(fc$n_closures, 4)
  expect_equal(fc$auxfile$start.time, notebook + fc$offset)
  expect_equal(fc$auxfile$start.time_notebook, notebook)
  fc2 <- find_clock_offset(traces, notebook, search = c(-300, 300), step = 5, plot = FALSE)
  expect_lte(abs(fc2$offset - 137), 5)
  expect_error(find_clock_offset(traces, notebook, search = c(20000, 30000), plot = FALSE), "never fall")
  expect_error(find_clock_offset(traces, as.character(notebook), plot = FALSE), "POSIXct")
})
