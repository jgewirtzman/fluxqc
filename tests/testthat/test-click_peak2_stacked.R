d <- example_closure()

test_that("resolve_gases handles auto, explicit and missing columns", {
  expect_equal(stemflux:::resolve_gases("auto", d), c("CO2dry_ppm", "CH4dry_ppb", "H2O_ppm"))
  expect_equal(stemflux:::resolve_gases(c("CH4dry_ppb", "CO2dry_ppm"), d), c("CH4dry_ppb", "CO2dry_ppm"))
  expect_error(stemflux:::resolve_gases("N2Odry_ppb", d), "not found")
  expect_error(stemflux:::resolve_gases("auto", data.frame(x = 1)), "no recognised")
})

test_that("stacked drawing works for 1, 2 and 3 panels on a non-interactive device", {
  pdf(NULL); on.exit(dev.off())
  for (g in list("CO2dry_ppm", c("CO2dry_ppm", "CH4dry_ppb"), c("CH4dry_ppb", "CO2dry_ppm", "H2O_ppm"))) {
    ys <- lapply(g, function(x) d[[x]])
    scr <- stemflux:::draw_stacked(d$POSIX.time, ys, g, stemflux:::stacked_ylims(ys, g, NULL), "title",
                                   v_blue = d$start.time[1], v_guide = d$start.time[1] + 600, time.axis = TRUE)
    expect_length(scr, length(g))
    close.screen(all.screens = TRUE)
  }
  ys <- list(d$CO2dry_ppm, d$CH4dry_ppb)
  scr <- stemflux:::draw_stacked(d$Etime, ys, c("CO2", "CH4"), stemflux:::stacked_ylims(ys, c("a", "b"), c(0, 1e5)),
                                 "t", v_red = c(0, 540), col = d$flag + 1)
  close.screen(all.screens = TRUE)
  expect_length(scr, 2)
})

test_that("plot.lim forms are honoured", {
  ys <- list(c(100, 200), c(1, 2))
  expect_equal(stemflux:::stacked_ylims(ys, c("a", "b"), NULL)[[1]], c(95, 205))
  expect_equal(stemflux:::stacked_ylims(ys, c("a", "b"), c(150, 5000))[[1]][1], 150)
  l <- stemflux:::stacked_ylims(ys, c("a", "b"), list(b = c(1.5, 10)))
  expect_equal(l[[1]], c(95, 205)); expect_equal(l[[2]][1], 1.5)
})

test_that("click_peak2_stacked with an empty seq returns an empty data.frame", {
  out <- click_peak2_stacked(list(d), gases = "auto", seq = integer(0), sleep = 0)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})
