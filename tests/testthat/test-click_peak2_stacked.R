test_that("stacked drawing works on a non-interactive device", {
  d <- example_closure()
  pdf(NULL)
  on.exit(dev.off())
  scr <- stemflux:::draw_stacked(d$POSIX.time, d$CO2dry_ppm, d$CH4dry_ppb,
                                 stemflux:::stacked_ylim(d$CO2dry_ppm, c(300, 5000)),
                                 stemflux:::stacked_ylim(d$CH4dry_ppb, c(1500, 1e5)),
                                 "CO2dry_ppm", "CH4dry_ppb", "title",
                                 v_blue = d$start.time[1], v_guide = d$start.time[1] + 600,
                                 time.axis = TRUE)
  expect_length(scr, 2)
  close.screen(all.screens = TRUE)
  # validation-style plot on Etime with red window lines
  scr <- stemflux:::draw_stacked(d$Etime, d$CO2dry_ppm, d$CH4dry_ppb,
                                 c(300, 1500), c(2000, 2500), "CO2", "CH4", "t",
                                 v_red = c(0, 540), col = d$flag + 1)
  close.screen(all.screens = TRUE)
  expect_length(scr, 2)
})

test_that("stacked_ylim clips to the plot limits", {
  expect_equal(stemflux:::stacked_ylim(c(100, 200), c(150, 5000))[1], 150)
  expect_equal(stemflux:::stacked_ylim(c(100, 200), c(0, 190))[2], 190)
  yl <- stemflux:::stacked_ylim(c(400, 500), c(300, 5000))
  expect_equal(yl, c(395, 505))
})

test_that("click_peak2_stacked with an empty seq returns an empty data.frame", {
  d <- example_closure()
  out <- click_peak2_stacked(list(d), seq = integer(0), sleep = 0)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})
