test_that("qc_c0_screen flags high starting concentrations", {
  res <- data.frame(UniqueID = paste0("m", 1:6), C0 = c(2100, 2050, 1980, 2300, 11700, NA),
                    camp = c("a", "a", "a", "a", "a", "b"))
  expect_message(out <- qc_c0_screen(res), "rejected = 1")
  expect_equal(out$qc_c0_pass, c(TRUE, TRUE, TRUE, TRUE, FALSE, NA))
  expect_match(out$qc_c0_reason[5], "not at ambient")
  expect_true(is.na(out$qc_c0_reason[1]))
  expect_equal(out$qc_c0_ratio[5], 11700 / median(c(2100, 2050, 1980, 2300, 11700)))
  thr <- attr(out, "qc_c0_threshold")
  expect_equal(thr$group, "all")
  # grouped medians
  out2 <- qc_c0_screen(res, group = "camp", verbose = FALSE)
  expect_equal(nrow(attr(out2, "qc_c0_threshold")), 2)
  expect_error(qc_c0_screen(res, C0.col = "nope"), "not found")
})
