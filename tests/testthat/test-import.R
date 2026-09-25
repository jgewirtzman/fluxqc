li7810 <- system.file("extdata", "LI7810", "LI7810.data", package = "goFlux")

write_fake_li7810 <- function() {
  f <- tempfile(fileext = ".data")
  hdr <- c("Model:\tLI-7810", "SN:\tTG10-00000", "Software Version:\t2.5.29",
           "Timestamp:\t2026-05-04 19:22:49", "Timezone:\tUTC",
           paste(c("DATAH", "SECONDS", "NANOSECONDS", "NDX", "DIAG", "REMARK", "DATE", "TIME",
                   "H2O", "CO2", "CH4"), collapse = "\t"),
           paste(c("DATAU", "secs", "nsecs", "index", "diag", "", "date", "time", "ppm", "ppm", "ppb"),
                 collapse = "\t"))
  row <- function(date, is) vapply(is, function(i)
    paste(c("DATA", format(1e9 + i, scientific = FALSE), 0, i, 0, '""', date,
            sprintf("10:00:%02d", i), 8000, 420 + i, 2000 + i), collapse = "\t"), character(1))
  writeLines(c(hdr, row("2026-05-01", 0:9), row("2026-05-02", 0:9), row("2026-05-03", 0:9)), f)
  f
}

test_that("import_li7810_subset keeps header and selected dates (awk and R paths)", {
  f <- write_fake_li7810()
  for (awk in c(TRUE, FALSE)) {
    out <- tempfile(fileext = ".data")
    expect_message(p <- import_li7810_subset(f, dates = c("2026-05-02", "2026-05-03"), outfile = out,
                                             import = FALSE, use.awk = awk), "20 data rows")
    l <- readLines(p)
    expect_length(l, 27)
    expect_match(l[1], "^Model")
    expect_false(any(grepl("2026-05-01", l)))
  }
  expect_warning(import_li7810_subset(f, dates = "1999-01-01", import = FALSE), "no DATA rows")
  expect_error(import_li7810_subset("nope.data", "2026-05-01"), "not found")
})

test_that("import_li7810_subset imports through goFlux", {
  skip_if(!nzchar(li7810), "goFlux example LI7810.data not installed")
  imp <- import_li7810_subset(li7810, dates = "2022-12-05", timezone = "UTC")
  expect_true(all(c("POSIX.time", "CO2dry_ppm", "CH4dry_ppb") %in% names(imp)))
  expect_gt(nrow(imp), 100)
  expect_true(file.exists(attr(imp, "subset.file")))
})

test_that("import2RData_flat imports a folder without writing RData", {
  skip_if(!nzchar(li7810), "goFlux example LI7810.data not installed")
  wd <- tempfile("wd_"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old))
  expect_message(imp <- import2RData_flat(dirname(li7810), "LI7810", pattern = "\\.data$",
                                          timezone = "UTC"), "imported")
  expect_false(dir.exists(file.path(wd, "RData")))
  expect_equal(unique(imp$source_file), "LI7810.data")
  expect_true(!is.unsorted(imp$POSIX.time))
  lst <- import2RData_flat(dirname(li7810), "LI7810", pattern = "\\.data$", timezone = "UTC", merge = FALSE)
  expect_type(lst, "list")
  expect_error(import2RData_flat(dirname(li7810), "NOPE"), "no importer")
  expect_error(import2RData_flat(tempdir(), "LI7810", pattern = "zzz$"), "no files")
})
