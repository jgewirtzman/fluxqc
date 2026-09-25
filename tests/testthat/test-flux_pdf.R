test_that("flux_pdf writes a PDF grouped by UniqueID", {
  skip_if_not_installed("ggplot2")
  make_plot <- function(uid, gas) {
    UniqueID <- uid
    df <- data.frame(x = 1:10, y = rnorm(10))
    ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point() + ggplot2::ggtitle(gas)
  }
  pl <- list(make_plot("A", "CO2"), make_plot("B", "CO2"), make_plot("A", "CH4"), make_plot("B", "CH4"))
  expect_equal(vapply(seq_along(pl), function(i) fluxqc:::plot_uid(pl[[i]], i), character(1)),
               c("A", "B", "A", "B"))
  out <- tempfile(fileext = ".pdf")
  expect_invisible(flux_pdf(pl, out))
  expect_true(file.exists(out))
  expect_gt(file.size(out), 1000)
  out2 <- tempfile(fileext = ".pdf")
  flux_pdf(pl, out2, group.by = "none")
  expect_true(file.exists(out2))
  expect_error(flux_pdf(list(1, 2), out), "ggplot")
  expect_error(flux_pdf(pl, outfile = 3), "character")
})

test_that("flux_pdf accepts real goFlux flux.plot output", {
  skip_if_not_installed("ggplot2")
  d <- example_closure()
  res <- suppressWarnings(goFlux::goFlux(d, "CH4dry_ppb"))
  best <- suppressWarnings(goFlux::best.flux(res))
  pl <- suppressWarnings(goFlux::flux.plot(best, d, "CH4dry_ppb"))
  expect_equal(fluxqc:::plot_uid(pl[[1]], 1), "emission_stem_semirigid")
  out <- tempfile(fileext = ".pdf")
  flux_pdf(pl, out)
  expect_true(file.exists(out))
})
