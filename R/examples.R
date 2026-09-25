#' Bundled example closures with real chamber geometry
#'
#' Six chamber closures shipped in `inst/extdata/`, returned in goFlux shape
#' (the structure of [goFlux::click.peak2()] output) together with an auxfile
#' that carries the real geometry of each closure. They are the teaching cases
#' of the companion guidelines paper (Figure 4 logic): the detection flags and
#' QC screens of [process_fluxes()] separate them.
#'
#' Closures `a`-`d` are the campaign's own *clicked* goFlux windows
#' (`semirigid_tree_manID_clicked.rds`, the manID output of the Yale-Myers
#' processing): each is shipped as the full `obs.win()` segment, i.e. the
#' surrounding context rows with `flag = 0` and the operator-clicked fitting
#' window with `flag = 1`, so the dead band actually applied is visible.
#' Geometry comes from the same rows. They were chosen by auditable criteria
#' against the campaign's final flux table (`flux_FINAL.csv`) and a linear fit
#' of CO2 on the flagged rows:
#'
#' | id | UniqueID | window (s after `start.time`, length) | CO2 slope (ppm/s), p | CH4 flux, MDF (nmol m-2 s-1), class in flux_FINAL | what it shows |
#' |---|---|---|---|---|---|
#' | `a` | `2021-05-18_5086_I_7.0` | +91 s, 765 s, n = 770 | 0.824, p < 1e-300 | 0.254, 0.00313, emission | high-signal valid closure: CO2 and CH4 both rise |
#' | `b` | `2020-09-23_5092_I_6` | +49 s, 668 s, n = 673 | 0.717, p < 1e-300 | -0.00118, 0.00363, below detection | valid closure, CH4 below detection |
#' | `c` | `2021-02-25_1885_I_5_1` | +64 s, 348 s, n = 351 | -0.0253, p = 1.4e-13 | 0.00181, 0.00752, below detection | failed closure: CO2 not rising (the only stem closure in the archive with CO2 not clearly positive and CH4 below detection) |
#' | `d` | `2020-12-09_5091_I_5` | +103 s, 583 s, n = 587 | -0.00029, p = 0.69 | -0.00463, 0.00461, uptake | ambiguous: CH4 detected (a marginal uptake, just past the MDF) over a flat CO2 tracer (the only closure with CH4 detected and CO2 flat) |
#'
#' The full diagnostics table is shipped as
#' `system.file("extdata", "stemflux_examples_diagnostics.csv", package = "stemflux")`.
#'
#' | id | UniqueID | what it shows | kind | source / geometry |
#' |---|---|---|---|---|
#' | `e` | `BF_CP40_2022-10-23_P06` | aquatic floating chamber with one bubble burst (5 consecutive 1 Hz samples at 42-46 s, ~675 ppb) inside a sealed closure | real | NASA BlueFlux campaign, Everglades (FL, USA), site CP40, ABB/LGR GLA131 (LGR2), 276 s; Area 324.3 cm2, Vtot 4.318 L from the campaign's chamber parameters, Tcham = mean air temperature of that site-day's water measurements (28.7 C); the trace is as segmented by the campaign's ebullition script (clock offset applied, 20 s trimmed at each end) |
#' | `f` | `SYN_uptake_from_a` | **synthetic** CH4 uptake (-0.08 ppb s-1) with positive CO2 (+0.6 ppm s-1), Gaussian noise at the measured Allan precision of closure `a` (2.2 ppb CH4, 0.8 ppm CO2), seed 42 | synthetic | geometry copied from `a`; 60 s dead band |
#'
#' @section Data release and citation:
#' The real closures are released with the package under its MIT license,
#' with the permission of the data owner. `a`-`d` come from the Yale-Myers
#' Forest tree-stem methane record (2020-2021; J. Gewirtzman and colleagues,
#' Yale School of the Environment); cite them as "Gewirtzman et al.,
#' Yale-Myers Forest stem-flux record 2020-2021, distributed with the
#' stemflux R package". `e` comes from the ground component of the NASA
#' Carbon Monitoring System BlueFlux field campaign (site CP40, Christian
#' Point, a ghost-forest mangrove site in south Florida; ABB/LGR GLA131); cite
#' Poulter, B., et al. (2023), Multi-scale observations of mangrove blue
#' carbon ecosystem fluxes: The NASA Carbon Monitoring System BlueFlux field
#' campaign, Environmental Research Letters 18, 075009,
#' \doi{10.1088/1748-9326/acdae6}, and the blueflux-ground repository (MIT
#' license). `f` is synthetic and carries no data rights.
#'
#' `a`-`d` carry the campaign's water column (`H2O_ppm`); `e` and `f` have
#' none and use 0. The `*_prec` columns hold the UGGA/GLA131 datasheet value
#' used by the source campaigns (`c(0.35, 0.9, 200)` for CO2 ppm, CH4 ppb,
#' H2O ppm).
#'
#' @param which character vector of example ids (`"a"` .. `"f"`) or
#'   `UniqueID`s; default all six.
#' @param dead.band `NULL` (default) keeps the shipped `flag` column (the
#'   clicked window for `a`-`d`, the whole trimmed trace for `e`, 60 s for
#'   `f`). A number re-flags every closure as `start.time + dead.band` to the
#'   end of its segment, ignoring the clicked windows.
#' @param prec manufacturer precision `c(CO2 ppm, CH4 ppb, H2O ppm)`.
#' @return list with `data` (one goFlux-shaped data.frame, all closures bound
#'   together, with `flag`, `Etime`, `start.time_corr`, `end.time_corr`,
#'   `obs.length_corr`) and `aux` (data.frame: `UniqueID`, `start.time`,
#'   `obs.length`, `Area`, `Vtot`, `Tcham`, `Pcham`, `start.time_corr`,
#'   `end.time_corr`, `example`, `kind`, `what`, `site`, `geometry_source`).
#' @examples
#' ex <- stemflux_examples()
#' ex$aux[, c("example", "UniqueID", "kind", "what")]
#' table(ex$data$UniqueID, ex$data$flag)
#' one <- stemflux_examples("c")$data
#' plot(one$Etime, one$CO2dry_ppm, pch = 16, cex = 0.4, col = one$flag + 1)
#' @seealso [example_closure()]
#' @export
stemflux_examples <- function(which = c("a", "b", "c", "d", "e", "f"), dead.band = NULL,
                              prec = c(0.35, 0.9, 200)) {
  f <- system.file("extdata", "stemflux_examples.csv", package = "stemflux")
  fa <- system.file("extdata", "stemflux_examples_aux.csv", package = "stemflux")
  if (f == "" || fa == "") stop("example files not found in the installed package")
  d <- read.csv(f, stringsAsFactors = FALSE)
  ax <- read.csv(fa, stringsAsFactors = FALSE)
  sel <- ax$example %in% which | ax$UniqueID %in% which
  if (!any(sel)) stop("no example matches: ", paste(which, collapse = ", "))
  ax <- ax[sel, ]
  for (cc in c("start.time", "start.time_corr", "end.time_corr")) ax[[cc]] <- as.POSIXct(ax[[cc]], tz = "UTC")
  d <- d[d$example %in% ax$example, ]
  d$POSIX.time <- as.POSIXct(d$POSIX.time, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  out <- list()
  for (i in seq_len(nrow(ax))) {
    a <- ax[i, ]
    x <- d[d$example == a$example, setdiff(names(d), "example")]
    x$start.time <- a$start.time; x$obs.length <- a$obs.length
    x$Area <- a$Area; x$Vtot <- a$Vtot; x$Tcham <- a$Tcham; x$Pcham <- a$Pcham
    x$CO2_prec <- prec[1]; x$CH4_prec <- prec[2]; x$H2O_prec <- prec[3]
    if (is.null(dead.band)) {
      fl <- x$flag
      x <- flag_window(x, a$start.time_corr, a$end.time_corr)
      x$flag <- fl                                   # keep the shipped (clicked) flags exactly
    } else {
      x <- flag_window(x, a$start.time + dead.band, max(x$POSIX.time) + 0.01)
    }
    out[[i]] <- x
  }
  data <- do.call(rbind, out); rownames(data) <- NULL
  aux <- ax[, c("UniqueID", "start.time", "obs.length", "Area", "Vtot", "Tcham", "Pcham",
                "start.time_corr", "end.time_corr", "example", "kind", "what", "site", "geometry_source")]
  rownames(aux) <- NULL
  list(data = data, aux = aux)
}

#' One real chamber closure, ready for goFlux
#'
#' Alias for `stemflux_examples("a")$data`: closure `2021-05-18_5086_I_7.0`,
#' a clean high-signal tree-stem closure (Yale-Myers Forest, ABB/LGR UGGA)
#' with the campaign's clicked window (765 s starting 91 s after the
#' recorded `start.time`, 770 flagged rows) and its obs.win context. See
#' [stemflux_examples()] for the geometry and the other bundled closures.
#'
#' @inheritParams stemflux_examples
#' @return goFlux-shaped data.frame for the one closure.
#' @examples
#' d <- example_closure()
#' nrow(d); sum(d$flag)
#' @export
example_closure <- function(dead.band = NULL, prec = c(0.35, 0.9, 200)) {
  stemflux_examples("a", dead.band = dead.band, prec = prec)$data
}
