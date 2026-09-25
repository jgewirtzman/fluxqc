#' Bundled example closures with real chamber geometry
#'
#' Six chamber closures shipped in `inst/extdata/`, returned in goFlux shape
#' (the structure of [goFlux::click.peak2()] output, with `flag`, `Etime`,
#' `start.time_corr`, ... already set) together with an auxfile that carries
#' the real geometry of each closure. They are the teaching cases of the
#' companion guidelines paper (Figure 4 logic): the detection flags and QC
#' screens of [process_fluxes()] separate them.
#'
#' | id | UniqueID | what it shows | kind | source / geometry |
#' |---|---|---|---|---|
#' | `a` | `2021-05-18_5086_I_7.0` | high-signal valid closure: CO2 and CH4 both rise | real | Yale-Myers Forest (CT, USA) tree stem, semi-rigid sleeve chamber, ABB/LGR UGGA, 1 Hz, 600 s; geometry from the campaign auxfile (`auxfile_goFlux_with_weather.csv`) |
#' | `b` | `2020-09-23_5092_I_6` | valid closure, CH4 below detection (healthy CO2, flat CH4) | real | same campaign |
#' | `c` | `2021-02-25_1641_WD_5` | failed closure: flat CO2 and CH4 (chamber not sealed, winter) | real | same campaign |
#' | `d` | `2020-09-21_3589_U_6` | ambiguous: CH4 rise over a flat CO2 tracer | real | same campaign |
#' | `e` | `BF_CP40_2022-10-23_P06` | aquatic floating chamber with ebullition (5 bubble events) | real | NASA BlueFlux campaign, Everglades (FL, USA), site CP40, ABB/LGR GLA131 (LGR2), 276 s; Area 324.3 cm2, Vtot 4.318 L (chamber + tubing + cell) from the campaign's chamber parameters, Tcham = mean air temperature of that site-day's water measurements (28.7 C); the trace is as segmented by the campaign's ebullition script (clock offset applied, 20 s trimmed at each end) |
#' | `f` | `SYN_uptake_from_a` | **synthetic** CH4 uptake (-0.08 ppb s-1) with positive CO2 (+0.6 ppm s-1), Gaussian noise at the measured Allan precision of closure `a` (2.2 ppb CH4, 0.8 ppm CO2), seed 42 | synthetic | geometry copied from `a` |
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
#' The archived traces hold dry mole fractions only, so `H2O_ppm` is 0 (no
#' water-vapour correction) except for `e`, where the record has none either
#' and 0 is used as well. Manufacturer precision written to the `*_prec`
#' columns is the UGGA/GLA131 datasheet value used by the source campaigns
#' (`c(0.35, 0.9, 200)` for CO2 ppm, CH4 ppb, H2O ppm).
#'
#' @param which character vector of example ids (`"a"` .. `"f"`) or
#'   `UniqueID`s; default all six.
#' @param dead.band seconds discarded at the start of every closure
#'   (`flag = 0`); for `e` the 20 s trim was already applied by the source
#'   script, so the default dead band is reduced to 0 for that closure unless
#'   `dead.band` is given explicitly as a named vector.
#' @param prec manufacturer precision `c(CO2 ppm, CH4 ppb, H2O ppm)`.
#' @return list with `data` (one goFlux-shaped data.frame, all closures bound
#'   together) and `aux` (data.frame: `UniqueID`, `start.time`, `obs.length`,
#'   `Area`, `Vtot`, `Tcham`, `Pcham`, `example`, `kind`, `what`, `site`,
#'   `geometry_source`).
#' @examples
#' ex <- stemflux_examples()
#' ex$aux[, c("example", "UniqueID", "kind", "what")]
#' table(ex$data$UniqueID)
#' one <- stemflux_examples("c")$data
#' plot(one$Etime, one$CO2dry_ppm, pch = 16, cex = 0.4)
#' @seealso [example_closure()]
#' @export
stemflux_examples <- function(which = c("a", "b", "c", "d", "e", "f"), dead.band = 60,
                              prec = c(0.35, 0.9, 200)) {
  f <- system.file("extdata", "stemflux_examples.csv", package = "stemflux")
  fa <- system.file("extdata", "stemflux_examples_aux.csv", package = "stemflux")
  if (f == "" || fa == "") stop("example files not found in the installed package")
  d <- read.csv(f, stringsAsFactors = FALSE)
  ax <- read.csv(fa, stringsAsFactors = FALSE)
  sel <- ax$example %in% which | ax$UniqueID %in% which
  if (!any(sel)) stop("no example matches: ", paste(which, collapse = ", "))
  ax <- ax[sel, ]
  ax$start.time <- as.POSIXct(ax$start.time, tz = "UTC")
  d <- d[d$example %in% ax$example, ]
  d$POSIX.time <- as.POSIXct(d$POSIX.time, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  out <- list()
  for (i in seq_len(nrow(ax))) {
    a <- ax[i, ]
    x <- d[d$example == a$example, setdiff(names(d), "example")]
    db <- if (!is.null(names(dead.band)) && a$example %in% names(dead.band)) dead.band[[a$example]]
          else if (a$example == "e" && is.null(names(dead.band))) 0 else dead.band[[1]]
    x$H2O_ppm <- 0
    x$start.time <- a$start.time; x$obs.length <- a$obs.length
    x$Area <- a$Area; x$Vtot <- a$Vtot; x$Tcham <- a$Tcham; x$Pcham <- a$Pcham
    x$CO2_prec <- prec[1]; x$CH4_prec <- prec[2]; x$H2O_prec <- prec[3]
    out[[i]] <- flag_window(x, a$start.time + db, max(x$POSIX.time) + 0.01)
  }
  data <- do.call(rbind, out); rownames(data) <- NULL
  aux <- ax[, c("UniqueID", "start.time", "obs.length", "Area", "Vtot", "Tcham", "Pcham",
                "example", "kind", "what", "site", "geometry_source")]
  rownames(aux) <- NULL
  list(data = data, aux = aux)
}

#' One real chamber closure, ready for goFlux
#'
#' Alias for `stemflux_examples("a")$data`: closure `2021-05-18_5086_I_7.0`,
#' a clean high-signal 600 s tree-stem closure (Yale-Myers Forest, ABB/LGR
#' UGGA). See [stemflux_examples()] for the geometry and the other bundled
#' closures.
#'
#' @inheritParams stemflux_examples
#' @return goFlux-shaped data.frame for the one closure.
#' @examples
#' d <- example_closure()
#' nrow(d); sum(d$flag)
#' @export
example_closure <- function(dead.band = 60, prec = c(0.35, 0.9, 200)) {
  stemflux_examples("a", dead.band = dead.band, prec = prec)$data
}
