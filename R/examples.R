#' Curated, anonymized example closures with real chamber geometry
#'
#' Ten real chamber closures from the author's campaigns, shipped in
#' `inst/extdata/` and returned in goFlux shape (the structure of
#' [goFlux::click.peak2()] output) together with an auxfile carrying the real
#' geometry of each. Every closure is the campaign's own goFlux manID output:
#' the full `obs.win()` segment (context rows `flag = 0`, operator-clicked
#' fitting window `flag = 1`; the ebullition closure has no clicked window and
#' is the campaign-segmented trace). They are the teaching cases of the
#' companion guidelines paper (Figure 4 logic).
#'
#' **Anonymization.** `UniqueID`s are descriptive keys; tree tags, plot codes,
#' species, operators, file names and dates are removed; the time axis is a
#' synthetic calendar (2000-01-01 + one day per closure) that preserves every
#' sampling interval; analyzer model, chamber family, geometry (`Area`,
#' `Vtot`, `Tcham`, `Pcham`), season and year are kept. The extraction is
#' documented, reproducibly from the author's drive, in
#' `system.file("extdata", "build_examples.R", package = "stemflux")`, and the
#' audit numbers below are in `stemflux_examples_diagnostics.csv` next to it.
#'
#' | key | shows | campaign type, season | analyzer, sampling | window (s), n | CO2 slope ppm/s, p | CH4 flux, MDF95 (nmol m-2 s-1), class |
#' |---|---|---|---|---|---|---|
#' | `emission_stem_semirigid` | clean emission (the worked-example closure) | upland tree stem, monthly survey, spring 2021 | ABB/LGR UGGA, 1 s | 765, 770 | 0.824, ~0 | 0.254, 0.0031, emission |
#' | `emission_upland_stem` | clean mid-range emission | upland tree stem, height campaign, summer 2021 | ABB/LGR UGGA, 5 s | 1209, 244 | 0.678, 1e-304 | 0.171, 0.018, emission |
#' | `below_detection_stem` | CO2 rises cleanly, CH4 flat | upland tree stem, height campaign, summer 2021 | ABB/LGR UGGA, 5 s | 1124, 227 | 0.628, 3e-256 | 0.0018, 0.019, below detection |
#' | `failed_closure_stem` | CO2 not rising (slightly falling), CH4 below detection | upland tree stem, monthly survey, winter 2021 | ABB/LGR UGGA, 1 s | 348, 351 | -0.0253, 1.4e-13 | 0.0018, 0.0075, below detection |
#' | `ambiguous_stem` | CH4 marginally detected (uptake just past MDF95) over a flat CO2 tracer | upland tree stem, monthly survey, winter 2020 | ABB/LGR UGGA, 1 s | 583, 587 | -0.0003, 0.69 | -0.0046, 0.0046, uptake |
#' | `uptake_stem` | strongest coherent stem uptake, CO2 rising, headspace drawn ~70 ppb below start | upland tree stem, monthly survey, summer 2020 | ABB/LGR UGGA, 1 s | 291, 293 | 0.730, ~0 | -0.176, 0.0078, uptake |
#' | `uptake_soil` | strong, clean soil CH4 uptake | forest soil collar, summer 2020 | ABB/LGR UGGA, 1 s | 387, 390 | 0.611, ~0 | -2.25, 0.060, uptake |
#' | `high_flux_wetland_stem` | large diffusive mangrove-stem emission | coastal wetland tree stem, dry season 2022 | Picarro G4301, 3 s | 472, 149 | 0.242, 2e-135 | 23.1, 0.0006, emission |
#' | `ebullition_floating` | one bubble burst (5 consecutive samples at 42-46 s, ~675 ppb) inside a sealed floating-chamber closure | open water, ghost-forest mangrove site, wet season 2022 | ABB/LGR GLA131, 1 s | 276, 278 | 0.114, ~0 | not processed for MDF |
#' | `emission_li7810_stem` | clean emission, third analyzer format; `Tcham` is the campaign default (24 C) | tree stem, Santa Fe 2026, spring | LI-COR LI-7810, 1 s | 367, 368 | 0.565, ~0 | 0.414, 0.0031, emission |
#'
#' `MDF95` and `class` come from each campaign's own processed table
#' (`flux_FINAL.csv`, `flux_after_tier1.csv`, the BlueFlux Picarro best-flux
#' table, the Santa Fe flux table); the CO2 slope is a linear fit on the
#' flagged rows made while building the set. The UGGA closures carry the
#' campaign's water column (`H2O_ppm`); the Picarro, GLA131 (ebullition) and
#' LI-7810 rows use the values present in their manID (0 where none was
#' logged). The `*_prec` columns hold `prec` (`c(0.35, 0.9, 200)` by default).
#'
#' @section Data release and citation:
#' The closures are released with the package under its MIT license, with
#' the permission of the data owner. Upland stem and soil closures: the
#' Yale-Myers Forest (CT, USA) chamber-flux record 2020-2021 (J. Gewirtzman
#' and colleagues, Yale School of the Environment); cite as "Gewirtzman et
#' al., Yale-Myers Forest chamber-flux record 2020-2021, distributed with the
#' stemflux R package". Wetland closures (`high_flux_wetland_stem`,
#' `ebullition_floating`): the ground component of the NASA Carbon Monitoring
#' System BlueFlux field campaign, south Florida; cite Poulter, B., et al.
#' (2023), Multi-scale observations of mangrove blue carbon ecosystem fluxes:
#' The NASA Carbon Monitoring System BlueFlux field campaign, Environmental
#' Research Letters 18, 075009, \doi{10.1088/1748-9326/acdae6}, and the
#' blueflux-ground repository (MIT license). `emission_li7810_stem`: the
#' author's Santa Fe 2026 campaign.
#'
#' @param which character vector of keys (see the table); default all ten.
#' @param dead.band `NULL` (default) keeps the shipped `flag` column (the
#'   clicked window). A number re-flags every closure as `start.time +
#'   dead.band` to the end of its segment, ignoring the clicked windows.
#' @param prec manufacturer precision `c(CO2 ppm, CH4 ppb, H2O ppm)` written
#'   to the `*_prec` columns.
#' @return list with `data` (one goFlux-shaped data.frame, all closures bound
#'   together, with `flag`, `Etime`, `start.time_corr`, `end.time_corr`,
#'   `obs.length_corr`) and `aux` (data.frame: `UniqueID`, `start.time`,
#'   `obs.length`, `Area`, `Vtot`, `Tcham`, `Pcham`, `start.time_corr`,
#'   `end.time_corr`, `kind`, `campaign_type`, `season_year`, `analyzer`,
#'   `chamber_family`).
#' @examples
#' ex <- stemflux_examples()
#' ex$aux[, c("UniqueID", "campaign_type", "analyzer")]
#' table(ex$data$UniqueID, ex$data$flag)
#' one <- stemflux_examples("failed_closure_stem")$data
#' plot(one$Etime, one$CO2dry_ppm, pch = 16, cex = 0.4, col = one$flag + 1)
#' @seealso [example_closure()]
#' @export
stemflux_examples <- function(which = NULL, dead.band = NULL, prec = c(0.35, 0.9, 200)) {
  f <- system.file("extdata", "stemflux_examples.csv", package = "stemflux")
  fa <- system.file("extdata", "stemflux_examples_aux.csv", package = "stemflux")
  if (f == "" || fa == "") stop("example files not found in the installed package")
  d <- read.csv(f, stringsAsFactors = FALSE)
  ax <- read.csv(fa, stringsAsFactors = FALSE)
  if (!is.null(which)) {
    miss <- setdiff(which, ax$UniqueID)
    if (length(miss)) stop("no example matches: ", paste(miss, collapse = ", "),
                           "; available: ", paste(ax$UniqueID, collapse = ", "))
    ax <- ax[match(which, ax$UniqueID), ]
  }
  for (cc in c("start.time", "start.time_corr", "end.time_corr")) ax[[cc]] <- as.POSIXct(ax[[cc]], tz = "UTC")
  d$POSIX.time <- as.POSIXct(d$POSIX.time, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  out <- list()
  for (i in seq_len(nrow(ax))) {
    a <- ax[i, ]
    x <- d[d$UniqueID == a$UniqueID, ]
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
  rownames(ax) <- NULL
  list(data = data, aux = ax)
}

#' One real chamber closure, ready for goFlux
#'
#' Alias for `stemflux_examples("emission_stem_semirigid")$data`: a clean
#' high-signal tree-stem closure (Yale-Myers Forest monthly survey, spring
#' 2021, ABB/LGR UGGA at 1 Hz) with the campaign's clicked window (765 s
#' starting 91 s after the recorded `start.time`, 770 flagged rows) and its
#' obs.win context. See [stemflux_examples()] for the geometry and the other
#' bundled closures.
#'
#' @inheritParams stemflux_examples
#' @return goFlux-shaped data.frame for the one closure.
#' @examples
#' d <- example_closure()
#' nrow(d); sum(d$flag)
#' @export
example_closure <- function(dead.band = NULL, prec = c(0.35, 0.9, 200)) {
  stemflux_examples("emission_stem_semirigid", dead.band = dead.band, prec = prec)$data
}
