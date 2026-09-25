# =============================================================================
# Ebullition add-on: a thin adapter over goFlux's own bubble detection and
# diffusion/ebullition partition (goFlux::goAquaFlux(), Camille Minaudo, goFlux
# >= 0.3.0). Nothing in this file is called by the core pipeline unless the
# user opts in with process_fluxes(ebullition = "flag" | "split") or calls
# detect_ebullition() directly.
# =============================================================================

#' Detect ebullition (bubble) events and partition fluxes with goFlux
#'
#' Ebullition add-on (opt-in). A thin wrapper around
#' [goFlux::goAquaFlux()], which since goFlux 0.3.0 (Camille Minaudo)
#' detects bubbling events per closure from the rolling dispersion of the
#' concentration series (adaptive median + k x MAD threshold on the rolling
#' variance of first differences), estimates each event's step magnitude with
#' a local step regression that separates the step from the diffusive trend
#' and, by default, from the transient peak and re-equilibration that follow
#' a bubble, then fits the diffusive flux with goFlux's own model selection on
#' the bubble-free part of the trace and reports the ebullitive flux from the
#' summed step magnitudes. Despite its name, `goAquaFlux()` needs only
#' `UniqueID`, `Etime`, `flag`, the gas column and `Vtot`/`Area`/`Pcham`/
#' `Tcham` (`offset`/`Vcham` are used only when `Vtot` is absent), so it
#' applies to any surface that bubbles: open water, saturated soil, flooded
#' vegetation.
#'
#' @details
#' **What this function adds** is only the goFlux-style bookkeeping: it
#' accepts a click.peak2-style data frame with several closures, forwards it
#' to `goAquaFlux()`, and returns the pieces in the shape the rest of fluxqc
#' uses. The former in-package detectors (a Hoffmann et al. 2017 style IQR
#' filter, a fixed threshold and a running-variance method) were removed once
#' goFlux provided the partition; they duplicated it with less careful step
#' estimation. Hoffmann, M., et al. (2017, Atmos. Meas. Tech. 10, 109-118,
#' \doi{10.5194/amt-10-109-2017}) remains the reference for the quartile-fence
#' family of filters; `goAquaFlux()`'s detector is different (rolling
#' dispersion with an adaptive MAD threshold) and should be cited through
#' goFlux (Rheault et al. 2024) and its changelog (Minaudo, 2026).
#'
#' **Limitations inherited from goAquaFlux.** The diffusive flux is fitted on
#' the observations *before* the first bubble (the window is restricted, not
#' de-ebulliated), so a closure whose first event comes very early leaves too
#' few points and `flux_diffusive`/`flux_total` are `NA` while
#' `flux_ebullition` is still reported; and the ebullitive rate assumes the
#' bubbles caught during a short closure represent the sampling period, which
#' for episodic ebullition is rarely true.
#'
#' @param dataframe click.peak2-style data.frame (`UniqueID`, `POSIX.time`,
#'   `flag`, `Etime`, the gas column and chamber columns).
#' @param gastype gas column, e.g. `"CH4dry_ppb"`.
#' @param H2O_col water column passed to `goAquaFlux()`; a zero column is
#'   added if absent.
#' @param ... further arguments to [goFlux::goAquaFlux()] (e.g.
#'   `bubble.window.size`, `bubble.method`, `bubble.args = list(k = 4)`,
#'   `prec`).
#' @return list of class `fluxqc_ebullition` with `events` (goFlux's bubble
#'   table: `UniqueID`, `start`, `end` (s since the first `Etime`),
#'   `magnitude` (gas units), `SE`, `slope`, `n_used`, `overshoot`, `tau`,
#'   ...; zero rows when nothing was detected), `summary` (per closure:
#'   `n_events`, `total_ebullitive` (summed magnitude), `first_bubble_time`,
#'   `flux_diffusive`, `flux_ebullitive`, `flux_total`, their SEs and
#'   `n_obs_diffusive`, straight from `goAquaFlux()$flux_summary`) and
#'   `diffusive` (goFlux's selected diffusive-model row per closure).
#' @references
#' Rheault, K., Christiansen, J. R., & Larsen, K. S. (2024). goFlux. Journal
#' of Open Source Software, 9(96), 6393. \doi{10.21105/joss.06393}
#' goFlux 0.3.0-0.4.0 changelog (2026): aquatic flux partitioning (diffusion +
#' ebullition) for floating chambers, by Camille Minaudo.
#' @family ebullition (add-on)
#' @seealso [goFlux::goAquaFlux()], [goFlux::flux.plot.aqua()], [process_fluxes()]
#' @examples
#' ex <- fluxqc_examples(c("ebullition_floating", "emission_stem_semirigid"))
#' eb <- detect_ebullition(ex$data, "CH4dry_ppb")
#' eb$summary
#' eb$events[, c("UniqueID", "start", "end", "magnitude", "SE")]
#' @export
detect_ebullition <- function(dataframe, gastype = "CH4dry_ppb", H2O_col = "H2O_ppm", ...) {
  d <- as.data.frame(dataframe)
  for (col in c("UniqueID", "Etime", "flag", gastype)) if (!col %in% names(d)) stop("column '", col, "' not found")
  if (!H2O_col %in% names(d)) d[[H2O_col]] <- 0
  r <- goFlux::goAquaFlux(d, gastype = gastype, H2O_col = H2O_col, ...)
  fs <- as.data.frame(r$flux_summary)
  ev <- if (is.null(r$bubbles)) data.frame() else as.data.frame(r$bubbles)
  if (nrow(ev)) ev <- ev[, c("UniqueID", setdiff(names(ev), "UniqueID"))]
  n_ev <- if (nrow(ev)) table(ev$UniqueID) else integer(0)
  tot <- if (nrow(ev)) tapply(ev$magnitude, ev$UniqueID, sum) else numeric(0)
  summary <- data.frame(UniqueID = fs$UniqueID,
                        n_events = as.integer(ifelse(is.na(n_ev[fs$UniqueID]), 0L, n_ev[fs$UniqueID])),
                        total_ebullitive = as.numeric(ifelse(is.na(tot[fs$UniqueID]), 0, tot[fs$UniqueID])),
                        first_bubble_time = fs$first_bubble_time,
                        flux_diffusive = fs$flux_diffusive, SE_diffusive = fs$SE_diffusive,
                        n_obs_diffusive = fs$n_obs.diffusion,
                        flux_ebullitive = fs$flux_ebullition, SE_ebullitive = fs$SE_ebullition,
                        flux_total = fs$flux_total, SE_total = fs$SE_total,
                        stringsAsFactors = FALSE)
  rownames(summary) <- NULL
  structure(list(events = ev, summary = summary,
                 diffusive = if (is.null(r$diffusive)) NULL else as.data.frame(r$diffusive),
                 gastype = gastype), class = "fluxqc_ebullition")
}

# Internal dispatcher used by process_fluxes() ONLY when ebullition != "none".
apply_ebullition <- function(d, gastype, ids, args, H2O_col) {
  eb <- do.call(detect_ebullition, c(list(dataframe = d[d$UniqueID %in% ids, ], gastype = gastype,
                                          H2O_col = H2O_col), args))
  list(events = eb$events, summary = eb$summary, diffusive = eb$diffusive)
}

# Internal: append the add-on's columns to the flux table (goFlux columns untouched)
append_ebullition_cols <- function(best, ebres, ids, mode, method) {
  best$ebullition_handling <- ifelse(best$UniqueID %in% ids, mode, "none")
  i <- match(best$UniqueID, ebres$summary$UniqueID)
  best$n_ebullition_events <- ifelse(is.na(i), 0L, as.integer(ebres$summary$n_events[i]))
  best$ebullition_flag <- best$n_ebullition_events > 0
  if (mode == "split") {
    pick <- function(col) ifelse(is.na(i), NA_real_, ebres$summary[[col]][i])
    best$flux_diffusive <- pick("flux_diffusive"); best$SE_diffusive <- pick("SE_diffusive")
    best$n_obs_diffusive <- pick("n_obs_diffusive")
    best$flux_ebullitive <- ifelse(is.na(i), 0, ebres$summary$flux_ebullitive[i])
    best$SE_ebullitive <- pick("SE_ebullitive")
    best$flux_total <- pick("flux_total"); best$SE_total <- pick("SE_total")
    best$first_bubble_time <- pick("first_bubble_time")
    best$ebullition_method <- ifelse(best$UniqueID %in% ids, method, NA_character_)
  }
  best
}
