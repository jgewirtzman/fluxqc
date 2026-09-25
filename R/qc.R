#' Toggleable physical quality screens for a flux table (never deletes rows)
#'
#' Applies a set of independent screens to a [goFlux::best.flux()] (or
#' [goFlux::goFlux()]) table and appends one logical column per screen
#' (`qc_<screen>`, `TRUE` = the screen fired), plus `qc_any` and a `qc_note`
#' string listing the screens that fired. Rows are retained and flagged only.
#' Screens that cannot be evaluated for a row (missing inputs)
#' give `NA`, which is *not* a clean bill of health.
#'
#' Each screen is a simple post-hoc toggle: `TRUE`, `FALSE` or a parameter
#' list in `screens`. None depends on what surface the chamber enclosed; the
#' only one to switch off deliberately is `co2_tracer`, for surfaces without
#' respiratory CO2 efflux (open water, lit foliage, dead wood).
#'
#' @section Screens:
#' \describe{
#'   \item{`c0` (`k = 1.5`)}{starting concentration `C0` greater than `k` times
#'     the median `C0` of its group (`group` column, e.g. campaign or day):
#'     the chamber was not flushed to ambient at closure or carried over a
#'     previous sample. The default `k` sits in the middle of a wide gap in the
#'     source campaigns' C0 distributions; check the sensitivity line printed
#'     with `verbose = TRUE`.}
#'   \item{`co2_tracer` (`min_slope = 0, require_detection = TRUE`)}{on live
#'     tissue the sealed chamber must accumulate CO2 from respiration; the
#'     screen fires when the CO2 slope is at or below `min_slope` (ppm s-1)
#'     *or*, with `require_detection`, when the CO2 rise is not distinguishable
#'     from noise: `best.flux <= MDF_emp` (or goFlux's `MDF`) in the `co2`
#'     table, or, from `traces`, a rise over the window (`slope x span`) of
#'     less than 2 sigma of the closure's Allan noise. Either way the chamber
#'     leaked, was not sealed, or was on dead/dormant tissue. The CO2 slope
#'     comes from, in order: the `co2` table (`LM.slope` by `UniqueID`),
#'     `traces` (linear fit of `CO2dry_ppm` on flagged rows), or
#'     `fluxes$LM.slope` when `fluxes` is itself the CO2 table. Turn it off
#'     (`co2_tracer = FALSE`) for surfaces without respiratory CO2 efflux (open
#'     water, lit foliage, dead wood).}
#'   \item{`convex` (`p = 0.05`)}{accelerating (convex) concentration change,
#'     which no chamber mechanism produces and which points at a leak, a
#'     disturbance or a delayed seal. With `traces`: a quadratic fit on the
#'     flagged rows whose curvature has the sign of the net linear slope and
#'     is significant at `p`. Without `traces`: `HM.k < 0` (only possible if
#'     `goFlux()` was run with `k.min < 0`), otherwise `NA`.}
#'   \item{`min_window` (`secs = 60`)}{the flagged window (`obs.length_corr`
#'     from `traces`, else `nb.obs` assuming 1 Hz) is shorter than `secs`.}
#'   \item{`ambient_start` (`sigma = 3, secs = 10`)}{the mean of the first
#'     `secs` seconds of the flagged window differs from the pre-closure
#'     ambient by more than `sigma` times the closure's Allan precision: the
#'     chamber was already enriched when the window started (late window
#'     start, or contamination). Ambient is the mean of the rows before the
#'     auxfile `start.time` (the shoulder from `obs.win()`), or, if that
#'     column is absent, of rows with `flag == 0` and `Etime < 0`. Needs
#'     `traces` with shoulder rows; `NA` otherwise.}
#'   \item{`noisy` (`k = 1.5`)}{the closure's own precision, MAD(dx)/sqrt(2) on
#'     its flagged rows ([precision_mad()]), exceeds `k` times the campaign
#'     (group) sigma that underlies `MDF_emp` (`sigma_emp` from
#'     [flag_detection()] if present, else MAD over the whole record per
#'     `group`). The per-closure sigma is a *diagnostic*, not an alternative
#'     MDF basis: a closure much noisier than its campaign was disturbed
#'     (moves, breaths, a leaking seal) or the analyzer was unstable. Note
#'     text: "closure sigma = X ppb = Yx campaign sigma (disturbed closure?)".
#'     Columns `qc_noisy_sigma`, `qc_noisy_ratio`. Needs `traces`.}
#'   \item{`leak` (`rate = NULL`)}{a bench-measured leak rate `rate` (fraction
#'     of the headspace exchanged per second, e.g. the decay constant of a
#'     spiked, sealed chamber) times the largest headspace-ambient difference
#'     during the closure gives the apparent slope a leak alone could produce;
#'     if that, times `flux.term`, is at least the row's `MDF_emp` (from
#'     [flag_detection()]; goFlux's `MDF` if `MDF_emp` is absent), the flux
#'     cannot be separated from leakage. Uses `traces` for the concentration
#'     range, else `abs(Ct - C0)`. Skipped when `rate` is `NULL`.}
#' }
#'
#' @param fluxes data.frame from [goFlux::best.flux()] / [goFlux::goFlux()],
#'   optionally after [flag_detection()] (needed for `leak`).
#' @param traces flagged concentration data the fluxes were computed from
#'   (click.peak2-style); optional but enables `convex`, `ambient_start`, the
#'   trace-based `co2_tracer` and `min_window` lengths.
#' @param screens named list; each element `TRUE` (defaults), `FALSE`/`NULL`
#'   (skip) or a list of parameters (merged with the defaults).
#' @param group name of a column in `fluxes` defining groups for the `c0`
#'   median (campaign, day, instrument); `NULL` uses all rows.
#' @param gastype gas column of `traces` (inferred when unique).
#' @param co2 optional CO2 flux table (`best.flux()` output on `CO2dry_ppm`)
#'   for the `co2_tracer` screen.
#' @param verbose print a per-screen summary.
#' @return `fluxes` with `qc_<screen>` columns, `qc_any` and `qc_note`. The
#'   parameters used are stored in the attribute `"qc_screens"`.
#' @seealso [qc_c0_screen()] (back-compatible wrapper), [flag_detection()]
#' @examples
#' d <- example_closure()
#' best <- goFlux::best.flux(goFlux::goFlux(d, "CH4dry_ppb"))
#' co2  <- goFlux::best.flux(goFlux::goFlux(d, "CO2dry_ppm"))
#' out <- qc_screens(best, traces = d, gastype = "CH4dry_ppb", co2 = co2,
#'                   screens = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE,
#'                                  min_window = list(secs = 120), ambient_start = TRUE,
#'                                  noisy = list(k = 1.5)))
#' # open water: no respiratory CO2, so switch the tracer screen off
#' # qc_screens(best_water, traces = water, screens = list(c0 = TRUE, co2_tracer = FALSE))
#' out[, grep("^qc_", names(out))]
#' @export
qc_screens <- function(fluxes, traces = NULL,
                       screens = list(c0 = list(k = 1.5), co2_tracer = list(min_slope = 0),
                                      convex = TRUE, min_window = list(secs = NULL),
                                      ambient_start = list(sigma = 3, secs = NULL),
                                      noisy = list(k = 1.5), leak = list(rate = NULL)),
                       group = NULL, gastype = NULL, co2 = NULL, verbose = TRUE) {
  fx <- as.data.frame(fluxes)
  if (!"UniqueID" %in% names(fx)) stop("'fluxes' needs a 'UniqueID' column")
  defaults <- list(c0 = list(k = 1.5), co2_tracer = list(min_slope = 0, require_detection = TRUE),
                   convex = list(p = 0.05),
                   min_window = list(secs = 60), ambient_start = list(sigma = 3, secs = 10),
                   noisy = list(k = 1.5), leak = list(rate = NULL))
  unknown <- setdiff(names(screens), names(defaults))
  if (length(unknown)) stop("unknown screen(s): ", paste(unknown, collapse = ", "))
  cfg <- list()
  for (s in names(defaults)) {
    v <- screens[[s]]
    if (is.null(v) || isFALSE(v)) next
    p <- defaults[[s]]
    if (is.list(v)) for (nm in names(v)) if (!is.null(v[[nm]])) p[[nm]] <- v[[nm]]
    cfg[[s]] <- p
  }
  n <- nrow(fx)
  uid <- as.character(fx$UniqueID)
  tr_flag <- NULL
  if (!is.null(traces)) {
    if (!"flag" %in% names(traces)) stop("'traces' needs a 'flag' column")
    tr_flag <- traces[!is.na(traces$flag) & traces$flag == 1, ]
    if (is.null(gastype)) {
      g <- intersect(goflux_gas_columns(), names(traces))
      gastype <- if (length(g) == 1) g else NULL
    }
  }
  per_uid <- function(f) {                       # apply f to each closure's flagged rows
    sp <- split(tr_flag, tr_flag$UniqueID)
    r <- vapply(sp, f, numeric(1))
    r[match(uid, names(sp))]
  }
  notes <- rep("", n)
  fired <- function(name, flag, note) {
    fx[[paste0("qc_", name)]] <<- flag
    hit <- !is.na(flag) & flag
    notes[hit] <<- ifelse(nzchar(notes[hit]), paste0(notes[hit], "; ", note), note)
    if (verbose) message(sprintf("qc %-14s fired %d, NA %d, of %d", name, sum(hit), sum(is.na(flag)), n))
  }

  # c0
  if (!is.null(cfg$c0)) {
    if (!"C0" %in% names(fx)) fired("c0", rep(NA, n), "")
    else {
      g <- if (is.null(group)) rep("all", n) else as.character(fx[[group]])
      med <- tapply(fx$C0, g, stats::median, na.rm = TRUE)
      ratio <- as.numeric(fx$C0) / as.numeric(med[g])
      fx$qc_c0_ratio <- ratio
      fired("c0", ifelse(is.finite(ratio), ratio > cfg$c0$k, NA),
            sprintf("C0 > %gx group median", cfg$c0$k))
      if (verbose) message("  c0 threshold sensitivity: ",
                           paste(sprintf("%.1fx->%d", c(1.2, 1.5, 2, 5, 10),
                                         vapply(c(1.2, 1.5, 2, 5, 10), function(m) sum(ratio > m, na.rm = TRUE), numeric(1))),
                                 collapse = "  "))
    }
  }
  # co2_tracer
  if (!is.null(cfg$co2_tracer)) {
    slope <- rep(NA_real_, n)
    if (!is.null(co2) && "LM.slope" %in% names(co2)) slope <- co2$LM.slope[match(uid, co2$UniqueID)]
    else if (!is.null(tr_flag) && "CO2dry_ppm" %in% names(tr_flag))
      slope <- per_uid(function(d) if (nrow(d) < 3) NA_real_ else unname(stats::coef(stats::lm(d$CO2dry_ppm ~ d$Etime))[2]))
    else if (identical(gastype, "CO2dry_ppm") && "LM.slope" %in% names(fx)) slope <- fx$LM.slope
    fx$qc_co2_slope <- slope
    flag <- ifelse(is.na(slope), NA, slope <= cfg$co2_tracer$min_slope)
    # noise-aware part: the CO2 rise must also exceed its own noise
    if (isTRUE(cfg$co2_tracer$require_detection)) {
      undetected <- rep(NA, n)
      if (!is.null(co2) && "UniqueID" %in% names(co2)) {
        j <- match(uid, co2$UniqueID)
        mdfc <- if ("MDF_emp" %in% names(co2)) co2$MDF_emp[j] else if ("MDF" %in% names(co2)) co2$MDF[j] else NA
        fl <- if ("best.flux" %in% names(co2)) co2$best.flux[j] else co2$LM.flux[j]
        undetected <- ifelse(is.na(mdfc) | is.na(fl), NA, fl <= mdfc)
      } else if (!is.null(tr_flag) && "CO2dry_ppm" %in% names(tr_flag)) {
        undetected <- per_uid(function(d) {
          if (nrow(d) < 3) return(NA_real_)
          s <- precision_allan(d$CO2dry_ppm); tt <- diff(range(d$Etime))
          sl <- unname(stats::coef(stats::lm(d$CO2dry_ppm ~ d$Etime))[2])
          as.numeric(sl * tt <= 2 * s)          # rise over the window below 2 sigma of the noise
        }) == 1
      }
      flag <- ifelse(is.na(flag), NA, flag | (!is.na(undetected) & undetected))
    }
    fired("co2_tracer", flag, sprintf("CO2 slope <= %g", cfg$co2_tracer$min_slope))
  }
  # convex
  if (!is.null(cfg$convex)) {
    if (!is.null(tr_flag) && !is.null(gastype)) {
      cv <- per_uid(function(d) {
        if (nrow(d) < 6) return(NA_real_)
        m <- stats::lm(d[[gastype]] ~ d$Etime + I(d$Etime^2))
        co <- summary(m)$coefficients
        if (nrow(co) < 3) return(NA_real_)
        net <- unname(stats::coef(stats::lm(d[[gastype]] ~ d$Etime))[2])   # sign of the overall change
        as.numeric(sign(co[3, 1]) == sign(net) && co[3, 4] < cfg$convex$p)
      })
      flag <- ifelse(is.na(cv), NA, cv == 1)
    } else if ("HM.k" %in% names(fx)) flag <- ifelse(is.na(fx$HM.k), NA, fx$HM.k < 0)
    else flag <- rep(NA, n)
    fired("convex", flag, "convex (accelerating) trace")
  }
  # min_window
  if (!is.null(cfg$min_window)) {
    len <- if (!is.null(tr_flag) && "obs.length_corr" %in% names(tr_flag))
      per_uid(function(d) d$obs.length_corr[1]) else if ("nb.obs" %in% names(fx)) fx$nb.obs else rep(NA_real_, n)
    fired("min_window", ifelse(is.na(len), NA, len < cfg$min_window$secs),
          sprintf("window < %g s", cfg$min_window$secs))
  }
  # ambient_start
  if (!is.null(cfg$ambient_start)) {
    if (is.null(traces) || is.null(gastype)) flag <- rep(NA, n)
    else {
      pre <- if ("start.time" %in% names(traces)) traces[traces$POSIX.time < traces$start.time, ]
             else traces[!is.na(traces$flag) & traces$flag == 0 & traces$Etime < 0, ]
      amb <- tapply(pre[[gastype]], pre$UniqueID, mean, na.rm = TRUE)
      dev <- per_uid(function(d) {
        a <- amb[unique(d$UniqueID)[1]]
        if (is.null(a) || is.na(a)) return(NA_real_)
        s <- precision_allan(d[[gastype]])
        first <- mean(d[[gastype]][d$Etime <= cfg$ambient_start$secs], na.rm = TRUE)
        abs(first - a) / s
      })
      fx$qc_ambient_dev_sigma <- dev
      flag <- ifelse(is.na(dev), NA, dev > cfg$ambient_start$sigma)
    }
    fired("ambient_start", flag, sprintf("start > %g sigma from ambient", cfg$ambient_start$sigma))
  }
  # noisy: per-closure sigma (diagnostic) vs the campaign sigma used for MDF_emp
  if (!is.null(cfg$noisy)) {
    if (is.null(tr_flag) || is.null(gastype)) flag <- rep(NA, n)
    else {
      s_clo <- per_uid(function(d) precision_mad(d[[gastype]]))
      s_camp <- if ("sigma_emp" %in% names(fx)) fx$sigma_emp
                else if (is.null(group)) rep(precision_mad(traces[[gastype]]), n)
                else {
                  gv <- as.character(fx[[group]]); tg <- gv[match(traces$UniqueID, fx$UniqueID)]
                  sg <- vapply(unique(gv), function(g) precision_mad(traces[[gastype]][!is.na(tg) & tg == g]), numeric(1))
                  unname(sg[gv])
                }
      ratio <- s_clo / s_camp
      fx$qc_noisy_sigma <- s_clo; fx$qc_noisy_ratio <- ratio
      flag <- ifelse(is.finite(ratio), ratio > cfg$noisy$k, NA)
      un <- if (grepl("ppb", gastype)) "ppb" else if (grepl("ppm", gastype)) "ppm" else ""
      note_v <- sprintf("closure sigma = %.2f %s = %.1fx campaign sigma (disturbed closure?)", s_clo, un, ratio)
    }
    if (exists("note_v", inherits = FALSE)) {
      hit <- !is.na(flag) & flag
      fx$qc_noisy <- flag
      notes[hit] <- ifelse(nzchar(notes[hit]), paste0(notes[hit], "; ", note_v[hit]), note_v[hit])
      if (verbose) message(sprintf("qc %-14s fired %d, NA %d, of %d", "noisy", sum(hit), sum(is.na(flag)), n))
      rm(note_v)
    } else fired("noisy", flag, "")
  }
  # leak
  if (!is.null(cfg$leak)) {
    if (is.null(cfg$leak$rate)) { if (verbose) message("qc leak           skipped (no bench leak rate)") }
    else {
      dC <- if (!is.null(tr_flag) && !is.null(gastype))
        per_uid(function(d) diff(range(d[[gastype]], na.rm = TRUE)))
      else if (all(c("C0", "Ct") %in% names(fx))) abs(fx$Ct - fx$C0) else rep(NA_real_, n)
      app <- cfg$leak$rate * dC * fx$flux.term
      fx$qc_leak_apparent_flux <- app
      mdfv <- if ("MDF_emp" %in% names(fx)) fx$MDF_emp else if ("MDF" %in% names(fx)) fx$MDF else rep(NA_real_, n)
      fired("leak", ifelse(is.na(app) | is.na(mdfv), NA, app >= mdfv), "leak-equivalent flux >= MDF")
    }
  }
  qc_cols <- grep("^qc_(c0|co2_tracer|convex|min_window|ambient_start|noisy|leak)$", names(fx), value = TRUE)
  if (length(qc_cols)) {
    m <- as.matrix(fx[, qc_cols, drop = FALSE])
    fx$qc_any <- apply(m, 1, function(r) if (all(is.na(r))) NA else any(r, na.rm = TRUE))
  } else fx$qc_any <- NA
  fx$qc_note <- ifelse(nzchar(notes), notes, NA_character_)
  attr(fx, "qc_screens") <- cfg
  fx
}

#' Chamber-contamination screen on the starting concentration (C0)
#'
#' Back-compatible wrapper around [qc_screens()] running only the `c0`
#' screen. Kept for scripts written against fluxqc 0.1.0.
#'
#' @inheritParams qc_screens
#' @param flux.result data.frame with a `C0` column.
#' @param k multiplier of the group median above which a measurement fails.
#' @param C0.col name of the starting-concentration column.
#' @return `flux.result` with `qc_c0_ratio`, `qc_c0_pass` (logical, `NA` if C0
#'   missing), `qc_c0_reason`, plus the `qc_c0` / `qc_any` / `qc_note`
#'   columns of [qc_screens()].
#' @examples
#' res <- data.frame(UniqueID = paste0("m", 1:6), C0 = c(2100, 2050, 1980, 2300, 11700, 2200))
#' qc_c0_screen(res)[, c("UniqueID", "C0", "qc_c0_ratio", "qc_c0_pass")]
#' @export
qc_c0_screen <- function(flux.result, k = 1.5, C0.col = "C0", group = NULL, verbose = TRUE) {
  if (!C0.col %in% names(flux.result)) stop("column '", C0.col, "' not found")
  fx <- flux.result
  if (C0.col != "C0") fx$C0 <- fx[[C0.col]]
  out <- qc_screens(fx, screens = list(c0 = list(k = k)), group = group, verbose = verbose)
  out$qc_c0_pass <- ifelse(is.na(out$qc_c0), NA, !out$qc_c0)
  out$qc_c0_reason <- ifelse(!is.na(out$qc_c0) & out$qc_c0,
                             sprintf("C0 %.0f = %.1fx group median; chamber not at ambient at closure",
                                     out$C0, out$qc_c0_ratio), NA_character_)
  out
}
