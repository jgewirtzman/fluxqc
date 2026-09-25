# =============================================================================
# build_examples.R -- how the bundled example closures were pulled and anonymized
#
# NOT run at install time. Documents, reproducibly from the author's drive, how
# inst/extdata/stemflux_examples.csv, stemflux_examples_aux.csv and
# stemflux_examples_diagnostics.csv were made. Each closure is the campaign's
# own goFlux manID output: the full obs.win() segment (context rows flag = 0,
# operator-clicked fitting window flag = 1) with the geometry from the same
# rows. Anonymization: UniqueIDs become descriptive keys; tree tags, plot codes,
# species, operators, file names and dates are dropped; the time axis is moved
# to a synthetic calendar (2000-01-01 + one day per closure) preserving every
# interval; analyzer model, chamber family, geometry (Area, Vtot, Tcham, Pcham),
# season and year (documentation only) are kept.
#
# Sources (author's drive, as of 2026-09):
#   YMF  = ~/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens
#   BF   = ~/My Drive/Research/Blueflux/blueflux-ground
#   SF   = ~/My Drive/Research/santafe-fluxes-2026
#   GP   = ~/Downloads/Guidelines Paper/figure-build/fig03_traces/out
#
# Continuity rule applied to every clicked window (columns max_gap_s,
# max_dCH4_over_MAD, max_dCO2_over_MAD in the diagnostics): no timestamp gap
# larger than twice the logging interval, and no single-sample |dCH4| beyond
# what Gaussian noise produces over the window (the 5 x MAD guideline, which for
# 300-800 one-second differences sits at the expected maximum of ~5 x MAD; the
# ebullition closure is exempt by design). Closures replaced for failing it:
# SF19_DBH_r1 (6.1 x MAD step) and Mar_22_T1_52_FLM30_stem (45 x MAD, 4 s gap).
# =============================================================================
suppressMessages(library(data.table))
YMF <- path.expand("~/My Drive/Research/YMF Tree Microbiomes & Methane/tree-methanogens")
BF  <- path.expand("~/My Drive/Research/Blueflux/blueflux-ground")
SF  <- path.expand("~/My Drive/Research/santafe-fluxes-2026")
GP  <- path.expand("~/Downloads/Guidelines Paper/figure-build/fig03_traces/out")
OUT <- path.expand("~/My Drive/Research/stemflux/inst/extdata")

# ---- selection ---------------------------------------------------------------
# key, source, original UniqueID, campaign type, season/year, analyzer, chamber family, MDF/class source
sel <- rbindlist(list(
  list("emission_stem_semirigid", "semirigid", "2021-05-18_5086_I_7.0",  "upland tree stem, monthly survey", "spring 2021", "ABB/LGR UGGA (GLA131)", "semi-rigid sleeve chamber", "flux_FINAL.csv"),
  list("emission_upland_stem",    "height",    "20210804_Lowland_BO1_200_6_165200", "upland tree stem, height campaign", "summer 2021", "ABB/LGR UGGA (GLA131)", "rigid static stem chamber", "flux_after_tier1.csv (MDF95)"),
  list("below_detection_stem",    "height",    "20210722_Ridge_BB11_50_6_94000", "upland tree stem, height campaign", "summer 2021", "ABB/LGR UGGA (GLA131)", "rigid static stem chamber", "flux_after_tier1.csv (MDF95)"),
  list("failed_closure_stem",     "semirigid", "2021-02-25_1885_I_5_1",  "upland tree stem, monthly survey", "winter 2021", "ABB/LGR UGGA (GLA131)", "semi-rigid sleeve chamber", "flux_FINAL.csv"),
  list("ambiguous_stem",          "semirigid", "2020-12-09_5091_I_5",    "upland tree stem, monthly survey", "winter 2020", "ABB/LGR UGGA (GLA131)", "semi-rigid sleeve chamber", "flux_FINAL.csv"),
  list("uptake_stem",             "semirigid", "2020-07-20_3325_U_3",    "upland tree stem, monthly survey", "summer 2020", "ABB/LGR UGGA (GLA131)", "semi-rigid sleeve chamber", "flux_FINAL.csv"),
  list("uptake_soil",             "soil",      "20200824_2-4-U_11",      "upland forest soil collar",        "summer 2020", "ABB/LGR UGGA (GLA131)", "soil collar chamber", "flux_after_tier1.csv (MDF95)"),
  list("high_flux_wetland_stem",  "bf_lgr3",   "Oct_22_84_FLM30_stem",   "mangrove tree stem, coastal wetland","wet season 2022 (October)", "ABB/LGR GLA131 (LGR3)", "rigid stem chamber", "CH4_best_flux_lgr3_results.csv (goFlux MDF)"),
  list("ebullition_floating",     "bf_ebull",  "LGR2_2022-10-23_CP40_P06","open water, floating chamber, ghost-forest mangrove site", "wet season 2022 (October)", "ABB/LGR GLA131", "floating chamber", "campaign ebullition outputs (no MDF)"),
  list("emission_li7810_stem",    "santafe",   "SF2_DBH2_r1",             "tree stem, Santa Fe 2026", "spring 2026 (May)", "LI-COR LI-7810", "rigid stem chamber", "santafe_2026_fluxes_with_mdf.csv (CH4_MDF)")))
setnames(sel, c("key", "src", "orig", "campaign_type", "season_year", "analyzer", "chamber_family", "class_source"))

# ---- readers: return the obs.win segment for one UniqueID in goFlux columns --
std <- function(x) {
  x <- as.data.frame(x)
  if (!"H2O_ppm" %in% names(x)) x$H2O_ppm <- 0
  x$POSIX.time <- as.POSIXct(sub("Z$", "", sub("T", " ", as.character(x$POSIX.time))), tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  for (cc in c("start.time", "start.time_corr", "end.time_corr")) x[[cc]] <- as.POSIXct(sub("Z$", "", sub("T", " ", as.character(x[[cc]]))), tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  x <- x[order(x$POSIX.time), ]
  x[, c("UniqueID", "POSIX.time", "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm", "start.time", "obs.length", "Area", "Vtot", "Tcham", "Pcham",
        "flag", "start.time_corr", "end.time_corr")]
}
read_src <- list(
  semirigid  = function(id) std(as.data.frame(readRDS(file.path(GP, "semirigid_tree_manID_clicked.rds")))[ , ][which(as.data.frame(readRDS(file.path(GP, "semirigid_tree_manID_clicked.rds")))$UniqueID == id), ]),
  height     = function(id) { m <- fread(file.path(YMF, "data/processed/flux/lgr_manual_identification_results_final.csv")); std(m[UniqueID == id]) },
  soil       = function(id) { m <- rbind(fread(file.path(YMF, "data/processed/flux/lgr_manual_identification_results_soil.csv")),
                                          fread(file.path(YMF, "data/processed/flux/lgr_manual_identification_results_december_soil.csv")), fill = TRUE); std(m[UniqueID == id]) },
  bf_picarro = function(id) { m <- fread(file.path(BF, "intermediate/results_trees/picarro_manual_identification_CH4_results.csv")); std(m[UniqueID == id]) },
  bf_lgr3    = function(id) { m <- fread(file.path(BF, "intermediate/results_trees/lgr3_manual_identification_results.csv")); std(m[UniqueID == id]) },
  santafe    = function(id) { e <- new.env(); load(file.path(SF, "RData/manID.RData"), envir = e); x <- e$manID[e$manID$UniqueID == id, ]; std(x) },
  bf_ebull   = function(id) {
    # trace as segmented by the campaign's detect_ebullition.R (clock offset applied, 20 s trimmed each end); no clicked window
    tr <- readRDS(file.path(BF, "output/ebullition/all_traces.rds"))[[id]]$trace
    w <- read.csv(file.path(BF, "output/data_products/soil_water_surface_fluxes_ORIGINAL.csv"))
    tc <- mean(w$air_temp[w$plot == "CP40" & w$date == "2022-10-23" & w$surface_type == "water"], na.rm = TRUE)   # site-day mean
    data.frame(UniqueID = id, POSIX.time = tr$datetime, CO2dry_ppm = tr$CO2_ppm, CH4dry_ppb = tr$CH4_ppm * 1000,
               H2O_ppm = if (all(is.na(tr$H2O_ppm))) 0 else tr$H2O_ppm, start.time = min(tr$datetime), obs.length = round(max(tr$elapsed_sec)),
               Area = 324.3, Vtot = 4.318, Tcham = round(tc, 2), Pcham = 101.325,          # goflux_reprocess_ebullition.R CHAMBER_PARAMS (LGR)
               flag = 1, start.time_corr = min(tr$datetime), end.time_corr = max(tr$datetime)) })

# ---- MDF / class lookups --------------------------------------------------------
fin <- fread(file.path(YMF, "outputs/data/flux_FINAL.csv")); t1 <- fread(file.path(YMF, "outputs/data/flux_after_tier1.csv"))
bfp <- fread(file.path(BF, "intermediate/results_trees/CH4_best_flux_picarro_results.csv"))
bfl <- fread(file.path(BF, "intermediate/results_trees/CH4_best_flux_lgr3_results.csv"))
sfr <- fread(file.path(SF, "results/santafe_2026_fluxes_with_mdf.csv"))
lookup <- function(src, id) switch(src,
  semirigid = fin[UniqueID == id, .(flux = best.flux, mdf = MDF, model = NA_character_)],
  height = , soil = t1[UniqueID == id, .(flux = best.flux, mdf = MDF95, model = model)],
  bf_picarro = bfp[UniqueID == id, .(flux = best.flux, mdf = MDF, model = model)],
  bf_lgr3 = bfl[UniqueID == id, .(flux = best.flux, mdf = MDF, model = model)],
  santafe = sfr[UniqueID == id, .(flux = CH4_best.flux, mdf = CH4_MDF, model = CH4_model)],
  bf_ebull = data.table(flux = NA_real_, mdf = NA_real_, model = NA_character_))

# ---- build ---------------------------------------------------------------------
T0 <- as.POSIXct("2000-01-01 00:00:00", tz = "UTC")
data <- list(); aux <- list(); diag <- list()
for (i in seq_len(nrow(sel))) {
  s <- sel[i]; x <- read_src[[s$src]](s$orig)
  stopifnot(nrow(x) > 0)
  shift <- as.numeric(difftime(T0 + (i - 1) * 86400 + 3600, x$start.time[1], units = "secs"))   # synthetic calendar
  fl <- x[x$flag == 1, ]; et <- as.numeric(difftime(fl$POSIX.time, x$start.time_corr[1], units = "secs"))
  a <- summary(lm(fl$CO2dry_ppm ~ et))$coefficients; b <- summary(lm(fl$CH4dry_ppb ~ et))
  lk <- lookup(s$src, s$orig)
  cls <- if (is.na(lk$flux)) NA_character_ else if (lk$flux > lk$mdf) "emission" else if (lk$flux < -lk$mdf) "uptake" else "below detection"
  data[[i]] <- data.frame(UniqueID = s$key, POSIX.time = format(x$POSIX.time + shift, "%Y-%m-%d %H:%M:%OS3"),
                          t_s = round(as.numeric(difftime(x$POSIX.time, x$start.time[1], units = "secs")), 3),
                          CO2dry_ppm = round(x$CO2dry_ppm, 4), CH4dry_ppb = round(x$CH4dry_ppb, 3), H2O_ppm = round(x$H2O_ppm, 1), flag = x$flag)
  aux[[i]] <- data.frame(UniqueID = s$key, start.time = format(x$start.time[1] + shift, "%Y-%m-%d %H:%M:%S"), obs.length = x$obs.length[1],
                         Area = x$Area[1], Vtot = x$Vtot[1], Tcham = x$Tcham[1], Pcham = x$Pcham[1],
                         start.time_corr = format(x$start.time_corr[1] + shift, "%Y-%m-%d %H:%M:%S"),
                         end.time_corr = format(x$end.time_corr[1] + shift, "%Y-%m-%d %H:%M:%S"),
                         kind = "real", campaign_type = s$campaign_type, season_year = s$season_year, analyzer = s$analyzer,
                         chamber_family = s$chamber_family)
  # continuity of the flagged window: timestamp gaps and single-sample steps relative to the window's own noise
  dtw <- diff(as.numeric(fl$POSIX.time)); dch <- abs(diff(fl$CH4dry_ppb)); dco <- abs(diff(fl$CO2dry_ppm))
  diag[[i]] <- data.frame(key = s$key, campaign_type = s$campaign_type, season_year = s$season_year, analyzer = s$analyzer,
                          chamber_family = s$chamber_family, sampling_interval_s = round(median(diff(et)), 1),
                          max_gap_s = round(max(dtw), 2), max_dCH4_over_MAD = round(max(dch) / mad(dch), 1),
                          max_dCO2_over_MAD = round(max(dco) / mad(dco), 1),
                          window_length_s = round(diff(range(et))), n_flagged = nrow(fl), n_segment = nrow(x),
                          co2_slope_ppm_s = signif(a[2, 1], 3), co2_p = signif(a[2, 4], 2), ch4_slope_ppb_s = signif(b$coefficients[2, 1], 3),
                          ch4_r2 = round(b$r.squared, 3), ch4_flux_nmol_m2_s = signif(lk$flux, 3), MDF95 = signif(lk$mdf, 3),
                          class = cls, class_source = s$class_source)
}
diag <- do.call(rbind, diag)
diag$verdict <- c("valid, detected: the worked-example closure (1 Hz, CH4 r2 0.997)",
                  "valid, detected: clean mid-range upland emission, 5 s logging",
                  "valid, CH4 below detection: CO2 rises cleanly, CH4 flat",
                  "failed closure: CO2 not rising (slightly falling), CH4 below detection; the only such stem closure in the monthly-survey archive",
                  "ambiguous: CH4 marginally detected (uptake just past MDF95) over a flat CO2 tracer; the only such closure in the archive",
                  "valid, detected uptake: strongest coherent negative stem closure, CO2 rising, headspace drawn ~70 ppb below start",
                  "valid, detected uptake: strong clean soil CH4 uptake",
                  "valid, detected: large diffusive mangrove-stem emission (continuous 1 s timestamps, CH4 steps 3.2 x MAD). No Picarro G4301 stem closure of the campaign passes the continuity test (irregular 3-7 s logging, single-sample steps >= 11 x MAD), so that format is not represented",
                  "ebullition: one bubble burst (5 consecutive samples at 42-46 s, ~675 ppb) inside a sealed floating-chamber closure; no clicked window (campaign-segmented trace)",
                  "valid, detected: clean LI-7810 stem emission, third analyzer format (continuous 1 s timestamps, CH4 steps 4.1 x MAD); Tcham is the campaign default (24 C) as no chamber temperature was logged")
write.csv(do.call(rbind, data), file.path(OUT, "stemflux_examples.csv"), row.names = FALSE)
write.csv(do.call(rbind, aux),  file.path(OUT, "stemflux_examples_aux.csv"), row.names = FALSE)
write.csv(diag, file.path(OUT, "stemflux_examples_diagnostics.csv"), row.names = FALSE)
options(width = 250); print(diag[, c("key", "season_year", "sampling_interval_s", "window_length_s", "n_flagged", "n_segment", "co2_slope_ppm_s", "co2_p", "ch4_flux_nmol_m2_s", "MDF95", "class")])
