# Contributing fluxqc work back to goFlux

goFlux (Qepanna/goFlux, v0.2.0, GPL-2) is distributed from GitHub, not CRAN:
its `DESCRIPTION` has no CRAN fields beyond the usual, the `README.md` tells
users to `remotes::install_github("Qepanna/goFlux")`, and `citation()` points
to the JOSS paper (doi:10.21105/joss.06393). The mentions of CRAN in the README
refer to *other* packages (HMR, flux, gasfluxes) that goFlux compares itself
to. So contributions go through GitHub issues and pull requests, and
`Remotes: Qepanna/goFlux` in fluxqc's `DESCRIPTION` is the right dependency
declaration.

Maintainer: Karelle Rheault <karh@ign.ku.dk> (cre, aut); co-author Klaus
Steenberg Larsen. Since 0.3.0 the aquatic/ebullition partition is by Camille
Minaudo (GitHub @camilleminaudo), the natural contact for anything about
`goAquaFlux()`/`find.bubbles()`; Katharina Jentzsch (@KatharinaJentzsch)
contributed the importer fixes (0.3.1-0.3.3, UniqueID generation in 0.3.2).

## (a) Bug reports, re-verified against goFlux 0.4.0 (2026-09-23; snapshot in /tmp/goflux040)

Status of the five issues first written against v0.2.0, plus one new one:

| # | issue | status in 0.4.0 |
|---|---|---|
| 1 | `flux2pdf()` class check breaks with ggplot2 >= 3.5 | **fixed** (0.4.0 `R/flux2pdf.R` uses `vapply(plot.list, inherits, logical(1), what = "ggplot")`); fluxqc keeps `flux_pdf()` only for the grouped-page layout and `plot_env` fallback |
| 2 | `import2RData()` writes into `getwd()/RData`; subfolders not searched | **partly**: `recursive = TRUE` now for DX4015, G2201i, G2508, G4301 and the GT5000 `.TXT` branch; still non-recursive for eosMX12, EGM5, GAIA, HT8850, N2OM1, N2Oi2, skyline, uCH4, uN2O and the GT5000 `.txt` branch; all 24 `import.*()` still hard-code `paste(getwd(), "RData")` for `save = TRUE` |
| 3 | `invokeRestart("muffleError")` in `import2RData()` | **open**: 23 occurrences in 0.4.0 `R/import2RData.R` |
| 4 | `import.LI7810(dates = )` | **open**: no `dates` argument; the patch below is rebased onto the 0.4.0 file and re-tested |
| 5 | `goFlux::import.LI7810()` fails on `match.call()[[1]]` | **fixed** (the `*_import` deprecation shims are gone from every importer in 0.4.0) |
| 6 | `crop.meas(max.obs.length = "aux")` errors | **new**, see below |

### Issue 2 (still open in 0.4.0): `import2RData()` output folder and inconsistent recursion

**Title:** import2RData(): output forced into getwd()/RData; `recursive`
differs between instruments

Every `import.*()` with `save = TRUE` does `RData_folder <- paste(getwd(),
"RData", sep = "/")` regardless of `path`, so a project that keeps raw data
and processed data in different places must `setwd()` around the call.
Suggest an `outdir` argument (default `file.path(path, "RData")`) or make
`save = FALSE` the default. The recursion into date-named subfolders was
added in 0.3.x for some instruments (DX4015, G2201i, G2508, G4301, GT5000
`.TXT`) but not for eosMX12, EGM5, GAIA, HT8850, N2OM1, N2Oi2, skyline, uCH4,
uN2O or the GT5000 `.txt` branch; zero-byte and non-data (`_b`, `_l`, `_p`)
files are still matched. Suggest `recursive = TRUE` everywhere, an optional
`pattern` argument and skipping zero-byte files. fluxqc's
`import2RData_flat()` shows the alternative.

### Issue 3 (still open in 0.4.0): "no restart 'muffleError' found"

**Title:** import2RData(): withCallingHandlers() calls a non-existent
'muffleError' restart, turning any per-file error into an uninformative
failure

`import2RData()` wraps each `import.*()` call in `withCallingHandlers(...,
error = function(e) { errs <<- c(errs, conditionMessage(e));
invokeRestart("muffleError") })`. `muffleWarning` and `muffleMessage` are
real restarts; there is no `muffleError` restart in R. As soon as one file
throws (a truncated file, a non-data `_b`/`_l` file matched by the pattern),
the handler itself errors with `no restart 'muffleError' found`, the original
message is lost and nothing is returned for the files that did import.
Reproduce: write `"this is not an analyzer file"` to `bad_f0001.txt` in a
folder and call `import2RData(path, "UGGA", "mdy", "UTC")`. Fix: `tryCatch()`
for the error branch (record `basename(file_list[i])` and the message,
return `NULL`), keeping `withCallingHandlers()` for warnings and messages;
report `errs` after the progress bar as the code already intends. The same
block appears 23 times, so the change is mechanical.

### Issue 4 (still open in 0.4.0): `import.LI7810(dates = NULL)` -- patch ready

`patches/import.LI7810-dates.patch` is a unified diff against the 0.4.0
`R/import.LI7810.R` (`patches/import.LI7810.R` is the whole modified file).
R only, no `awk`. With `dates = NULL` the function is byte-for-byte the 0.4.0
one (`all.equal()` on the bundled `LI7810.data` is `TRUE`, 330 rows); with
`dates = "2022-12-05"` it returns the same 330 rows; a non-matching date stops
with a clear message. Apply with `git apply
inst/upstream/patches/import.LI7810-dates.patch` in a goFlux clone.

**PR title:** import.LI7810(): add `dates` argument to subset campaign-long
.data files before parsing

**PR description:** LI-7810 `.data` files grow to 100-200 MB over a
campaign (~1.5 million rows); `import.LI7810()` reads and time-parses every
row even when only a few field days are needed. This PR adds `dates = NULL`:
when a character vector of dates (as written in the file's `DATE` column) is
given, the file is read with `readLines()`, the `DATAH` header row is kept
and only the `DATA` rows whose `DATE` field is in `dates` are passed to
`read.delim(text = ...)`; everything after that is unchanged. Roxygen
`@param` added. fluxqc's `import_li7810_subset()` forwards to this argument
when it is present and falls back to an external subset otherwise.

### Issue 6 (new in 0.4.0): `crop.meas(max.obs.length = "aux")` fails with a NULL join

**Title:** crop.meas(): `max.obs.length = "aux"` errors ("`x` and `y` must
share the same src ... `y` is NULL")

In 0.4.0 `R/crop.meas.R`, `unique_obs.length` is initialised to `NULL` in
the "Assign NULL to variables without binding" line (line 162) and, unlike
`unique_deadband` (line 179) and `unique_crop_start` (line 185), is never
assigned from the auxfile before `meta <- left_join(meta, unique_obs.length,
by = "UniqueID")` (line 199), so the per-UniqueID path documented for
`max.obs.length` always fails:

```r
library(goFlux); data(manID.UGGA)
aux <- data.frame(UniqueID = unique(manID.UGGA$UniqueID), max.obs.length = 120)
crop.meas(manID.UGGA, auxfile = aux, max.obs.length = "aux")
#> Error in `auto_copy()`: `x` and `y` must share the same src. ... `y` is `NULL`.
crop.meas(manID.UGGA, max.obs.length = 120)   # the scalar form works
```

**Fix:** add `unique_obs.length <- select(auxfile, UniqueID, max.obs.length)`
in the `max.obs.length == "aux"` branch, mirroring the `deadband` and
`crop_start`/`crop_end` branches. fluxqc's `windows_from_table()` works
around it by calling `crop.meas()` once per closure with scalar arguments.

## (b) What belongs upstream and what stays in fluxqc

| fluxqc function | Upstream candidate? | Proposal |
|---|---|---|
| `click_peak2_stacked()` | **Yes** | `click.peak2(..., gases = NULL)`: when more than one gas is given, draw one `split.screen()` panel per gas sharing the time axis and identify on the top one. Output is unchanged. Small diff, no new dependency. |
| `precision_allan()`, `precision_by_id()` | **Yes** | A `prec.empirical()` helper (or a `prec = "empirical"` option in `goFlux()`) that fills the `*_prec` columns from `SD(diff(x))/sqrt(2)` of the flagged window, so `MDF` and `flux.plot()` use measured noise. Document that datasheet precision underestimates noise on field-worn analyzers (2.5x on the example closure). |
| `precision_mad()` | Maybe | Same helper with `method = "mad"`; document as the robust whole-record option. |
| `flag_detection()`, `mdf_multi()` | **Yes**, partially | `best.flux()` already applies an `MDF` criterion; propose an `MDF.method` argument (`"prec"`, `"two_sigma"`, `"christiansen"`, `"wassmann"`) and a `det_class` column (`emission` / `uptake` / `below detection`). Keep the full precision x method table in fluxqc. |
| `flux_pdf()` | No longer needed | Issue 1 is fixed in 0.4.0; fluxqc keeps `flux_pdf()` for its grouped pages only. |
| `import2RData_flat()` | **Yes** (as a fix) | Issues 2 and 3; the function itself need not go upstream if `import2RData()` is fixed. |
| `import_li7810_subset()` | **Yes** (as an argument) | Issue 4. |
| `windows_from_table()`, `flag_window()` | No longer needed | goFlux 0.4.0's `crop.meas()` / `auto.deadband()` cover cropping; `windows_from_table()` is now a wrapper (plus the Issue 6 fix). |
| `auto_id_rise()`, `find_rise()` | Maybe later | goFlux has `autoID()` (LI-8200) and now `auto.deadband()`; a generic rise-based window finder could sit beside them. |
| `qc_screens()`, `qc_c0_screen()` | Stay | Campaign-level QC policy, not a flux-calculation step. |
| `process_fluxes()`, `write_outputs()` | Stay | Opinionated wrapper around goFlux; goFlux keeps its modular design. |
| `find_clock_offset()` | Stay | Analyzer-independent utility; could be linked from goFlux's documentation. |
| Rescue scripts in `inst/rescue/` (`detect_clock_streams()`, `separate_window_streams()`, `obs_win_deinterleave()`, `click_flux_2d()`, `write_stream_files()`) | **No** | One-off tools for a specific data-damage problem; no longer exported from fluxqc and not upstream candidates. |
| Ebullition add-on (`detect_ebullition()`) | **No PR needed** | goFlux >= 0.3.0 ships `goAquaFlux()` (C. Minaudo); fluxqc's add-on is a thin adapter over it. Two suggestions could go to Minaudo as issues: export `find.bubbles()` (0.4.0 keeps it internal), and offer a de-ebulliated (step-subtracted) diffusive fit as an alternative to the pre-first-bubble window, which leaves `flux_diffusive = NA` when the first bubble comes early. |
| `fluxqc_examples()`, `example_closure()` | Stay (data offer) | goFlux has its own example data; the four Yale-Myers teaching closures could be offered to goFlux as additional example data if wanted. |

## (c) Practical steps

1. Fork `Qepanna/goFlux` to `jongewirtzman/goFlux` and clone it. Add the
   upstream remote (`git remote add upstream https://github.com/Qepanna/goFlux`).
2. Open the open issues (2, 3, 4, 6) on `Qepanna/goFlux` first, so each PR can
   say "Fixes #N" (1 and 5 are already fixed in 0.4.0). Bug-fix PRs first,
   one branch per fix:
   - `fix/crop-meas-max-obs-length-aux` (Issue 6; one line in `R/crop.meas.R`)
   - `fix/import2RData-error-handling` (Issue 3; `R/import2RData.R`)
   - `feat/import2RData-outdir-recursive` (Issue 2)
   - `feat/import-LI7810-dates` (Issue 4; patch rebased onto 0.4.0 in `inst/upstream/patches/`)
3. Feature PRs, each after the corresponding issue/discussion:
   - `feat/click-peak2-stacked` porting `click_peak2_stacked()` into
     `click.peak2()` behind a `gases` argument;
   - `feat/empirical-prec` adding the precision helper and a `prec =
     "empirical"` path;
   - `feat/best-flux-mdf-method` adding `MDF.method` and `det_class` to
     `best.flux()`.
   Each PR: roxygen docs, an entry in `NEWS.md` (goFlux keeps one), a small
   test (goFlux has no `tests/` folder yet, so the first PR would add
   `tests/testthat` with a note in the description), and `R CMD check` clean.
   Keep goFlux's coding style (dplyr pipelines, `.` placeholders) even where
   fluxqc uses base R.
4. Email Karelle Rheault (karh@ign.ku.dk), copying Camille Minaudo for the
   ebullition points, once the issues are open:
   introduce the chamber-flux use case (Yale-Myers, Harvard Forest, BlueFlux,
   Santa Fe, BCI campaigns, all processed with goFlux), point to
   the fluxqc repository and the issues, offer the PRs in the order above,
   and propose co-authoring a short package-update note (e.g. a JOSS
   "software update" or a methods note) covering stacked window selection,
   empirical precision and multi-method detection limits, with her as
   corresponding author of the package. Ask which of the feature PRs she
   would want in goFlux versus linked from its documentation.
5. After merges: bump the goFlux dependency again, remove `click_peak2_stacked()`
   with a deprecation wrapper if it is merged, and update the function map in
   `README.md`. (Already done for 0.4.0: dependency bumped, ebullition and
   window cropping delegated to goFlux.)
