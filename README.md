# fluxqc

Extensions to [goFlux](https://github.com/Qepanna/goFlux) (Rheault et al.
2024) for chamber fluxes from any surface (soil, water, tree stems, leaves;
diffusive or ebullitive): multi-gas stacked window
selection, non-interactive window detection, empirical per-measurement
precision, minimum detectable flux by several published methods with detection
flags, toggleable physical QC screens, a one-call pipeline with explicit
outputs, a `ggplot2 >= 3.5`-safe PDF exporter, import helpers for very large
LI-7810 files, clock utilities, and an opt-in ebullition add-on.

## Status: staging ground for goFlux

fluxqc is not meant to live for ever as a separate package. Everything in it
that is general is being folded into goFlux itself, developed as feature
branches on a local clone of [Qepanna/goFlux](https://github.com/Qepanna/goFlux)
(each with a `PR_<branch>.md` describing the change) and submitted upstream
as pull requests. As each one is merged and released, the corresponding
fluxqc function is retired in favour of the goFlux one:

| goFlux branch | goFlux function(s) | replaces in fluxqc |
|---|---|---|
| `feat/empirical-precision-mdf` | `empirical.prec()`, `conf =` in `MDF()`, `goFlux()`, `best.flux()` | `precision_mad()`, `precision_allan()`, `precision_by_id()`, the `wassmann` MDF in `flag_detection()` |
| `feat/qc-flags` | `qc.flags()`, `co2.tracer()` | `qc_screens()`, `qc_c0_screen()` |
| `feat/click-peak2-stacked` | `click.peak2(gases = )` | `click_peak2_stacked()` |
| `feat/import-li7810-dates` | `import.LI7810(dates = )` | `import_li7810_subset()` shim |
| `fix/crop-meas-aux-obslength` | `crop.meas(max.obs.length = "aux")` | per-closure workaround in `windows_from_table()` |
| `feat/aqua-diffusive-deebulliated` | `goAquaFlux(diffusion.window = "deebulliated")` | nothing (new capability; `detect_ebullition()` passes it through once released) |

What stays here for now: `process_fluxes()`/`write_outputs()` (the one-call
pipeline), `mdf_multi()` (method comparison), `flux_pdf()`, `find_clock_offset()`,
`auto_id_rise()` (rise detection is a different approach from goFlux's
`autoID()` window construction and does not fit it cleanly), the bundled
example closures and the rescue tools in `inst/rescue/`.

Everything keeps goFlux's column conventions (`UniqueID`, `start.time`,
`obs.length`, `Vtot`, `Area`, `Pcham`, `Tcham`, `flag`, `Etime`,
`start.time_corr`, `flux.term`, `nb.obs`, `MDF`, `prec`, `best.flux`, ...) so
the functions drop into an existing goFlux pipeline.

## Installation

goFlux is on GitHub, not CRAN, so install it first.

```r
# install.packages("remotes")
remotes::install_github("Qepanna/goFlux")   # >= 0.4.0
remotes::install_github("jongewirtzman/fluxqc", build_vignettes = TRUE)
```

## Quick start

```r
library(goFlux)
library(fluxqc)

d <- example_closure()                   # one real UGGA closure (a tree stem) with its clicked window

res <- process_fluxes(d, gastype = "CH4dry_ppb",                  # defaults: mad / wassmann / 0.95
                      qc = list(c0 = TRUE, convex = TRUE, min_window = list(secs = 120)))
res
#> fluxqc result: 1 closures, CH4dry_ppb
#>   MDF_emp: wassmann | mad | 0.95 | detected: 1 | qc_any: 0
res$fluxes[, c("best.flux", "model", "MDF", "sigma_emp", "MDF_emp", "MDF_emp_method", "det_class_emp", "qc_any")]
#>   best.flux model         MDF sigma_emp    MDF_emp        MDF_emp_method det_class_emp qc_any
#> 1  0.253722    HM 0.000826   2.2069 0.00378 wassmann | mad | 0.95      emission  FALSE
write_outputs(res, "results/CH4")        # fluxes.csv, settings.json, flux_plots.pdf
```

goFlux's columns (`best.flux`, `MDF`, `quality.check`, ...) come back exactly
as `best.flux()` produced them; fluxqc only appends (`sigma_emp`,
`MDF_emp`, `below_MDF_emp`, `MDF_emp_method`, `detected_emp`,
`det_class_emp`, `qc_*`). Multi-method comparison columns are available with
`extra = list(precision = ..., mdf = ..., conf = ...)` or `extra = "all"`.

```r
```

A full campaign:

```r
licor <- import_li7810_subset("TG10-01503-2026-05-01T000000.data",
                              dates = c("2026-05-22", "2026-05-23"), timezone = "UTC")
ow    <- obs.win(licor, auxfile = aux, obs.length = 600, shoulder = 120)
manID <- auto_id_rise(ow, plot.dir = "plots/auto_id")     # or click_peak2_stacked(ow, gases = "auto")
CH4   <- process_fluxes(manID, aux = aux, gastype = "CH4dry_ppb", group = "campaign")
write_outputs(CH4, "results/CH4")
```

Vignettes: `vignette("worked-example", "fluxqc")` (every intermediate number
of one closure, and how to choose precision and MDF), `vignette("pipeline",
"fluxqc")` (the end-to-end pipeline and where every output column comes
from), `vignette("ebullition", "fluxqc")` (the opt-in add-on).

## Choosing precision and MDF

`flag_detection()` (called by `process_fluxes()`) takes one precision
estimator, one MDF method and one two-sided confidence level and appends
`sigma_emp`, `MDF_emp`, `below_MDF_emp`, `MDF_emp_method`. The default is the
companion data-processing manuscript's method: `precision = "mad"` (sigma =
MAD(dx)/sqrt 2 over the entire 1 Hz record of each field period, `group`),
`mdf = "wassmann"` (MDF = z sigma / t flux.term) and `conf = 0.95` two-sided
(z = qnorm(0.975) = 1.96, i.e. the 2-sigma criterion; the companion
manuscript itself used `conf = 0.90`, z = 1.645).

| precision | formula | over |
|---|---|---|
| `mad` (default) | MAD(dx) x 1.4826 / sqrt 2 | whole record per `group` |
| `allan` | SD(dx) / sqrt 2 | flagged rows per closure |
| `rolling` | quietest-window residual SD | flagged rows per closure (comparison only) |
| `datasheet` | goFlux `prec` | -- |
| `custom` | user `sigma` (bench test, Box 2) | -- |

| MDF (`dC_min / t x flux.term`) | `dC_min` | origin |
|---|---|---|
| `wassmann` (default) | z sigma_global, z two-sided | Wassmann et al. 2018 |
| `goflux` | `prec` | goFlux `MDF()`, which cites Christiansen et al. 2015 |
| `christiansen` | sigma_i x 3 x t_crit(nb.obs - 2) | as implemented in the author's scripts; the x3 is their safety factor |
| `two_sigma` | 2 sigma_i | guidelines paper Table 2.5 illustration |

The per-closure sigma is a QC diagnostic, not an MDF basis: `qc_screens(noisy
= list(k = 1.5))` flags closures whose own MAD sigma exceeds k x the campaign
sigma. Comparison columns (`sigma_<precision>`, `MDF_<mdf>_<precision>[_<conf>]`)
are appended only with `extra = list(...)` or `extra = "all"`. The `prec`
argument of `goFlux::import.*()` feeds only goFlux's own `MDF` column;
fluxqc's estimators are applied after `best.flux()`. Nothing is deleted at
any step, only flagged.

## Analyzer independence

The core functions consume goFlux's standard columns (`UniqueID`,
`POSIX.time`, `flag`, `Etime`, the gas columns `CO2dry_ppm`, `CH4dry_ppb`,
`N2Odry_ppb`, `H2O_ppm`, ..., and the `best.flux()` output), so they work with
any goFlux importer: UGGA/LGR, LI-7810/7820/8x00, Picarro G2xxx/G4301, GAIA,
and the rest. Gas columns are always named through the `gastype` / `gases` /
`detect.gas` arguments; the only goFlux column name used implicitly is
`CO2dry_ppm` in the `co2_tracer` screen, because that screen *is* about CO2.
The analyzer-specific pieces are `import_li7810_subset()` (LI-7810 `.data`
files), the shipped example closures (UGGA / GLA131) and the unexported
rescue scripts (raw LGR text-file layout). This was verified by grepping the
core sources for gas-column and file-format literals.

## Function map

| Capability | fluxqc function | Nearest goFlux function |
|---|---|---|
| **Core path** | `process_fluxes()`, `write_outputs()` | `goFlux()` + `best.flux()` |
| Manual window selection, any gases stacked on one time axis | `click_peak2_stacked(gases = ...)` | `click.peak2()` |
| Non-interactive window from the sustained CO2 rise, notebook fallback, diagnostic PNG per closure | `auto_id_rise()`, `find_rise()` | `autoID()` (LI-8200 only) / `click.peak2()` |
| manID from a table of explicit start/end times (wrapper over `crop.meas()`) | `windows_from_table()`, `flag_window()` | `crop.meas()`, `auto.deadband()` (goFlux >= 0.4.0) |
| Per-measurement precision, first-differenced (Allan, tau = 1 s) | `precision_allan()`, `precision_by_id()` | `prec` argument of `goFlux()` |
| Robust whole-record precision, MAD(dx)/sqrt(2) | `precision_mad()` | -- |
| Quietest-window precision scan | `precision_rolling()` | -- |
| Empirical precision + MDF appended as `sigma_emp`, `MDF_emp`, `below_MDF_emp`, `MDF_emp_method` (comparison columns via `extra`) | `flag_detection()`, `mdf_multi()` | `MDF()`, `best.flux()` (`MDF` criterion) |
| (V/A)(P/RT) with goFlux unit conventions | `flux_term()` | `flux.term()` (unexported) |
| Toggleable physical QC screens (C0, CO2 tracer, convex, min window, ambient start, noisy closure, leak), flags only; switch `co2_tracer` off for surfaces without respiratory CO2 | `qc_screens()` (`qc_c0_screen()` wrapper) | -- |
| Multi-page PDF of `flux.plot()` output, ggplot2 >= 3.5 safe | `flux_pdf()` | `flux2pdf()` |
| Date subset of huge LI-7810 `.data` files (shim: forwards to `import.LI7810(dates=)` once the upstream patch in `inst/upstream/patches/` is merged) | `import_li7810_subset()` | `import.LI7810()` |
| Recursive per-file import without the `muffleError` wrapper or `getwd()/RData` writes | `import2RData_flat()` | `import2RData()` |
| Utility: notebook-to-analyzer clock offset by cross-correlating notebook starts with rise onsets | `find_clock_offset()` | -- |
| **Ebullition add-on (opt-in)**: thin adapter over goFlux's bubble detection and diffusion/ebullition partition, appended to the flux table | `detect_ebullition()`; `process_fluxes(ebullition = "flag"|"split")` | `goAquaFlux()` (goFlux >= 0.3.0, C. Minaudo) |
| Six bundled closures with real geometry (see Data) | `fluxqc_examples()`, `example_closure()` | `manID.UGGA` (example data) |
| One-off rescue scripts for LGR clock resets and interleaved streams (unexported: `detect_clock_streams()`, `write_stream_files()`, `separate_window_streams()`, `obs_win_deinterleave()`, `click_flux_2d()`) | `system.file("rescue", package = "fluxqc")` | -- |

## Data

`fluxqc_examples()` returns a curated, anonymized, cross-campaign set of
ten real closures in goFlux shape (the campaign's clicked window with its
obs.win context) plus an auxfile with the real geometry; `example_closure()`
is the alias for `emission_stem_semirigid`. Identifiers, dates, tree tags and
plot codes are removed (the time axis is a synthetic calendar preserving every
interval); analyzer, chamber family, geometry, season and year are kept.
`inst/extdata/build_examples.R` documents exactly how each was pulled, and
`fluxqc_examples_diagnostics.csv` holds the audit numbers, including a
continuity check of every clicked window (no timestamp gaps; single-sample
CH4 steps within the noise; a Picarro G4301 closure was rejected on it, so
that format is not represented).

| key | shows | campaign / season | analyzer | CO2 slope (ppm/s) | CH4 flux, MDF95 (nmol m-2 s-1) |
|---|---|---|---|---|---|
| `emission_stem_semirigid` | clean emission, the worked-example closure | upland stem, monthly survey, spring 2021 | UGGA 1 Hz | 0.82 | 0.254, 0.0031 |
| `emission_upland_stem` | clean mid-range upland emission | upland stem, height campaign, summer 2021 | UGGA 5 s | 0.68 | 0.171, 0.018 |
| `below_detection_stem` | CO2 rising, CH4 below detection | upland stem, height campaign, summer 2021 | UGGA 5 s | 0.63 | 0.0018, 0.019 |
| `failed_closure_stem` | CO2 not rising, CH4 below detection | upland stem, monthly survey, winter 2021 | UGGA 1 Hz | -0.025 | 0.0018, 0.0075 |
| `ambiguous_stem` | CH4 marginally detected over a flat CO2 tracer | upland stem, monthly survey, winter 2020 | UGGA 1 Hz | -0.0003 (p 0.69) | -0.0046, 0.0046 |
| `uptake_stem` | strongest coherent stem uptake, CO2 rising | upland stem, monthly survey, summer 2020 | UGGA 1 Hz | 0.73 | -0.176, 0.0078 |
| `uptake_soil` | strong clean soil uptake | forest soil collar, summer 2020 | UGGA 1 Hz | 0.61 | -2.25, 0.060 |
| `high_flux_wetland_stem` | large diffusive mangrove-stem emission | coastal wetland stem, wet season 2022 | GLA131 1 Hz | 0.18 | 29.3, 0.021 |
| `ebullition_floating` | one bubble burst inside a sealed floating-chamber closure | open water, ghost-forest mangrove site, wet season 2022 | GLA131 1 Hz | 0.11 | -- |
| `emission_li7810_stem` | clean emission, second analyzer format | tree stem, Santa Fe 2026, spring | LI-7810 1 Hz | 0.61 | 0.653, 0.0028 |

The real closures are released with the package under its MIT license with
the data owner's permission. Upland stem and soil closures: Yale-Myers Forest
(CT, USA) record 2020-2021 (Gewirtzman et al.; cite as "Gewirtzman et al.,
Yale-Myers Forest chamber-flux record 2020-2021, distributed with the
fluxqc R package"). Wetland closures: NASA Carbon Monitoring System BlueFlux
field campaign, ground component (Poulter et al. 2023, *Environmental Research
Letters* 18, 075009, https://doi.org/10.1088/1748-9326/acdae6;
blueflux-ground repository, MIT). LI-7810 closure: the author's Santa Fe 2026
campaign.

## Contributing upstream

`inst/upstream/UPSTREAM.md` holds ready-to-paste goFlux issue reports (five
bugs found while writing these tools) and a plan for which functions to
propose as goFlux pull requests.

## Citation

See `citation("fluxqc")`, `CITATION.cff` and `.zenodo.json`. goFlux itself
should be cited as Rheault, K., Christiansen, J. R., & Larsen, K. S. (2024).
goFlux: A user-friendly way to calculate GHG fluxes yourself, regardless of
user experience. *Journal of Open Source Software*, 9(96), 6393.
https://doi.org/10.21105/joss.06393

## License

MIT, Jonathan Gewirtzman.
