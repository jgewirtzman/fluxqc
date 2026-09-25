# stemflux

Extensions to [goFlux](https://github.com/Qepanna/goFlux) (Rheault et al.
2024) for tree-stem, soil and water chamber fluxes: multi-gas stacked window
selection, non-interactive window detection, empirical per-measurement
precision, minimum detectable flux by several published methods with detection
flags, toggleable physical QC screens, a one-call pipeline with explicit
outputs, a `ggplot2 >= 3.5`-safe PDF exporter, import helpers for very large
LI-7810 files, clock utilities, and an opt-in ebullition add-on.

Everything keeps goFlux's column conventions (`UniqueID`, `start.time`,
`obs.length`, `Vtot`, `Area`, `Pcham`, `Tcham`, `flag`, `Etime`,
`start.time_corr`, `flux.term`, `nb.obs`, `MDF`, `prec`, `best.flux`, ...) so
the functions drop into an existing goFlux pipeline.

## Installation

goFlux is on GitHub, not CRAN, so install it first.

```r
# install.packages("remotes")
remotes::install_github("Qepanna/goFlux")
remotes::install_github("jongewirtzman/stemflux", build_vignettes = TRUE)
```

## Quick start

```r
library(goFlux)
library(stemflux)

d <- example_closure(dead.band = 60)     # one real 600 s UGGA tree-stem closure

res <- process_fluxes(d, gastype = "CH4dry_ppb",                  # defaults: allan / christiansen / 0.95
                      qc = list(c0 = TRUE, convex = TRUE, min_window = list(secs = 120)))
res
#> stemflux result: 1 closures, CH4dry_ppb
#>   MDF_emp: christiansen | allan | 0.95 | detected: 1 | qc_any: 0
res$fluxes[, c("best.flux", "model", "MDF", "sigma_emp", "MDF_emp", "MDF_emp_method", "det_class_emp", "qc_any")]
#>   best.flux model        MDF sigma_emp   MDF_emp              MDF_emp_method det_class_emp qc_any
#> 1  0.255676    HM 0.00113087    2.2253 0.0137646 christiansen | allan | 0.95      emission  FALSE
write_outputs(res, "results/CH4")        # fluxes.csv, settings.json, flux_plots.pdf
```

goFlux's columns (`best.flux`, `MDF`, `quality.check`, ...) come back exactly
as `best.flux()` produced them; stemflux only appends (`sigma_emp`,
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

Vignettes: `vignette("worked-example", "stemflux")` (every intermediate number
of one closure, and how to choose precision and MDF), `vignette("pipeline",
"stemflux")` (the end-to-end pipeline and where every output column comes
from), `vignette("ebullition", "stemflux")` (the opt-in add-on).

## Choosing precision and MDF

`flag_detection()` (called by `process_fluxes()`) takes one precision
estimator, one MDF method and one confidence level (defaults `"allan"`,
`"christiansen"`, `0.95`) and appends `sigma_emp`, `MDF_emp`,
`below_MDF_emp`, `MDF_emp_method`. The multi-method comparison
(`sigma_<precision>`, `MDF_<mdf>_<precision>[_<conf>]`) is appended only
with `extra = list(...)` or `extra = "all"`.

| precision | use when |
|---|---|
| bench sigma (Box 2 of the guidelines paper), passed as `sigma =` (`custom`) | you ran a sealed/standard-gas test on the same analyzer |
| `allan` (per closure, SD(diff x)/sqrt 2) | default; tracks the analyzer's state on the day |
| `mad` (record-wide, MAD(diff x)/sqrt 2) | records with breaths, moves, chamber swaps |
| `rolling` (quietest windows) | comparison only; biased low |
| `datasheet` (goFlux `prec`) | comparability with goFlux only; too optimistic on field analyzers |

| question | MDF method |
|---|---|
| means / budgets / upscaling | keep every flux, flag with `two_sigma` or `christiansen`; never delete |
| is this flux non-zero? | `christiansen` at a stated confidence (per-closure sigma, t-quantile) |
| cross-study comparison | report all four (`goflux`, `wassmann`, `christiansen`, `two_sigma`) |

The `prec` argument of `goFlux::import.*()` feeds only goFlux's own `MDF`
column; stemflux's estimators are applied after `best.flux()`. Nothing is
deleted at any step, only flagged.

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

| Capability | stemflux function | Nearest goFlux function |
|---|---|---|
| **Core path** | `process_fluxes()`, `write_outputs()` | `goFlux()` + `best.flux()` |
| Manual window selection, any gases stacked on one time axis | `click_peak2_stacked(gases = ...)` | `click.peak2()` |
| Non-interactive window from the sustained CO2 rise, notebook fallback, diagnostic PNG per closure | `auto_id_rise()`, `find_rise()` | `autoID()` (LI-8200 only) / `click.peak2()` |
| manID from a table of explicit start/end times | `windows_from_table()`, `flag_window()` | `click.peak2()` |
| Per-measurement precision, first-differenced (Allan, tau = 1 s) | `precision_allan()`, `precision_by_id()` | `prec` argument of `goFlux()` |
| Robust whole-record precision, MAD(dx)/sqrt(2) | `precision_mad()` | -- |
| Quietest-window precision scan | `precision_rolling()` | -- |
| Empirical precision + MDF appended as `sigma_emp`, `MDF_emp`, `below_MDF_emp`, `MDF_emp_method` (comparison columns via `extra`) | `flag_detection()`, `mdf_multi()` | `MDF()`, `best.flux()` (`MDF` criterion) |
| (V/A)(P/RT) with goFlux unit conventions | `flux_term()` | `flux.term()` (unexported) |
| Toggleable physical QC screens (C0, CO2 tracer, convex, min window, ambient start, leak), flags only | `qc_screens()` (`qc_c0_screen()` wrapper) | -- |
| Multi-page PDF of `flux.plot()` output, ggplot2 >= 3.5 safe | `flux_pdf()` | `flux2pdf()` |
| Date subset of huge LI-7810 `.data` files (shim: forwards to `import.LI7810(dates=)` once the upstream patch in `inst/upstream/patches/` is merged) | `import_li7810_subset()` | `import.LI7810()` |
| Recursive per-file import without the `muffleError` wrapper or `getwd()/RData` writes | `import2RData_flat()` | `import2RData()` |
| Utility: notebook-to-analyzer clock offset by cross-correlating notebook starts with rise onsets | `find_clock_offset()` | -- |
| **Ebullition add-on (opt-in)**: IQR (Hoffmann 2017) / threshold / running-variance event detection, de-ebulliated goFlux input, ebullitive flux | `detect_ebullition()`, `split_ebullition()`, `merge_ebullition()`; `process_fluxes(ebullition = "flag"|"split")` | -- (goFlux has none) |
| Six bundled closures with real geometry (see Data) | `stemflux_examples()`, `example_closure()` | `manID.UGGA` (example data) |
| One-off rescue scripts for LGR clock resets and interleaved streams (unexported: `detect_clock_streams()`, `write_stream_files()`, `separate_window_streams()`, `obs_win_deinterleave()`, `click_flux_2d()`) | `system.file("rescue", package = "stemflux")` | -- |

## Data

`stemflux_examples()` returns six bundled closures in goFlux shape with an
auxfile holding the real geometry of each (`example_closure()` is the alias
for `a`):

| id | UniqueID | shows | kind | geometry source |
|---|---|---|---|---|
| a | `2021-05-18_5086_I_7.0` | high-signal valid closure (CO2 and CH4 rise) | real | Yale-Myers campaign auxfile |
| b | `2020-09-23_5092_I_6` | valid closure, CH4 below detection | real | same |
| c | `2021-02-25_1641_WD_5` | failed closure (flat CO2 and CH4) | real | same |
| d | `2020-09-21_3589_U_6` | ambiguous: CH4 rise over flat CO2 | real | same |
| e | `BF_CP40_2022-10-23_P06` | floating chamber on water with 5 ebullition events | real | BlueFlux chamber parameters (Area 324.3 cm2, Vtot 4.318 L); Tcham = site-day mean air temperature |
| f | `SYN_uptake_from_a` | synthetic CH4 uptake built from the measured noise of a | synthetic | copied from a |

The real closures are released with the package under its MIT license with
the data owner's permission. `a`-`d`: Yale-Myers Forest (CT, USA) tree-stem
record 2020-2021, ABB/LGR UGGA, semi-rigid sleeve chambers; cite as
"Gewirtzman et al., Yale-Myers Forest stem-flux record 2020-2021, distributed
with the stemflux R package". `e`: NASA Carbon Monitoring System BlueFlux
field campaign, ground component, site CP40 (Christian Point, ghost-forest
mangrove, south Florida), ABB/LGR GLA131; cite Poulter et al. (2023),
*Environmental Research Letters* 18, 075009,
https://doi.org/10.1088/1748-9326/acdae6 (blueflux-ground repository, MIT).

## Contributing upstream

`inst/upstream/UPSTREAM.md` holds ready-to-paste goFlux issue reports (five
bugs found while writing these tools) and a plan for which functions to
propose as goFlux pull requests.

## Citation

See `citation("stemflux")`, `CITATION.cff` and `.zenodo.json`. goFlux itself
should be cited as Rheault, K., Christiansen, J. R., & Larsen, K. S. (2024).
goFlux: A user-friendly way to calculate GHG fluxes yourself, regardless of
user experience. *Journal of Open Source Software*, 9(96), 6393.
https://doi.org/10.21105/joss.06393

## License

MIT, Jonathan Gewirtzman.
