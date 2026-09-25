# stemflux

Extensions to [goFlux](https://github.com/Qepanna/goFlux) (Rheault et al.
2024) for tree-stem, soil and water chamber fluxes: stacked-panel window
selection, non-interactive window detection, empirical per-measurement
precision, multi-method minimum detectable flux with detection flags, a
starting-concentration contamination screen, a `ggplot2 >= 3.5`-safe PDF
exporter, analyzer clock-reset handling, ebullition detection/partitioning,
and import helpers for very large LI-7810 files.

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

# one real closure ships with the package (600 s, UGGA, tree stem)
d <- example_closure(dead.band = 60)

CH4 <- best.flux(goFlux(d, "CH4dry_ppb"))
CO2 <- best.flux(goFlux(d, "CO2dry_ppm"))

sig <- precision_by_id(d, "CH4dry_ppb")     # sigma = SD(diff x) / sqrt(2), per closure
CH4 <- flag_detection(CH4, sigma = sig)     # MDF by 8 methods + below_* flags + det_class
CH4 <- qc_c0_screen(CH4, k = 1.5)           # chamber at ambient at closure?

CH4[, c("best.flux", "model", "sigma", "MDF_goflux", "MDF_2sigma", "det_class", "qc_c0_pass")]
#>   best.flux model  sigma  MDF_goflux MDF_2sigma det_class qc_c0_pass
#> 1  0.255676    HM 2.2253 0.00112621 0.00556930  emission       TRUE
```

A full campaign:

```r
licor <- import_li7810_subset("TG10-01503-2026-05-01T000000.data",
                              dates = c("2026-05-22", "2026-05-23"), timezone = "UTC")
ow    <- obs.win(licor, auxfile = aux, obs.length = 600, shoulder = 120)
manID <- auto_id_rise(ow, plot.dir = "plots/auto_id")          # or click_peak2_stacked(ow)
CH4   <- flag_detection(best.flux(goFlux(manID, "CH4dry_ppb")),
                        sigma = precision_by_id(manID, "CH4dry_ppb"))
flux_pdf(flux.plot(CH4, manID, "CH4dry_ppb"), "plots/CH4.pdf")
```

See `vignette("worked-example", "stemflux")` for every intermediate number of
one closure and `vignette("pipeline", "stemflux")` for the recommended
end-to-end pipeline.

## Function map

| Capability | stemflux function | Nearest goFlux function |
|---|---|---|
| Manual window selection with CO2 and CH4 stacked on one time axis | `click_peak2_stacked()` | `click.peak2()` |
| Non-interactive window from the sustained CO2 rise, notebook fallback, diagnostic PNG per closure | `auto_id_rise()`, `find_rise()` | `autoID()` (LI-8200 only) / `click.peak2()` |
| manID from a table of explicit start/end times | `windows_from_table()`, `flag_window()` | `click.peak2()` |
| Per-measurement precision, first-differenced (Allan, tau = 1 s) | `precision_allan()`, `precision_by_id()` | `prec` argument of `goFlux()` |
| Robust whole-record precision, MAD(dx)/sqrt(2) | `precision_mad()` | -- |
| Quietest-window precision scan of a continuous record | `precision_rolling()` | -- |
| MDF by goFlux / Wassmann 2018 / Christiansen 2015 / 2-sigma methods | `mdf_multi()` | `MDF()` |
| Append MDFs, below-detection flags and emission/uptake/below class to `best.flux()` output | `flag_detection()` | `best.flux()` (`MDF` criterion) |
| (V/A)(P/RT) with goFlux unit conventions | `flux_term()` | `flux.term()` (unexported) |
| Chamber contamination screen on starting concentration | `qc_c0_screen()` | -- |
| Multi-page PDF of `flux.plot()` output, ggplot2 >= 3.5 safe | `flux_pdf()` | `flux2pdf()` |
| Detect analyzer clock resets from raw write order | `detect_clock_streams()`, `write_stream_files()` | -- |
| De-interleave two time streams per observation window | `separate_window_streams()`, `obs_win_deinterleave()` | `obs.win()` |
| 2-D (time x concentration) box selection with immediate flux | `click_flux_2d()` | `click.peak2()` + `goFlux()` |
| Ebullition jump detection, de-ebulliated trace | `detect_ebullition()` | -- |
| Diffusive / ebullitive partition | `split_ebullition()` | -- |
| awk/date subset of huge LI-7810 `.data` files before import | `import_li7810_subset()` | `import.LI7810()` |
| Recursive per-file import without the `muffleError` progress-bar wrapper or `getwd()/RData` writes | `import2RData_flat()` | `import2RData()` |
| One real closure with auxfile geometry | `example_closure()` | `manID.UGGA` (example data) |

## Contributing upstream

`inst/upstream/UPSTREAM.md` holds ready-to-paste goFlux issue reports (four
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
