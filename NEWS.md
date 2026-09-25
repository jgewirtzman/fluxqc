# stemflux 0.2.0

Revision after author review.

* `click_peak2_stacked(gases = ...)` stacks any number and order of gas
  columns (`"auto"` = every recognised goFlux gas present); `auto_id_rise()`
  gained `detect.gas` and `gases` for its diagnostic plots. Replaces the
  `gastype` / `gastype2` arguments.
* `flag_detection(flux.result, traces, precision = "mad", mdf = "wassmann",
  conf = 0.95, group = NULL, extra = NULL)` (defaults = the companion
  manuscript's precision and MDF form, record-wide MAD sigma per group, at
  the stated 95 percent two-sided confidence, z = qnorm(0.975) = 1.96 (2 sigma);
  quantiles are two-sided) appends exactly `sigma_emp`,
  `MDF_emp`, `below_MDF_emp`, `MDF_emp_method` (plus `detected_emp`,
  `det_class_emp`). The multi-method comparison (`sigma_<precision>`,
  `MDF_<mdf>_<precision>[_<conf>]`) is appended only with `extra = list(...)`
  or `extra = "all"`. goFlux's columns (`MDF`, `best.flux`, ...) are never
  modified or duplicated; a test checks that `goFlux::flux.plot()` runs
  unchanged on `process_fluxes()$fluxes`.
* `detect_ebullition(merge.consecutive = TRUE)`: a run of consecutive flagged
  samples is one event (bubble bursts span several 1 Hz samples).
* `stemflux_examples()`: ten anonymized real closures across campaigns and
  two analyzer families (ABB/LGR UGGA/GLA131, LI-COR LI-7810), each the
  campaign's clicked goFlux window with obs.win context and real geometry,
  chosen by audited criteria (`inst/extdata/build_examples.R`, one BlueFlux ebullition closure, one synthetic
  uptake), released with the package; `example_closure()` is the alias for
  `a`. Data files `inst/extdata/stemflux_examples*.csv` replace
  `example_closure.csv`.
* `detect_clock_streams()` un-exported (moved to `inst/rescue/`): it depends
  on the raw LGR file layout and a one-campaign failure mode.
* `import_li7810_subset()` is now a shim: forwards to `import.LI7810(dates=)`
  when the installed goFlux has that argument (upstream patch in
  `inst/upstream/patches/import.LI7810-dates.patch`), else subsets locally.
* Provenance statement for the ebullition add-on (no code copied from the
  Hoffmann et al. 2017 script or fluxCH4) in the vignette and `inst/CITATION`;
  "Analyzer independence" section in the README.
* New `qc_screens()` with toggleable screens `c0`, `co2_tracer`, `convex`,
  `min_window`, `ambient_start`, `leak`, each a `qc_<screen>` flag column plus
  `qc_any` and `qc_note`; rows are never deleted. `qc_c0_screen()` is a thin
  wrapper.
* New `process_fluxes()` (goFlux -> best.flux -> flag_detection -> qc_screens
  in one call, `settings` recorded) and `write_outputs()` (`fluxes.csv`,
  `settings.json`, `flux_plots.pdf`, and `ebullition_events.csv` when the
  add-on ran).
* Ebullition is now an opt-in add-on: `detect_ebullition()` with
  `method = "iqr"` (Hoffmann et al. 2017 quartile fences, default),
  `"threshold"` and `"runvar"`, working per `UniqueID` on flagged data;
  `split_ebullition()` returns a de-ebulliated goFlux input plus the
  ebullitive flux; `merge_ebullition()` appends totals;
  `process_fluxes(ebullition = "flag" | "split", ebullition_ids = ...)`
  dispatches to it. Documented in `vignette("ebullition")` with a survey of
  how other packages handle ebullition.
* New `find_clock_offset()`: constant notebook-to-analyzer clock offset by
  cross-correlating notebook start times with concentration-rise onsets.
* `separate_window_streams()`, `obs_win_deinterleave()`, `click_flux_2d()`
  and `write_stream_files()` are no longer exported; they live as plain
  scripts in `inst/rescue/` with a README. `detect_clock_streams()` stays.
* `UPSTREAM.md`: fifth goFlux issue (namespace-qualified `import.*()` calls
  fail); rescue tools removed from the PR candidates.

# stemflux 0.1.0

First release. Consolidates helper code written between 2023 and 2026 while
processing chamber fluxes with goFlux across the Yale-Myers, Harvard Forest,
BlueFlux (Everglades), Santa Fe and Barro Colorado campaigns.

* Window selection: `click_peak2_stacked()`, `auto_id_rise()`, `find_rise()`,
  `windows_from_table()`, `flag_window()`.
* Precision: `precision_allan()`, `precision_mad()`, `precision_rolling()`,
  `precision_by_id()`.
* Detection limits: `mdf_multi()`, `flag_detection()`, `flux_term()`.
* QC and plots: `qc_c0_screen()`, `flux_pdf()` (works with ggplot2 >= 3.5).
* Clock resets and interleaved streams: `detect_clock_streams()`,
  `write_stream_files()`, `separate_window_streams()`,
  `obs_win_deinterleave()`, `click_flux_2d()`.
* Ebullition: `detect_ebullition()`, `split_ebullition()`.
* Import: `import_li7810_subset()`, `import2RData_flat()`.
* Data: `example_closure()` (one real UGGA tree-stem closure with auxfile
  geometry).
* Vignettes: "From raw trace to reported flux" (fully evaluated) and "The
  recommended stem-flux pipeline".
* `inst/upstream/UPSTREAM.md`: goFlux bug reports and a contribution plan.
