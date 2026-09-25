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
