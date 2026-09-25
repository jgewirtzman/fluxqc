# Rescue scripts (not exported)

These functions were one-off rescue tools written for the BlueFlux
(Everglades) ground campaign, where three ABB/LGR GLA131 analyzers produced
files with mid-file clock resets. After a reset, rows written later carry
earlier time stamps, so a time-sorted import interleaves two streams and the
observation windows contain two overlapping traces. They were moved out of
the stemflux namespace in 0.2.0 because they solve a specific data-damage
problem and are not part of a normal pipeline. They are kept here, tested by
hand, for anyone who meets the same problem.

Source them with

```r
source(system.file("rescue", "rescue_tools.R", package = "stemflux"))
```

| Function | Purpose |
|---|---|
| `write_stream_files(streams, filename, dir)` | write the "ahead" and "behind" streams found by `stemflux::detect_clock_streams()` to separate files so each can be imported with `goFlux::import.UGGA()` |
| `separate_window_streams(window_data, gas_col, threshold)` | label two interleaved streams within one short window by switching on jumps larger than `threshold` |
| `obs_win_deinterleave(raw_data, auxfile, ...)` | `obs.win()` replacement that keeps the smoother stream per window and adds `Etime` and chamber columns |
| `click_flux_2d(raw_data, auxfile, ...)` | interactive: draw a time x concentration box around the good points and run `goFlux::goFlux()` on them at once |

The exported `stemflux::detect_clock_streams()` (reset detection from raw
write order) and `stemflux::find_clock_offset()` (notebook-to-analyzer clock
offset) are the general-purpose parts of this work.
