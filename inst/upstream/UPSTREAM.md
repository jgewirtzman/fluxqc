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
Steenberg Larsen.

## (a) Bug reports, ready to paste as GitHub issues

### Issue 1: `flux2pdf()` class check fails with ggplot2 >= 3.5 (and is now disabled)

**Title:** flux2pdf(): plot.list validation breaks with ggplot2 >= 3.5; check
is commented out in v0.2.0

**Description**

`flux2pdf()` originally validated its input with

```r
if (length(grep("ggplot", sapply(plot.list, class)[1, ])) != length(plot.list))
  stop("all elements of 'plot.list' must be of class 'gg, ggplot'")
```

`sapply(plot.list, class)` returns a character *matrix* only while every
element has exactly the same number of classes, which was true when every
ggplot object was `c("gg", "ggplot")`. The class vector changed in ggplot2 3.5
and again in the S7-based 4.0 series, where it has five elements
(`ggplot2::ggplot`, `ggplot`, `ggplot2::gg`, `S7_object`, `gg`). Any list in
which elements differ in class length (a `patchwork` or a plain `gg` mixed in,
or an accidental non-plot) makes `sapply()` return a *list*, and `[1, ]` then
errors with "incorrect number of dimensions" instead of the intended message.
In v0.2.0 the check is commented out entirely ("Re-occuring error even when
all elements are ggplot objects..."), so a non-ggplot element now fails later
inside `marrangeGrob()` with an unhelpful message. Downstream code that copied
the original check (our `07_quality_plots.R`) broke outright when ggplot2 was
upgraded.

**Minimal reproducible example**

```r
library(ggplot2)   # 4.0.2
p <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point()
class(p)
#> [1] "ggplot2::ggplot" "ggplot"  "ggplot2::gg"  "S7_object"  "gg"
sapply(list(p, p), class)[1, ]   # a matrix; happens to contain "ggplot"
sapply(list(p, 1), class)[1, ]   # Error in ...[1, ]: incorrect number of dimensions
```

**Proposed fix** (one line, no new dependency):

```r
if (!all(vapply(plot.list, inherits, logical(1), what = "ggplot")))
  stop("all elements of 'plot.list' must be ggplot objects")
```

fluxqc's `flux_pdf()` implements this and also reads the `UniqueID` from
`plot_env` defensively, falling back to the plot index.

### Issue 2: `import2RData()` writes into `getwd()/RData` and (for some instruments) does not search subfolders

**Title:** import2RData(): output is forced into getwd()/RData; UGGA/LI-COR
file search is not consistent across instruments

**Description**

1. Every `import.*()` function with `save = TRUE` does
   `RData_folder <- paste(getwd(), "RData", sep = "/")` and creates that
   folder in the *current working directory*, regardless of `path`. There is
   no argument to choose the output folder, so a project that keeps raw data
   on one drive and processed data on another has to `setwd()` around the
   call, and running the import from a package or an `Rscript` in a different
   directory scatters `RData` folders. Suggest an `outdir = file.path(path,
   "RData")` argument (default keeps the current behaviour) or, better, make
   `save = FALSE` the default and return the merged data.
2. The file search differs by instrument: e.g. DX4015 uses
   `list.files(path, pattern = "\\.TXT", recursive = TRUE)` while the last
   branch (`\\.csv`) is non-recursive. Field campaigns store raw files in
   date-named subfolders (`LGR1/2023-06-29/micro_..._f0001.txt`), and the
   UGGA/LGR folders also contain `_b`, `_l`, `_p` and zero-byte files that
   match the pattern but are not concentration data; these currently have to
   be staged into a flat folder by hand (see `01_prepare_raw_data.R` in the
   whole_tree_flux project). Suggest `recursive = TRUE` everywhere plus an
   optional `pattern` argument and skipping zero-byte files.

**Minimal reproducible example**

```r
library(goFlux)
tmp <- tempfile(); dir.create(tmp)
setwd(tmp)
import2RData(path = system.file("extdata", "LI7810", package = "goFlux"),
             instrument = "LI7810", date.format = "ymd", timezone = "UTC")
list.files(tmp)      # "RData" appears here, not next to the data
```

fluxqc's `import2RData_flat()` shows the alternative: recursive search with
a pattern, zero-byte filtering, per-file import with `save = FALSE`, nothing
written to disk.

### Issue 3: `import2RData()` dies with "no restart 'muffleError' found"

**Title:** import2RData(): withCallingHandlers() calls a non-existent
'muffleError' restart, turning any per-file error into an uninformative
failure

**Description**

`import2RData()` wraps each `import.*()` call in

```r
withCallingHandlers(import.XXX(...),
  error = function(e) { errs <<- c(errs, conditionMessage(e)); invokeRestart("muffleError") },
  warning = function(w) { ...; invokeRestart("muffleWarning") },
  message = function(m) { ...; invokeRestart("muffleMessage") })
```

`muffleWarning` and `muffleMessage` are real restarts established by
`warning()` and `message()`, but there is no `muffleError` restart in R. As
soon as any single file throws (a truncated file, an unexpected column, a
non-data `_b`/`_l` file matched by the pattern), the handler runs
`invokeRestart("muffleError")`, which itself errors with
`no restart 'muffleError' found`, the progress bar stops, the original error
message is lost, and nothing is returned for the files that *did* import.
Observed with UGGA files from Yale-Myers Forest (whole_tree_flux
`ymf_02_import.R` documents the workaround: call `import.UGGA()` file by
file).

**Minimal reproducible example**

```r
library(goFlux)
tmp <- tempfile(); dir.create(tmp)
writeLines("this is not an analyzer file", file.path(tmp, "bad_f0001.txt"))
import2RData(path = tmp, instrument = "UGGA", date.format = "mdy", timezone = "UTC")
#> Error in invokeRestart("muffleError") : no restart 'muffleError' found
```

**Proposed fix:** use `tryCatch()` for the error branch (there is nothing to
"continue" from after an error) and keep `withCallingHandlers()` only for
warnings and messages:

```r
res <- tryCatch(
  withCallingHandlers(import.XXX(...),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }),
  error = function(e) { errs <<- c(errs, paste(basename(file_list[i]), conditionMessage(e))); NULL })
```

and report `errs` after the progress bar, as the code already intends.

### Issue 4: `import.LI7810()` has no way to read a date subset of a huge `.data` file

**Title:** import.LI7810(): add a `dates` (or `start`/`end`) argument to skip
rows before parsing; campaign-long files are 100+ MB

**Description**

The LI-7810 writes one continuous tab-delimited `.data` file per power cycle;
over a multi-week campaign this reaches 150-200 MB (e.g.
`TG10-01503-2026-05-01T000000.data`, 163 MB, ~1.5 million rows at 1 Hz).
`import.LI7810()` reads the whole file with `read.delim()` and then builds
POSIX times for every row, which takes minutes and several GB of memory when
only a handful of field days are needed. The `DATE` column (7th field of
`DATA` rows) makes a cheap pre-filter possible.

**Minimal reproducible example** (timing on the goFlux example file scaled up)

```r
f <- system.file("extdata", "LI7810", "LI7810.data", package = "goFlux")
big <- tempfile(fileext = ".data")
l <- readLines(f); writeLines(c(l[1:7], rep(l[-(1:7)], 2000)), big)   # ~ 100 MB
system.time(import.LI7810(big, timezone = "UTC"))
```

**Proposed fix:** an optional `dates = NULL` argument; when supplied, either
(a) pre-filter with `awk`/`grep` when available, or (b) read with
`data.table::fread(select = ...)` and drop non-matching `DATE` values before
parsing times. fluxqc's `import_li7810_subset()` implements (a) with an R
fallback and then calls `import.LI7810()` on the subset.

### Issue 5 (found while testing fluxqc): `import.*()` fail when called with `goFlux::`

**Title:** import.LI7810() and siblings error with "the condition has length
> 1" when called as goFlux::import.LI7810()

**Description**

Each importer starts with a deprecation shim of the form

```r
if (as.character(match.call()[[1]]) == "LI7810_import") { warning(...) }
```

When the function is called with a namespace prefix, `match.call()[[1]]` is
the call `goFlux::import.LI7810`, whose `as.character()` has length 3
(`"::"`, `"goFlux"`, `"import.LI7810"`), and R >= 4.2 turns the
length-3 condition into an error. The same happens when the importer is
passed as a function object to `do.call()`. Only a bare, attached call works.
Packages that depend on goFlux (like fluxqc) cannot use the recommended
`pkg::fun()` style.

**Minimal reproducible example**

```r
f <- system.file("extdata", "LI7810", "LI7810.data", package = "goFlux")
goFlux::import.LI7810(f, timezone = "UTC")
#> Error in if (as.character(match.call()[[1]]) == "LI7810_import") { :
#>   the condition has length > 1
library(goFlux); import.LI7810(f, timezone = "UTC")   # works
```

**Proposed fix:** `if (identical(deparse(match.call()[[1]]), "LI7810_import"))`
or, since the old names were removed a version ago, drop the shim.


### Ready-made patch for Issue 4: `import.LI7810(dates = NULL)`

`patches/import.LI7810-dates.patch` is a unified diff against
`R/import.LI7810.R` of the goFlux snapshot (v0.2.0, `~/Downloads/goFlux-master`,
RemoteSha 9c0b316 of the installed copy); `patches/import.LI7810.R` is the
whole modified file. R only, no `awk`. Apply with
`git apply inst/upstream/patches/import.LI7810-dates.patch` in a goFlux clone.

**PR title:** import.LI7810(): add `dates` argument to subset campaign-long
.data files before parsing (Fixes #<Issue 4>)

**PR description:**

> LI-7810 `.data` files grow to 100-200 MB over a campaign. `import.LI7810()`
> reads and time-parses every row even when only a few field days are needed.
> This PR adds `dates = NULL`: when a character vector of dates (as written in
> the file's `DATE` column) is given, the file is read with `readLines()`, the
> `DATAH` header row is kept and only the `DATA` rows whose `DATE` field is in
> `dates` are passed to `read.delim(text = ...)`; everything after that is
> unchanged. With `dates = NULL` the behaviour is byte-for-byte the current one
> (verified: `all.equal()` on the bundled `LI7810.data` is `TRUE`). A
> non-matching `dates` stops with a clear message. Roxygen `@param` added.
> Tested on the bundled example (330 rows) and on a 163 MB campaign file
> (six field days out of ~1.5 million rows) through the fluxqc shim
> `import_li7810_subset()`, which forwards to this argument when present.

Issue 5 (namespace-qualified calls) is a one-line change per importer
(`identical(deparse(match.call()[[1]]), "LI7810_import")` or dropping the
shim) and is best sent as its own tiny PR touching every `import.*.R`; it is
not bundled here to keep the `dates` PR single-purpose. Issue 3
(`muffleError`) needs the same restructuring in each of the ~20
`withCallingHandlers()` blocks of `import2RData.R` (wrap in `tryCatch()`,
keep `withCallingHandlers()` for warnings and messages), which is mechanical
but not small; the sketch in Issue 3 is the intended shape.

## (b) What belongs upstream and what stays in fluxqc

| fluxqc function | Upstream candidate? | Proposal |
|---|---|---|
| `click_peak2_stacked()` | **Yes** | `click.peak2(..., gases = NULL)`: when more than one gas is given, draw one `split.screen()` panel per gas sharing the time axis and identify on the top one. Output is unchanged. Small diff, no new dependency. |
| `precision_allan()`, `precision_by_id()` | **Yes** | A `prec.empirical()` helper (or a `prec = "empirical"` option in `goFlux()`) that fills the `*_prec` columns from `SD(diff(x))/sqrt(2)` of the flagged window, so `MDF` and `flux.plot()` use measured noise. Document that datasheet precision underestimates noise on field-worn analyzers (2.5x on the example closure). |
| `precision_mad()` | Maybe | Same helper with `method = "mad"`; document as the robust whole-record option. |
| `flag_detection()`, `mdf_multi()` | **Yes**, partially | `best.flux()` already applies an `MDF` criterion; propose an `MDF.method` argument (`"prec"`, `"two_sigma"`, `"christiansen"`, `"wassmann"`) and a `det_class` column (`emission` / `uptake` / `below detection`). Keep the full precision x method table in fluxqc. |
| `flux_pdf()` | **Yes** | Just the `inherits()` fix and the defensive `plot_env` read (Issue 1). |
| `import2RData_flat()` | **Yes** (as a fix) | Issues 2 and 3; the function itself need not go upstream if `import2RData()` is fixed. |
| `import_li7810_subset()` | **Yes** (as an argument) | Issue 4. |
| `auto_id_rise()`, `find_rise()`, `windows_from_table()`, `flag_window()` | Maybe later | goFlux has `autoID()` for LI-8200 only; a generic rise-based window finder could sit beside it. Propose after it has been used on more campaigns; keep in fluxqc for now. |
| `qc_screens()`, `qc_c0_screen()` | Stay | Campaign-level QC policy, not a flux-calculation step. |
| `process_fluxes()`, `write_outputs()` | Stay | Opinionated wrapper around goFlux; goFlux keeps its modular design. |
| `find_clock_offset()` | Stay | Analyzer-independent utility; could be linked from goFlux's documentation. |
| Rescue scripts in `inst/rescue/` (`detect_clock_streams()`, `separate_window_streams()`, `obs_win_deinterleave()`, `click_flux_2d()`, `write_stream_files()`) | **No** | One-off tools for a specific data-damage problem; no longer exported from fluxqc and not upstream candidates. |
| Ebullition add-on (`detect_ebullition()`, `split_ebullition()`, `merge_ebullition()`) | Stay | Out of goFlux's scope (it fits one smooth model per closure, and its source has no ebullition handling); the add-on produces a goFlux-shaped input so no goFlux change is needed. |
| `fluxqc_examples()`, `example_closure()` | Stay (data offer) | goFlux has its own example data; the four Yale-Myers teaching closures could be offered to goFlux as additional example data if wanted. |

## (c) Practical steps

1. Fork `Qepanna/goFlux` to `jongewirtzman/goFlux` and clone it. Add the
   upstream remote (`git remote add upstream https://github.com/Qepanna/goFlux`).
2. Open the five issues above on `Qepanna/goFlux` first, so each PR can say
   "Fixes #N". Bug-fix PRs are more likely to be merged quickly than feature
   PRs, so send them first, one branch per fix:
   - `fix/flux2pdf-class-check` (Issue 1; touches `R/flux2pdf.R` only)
   - `fix/import2RData-error-handling` (Issue 3; `R/import2RData.R`)
   - `fix/import-namespace-call` (Issue 5; one line per importer)
   - `feat/import2RData-outdir-recursive` (Issue 2)
   - `feat/import-LI7810-dates` (Issue 4; patch ready in `inst/upstream/patches/`)
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
4. Email Karelle Rheault (karh@ign.ku.dk) once the issues are open:
   introduce the chamber-flux use case (Yale-Myers, Harvard Forest, BlueFlux,
   Santa Fe, BCI campaigns, all processed with goFlux), point to
   the fluxqc repository and the issues, offer the PRs in the order above,
   and propose co-authoring a short package-update note (e.g. a JOSS
   "software update" or a methods note) covering stacked window selection,
   empirical precision and multi-method detection limits, with her as
   corresponding author of the package. Ask which of the feature PRs she
   would want in goFlux versus linked from its documentation.
5. After merges: bump fluxqc to depend on the new goFlux version, remove the
   duplicated functions (`flux_pdf()`, `click_peak2_stacked()`) with
   deprecation wrappers, and update the function map in `README.md`.
