# Turn the console log off or on

Turns the package's console log off (`quiet = TRUE`) or back on. The
package writes a trace of INFO lines and above to
`file.path(tempdir(), "spatialkit_model_log.log")`, a file that lasts as
long as the R session, and echoes WARN lines and above to the console;
this function changes only the console echo.
[`logger::log_threshold()`](https://daroczig.github.io/logger/reference/log_threshold.html)
changes the log file's threshold by default (`index = 1`), not the
console echo's, so use this function for the echo.

## Usage

``` r
spatialkit_quiet(quiet = TRUE)
```

## Arguments

- quiet:

  `TRUE` (default) silences the console echo; `FALSE` restores the
  package default, WARN+. A logger threshold
  ([`logger::ERROR`](https://daroczig.github.io/logger/reference/log_levels.html),
  or the value a previous call returned) sets that level instead, which
  is how to put back exactly what was in force rather than the default:
  `old <- spatialkit_quiet(); spatialkit_quiet(old)`.

## Value

Invisibly, the threshold that was in force before the change, a logger
level that can be passed back as `quiet`.

## Details

These are log records, not R conditions:
[`suppressWarnings()`](https://rdrr.io/r/base/warning.html) and
`tryCatch(warning = )` do not see them. Conditions the package raises as
R warnings are unaffected by this function (see "Warnings and log lines"
in
[spatialkit-conditions](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)).

While a document is being knitted (R Markdown, Quarto, a pkgdown
article) the console echo is also sent as an R message, because knitr
does not capture what is written to the console's error stream and the
cautions would otherwise be missing from the output. They appear as
`## WARN [...]` lines, and the chunk option `message = FALSE`, like
[`suppressMessages()`](https://rdrr.io/r/base/message.html), keeps them
out of the document. A line the package also raises as an R warning is
not repeated, since the document shows the warning. Outside knitr
nothing changes.

## See also

Other options and caches:
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md),
[`clear_grid_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md),
[`spatialkit-conditions`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md),
[`spatialkit-options`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md)

## Examples

``` r
old <- spatialkit_quiet()      # console echo off
spatialkit_quiet(old)          # back to whatever it was
spatialkit_quiet(FALSE)        # or back to the WARN+ default
```
