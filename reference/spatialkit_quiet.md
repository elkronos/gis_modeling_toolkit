# Quieten (or restore) spatialkit's console log

The package writes an INFO+ trace to a session temp file (logger
appender index 1) and echoes WARN+ to the console (index 2). Both
[`logger::log_appender()`](https://daroczig.github.io/logger/reference/log_appender.html)
and
[`logger::log_threshold()`](https://daroczig.github.io/logger/reference/log_threshold.html)
default to `index = 1`, so the obvious two-line recipe silences the file
and leaves the console untouched, which is the opposite of what anyone
wants. This helper names the right index.

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
real R warnings are unaffected by this function.

## Examples

``` r
old <- spatialkit_quiet()      # console echo off
spatialkit_quiet(old)          # back to whatever it was
spatialkit_quiet(FALSE)        # or back to the WARN+ default
```
