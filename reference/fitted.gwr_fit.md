# In-sample fitted values from a GWR fit

Reads the fitted values out of the GWmodel result, which stores them
under one of several names depending on version and entry point; the
extraction falls back through the local coefficients and the residuals
when no direct column is present. These are **in-sample** values (each
observation was inside its own bandwidth window), so
[`summary()`](https://rdrr.io/r/base/summary.html) on a `gwr_fit`
reports an optimistic fit. Use
[`cv_gwr`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
for a spatially blocked estimate.

## Usage

``` r
# S3 method for class 'gwr_fit'
fitted(object, ...)
```

## Arguments

- object:

  A `gwr_fit`.

- ...:

  Ignored.

## Value

Numeric vector of length `object$n` (`NA` where extraction failed).
