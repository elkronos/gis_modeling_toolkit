# In-sample residuals from a GWR fit

Observed response minus
[`fitted.gwr_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
so these are **in-sample** residuals.

## Usage

``` r
# S3 method for class 'gwr_fit'
residuals(object, ...)
```

## Arguments

- object:

  A `gwr_fit`.

- ...:

  Ignored.

## Value

Numeric vector of length `object$n`.
