# Out-of-bag predictions from a random forest fit

Returns out-of-bag predictions rather than in-sample ones. See
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
for why, and for what it means for
[`summary()`](https://rdrr.io/r/base/summary.html).

## Usage

``` r
# S3 method for class 'rf_fit'
fitted(object, ...)
```

## Arguments

- object:

  An `rf_fit`.

- ...:

  Ignored.

## Value

Numeric vector of length `object$n`.
