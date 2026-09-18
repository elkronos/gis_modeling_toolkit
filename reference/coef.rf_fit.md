# Coefficients are undefined for a random forest

Consistent with
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md)
and
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
which also error rather than returning `NULL` when they cannot supply
coefficients – see the [`coef()`](https://rdrr.io/r/stats/coef.html)
contract in
[`new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).

## Usage

``` r
# S3 method for class 'rf_fit'
coef(object, ...)
```

## Arguments

- object:

  An `rf_fit`.

- ...:

  Ignored.

## Value

Never returns; always signals an error.
