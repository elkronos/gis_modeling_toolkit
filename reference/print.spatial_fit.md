# Print a fitted spatial model

Shows the one-screen summary of a `gwr_fit`, a `bayesian_fit` or a
custom subclass: backend, formula, number of observations, CRS, and the
few backend-specific numbers worth seeing immediately (GWR bandwidth, GP
basis size). An `rf_fit` has its own method (see
[`print.rf_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md)),
which shows the same header plus the forest settings. It is what you get
by typing the object's name, and the quickest way to confirm a fit used
the data, predictors and CRS you meant. For fit quality use
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
or [`summary()`](https://rdrr.io/r/base/summary.html) instead. Nothing
printed here is an out-of-sample score.

## Usage

``` r
# S3 method for class 'spatial_fit'
print(x, ...)
```

## Arguments

- x:

  A `spatial_fit` object.

- ...:

  Ignored.

## Value

`x`, invisibly (called for its side effect).
