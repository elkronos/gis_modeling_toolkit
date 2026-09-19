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

## See also

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md),
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md),
[`fitted.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md),
[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md),
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md),
[`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md),
[`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
