# Extract Bayesian model fixed-effect summaries

Returns the posterior summary of the global (non-spatial) regression
terms: estimate, error and credible interval per predictor, as
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
reports them. Reach for it to read the average effect of a predictor
with its uncertainty attached (the Bayesian counterpart to a coefficient
table), remembering that the Gaussian-process term has already absorbed
the spatially structured part of the signal, so these are effects net of
location.

## Usage

``` r
# S3 method for class 'bayesian_fit'
coef(object, ...)
```

## Arguments

- object:

  A `bayesian_fit` object.

- ...:

  Ignored.

## Value

A matrix of fixed-effect posterior summaries, as returned by
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html).
Never `NULL`: a missing 'brms' or a failing `fixef()` call errors,
following the [`coef()`](https://rdrr.io/r/stats/coef.html) contract
described in
[`new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).

## See also

Other methods on a fitted model:
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
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
