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
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html), on
the fitted scale (per standard deviation of each predictor under
`standardize_predictors = TRUE`; see above). Never `NULL`: a missing
'brms' or a failing `fixef()` call errors, following the
[`coef()`](https://rdrr.io/r/stats/coef.html) contract described in
[`new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).

## Standardised predictors

The summaries are on the scale the model was fitted on. A fit made with
`standardize_predictors = TRUE` was fitted on centred and scaled numeric
predictors, so each slope is the change in the linear predictor per
*standard deviation* of its predictor and the intercept is its value at
the predictor *means*, not the raw-unit numbers
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) reports on the same
formula. Nothing on the returned matrix says so;
[`print()`](https://rdrr.io/r/base/print.html) on the fit does, and the
centre and scale of each predictor are in
`object$info$predictor_scaling`. To put a slope back in raw units divide
its `Estimate`, `Est.Error` and interval bounds by that predictor's
`scale`. The intercept's `Estimate` follows by linearity (subtract each
raw-unit slope times its predictor's `center`), but its `Est.Error` and
interval depend on the posterior covariance of the coefficients:
transform the draws from `brms::as_draws_df(object$engine)` for those,
or refit without standardising.

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
