# Predict from a Bayesian spatial GP model

Applies the same newdata preparation pipeline as
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md):
non-point geometries are coerced to points, the data is projected to the
CRS used during fitting (via
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)),
and rows with missing or non-finite values are dropped. Coordinate
scaling and predictor standardisation stored at fit time are then
applied before delegating to
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
or
[`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

## Usage

``` r
# S3 method for class 'bayesian_fit'
predict(
  object,
  newdata = NULL,
  summary = c("mean", "median"),
  type = c("epred", "predict"),
  draws = FALSE,
  ...
)
```

## Arguments

- object:

  A `bayesian_fit` object.

- newdata:

  An sf object with the same predictors. The response variable need not
  be present (true out-of-sample prediction is supported). NULL = fitted
  values.

- summary:

  "mean" (default) or "median" over posterior draws.

- type:

  "epred" (default) for expected predictions (no obs noise), or
  "predict" for full posterior predictive draws (includes obs noise).

- draws:

  If TRUE, return the full posterior draw matrix instead of a point
  summary. Default FALSE.

- ...:

  Ignored.

## Value

Numeric vector of length `nrow(newdata)`, or a `n_draws x nrow(newdata)`
matrix when `draws = TRUE` (a 1-row all-`NA` matrix if the posterior
draw fails). With `newdata = NULL` the cached
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) values are
returned only for the default `summary = "mean"`, `type = "epred"`,
`draws = FALSE` combination; any other combination is recomputed against
the training data, because the cache holds epred column means and
nothing else.

## The GP boundary is pinned

brms 2.x does not store the Hilbert-space boundary \\L\\ in a fitted GP
basis, so `brms:::.data_gp()` recomputes it from whatever rows
[`predict()`](https://rdrr.io/r/stats/predict.html) is handed, which
moved every eigenfunction of the approximation with the newdata bounding
box while the fitted basis coefficients stayed put. Two synthetic rows
at the training coordinate extrema are therefore appended before the
posterior draw and dropped from the result, reproducing the boundary the
model was fitted with, so chunked, fold-wise and single-call predictions
agree.

That is exact only for `newdata` **inside** the training coordinate
envelope. Beyond it the boundary has to grow whatever is done, so
predictions there are extrapolation from a basis that was not built for
them *and* depend on which other rows share the call, including on
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)'s
`chunk_size`. A notice is written to the log (not raised as a warning)
when it happens.

## See also

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md),
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md),
[`fitted.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md),
[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md),
[`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md),
[`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
