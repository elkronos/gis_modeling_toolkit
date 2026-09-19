# Predict from a random forest fit

With `newdata = NULL` this returns **out-of-bag** predictions, not
in-sample ones. See
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md).

## Usage

``` r
# S3 method for class 'rf_fit'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  An `rf_fit`.

- newdata:

  Optional sf object carrying the same predictors. It is transformed to
  the CRS used at fitting time first, so a forest that includes the
  coordinates is not fed a different coordinate system. Categorical
  predictors must not carry a level the forest was never grown with,
  meaning a level with **no training rows**, not merely one absent from
  [`levels()`](https://rdrr.io/r/base/levels.html): an ordinary subset,
  or a spatial-CV fold that holds out a whole class, keeps the unused
  level while the forest has no split for it. An unseen level is an
  error, not a guess. A predictor that was numeric or logical when the
  forest was grown must also arrive numeric or logical: ranger would
  otherwise factor-code a character column and apply the numeric split
  thresholds to the codes, predicting confidently from nonsense, so a
  character column is refused instead. (ranger sees a logical as 0/1, so
  either form is accepted for one.)

- ...:

  Passed to `ranger`'s predict method. Arguments that make `ranger`
  return a matrix (`predict.all = TRUE`, `type = "quantiles"`,
  `type = "se"` with `predict.all`) are rejected, because this method's
  contract is one number per row of `newdata`. Call
  `predict(fit$engine, data = ...)` directly for those. `seed` defaults
  to a constant: an unset `seed` makes `ranger` draw one uniform from
  the global RNG stream per call, so the number of
  [`predict()`](https://rdrr.io/r/stats/predict.html) calls a script
  happens to make (via
  [`predict_surface`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)'s
  `chunk_size`, say) would otherwise shift every later random draw. It
  does not affect a regression forest's predictions; pass your own if
  you need one.

## Value

Numeric vector, aligned to `nrow(newdata)` with `NA` for rows dropped as
incomplete.

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
[`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
