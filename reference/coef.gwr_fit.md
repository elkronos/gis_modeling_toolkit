# Extract GWR local coefficients

Returns the whole surface of coefficients, one row per observation and
one column per term, against the single global vector
[`coef()`](https://rdrr.io/r/stats/coef.html) returns for an `lm`. That
table is the point of fitting a GWR at all: inspect the spread of a
predictor's column to see where, and by how much, its relationship with
the response changes across the study area, and join it back to
`object$data_sf` to map it. Use
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)
for a quick look at that map.

## Usage

``` r
# S3 method for class 'gwr_fit'
coef(object, ...)
```

## Arguments

- object:

  A `gwr_fit` object.

- ...:

  Ignored.

## Value

A data.frame of local coefficient estimates: one row per observation,
one column per model term. Never `NULL`: when the engine carries no
`SDF` component this errors, following the
[`coef()`](https://rdrr.io/r/stats/coef.html) contract described in
[`new_spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).

## What is and is not returned

Only the model terms: the intercept and one column per predictor.
GWmodel's `SDF` data slot carries a good deal more alongside them
(standard errors, t-values, the observed response, the fitted values,
the residuals, `Local_R2`): 15 columns for a two-predictor fit, of which
3 are coefficients. Returning the whole slot would have made
`coef(fit)$Local_R2` and `coef(fit)$a_SE` read like coefficients and
`ncol(coef(fit))` a meaningless number. Reach for `object$engine$SDF`
when you want the rest; it is the unmodified GWmodel object.

If the model terms cannot be located in the `SDF` (a GWmodel that names
its coefficient columns differently), the whole slot is returned with a
warning saying so. The call neither errors nor returns a silently short
table.

## See also

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
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
