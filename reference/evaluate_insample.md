# Compute in-sample (or out-of-sample) metrics for fitted spatial models

Accepts a single `spatial_fit` object or a named list of them. Does NOT
refit. Uses [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) for
in-sample and [`predict()`](https://rdrr.io/r/stats/predict.html) for
new data.

## Usage

``` r
evaluate_insample(fits, newdata = NULL, ...)
```

## Arguments

- fits:

  A `spatial_fit` object, or a named list of them (e.g.
  `list(GWR = gwr_obj, Bayesian = bayes_obj)`). The names are used as
  the model labels and every element must have one; an unnamed list is
  an error, and so are duplicated names. `model` is the key the
  comparison table is assembled on, so two fits sharing a name cannot be
  told apart in the output.

- newdata:

  Optional sf object for out-of-sample evaluation. Must contain the
  response variable and all predictors. If NULL, in-sample metrics are
  computed.

- ...:

  Extra arguments passed to predict().

## Value

A data.frame with one row per model and columns for model name, all
regression metrics, and `metric_basis`: what the row's metrics were
computed on, `"in-sample"` (fitted values), `"out-of-bag"` (an
`rf_fit`'s fitted values, see "What the metrics are computed on") or
`"newdata"`. Rows with different bases do not compare like for like. An
element that is not a `spatial_fit` is skipped, with a logged warning,
and has no row; a list in which no element is a `spatial_fit` is an
error.

## What the metrics are computed on

With `newdata = NULL` the metrics come from `fitted(object)`. That is
**in-sample** for a `gwr_fit` or a `bayesian_fit`, but **out-of-bag**
for an `rf_fit`, whose
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) method returns
out-of-bag predictions (see
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
The data.frame
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
returns carries no label distinguishing the two, so check
`object$info$fitted_are_oob` before comparing numbers across backends;
`evaluate_insample()` and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
record it per model in a `metric_basis` column.
[`compare_models_cv`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
scores every backend the same way.

\\R^2\\ is \\1 - RSS/TSS\\ with the total sum of squares taken about the
mean of the response the model was *fitted* to. In sample that is the
ordinary \\R^2\\. With `newdata` it is out-of-sample \\R^2\\, the
convention every `cv_*()` function uses: the model is measured against
the prediction it had to beat, the training mean, not against the new
rows' own mean, which it could not have known. It is below 0 when the
model predicts the new rows worse than the training mean does, and it is
`NA` when the response does not vary about that baseline by more than
rounding error (100 machine epsilons of its magnitude, whatever its
units).

## Percentage errors on responses with zeros

`MAPE` divides by the observed value and `SMAPE` by \\\|y\| +
\|\hat{y}\|\\, so neither is defined where its denominator is zero.
Neither returns `Inf` or `NaN`. Both are averaged over the rows whose
denominator is non-zero, and are `NA` when no row qualifies. Non-zero is
judged at the scale of the data: a denominator no larger than 100
machine epsilons times the largest one counts as zero, so the rule does
not depend on the units of the response. The `n_MAPE` and `n_SMAPE`
columns record how many rows that was; the `n` column counts finite
observation/prediction pairs. Read a percentage error next to its count:
when `n_MAPE < n`, `MAPE` is an average over a subset of the data,
whatever its value.

This bites on any response taking exact zeros: counts, rainfall,
abundance, claim amounts. On a zero-inflated response with 62 zeros out
of 120, `MAPE` is an average over the 58 non-zero rows, which
`n_MAPE = 58` now says. `SMAPE` fails differently and more subtly: it
drops the rows where observation and prediction are both near zero
(which on a well-fitted zero-inflated model are the rows it got
*right*), so it averages the harder rows only and reads worse than the
fit deserves; `n_SMAPE` shows how many rows it kept, and the count is
only a label, not a repair.

`RMSE`, `MAE` and \\R^2\\ use every finite row and are unaffected;
prefer them whenever the response can be zero. For a Bayesian fit,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
additionally reports CRPS and interval coverage, which are proper
scoring rules and have no such failure mode.

## See also

Other model evaluation:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  pts <- st_as_sf(
    data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
               a = rnorm(60)),
    coords = c("x", "y"), crs = 32632
  )
  pts$z <- 2 * pts$a + rnorm(60, 0, 0.3)
  # Fit on 40 rows and keep 20 back, so the second call really is out of
  # sample; scoring the training rows again would only re-read the fit.
  fit <- fit_rf_model(pts[1:40, ], "z", "a", num_trees = 50, seed = 1)
  print(evaluate_insample(fit))                    # in-sample (out-of-bag for RF)
  evaluate_insample(fit, newdata = pts[41:60, ])   # on the 20 held-out rows
}
#>    model  n     RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_MAPE
#> 1 rf_fit 40 0.409009 0.3236812 86.57837 50.87364 0.9541711     NA     40
#>   n_SMAPE metric_basis
#> 1      40   out-of-bag
#>    model  n      RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_MAPE
#> 1 rf_fit 20 0.5080838 0.3836504 158.4334 63.30707 0.9022065     NA     20
#>   n_SMAPE metric_basis
#> 1      20      newdata
```
