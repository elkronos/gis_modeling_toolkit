# Compute goodness-of-fit metrics for a fitted model

Reports RMSE, MAE, MAPE, SMAPE, \\R^2\\ and adjusted \\R^2\\ for any
`spatial_fit` as a one-row data frame. With `newdata` it scores rows you
held out yourself. Without it, it reads the fit's own fitted values:
in-sample for a GWR or a Bayesian GP fit, out-of-bag for a random
forest. Those numbers do not compare across backends, and the data frame
carries no label saying which kind it holds.

## Usage

``` r
model_metrics(object, ...)

# S3 method for class 'spatial_fit'
model_metrics(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `spatial_fit` object.

- ...:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html) when
  `newdata` is given. Without `newdata` nothing reads them, so any
  argument here is an error, which catches a misspelt `newdata`. The
  metrics are computed from one predicted value per row, so
  `draws = TRUE` (a Bayesian GP fit) and `se.fit = TRUE` (a GWR fit),
  with which [`predict()`](https://rdrr.io/r/stats/predict.html) returns
  a matrix or a list, are refused.

- newdata:

  Optional sf object of held-out rows, carrying the response and the
  predictors. If `NULL`, the fitted values are used (in-sample, or
  out-of-bag for an `rf_fit`; see "What the metrics are computed on").

## Value

A one-row data.frame with `n`, `RMSE`, `MAE`, `MAPE`, `SMAPE`, `R2`,
`Adj_R2`, `n_MAPE` and `n_SMAPE` (the last two are the rows each
percentage error was averaged over). `Adj_R2` is always `NA`: a GWR's
effective number of parameters far exceeds the global predictor count
and a GP model has no simple `p`. A character or factor response is an
error; a logical response is scored as 0/1.

## Details

Use this function for one fit.
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
runs it over a list of fits and labels each row with its basis, and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
adds information criteria and residual Moran's I to that table. None of
the three refits a model. An in-sample score is not a cross-validation
score: a GWR can reach a near-perfect in-sample \\R^2\\ at a small
bandwidth without predicting anything. For a figure you can report, use
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
or
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md).

## What the metrics are computed on

With `newdata = NULL` the metrics come from `fitted(object)`. That is
in-sample for a GWR or a Bayesian GP fit, and out-of-bag for a random
forest, whose [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
method returns out-of-bag predictions (see
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
Check `object$info$fitted_are_oob` before comparing numbers across
backends;
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
and
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md)
record it per model in a `metric_basis` column, and
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
scores every backend the same way.

\\R^2\\ is \\1 - RSS/TSS\\ with the total sum of squares taken about the
mean of the response the model was fitted to. In sample that is the
ordinary \\R^2\\. With `newdata` it is the out-of-sample \\R^2\\ every
`cv_*()` function uses: the model is measured against the training mean,
the prediction it had to beat, and not against the new rows' own mean,
which it could not have known. It is below 0 when the model predicts the
new rows worse than the training mean does, and `NA` when the response
does not vary about that mean beyond the rounding of its values,
whatever its units and origin.

A GWR fit whose response
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
took about a centre (`fit$engine$centres$response`, a response at a
large offset) is scored on the response and the fitted values or
predictions less that centre, which keep digits that rounding at the
response's size would lose: `RMSE`, `MAE` and \\R^2\\ are those of the
same values at 0, and `MAPE` and `SMAPE` divide the same errors by the
values in the response's units.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) still return values
in the response's units.

## Percentage errors on responses with zeros

`MAPE` divides by the observed value and `SMAPE` by \\\|y\| +
\|\hat{y}\|\\, so neither is defined where its denominator is zero.
Neither returns `Inf` or `NaN`. Both are averaged over the rows whose
denominator is non-zero, and are `NA` when no row qualifies. Non-zero is
judged at the scale of the data: a denominator no larger than 100
machine epsilons times the largest one counts as zero, so the rule does
not depend on the units of the response. The `n_MAPE` and `n_SMAPE`
columns record how many rows that was; the `n` column (`n_pred` in a
cross-validation result) counts finite observation/prediction pairs.
Read a percentage error next to that count: when `n_MAPE` is smaller,
`MAPE` is an average over a subset of the data, whatever its value.

Any response taking exact zeros is affected: counts, rainfall,
abundance, claim amounts. On a zero-inflated response with 62 zeros out
of 120, `MAPE` averages the 58 non-zero rows and `n_MAPE` is 58. `SMAPE`
drops the rows where observation and prediction are both near zero,
which on a well-fitted zero-inflated model are the rows it got right, so
it averages the harder rows only and reads worse than the fit deserves.
`n_SMAPE` shows how many rows it kept; it does not correct the average.

`RMSE`, `MAE` and \\R^2\\ use every finite row; prefer them whenever the
response can be zero. For a Bayesian GP fit,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
also reports CRPS and interval coverage, which are proper scoring rules
and are defined at zero.

## Which metrics survive a non-Gaussian response

RMSE and MAE are defined for any numeric response and are what to read
for a count, a rate or a bounded outcome. MAPE and SMAPE assume a
response that is rarely zero (see "Percentage errors on responses with
zeros"), and \\R^2\\ and adjusted \\R^2\\ compare residual variance to
total variance, which suits a Gaussian response and fits loosely
anything whose variance tracks its mean. On a Poisson or zero-inflated
response read the four as a rough summary rather than a score.

For the Bayesian GP,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
also reports CRPS and interval coverage at its `coverage_levels` (50, 80
and 95 percent by default). Both are proper scoring rules computed from
posterior draws, so they are meaningful for any `family` that predicts
one number per row (a count, a rate, a binary or bounded outcome), and
they are the numbers to compare when the response is not Gaussian. A
categorical or ordinal family predicts a probability per response
category instead, so
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
refuses one before fitting anything. When every fold fails, the
`fold_metrics` frame
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
returns carries the `CRPS` column but not the `coverage_*` columns, so
code that reads those columns must allow for their absence.

## See also

Other in-sample checks:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
             a = rnorm(60)),
  coords = c("x", "y"), crs = 32632
)
# An offset keeps the response away from zero, where MAPE is undefined.
pts$z <- 10 + 2 * pts$a + rnorm(60, 0, 0.3)
# Fit on 40 rows and keep 20 back, so the second call is out of sample.
fit <- fit_rf_model(pts[1:40, ], "z", "a", num_trees = 50, seed = 1)
model_metrics(fit)                           # out-of-bag for a forest
#>    n     RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_MAPE n_SMAPE
#> 1 40 0.409009 0.3236812 3.425689 3.390414 0.9541711     NA     40      40
model_metrics(fit, newdata = pts[41:60, ])   # on the 20 held-out rows
#>    n      RMSE       MAE    MAPE    SMAPE        R2 Adj_R2 n_MAPE n_SMAPE
#> 1 20 0.5080838 0.3836504 3.56621 3.649174 0.9022065     NA     20      20
```
