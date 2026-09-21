# Compute goodness-of-fit metrics for a spatial model

Reports RMSE, MAE, MAPE, SMAPE, \\R^2\\ and adjusted \\R^2\\ for any
`spatial_fit`, in one row and on one scale, so that fits from different
backends can be read side by side. Reach for it to score a model on data
you hold out yourself (pass it as `newdata`), or to get a quick
in-sample reading of how closely a fit tracks its training data.

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

  Additional arguments passed to predict().

- newdata:

  Optional sf object for out-of-sample evaluation. If NULL, fitted
  values are used (in-sample, or out-of-bag for an `rf_fit`; see above).

## Value

A data.frame with n, RMSE, MAE, MAPE, SMAPE, R2, Adj_R2, n_MAPE and
n_SMAPE (the last two are the rows each percentage error was averaged
over; see "Percentage errors on responses with zeros"). `Adj_R2` is
always `NA`: GWR's effective parameter count far exceeds the global
predictor count and a GP model has no simple `p`, so it is deliberately
suppressed. A non-numeric response is an error: a character or factor
response cannot be scored, and used to come back as `n = 0` with every
metric `NA`; a logical response is treated as 0/1.

## Details

It is not a substitute for cross-validation. With `newdata = NULL` the
numbers are in-sample for a `gwr_fit` or `bayesian_fit`, and a GWR can
reach a near-perfect in-sample \\R^2\\ at a small bandwidth without
predicting anything. For a figure you can report, use
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
or
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md).

## What the metrics are computed on

With `newdata = NULL` the metrics come from `fitted(object)`. That is
**in-sample** for a `gwr_fit` or a `bayesian_fit`, but **out-of-bag**
for an `rf_fit`, whose
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) method returns
out-of-bag predictions (see
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
The returned data.frame carries no label distinguishing the two, so
check `object$info$fitted_are_oob` before comparing numbers across
backends, or use
[`compare_models_cv`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
which scores every backend the same way.

## Percentage errors on responses with zeros

`MAPE` divides by the observed value and `SMAPE` by \\\|y\| +
\|\hat{y}\|\\, so neither is defined where its denominator is zero.
Neither returns `Inf` or `NaN`. Both are averaged over the rows whose
denominator is non-zero, and are `NA` when no row qualifies. The
`n_MAPE` and `n_SMAPE` columns record how many rows that was; the `n`
column counts finite observation/prediction pairs. Read a percentage
error next to its count: when `n_MAPE < n`, `MAPE` is an average over a
subset of the data, whatever its value.

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

## Which metrics survive a non-Gaussian response

RMSE and MAE are defined for any numeric response and are what to read
for a count, a rate or a bounded outcome. MAPE and SMAPE assume a
response that is rarely zero (see the previous section), and R-squared
and adjusted R-squared compare residual variance to total variance,
which is the right comparison for a Gaussian response and a loose one
for anything whose variance tracks its mean. None of the four is wrong
to compute; each is Gaussian-shaped thinking, and on a Poisson or
zero-inflated response should be read as a rough summary rather than a
score.

For the Bayesian backend,
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
additionally reports CRPS and interval coverage at 50, 80 and 95
percent. Both are proper scoring rules computed from posterior draws, so
they are meaningful for any `family` the backend accepts, and they are
the numbers to compare when the response is not Gaussian. When every
fold fails, the `fold_metrics` frame
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
returns carries the CRPS column but not the `coverage_*` columns, so
code that reads those columns must tolerate their absence.

## See also

Other model evaluation:
[`compare_models()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models.md),
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md),
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md),
[`fitted.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md),
[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
[`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md),
[`predict.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.gwr_fit.md),
[`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md),
[`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)

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
  print(model_metrics(fit))                    # in-sample (out-of-bag for RF)
  model_metrics(fit, newdata = pts[41:60, ])   # on the 20 held-out rows
}
#>    n     RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_MAPE n_SMAPE
#> 1 40 0.409009 0.3236812 86.57837 50.87364 0.9541711     NA     40      40
#>    n      RMSE       MAE     MAPE    SMAPE        R2 Adj_R2 n_MAPE n_SMAPE
#> 1 20 0.5080838 0.3836504 158.4334 63.30707 0.8999583     NA     20      20
```
