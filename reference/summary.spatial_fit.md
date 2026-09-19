# Summarise a fitted spatial model

Computes goodness-of-fit metrics from `fitted(object)` against the
observed response. A non-numeric response is an error: a character or
factor response cannot be scored, and used to come back as `n = 0` with
every metric `NA`; a logical response is treated as 0/1.

## Usage

``` r
# S3 method for class 'spatial_fit'
summary(object, ...)
```

## Arguments

- object:

  A `spatial_fit` object.

- ...:

  Ignored.

## Value

An object of class `summary.spatial_fit`: a list with `class`,
`formula`, `n`, `response_var`, `predictor_vars`, `info` and `in_sample`
(the metric data.frame, out-of-bag for an `rf_fit`).

## What the metrics are computed on

For a `gwr_fit` or a `bayesian_fit` these are **in-sample** metrics:
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns values
computed at the training locations from the model that saw them. For an
`rf_fit` they are **out-of-bag**, because
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md)
returns out-of-bag predictions in place of in-sample ones (see
[`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
The two are not comparable, and
[`print()`](https://rdrr.io/r/base/print.html) on the result labels
which one it is holding, driven by `$info$fitted_are_oob`. Use
[`compare_models_cv`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
to compare backends.

Adjusted R-squared is suppressed: GWR's effective parameter count far
exceeds the global predictor count, and a GP model has no simple `p`.

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
