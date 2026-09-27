# Side-by-side comparison of fitted spatial models

Takes a named list of already-fit `spatial_fit` objects and produces a
tidy comparison table including in-sample (for a forest, out-of-bag)
metrics and model-specific information criteria (AICc, LOOIC).

## Usage

``` r
compare_models(fits, newdata = NULL, ...)
```

## Arguments

- fits:

  A named list of `spatial_fit` objects. Names must be unique; see
  [`evaluate_insample`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md).

- newdata:

  Optional sf for out-of-sample evaluation.

- ...:

  Extra arguments passed to predict().

## Value

A data.frame comparing all models. Its `metric_basis` column says what
each row's metrics were computed on (see
[`evaluate_insample`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md));
a table that mixes `"out-of-bag"` and `"in-sample"` rows does not rank
the models, and says so in the log. `AICc` (GWR) and `LOOIC` (Bayesian)
are sums over the rows a model was fitted to, so each column is set to
`NA`, with a warning, when the fits carrying it were fitted to different
rows (a predictor with missing values drops rows, for example) or to
different responses (a transformed response on the same rows), and the
warning says which. `convergence_ok` is `TRUE` or `FALSE` for a Bayesian
fit whose convergence was checked (see
[`fit_bayesian_spatial_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md))
and `NA` otherwise; a fit that did not converge is ranked like the
others, so it also raises a warning. Alongside the metrics it carries
`resid_morans_I`, `resid_morans_z`, `resid_morans_p` and
`resid_morans_null`, the last of which names the null
[`residual_morans_i`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
scored each model against, since that choice is per-fit and governs how
much the p-value is worth. A significant p-value is noted in the log
(not raised as an R warning): positive autocorrelation as structure the
model may have missed, negative (`resid_morans_z < 0`) as the
alternating residuals of a model that tracks its data closely. Read the
caveats in
[`?residual_morans_i`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
before treating silence as evidence of no residual structure.

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
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)
and `compare_models()` record it per model in a `metric_basis` column.
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
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  pts <- st_as_sf(
    data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
               a = rnorm(60), b = rnorm(60)),
    coords = c("x", "y"), crs = 32632
  )
  # An offset keeps the response away from zero, where MAPE is undefined.
  pts$z <- 10 + 2 * pts$a - pts$b + rnorm(60, 0, 0.3)
  # Two predictor sets for the same learner: the model that leaves `b` out
  # against the one that has it.  These are in-sample (out-of-bag) numbers;
  # compare_models_cv() scores the same question on spatial folds.
  fits <- list(a_only = fit_rf_model(pts, "z", "a", num_trees = 100, seed = 1),
               a_and_b = fit_rf_model(pts, "z", c("a", "b"), num_trees = 100, seed = 1))
  compare_models(fits)
}
#>     model  n      RMSE       MAE      MAPE    SMAPE        R2 Adj_R2 n_MAPE
#> 1  a_only 60 1.2209445 0.9721158 10.181689 9.914608 0.6804065     NA     60
#> 2 a_and_b 60 0.7698265 0.5989365  6.332639 6.133211 0.8729450     NA     60
#>   n_SMAPE metric_basis AICc LOOIC bandwidth_is_fallback convergence_ok
#> 1      60   out-of-bag   NA    NA                    NA             NA
#> 2      60   out-of-bag   NA    NA                    NA             NA
#>   resid_morans_I resid_morans_z resid_morans_p resid_morans_null
#> 1    -0.04684408     -0.5291138      0.5967265     randomisation
#> 2    -0.06447357     -0.8456951      0.3977229     randomisation
```
