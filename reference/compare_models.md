# Compare fitted models in sample

The table from
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
plus AICc (GWR), LOOIC (Bayesian GP), convergence and GP-adequacy flags
and the residual Moran's I of each fit. No model is refitted. For scores
on held-out blocks use
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md).

## Usage

``` r
compare_models(fits, newdata = NULL, ...)
```

## Arguments

- fits:

  A `spatial_fit`, or a named list of them. The names label the rows and
  must be unique (see
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md));
  a single fit is labelled by its class.

- newdata:

  Optional sf object of held-out rows, carrying the response and the
  predictors.

- ...:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html) when
  `newdata` is given. Without `newdata` nothing reads them, so any
  argument here is an error, which catches a misspelt `newdata`.

## Value

A data.frame with one row per model and these columns:

- `model` and the metric columns:

  As in
  [`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md).

- `metric_basis`:

  What each row's metrics were computed on (see
  [`evaluate_insample`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md)).
  A table that mixes `"out-of-bag"` and `"in-sample"` rows does not rank
  the models, and says so in the log.

- `AICc` (GWR) and `LOOIC` (Bayesian):

  Sums over the rows a model was fitted to, so each column is set to
  `NA`, with a warning, when the fits carrying it were fitted to
  different rows (a predictor with missing values drops rows, for
  example) or to different responses (a transformed response on the same
  rows), and the warning says which.

- `convergence_ok`:

  `TRUE` or `FALSE` for a Bayesian fit whose convergence was checked
  (see
  [`fit_bayesian_spatial_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md))
  and `NA` otherwise. A fit that did not converge is ranked like the
  others, so it also raises a warning.

- `gp_adequacy`:

  What the posterior of a Bayesian fit said about its Gaussian-process
  basis (`"ok"`, `"below_resolution"`, `"unresolved_short_range"`,
  `"no_signal_above_floor"` or `"not_checked"`; see "What the basis
  resolves" in
  [`fit_bayesian_spatial_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md))
  and `NA` for the other backends. A basis found too coarse scores like
  one that was adequate, so read the column beside the metrics.

- `resid_morans_I`, `resid_morans_z`, `resid_morans_p`:

  Moran's I of the residuals, from
  [`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md).

- `resid_morans_null`:

  The null each model was scored against, since that choice is per fit
  and governs how much the p-value is worth.

- `resid_morans_basis`:

  Which residuals the statistic was computed on: always those of the
  rows the model was fitted to, also when `newdata` is given, so the
  column reads `"in-sample"`, or `"out-of-bag"` for a forest, and never
  `"newdata"`. Beside `metric_basis = "newdata"` the Moran columns
  describe the training residuals, not the prediction errors on
  `newdata`, and the two can differ in sign.

A significant p-value is noted in the log (not raised as an R warning):
positive autocorrelation as structure the model may have missed,
negative (`resid_morans_z < 0`) as the alternating residuals of a model
that tracks its data closely. Read the caveats in
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)
before treating silence as evidence of no residual structure.

## Details

The metrics are those of
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md):
see its sections "What the metrics are computed on" and "Percentage
errors on responses with zeros".

## See also

Other in-sample checks:
[`evaluate_insample()`](https://elkronos.github.io/gis_modeling_toolkit/reference/evaluate_insample.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`print.morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.morans_i.md),
[`residual_morans_i()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residual_morans_i.md)

## Examples

``` r
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
# against the one that has it.  These are out-of-bag numbers;
# compare_models_cv() scores the same question on spatial folds.
fits <- list(
  a_only  = fit_rf_model(pts, "z", "a", num_trees = 100, seed = 1),
  a_and_b = fit_rf_model(pts, "z", c("a", "b"), num_trees = 100, seed = 1))
cmp <- compare_models(fits)
cmp[, c("model", "RMSE", "R2", "metric_basis", "resid_morans_I",
        "resid_morans_p")]
#>     model      RMSE        R2 metric_basis resid_morans_I resid_morans_p
#> 1  a_only 1.2209445 0.6804065   out-of-bag    -0.04684408      0.5967265
#> 2 a_and_b 0.7698265 0.8729450   out-of-bag    -0.06447357      0.3977229
```
