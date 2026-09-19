# Cross-validate a random forest with spatial folds

Thin wrapper over
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
that refits a `ranger` forest on each training fold. Unlike the
out-of-bag error, this holds out whole spatial blocks, so neighbours of
a held-out point are not sitting in the training set.

## Usage

``` r
cv_rf(
  data_sf,
  response_var,
  predictor_vars,
  folds = NULL,
  k = 5,
  seed = 123,
  parallel = FALSE,
  block_size = NULL,
  auto_range = FALSE,
  boundary = NULL,
  pointize = "auto",
  metrics = NULL,
  ...
)
```

## Arguments

- data_sf:

  An sf object.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- folds:

  Optional fold definitions, in any of three shapes: a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  return value; a bare list of `list(train =, test =)` pairs of
  `..row_id` values; or a vector of fold labels, one per row, which
  becomes leave-that-label-out splits. The label vector is how folds
  built by another package are used here, since
  `blockCV::cv_spatial()`'s `$folds_list` holds two *unnamed* vectors
  per fold and is refused by name. That function returns a label vector
  as `$folds_ids`. Train and test must be disjoint: a fold that trains
  on its own test rows is not a cross-validation split and is refused
  with an error. IDs naming no row in the prepared data are dropped with
  a logged count.

- k:

  Number of folds when `folds` is `NULL`. Default 5.

- seed:

  RNG seed. Default 123. It seeds fold construction **and**, through a
  per-fold draw, each fold's forest, so two different seeds give
  different results even on identical `folds`. To grow every fold's
  forest from one fixed ranger seed instead, call
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  with your own `fit_fn` wrapping `fit_rf_model(seed = )`.

- parallel:

  Passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  Default `FALSE`. Under forked workers each fold's forest runs
  single-threaded unless `num_threads` is passed explicitly, so
  `parallel = 4` means four threads in total, where it would otherwise
  mean four times the session's `mc.cores`.

- block_size, auto_range, boundary:

  Passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).

- pointize:

  How non-POINT geometry is reduced to a point before fitting; passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
  Default `"auto"`.

- metrics:

  Optional scoring function of your own, passed to
  [`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md);
  see **Your own metrics** there.

- ...:

  Passed to
  [`fit_rf_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
  on every fold: `num_trees`, `mtry`, `importance`, `include_coords` and
  so on. `data_sf`, `response_var`, `predictor_vars` and
  `.already_prepped` are set by this function and must not be passed
  here (every fold would fail with "matched by multiple actual
  arguments"). A `seed` given here overrides the per-fold draw described
  above.

## Value

The
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
result.

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

## See also

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 200
  dat <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$z <- 2 * dat$a + rnorm(n, 0, 0.3)
  cv_rf(dat, "z", "a", k = 4)$overall
}
#> cv_rf(): no folds supplied -- using spatial block k-fold CV (k=4).
#>       RMSE       MAE     MAPE   SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 0.417048 0.3264843 83.10643 42.3203 0.9587525     NA    200    200     200
```
