# K-fold cross-validation for GWR

Refits a geographically weighted regression from scratch on each
training fold and scores it on the held-out fold, so the reported error
is what the model achieves at locations it did not see. Reach for it
whenever you need a defensible accuracy figure for a GWR: the in-sample
\\R^2\\ that
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
reports on a `gwr_fit` is close to meaningless, because a local
regression with a small bandwidth can track the training points almost
exactly. Bandwidth is re-selected per fold unless you fix it with
`bandwidth`, which keeps the selection itself inside the
cross-validation, with no tuning on the full data first.

## Usage

``` r
cv_gwr(
  data_sf,
  response_var,
  predictor_vars,
  folds = NULL,
  k = 5,
  seed = 123,
  adaptive = TRUE,
  bandwidth = NULL,
  kernel = c("bisquare", "gaussian", "tricube", "boxcar", "exponential"),
  boundary = NULL,
  pointize = "auto",
  block_size = NULL,
  auto_range = FALSE,
  parallel = FALSE,
  metrics = NULL
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
  built by another package are used here (`blockCV::cv_spatial()`
  returns one as `$folds_ids`), since its `$folds_list` holds two
  *unnamed* vectors per fold and is refused by name. Train and test must
  be disjoint: a fold that trains on its own test rows is not a
  cross-validation split and is refused with an error. IDs naming no row
  in the prepared data are dropped with a logged count (expected when
  rows were removed for missing values; a sign the folds came from other
  data when they were not).

- k:

  Number of folds. Default 5.

- seed:

  RNG seed. Default 123.

- adaptive:

  Logical; use adaptive bandwidth. Default TRUE.

- bandwidth:

  Optional bandwidth, applied to every fold. For `adaptive = TRUE` an
  integer number of neighbours; for `adaptive = FALSE` a distance in the
  units of the **projected** CRS the folds are fitted in (geographic
  input is projected first, so a value in degrees is read as metres).
  `NULL` (default) selects a bandwidth per fold with
  [`GWmodel::bw.gwr()`](https://rdrr.io/pkg/GWmodel/man/bw.gwr.html),
  which is reported in `fold_metrics$bandwidth`. See
  [`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md).

- kernel:

  Kernel function type.

- boundary:

  Optional polygonal sf/sfc for CRS alignment.

- pointize:

  Geometry coercion strategy.

- block_size:

  Optional minimum block edge length for spatial CV blocks (projected
  CRS units). Blocks are then at least as large as the spatial
  autocorrelation range.

- auto_range:

  Logical. If `TRUE` and `folds` is `NULL`, estimate the autocorrelation
  range and use it as the minimum block size. Default `FALSE`.

- parallel:

  Logical or positive integer. If `TRUE`, auto-detect the number of
  cores and fit folds in parallel via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (macOS / Linux; falls back to sequential on Windows). If an integer \>
  1, use that many cores. Default `FALSE` (sequential).

- metrics:

  Optional scoring function of your own, a `function(y, yhat)` returning
  a named numeric vector; its names become columns of `fold_metrics`
  (per fold) and `overall` (pooled) beside the built-in ones. See **Your
  own metrics** on
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  for the contract.

## Value

A list with `overall`, `fold_metrics`, `predictions`, `folds`,
`n_folds_attempted`, `n_folds_succeeded`, `fold_status`, `orphan_rows`,
`n_unknown_ids`, `n_dropped`, `formula` and `adaptive`. The two fold
counts make a run where every fold failed visible in the return value
itself, beyond the warning, since `overall` is a well-formed all-`NA`
row either way, and `fold_status` (one row per fold: `fold`, `status`,
`message`) says why each missing fold is missing. See
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
for the five statuses and for `orphan_rows`. `Adj_R2` is `NA` in both
`overall` and `fold_metrics`: the pooled predictions have no single
parameter count, and a GWR's effective parameter count is not its
predictor count (see
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)).

## Details

Folds default to spatial blocks
([`make_folds`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)`(method = "block_kfold")`),
not random ones. With autocorrelated data a random split leaves a
held-out point's neighbours in the training set and the score comes back
flattering. Use
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
for the same treatment of a Bayesian GP model,
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for a forest, and
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
to score several backends on one set of folds.

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
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("GWmodel", quietly = TRUE) &&
    requireNamespace("sp", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 60
  dat <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
    2 * dat$elev + rnorm(n)
  cv <- cv_gwr(dat, "price", "elev", k = 3, bandwidth = 30)
  print(cv$overall)       # pooled over the held-out rows of every fold
  cv$fold_metrics         # and fold by fold
}
#> cv_gwr(): no folds supplied -- using spatial block k-fold CV (k=3).
#>       RMSE      MAE     MAPE    SMAPE       R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 1.970881 1.640318 11.65992 11.14015 0.620296     NA     60     60      60
#>   fold n_train n_test n_pred    RMSE      MAE      MAPE     SMAPE        R2
#> 1    1      39     21     21 1.34471 1.167010  7.532286  7.606738 0.6918865
#> 2    2      41     19     19 2.31357 1.985458 14.905985 13.676105 0.5117606
#> 3    3      40     20     20 2.16090 1.809408 12.910182 12.441061 0.6691330
#>   Adj_R2 n_MAPE n_SMAPE bandwidth
#> 1     NA     21      21        30
#> 2     NA     19      19        30
#> 3     NA     20      20        30
```
