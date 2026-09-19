# Model-agnostic spatial cross-validation

Run K-fold CV for any model that returns a `spatial_fit` object. This is
the extensibility point: to plug in a new model type, supply a
`fit_fn(train_sf)` that returns a `spatial_fit`.

## Usage

``` r
cv_spatial(
  data_sf,
  response_var,
  predictor_vars,
  fit_fn,
  folds = NULL,
  k = 5,
  seed = 123,
  boundary = NULL,
  pointize = "auto",
  predict_args = list(),
  fold_info_fn = NULL,
  p = NULL,
  block_size = NULL,
  auto_range = FALSE,
  parallel = FALSE,
  metrics = NULL,
  .caller = "cv_spatial"
)
```

## Arguments

- data_sf:

  An sf object.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- fit_fn:

  A function of one argument, the training slice of `data_sf`, returning
  a `spatial_fit` built with
  [`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md).
  It is called once per fold on the training rows only, so anything done
  inside it (scaling, tuning, an inner variable sweep) is already nested
  and leak-free. The `subclass` it stamps must have a
  `predict.<subclass>()` method registered, because that is how each
  fold is scored.

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
  in the prepared data are dropped with a logged count. Built via
  `block_kfold` when `NULL`.

- k:

  Number of folds.

- seed:

  RNG seed.

- boundary:

  Optional boundary for fold construction.

- pointize:

  Geometry coercion strategy.

- predict_args:

  Extra arguments for predict().

- fold_info_fn:

  Optional `function(fit, test_sf, y, yhat)` returning a named list of
  per-fold extras (a bandwidth, a tuning value, anything read off the
  fitted object), added as columns of `fold_metrics`. It sees the fit
  and the held-out layer, which `metrics` does not; it is applied per
  fold only, and its values are not pooled. An element `..per_row` that
  is a data frame with one row per held-out observation is spliced into
  `predictions` instead.

- p:

  Number of predictors for Adj R² (NULL to skip). Only meaningful for
  models with a fixed global parameter count; pass NULL for models with
  spatially varying coefficients (e.g. GWR).

- block_size:

  Optional minimum block edge length for spatial CV blocks (projected
  CRS units).

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

  Optional scoring function of your own; see **Your own metrics** below.
  Default `NULL`: the built-in metrics only.

- .caller:

  Internal. The name the messages carry, so a wrapper such as
  [`cv_rf`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  reports itself in place of `cv_spatial()`.

## Value

A list with `overall`, `fold_metrics`, `predictions`, `folds`, the two
fold counts `n_folds_attempted` and `n_folds_succeeded`, and four
elements that say what happened to the difference between them and to
the rows: `fold_status`, `orphan_rows`, `n_unknown_ids` and `n_dropped`.
The counts are reported deliberately: a `fit_fn` that fails on every
fold otherwise looks like a successful run that happened to score `NA`,
so compare them before trusting `overall`. `fold_status` is a data.frame
with one row per fold supplied (`fold`, `status` and `message`), where
`status` is `"ok"`; `"error"` (the fit or its
[`predict()`](https://rdrr.io/r/stats/predict.html) threw; `message` is
the error text); `"skipped"` (nothing scorable: too few matched rows, a
prediction of the wrong length, or no finite observed/predicted pair);
`"dropped"` (an empty test set, or fewer than two training rows, once
incomplete rows were removed, so the fold never reached the fitter); or
`"worker_error"` (a parallel worker died). Every fold missing from
`fold_metrics` has its reason there, which matters most when the console
output of a long run is gone. `orphan_rows` holds the `..row_id`s of
rows in the data that no fold names (they enter no training set and are
never scored; non-empty only when the folds were built on a different or
subsetted layer), and `n_unknown_ids` counts the distinct row IDs the
folds name that the data does not have. Each such row is named by every
fold, once as a test row and once in each other fold's training set, and
this counts the row, not the mentions (expected when rows were removed
for missing values). `n_dropped` is how many rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry before any
fold was fitted. The `fold` column of `fold_metrics`, `predictions` and
`fold_status` carries the fold's index in the `folds` object that was
supplied, so it lines up with `make_folds()$assignment$fold` even when
some folds were unusable and dropped. `overall$Adj_R2` is always `NA`:
the pooled out-of-sample predictions come from `k` separately fitted
models and have no single parameter count to adjust for. The per-fold
`fold_metrics$Adj_R2` carries the adjusted value when `p` is supplied,
and is `NA` otherwise.

## Name collision with blockCV

The blockCV package exports a function of the same name that does the
opposite job: `blockCV::cv_spatial()` *builds* spatial folds, where this
function *runs* a cross-validation over folds it is given. With both
packages attached, whichever was attached last masks the other;
`spatialkit::cv_spatial()` always resolves to this one. The two
cooperate: `blockCV::cv_spatial()` returns its fold assignment as
`$folds_ids`, a vector of fold labels, and that vector is accepted
directly as the `folds` argument here and in every other `cv_*()`
function. Fold construction is blockCV's home ground and this package
does not try to match its breadth there; what this package adds is the
path from irregular points through data-drawn regions to a
cross-validated, compared model.

## Your own metrics

The built-in columns (`RMSE`, `MAE`, `MAPE`, `SMAPE`, `R2`, `Adj_R2`)
are the Gaussian regression set, and the section below says which of
them survive a count or a bounded response. `metrics` is the way to
score what they cannot: a `function(y, yhat)` that returns a named
numeric vector (a named list of scalars, or a one-row data frame, also
serve), for example a Poisson deviance, a log score on a probability, or
a loss with your own weights. It is applied twice, in the same way the
built-in metrics are: once per fold, to that fold's held-out rows, so
each name becomes a column of `fold_metrics`; and once to the pooled
out-of-sample predictions of every fold, so each name becomes a column
of `overall`. Only the pairs the built-in metrics use reach the function
(both `y` and `yhat` finite), so its columns describe the same rows as
`RMSE`, and `n_pred` counts them.

The contract: every element named, names unique and not one of the
built-in column names, one number per name. Anything else is an error,
because a scoring function that returns the wrong shape is a mistake to
surface instead of a fold to skip. A function that *throws* on a fold is
logged and its columns are `NA` for that fold (and for `overall`, if it
throws on the pooled predictions); a fold is never dropped for it. When
no fold produced a prediction the empty `fold_metrics` frame still
carries the function's columns, typed, provided the function can be
called on zero-length input.

`fold_info_fn` is the per-fold half of the same mechanism, with access
to the fitted object and the held-out layer; `metrics` sees only the two
vectors but is also pooled.
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
hands one `metrics` to every backend, so the columns are comparable
across the rows of its `overall`.

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

[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)
for the constructor a `fit_fn` must use;
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
and
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
for the built-in backends, which are thin wrappers over this function.

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 80
site <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
site$price <- 10 + 0.01 * (st_coordinates(site)[, 1] - 5e5) +
  2 * site$elev + rnorm(n)

# 1. A fit_fn returning a spatial_fit of your own subclass.
lm_fit <- function(train_sf) {
  new_spatial_fit(
    subclass       = "lm_fit",
    engine         = lm(price ~ elev, st_drop_geometry(train_sf)),
    formula        = price ~ elev,
    response_var   = "price",
    predictor_vars = "elev",
    data_sf        = train_sf
  )
}

# 2. The predict() method cv_spatial() scores each fold with. Without it
#    every fold fails and `overall` comes back all-NA.
predict.lm_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data_sf
  as.numeric(stats::predict(object$engine, st_drop_geometry(newdata)))
}
registerS3method("predict", "lm_fit", predict.lm_fit)

cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1)
#> cv_spatial(): no folds supplied -- using spatial block k-fold CV (k=3).
cv$overall
#>       RMSE      MAE     MAPE    SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 3.472328 3.009688 21.29953 20.27029 0.2341257     NA     80     80      80
# Compare these before trusting the metrics above; fold_status says why
# any fold is missing from fold_metrics.
c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
#> attempted succeeded 
#>         3         3 
cv$fold_status
#>   fold status message
#> 1    1     ok        
#> 2    2     ok        
#> 3    3     ok        

# 3. A metric of your own beside the built-in ones: per fold and pooled.
med_ae <- function(y, yhat) c(MedAE = stats::median(abs(y - yhat)))
cv2 <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1,
                  metrics = med_ae)
#> cv_spatial(): no folds supplied -- using spatial block k-fold CV (k=3).
cv2$overall$MedAE
#> [1] 2.981173
cv2$fold_metrics[, c("fold", "RMSE", "MedAE")]
#>   fold     RMSE    MedAE
#> 1    1 3.793601 3.536960
#> 2    2 3.198481 2.590833
#> 3    3 3.359840 2.678648
```
