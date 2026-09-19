# K-fold cross-validation for the Bayesian spatial model

Refits the Gaussian-process model of
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
on each training fold and scores it on the held-out fold. Beyond the
point-prediction metrics the other CV wrappers report, this one scores
the whole predictive *distribution*: `predictive_coverage` says what
fraction of held-out observations fell inside the 50/80/95\\ calibration
together. That is the reason to reach for it. A Bayesian model is
usually chosen for its uncertainty, and only held-out coverage shows
whether those intervals are honest at locations the model has not seen.

## Usage

``` r
cv_bayes(
  data_sf,
  response_var,
  predictor_vars,
  folds = NULL,
  k = 5,
  seed = 123,
  boundary = NULL,
  pointize = "auto",
  fit_args = list(),
  summary = c("mean", "median"),
  compute_pred_intervals = TRUE,
  coverage_levels = c(0.5, 0.8, 0.95),
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

  RNG seed. Default 123. It seeds fold construction **and**, through a
  per-fold draw, each fold's Stan sampler, so two seeds give different
  posteriors even on identical `folds`. A `seed` in `fit_args` overrides
  the per-fold draw with one fixed sampler seed.

- boundary:

  Optional polygonal sf/sfc for CRS alignment.

- pointize:

  Geometry coercion strategy.

- fit_args:

  Named list of extra arguments for fit_bayesian_spatial_model(). A
  user-supplied `gp_k` is respected in every fold; when omitted, the GP
  rank is auto-selected per training fold. `compute_loo`, `boundary`,
  and `pointize` are always overridden by the CV internals.

- summary:

  "mean" or "median" for posterior predictions.

- compute_pred_intervals:

  Logical; compute predictive intervals.

- coverage_levels:

  Numeric vector of coverage levels.

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
  1, use that many cores. Default `FALSE` (sequential). Bayesian folds
  with full MCMC runs are the primary beneficiary of this option.

- metrics:

  Optional scoring function of your own, a `function(y, yhat)` returning
  a named numeric vector; its names become columns of `fold_metrics`
  (per fold) and `overall` (pooled) beside the built-in ones. See **Your
  own metrics** on
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  for the contract. `yhat` is the posterior predictive `summary` (mean
  or median) of each held-out row; a score that needs the draws belongs
  in `predictive_coverage` and `CRPS`, which this function computes
  itself.

## Value

A list with `overall`, `fold_metrics`, `predictions`, `folds`,
`n_folds_attempted`, `n_folds_succeeded`, `fold_status`, `orphan_rows`,
`n_unknown_ids`, `n_dropped`, `formula` and `predictive_coverage`. The
two fold counts make a run where every fold failed visible in the return
value itself, beyond the warning, and `fold_status` (one row per fold:
`fold`, `status`, `message`) keeps the reason each missing fold is
missing, including the error text of a fold whose sampler failed, where
a long run's console output would not. See
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
for the five statuses and for `orphan_rows`. `predictions` carries,
beyond the columns its siblings share, `yhat_sd`: the posterior
predictive standard deviation of each held-out row, from the same draws
that give the coverage below (`NA` when `compute_pred_intervals = FALSE`
or the draws failed for that fold). `overall$Adj_R2` is always `NA`, as
for every `cv_*()`: see
[`cv_spatial`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
The `predictive_coverage` entries (one per `coverage_levels` value, plus
`mean_CRPS`) are averages across folds **weighted by each fold's
`n_pred`**, because the per-fold values in `fold_metrics` are themselves
means over that fold's test rows; an unweighted average would not be the
pooled quantity when fold sizes differ, which for spatially blocked
folds they routinely do.

## Details

It is the most expensive wrapper in the package by a wide margin: every
fold is a full MCMC run. Use few folds, and `parallel = TRUE` if you
have the cores. For a cheap first pass on the same question,
cross-validate a forest with
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
and come back here once the predictor set has settled.

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
`cv_bayes()` additionally reports CRPS and interval coverage, which are
proper scoring rules and have no such failure mode.

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

For the Bayesian backend, `cv_bayes()` additionally reports CRPS and
interval coverage at 50, 80 and 95 percent. Both are proper scoring
rules computed from posterior draws, so they are meaningful for any
`family` the backend accepts, and they are the numbers to compare when
the response is not Gaussian. When every fold fails, the `fold_metrics`
frame `cv_bayes()` returns carries the CRPS column but not the
`coverage_*` columns, so code that reads those columns must tolerate
their absence.

## See also

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run: fits with Stan, which needs a working C++ toolchain and takes
# minutes of MCMC -- both outside what an example may assume.
if (requireNamespace("brms", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 60
  dat <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
  cv <- cv_bayes(dat, "price", "elev", k = 2,
                 fit_args = list(chains = 2, iter = 500))
  cv$overall
  cv$predictive_coverage  # coverage at 50/80/95% plus mean CRPS
}
} # }
```
