# Cross-validate a Bayesian GP with spatial folds

Refits the Gaussian-process model of
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
on each training fold and scores it on the held-out fold. Besides the
point-prediction metrics the other `cv_*()` functions report, it scores
the whole predictive distribution: `predictive_coverage` gives the share
of held-out observations inside the 50, 80 and 95 percent intervals, and
`mean_CRPS` rates sharpness and calibration together. A Bayesian model
is often chosen for its uncertainty, and held-out coverage shows whether
those intervals are honest at locations the model has not seen.

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

  An sf layer with the response and the predictors.

- response_var:

  Name of the response column.

- predictor_vars:

  Character vector of predictor column names.

- folds:

  Optional fold definitions, in one of three shapes: a
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  result; a bare list of `list(train =, test =)` pairs of `..row_id`
  values; or a vector of fold labels, one per row, which becomes
  leave-that-label-out splits. The label vector is how folds built by
  another package are used here: `blockCV::cv_spatial()` returns one as
  `$folds_ids`, while its `$folds_list` holds two unnamed vectors per
  fold and is refused by name. Train and test must be disjoint; a fold
  that trains on its own test rows is not a cross-validation split and
  is refused with an error. IDs naming no row of the prepared data are
  dropped with a logged count (expected when rows were removed for
  missing values; a sign the folds came from other data when they were
  not). When `NULL`, the folds are `make_folds(method = "block_kfold")`,
  built with `k`, `seed`, `boundary`, `block_size` and `auto_range`.
  With `folds = NULL` the score estimates the error at the distance the
  block grid built here puts between the held-out points and their
  training data, which overstates the error of a map where the samples
  cover the area it predicts over. The excess is how much the model's
  error grows from the map's distances to the data to the block grid's,
  and is small for a model whose error changes little with the distance
  from the data.
  [`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md)
  with `prediction_points` compares the held-out distances with those of
  the locations the model will predict at, and
  `make_folds(method = "knndm", prediction_points = )` builds folds
  matched to them.

- k:

  Number of folds to build when `folds` is `NULL`; not read when `folds`
  is supplied. Default 5.

- seed:

  Seed for building the folds when `folds` is `NULL` and, through one
  seed drawn per fold, for each fold's Stan sampler, so two seeds give
  different posteriors even on identical `folds`. A `seed` in `fit_args`
  overrides the per-fold draw with one fixed sampler seed. Default 123.

- boundary:

  Optional polygon layer (sf or sfc). Data without a CRS are aligned to
  it, and the blocks built when `folds` is `NULL` are clipped to it (see
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
  and
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)).

- pointize:

  How a feature that is not a point is reduced to one before fitting: a
  mode of
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).
  Default `"auto"`.

- fit_args:

  Named list of extra arguments for
  [`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md).
  A `gp_k` given here is used in every fold; without one, the GP rank is
  chosen in each training fold. `compute_loo`, `boundary` and `pointize`
  are always set by the cross-validation. A categorical or ordinal
  `family`
  ([`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
  `cumulative`, `sratio`, `cratio`, `acat`) is refused before anything
  is fitted: its prediction is a probability per response category, and
  every score here needs one number per row.

- summary:

  `"mean"` (default) or `"median"`: the summary of the posterior
  predictive draws used as each held-out row's prediction.

- compute_pred_intervals:

  Logical, default `TRUE`: draw the posterior predictive distribution of
  the held-out rows. `FALSE` skips the draws and saves time, and leaves
  the `coverage_*` columns, `CRPS`, `n_draws`, `predictive_coverage` and
  `predictions$yhat_sd` `NA`.

- coverage_levels:

  Numeric vector of the nominal coverage levels to score, as proportions
  strictly between 0 and 1 (`0.95`, not `95`), each given once; anything
  else is an error. Each becomes a column `coverage_<percent>` of
  `fold_metrics`, named at full precision (`0.975` gives
  `coverage_97.5`).

- block_size:

  Optional minimum block edge length, in the units of the CRS the folds
  are built in, for the blocks built when `folds` is `NULL`; see
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md).
  It is a size you choose and is compared with the autocorrelation range
  only when one can be estimated, in which case a smaller size is warned
  about. Beside supplied `folds` it has nothing to size and raises a
  warning of class `spatialkit_fold_args_ignored`.

- auto_range:

  Logical, default `FALSE`. `TRUE` with `folds = NULL` estimates the
  autocorrelation range and uses it as the minimum block size. `TRUE`
  beside supplied `folds` raises the same warning as `block_size`.

- parallel:

  Logical or positive integer. `TRUE` detects the number of cores and
  fits the folds in parallel with
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (macOS and Linux; sequential on Windows); an integer above 1 uses that
  many cores. Default `FALSE` (sequential). Mind the memory: every fold
  compiles its own Stan model, and one compilation can take several GB,
  so `parallel = n` runs `n` of them at once. A compiler killed for lack
  of memory fails its fold with rstan's `"invalid connection"` error,
  which `fold_status` records; use fewer cores if you see it.

- metrics:

  Optional scoring function of your own, a `function(y, yhat)` returning
  a named numeric vector; its names become columns of `fold_metrics`
  (per fold) and `overall` (pooled) beside the built-in ones. See the
  section "Your own metrics" in
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  for the contract. `yhat` is the posterior predictive `summary` (mean
  or median) of each held-out row; a score that needs the draws belongs
  in `predictive_coverage` and `CRPS`, which this function computes
  itself.

## Value

A list with the elements
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
returns (`overall`, `fold_metrics`, `predictions`, `folds`,
`n_folds_attempted`, `n_folds_succeeded`, `fold_status`, `orphan_rows`,
`n_unknown_ids` and `n_dropped`; see its Value for each, and for the
five fold statuses), plus:

- `formula`:

  The model formula, as text.

- `predictive_coverage`:

  One entry per `coverage_levels` value, plus `mean_CRPS`: averages
  across folds weighted by each fold's `n_pred`. The per-fold values in
  `fold_metrics` are means over that fold's test rows, so an unweighted
  average would not be the pooled quantity when fold sizes differ, as
  they do for spatial blocks.

- `coverage_levels`:

  The nominal levels, named by their `coverage_*` column.

`fold_status` keeps the error text of a fold whose sampler failed.
`predictions` also carries `yhat_sd`, the posterior predictive standard
deviation of each held-out row, from the same draws that give the
coverage (`NA` when `compute_pred_intervals = FALSE` or the draws failed
for that fold). `overall$Adj_R2` is always `NA`, as for every `cv_*()`.

`fold_metrics` also carries `gp_k`, `gp_n_basis`, `n_draws`, `CRPS`, the
`coverage_*` columns, `convergence_ok` and `gp_adequacy`.
`convergence_ok` is `TRUE` or `FALSE` as
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
judged that fold's sampler (divergent transitions, R-hat, effective
sample size and E-BFMI; see its `check_convergence`), and `NA` when
`fit_args` sets `check_convergence = FALSE`. A fold that did not
converge is scored like the others, so a run with any `FALSE` raises one
warning of class `spatialkit_convergence` naming those folds and the
checks each failed. `gp_adequacy` is what that fold's posterior said
about its Gaussian-process basis (`"ok"`, `"below_resolution"`,
`"unresolved_short_range"`, `"no_signal_above_floor"` or
`"not_checked"`; see "What the basis resolves" in
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).
A fold whose basis was found too coarse is scored like the others too;
see the section "Warnings from the folds".

## Details

Every fold is a full MCMC run, so this is the slowest of the `cv_*()`
functions. Use few folds, and `parallel = TRUE` if you have the cores.
For a quicker first pass on the same question, cross-validate a random
forest with
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
and come back here once the predictor set has settled.

## Warnings from the folds

[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
warns about a fit under four condition classes:
`spatialkit_gp_basis_size`, `spatialkit_gp_resolution`,
`spatialkit_convergence` and `spatialkit_gp_adequacy` (see its section
"Warnings"). Each fold is a fit, and none of those warnings says which
fold it came from, so they are collected as the folds run and each class
is raised once after the last fold, under this function's name and with
the same class. Each warning names the folds and carries `folds` (their
labels) and `fold_messages` (every fold's own text, named by fold) as
fields of the condition. The convergence warning covers the folds that
were scored and gives the checks each failed, also as the field
`failed`; the other three quote the first fold's own text. Each fold's
own line is in the log. Warnings that Stan or brms raise themselves are
not collected and are passed on.

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
response can be zero. For a Bayesian GP fit, `cv_bayes()` also reports
CRPS and interval coverage, which are proper scoring rules and are
defined at zero.

## Which metrics survive a non-Gaussian response

RMSE and MAE are defined for any numeric response and are what to read
for a count, a rate or a bounded outcome. MAPE and SMAPE assume a
response that is rarely zero (see "Percentage errors on responses with
zeros"), and \\R^2\\ and adjusted \\R^2\\ compare residual variance to
total variance, which suits a Gaussian response and fits loosely
anything whose variance tracks its mean. On a Poisson or zero-inflated
response read the four as a rough summary rather than a score.

For the Bayesian GP, `cv_bayes()` also reports CRPS and interval
coverage at its `coverage_levels` (50, 80 and 95 percent by default).
Both are proper scoring rules computed from posterior draws, so they are
meaningful for any `family` that predicts one number per row (a count, a
rate, a binary or bounded outcome), and they are the numbers to compare
when the response is not Gaussian. A categorical or ordinal family
predicts a probability per response category instead, so `cv_bayes()`
refuses one before fitting anything. When every fold fails, the
`fold_metrics` frame `cv_bayes()` returns carries the `CRPS` column but
not the `coverage_*` columns, so code that reads those columns must
allow for their absence.

## See also

Other cross-validation:
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Not run: fits with Stan, which needs a working C++ toolchain and takes
# minutes of MCMC.
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
# Two short chains per fold keep this to a few minutes and leave the
# posterior rough, so the intervals below will run narrow; a run to
# report uses the defaults of brms (chains = 4, iter = 2000).
cv <- cv_bayes(dat, "price", "elev", k = 2,
               fit_args = list(chains = 2, iter = 1000))
cv$overall
cv$predictive_coverage  # coverage at 50, 80 and 95 percent, and mean CRPS
} # }
```
