# Cross-validate at a ladder of block sizes

A blocked cross-validation with blocks smaller than the autocorrelation
range leaks: every held-out point has a near-identical neighbour in the
training set, and the score is optimistic in proportion. The single
number a `cv_*()` call returns cannot show this. This runs the same
cross-validation at a ladder of block sizes, and by default once with
random folds as the fully leaky reference. It returns the metric at
each, with the estimated autocorrelation range alongside so that the
curve can be read against it: it rises as the blocks pass the range and
then plateaus, and the height of the rise is how much the random-fold
number overstated the model.

## Usage

``` r
cv_block_size_sweep(
  data_sf,
  response_var,
  predictor_vars,
  fit_fn,
  block_sizes = NULL,
  n_sizes = 6L,
  k = 5L,
  metric = "RMSE",
  include_random = TRUE,
  max_fits = 60L,
  sac = NULL,
  seed = 123L,
  quiet = FALSE,
  ...
)
```

## Arguments

- data_sf:

  An sf object with the response and predictors.

- response_var, predictor_vars:

  Column names.

- fit_fn:

  A `function(train_sf)` returning a `spatial_fit`, as for
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md);
  see the example for wrapping a built-in backend.

- block_sizes:

  Optional numeric vector of block edge lengths to sweep, in the CRS
  units the folds are built in. Default `NULL`: the ladder described
  above.

- n_sizes:

  Number of sizes in the default ladder. Default 6.

- k:

  Folds per cross-validation. Default 5.

- metric:

  Which column of `overall` to read. Default `"RMSE"`.

- include_random:

  Also cross-validate with random folds, as the leaky reference. Default
  `TRUE`.

- max_fits:

  The fit budget; see above. Default 60.

- sac:

  Optional `sac_range` from
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
  to mark on the curve. Default `NULL`: estimated here from the
  response, detrended on `predictor_vars`, when gstat is installed.

- seed:

  Seed for the fold construction at every size.

- quiet:

  Suppress the progress messages. Default `FALSE`.

- ...:

  Passed to
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  at every size (`predict_args`, `p`, `parallel`, `metrics`, ...).
  `block_size`, `folds`, `k`, `seed` and `auto_range` are set here and
  cannot be passed.

## Value

A data.frame of class `"block_size_sweep"` with one row per
cross-validation: `block_size` (`NA` for the random reference),
`method`, `blocks_used`, `k` (the folds actually built),
`n_folds_succeeded`, `value` (the pooled metric), `fold_min`, `fold_max`
and `fold_sd` (its spread across folds). Attributes: `metric`,
`sac_range` (the effective range, or `NA`), `crs`, `n_fits`, and
`results`, the full
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
result at every size.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it.

## The fit budget

Each block size is a full cross-validation, so the cost is
`length(block_sizes) * k` fits, plus `k` for the random reference.
`max_fits` caps that (default 60: six sizes at `k = 5`, plus the
reference). A sweep that would run past the cap refuses to start, naming
the number of fits it would have needed. Raise `max_fits` deliberately;
the RF example below takes seconds, a Bayesian `fit_fn` takes minutes
per fit.

## The ladder

When `block_sizes` is `NULL`, `n_sizes` values are log-spaced from a
twenty-fifth to a half of the shorter side of the data's extent, and any
size at which the grid would hold fewer than `k` blocks is dropped, so
every point on the curve is a `k`-fold cross-validation of the same
shape. Sizes are in the units of the CRS the folds are built in
([`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)'s
`params$crs`, metres for geographic input), and the returned table
records that CRS.

## See also

[`plot.block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.block_size_sweep.md).

Other cross-validation:
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gwr_model_selection.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)

## Examples

``` r
if (requireNamespace("ranger", quietly = TRUE) &&
    requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(dist(cbind(x, y)))
  field <- as.numeric(t(chol(exp(-d / 100) + diag(1e-6, n))) %*% rnorm(n))
  dat <- st_as_sf(data.frame(x = x, y = y, a = rnorm(n)), coords = c("x", "y"),
                  crs = 32632)
  dat$z <- field + 0.5 * dat$a + rnorm(n, 0, 0.2)
  rf_fn <- function(train_sf)
    fit_rf_model(train_sf, "z", "a", include_coords = TRUE, num_trees = 100, seed = 1)
  sw <- cv_block_size_sweep(dat, "z", "a", fit_fn = rf_fn, k = 4, n_sizes = 4,
                            quiet = TRUE)
  sw
  if (requireNamespace("ggplot2", quietly = TRUE)) plot(sw)
}
```
