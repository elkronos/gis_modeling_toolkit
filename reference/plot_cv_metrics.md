# Plot one cross-validation metric fold by fold

One pooled RMSE can come from the same error in every fold or from small
errors in most folds and a large one in a single fold (a model that
works, and a model that fails in one block), and the pooled number
cannot tell the two apart. This draws the metric of each fold as a
point, sized by the number of held-out predictions the fold contributed,
with the pooled value from `overall` as a horizontal line, so the spread
behind the number is visible. A
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
result draws one panel per model on a shared scale, folds aligned, which
is the comparison the shared fold set was built for.

## Usage

``` r
plot_cv_metrics(cv, metric = "RMSE", ...)
```

## Arguments

- cv:

  The list returned by
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
  [`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  [`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md)
  or
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md).

- metric:

  Character(1) naming a column of `fold_metrics` (or of `by_fold`).
  Default `"RMSE"`.

- ...:

  Not used.

## Value

A `ggplot` object.

## Details

Any column of `fold_metrics` can be drawn, including a backend's extras
(`bandwidth`, `CRPS`, `coverage_95`) and columns a `metrics` function
added. A column that is `NA` in every fold is refused with a message
saying why rather than drawn as an empty panel: `Adj_R2` is `NA` for
every backend unless `p` was passed to
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
by design. The pooled line is drawn from `overall` when it carries the
metric, from `predictive_coverage` for
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)'s
coverage and CRPS columns, and not at all for a per-fold extra that has
no pooled counterpart (a bandwidth) or for a count (`n_pred`, `n_MAPE`,
`n_SMAPE`), whose `overall` value is the total over the folds; the
caption says which. A model with no finite per-fold value gets no panel
and no pooled line, and the caption names it.

## See also

[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
and
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
which make the results.

Other cross-validation:
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_rf()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_rf.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 150
dat <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
dat$z <- 2 * dat$a + 0.003 * (st_coordinates(dat)[, 1] - 5e5) +
  rnorm(n, 0, 0.5)
cv <- cv_rf(dat, "z", "a", k = 5, num_trees = 100)
#> cv_rf(): no folds supplied -- using spatial block k-fold CV (k=5).
plot_cv_metrics(cv, "RMSE")
```
