# Plot the interval calibration of a Bayesian cross-validation

[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
scores every fold's posterior predictive intervals at each of
`coverage_levels`: the share of held-out observations the interval of
that nominal level contained. One coverage number cannot show a pattern;
the pairs can. This draws observed coverage against nominal with the
diagonal, one point per level for the fold-weighted pooled value (from
`predictive_coverage`) and one faint point per fold, so systematic
over-confidence (points below the line) or intervals wider than they
need to be (above it) are read at a glance. Three levels is a thin
curve; pass `coverage_levels = seq(0.1, 0.9, by = 0.1)` to
[`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md)
for a full one — the levels are read off the column names, so whatever
was computed is drawn.

## Usage

``` r
plot_calibration(cv, ...)
```

## Arguments

- cv:

  The list returned by
  [`cv_bayes()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_bayes.md),
  or a
  [`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md)
  result that ran the Bayesian backend (its `$bayes_cv` is used).

- ...:

  Ignored.

## Value

A `ggplot` object.

## See also

Other plotting:
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md),
[`plot.block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.block_size_sweep.md),
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md),
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md),
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md),
[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)

## Examples

``` r
# \donttest{
if (requireNamespace("brms", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  # cv <- cv_bayes(dat, "z", "a", k = 3, coverage_levels = seq(0.1, 0.9, 0.2))
  # plot_calibration(cv)
}
#> NULL
# }
```
