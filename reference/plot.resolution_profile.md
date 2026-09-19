# Plot a resolution profile

The criteria of a
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
against the number of cells, one panel per criterion on a shared x axis,
with the level each criterion selects marked and the region over which
it is within `tol` of its optimum shaded. That shaded band is the flat
region
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reports. A criterion whose optimum sits at the support ceiling or the
range floor is captioned as such, because there the bound is choosing,
not the criterion.

## Usage

``` r
# S3 method for class 'resolution_profile'
plot(x, criteria = NULL, tol = 0.02, ...)
```

## Arguments

- x:

  A `resolution_profile`.

- criteria:

  Character vector of criteria to draw, any of `"cp"`, `"reliability"`,
  `"elbow"`, `"moran_z"`, `"wss"`. Default: every one of the four
  selectable criteria that is finite at some level. `"wss"` is the raw
  curve behind `"elbow"` and has no selected point.

- tol:

  Passed to
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  for the flat region.

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
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md),
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md),
[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(sf)
  set.seed(3)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  xy$z <- sin(xy$x / 200) + cos(xy$y / 250) + rnorm(n, sd = 0.3)
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, response_var = "z", n_levels = 10)
  plot(prof)
  plot(prof, criteria = c("cp", "wss"))
}
```
