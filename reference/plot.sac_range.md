# Plot an estimated spatial autocorrelation range

Draws the empirical variogram that
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
attaches to its result, with the fitted model and the effective range
overlaid where a range was identified, and a subtitle saying why not
where it was not. A variogram that never reaches a sill, or that has no
spatial structure at the lags resolved, is the single most useful thing
to *see* when a range comes back `NA`, so those cases are drawn rather
than refused.

## Usage

``` r
# S3 method for class 'sac_range'
plot(x, ...)
```

## Arguments

- x:

  An object of class `sac_range`, as returned by
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md).

- ...:

  Ignored.

## Value

A `ggplot` object.

## Details

Nothing is recomputed: the plot reads the `variogram`,
`variogram_model`, `crs`, `directional` and `anisotropy_used` attributes
the estimate already carries. The distance axis is in the units of the
CRS the variogram was actually fitted in, which
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
may have chosen itself for lon/lat input.

## See also

[`plot.spatial_fit`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)`(type = "variogram")`,
which draws the same picture for a fitted model's residuals.

Other plotting:
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md),
[`plot.block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.block_size_sweep.md),
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md),
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
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
  n <- 150
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  xy$z <- sin(xy$x / 150) + rnorm(n, sd = 0.3)
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  r <- estimate_sac_range(pts, response_var = "z")
  plot(r)
}
```
