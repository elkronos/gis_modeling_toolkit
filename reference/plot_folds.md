# Map a cross-validation fold scheme

Shows which fold each observation belongs to. This is the fastest way to
see whether spatial blocks are actually separating the data, or whether
the blocks are smaller than the autocorrelation range and therefore
leaking. For `"block_kfold"` folds the block outlines are drawn too,
from `folds$params$blocks`, so a fold can be seen to be one region or
several and an empty block can be seen to be empty.

## Usage

``` r
plot_folds(folds, points_sf, boundary = NULL, blocks = TRUE)
```

## Arguments

- folds:

  A list returned by
  [`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md).

- points_sf:

  The `sf` layer the folds were built from.

- boundary:

  Optional polygonal `sf`/`sfc` to draw underneath.

- blocks:

  Logical; draw the block polygons when `folds` carries them. Default
  `TRUE`.

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
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 80
  pts <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000)),
    coords = c("x", "y"), crs = 32632
  )
  f <- make_folds(pts, k = 5, method = "block_kfold", block_size = 300)
  plot_folds(f, pts)
}
```
