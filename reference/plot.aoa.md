# Plot the dissimilarity distribution behind an area of applicability

`n_outside` says how many prediction locations fall outside the area of
applicability; it does not say whether the rest sit comfortably inside
or crowd against the threshold, nor how far outside the outsiders are.
This draws the dissimilarity index of the prediction locations against
that of the cross-validated training data, with the threshold marked, so
the prediction set can be read as mostly inside, marginal or largely
outside. The training curve is the reference the threshold was derived
from: its upper tail ends where the threshold is (or below it, when the
outlier fence removed the tail).

## Usage

``` r
# S3 method for class 'aoa'
plot(x, type = c("ecdf", "histogram"), ...)
```

## Arguments

- x:

  An `aoa` object from
  [`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md).

- type:

  `"ecdf"` (default), the two empirical distribution functions on one
  axis, or `"histogram"`, the prediction DI as bars with the training DI
  as an outline.

- ...:

  Ignored.

## Value

A `ggplot` object.

## See also

Other plotting:
[`plot.block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.block_size_sweep.md),
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md),
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md),
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md),
[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(sf)
  set.seed(2)
  n <- 200
  train <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  train$z <- train$a - train$b + rnorm(n, 0, 0.3)
  # Prediction locations whose predictor `a` drifts beyond the training range.
  new <- st_as_sf(
    data.frame(x = 5e5 + runif(100, 0, 1000), y = 5e6 + runif(100, 0, 1000),
               a = rnorm(100, mean = 2), b = rnorm(100)),
    coords = c("x", "y"), crs = 32632)
  aoa <- area_of_applicability(new, train_sf = train, predictor_vars = c("a", "b"))
  plot(aoa)
  plot(aoa, type = "histogram")
}
```
