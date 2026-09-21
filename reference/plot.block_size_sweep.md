# Plot cross-validation error against block size

The curve from
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md):
the pooled metric at each block size, the fold-to-fold range as a band,
the random-fold reference as a dashed line, and the estimated
autocorrelation range as a vertical marker. Blocks smaller than the
range leak, so the curve rises from the reference towards the range and
plateaus beyond it; the height of the rise is what the random-fold
number overstated.

## Usage

``` r
# S3 method for class 'block_size_sweep'
plot(x, ...)
```

## Arguments

- x:

  A `block_size_sweep`.

- ...:

  Ignored.

## Value

A `ggplot` object.

## See also

Other plotting:
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md),
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
if (requireNamespace("ranger", quietly = TRUE) &&
    requireNamespace("gstat", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 200
  x <- 5e5 + runif(n, 0, 1000); y <- 5e6 + runif(n, 0, 1000)
  d <- as.matrix(dist(cbind(x, y)))
  field <- as.numeric(t(chol(exp(-d / 100) + diag(1e-6, n))) %*% rnorm(n))
  dat <- st_as_sf(data.frame(x = x, y = y, a = rnorm(n)), coords = c("x", "y"),
                  crs = 32632)
  dat$z <- field + 0.5 * dat$a + rnorm(n, 0, 0.2)
  rf_fn <- function(train_sf)
    fit_rf_model(train_sf, "z", "a", include_coords = TRUE, num_trees = 100,
                 seed = 1)
  sw <- cv_block_size_sweep(dat, "z", "a", fit_fn = rf_fn, k = 4, n_sizes = 4,
                            quiet = TRUE)
  # Error rises from the random-fold reference (dashed) as the blocks grow
  # towards the estimated autocorrelation range (vertical marker).  The
  # height of that rise is what random folds were hiding; with four sizes the
  # ladder stops near the range, so the plateau beyond it is not drawn here.
  plot(sw)
}
```
