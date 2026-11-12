# Plot an estimated spatial autocorrelation range

Draws the empirical variogram that
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
attaches to its result, as a ggplot2 object, with the fitted model and
the effective range overlaid where a range was identified, and a
subtitle that gives the reason where it was not. A variogram that never
reaches a sill, or that has no spatial structure at the lags resolved,
is drawn as well, since the curve shows why the range came back `NA`.

## Usage

``` r
# S3 method for class 'sac_range'
plot(x, ...)
```

## Arguments

- x:

  A `sac_range` object, as returned by
  [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md).

- ...:

  Ignored.

## Value

A `ggplot` object.

## Details

Nothing is recomputed: the plot reads the `variogram`,
`variogram_model`, `crs` and `rejected_reason` attributes the estimate
carries. The variogram is the one over all point pairs. The distance
axis is in the units of the CRS the variogram was fitted in, which
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md)
chooses itself for lon/lat input.

## See also

[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
which makes the object, and
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)
with `type = "variogram"`, which draws the same picture for a fitted
model's residuals.

Other spatial autocorrelation:
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`sac_nugget()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sac_nugget.md)

## Examples

``` r
library(sf)
# An exponential field with range parameter 150 (effective range about
# 450 m) and a nugget of half the sill, so the fitted model, the nugget and
# the range line all have something to show.
set.seed(3)
n <- 250
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 150) + diag(0.5, n))) %*% rnorm(n))
pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
r <- estimate_sac_range(pts, response_var = "z")
plot(r)


# A field whose variogram never reaches a sill: the plot still draws it,
# and the subtitle says why no range is marked.
xy$trend <- sin(xy$x / 400) + rnorm(n, sd = 0.2)
r2 <- estimate_sac_range(st_as_sf(xy, coords = c("x", "y"), crs = 32632),
                         response_var = "trend")
plot(r2)
```
