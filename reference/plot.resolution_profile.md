# Plot a resolution profile

The criteria of a
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
against the number of cells, one panel per criterion on a shared x axis,
with the level each criterion selects marked and the region over which
it is within `tol` of its optimum shaded. The shading marks the flat
region
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reports. A criterion whose optimum sits at an end of the levels it was
scored at is captioned with which bound that is (the support ceiling,
the range floor, the ladder's own end, or the first level the criterion
is computable at), because there the bound is choosing, not the
criterion. A criterion that is `NA` at every level is left out with a
log line rather than drawn empty.

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
  `"elbow"`, `"loo"`, `"moran_z"`, `"wss"`. Default: the first five,
  each where it is finite at some level. `"cp"`, `"reliability"`,
  `"elbow"` and `"loo"` (the leave-one-out error of the estimate the
  profile was built for; a profile built without `estimate` has no such
  column, and a column of that name added to one is drawn only when
  named here) are count criteria and get a selected level and a shaded
  band. `"moran_z"` is a diagnostic of the cell means and not a count
  rule (see
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)),
  and `"wss"` is the raw curve behind `"elbow"`: both are drawn as
  curves with no selected point.

- tol:

  Passed to
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  for the flat region.

- ...:

  Not used.

## Value

A `ggplot` object.

## Details

A profile whose ladder is one size (`obs_per_cell` set to one number, or
the default profile) has one level, set by the size. Each criterion that
is finite there is drawn as one point in its panel, with no selected
level and no band, and the caption says that the size set the count.

## See also

[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
which makes the object, and
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md).

Other resolution:
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
library(sf)
# The same field as ?resolution_profile: an exponential covariance with
# range parameter 200 and a nugget of 0.6 on a unit sill.
set.seed(4)
n <- 400
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
prof <- resolution_profile(pts, response_var = "z", n_levels = 12)
# Two of the criteria it scored (plot(prof) draws every one, a panel each)
plot(prof, criteria = c("cp", "reliability"))

plot(prof, criteria = c("cp", "wss")) # Cp beside the raw WSS curve
```
