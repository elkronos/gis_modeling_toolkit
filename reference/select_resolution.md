# Read a level, and the region over which it is not distinguishable, off a profile

Picks the level a criterion prefers, together with the *flat region*:
every level whose criterion value is within `tol` of the optimum. On the
criteria this package computes the flat region is routinely wide: the
reliability curve is flat to within 2 percent over a factor of 3–6 in
the number of cells, and \\C_p\\ on a smooth field descends to the
support ceiling. The region is the answer, and the argmin only a point
in it. When the optimum sits at the ladder's ceiling or floor the result
says so, because a bound is then doing the choosing rather than the
criterion (see
[`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
for what each criterion measures and how it behaved on simulated
fields).

## Usage

``` r
select_resolution(
  profile,
  criterion = c("cp", "reliability", "elbow", "moran_z"),
  tol = 0.02
)
```

## Arguments

- profile:

  A `resolution_profile`.

- criterion:

  Which column decides: `"cp"` (minimised), `"reliability"` (maximised),
  `"elbow"` (maximised) or `"moran_z"` (\\\|z\|\\ minimised).

- tol:

  Width of the flat region. For `cp` and `reliability` it is relative to
  the optimum's value (`0.02` keeps levels within 2 percent of it); for
  `elbow` and `moran_z`, whose optimum can be zero, it is relative to
  the criterion's range over the ladder.

## Value

A list of class `resolution_selection` with `best` (the level), `flat`
(the levels in the flat region, ascending), `criterion`, `value` (the
optimum), `at_ceiling` and `at_floor` (logical: the optimum is the last
or first level of the ladder), `n_levels` and `values` (the criterion at
every level).

## See also

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  set.seed(2)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 100) + diag(0.3, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, response_var = "z", n_levels = 12)

  sel <- select_resolution(prof, criterion = "reliability")
  sel                      # the level, and the flat region around it
  sel$flat                 # every level within `tol` of the optimum
  sel$at_ceiling           # TRUE would mean the ladder, not the criterion, chose

  # A different criterion can prefer a different level while agreeing on the
  # region: the flat region is the answer, the argmin a point in it.
  select_resolution(prof, criterion = "cp")$flat
}
#> [1] 37 39 41
```
