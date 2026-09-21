# Read a level, and the region over which it is not distinguishable, off a profile

Picks the level a criterion prefers, together with the *flat region*:
every level whose criterion value is within `tol` of the optimum. On the
criteria this package computes the flat region is routinely wide: the
reliability curve is flat to within 2 percent over a factor of 3–6 in
the number of cells, and \\C_p\\ on a smooth field descends to the
support ceiling. The region is the answer, and the argmin only a point
in it. When the optimum sits at an end of the levels the criterion was
scored at, the result says so and names the bound, because a bound is
then doing the choosing rather than the criterion (see
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
(the levels in the flat region, ascending; a set, which can skip a
rung), `criterion`, `value` (the optimum), `at_ceiling` and `at_floor`
(logical: the optimum is the last or first of the levels this criterion
was scored at, which for `moran_z` starts above nine cells), `edge`
(which bound that is, in words: the support ceiling, the range floor,
the ladder's own end, or the first or last level the criterion is
computable at; `NA` for an interior optimum), `n_levels` and `values`
(the criterion at every level, `NA` where it could not be computed).

## See also

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  # An exponential field with range parameter 200 (true effective range
  # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill: enough
  # noise for Mallows' Cp to have an interior optimum rather than descend
  # to the ceiling.
  set.seed(2)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, response_var = "z", n_levels = 12)

  sel <- select_resolution(prof, criterion = "cp")
  print(sel)               # the level, and the flat region around it
  print(sel$flat)          # every level within `tol` of the optimum
  print(sel$edge)          # NA here: the optimum is interior

  # Reliability prefers coarse cells and here runs into the floor the
  # autocorrelation range sets, which the result says in words.
  rel <- select_resolution(prof, criterion = "reliability")
  print(rel)
  c(at_floor = rel$at_floor, edge = rel$edge)
}
#> Resolution by cp: 37 cells
#>   flat region : 31 to 44 (3 of 12 levels)
#> [1] 31 37 44
#> [1] NA
#> Resolution by reliability: 6 cells
#>   flat region : 6 to 9 (3 of 12 levels)
#>   note        : the optimum is the range floor (area / range^2); the bound is
#>                 choosing, not the criterion. Fewer cells would be wider than
#>                 the range and average over more than one patch of the field.
#>                           at_floor                               edge 
#>                             "TRUE" "the range floor (area / range^2)" 
```
