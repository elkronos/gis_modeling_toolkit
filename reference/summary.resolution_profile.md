# Every criterion's pick, side by side

[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reads one criterion at a time. This puts the count criteria in one
table: the level each prefers, the flat region around it, and whether a
ladder bound is doing the choosing rather than the criterion. The
closing line gives the levels that lie in *every* flat region, the cell
counts no criterion objects to.

## Usage

``` r
# S3 method for class 'resolution_profile'
summary(object, criteria = NULL, tol = 0.02, ...)
```

## Arguments

- object:

  A
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

- criteria:

  Character vector, any of `"cp"`, `"reliability"`, `"elbow"`,
  `"moran_z"`, `"loo"`. Default: every one of `"cp"`, `"reliability"`,
  `"elbow"` and `"loo"` that is finite at some level (`loo` is a column
  of a profile built with `estimate`; a column of that name on a profile
  built without one is read only when named here). `"moran_z"`, a
  diagnostic, is read only when named here. `"size"` and `"auto"` are
  not criteria of the table: a size sets the count and has no flat
  region to compare.

- tol:

  Passed to
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  for the flat region. Default 0.02.

- ...:

  Ignored.

## Value

A data.frame of class `resolution_summary`, one row per criterion in the
order given, with columns `criterion`, `best`, `flat_min`, `flat_max`,
`n_flat`, `value` (the optimum itself, on that criterion's own scale and
so not comparable across rows, which is why the print method leaves it
out), `at_floor`, `at_ceiling`, `edge` (which bound an edge optimum sits
on, as
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reports it; `NA` when interior) and `structured` (as
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
reports it: `FALSE` for a `cp` that one cell matches). Attributes:
`bands` (a named list holding each criterion's flat region in full,
since `flat_min` and `flat_max` are only its ends and the region can
have holes in it), `common` (the levels in every flat region of the
count criteria, `moran_z` left out; an integer vector that is empty when
the regions do not overlap), `scored` (the levels each criterion
returned a finite value at, which is what decides whether a gap in a
band is a rejection or a level nothing was computed at), `ladder` (every
level on the profile), `n_levels`, `tol` and `variable` (what the
criteria were scored on). The print method recomputes the closing
comparison from `bands`, so a row subset of the result reads honestly.

## Details

The criteria are not peers, and the table should be read with what each
returns in practice in mind (see
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)):
`reliability` is read at or above the range floor and usually returns
it, `elbow` is absent unless the points cluster, and `cp`'s optimum is
often the support ceiling. `moran_z` is a diagnostic of the cell means
and not a count rule, so it is not in the default set; named in
`criteria` it gets a row, and it is left out of the closing comparison.

A flat region is a set, not an interval. The criterion curves are not
monotone, so a region can skip a rung of the ladder, and the table
prints what the criterion actually accepts (`"26, 31"`, not
`"26 to 31"`) rather than a range that would quietly include the levels
it rejected.

That set is often empty, and an empty one is a result rather than a
failure. The criteria answer different questions: how well the cells
represent the field (\\C_p\\), whether the cell values are
distinguishable from noise (`reliability`), and where the within-cluster
sum of squares bends (`elbow`). A field with no single right resolution
shows up here as disjoint bands, and the spread between the picks is
printed for the same reason.

Nothing in the table is a decision procedure. Each flat region is
routinely wide, the choice within it belongs to the analyst, and
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md)
draws the curves the bands were read from.

A profile whose ladder is one size (`obs_per_cell` set to one number, or
the default profile) has one level, and the size set it. The table then
has one row per criterion that is finite at that level, and the line
under it says that the count is set by the size and that no criterion
chooses among levels. The result carries the size as
`attr(, "obs_per_cell")`.

## See also

[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
for one criterion, with the per-level values attached;
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md)
for the curves behind these bands.

Other aggregation:
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`kriging_adequacy()`](https://elkronos.github.io/gis_modeling_toolkit/reference/kriging_adequacy.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)

## Examples

``` r
if (requireNamespace("gstat", quietly = TRUE)) {
  library(sf)
  # An exponential field with range parameter 200 (true effective range
  # 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
  set.seed(4)
  n <- 400
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  D  <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
  pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  prof <- resolution_profile(pts, response_var = "z", n_levels = 12)

  # print() is explicit because only the last value of a braced block is
  # shown, and the table is the thing worth seeing here.
  print(summary(prof))
  print(attr(summary(prof), "common")) # the levels all of them accept, if any
  attr(summary(prof), "bands")         # each criterion's region in full
}
#> Resolution picks: 3 criteria over 12 levels (2 to 44 cells)
#> 
#>   cells: density-weighted k-means++ starts (power 1, 8 neighbours)
#> 
#>    criterion best flat region levels in band
#>           cp   44          44              1
#>  reliability    4      4 to 6              3
#>          loo   44          44              1
#> 
#>   cp, loo: the optimum is the support ceiling (n / min_cell_n).
#>   reliability: the optimum is the range floor (area / range^2).
#>   There the bound is choosing, not the criterion.
#> 
#>   picks span 4 to 44 cells (11.0x)
#>   no level is in every flat region: the criteria disagree over the
#>   whole ladder. plot() draws the curves they were read from.
#> integer(0)
#> $cp
#> [1] 44
#> 
#> $reliability
#> [1] 4 5 6
#> 
#> $loo
#> [1] 44
#> 
```
