# Compare every criterion's pick on a resolution profile

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
reports it: `FALSE` for a `cp` that one cell matches, and where the
cells beat one cell but the F test does not tell their means from one
mean or, for `cp`, the profile's variogram found no spatial structure).
Attributes: `bands` (a named list holding each criterion's flat region
in full, since `flat_min` and `flat_max` are only its ends and the
region can have holes in it), `common` (the levels in every flat region
of the count criteria, `moran_z` left out; an integer vector that is
empty when the regions do not overlap), `scored` (the levels each
criterion returned a finite value at, which is what decides whether a
gap in a flat region is a rejection or a level nothing was computed at),
`ladder` (every level on the profile), `n_levels`, `tol` and `variable`
(what the criteria were scored on), and, only when some row has
`structured` `FALSE` although its cells beat one cell by more than
`tol`, `not_identified` (a character vector named by criterion that says
why: the F test, or for `cp` the profile's variogram). The print method
recomputes the closing comparison from `bands`, so a row subset of the
result prints its own comparison.

## Details

Read the table with what each criterion returns in mind (see
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)):
`reliability` is read at or above the range floor and often returns it,
`elbow` is absent unless the points cluster, and `cp`'s optimum is often
the support ceiling. `moran_z` is a diagnostic of the cell means and not
a count rule, so it is not in the default set; named in `criteria` it
gets a row, and it is left out of the closing comparison.

The criterion curves are not monotone, so a flat region can skip rungs
of the ladder, and the table lists the levels it accepts (`"26, 31"`)
rather than a range (`"26 to 31"`) that would include the levels it
rejected.

The levels common to every flat region are often none. The criteria
answer different questions: how well the cells represent the field
(\\C_p\\), whether the cell values are distinguishable from noise
(`reliability`), and where the within-cluster sum of squares bends
(`elbow`). An empty set means the criteria favour different counts, and
the print gives the spread between the picks.

A flat region can be wide, and the choice within it is the analyst's;
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md)
draws the curves the regions were read from.

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
for the curves behind these flat regions.

Other resolution:
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)

## Examples

``` r
library(sf)
# An exponential field with range parameter 200 (true effective range
# 600 m) on a 1 km square, with a nugget of 0.6 on a unit sill.
set.seed(4)
n <- 400
xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
D  <- as.matrix(dist(xy))
xy$z <- as.numeric(t(chol(exp(-D / 200) + diag(0.6, n))) %*% rnorm(n))
pts <- st_as_sf(xy, coords = c("x", "y"), crs = 32632)
prof <- resolution_profile(pts, response_var = "z", n_levels = 12,
                           estimate = NULL)

s <- summary(prof)
s
#> Resolution picks: 2 criteria over 12 levels (2 to 44 cells)
#> 
#>   cells: density-weighted k-means++ starts (power 1, 8 neighbours)
#> 
#>    criterion best  flat region levels in flat region
#>           cp   19 19 to 25, 44                     3
#>  reliability    4       4 to 6                     3
#> 
#>   reliability: the optimum is the range floor (area / range^2).
#>   There the bound is choosing, not the criterion.
#> 
#>   picks span 4 to 19 cells (4.8x)
#>   no level is in every flat region: the criteria disagree over the
#>   whole ladder. plot() draws the curves they were read from.
attr(s, "common")   # the levels every flat region accepts, if any
#> integer(0)
attr(s, "bands")    # each criterion's flat region in full
#> $cp
#> [1] 19 25 44
#> 
#> $reliability
#> [1] 4 5 6
#> 
```
