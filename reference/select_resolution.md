# Choose a cell count from a resolution profile

Picks the level a criterion prefers, together with the *flat region*:
every level whose criterion value is within `tol` of the optimum. Report
the flat region with the count. On the default profile, one size scored
for pooled values, the default criterion is that size: it returns the
count the size gives, with no optimum and no flat region (see "A profile
built with an estimate"); what follows is about the criteria read on a
ladder of counts. When the optimum sits at an end of the levels the
criterion was scored at, the result says so and names the bound, because
a bound is then doing the choosing rather than the criterion (see
[`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
for what each criterion measures and what it returns in practice). The
functions that take the result as a cell count raise that as a warning
of class `spatialkit_bound_chosen`.

## Usage

``` r
select_resolution(
  profile,
  criterion = c("auto", "cp", "reliability", "elbow", "moran_z", "loo", "size"),
  tol = 0.02,
  range_floor = TRUE,
  diagnostic = FALSE,
  obs_per_cell = NULL
)
```

## Arguments

- profile:

  A `resolution_profile`.

- criterion:

  Which rule sets the count: `"cp"` (minimised), `"reliability"`
  (maximised at or above the range floor), `"elbow"` (maximised),
  `"loo"` (the leave-one-out error of the profile's estimate,
  minimised), `"size"` (the count `obs_per_cell` gives), `"auto"`
  (`"size"` for a profile of pooled values, `"loo"` for one of plain
  means, `"cp"` for a profile without an estimate) or, with
  `diagnostic = TRUE`, `"moran_z"` (\\\|z\|\\ minimised). Default
  `"auto"`.

- tol:

  Width of the flat region. For `cp`, `reliability` and `loo` it is
  relative to the optimum's value (`0.02` keeps levels within 2 percent
  of it); for `elbow` and `moran_z`, whose optimum can be zero, it is
  relative to the criterion's range over the ladder. Under `"size"`
  there is no flat region, and `tol` is the margin by which the cells
  have to beat one mean for `structured` and nothing else.

- range_floor:

  For `criterion = "reliability"` only. `TRUE` (default) reads the
  criterion at the levels at or above the range floor in
  `attr(profile, "bounds")`. It applies when the layer supports the
  floor, some scored level reaches it, and the profile built its own
  ladder (a ladder set with `levels` is read as given). `FALSE` reads
  every level.

- diagnostic:

  `criterion = "moran_z"` is refused, with an error of class
  `spatialkit_diagnostic_criterion` that says what the column measures,
  unless this is `TRUE`: `moran_z` is a diagnostic of the cell means,
  not a count rule. Default `FALSE`.

- obs_per_cell:

  For `criterion = "size"`: the mean number of rows per cell wanted, a
  single number of at least 1 (fractions are allowed). `NULL` (default)
  is the profile's own size when it was built at exactly one
  ([`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)`(obs_per_cell = )`),
  and 4.5 otherwise. Given with a criterion that is not `"size"` it is
  an error.

## Value

A list of class `resolution_selection` with these elements:

- `best`:

  The level.

- `flat`:

  The levels in the flat region, ascending: a set, which can skip a
  rung.

- `criterion`:

  The criterion applied, never `"auto"`.

- `value`:

  The optimum.

- `at_ceiling`, `at_floor`:

  Whether the optimum is the last or first of the levels this criterion
  was read at, which for `moran_z` starts above nine cells and for
  `reliability` at the range floor; both `NA` under `"size"`, which
  searches for no optimum.

- `edge`:

  Which bound that is, in words: the support ceiling, the subsample's
  ceiling, the range floor, the ladder's own end, the finest or the
  coarsest size offered on a ladder `obs_per_cell` set, the only level
  scored, or the first or last level the criterion is computable at;
  `NA` for an interior optimum and under `"size"`.

- `n_levels`:

  The levels on the profile.

- `values`:

  The criterion at every level, `NA` where it could not be computed.

- `structured`:

  For `cp`: `FALSE` when `cp` at one cell is not above the optimum by
  more than `tol` (below it included), when the cell means at the
  optimum differ from one mean by no more than chance allows, or when
  the variogram found no spatial structure; `TRUE` when the optimum
  beats one cell on both counts (see "What each criterion returns"). For
  `loo` and `size` on a profile built with `estimate`: `FALSE` when the
  leave-one-out error of one mean for all scored rows is within `tol` of
  the error at `best`, or below it, or when the cell means at `best`
  differ from one mean by no more than chance allows; `TRUE` when the
  cells beat one mean on both counts. `TRUE` for an `elbow`; `NA` for
  `reliability` and `moran_z`, and for a profile that carries no
  one-cell value.

- `one_cell`:

  `cp` at one cell under `cp`; the leave-one-out error of one mean for
  all scored rows under `loo` and `size`; `NA` otherwise.

- `structure_p`:

  The p-value of the F test of the cells at `best` against one mean,
  multiplied by the number of levels searched and capped at 1 (a size
  searches one level); `NA` where `structured` is not read off `cp` or
  the leave-one-out error.

- `range_floor`:

  The floor `reliability` was held to, or `NA`.

- `note`:

  What the flags mean, in words: a character vector, one note an
  element; `NULL` when there is nothing to say. The notes written under
  `"loo"` and `"size"`, and for a profile built with `estimate` or with
  one size, carry names: `size` (the size that set the count), `held` (a
  size the profile held to the most cells k-means can place),
  `structure` (why the count is not identified), `finest` (the lowest
  error of pooled values is at the finest level scored), `estimate`
  (`cp` scores plain means, on a profile of pooled values) and `ladder`
  (the ladder is one size). A note of `cp` or `moran_z` about the
  criterion itself has no name. The warning of class
  `spatialkit_no_structure` that the functions taking the selection as a
  count raise carries every note but those named `size`, `held`,
  `finest`, `estimate` and `ladder`.

- `seeds`:

  An sf POINT layer, with `seed_id`, of the centres of the partition the
  profile scored at `best`, in the CRS it was computed in; `NULL` for a
  profile that does not carry them. See "Seeds" below.

- `partition`:

  For a profile built with `partition = "density"` or `init = "random"`:
  `kind`, `init`, `density_power` and `density_k`. The seeds carry it as
  `attr(, "partition")`; the cells are those the profile built, whatever
  kind built them. It stays the last element.

For `"loo"` and `"size"`, and for any criterion on a profile built with
`estimate`, the list has these elements after `seeds`:

- `estimate`:

  The profile's, `"mean"` or `"pooled"`; `NA` for a profile without one.

- `obs_per_cell`:

  Rows per cell at `best`: the rows a size divides over `best`; `NA` for
  an object that keeps no row count.

- `obs_per_cell_asked`:

  The size asked for; `NA` unless the criterion is `"size"`.

- `scored_per_cell`:

  The profile's column at `best`: the rows the leave-one-out error read
  over the cells that hold one. It is not `obs_per_cell`: it leaves out
  rows the error did not read and cells that hold none of them, as under
  `select_on = "split"`.

- `loo`, `loo_one_cell`:

  The leave-one-out error at `best` and that of one mean for all scored
  rows.

- `share_alone`, `own_weight`:

  The profile's columns at `best`.

- `count_rule`:

  A list of class `spatialkit_count_rule` that records how the count was
  set: `rule` (the criterion applied), `estimate` (the estimate the rule
  was computed for: the profile's for `"size"` and `"loo"`, `"mean"` for
  `"cp"`, `NA` otherwise), `cells`, `obs_per_cell`,
  `obs_per_cell_asked`, `response`, `n` (the rows the size divides),
  `n_points`, `n_fit` (the rows the centres were fitted to; fewer than
  `n_points` for a subsample), `seed`, `features`, `scored` (what the
  profile scored at `best` for its estimate: `estimate`, `n`, `cells`,
  `loo`, `loo_se`, `loo_one_cell`, `share_alone`, `own_weight`,
  `structured`, `rows` and `centres`, the matrix of the centres the rows
  were scored in, in `crs` and in the profile's order; `NULL` without an
  estimate), `partition`, `crs` and `scored_cells` (`NA`: the record
  stands for a count).

A `"cp"` selection off a profile without an estimate has `count_rule` as
well, after `seeds`, and nothing else of this list: its `estimate` is
`"mean"`, since \\C_p\\ scores plain means, and `scored` is `NULL`.
Under `"loo"` and `"size"` `one_cell` holds the leave-one-out error of
one mean and `range_floor` is `NA`.

`seeds` carries the same record as `attr(, "count_rule")` with
`scored_cells = TRUE`, since they are the centres that were scored, or
`FALSE` when they were refitted to every point, since they are then
other cells built from the count. The record travels on with the seeds,
the cells and the assigned layer to
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
which warns when it reports another cell value than the count was chosen
for (class `spatialkit_estimate_mismatch`) and when the cells or rows it
is given are not the ones the count was scored on
(`spatialkit_cells_differ`).

Seeds: the Voronoi cells of the seeds are the cells that were scored;
[`get_voronoi_seeds`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)`(method = "kmeans", n = <this>)`
and
[`voronoi_seeds_kmeans`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)`(k = <this>)`
return them. On a profile without an estimate whose centres were fitted
to a subsample (`n_sample` below `n` in `attr(profile, "bounds")`) the
seeds are those centres refitted to every point: the kind's own
iteration (Hartigan and Wong's for `"kmeans"`, Lloyd's on every point's
density weights for `"density"`) run from the centres scored, with the
profile's iteration limit and no random draw. Where Hartigan and Wong's
algorithm refuses the centres, Lloyd's iteration with every weight one
runs in its place, and a log line says so. Such seeds carry
`attr(, "refit")`: `rows` (the points), `subsample` (the points the
profile fitted to), `kind`, `iter`, `iter_max`, `converged`,
`wss_scored` (the within-cell sum of squares of every point in the cells
of the centres scored, each cell about its own mean) and `wss` (the same
of the refitted cells). Their cells are then not the cells scored, and
the count is still the one the subsample gave; the centres scored are
`attr(profile, "centres")[[as.character(best)]]`. The refit reads every
point's coordinates from `attr(profile, "refit_rows")`; on a profile
without that attribute the centres scored are the seeds, and a log line
says so. The refit is not run either, and the centres scored are the
seeds, when it would pass the work limits of a profile's own fits (see
"Conditions" in
[`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md));
a log line says so. A profile with an estimate scores it on every row of
the layer in the cells of the centres fitted to the subsample, so its
seeds are those centres, whose cells are the cells scored.
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a profile
read no seeds and refit nothing.

## What each criterion returns

- `cp`:

  The level that minimises Mallows' \\C_p\\: the count at which
  re-measuring the sampled sites, each predicted by its cell's mean, has
  the least estimated error
  ([`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  says what that leaves out for a new location). On a smooth response
  with little noise that is the support ceiling. The result also says
  whether the optimum beats one cell for the whole layer
  (`attr(profile, "one_cell")`), in `structured`. It is `TRUE` when two
  things hold: \\C_p\\ at one cell is above the optimum by more than
  `tol`, and the cell means at the optimum differ from one mean by more
  than chance allows. The second is the F test of the cells against one
  mean, \\(TSS - RSS) / (L - 1)\\ over \\RSS / (m - L)\\, at 5 percent
  after a Bonferroni adjustment for the levels the optimum was searched
  over (`structure_p` is the adjusted p-value). The cells are fitted to
  the coordinates alone, so the test holds its level on a response with
  no spatial structure whatever \\C_p\\'s noise variance is, and it is
  what decides when that variance is 0: \\C_p\\ is then \\RSS / n\\,
  which falls with every finer level by chance alone, and the first
  comparison says nothing. When either fails, `structured` is `FALSE`:
  one mean for the whole layer predicts as well as any level scored, so
  the response shows no spatial structure the cells can use and the
  count is not identified. `structured` is also `FALSE` when the
  profile's variogram was refused for showing no spatial structure.
  `best` is still the argmin. The F distribution is exact for a Gaussian
  response. The test has little power on a small layer, so a `FALSE`
  there is weak evidence that there is no structure.

- `reliability`:

  The level of highest reliability at or above the range floor the
  profile reports (when the layer supports the floor and the profile
  built its own ladder). The reliability curve often peaks below the
  floor, and the floor is a level of every such ladder, so the result is
  often the floor, the coarsest count the range allows, with `edge`
  naming it (on a profile whose k-means failed at that level, the first
  level above it). `range_floor = FALSE` reads every level and returns
  the curve's own maximum, which can be a count below the floor.

- `elbow`:

  The number of clusters of the points, when they cluster; an error
  otherwise, because the column is `NA` at every level unless the WSS
  curve bends.

- `moran_z`:

  Refused unless `diagnostic = TRUE`. It is the level at which
  neighbouring cell means (or, with predictors, their regression
  residuals) are least autocorrelated. \\\|z\|\\ is near 0 at 10 to 20
  cells whatever the data and grows as cells get finer than the range,
  so its minimum is usually the first level it can be computed at, or
  one near it, and it says nothing about how well the cells represent
  the response.

- `loo`:

  The level with the lowest leave-one-out error of the estimate the
  profile was built for
  ([`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)`(estimate = )`);
  of levels with the same error, the one with the fewest cells. The flat
  region is every level whose error is within `tol` of the lowest. The
  column needs no variogram. `structured` is read as for `cp`, with the
  leave-one-out error in place of \\C_p\\: it is `TRUE` when the error
  of one mean for all scored rows is above the optimum by more than
  `tol` and the F test of the cells at the optimum against one mean,
  adjusted for the levels searched, has \\p \< 0.05\\; `one_cell` then
  holds the error of one mean. Under this criterion the count depends on
  the values of the response. When the profile's estimate is the pooled
  value and the lowest error is at the finest level scored, a note says
  that the end of the ladder chose the count. For pooled values the
  error tends to keep falling to the finest level (see "How this was
  checked" in
  [`resolution_profile`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)).
  On a ladder `obs_per_cell` set, an optimum at an end has the `edge`
  `"the finest size offered (obs_per_cell)"` or
  `"the coarsest size offered (obs_per_cell)"`; a finest level the
  profile held to the most cells k-means can place keeps the label of
  the support ceiling.

- `size`:

  The count a size gives, `max(2, round(n / obs_per_cell))`, with `n`
  the rows the size divides: the rows with a response, or every point
  for a profile without one. No optimum is searched for. The level
  returned is the first of three that applies: the level the profile
  built for that size; the count itself, when the profile has a
  partition at it; otherwise the built level nearest to the count on the
  log scale (of two as near, the one with fewer cells), with a warning
  of class `spatialkit_size_not_reached`, which says so when the count
  is more than k-means can place on the layer. It reads any profile, one
  built without a response or without an estimate included, but not a
  column subset, which has lost the row count. `flat` is the level,
  `value` the rows per cell it leaves (`n / best`), and `at_ceiling`,
  `at_floor` and `edge` are `NA`: no bound chose anything, so the
  functions that take the selection as a count raise no
  `spatialkit_bound_chosen`. Under this rule neither the cells, which
  are fitted to the coordinates, nor their number, which is a count of
  rows, depend on the values of the response. A size does not ask
  whether the response has structure the cells can use, so `structured`
  asks it separately, on a profile built with an `estimate`. It is
  `FALSE` when the leave-one-out error of one mean for all scored rows
  is within `tol` of the error at the level, or below it, or when the F
  test of the cells at the level against one mean has \\p \ge 0.05\\.
  The count stands either way, since a size set it and no criterion:
  `structured = FALSE` says that the cell values will say little that
  one mean for the layer does not. The functions that take the selection
  as a count then raise `spatialkit_no_structure`. On a profile without
  an estimate `structured` is `NA`.

- `auto`:

  The rule that goes with the profile's estimate: `"size"` for a profile
  built with `estimate = "pooled"`, `"loo"` for `estimate = "mean"`, and
  `"cp"` for a profile without an estimate. A column subset of a profile
  has lost the record of its estimate and is read by `"cp"`. `criterion`
  in the result is the criterion applied, never `"auto"`.

## A profile built with an estimate

This function's default criterion is `"auto"`: `"size"` for a profile of
pooled values, `"loo"` for one of plain means and `"cp"` for a profile
without an estimate. A profile built with `estimate` that is passed
straight to
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
or
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
as the count is read the same way, so `select_resolution(profile)` and
the profile itself give the same count. A profile of pooled values built
with a ladder argument (`n_levels` and the others) is read by the size
as well: the default returns the level nearest the count the size gives,
on the automatic ladder its top level, with
`spatialkit_size_not_reached`. `criterion = "cp"` reads \\C_p\\ on a
profile that has it. `cp` is computed on such a profile as on any other
and scores plain means; when the profile's estimate is the pooled value
the result's note says so.

`"cp"` and `"reliability"` cannot be read off a profile whose ladder
`obs_per_cell` set under `sac = "auto"`. Such a profile fits no
variogram, both columns are `NA` at every level, and the call is an
error that names the ways to get them: build the profile with
`sac = NULL`, with a `sac` of your own, or with `obs_per_cell = NULL`.
The default profile is such a profile. On a ladder of one size that has
a variogram either criterion returns the one level, with `edge`
`"the only level scored"`.

For a ladder of one size,
[`summary()`](https://rdrr.io/r/base/summary.html) of the profile lists
each criterion that is finite at the level and says that the size set
the count, and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws each as one point, with no selected level and no flat region. When
none is finite, as on a constant response, where `loo` is `NA` too, or
when no row of the selection half has a response
(`select_on = "split"`),
[`summary()`](https://rdrr.io/r/base/summary.html) stops with an error
that says so. In that last case `"cp"` and `"moran_z"`, and `"loo"` with
an `estimate`, are errors that give the same reason, except `"cp"` on a
size ladder that fitted no variogram, which gives the size ladder's
reason, as above.

## See also

Other resolution:
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`summary.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.resolution_profile.md)

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

# The default profile is one size scored for pooled values, and the
# default criterion, "auto", reads it by that size.
select_resolution(resolution_profile(pts, response_var = "z"))
#> Resolution by size: 89 cells
#>   size        : 4.5 rows per cell, the default (obs_per_cell); 400 rows / 89 cells = 4.49
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
#>   estimate    : pooled toward the 8 nearest cells
#>   at this size: 1.0% of the rows are alone in their cell; median own weight 0.40
#>   loo         : 0.903 (one cell for the layer: 1.23)
#>   note        : the count is set by obs_per_cell = 4.5, the default; no
#>                 criterion chose it.

# A ladder of counts for plain means, read by Mallows' Cp.
prof <- resolution_profile(pts, response_var = "z", n_levels = 12,
                           estimate = NULL)
sel <- select_resolution(prof, criterion = "cp")
sel               # the level, and the flat region around it
#> Resolution by cp: 19 cells
#>   flat region : 19 to 25, 44 (3 of 12 levels)
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
sel$flat          # every level within `tol` of the optimum
#> [1] 19 25 44
sel$edge          # the bound the optimum sits on; NA when interior
#> [1] NA
sel$structured    # FALSE when one cell for the layer does as well
#> [1] TRUE

# Reliability prefers coarse cells: it is read at or above the floor the
# autocorrelation range sets, and the result says in words when the
# floor is choosing.
rel <- select_resolution(prof, criterion = "reliability")
rel
#> Resolution by reliability: 4 cells
#>   flat region : 4 to 6 (3 of 12 levels)
#>   cells       : density-weighted k-means++ starts (power 1, 8 neighbours)
#>   note        : the optimum is the range floor (area / range^2); the bound is
#>                 choosing, not the criterion. Fewer cells would be wider than
#>                 the range and average over more than one patch of the field.
rel$at_floor
#> [1] TRUE
rel$edge
#> [1] "the range floor (area / range^2)"
```
