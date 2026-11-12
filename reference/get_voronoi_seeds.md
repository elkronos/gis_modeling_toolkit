# Place seeds for a Voronoi tessellation

Creates an sf POINT layer of seed locations, one per cell of the Voronoi
tessellation built on them: k-means centres of a cloud of points, points
drawn uniformly inside a boundary, or points you provide.

## Usage

``` r
get_voronoi_seeds(
  boundary = NULL,
  method = c("kmeans", "random", "provided"),
  n = NULL,
  seeds = NULL,
  sample_points = NULL,
  kmeans_nstart = NULL,
  kmeans_iter = 100,
  set_seed = NULL,
  partition = c("density", "kmeans"),
  init = NULL,
  density_power = 1,
  density_k = 8
)
```

## Arguments

- boundary:

  Optional polygonal sf object defining the sampling area. Rows without
  a location, rows the repair leaves without a polygon, vertices at a
  latitude beyond 90 and mixed dimensions are handled as "Boundaries" in
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  says; "Sampling inside the boundary" gives how points are drawn in it.

- method:

  How to place the seeds: `"kmeans"` (default; cluster centres of
  `sample_points`, or of points drawn in `boundary`), `"random"`
  (uniform in `boundary`) or `"provided"` (the points in `seeds`).

- n:

  Integer; number of seeds to return. Required for `method = "kmeans"`
  and `method = "random"`. Ignored for `method = "provided"`, where
  every row of `seeds` is returned; a mismatch between `n` and
  `nrow(seeds)` is reported as a warning.

  For `method = "kmeans"` it is an upper bound only: k-means cannot
  produce more centres than there are distinct positions in the sampling
  cloud, nor as many centres as there are rows. When `n` exceeds either
  ceiling it is clamped, with a warning naming the count used. Asking
  for as many seeds as there are points (`n = nrow(sample_points)`)
  gives `nrow(sample_points) - 1`. Check
  [`nrow()`](https://rdrr.io/r/base/nrow.html) on the result; do not
  assume `n`.

  Besides a number, `n` accepts what the level-selection step returned:
  the integer vector of ranked candidates from
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  (its first element is used; its `NA` is refused with a note on where
  it came from), a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result (its `$best`), or a
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).
  A profile built without an `estimate` is read with
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  by `cp`. A profile built with `estimate` is read by
  `criterion = "auto"`, the rule that goes with its estimate: a size for
  pooled values and the lowest leave-one-out error for plain means. The
  size is the profile's own when it was built at one, and 4.5 rows per
  cell otherwise, and the count is the profile's level for that size, or
  the built level nearest to `max(2, round(n / obs_per_cell))` when the
  profile has none for it (with a warning of class
  `"spatialkit_size_not_reached"`). `select_resolution(profile)` with no
  `criterion` reads such a profile the same way, since its default
  criterion is `"auto"`. A profile built with `obs_per_cell` and without
  an `estimate` is read by `cp` like any profile without one; pass
  `select_resolution(profile, criterion = "size")` for the count its
  size gives. The result then carries `attr(, "n_from")` saying which. A
  count read off a selection or a profile comes with what the selection
  found: a warning of class `"spatialkit_bound_chosen"` when the optimum
  sits on a bound of the ladder, so the bound chose the count, and one
  of class `"spatialkit_no_structure"` when one cell for the whole layer
  scores as well. A count set by a size has no optimum, so no bound can
  choose it; its `"spatialkit_no_structure"` warning says that the count
  stands and that the cell values will say little that one mean does
  not. A profile without an estimate and with no `cp` to read is an
  error (class `"spatialkit_no_default_criterion"` when another
  criterion is finite); see
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md).

  With `method = "kmeans"`, a selection or profile also carries the
  centres of the partition the profile scored at that count, and those
  are returned instead of a new k-means: a k-means partition, weighted
  or not, is the Voronoi partition of its centres, so their cells are
  the cells the criteria judged. On a profile without an estimate whose
  centres were fitted to a subsample,
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  refits them to every point first, and those centres are returned, with
  `attr(, "refit")`: their cells are near the cells judged but not the
  same. `sample_points`, when given, is then only assigned to them (for
  `attr(, "kmeans")`) and gives the cover, and `boundary` only sets the
  CRS and the cover; `kmeans_nstart`, `kmeans_iter` and `set_seed` are
  not used, and neither are `partition`, `init` and the density
  settings: the centres are of the kind the profile was built with, and
  naming another here is a warning of class
  `"spatialkit_partition_unused"`. A call names what it supplies: a
  `partition` or an `init` passed in the call; with
  `partition = "density"` passed, both density settings at their values,
  defaults included; without it, a density setting other than its
  default. A call that supplies none of the four is told nothing,
  whichever kind and settings the profile was built with. Pass the count
  as a number (`n = sel$best`) for a fresh k-means instead.

- seeds:

  sf POINT object of user-provided seeds (method = "provided"). Methods
  `"kmeans"` and `"random"` do not read it, and ignore it with a
  warning; so does a call that passes `seed =`, which R matches to
  `seeds` (the seed of the random stream is `set_seed`). Seeds in
  another CRS than `boundary` are returned in the boundary's CRS; seeds
  that mix XY and XYZ (or XYM) points are then returned in XY, with a
  warning of class `"spatialkit_zm_dropped"`, since sf reads a layer in
  the dimension of its first geometry.

- sample_points:

  Optional sf POINT cloud for k-means clustering. Method `"random"`
  places its seeds in `boundary` and does not read it; it is ignored
  there with a warning. Only the first two coordinate columns are
  clustered, so a Z or M dimension does not join the distance
  calculation and dominate it; rows without a location (an empty point,
  a coordinate that is not finite or lies beyond 1e150, or in a
  geographic CRS a latitude beyond 90 degrees, read as given before any
  move) are dropped with a warning, so they never reach
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html), which fails
  on them without naming a cause. A lon/lat cloud, or one with no CRS
  whose coordinates look like lon/lat (the heuristic
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  applies, with its warning), is projected before clustering. A cloud
  that mixes XY and XYZ (or XYM) points is read without Z and M before
  it is projected or clustered, since sf reads a layer in the dimension
  of its first geometry.

- kmeans_nstart:

  The number of starts, of which the best is kept. `NULL` (default) uses
  3 starts for `partition = "density"` and 10 for
  `partition = "kmeans"`. Under `partition = "kmeans"` from random
  starts it is `nstart` of
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html). A fresh
  partition is not the one
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  scored a count on; pass the selection or profile itself as `n` for
  that one (see `n`). Under `partition = "density"`, and with
  `init = "kmeans++"`, it has to be `NULL` or a positive whole number,
  so a number that is not whole, which `partition = "kmeans"` takes,
  stops a call that names no kind. It is not read by `method = "random"`
  or `"provided"`.

- kmeans_iter:

  Integer; `iter.max` for
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html). Default 100
  for both kinds
  ([`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  uses 50 for plain k-means by default). Under `partition = "density"`
  it is the most steps of Lloyd's iteration a start takes. Under
  `partition = "density"`, and with `init = "kmeans++"`, it has to be a
  positive whole number, and a kept start that was cut off at it is a
  warning of class `"spatialkit_not_converged"` (`partition = "kmeans"`
  from random starts keeps the warning
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) gives).

- set_seed:

  Optional integer seed for the random draws. Default `NULL`: the draws
  come from the session's random-number stream, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes it reproducible. A number gives the same seeds every time and
  leaves the caller's stream as it was.
  ([`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
  defaults to 456 instead.)

- partition:

  For `method = "kmeans"`: what the seeds minimise. `"density"`
  (default) is the sum of squared distances from the points of the cloud
  to their nearest seed, each point weighted by its local density to the
  power `density_power`: compact cells, many small ones where the points
  are dense and a few large ones where they are sparse. At the default
  power they hold about the same number of points each once `n` is above
  the number of dense places (clusters) in the cloud; below that the
  counts are about as uneven as under `"kmeans"`. `"kmeans"` is the same
  sum unweighted: compact cells, with several times the points in a cell
  of a dense place that a cell of a sparse place holds. It is
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) from 10
  random starts on the rows as they come, with the coordinates centred,
  divided by the diagonal of their bounding box and rounded to 2^-30 of
  it, so that the seeds do not move with the unit or the origin; where
  the locations span less than a kilometre a change of unit or origin
  can still move them (see "The partition" in
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)).
  The density kind follows the density of `sample_points`; a cloud this
  function draws inside `boundary` is uniform and has none to follow.
  Without the FNN package the density kind finds its weights by
  comparing every pair of distinct locations, so its time grows with the
  square of their number; the seeds are the same, and nothing stops the
  call. See "The partition" in
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  for what each kind does and how.

- init:

  For `method = "kmeans"`: how each start is drawn, `"random"` (`n`
  distinct locations at random, as
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) draws them)
  or `"kmeans++"` (the seeds one at a time, each in proportion to its
  squared distance from those already drawn, times its weight under
  `partition = "density"`). `NULL` (default) is each kind's own rule:
  `"random"` for `partition = "kmeans"` and `"kmeans++"` for
  `partition = "density"`. Under `partition = "kmeans"` from random
  starts the cloud is clustered in the order its rows come in; with
  `init = "kmeans++"`, and under `partition = "density"`, the rows are
  put in coordinate order first, so the seeds do not depend on their
  order.

- density_power:

  For `partition = "density"`: the power the local density is raised to
  for a point's weight, from 0 to 2. Default 1, which evens out the cell
  counts once the cells outnumber the dense places; 0 weights every
  point alike. Setting it under `partition = "kmeans"` is a warning.

- density_k:

  For `partition = "density"`: the number of nearest neighbours the
  local density is read from. Default 8.

## Value

An sf POINT object with seed_id and method columns. With
`method = "kmeans"` it also carries `attr(, "kmeans")`, the run behind
the seeds: `cluster` (the `seed_id` each clustered cloud point was
assigned to), `rows` (those points' positions in the cloud, since rows
with unusable coordinates are dropped first), `size` (points per seed),
`withinss` and `tot_withinss` (the within-cluster sums of squares),
`iter` and `nstart`. For the scored centres of a selection or profile
(see `n`) it describes `sample_points` assigned to them, each point to
its nearest seed, with `iter` `NA` and `nstart` the profile's restarts,
and is absent when no `sample_points` were given. A point the same
distance from two seeds (squared distances within one part in ten
million) goes to the lower `seed_id`, the rule the profile scores its
cells by, so the record is the same in another unit or with another
origin. Under `partition = "density"` `withinss` and `tot_withinss`
count every point once, about its cell's mean. A run this call makes
with `partition = "density"` or `init = "kmeans++"` also records
`converged`, whether the kept start stopped before `kmeans_iter`, and
under `partition = "density"` `objective`, the weighted sum the seeds
minimise; scored centres have neither. Seeds placed with
`partition = "density"` or `init = "kmeans++"`, and the scored centres
of a profile of the density kind or from random starts, carry
`attr(, "partition")`: a list with `kind`, `init`, `density_power` and
`density_k`, and for seeds this call placed `iter`, `iter_max` and
`converged` (the iterations the kept start took, the limit, and whether
it stopped before the limit).
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
returns the list as `params$partition`. Seeds with no such record are
plain k-means: from random starts when this call placed them, or the
scored centres of a k-means profile from k-means++ starts. Every method
also attaches `attr(, "cover")` when the call says what ground the seeds
stand for (see "Details"): an sfc holding the boundary's polygons or the
hull of the clustered points, in its own CRS.

When `n` is a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result that carries `count_rule`, or a profile built with an `estimate`
or read by `cp`, the seeds carry that record as `attr(, "count_rule")`:
how the count was set, for which cell value, and what the profile scored
at it. Its `scored_cells` is `TRUE` for the scored centres
`method = "kmeans"` hands on, and `FALSE` for centres refitted to every
point (see `n`) and for seeds this call places itself at that count: a
k-means of its own when the selection holds no centres at the count, and
every `method = "random"` seeding. `method = "provided"` does not read
`n` and writes no record; seeds that already carry one (`sel$seeds`)
keep it, marked `FALSE` when their number is not the record's count.
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
hand the record on, and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
checks its result against it. A count given as a number, or read off a
profile without an estimate by `cp` or `elbow`, has no record.

## Details

Seeds lie inside the ground they stand for (k-means seeds inside the
hull of the points they were fitted to, random seeds inside the
boundary), so Voronoi cells clipped to the seeds' own hull would leave
its edge uncovered. The result therefore carries that ground as
`attr(, "cover")`: the `boundary` when one was given, otherwise the
convex hull of `sample_points` (for lon/lat points, the hull of their
coordinates in degrees, with the longitudes of points on both sides of
180 degrees taken from 0 to 360, -170 as 190), otherwise (for the scored
centres of a
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
or a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result passed as `n`) the hull of the points the profile was computed
on.
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
and
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
derive their default boundary from it, so
`build_tessellation(seeds, method = "voronoi")` holds every observation
without a `boundary`. Seeds with no cover are tessellated over their own
hull unless a `boundary` or `cover` is passed then:
`method = "provided"` seeds given with neither `boundary` nor
`sample_points` have none, and no warning is raised for them, since they
may be the observation sites themselves. The attribute survives
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
and row subsetting; [`rbind()`](https://rdrr.io/r/base/cbind.html) drops
it.

A k-means seed is a mean of the points of its cluster, weighted under
the density kind, so it can lie outside a `boundary` that is not convex,
in a bay or a bend of a coast, even when every point lies inside.
Tessellated within that boundary, its cell is built if any of it lies
inside: the part of the study area nearer the seed than any other seed.
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
gives the seed that cell in `index`, as it gives every seed the cell it
generates, so `seeds$cell_id <- tess$index` holds each seed's cell for
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md).
A seed with every part of the study area nearer another seed has no
cell, `NA` in `index`, and a warning of class `spatialkit_unassigned`
gives the number.

## Sampling inside the boundary

Seeds drawn by `method = "random"`, and the cloud `method = "kmeans"`
draws when no `sample_points` are given, are sampled inside `boundary`.
Its rows are repaired with
[`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html)
before they are unioned. In longitude and latitude the union runs with
s2, and a row s2 still refuses (crossing edges, a vertex listed twice)
is repaired again, with crossing edges split as GEOS splits them in a
projected CRS; a row s2 refuses even then stops the call with an error
of class `"spatialkit_cells_degenerate"` that names it and says to
transform `boundary` to a projected CRS.

In longitude and latitude, a boundary whose bounding box is over 10,000
times its area (small pieces far apart), which sf does not sample whole,
is sampled piece by piece, the points shared among the pieces in
proportion to their areas; a piece whose own box is over 10,000 times
its area stops the call with an error that says to transform `boundary`
to a projected CRS. A boundary whose bounding box has no height or no
width (three vertices on one parallel) is refused. In longitude and
latitude s2 reads the edges as great circles, which leave the bounding
box read from the coordinates, and sf (without lwgeom) draws in that
box: a boundary whose box has less than half its area on the sphere
(three vertices on or near one parallel) is refused too. That error, and
the error for a box with no height whose vertices lie on one parallel
other than the equator, say to transform `boundary` to a projected CRS
in which a parallel is a curve, such as the one
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
picks for it (a polar one, such as EPSG:3413 or 3031, for a ring around
a pole), where the region has an area and lies in its box.

sf draws the points in the bounding box and keeps those inside, so it
draws the number asked for times the ratio of the box to the area. In
any CRS, a boundary that would take more than 10^6 points drawn in its
box, or 10 times the points asked for if more (small pieces far apart, a
thin strip at a slant), is sampled piece by piece in the same way (or
whole, in draws of at most 10^6 points in its box, where that takes less
time), and a piece that would still take more is drawn in such draws.
Where the points drawn in all would exceed 2e7 (or 10 times the points
asked for, if more), the call stops before any draw with an error that
gives the ratio and says to pass `sample_points` (`method = "kmeans"`,
whose cloud is `max(2000, 50 * n)` points) or to ask for fewer seeds.

## See also

Other cell building:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
get_voronoi_seeds(bnd, method = "random", n = 5, set_seed = 1)
#> Simple feature collection with 5 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 20.16819 ymin: 6.178627 xmax: 90.82078 ymax: 94.46753
#> Projected CRS: WGS 84 / UTM zone 32N
#>   seed_id method                  geometry
#> 1       1 random POINT (26.55087 89.83897)
#> 2       2 random POINT (37.21239 94.46753)
#> 3       3 random POINT (57.28534 66.07978)
#> 4       4 random  POINT (90.82078 62.9114)
#> 5       5 random POINT (20.16819 6.178627)

# Seeds from a cloud with one dense place: density-weighted k-means (the
# default, 3 starts) against plain k-means.  The fewest and most points a
# seed's cell holds under each kind, with more seeds than dense places:
set.seed(2)
cloud <- st_as_sf(data.frame(x = c(rnorm(150, 30, 5), runif(50, 0, 100)),
                             y = c(rnorm(150, 30, 5), runif(50, 0, 100))),
                  coords = c("x", "y"), crs = 32632)
dens <- get_voronoi_seeds(bnd, method = "kmeans", n = 20,
                          sample_points = cloud, set_seed = 1)
km   <- get_voronoi_seeds(bnd, method = "kmeans", n = 20,
                          sample_points = cloud, set_seed = 1,
                          partition = "kmeans")
sizes <- list(density = attr(dens, "kmeans")$size,
              kmeans  = attr(km, "kmeans")$size)
sapply(sizes, range)   # fewest and most points per seed
#>      density kmeans
#> [1,]       6      1
#> [2,]      13     32
```
