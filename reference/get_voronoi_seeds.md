# Generate seed points for Voronoi tessellation

Creates an sf POINT layer of "seed" locations. Multiple strategies are
supported: user-provided points, uniform random sampling within a
boundary, or k-means clustering of a sampling cloud.

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

  Optional polygonal sf object defining the sampling area.

- method:

  One of "kmeans", "random", "provided".

- n:

  Integer; number of seeds to return. Required for `method = "kmeans"`
  and `method = "random"`. **Ignored** for `method = "provided"`, where
  every row of `seeds` is returned; a mismatch between `n` and
  `nrow(seeds)` is reported as a warning.

  For `method = "kmeans"` it is an upper bound only: k-means cannot
  produce more centres than there are distinct positions in the sampling
  cloud, nor as many centres as there are rows. When `n` exceeds either
  ceiling it is clamped, with a warning naming the count actually used.
  `n = nrow(sample_points)` is the common case, and yields
  `nrow(sample_points) - 1` seeds. Check
  [`nrow()`](https://rdrr.io/r/base/nrow.html) on the result; do not
  assume `n`.

  Besides a number, `n` accepts what the level-selection step returned:
  the integer vector of ranked candidates from
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  (its first element is used), a
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
  the cells the criteria judged. `sample_points`, when given, is then
  only assigned to them (for `attr(, "kmeans")`) and gives the cover,
  and `boundary` only sets the CRS and the cover (see "Details");
  `kmeans_nstart`, `kmeans_iter` and `set_seed` are not used, and
  neither are `partition`, `init` and the density settings: the centres
  are of the kind the profile was built with, and naming another here is
  a warning of class `"spatialkit_partition_unused"`. A call names what
  it supplies: a `partition` or an `init` passed in the call; with
  `partition = "density"` passed, both density settings at their values,
  defaults included; without it, a density setting other than its
  default. A call that supplies none of the four is told nothing,
  whichever kind and settings the profile was built with. Pass the count
  as a number (`n = sel$best`) for a fresh k-means instead.

- seeds:

  sf POINT object of user-provided seeds (method = "provided").

- sample_points:

  Optional sf POINT cloud for k-means clustering. Only the first two
  coordinate columns are clustered, so a Z or M dimension does not join
  the distance calculation and dominate it; rows with empty or
  non-finite coordinates are dropped with a warning, so they never reach
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html), which fails
  on them without naming a cause. A lon/lat cloud, or one with no CRS
  whose coordinates look like lon/lat (the heuristic
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  applies, with its warning), is projected before clustering.

- kmeans_nstart:

  The number of starts, of which the best is kept. `NULL` (default) is
  the kind's own number: 3 for `partition = "density"` and 10 for
  `partition = "kmeans"`, the default of earlier versions.
  `kmeans_nstart = 10` gives the density kind the starts it had when 10
  was the default. Under `partition = "kmeans"` from random starts it is
  `nstart` of [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).
  A fresh partition is not the one
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  scored a count on; pass the selection or profile itself as `n` for
  that one (see `n`). Under `partition = "density"`, and with
  `init = "kmeans++"`, it has to be `NULL` or a positive whole number,
  so a number that is not whole, which `partition = "kmeans"` takes,
  stops a call that names no kind. It is not read by `method = "random"`
  or `"provided"`.

- kmeans_iter:

  Integer; iter.max for kmeans(). Default 100. Under
  `partition = "density"` it is the most steps of Lloyd's iteration a
  start takes. Under `partition = "density"`, and with
  `init = "kmeans++"`, it has to be a positive whole number, and a kept
  start that was cut off at it is a warning of class
  `"spatialkit_not_converged"` (`partition = "kmeans"` from random
  starts keeps the warning
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) gives).

- set_seed:

  Optional integer RNG seed.

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
  of a dense place that a cell of a sparse place holds. It is the run of
  earlier versions, whose default it was:
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) from 10
  random starts on the rows as they come. The density kind follows the
  density of `sample_points`; a cloud this function draws inside
  `boundary` is uniform and has none to follow. Without the 'FNN'
  package the density kind finds its weights by comparing every pair of
  distinct locations, so its time grows with the square of their number;
  the seeds are the same, and nothing stops the call. See "The
  partition" in
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  for what each kind does and how.

- init:

  For `method = "kmeans"`: how each start is drawn, `"random"` (`n`
  distinct locations at random, as
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) draws them)
  or `"kmeans++"` (the seeds one at a time, each in proportion to its
  squared distance from those already drawn, times its weight under
  `partition = "density"`). `NULL` (default) is each kind's own rule:
  `"random"` for `partition = "kmeans"`, the run this function has
  always made, and `"kmeans++"` for `partition = "density"`. Under
  `partition = "kmeans"` from random starts the cloud is clustered in
  the order its rows come in; with `init = "kmeans++"`, and under
  `partition = "density"`, the rows are put in coordinate order first,
  so the seeds do not depend on their order.

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
and is absent when no `sample_points` were given. Under
`partition = "density"` `withinss` and `tot_withinss` count every point
once, about its cell's mean. A run this call makes with
`partition = "density"` or `init = "kmeans++"` also records `converged`,
whether the kept start stopped before `kmeans_iter`, and under
`partition = "density"` `objective`, the weighted sum the seeds
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
result that carries `count_rule`, or a profile built with an `estimate`,
the seeds carry that record as `attr(, "count_rule")`: how the count was
set, for which cell value, and what the profile scored at it. Its
`scored_cells` is `TRUE` for the scored centres `method = "kmeans"`
hands on, and `FALSE` for seeds this call places itself at that count: a
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
convex hull of `sample_points`, otherwise (for the scored centres of a
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
inside, but the seed itself lies in no cell:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
warns that it falls in no cell (class `spatialkit_unassigned`) and gives
it `NA` in `index`. `seeds` in
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md)
says how to give such a seed its cell's ID.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
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

# Seeds from a cloud with one dense place: density-weighted k-means by
# default, from 3 starts; partition = "kmeans" gives the k-means run of
# earlier versions.  The points each seed's cell holds, under each kind,
# with more seeds than dense places:
set.seed(2)
cloud <- st_as_sf(data.frame(x = c(rnorm(150, 30, 5), runif(50, 0, 100)),
                             y = c(rnorm(150, 30, 5), runif(50, 0, 100))),
                  coords = c("x", "y"), crs = 32632)
dens <- get_voronoi_seeds(bnd, method = "kmeans", n = 20,
                          sample_points = cloud, set_seed = 1)
km   <- get_voronoi_seeds(bnd, method = "kmeans", n = 20,
                          sample_points = cloud, set_seed = 1,
                          partition = "kmeans")
rbind(density = sort(attr(dens, "kmeans")$size),
      kmeans  = sort(attr(km, "kmeans")$size))   # points per seed
#>         [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#> density    6    8    8    8    8    9    9    9   10    10    11    11    11
#> kmeans     1    2    2    3    5    5    5    6    6     7     7     8     9
#>         [,14] [,15] [,16] [,17] [,18] [,19] [,20]
#> density    11    11    11    12    12    12    13
#> kmeans     10    10    14    18    24    26    32
```
