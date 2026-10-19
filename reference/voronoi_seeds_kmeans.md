# K-means seed generation from point coordinates

Places `k` seed points at k-means cluster centres of the observed
coordinates, so seeds (and the Voronoi cells built from them) follow the
sampling density: clusters of observations attract seeds, empty ground
gets none. Reach for this when you want cells that follow the data, so
that counts per cell vary far less than on a fixed grid over clustered
points (the density kind, the default, evens them out once there are
more cells than clusters of points; plain k-means,
`partition = "kmeans"`, minimises the spread of points around each
centre and does not equalise them), which is what keeps per-cell
aggregates in
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
from resting on one or two observations. Use
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
instead when you want coverage of the study area rather than of the
data, and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
to pick between them by name.

## Usage

``` r
voronoi_seeds_kmeans(
  points_sf,
  k,
  set_seed = 456,
  nstart = NULL,
  partition = c("density", "kmeans"),
  init = NULL,
  density_power = 1,
  density_k = 8,
  iter_max = NULL
)
```

## Arguments

- points_sf:

  An sf object with POINT geometries.

- k:

  Integer; requested number of clusters, treated as an upper bound only.
  It is clamped, with a warning, to whichever is smaller of the number
  of distinct point positions and `nrow(points_sf) - 1`, because k-means
  can produce neither more centres than there are distinct points nor as
  many centres as there are rows. Check
  [`nrow()`](https://rdrr.io/r/base/nrow.html) on the result. A
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result or a
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  is also accepted: its scored centres are returned when it carries
  them, whatever `partition`, `init` and density settings the call names
  (see `n` in
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  for which arguments a call is taken to name), and otherwise its count
  is used. A profile is read as
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  reads one passed as `n`: by `cp` when it was built without an
  `estimate`, and by `select_resolution(criterion = "auto")` when it was
  built with one (a size for pooled values, the lowest leave-one-out
  error for plain means). Such a count raises the warnings and the error
  described under `n` in
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  (classes `"spatialkit_bound_chosen"`, `"spatialkit_no_structure"` and
  `"spatialkit_no_default_criterion"`).

- set_seed:

  Optional integer RNG seed. Default 456, so a call gives the same seeds
  every time whatever the session's random-number state; an outer
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not change
  them, and the caller's random-number stream is left as it was. Pass
  `NULL` to draw the k-means starts from the session's stream instead
  (the default of
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).

- nstart:

  The number of starts, of which the best is kept. `NULL` (default) is
  the kind's own number: 3 for `partition = "density"` and 10 for
  `partition = "kmeans"`, the default of earlier versions. A number is
  run under either kind and either `init`: `nstart = 10` gives the
  density kind the starts it had when 10 was the default.

- partition:

  What the seeds minimise. `"density"` (default): the sum of squared
  distances from the points to their nearest seed, each point weighted
  by its local density to the power `density_power`; compact cells that,
  at the default power, hold about the same number of points each once
  `k` is above the number of dense places, and counts about as uneven as
  `"kmeans"` gives below that. `"kmeans"`: the same sum unweighted;
  compact cells whose counts are several times higher in dense places
  than in sparse ones. It is the run of earlier versions, whose default
  it was: [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) from
  10 random starts on the rows as they come. Without the 'FNN' package
  the density kind finds its weights by comparing every pair of distinct
  locations, so its time grows with the square of their number; the
  seeds are the same, and nothing stops the call. See "Details".

- init:

  How each start is drawn: `"random"` (`k` distinct locations at random,
  as [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) draws
  them) or `"kmeans++"` (the seeds one at a time, each in proportion to
  its squared distance from those already drawn, times its weight under
  `partition = "density"`). `NULL` (default) is each kind's own rule:
  `"random"` for `partition = "kmeans"`, the run this function has
  always made, and `"kmeans++"` for `partition = "density"`. Under
  `partition = "kmeans"` from random starts the points are clustered in
  the order their rows come in; with `init = "kmeans++"`, and under
  `partition = "density"`, the rows are put in coordinate order first,
  so the seeds do not depend on their order.

- density_power:

  For `partition = "density"`: the power the local density is raised to
  for a point's weight, from 0 to 2. Default 1, which evens out the cell
  counts once the cells outnumber the dense places; 0 weights every
  point alike, and values above 1 give the cells of dense places fewer
  points than those of sparse places. Setting it under
  `partition = "kmeans"` is a warning.

- density_k:

  For `partition = "density"`: the number of nearest neighbours the
  local density is read from. Default 8.

- iter_max:

  The most iterations a start runs for: a single positive whole number,
  or `NULL` (default) for the kind's own limit, 50 for
  `partition = "kmeans"` and 100 steps of Lloyd's iteration for
  `partition = "density"`. Under `partition = "density"`, and with
  `init = "kmeans++"`, a kept start that was cut off at the limit is a
  warning of class `"spatialkit_not_converged"`; `partition = "kmeans"`
  from random starts keeps the warning
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) gives.

## Value

An sf object of **at most** `k` cluster-centre POINTs (fewer when `k`
exceeds the number of distinct positions), with `seed_id` and
`method = "kmeans"` columns matching
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
and `attr(, "cover")`, the convex hull of `points_sf` (an sfc in its
CRS). Seeds placed with `partition = "density"` or `init = "kmeans++"`,
and the scored centres of a profile of the density kind or from random
starts, also carry `attr(, "partition")`: a list with `kind`, `init`,
`density_power` and `density_k`, and for seeds this call placed `iter`,
`iter_max` and `converged` (the iterations the kept start took, the
limit, and whether it stopped before the limit).
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
returns the list as `params$partition`. Seeds with no such record are
plain k-means: from random starts when this call placed them, or the
scored centres of a k-means profile from k-means++ starts. When `k` is a
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result that carries `count_rule`, or a profile built with an `estimate`,
the seeds carry that record as `attr(, "count_rule")`, with
`scored_cells` `TRUE` for the scored centres and `FALSE` for centres
this call fits itself at that count (see
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).

## Details

Lon/lat input is projected first so the k-means distances are metric and
not degrees; so is input with no CRS whose coordinates look like lon/lat
(the heuristic
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
applies, with its warning), and the seeds come back in the input's own
coordinates. Rows with empty or non-finite coordinates are dropped with
a warning, and `k` is clamped to the number of distinct positions.

A count given as a number gets a fresh partition, which is not the
partition
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
or
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
scored that count on: those fit their own run, the profile to every row
or to a subsample, and keep the best of their own starts. Pass the
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result or the profile itself as `k` to get the scored one: the centres
of the scored partition are then returned (in the CRS of `points_sf`,
with `attr(, "n_from")`), and no k-means is run.

`partition` says what the seeds minimise. `"kmeans"` is the sum of
squared distances from the points to their nearest seed, from `nstart`
random starts by
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) (Hartigan-Wong)
unless `init` says otherwise; it was the default of earlier versions.
The seeds follow the square root of the point density, so the cells are
compact, smaller where the points are dense, and a cell in a dense place
still holds several times the points of a cell in a sparse place.
`"density"` weights each point by its local density (the number of
points in the disc that reaches its `density_k`-th nearest neighbour,
over the disc's area) to the power `density_power`, and is the default.
At the default power of 1 the seeds follow the point density itself, and
the cells are compact: many small cells where the points are dense, a
few large ones where they are sparse. They hold about the same number of
points each once `k` is above the number of dense places (clusters of
points). With fewer cells than that a cell takes in whole clusters,
large and small, and the counts are about as uneven as under `"kmeans"`.
The density kind is the best of `nstart` starts, each refined by at most
100 steps of Lloyd's iteration (`iter_max` sets another limit, and a
kept start that was cut off at the limit is a warning of class
`"spatialkit_not_converged"`); it does not depend on the order of the
rows or on the unit of the coordinates, and points at one location share
a cell. A seed can end with no point at counts from about a fifth of the
points up, mostly from `init = "random"`; that is a warning of class
`"spatialkit_empty_cells"`. "The partition" in
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
gives the details and a simulated comparison of the two kinds.

The centres lie inside the points, so their own hull is smaller than the
data's. The result carries the convex hull of `points_sf` as
`attr(, "cover")`, and
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
clips the cells to that hull (grown by 2 percent of its diagonal) when
no `boundary` is passed, so every point of `points_sf` falls in a cell.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(100, 0, 1000), y = 5e6 + runif(100, 0, 1000)),
  coords = c("x", "y"), crs = 32632
)
seeds <- voronoi_seeds_kmeans(pts, k = 8)   # density-weighted, 3 starts
nrow(seeds)   # at most 8: one seed per non-empty cluster
#> [1] 8
seeds         # the cluster centres, as an sf POINT layer in the points' CRS
#> Simple feature collection with 8 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 500197.6 ymin: 5000135 xmax: 500802.6 ymax: 5000911
#> Projected CRS: WGS 84 / UTM zone 32N
#>                   geometry seed_id method
#> 1 POINT (500596.5 5000605)       1 kmeans
#> 2 POINT (500371.1 5000409)       2 kmeans
#> 3 POINT (500802.6 5000255)       3 kmeans
#> 4 POINT (500316.8 5000797)       4 kmeans
#> 5 POINT (500197.6 5000577)       5 kmeans
#> 6   POINT (500799 5000911)       6 kmeans
#> 7 POINT (500785.6 5000620)       7 kmeans
#> 8 POINT (500513.3 5000135)       8 kmeans
# The k-means seeds of earlier versions: 10 random starts.
voronoi_seeds_kmeans(pts, k = 8, partition = "kmeans")
#> Simple feature collection with 8 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 500102 ymin: 5000143 xmax: 500813 ymax: 5000897
#> Projected CRS: WGS 84 / UTM zone 32N
#>                   geometry seed_id method
#> 1   POINT (500145 5000616)       1 kmeans
#> 2   POINT (500813 5000237)       2 kmeans
#> 3 POINT (500786.3 5000897)       3 kmeans
#> 4   POINT (500102 5000222)       4 kmeans
#> 5 POINT (500354.4 5000832)       5 kmeans
#> 6 POINT (500377.9 5000474)       6 kmeans
#> 7 POINT (500502.9 5000143)       7 kmeans
#> 8   POINT (500756 5000585)       8 kmeans
```
