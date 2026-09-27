# K-means seed generation from point coordinates

Places `k` seed points at k-means cluster centres of the observed
coordinates, so seeds (and the Voronoi cells built from them) follow the
sampling density: clusters of observations attract seeds, empty ground
gets none. Reach for this when you want cells that follow the data, so
that counts per cell vary far less than on a fixed grid over clustered
points (k-means does not equalise them, it minimises the spread of
points around each centre), which is what keeps per-cell aggregates in
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
from resting on one or two observations. Use
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)
instead when you want coverage of the study area rather than of the
data, and
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
to pick between them by name.

## Usage

``` r
voronoi_seeds_kmeans(points_sf, k, set_seed = 456, nstart = 10)
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
  them, and otherwise its count is used.

- set_seed:

  Optional integer RNG seed. Default 456, so a call gives the same seeds
  every time whatever the session's random-number state; an outer
  [`set.seed()`](https://rdrr.io/r/base/Random.html) does not change
  them, and the caller's random-number stream is left as it was. Pass
  `NULL` to draw the k-means starts from the session's stream instead
  (the default of
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)).

- nstart:

  Number of random starts for
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html); the best is
  kept. Default 10.

## Value

An sf object of **at most** `k` cluster-centre POINTs (fewer when `k`
exceeds the number of distinct positions), with `seed_id` and
`method = "kmeans"` columns matching
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md).

## Details

Lon/lat input is projected first so the k-means distances are metric and
not degrees; so is input with no CRS whose coordinates look like lon/lat
(the heuristic
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
applies, with its warning), and the seeds come back in the input's own
coordinates. Rows with empty or non-finite coordinates are dropped with
a warning, and `k` is clamped to the number of distinct positions.

The partition is
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) (Hartigan-Wong)
with `nstart` random starts.
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
and
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
score each count on a different run, by default the best of 25 k-means++
restarts, which usually reaches a lower within-cluster sum of squares;
the seeds for a count passed as a number are therefore not the partition
that count was scored on. Raising `nstart` narrows the gap but does not
close it. Pass the
[`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
result or the profile itself as `k` to close it: the centres of the
scored partition are then returned (in the CRS of `points_sf`, with
`attr(, "n_from")`), and no k-means is run.

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
seeds <- voronoi_seeds_kmeans(pts, k = 8)
nrow(seeds)   # at most 8: one seed per non-empty cluster
#> [1] 8
seeds         # the cluster centres, as an sf POINT layer in the points' CRS
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
