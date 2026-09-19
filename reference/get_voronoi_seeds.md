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
  kmeans_nstart = 10,
  kmeans_iter = 100,
  set_seed = NULL
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
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  (read with
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  at its default criterion). The result then carries `attr(, "n_from")`
  saying which.

- seeds:

  sf POINT object of user-provided seeds (method = "provided").

- sample_points:

  Optional sf POINT cloud for k-means clustering. Only the first two
  coordinate columns are clustered, so a Z or M dimension does not join
  the distance calculation and dominate it; rows with empty or
  non-finite coordinates are dropped with a warning, so they never reach
  [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html), which fails
  on them without naming a cause. A lon/lat cloud is projected before
  clustering.

- kmeans_nstart:

  Integer; nstart for kmeans(). Default 10.

- kmeans_iter:

  Integer; iter.max for kmeans(). Default 100.

- set_seed:

  Optional integer RNG seed.

## Value

An sf POINT object with seed_id and method columns. With
`method = "kmeans"` it also carries `attr(, "kmeans")`, the run behind
the seeds: `cluster` (the `seed_id` each clustered cloud point was
assigned to), `rows` (those points' positions in the cloud, since rows
with unusable coordinates are dropped first), `size` (points per seed),
`withinss` and `tot_withinss` (the within-cluster sums of squares),
`iter` and `nstart`.

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
```
