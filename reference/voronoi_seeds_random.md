# Random seed generation within a polygonal boundary

Draws `k` seed points uniformly at random inside `boundary`, ignoring
where the observations are. Reach for this when the cells should cover
the study area evenly (so that sparsely sampled ground still gets its
own cells and is visibly under-sampled in the results) instead of
concentrating resolution where the data already are, which is what
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)
does. It is also the honest choice for a null or sensitivity comparison:
re-running an analysis over several random seedings shows how much of a
result depends on one particular tessellation.

## Usage

``` r
voronoi_seeds_random(boundary, k, set_seed = 456)
```

## Arguments

- boundary:

  An sf or sfc polygonal object.

- k:

  Integer; number of random seeds.

- set_seed:

  Integer RNG seed. Default 456.

## Value

An sf object of **at most** `k` random POINTs (rejection sampling inside
an awkward geometry can fall short of `k`, which is warned about), with
`seed_id` and `method = "random"` columns matching
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md).

## Details

Sampling is by rejection inside the polygon, so an awkward geometry can
return fewer than `k` seeds; that shortfall is warned about rather than
silently padded.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
seeds <- voronoi_seeds_random(bnd, k = 10)
nrow(seeds)   # at most 10: a seed that lands outside the boundary is dropped
#> [1] 10
seeds
#> Simple feature collection with 10 features and 2 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 8.243274 ymin: 21.79086 xmax: 85.21335 ymax: 84.31172
#> Projected CRS: WGS 84 / UTM zone 32N
#>                     geometry seed_id method
#> 1   POINT (8.95516 37.29459)       1 random
#> 2  POINT (21.05123 21.79086)       2 random
#> 3   POINT (73.29553 75.5105)       3 random
#> 4  POINT (85.21335 82.16811)       4 random
#> 5  POINT (78.83979 59.89182)       5 random
#> 6    POINT (33.196 65.10336)       6 random
#> 7  POINT (8.243274 84.31172)       7 random
#> 8  POINT (28.55269 45.32381)       8 random
#> 9  POINT (23.75033 71.67571)       9 random
#> 10 POINT (38.52362 29.12222)      10 random
```
