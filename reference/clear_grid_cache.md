# Clear the in-session grid cache

Empties the cache that
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md)
fills, and returns the number of grids removed.

## Usage

``` r
clear_grid_cache(cache_env = .gmt_cache)
```

## Arguments

- cache_env:

  The cache to clear. The default is the package's own; pass an
  environment only if you passed the same one to
  `create_grid_polygons_cached(cache_env = )`.

## Value

Invisibly, the number of entries removed.

## See also

[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
which fills the cache.

Other options and caches:
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md),
[`spatialkit-conditions`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md),
[`spatialkit-options`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-options.md),
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
invisible(clear_grid_cache())   # start from an empty cache
g1 <- create_grid_polygons_cached(bnd, target_cells = 9)
g2 <- create_grid_polygons_cached(bnd, target_cells = 9)   # cache hit
print(clear_grid_cache())   # 1: the grid cached above
#> [1] 1
print(clear_grid_cache())   # 0: already empty
#> [1] 0
```
