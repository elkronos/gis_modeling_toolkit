# Clear the in-session grid cache

Removes all memoized grid results from the internal cache environment.

## Usage

``` r
clear_grid_cache(cache_env = .gmt_cache)
```

## Arguments

- cache_env:

  Environment to clear. Default .gmt_cache.

## Value

Invisibly, the number of entries removed.

## See also

Other package options and caches:
[`clear_fitted_cache()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_fitted_cache.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
g1 <- create_grid_polygons_cached(bnd, target_cells = 9)
g2 <- create_grid_polygons_cached(bnd, target_cells = 9)   # cache hit
clear_grid_cache()                                          # entries removed
```
