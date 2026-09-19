# Create and cache grid polygons over a boundary

Builds a grid via
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
and memoizes the result so repeated calls with the same inputs return
instantly.

## Usage

``` r
create_grid_polygons_cached(
  boundary,
  target_cells,
  type = c("square", "hex"),
  ...,
  cache_env = .gmt_cache,
  max_entries = 50L
)
```

## Arguments

- boundary:

  An sf or sfc polygonal object.

- target_cells:

  Approximate desired number of cells.

- type:

  Grid type: `"square"` (the default) or `"hex"`, matching
  [`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md).

- ...:

  Additional arguments forwarded to create_grid_polygons().

- cache_env:

  Environment for memoized grids. Default .gmt_cache.

- max_entries:

  Maximum number of grids the cache holds. Default 50. Once full, adding
  a grid evicts the one added earliest, so a loop over many boundaries
  holds at most this many grids (about 2 MB per 2,500-cell grid) rather
  than every grid it ever built for the life of the session.
  [`clear_grid_cache`](https://elkronos.github.io/gis_modeling_toolkit/reference/clear_grid_cache.md)
  empties it outright.

## Value

An sf data frame with a stable poly_id column. The rows are re-ordered
and re-numbered by
[`ensure_stable_poly_id`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
which
[`create_grid_polygons`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
does not do: the same cell therefore carries a different `poly_id`
depending on which of the two builders produced it. Use one builder
throughout an analysis; joining a summary keyed on IDs from one onto
geometries from the other draws the values on the wrong polygons.

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
g <- create_grid_polygons_cached(bnd, target_cells = 16, type = "hex")
nrow(g)
#> [1] 27
head(g$poly_id)   # stable IDs from ensure_stable_poly_id()
#> [1] 1 2 3 4 5 6
```
