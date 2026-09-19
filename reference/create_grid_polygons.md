# Create square or hexagonal grid polygons over a boundary

Lays a regular grid of equal-area cells over `boundary` and clips it to
that boundary. Reach for this rather than
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
when cell size should be a decision you make, instead of one dictated by
where the observations happen to be. That is the case when you need
per-cell rates comparable across the map, or a resolution that stays
fixed as the sample grows. Hexagons (`type = "hex"`) avoid the
axis-aligned artefacts of squares and give every cell the same distance
to all six neighbours, which matters for anything that reads
neighbourhoods.

## Usage

``` r
create_grid_polygons(
  boundary,
  target_cells = NULL,
  type = c("square", "hex"),
  cellsize = NULL,
  n = NULL,
  clip = TRUE,
  crs = NULL,
  quiet = FALSE,
  max_cells = 1e+06
)
```

## Arguments

- boundary:

  Polygonal sf or sfc object.

- target_cells:

  Optional approximate desired number of cells. The cell *size* is
  derived from it as `sqrt(area / target_cells)`, so square grids get
  square cells; for hex grids the count is adjusted for hexagonal
  packing density. The word "approximate" is load bearing: a grid of
  square cells over an elongated bounding box needs more of them than a
  grid of rectangles would (a 1000 x 1 strip at `target_cells = 9`
  yields cells of side 10.5 and about 95 of them), and clipping to an
  irregular boundary moves the count again. Pass `cellsize` when the
  count matters more than the shape.

- type:

  Grid type: `"square"` (the default) or `"hex"`.

- cellsize:

  Optional numeric cell size (length 1 or 2), in the units of the
  working CRS. Takes precedence over `n`: if both are supplied,
  `cellsize` is used, `n` is ignored and a warning is logged. Supply
  exactly one of `target_cells`, `cellsize` and `n`. For `type = "hex"`
  a hexagon is defined by a single edge-to-edge distance, so only
  `cellsize[1]` is used and a differing `cellsize[2]` is ignored with a
  logged warning.

- n:

  Optional grid resolution (integer, length 1 or 2) giving the number of
  columns and rows to divide the boundary's bounding box into; the cell
  size is derived from it.
  [`sf::st_make_grid()`](https://r-spatial.github.io/sf/reference/st_make_grid.html)
  derives hexagon placement from `cellsize` alone, so for `type = "hex"`
  `n` does not set the number of cells, although it does change the
  grid, because the `cellsize` derived from it is what the hexagons are
  built with. Ignored (with a logged warning) when `cellsize` is also
  supplied. Passing both would otherwise truncate the grid to `n[1]` x
  `n[2]` cells anchored at the bounding-box corner, covering only part
  of the boundary.

- clip:

  Logical; clip grid to boundary.

- crs:

  Optional target CRS. When `NULL` (default) the boundary is projected
  with
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
  which changes the CRS of the returned grid; a message reports this
  unless `quiet = TRUE`.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

- max_cells:

  Upper bound on the number of cells the grid may have, estimated from
  the boundary's bounding box before anything is built. Default `1e6`. A
  `cellsize` in the wrong units (metres on a boundary in kilometres,
  say) asks for a grid that cannot be built, and this refuses it with a
  message before memory is exhausted. Set to `Inf` to disable.

## Value

An sf polygon layer with poly_id column.

## Details

Size the grid with exactly one of `target_cells` (roughly how many cells
you want, the package derives the rest), `cellsize` (a fixed edge length
in CRS units) or `n` (a fixed number of columns and rows). See
`@param cellsize` for what happens when more than one is given.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md),
[`ensure_stable_poly_id()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_stable_poly_id.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md),
[`voronoi_seeds_kmeans()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_kmeans.md),
[`voronoi_seeds_random()`](https://elkronos.github.io/gis_modeling_toolkit/reference/voronoi_seeds_random.md)

## Examples

``` r
library(sf)
bnd <- st_sf(geometry = st_sfc(st_polygon(list(rbind(
  c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)
))), crs = 32632))
grid_sq  <- create_grid_polygons(bnd, target_cells = 100, type = "square")
grid_hex <- create_grid_polygons(bnd, target_cells = 100, type = "hex")
nrow(grid_sq)
#> [1] 100
# Hex counts run above target because clipping keeps every hexagon that
# merely overhangs the boundary; the inflation is proportionally larger
# at small target_cells.
nrow(grid_hex)
#> [1] 114
```
