# Build a tessellation (Voronoi, Delaunay triangles, hex grid, or square grid)

The single entry point for turning a point pattern into analysis
regions, and the first step of the package's pipeline. It wraps the four
tessellation methods behind one interface that handles CRS projection,
clipping and stable cell identifiers consistently, and returns the cell
layer together with the point-to-cell index that
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
and
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
consume. Prefer it to the individual constructors whenever you might
want to compare methods: the return shape does not change with `method`,
so swapping `"voronoi"` for `"hex"` costs one argument.

## Usage

``` r
build_tessellation(
  points_sf,
  boundary = NULL,
  method = c("voronoi", "triangles", "hex", "square"),
  approx_n_cells = NULL,
  cellsize = NULL,
  expand = 0,
  clip = TRUE,
  keep_duplicates = FALSE,
  crs = NULL,
  quiet = FALSE
)
```

## Arguments

- points_sf:

  An sf object with POINT/MULTIPOINT geometry.

- boundary:

  Polygonal sf/sfc study area. **Required** for `method = "hex"` and
  `method = "square"`, which have no extent of their own to lay a grid
  over and error without it; supply the study-area polygon, or build one
  from the points with
  [`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md).
  **Optional** for `method = "voronoi"` and `method = "triangles"`,
  which derive their extent from the points themselves and use
  `boundary` only to clip the result when `clip = TRUE`.

- method:

  One of "voronoi", "triangles", "hex", "square".

- approx_n_cells:

  Approximate number of cells. Read by `method = "hex"` and `"square"`
  only: `"voronoi"` grows one cell per input point and `"triangles"` one
  triangle per neighbouring triple, so neither has a count to set, and
  both warn that the argument was ignored. For a Voronoi cell count,
  place the seeds with
  [`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md)
  and tessellate those. For hex grids the target is adjusted for packing
  density; the actual count after clipping to an irregular boundary may
  differ noticeably. Besides a number, this accepts what the
  level-selection step returned: the integer vector of ranked candidates
  from
  [`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
  (its first element is used), a
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  result (its `$best`), or a
  [`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md)
  (read with
  [`select_resolution()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_resolution.md)
  at its default criterion). The count used is returned as
  `params$approx_n_cells` and where it came from as
  `params$approx_n_cells_from` (`NULL` for a plain number).

- cellsize:

  Numeric cell size, in the units of the working CRS. Read by
  `method = "hex"` and `"square"` only; the other two methods warn that
  it was ignored. When both `cellsize` and `approx_n_cells` are given,
  `cellsize` wins and `approx_n_cells` is ignored with a logged warning;
  supply one or the other.

- expand:

  Buffer distance for the Voronoi envelope. Applied by
  `method = "voronoi"` only; the `"hex"`, `"square"` and `"triangles"`
  methods ignore it (the value you passed is still echoed back in
  `params$expand`).

- clip:

  Logical; clip to boundary.

- keep_duplicates:

  Logical; keep duplicate points.

- crs:

  Optional target CRS.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

## Value

A list with components:

- `cells`:

  An sf polygon layer, one row per cell. It always carries a `cell_id`
  column; the `"hex"` and `"square"` methods additionally carry
  `poly_id`, which holds the same values.

- `index`:

  Integer vector of `cell_id` values, one per row of `points_sf`, and
  `NA` for a point that falls inside no cell, that is, one outside the
  study area. Only a point within a thousandth of the median cell width
  of a cell is snapped to it. That covers points sitting exactly on a
  shared edge, and leaves points outside the study area as `NA`. A
  summary built from `index` therefore counts only the points the
  tessellation actually covers.

- `boundary`:

  The boundary used (possibly derived and/or reprojected).

- `method`:

  The method actually used.

- `params`:

  The parameters the tessellation was built with, plus `snapped`, the
  record of that nearest-cell repair: a list with `n`, `which` (row
  positions in `points_sf`) and `distance` (how far outside every cell
  each sat, in CRS units).

## Details

Which method to reach for. `"voronoi"` gives one cell per point, so
resolution follows sampling density. That is the choice when the
observations themselves define the regions. `"hex"` and `"square"` give
equal-area cells on a fixed grid, so cell size is a decision you make,
and the data does not make it for you; hexagons avoid the axis-aligned
artefacts of squares and have uniform neighbour distances. `"triangles"`
returns the Delaunay triangulation, useful for interpolation and
adjacency work; it is not meant as an aggregation unit.
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md)
will suggest a cell count from the spatial structure of the data.

## See also

Other tessellation:
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
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
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(20, 0, 100), y = 5e6 + runif(20, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
tess <- build_tessellation(pts, method = "voronoi", quiet = TRUE)
tess$cells
#> Simple feature collection with 20 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 500003.6 ymin: 4999999 xmax: 500101.8 ymax: 5000096
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>                          geometry cell_id
#> 1  POLYGON ((500016.9 5000038,...       1
#> 2  POLYGON ((500021.3 5000077,...       2
#> 3  POLYGON ((500016.9 5000038,...       3
#> 4  POLYGON ((500007.8 5000047,...       4
#> 5  POLYGON ((500044.6 5000090,...       5
#> 6  POLYGON ((500021.3 5000077,...       6
#> 7  POLYGON ((500043.4 5000044,...       7
#> 8  POLYGON ((500043.6 5000044,...       8
#> 9  POLYGON ((500053.2 5000027,...       9
#> 10 POLYGON ((500070.6 5000087,...      10
```
