# Create Voronoi polygons from points with CRS handling and optional clipping

Assigns every location in the study area to its nearest input point,
giving one cell per point. This is the tessellation to reach for when
the observations themselves define the regions of interest (sampling
sites, monitoring stations, service points), because cell size then
adapts to sampling density instead of being imposed by a fixed grid:
dense areas get small cells and sparse areas large ones. Prefer
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
instead when you need equal-area cells or a resolution independent of
where the data happen to be.

## Usage

``` r
create_voronoi_polygons(
  points_sf,
  boundary = NULL,
  expand = 0,
  clip = TRUE,
  keep_duplicates = FALSE,
  crs = NULL,
  quiet = FALSE
)
```

## Arguments

- points_sf:

  An sf object with POINT/MULTIPOINT geometries.

- boundary:

  Optional polygonal sf object.

- expand:

  Numeric; absolute buffer distance for the envelope.

- clip:

  Logical; intersect cells with boundary.

- keep_duplicates:

  Logical; keep coincident points for graph construction.

- crs:

  Optional target CRS.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

## Value

A list with `cells`, `index`, `boundary`, `method` and `params`. `index`
holds one `cell_id` per row of `points_sf`, and `NA` for a point that
falls outside every cell, which means outside the study area, so a
summary built from it counts only the points the tessellation actually
covers.

## Details

The heavy lifting is
[`sf::st_voronoi()`](https://r-spatial.github.io/sf/reference/geos_unary.html).
What this adds is the surrounding bookkeeping: projecting lon/lat input,
building and buffering an envelope so edge cells are bounded, clipping
to `boundary`, restoring the point-to-cell correspondence that
[`st_voronoi()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
scrambles, and stamping stable `cell_id` values.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`create_grid_polygons_cached()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons_cached.md),
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
  data.frame(x = 5e5 + runif(15, 0, 100), y = 5e6 + runif(15, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
res <- create_voronoi_polygons(pts, quiet = TRUE)
res$cells   # one polygon per unique point, with stable cell_id
#> Simple feature collection with 15 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 500003.5 ymin: 4999999 xmax: 500097.1 ymax: 5000102
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>                          geometry cell_id
#> 1  POLYGON ((500024.2 5000020,...       1
#> 2  POLYGON ((500030.5 5000080,...       2
#> 3  POLYGON ((500007.1 5000040,...       3
#> 4  POLYGON ((500024.2 5000020,...       4
#> 5  POLYGON ((500048.1 5000053,...       5
#> 6  POLYGON ((500025 5000064, 5...       6
#> 7  POLYGON ((500055.9 5000081,...       7
#> 8  POLYGON ((500037.4 5000019,...       8
#> 9  POLYGON ((500071.5 5000085,...       9
#> 10 POLYGON ((500044.6 5000033,...      10
res$index   # cell_id assignment for each input point
#>  [1]  5  6  9 15  2 13 14 11  8  1  3  4 10  7 12
```
