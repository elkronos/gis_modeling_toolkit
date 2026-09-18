# Create Voronoi polygons from points with robust CRS and optional clipping

Assigns every location in the study area to its nearest input point,
giving one cell per point. This is the tessellation to reach for when
the observations themselves define the regions of interest — sampling
sites, monitoring stations, service points — because cell size then
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
falls outside every cell – outside the study area, in other words – so a
summary built from it counts only the points the tessellation actually
covers.

## Details

The heavy lifting is
[`sf::st_voronoi()`](https://r-spatial.github.io/sf/reference/geos_unary.html);
what this adds is the surrounding bookkeeping — projecting lon/lat
input, building and buffering an envelope so edge cells are bounded,
clipping to `boundary`, restoring the point-to-cell correspondence that
[`st_voronoi()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
scrambles, and stamping stable `cell_id` values.

## See also

Other tessellation:
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md),
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md),
[`get_voronoi_seeds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/get_voronoi_seeds.md),
[`plot_tessellation_map()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_tessellation_map.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = runif(15, 0, 100), y = runif(15, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
res <- create_voronoi_polygons(pts, quiet = TRUE)
res$cells   # one polygon per unique point, with stable cell_id
#> Simple feature collection with 15 features and 1 field
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 3.542734 ymin: -1.296593 xmax: 97.10325 ymax: 101.8263
#> Projected CRS: WGS 84 / UTM zone 32N
#> First 10 features:
#>                          geometry cell_id
#> 1  POLYGON ((24.18211 19.5762,...       1
#> 2  POLYGON ((30.51111 79.94056...       2
#> 3  POLYGON ((7.149734 40.23226...       3
#> 4  POLYGON ((24.18211 19.5762,...       4
#> 5  POLYGON ((48.09213 52.90715...       5
#> 6  POLYGON ((24.95868 64.12208...       6
#> 7  POLYGON ((55.8777 80.67853,...       7
#> 8  POLYGON ((37.41231 18.53199...       8
#> 9  POLYGON ((71.52246 84.72244...       9
#> 10 POLYGON ((44.60761 32.96868...      10
res$index   # cell_id assignment for each input point
#>  [1]  5  6  9 15  2 13 14 11  8  1  3  4 10  7 12
```
