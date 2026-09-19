# Build a polygonal clip target from points and/or a boundary

Resolves the single polygon that every tessellation method clips
against. With a `boundary` it is that boundary (optionally buffered by
`expand`); without one it is the convex hull of `points_sf`, again
optionally buffered. Reach for it when you want to see or reuse the
exact clip target
[`build_tessellation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/build_tessellation.md)
will apply, for instance to check that a study-area polygon actually
contains the observations before tessellating, or to pass the same
envelope to
[`create_voronoi_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_voronoi_polygons.md)
and
[`create_grid_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/create_grid_polygons.md)
so that two tessellations of one dataset cover identical ground.

## Usage

``` r
clip_target_for(points_sf, boundary = NULL, expand = 0, quiet = FALSE)
```

## Arguments

- points_sf:

  An sf object with POINT/MULTIPOINT geometry.

- boundary:

  Optional polygonal sf object.

- expand:

  Numeric expansion distance or fraction (0–1 = fraction of extent).
  Absolute values are expressed in the units of the CRS the clip target
  is built in. Because
  [`sf::st_buffer()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
  interprets `dist` as **metres** for lon/lat data while the
  fraction-of-extent form is derived from a bounding box measured in
  **degrees**, lon/lat input is projected to a local projected CRS first
  (see `@return`), so both forms agree.

- quiet:

  Logical; suppress this function's progress
  [`message()`](https://rdrr.io/r/base/message.html)s. It does not
  silence R warnings, nor the package's console log echo (see
  [`spatialkit_quiet`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
  for that). Default `FALSE`.

## Value

An sf polygon layer representing the clip target. For lon/lat input the
layer is returned in the automatically selected local projected CRS, not
the input CRS; a message reports this unless `quiet = TRUE`.

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(30, 0, 100), y = 5e6 + runif(30, 0, 100)),
  coords = c("x", "y"), crs = 32632
)
# No boundary: the convex hull, expanded by 10% of the extent
hull <- clip_target_for(pts, expand = 0.1, quiet = TRUE)
st_area(hull)
#> 12054.18 [m^2]
```
