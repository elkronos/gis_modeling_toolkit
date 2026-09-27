# Coerce arbitrary geometries to representative points

Converts the geometry column of an sf object to POINTs using one of
several strategies.

## Usage

``` r
coerce_to_points(
  x,
  mode = c("auto", "centroid", "point_on_surface", "surface", "line_midpoint",
    "bbox_center"),
  tmp_project = TRUE
)
```

## Arguments

- x:

  An sf object.

- mode:

  One of "auto", "centroid", "point_on_surface", "surface",
  "line_midpoint", "bbox_center".

- tmp_project:

  Logical; temporarily project for line-based midpoints. When `x` has no
  CRS and the lon/lat heuristic of
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  takes its coordinates for degrees (inside the lon/lat envelope and
  more than one unit across, or with decimal-degree precision), that
  temporary projection interprets them as EPSG:4326 (with a warning) and
  the midpoints returned are geodesic ones brought back to the input's
  numbers, not planar midpoints. Set the CRS, or pass
  `tmp_project = FALSE`, for planar data.

## Value

An sf object with geometry coerced to POINTs, row for row with `x`; an
empty input geometry gives an empty POINT.

## Details

The result has one row per row of `x`, in the same order. An EMPTY
geometry of any type, lines included, becomes an EMPTY POINT in its own
row;
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
and
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
then drop such rows, as they drop any other empty geometry. Empty lines
are never handed to
[`sf::st_line_sample()`](https://r-spatial.github.io/sf/reference/st_line_sample.html):
it yields no midpoint for them, which would misalign the result, and
with sf 1.0.x an empty MULTILINESTRING (or an empty part of one) crashed
the R session. An empty part inside a non-empty feature is ignored, so
the feature gets the point its other parts give; GEOS's interior point,
used by `"point_on_surface"` and by the temporary projection's choice of
CRS, segfaulted on an empty line part too.

## See also

Other spatial data preparation:
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

## Examples

``` r
library(sf)
poly <- st_sf(
  id = 1,
  geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)
  ))), crs = 32632)
)
coerce_to_points(poly, "auto")  # interior representative point
#> Simple feature collection with 1 feature and 1 field
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1 ymin: 1 xmax: 1 ymax: 1
#> Projected CRS: WGS 84 / UTM zone 32N
#>   id    geometry
#> 1  1 POINT (1 1)
```
