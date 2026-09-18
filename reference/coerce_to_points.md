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
  CRS and its coordinates fall inside the lon/lat envelope, that
  temporary projection interprets them as EPSG:4326 (with a warning) and
  the midpoints returned are geodesic ones brought back to the input's
  numbers, not planar midpoints. Set the CRS, or pass
  `tmp_project = FALSE`, for planar data.

## Value

An sf object with geometry coerced to POINTs.

## Details

LINESTRING midpoints are sampled with
[`sf::st_line_sample()`](https://r-spatial.github.io/sf/reference/st_line_sample.html),
which yields no point for an EMPTY LINESTRING. Rather than silently
misaligning the result (or letting sf crash), such input raises an
error; drop empty geometries first with `x <- x[!sf::st_is_empty(x), ]`.

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
