# Reduce each geometry to one representative point

Replaces the geometry of each row of an sf object with a point, by the
rule `mode` names. The result has one row per row of `x`, in the same
order.

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

  How each geometry becomes a point: `"auto"` (default), `"centroid"`,
  `"point_on_surface"`, `"surface"` (an alias for `"point_on_surface"`),
  `"line_midpoint"` or `"bbox_center"`. See "Modes".

- tmp_project:

  Logical; take line midpoints in a temporary projection. Default
  `TRUE`. With `FALSE`, lines in longitude and latitude get their
  centroid. When `x` has no CRS and the lon/lat test of
  [`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  takes its coordinates for degrees, the temporary projection reads them
  as EPSG:4326 (with a warning), and the midpoints returned are geodesic
  ones brought back to the input's numbers, not planar midpoints. Set
  the CRS, or pass `tmp_project = FALSE`, for planar data.

## Value

An sf object with its geometry reduced to POINTs, row for row with `x`;
an empty input geometry gives an empty POINT.

## Modes

- `"auto"` (default):

  POINT rows are kept. POLYGON and MULTIPOLYGON rows get a point on the
  surface. A LINESTRING gets its midpoint, and a MULTILINESTRING the
  midpoint of its longest part (see "Lines"). MULTIPOINT rows and rows
  of other types, such as geometry collections, get their centroid.

- `"centroid"`:

  The centroid of every row.

- `"point_on_surface"` (alias `"surface"`):

  A point that lies on each geometry, taken with GEOS.

- `"line_midpoint"`:

  The midpoint of each LINESTRING, taken in a temporary projection (with
  `tmp_project = FALSE` on lon/lat data, its centroid). Rows of other
  types (points, polygons and their MULTI forms) get their centroid. A
  MULTILINESTRING or GEOMETRYCOLLECTION row is an error: cast it to
  LINESTRING first, or use `"auto"`, which takes the midpoint of a
  MULTILINESTRING's longest part.

- `"bbox_center"`:

  The centre of each row's bounding box.

## Empty geometries and rows without a location

An empty geometry of any type, lines included, becomes an empty POINT in
its own row;
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
and
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)
then drop such rows, as they drop any other empty geometry. An empty
part inside a non-empty feature is ignored, so the feature gets the
point its other parts give.

A geometry with no location (see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
"Rows without a location") gets an empty POINT in every mode, as an
empty geometry does. So does a geometry of a layer in a geographic CRS
with a vertex at a latitude beyond 90 degrees, and a warning of class
`spatialkit_lat_beyond_90` gives its rows of `x` (`rows`).

## Lines

Modes `"auto"` and `"line_midpoint"` take the midpoint of a line in a
temporary projection (`tmp_project`), the one
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
chooses for the lines. A line that projection returns empty or without
some of its vertices is not sampled: the call stops with an error of
class `spatialkit_features_degenerate` (see "Errors").

In mode `"auto"` a MULTILINESTRING gets the midpoint of its longest
part. When several parts share the greatest length (to one part in ten
million, measured in the projection used), their midpoints are taken and
the one lowest in x, then in y, is kept, so the point does not depend on
the order in which the feature lists its parts. x values within one part
in ten million of the diagonal of the midpoints' bounding box count as
one, so y decides between midpoints whose x differ only in their last
bits, and the point does not move with the origin of the coordinates.

## Longitude and latitude

Modes `"auto"` (for polygons) and `"point_on_surface"` take the point
with GEOS, which reads a layer in longitude and latitude on the plane in
degrees, under either
[`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
setting. sf's warning that such a point "may not give correct results
for longitude/latitude data" is not passed on. A polygon that spans more
than 180 degrees of longitude (one written across longitude 180, from
179.5 to -179.5) GEOS reads the long way round the globe, so its point
is tested on the polygon read on the sphere with s2; when it does not
lie on it, the point is taken again with the longitudes below 0 read
from 180 to 360, and brought back into -180 to 180. A polygon whose
point lies on it neither way, or that s2 cannot read, keeps GEOS's
point, and a warning of class `spatialkit_surface_point_off` gives its
rows of `x` (`rows`). When `x` has no CRS and the lon/lat test of
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
takes its coordinates for degrees, its polygons that span more than 180
units are read so too, as EPSG:4326, and a warning names the rows whose
point that reading moves; the points keep the missing CRS.

Mode `"bbox_center"`, and the centroids that GEOS takes on the plane in
degrees (a lon/lat layer with s2 off) under `"centroid"`, for the rows
other than LINESTRINGs under `"line_midpoint"`, and for MULTIPOINTs and
collections under `"auto"`, read a geometry written across longitude 180
the same way: one that spans more than 180 degrees of longitude as
written, and no more than 180 with the longitudes below 0 read from 180
to 360, gets the centre of that box, or GEOS's centroid of the geometry
so read, brought back into -180 to 180. A layer without a CRS that the
lon/lat test takes for degrees is read so too, with the same warning. A
geometry wider than 180 degrees in both readings, such as one around a
pole, keeps the centre of its box on the plane, or GEOS's planar
centroid. The centroids s2 takes on the sphere are not changed.

With s2 on, sf takes the centroids of a lon/lat layer with s2 (under
`"centroid"`, for the rows other than LINESTRINGs under
`"line_midpoint"`, and for collections under `"auto"`), and s2 refuses
some polygons GEOS accepts. A polygon s2 refuses because a vertex
repeats the one before it, on its own or inside a collection, gets the
centroid of the same ring without the repeat; the ring's points are
unchanged. A polygon s2 refuses for another reason (edges that cross, a
figure eight) is an error (see "Errors"). With s2 off, the centroids are
GEOS's, planar in degrees.

## Z and M coordinates

A layer whose geometries do not all have the same dimensions (XY beside
XYZ or XYM, as [`rbind()`](https://rdrr.io/r/base/cbind.html) of an XY
and an XYZ layer gives) is read in x and y, and its points come back in
XY with a warning of class `spatialkit_zm_dropped` whose `rows` are the
rows that carried Z or M (see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
"Mixed dimensions, Z and M"). A layer in one dimension, XY or XYZ
throughout, is read as it is, and so is a layer whose only geometries of
another dimension than the first are empty and not points.

GEOS does not read M coordinates. When the mode reduces a geometry that
carries M (XYM or XYZM) to a point (any such geometry under
`"centroid"`, `"point_on_surface"` and `"line_midpoint"`, any that is
not a POINT under `"auto"`), the layer is reduced in x and y: Z and M
are dropped from every geometry, the points come back in XY, and a
warning of class `spatialkit_zm_dropped` says so, with the rows whose
geometry carries M as its `rows`. Empty geometries are not counted. A
POINT layer that carries M is returned as it is under `"auto"`. Mode
`"bbox_center"` reads the box of each geometry and gives XY points for
any layer, with no warning about Z or M.

## Errors

Each of these is an error of class `spatialkit_features_degenerate` that
gives the rows of `x` (`rows` on the condition):

- A polygon s2 refuses:

  In longitude and latitude with s2 on, a polygon whose centroid s2
  cannot take even without its repeated vertices. The message quotes s2
  and says to transform `x` to a projected CRS or to use another `mode`.

- A centroid outside its feature:

  GEOS, which takes the centroids in a projected CRS, without a CRS, and
  in degrees with s2 off, weights each part of a ring by its signed
  area. Where the edges of a ring cross (a bowtie), or a hole lies
  outside its shell, parts of the feature count with opposite signs, and
  the point can lie far outside the feature. A centroid outside the
  bounding box of a feature GEOS reports invalid is refused: the
  condition carries `site = "geos_centroid"` and GEOS's reason for each
  row (`reasons`), and the message says to pass
  `mode = "point_on_surface"` or to repair the features with
  [`sf::st_make_valid()`](https://r-spatial.github.io/sf/reference/valid.html).
  A centroid inside the box is returned as GEOS gives it.

- A line the temporary projection cannot place:

  A line that the projection returns empty or without some of its
  vertices. The condition carries `site = "projection"`, and the message
  says to pass `mode = "centroid"` or `tmp_project = FALSE`, which take
  the point in longitude and latitude, without the projection.

The functions that reduce a layer they are given to points this way
([`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md),
[`fold_separation()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fold_separation.md),
[`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
[`resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/resolution_profile.md),
[`determine_optimal_levels()`](https://elkronos.github.io/gis_modeling_toolkit/reference/determine_optimal_levels.md),
[`cell_values_at()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cell_values_at.md),
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md),
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md),
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md),
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md))
raise each of these errors in their own name, with the same class and
fields, naming their argument and its rows.
[spatialkit-conditions](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-conditions.md)
lists every class.

## See also

[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
for how input geometry and CRSs are read.

Other data preparation:
[`[.spatialkit_rows()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
[`spatialkit-geometry`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)

## Examples

``` r
library(sf)
poly <- st_sf(
  id = 1,
  geometry = st_sfc(st_polygon(list(rbind(
    c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)
  ))), crs = 32632)
)
coerce_to_points(poly, "auto")  # a point on the polygon's surface
#> Simple feature collection with 1 feature and 1 field
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 1 ymin: 1 xmax: 1 ymax: 1
#> Projected CRS: WGS 84 / UTM zone 32N
#>   id    geometry
#> 1  1 POINT (1 1)
```
