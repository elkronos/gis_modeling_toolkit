# Project an sf object to a CRS suited to its extent

Returns a spatial object in a projected CRS, so that distances and areas
computed on it are close to true. Without `target_crs`, an object that
already has a projected CRS is returned as it is, and lon/lat input is
moved to a CRS chosen from its extent; with `target_crs`, the object is
moved to that CRS.

## Usage

``` r
ensure_projected(x, target_crs = NULL, purpose = c("distance", "area"))
```

## Arguments

- x:

  An sf or sfc object (other objects are returned unchanged).

- target_crs:

  Optional target CRS (sf object, integer EPSG, or crs). Must resolve to
  a usable CRS via
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html);
  an unusable value (one that resolves to `NA_crs_`) raises an error, so
  `x` is never left unprojected without a word.

- purpose:

  `"distance"` (default): the candidate that distorts pairwise distances
  least, which is what ranges, block sizes, bandwidths and length-scales
  read off the coordinates. `"area"`: an equal-area projection, for
  densities or rates per cell. For lon/lat input the choice is then made
  among a Lambert azimuthal centred on the data and an Albers conic
  (where its parallels do not degenerate or pass a pole), whichever
  distorts distances less; a UTM zone, which is not equal-area, never
  enters that comparison, and global coverage gets Equal Earth in place
  of Web Mercator. For input already projected, the area distortion is
  measured and a warning is logged above 1 percent; see "Equal-area
  check". Ignored when `target_crs` is given.

## Value

`x` in `target_crs` when one is given. Otherwise `x` in the projected
CRS chosen here, as it came when it already had a projected CRS, or
without a CRS when it had none and the lon/lat test declined it (see
"Missing CRS"). Two attributes record how the CRS was reached:

- `crs_choice`:

  Present when a projection was chosen here (lon/lat input, no
  `target_crs`): a data.frame with one row per projection considered,
  holding `name`, `crs` (its definition), the measured worst-case
  `distance_error` (relative, over sampled pairs; `NA` where it could
  not be measured) and `chosen`. The figure is measured for the single
  candidate on the paths where no comparison runs (a UTM zone on a local
  extent, the equal-area projection chosen for a layer straddling the
  antimeridian). `NULL` when no local projection was chosen here: input
  that already had a projected CRS, a `target_crs` you supplied, or an
  extent no local projection fits, which falls back to Web Mercator or
  Equal Earth.

- `crs_assumed`:

  On input without a CRS: `"EPSG:4326"` when the lon/lat test took the
  coordinates for degrees, `"none"` when it declined. The attribute is
  also read on the way in: an object already carrying `"none"` is
  returned as it is, without the test, which is how a
  [`predict()`](https://rdrr.io/r/stats/predict.html) method replays a
  fit's decision so that a subset of the training rows is not judged
  differently from the whole.

## Details

Without `target_crs`, only geographic (lon/lat) input is transformed,
and the CRS chosen depends on the extent of the data:

- Local extents:

  The UTM zone containing the data's centre (EPSG:326xx north of the
  equator, EPSG:327xx south). Distances and areas are close to true over
  a few degrees of longitude, which is the case this package is usually
  in. The centre is the centroid on the sphere, computed with s2
  whatever
  [`sf::sf_use_s2()`](https://r-spatial.github.io/sf/reference/s2.html)
  is set to, so data near a zone edge get the same zone in every
  session.

- Wide extents:

  Once the data reach well beyond the roughly 3 degrees a UTM zone is
  designed for, a single zone can distort distances by several percent,
  and that error propagates into variogram ranges, block sizes, GWR
  bandwidths and GP length-scales. The projection used is then the one
  with the smallest measured distortion among the zone, a Lambert
  azimuthal equal-area centred on the data and (where its standard
  parallels do not degenerate or pass a pole) an Albers conic. Each is
  scored by projecting representative points of the data and comparing
  planar with geodesic distances over the pairs at least 1 mm apart. A
  non-POINT layer is reduced to one point per feature, plus the vertices
  of its outline when fewer than 40 distinct points remain, so a single
  study-area polygon is scored too, however many rows repeat it; copies
  of one location count once, and rows without a finite location are
  left out. A candidate must place every location of the layer to be
  used: one under which a vertex of the layer has no coordinates (a UTM
  zone about 90 degrees or more from its central meridian) gets
  `distance_error` `NA` and is not chosen, and the log line names it. A
  vertex at a latitude beyond 90, which no CRS places, takes no part in
  this test. When no candidate places every vertex, Web Mercator (Equal
  Earth for `purpose = "area"`) is used. The choice, its distance error
  and the zone's, and this argument are logged (see
  [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md));
  they are not R warnings, so `tryCatch(warning = )` does not see them.

- Antimeridian:

  Data straddling longitude 180 have a bounding box wider than a
  hemisphere. The wrap is detected from the coordinates (one very large
  gap in the sorted longitudes) and an equal-area projection centred on
  the true extent is used.

- Around a pole:

  Data spanning more than 180 degrees of longitude with no such gap
  surround a pole. When every point also lies on one side of the equator
  (Antarctic stations, a pan-Arctic network), a Lambert azimuthal
  equal-area centred on that pole is used, provided it measures a
  smaller distance error than the global fallback. Web Mercator splits
  such a layer at longitude 180 and stretches it towards the pole. Only
  the coverage left over, spanning both hemispheres or a low-latitude
  belt the polar projection fits worse, falls back to EPSG:3857 (Equal
  Earth for `purpose = "area"`).

- Missing CRS:

  A layer without a CRS is read by the lon/lat test of
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
  ("Layers without a CRS"). With no `target_crs`, coordinates the test
  takes for degrees are read as EPSG:4326, with a warning, and the rules
  above then apply; coordinates it declines are left as they are. With
  `target_crs`, coordinates the test takes for degrees are read as
  EPSG:4326 and reprojected to the target, with a warning, and any
  others have the target stamped on without reprojection. A stamp is a
  relabel, logged only, so check that the coordinates are in that CRS.
  Set the CRS on `x` to avoid either. A `target_crs` in longitude and
  latitude is not stamped on coordinates that cannot be in it: an x
  outside -180 to 360 or a y beyond 90 north or south, in the box of the
  coordinates that have a location, stops the call with an error that
  gives the box.

- Rows without a location:

  Such rows (see
  [spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
  "Rows without a location") take no part in the lon/lat test or the
  choice of a projection. The lon/lat test reads the coordinates that
  have a location vertex by vertex, so a polygon with one vertex beyond
  1e150 is read on its other vertices. A layer in which no row has a
  location, or each row of which is empty, has no location or has a
  vertex at a latitude beyond 90, gets Web Mercator (Equal Earth for
  `purpose = "area"`). Moved to another CRS, a row without a location
  comes back empty, in its own geometry type. A row that has a location
  and that the move cannot place in full comes back empty, or without
  the vertices the CRS cannot place, with a warning of class
  `spatialkit_rows_unplaced` that names the rows (`rows` on the
  condition). With the projection chosen here that happens only at a
  coordinate no CRS places, such as a latitude beyond 90; with a
  `target_crs`, wherever that CRS cannot place the layer.

`target_crs` overrides all of this. Pass it whenever you need a
specific, reproducible projection: comparing runs, matching an existing
layer, or fixing the units that
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)'s
`block_size` will be interpreted in.

A layer that mixes XY and XYZ (or XYM) geometries is returned in XY when
it is transformed, with a warning of class `spatialkit_zm_dropped`, and
the projection is chosen from the layer in XY; see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
"Mixed dimensions, Z and M". A layer that is not transformed is returned
with its own geometry.

## Equal-area check

Under `purpose = "area"`, input that already has a projected CRS is
returned as it is, but its area distortion over the extent is measured
and logged as a warning when it exceeds 1 percent. The measure is the
spread of planar-to-geodesic area ratios over probe polygons. The probes
are the layer's distinct polygons (exact copies of a polygon count
once). A 6 by 6 grid over the bounding box is used instead for a layer
without polygons, a layer of one distinct polygon, which gives one ratio
and no spread, and a layer with a polygon s2 cannot measure in degrees
(edges that cross, a spike) or that
[`sf::st_is_valid()`](https://r-spatial.github.io/sf/reference/valid.html)
rejects in its CRS (a hole outside its shell), whose planar area is not
the area it encloses. When the grid gives no figure either (a single
point, a line along one axis), a warning says that the distortion could
not be measured. An equal-area projection still shows some spread, from
the sphere the geodesic areas are computed on against the ellipsoid the
projection uses.
[`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
applies the same check before it computes a density.

## See also

[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
for how input geometry and CRSs are read.

Other data preparation:
[`[.spatialkit_rows()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md),
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
[`spatialkit-geometry`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)

## Examples

``` r
library(sf)
pts_ll <- st_as_sf(
  data.frame(lon = c(9.1, 9.2), lat = c(48.7, 48.8)),
  coords = c("lon", "lat"), crs = 4326
)
# A local extent gets the containing UTM zone; the zone's measured
# distance error over the extent travels with the result.
st_crs(ensure_projected(pts_ll))$epsg  # 32632
#> [1] 32632
attr(ensure_projected(pts_ll), "crs_choice")
#>          name        crs distance_error chosen
#> 1 UTM zone 32 EPSG:32632   0.0003984524   TRUE

# A continental extent is scored against the zone and may get an
# equal-area projection instead.  The choice is logged with its distance
# error and the zone's, not raised as an R warning (see ?spatialkit_quiet);
# crs_choice holds every candidate scored.
wide <- st_as_sf(
  data.frame(lon = c(-120, -70), lat = c(30, 48)),
  coords = c("lon", "lat"), crs = 4326
)
wide_p <- ensure_projected(wide)
st_crs(wide_p)$proj4string
#> [1] "+proj=aea +lat_0=41.729902 +lon_0=-98.42263 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
attr(wide_p, "crs_choice")[, c("name", "distance_error", "chosen")]
#>                                                      name distance_error chosen
#> 1                                             UTM zone 14   0.0175763956  FALSE
#> 2   Lambert azimuthal equal-area centred on (-98.4, 41.7)   0.0055662736  FALSE
#> 3 Albers equal-area (lat_1=33.0, lat_2=45.0, lon_0=-98.4)   0.0001446928   TRUE

# target_crs overrides the choice.
st_crs(ensure_projected(pts_ll, target_crs = 3035))$epsg  # 3035
#> [1] 3035

# For densities per cell the CRS has to be equal-area.  Two candidates are
# scored, an Albers conic and a Lambert azimuthal, both centred on the data,
# and whichever distorts distance less over the extent is used: here Albers.
st_crs(ensure_projected(pts_ll, purpose = "area"))$proj4string
#> [1] "+proj=aea +lat_0=48.750011 +lon_0=9.14995 +lat_1=48.716667 +lat_2=48.783333 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
```
