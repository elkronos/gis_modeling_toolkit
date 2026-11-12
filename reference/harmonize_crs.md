# Bring two spatial objects to one CRS

Brings two sf or sfc objects to one CRS: the CRS of `a` or `b`, as
`prefer` says, or `target_crs` when it is given. Without `target_crs`,
two objects that both lack a CRS are returned unchanged.

## Usage

``` r
harmonize_crs(
  a,
  b,
  prefer = c("a", "b"),
  target_crs = NULL,
  on_transform_error = c("stop", "set_crs")
)
```

## Arguments

- a, b:

  Objects of class sf or sfc.

- prefer:

  Which object's CRS both are brought to when they differ: `"a"`
  (default) or `"b"`. Ignored when `target_crs` is given.

- target_crs:

  Optional target CRS to apply to both: anything
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  accepts, including an sf or sfc object, whose CRS is used.

- on_transform_error:

  What to do when
  [`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
  fails: `"stop"` (default) raises an error; `"set_crs"` relabels the
  CRS with
  [`sf::st_set_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  without moving the coordinates. Use `"set_crs"` only when the
  coordinates are already in the target CRS.

## Value

A named list with components `a` and `b`, each of the class it came in.
A layer that is not transformed is returned as it came.

## Details

When one input has no CRS, the lon/lat test that
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
uses decides what happens (see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
"Layers without a CRS"), so both functions place the same data in the
same place: coordinates the test takes for degrees are read as EPSG:4326
and reprojected to the other object's CRS (or `target_crs`); other
coordinates are stamped with that CRS by
[`sf::st_set_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html),
which relabels them without moving them. Either way a warning says so. A
CRS in longitude and latitude is not stamped on coordinates that cannot
be in it: a layer with an x outside -180 to 360 or a y beyond 90 north
or south, read on the coordinates that have a location, stops the call
with an error that names `a` or `b`, its box and the target CRS.

A layer that is transformed is read as
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
describes: a layer that mixes XY and XYZ (or XYM) geometries comes back
in XY, with a warning of class `spatialkit_zm_dropped`, and a row
without a location comes back empty, in its own geometry type. So does a
row of a projected layer moved to a CRS in degrees with a coordinate
beyond 1e10 in its CRS's own unit, further out than any projected CRS
places a point of the Earth. A row that has a location and that the
target CRS cannot place in full, such as a point far outside a UTM zone
or a ring with a vertex at a latitude beyond 90, comes back empty, or
without the vertices the CRS cannot place, and a warning of class
`spatialkit_rows_unplaced` names such rows of `a` or `b` (`rows` on the
condition).

## See also

[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
for how input geometry and CRSs are read.

Other data preparation:
[`[.spatialkit_rows()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md),
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md),
[`spatialkit-geometry`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)

## Examples

``` r
library(sf)
a <- st_as_sf(data.frame(x = c(500000, 500100), y = c(4000000, 4000100)),
              coords = c("x", "y"), crs = 32632)
b <- st_transform(a, 4326)                 # same points, lon/lat
h <- harmonize_crs(a, b)                    # b is brought into a's CRS
c(a = st_crs(h$a)$epsg, b = st_crs(h$b)$epsg)
#>     a     b 
#> 32632 32632 
st_crs(h$a) == st_crs(h$b)
#> [1] TRUE
```
