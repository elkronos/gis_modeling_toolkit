# Ensure an object has a projected CRS (with sensible defaults)

Coerces spatial objects to a projected coordinate reference system
suitable for distance/area calculations.

## Usage

``` r
ensure_projected(x, target_crs = NULL, purpose = c("distance", "area"))
```

## Arguments

- x:

  An sf or sfc object (other objects returned unchanged).

- target_crs:

  Optional target CRS (sf object, integer EPSG, or crs). Must resolve to
  a usable CRS via
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html);
  an unusable value (one that resolves to `NA_crs_`) raises an error, so
  `x` is never left silently unprojected.

- purpose:

  Which property the projection is for. `"distance"` (the default, and
  everything above): the candidate that distorts pairwise distances
  least, which is what ranges, block sizes, bandwidths and length-scales
  read off the coordinates. `"area"`: densities or rates per cell are
  going to be computed, so the CRS must be equal-area. For lon/lat input
  the choice is then made among equal-area projections only: a Lambert
  azimuthal centred on the data, or an Albers conic where its parallels
  do not degenerate, whichever distorts distances less. A UTM zone
  (conformal, not equal-area) never enters that comparison, and global
  coverage gets Equal Earth in place of Web Mercator. Already-projected
  input is still returned untouched, but its area distortion over the
  extent is measured (the spread of planar-to-geodesic area ratios over
  probe polygons) and logged as a warning when it exceeds 1 percent.
  Measured: a zone's own width edge to edge, 0.25 percent; the
  conterminous United States forced into one zone, 14 percent; Web
  Mercator over 2.5 degrees of latitude at 48N, 4 percent; an equal-area
  projection, a few tenths of a percent, which is the sphere the
  geodesic areas are computed on against the ellipsoid the projection
  uses.
  [`summarize_by_cell()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summarize_by_cell.md)
  applies the same measurement before it computes a density. Ignored
  when `target_crs` is given.

## Value

x, potentially with a new projected CRS. When a projection was chosen
here (lon/lat input, no `target_crs`) the result carries
`attr(x, "crs_choice")`: a data.frame with one row per projection
considered, holding `name`, `crs` (its definition), the measured
worst-case `distance_error` (relative, over sampled pairs; `NA` where it
could not be measured) and `chosen`. The figure the log line quotes for
the winner is therefore recoverable for every candidate, and is measured
for the single candidate on the paths where no comparison runs (a UTM
zone on a local extent, the equal-area projection chosen for a layer
straddling the antimeridian). It is `NULL` exactly when no local
projection was chosen here: input that already carried a projected CRS,
a `target_crs` you supplied, or an extent no local projection fits,
which falls back to Web Mercator or Equal Earth. CRS-less input
additionally carries `attr(x, "crs_assumed")`: `"EPSG:4326"` when the
lon/lat heuristic fired, `"none"` when it declined. That attribute is
also read on the way IN. An object already carrying `"none"` is returned
untouched, with the heuristic skipped, which is how a
[`predict()`](https://rdrr.io/r/stats/predict.html) method replays a
fit's negative decision so that a subset of the training rows is not
judged differently from the whole.

## Details

An object that already has a projected CRS is returned untouched. Only
geographic (lon/lat) input is transformed, and the CRS chosen depends on
the extent of the data. It is **not** always UTM:

- Local extents:

  The UTM zone containing the data's centre (EPSG:326xx north of the
  equator, EPSG:327xx south). Distances and areas are close to true over
  a few degrees of longitude, which is the case this package is usually
  in.

- Wide extents:

  Once the data reach well beyond the roughly 3 degrees a UTM zone is
  designed for, a single zone can distort distances by several percent,
  and that error propagates straight into variogram ranges, block sizes,
  GWR bandwidths and GP length-scales. Which projection is actually best
  is then **measured, not assumed**: the zone, a Lambert azimuthal
  equal-area centred on the data and (where its standard parallels do
  not degenerate) an Albers conic are each scored by projecting
  representative points of the data (a non-POINT layer is reduced to
  points first) and comparing planar with geodesic pairwise distances,
  and the one that distorts least is used. The choice, both error
  figures and this argument are **logged** (see the logging note under
  [`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md));
  they are not R warnings, so `tryCatch(warning = )` does not see them.

- Antimeridian:

  Data straddling ±180° have a bounding box wider than a hemisphere. The
  wrap is detected from the coordinates (one very large gap in the
  sorted longitudes) and an equal-area projection centred on the true
  extent is used. Only truly global coverage falls back to EPSG:3857.

- Missing CRS:

  With no `target_crs`, a bounding box that looks like lon/lat means
  EPSG:4326 is assumed (a real warning) and the rules above then apply;
  coordinates the heuristic declines are left exactly as they are. With
  `target_crs` supplied there is no source CRS to reproject from, so the
  same heuristic decides between two outcomes: lon/lat-looking
  coordinates are read as EPSG:4326 and reprojected to the target (a
  real warning), and anything else has the target **stamped on without
  reprojection**. That is a relabel, logged only, so verify the
  coordinates really are in that CRS. Set the CRS explicitly to suppress
  either.

`target_crs` overrides all of this. Pass it whenever you need a
specific, reproducible projection: comparing runs, matching an existing
layer, or fixing the units that
[`make_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/make_folds.md)'s
`block_size` will be interpreted in.

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

# A continental extent is scored against the zone and may get an equal-area
# projection instead; the choice and both error figures are LOGGED, not
# warned -- see Details and ?spatialkit_quiet.
wide <- st_as_sf(
  data.frame(lon = c(-120, -70), lat = c(30, 48)),
  coords = c("lon", "lat"), crs = 4326
)
st_crs(ensure_projected(wide))$proj4string
#> [1] "+proj=aea +lat_0=41.729902 +lon_0=-98.42263 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# target_crs overrides the choice entirely.
st_crs(ensure_projected(pts_ll, target_crs = 3035))$epsg  # 3035
#> [1] 3035

# For densities per cell the CRS has to be equal-area: a Lambert azimuthal
# centred on the data rather than the UTM zone.
st_crs(ensure_projected(pts_ll, purpose = "area"))$proj4string
#> [1] "+proj=aea +lat_0=48.750011 +lon_0=9.14995 +lat_1=48.716667 +lat_2=48.783333 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
```
