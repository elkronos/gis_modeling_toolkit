# Harmonize CRS between two spatial objects

Aligns two sf objects to a common CRS.

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

  Which object's CRS to keep ("a" or "b").

- target_crs:

  Optional target CRS to apply to both.

- on_transform_error:

  What to do when st_transform() fails: `"stop"` (default) raises an
  error immediately; `"set_crs"` falls back to st_set_crs() (UNSAFE:
  coordinates are NOT reprojected, only the CRS label is overwritten).
  The `"set_crs"` option exists only for rare edge cases where you are
  certain the coordinates already match the target CRS definition.

## Value

A named list with components a and b.

## Details

When one input carries no CRS, the same lon/lat heuristic
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
uses decides what happens, so both entry points place identical data in
the same place: coordinates that look like degrees are taken as
EPSG:4326 and **reprojected** to the other object's CRS (or
`target_crs`); coordinates that do not are *stamped* with
[`sf::st_set_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html),
which relabels without moving them. Either way the assumption is
announced with a warning.

## See also

Other spatial data preparation:
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

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
