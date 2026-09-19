# Subset a layer that carries a row record

[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
return a layer with an attribute recording what happened to its rows
(`"dropped"` and `"ties"` respectively). Those records describe the rows
the layer was built with, and `[` on an `sf` object copies attributes
through unchanged, which would leave a subset reporting its parent's
numbers with row positions that no longer resolve. Subsetting therefore
returns a plain layer with the record removed; read the record from the
layer the function returned, before subsetting it.

## Usage

``` r
# S3 method for class 'spatialkit_rows'
x[...]
```

## Arguments

- x:

  A layer returned by
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  or
  [`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md).

- ...:

  Passed to the underlying `sf` or data frame method.

## Value

The subset, without the row records and without this class. A subset
that is not a data frame (a single column taken with `drop = TRUE`) is
returned unchanged.

## Examples

``` r
library(sf)
dat <- st_as_sf(
  data.frame(x = 1:5, y = 5:1,
             resp = c(1, 2, NA, 4, 5), pred = c(1, 2, 3, 4, Inf)),
  coords = c("x", "y"), crs = 32632
)
clean <- prep_model_data(dat, "resp", "pred")
attr(clean, "dropped")$n          # 2
#> [1] 2
attr(clean[1:2, ], "dropped")     # NULL -- the record does not follow
#> NULL
```
