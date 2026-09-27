# Subset a layer that carries a row record

[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md)
return a layer with an attribute recording what happened to its rows
(`"dropped"` and `"ties"` respectively). Those records describe the rows
the layer was built with, and `[` on an `sf` object copies attributes
through unchanged, which would leave a subset reporting its parent's
numbers for a different set of rows. Subsetting therefore returns a
plain layer with the record removed, and so do the dplyr verbs that
select or reorder rows
([`filter()`](https://rdrr.io/r/stats/filter.html), `slice()`,
`arrange()`, `distinct()`); read the record from the layer the function
returned, before subsetting it. Binding such layers
([`rbind()`](https://rdrr.io/r/base/cbind.html),
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html))
likewise returns a plain `sf` layer.

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

## Details

Each record carries `n_rows`, the number of rows it was made for.
[`sf::st_drop_geometry()`](https://r-spatial.github.io/sf/reference/st_geometry.html)
keeps the rows, and with them the record: it returns a data frame of
class `c("spatialkit_rows", "data.frame")`. Binding such data frames
with [`rbind()`](https://rdrr.io/r/base/cbind.html) or
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
keeps the first one's record and class, so the record then describes
only the first input's rows: its `n_rows` no longer equals
[`nrow()`](https://rdrr.io/r/base/nrow.html) of the result. The
package's own readers ignore a record whose `n_rows` does not match;
when reading `attr(x, "dropped")` or `attr(x, "ties")` yourself from a
layer that has been through such steps, check it the same way.

## See also

[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
and
[`assign_features_to_polygons()`](https://elkronos.github.io/gis_modeling_toolkit/reference/assign_features_to_polygons.md),
the two functions that attach the records this method removes.

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
