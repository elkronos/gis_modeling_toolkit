# Prepare an sf layer for model fitting

Reduces the layer to 2-D points in a projected CRS and drops the rows a
model cannot use: rows with a missing or non-finite value in the
response or a predictor, and rows without a location. Every fitter and
every `cv_*()` function prepares its data here.

## Usage

``` r
prep_model_data(
  data_sf,
  response_var,
  predictor_vars,
  boundary = NULL,
  pointize = c("auto", "surface", "point_on_surface", "centroid", "line_midpoint",
    "bbox_center"),
  require_response = TRUE
)
```

## Arguments

- data_sf:

  An sf object. Rows without a location are dropped (see Details).

- response_var:

  Name of the response column.

- predictor_vars:

  Predictor column names. May be `character(0)` for an intercept-only
  model
  ([`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
  supports one; the GWR and random-forest backends do not and reject it
  themselves).

- boundary:

  Optional sf/sfc whose CRS the data are projected to. It is not a clip:
  rows that fall outside it are kept. Subset `data_sf` yourself (for
  example with
  [`sf::st_filter()`](https://r-spatial.github.io/sf/reference/st_join.html))
  to restrict the rows. One that mixes XY and XYZ (or XYM) geometries is
  read in x and y.

- pointize:

  How non-point geometries become points, passed to
  [`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md)
  as its `mode`: `"auto"` (default), `"centroid"`, `"point_on_surface"`,
  `"surface"` (an alias for `"point_on_surface"`), `"line_midpoint"` or
  `"bbox_center"`. Its errors, such as the one of class
  `spatialkit_features_degenerate` for a polygon s2 cannot take a
  centroid of, are raised in this function's name and give the rows of
  `data_sf`.

- require_response:

  Logical; if `FALSE` the response column is not required to be present
  (useful for out-of-sample prediction where the response is unknown).
  Default `TRUE`.

## Value

An sf object with 2-D (XY) POINT geometry, without the rows that were
dropped. The attribute `"dropped"` records what was removed and two
marks on the rows that were kept, in a list with these fields:

- `n`:

  Rows dropped. Every fit stores it as `$info$n_dropped`.

- `n_geometry`:

  How many of them were dropped for having no location.

- `which`:

  Their positions in `data_sf`.

- `row_id`:

  Their `..row_id` values when the layer carries that column, else
  `NULL`.

- `reason`:

  One per dropped row: `"geometry"`, `"missing"` or `"non_finite"`, in
  that order of precedence when several apply.

- `by_column`:

  For each modelling column, and for `"(geometry)"`, how many of the
  dropped rows it is missing or non-finite in; a row can count under
  several.

- `response_std_diff`:

  The mean response of the dropped rows minus that of the kept rows, in
  standard deviations of the response over both; `NA` when the response
  is not numeric, is itself missing in every dropped row, or takes one
  value up to the rounding of its values.

- `n_rows`:

  The number of rows returned, which the record was made for.

- `constant_predictors`:

  The predictors that take one value in every returned row
  (`character(0)` when none does).

- `n_duplicate_locations`:

  The returned rows whose coordinates repeat those of an earlier row.

## Details

Non-POINT geometries (MULTIPOINT included) are reduced to representative
points by
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
so that coordinates line up one row per observation. Any Z or M
coordinate (POINT Z from a GPS, a GeoPackage or KML) is dropped, because
every backend works in 2-D distance.

A row has no location when its geometry is empty, a coordinate is not
finite or lies beyond 1e150, or, in a geographic CRS, a vertex lies at a
latitude beyond 90 degrees (see
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md),
"Rows without a location"). The 1e150 bound is read on the coordinates
as given and again after the projection, because a move between CRSs can
bring a coordinate beyond it back finite, and that point is not a place
the caller gave.

The CRS is projected whenever one can be established
([`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)).
A layer without a CRS is read by the lon/lat test of
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)
("Layers without a CRS"): coordinates the test takes for degrees are
read as EPSG:4326 and projected, with a warning, and others are passed
through as they are. A small planar survey inside the lon/lat envelope
is read as degrees, so set the CRS on `data_sf` if the data are planar.

The response may not appear in `predictor_vars`. Using it as its own
predictor is leakage no backend catches (an out-of-bag \\R^2\\ near 1 in
the random forest, a silently reduced design matrix in GWR, duplicated
rows and a phantom `<none>` entry in the GWR selection table), so it is
refused here.

Column names must be syntactically valid R names (`make.names(x) == x`).
Every backend builds a model formula from these names, and a name R
parses as an expression would fit another model than the one requested
while the fit still recorded the name you asked for: `"B5-B4"` is
`B5 - B4`, and `"log(a)"` is a function call. Rename the column (for
example with [`make.names()`](https://rdrr.io/r/base/make.names.html))
before fitting.

The `"dropped"` record describes the rows this call returned and does
not survive subsetting: `clean[i, ]`, like
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
`slice()` or `arrange()` of it, is a plain layer with no `"dropped"`
attribute, and a fit given such a subset with `.already_prepped = TRUE`
reports `n_dropped = 0` even when the parent layer dropped rows.
[`sf::st_drop_geometry()`](https://r-spatial.github.io/sf/reference/st_geometry.html)
keeps the record, since the rows are the same; see
[`[.spatialkit_rows`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md)
for what binding such data frames does.

## Warnings

When preparing data for a fit (`require_response = TRUE`) three
conditions are raised as R warnings, each with a class so it can be
caught by name with
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) or
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html):

- `spatialkit_rows_dropped`:

  Rows were dropped. The message gives the number and share of rows, the
  count by reason and by column, and, for a numeric response, how far
  the mean response of the dropped rows sits from that of the kept rows.
  A fit on the rows that remain describes those rows only; when the
  difference is large the rows are not missing at random and the fit is
  biased towards the complete ones.

- `spatialkit_response_imprecise`:

  The response lies so far from zero that doubles there are spaced more
  than a millionth of its standard deviation apart (a standard deviation
  of 1 at about 1e10, \\2^{33}\\, or beyond). Its values are rounded to
  that spacing, and so is every mean and prediction computed at that
  size: the fit, its scores and any choice made on them read that
  rounding. The rows are kept; subtract a constant from the response
  first.

- `spatialkit_constant_predictor`:

  A predictor takes a single value in every row that remains. Its
  coefficient cannot be told from the intercept: a GWR window is
  rank-deficient, a Bayesian slope is its prior, and a forest ignores
  it. The predictor is kept; drop it from `predictor_vars`.

The fitters and the `cv_*()` functions pass these warnings on.
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
and
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
run several cross-validations of one layer, each of which prepares it,
and raise `spatialkit_rows_dropped` and `spatialkit_constant_predictor`
once per distinct text for the run and `spatialkit_response_imprecise`
once per response column;
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
raises the last once as well when its `fit_fn` prepares each fold's rows
again.

With `require_response = FALSE` (how every
[`predict()`](https://rdrr.io/r/stats/predict.html) method prepares
`newdata`) none is raised: the dropped rows are written to the log and
to the record, and the prediction for each comes back `NA`, so nothing
is lost without a trace. Repeated locations are legitimate (repeat
visits to one site) and are counted in the record and the log, never
warned about.

## See also

Other data preparation:
[`[.spatialkit_rows()`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md),
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md),
[`spatialkit-geometry`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)

## Examples

``` r
library(sf)
dat <- st_as_sf(
  data.frame(x = 1:5, y = 5:1,
             resp = c(1, 2, NA, 4, 5),
             pred = c(1, 2, 3, 4, Inf)),
  coords = c("x", "y"), crs = 32632
)
# Drops rows 3 (NA) and 5 (Inf), and warns that it did, with the counts.
clean <- suppressWarnings(prep_model_data(dat, "resp", "pred"))
attr(clean, "dropped")
#> $n
#> [1] 2
#> 
#> $n_geometry
#> [1] 0
#> 
#> $which
#> [1] 3 5
#> 
#> $row_id
#> NULL
#> 
#> $reason
#> [1] "missing"    "non_finite"
#> 
#> $by_column
#> resp pred 
#>    1    1 
#> 
#> $response_std_diff
#> [1] 1.460593
#> 
#> $constant_predictors
#> character(0)
#> 
#> $n_duplicate_locations
#> [1] 0
#> 
#> $n_rows
#> [1] 3
#> 
```
