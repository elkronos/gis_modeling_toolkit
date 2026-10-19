# Prepare and sanitize an sf dataset for spatial modeling

Ensures point geometry, projected CRS, and removes rows with missing or
non-finite values in modeling columns, *including* rows whose geometry
is empty or whose coordinates are not finite, which no model backend can
use. All non-POINT geometries (including MULTIPOINT) are coerced to
representative points via
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
so downstream coordinate extraction always aligns one row per
observation. Any Z or M coordinate (POINT Z from a GPS, a GeoPackage or
KML) is dropped, because every backend works in 2-D map distance.

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

  An sf object.

- response_var:

  Response variable column name.

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
  to restrict the rows.

- pointize:

  Strategy for non-point geometry coercion, passed to
  [`coerce_to_points`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md).
  One of `"auto"` (default), `"centroid"`, `"point_on_surface"`,
  `"surface"` (an alias for `"point_on_surface"`), `"line_midpoint"` or
  `"bbox_center"`.

- require_response:

  Logical; if FALSE the response column is not required to be present
  (useful for out-of-sample prediction where the response is unknown).
  Default TRUE.

## Value

An sf object with 2-D (XY) POINT geometry, cleaned of rows carrying
missing or non-finite values in the modelling columns or in the
coordinates. What was removed is recorded on the attribute `"dropped"`,
a list with `n` (rows dropped), `n_geometry` (how many of them for an
empty or non-finite geometry), `which` (their positions in `data_sf`),
`row_id` (their `..row_id` values when the layer carries that column,
else `NULL`) and `reason` (one per dropped row: `"geometry"`,
`"missing"` or `"non_finite"`, in that order of precedence when several
apply), `by_column` (for each modelling column, and for `"(geometry)"`,
how many of the dropped rows it is missing or non-finite in; a row can
count under several), `response_std_diff` (the mean response of the
dropped rows minus that of the kept rows, in standard deviations of the
response over both; `NA` when the response is not numeric, or is itself
missing in every dropped row), plus `n_rows`, the number of rows
returned, which the record was made for. The same record carries two
marks on the rows that were kept: `constant_predictors`, the names of
the predictors that take one value in every returned row (`character(0)`
when none does), and `n_duplicate_locations`, the number of returned
rows whose coordinates repeat those of an earlier row. Every fit stores
`n` as `$info$n_dropped`. The record describes the rows this call
returned and does not survive subsetting: `clean[i, ]`, like
[`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html),
`slice()` or `arrange()` of it, is a plain layer with no `"dropped"`
attribute, and a fit given such a subset with `.already_prepped = TRUE`
reports `n_dropped = 0` even when the parent layer dropped rows.
[`sf::st_drop_geometry()`](https://r-spatial.github.io/sf/reference/st_geometry.html)
keeps the record, since the rows are the same; see
[`[.spatialkit_rows`](https://elkronos.github.io/gis_modeling_toolkit/reference/sub-.spatialkit_rows.md)
for what binding such data frames does. The CRS is projected whenever
one can be established. A CRS-less layer is decided by the lon/lat
heuristic (see
[`ensure_projected`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)):
if its bounding box fits the lon/lat envelope *and* it either spans more
than one unit on some axis or carries decimal-degree-like precision, it
is read as EPSG:4326 and projected, with a warning. A small planar
survey inside that envelope is included in that, deliberately; only
coordinates the heuristic declines are passed through as-is. Set the CRS
on `data_sf` if the data are planar.

## Details

The response may not appear in `predictor_vars`. Using it as its own
predictor is leakage no backend catches: an out-of-bag R^2 near 1 in the
random forest, a silently reduced design matrix in GWR, duplicated rows
and a phantom `<none>` entry in the GWR selection table. It is refused
here.

Column names must be syntactically valid R names (`make.names(x) == x`).
Every backend builds a model formula from these names. A name R parses
as an expression would fit a different model from the one requested
while the fit object still recorded the name you asked for: `"B5-B4"` is
`B5 - B4`, and `"log(a)"` is a function call. Rename the column (for
example with [`make.names()`](https://rdrr.io/r/base/make.names.html))
before fitting.

## Warnings

When preparing data for a fit (`require_response = TRUE`) two conditions
are raised as R warnings, each with a class so it can be caught by name
with [`tryCatch()`](https://rdrr.io/r/base/conditions.html) or
[`withCallingHandlers()`](https://rdrr.io/r/base/conditions.html):

- `spatialkit_rows_dropped`:

  Rows were dropped. The message gives the number and share of rows, the
  count by reason and by column, and, for a numeric response, how far
  the mean response of the dropped rows sits from that of the kept rows.
  A fit on the rows that remain describes those rows only; when the
  difference is large the rows are not missing at random and the fit is
  biased towards the complete ones.

- `spatialkit_constant_predictor`:

  A predictor takes a single value in every row that remains. Its
  coefficient cannot be told from the intercept: a GWR window is
  rank-deficient, a Bayesian slope is its prior, and a forest ignores
  it. The predictor is kept; drop it from `predictor_vars`.

Every fitter and every `cv_*()` function prepares its data here and
passes these on.
[`compare_models_cv()`](https://elkronos.github.io/gis_modeling_toolkit/reference/compare_models_cv.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
and
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
run several cross-validations of one layer, each of which prepares it,
and raise a warning of either class once per distinct text for the run.

With `require_response = FALSE` (how every
[`predict()`](https://rdrr.io/r/stats/predict.html) method prepares
`newdata`) neither is raised: the dropped rows are written to the log
and to the record, and the prediction for each comes back `NA`, so
nothing is lost without a trace. Repeated locations are legitimate
(repeat visits to one site) and are counted in the record and the log,
never warned about.

## See also

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)

Other spatial data preparation:
[`clip_target_for()`](https://elkronos.github.io/gis_modeling_toolkit/reference/clip_target_for.md),
[`coerce_to_points()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coerce_to_points.md),
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md),
[`harmonize_crs()`](https://elkronos.github.io/gis_modeling_toolkit/reference/harmonize_crs.md)

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
