# Fit a Geographically Weighted Regression (GWR) via GWmodel

Fits a GWR using GWmodel on an sf dataset with either adaptive or fixed
bandwidth.

## Usage

``` r
fit_gwr_model(
  data_sf,
  response_var,
  predictor_vars,
  adaptive = TRUE,
  bandwidth = NULL,
  kernel = c("bisquare", "gaussian", "tricube", "boxcar", "exponential"),
  .already_prepped = FALSE
)
```

## Arguments

- data_sf:

  An sf object with response, predictors, and geometries.

- response_var:

  Response column name.

- predictor_vars:

  Predictor column names.

- adaptive:

  Logical; use adaptive bandwidth. Default TRUE. When TRUE, bandwidth is
  an integer number of nearest neighbours. When FALSE, bandwidth is a
  fixed distance in CRS units.

- bandwidth:

  Optional numeric bandwidth value. For adaptive mode this is an integer
  (number of neighbours); for fixed mode a distance in the units of the
  **projected** CRS the fit runs in.
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  projects geographic input before the bandwidth is used, so 0.2
  supplied for lon/lat data is 0.2 metres, not 0.2 degrees. Read the CRS
  off `sf::st_crs(fit$data_sf)`, or pass `target_crs` to
  [`ensure_projected`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)
  beforehand to fix the units yourself. If NULL (default), bandwidth is
  selected automatically via
  [`GWmodel::bw.gwr()`](https://rdrr.io/pkg/GWmodel/man/bw.gwr.html).
  With `adaptive = FALSE`, a bandwidth smaller than a ten-thousandth of
  the data's extent raises a warning naming the extent and the CRS the
  fit runs in: every local window is then likely to be empty, which used
  to produce a fit whose coefficients were all `NaN` with nothing raised
  anywhere.

- kernel:

  Kernel function type. One of "bisquare" (default), "gaussian",
  "tricube", "boxcar", "exponential".

- .already_prepped:

  Logical (internal). If `TRUE`, skip the
  [`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
  call because the caller has already projected, coerced, and filtered
  the data. Used by the CV internals to avoid a redundant second pass on
  every fold. End users should leave this at the default `FALSE`.

## Value

A `gwr_fit` object (inherits from `spatial_fit`). Supports
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md).
Model-specific metadata lives in `$info`: bandwidth, adaptive, kernel,
AICc, `bandwidth_is_fallback` (`TRUE` when automatic selection failed
and the arbitrary fallback was used), `condition_index`,
`local_collinearity`, `n_local_collinear`, `n_local_singular`,
`nonfinite_coef` (a logical matrix, one row per observation and one
column per term, with `Intercept` first, `TRUE` where the local
coefficient came back non-finite, so the count in `n_local_singular` can
be placed) and `n_dropped` (the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`). The raw GWmodel result is in
`$engine`.

## Collinearity diagnostics

The function computes the **scaled condition index** of the design and
warns when it exceeds 30, the conventional threshold, which Wheeler &
Tiefelsdorf (2005) carry over to the local designs of GWR. The index is
the ratio of the largest to the smallest singular value after each
column is scaled to unit length (Belsley, Kuh & Welsch 1980). Scaling
makes the index independent of the predictors' units;
[`kappa()`](https://rdrr.io/r/base/kappa.html) on the raw matrix is not,
and a threshold on it is a threshold on nothing in particular. A
**global** index is computed on the full design (intercept plus
predictors). In addition, a **local** spot-check is performed at up to
30 locations: every location when there are 30 or fewer, otherwise 30
spread evenly over the extent (evenly spaced ranks of the observations
ordered by x, then y), so the diagnostic is reproducible, draws no
random numbers, does not depend on the row order of the data, and the
count is not configurable. For each sampled point the nearest neighbours
within the bandwidth window (the bandwidth the model is actually fitted
with, not a stand-in) are selected and the condition number of that
local design sub-matrix is evaluated. That sub-matrix is the predictors
**plus an intercept column**, matching the design GWmodel fits, and is
unweighted; the global condition number is computed on the predictors
alone, so the two numbers are not directly comparable. An indicator that
is constant inside a window is collinear with the intercept and with
nothing else, which is why the intercept has to be there. A non-finite
condition number counts as extreme:
[`kappa()`](https://rdrr.io/r/base/kappa.html) returns `Inf` for an
exactly singular design, which is the worst case, not an exempt one.

A warning is issued whenever **any** sampled location has a singular or
near-singular local design; the wording reports a percentage when more
than 25\\ R warnings, not log lines.

After the fit, the local coefficient surfaces are scanned and a further
warning counts local regressions that came back non-finite. Their
windows were singular.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
all drop those rows, so when this warning fires the metrics describe
only the part of the study area that fitted.

Because the local spot-check examines only a subset of locations, it may
not detect every problematic neighbourhood. Users working with highly
clustered data or near-collinear predictors should consider a full
local-collinearity audit as a post-fit diagnostic.

## See also

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

## Examples

``` r
if (requireNamespace("GWmodel", quietly = TRUE) &&
    requireNamespace("sp", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 60
  dat <- st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  dat$price <- 10 + 0.01 * (st_coordinates(dat)[, 1] - 5e5) +
    2 * dat$elev + rnorm(n)
  fit <- fit_gwr_model(dat, "price", "elev", bandwidth = 30)
  summary(fit)
  head(predict(fit, newdata = dat))   # newdata is re-projected if needed
}
#> [1] 16.47047 13.69566 17.02994 18.13618 11.49180 18.37175
```
