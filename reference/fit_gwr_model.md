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

  Predictor column names (numeric columns; a name given twice counts
  once).

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
  anywhere. With `adaptive = TRUE` the count is rounded, and one too
  small for the model is raised, with a warning, to the smallest that
  gives every local regression more points of non-zero weight than
  parameters: the number of predictors plus 3 for the bisquare and
  tricube kernels (which give the farthest neighbour in a window weight
  0), plus 2 for the others. That floor is enough unless several
  neighbours tie at the kernel's edge (a regular grid), which leaves a
  window fewer weighted points; then use a larger bandwidth. A count
  above the number of observations is capped at it, with a warning (a
  distance meant for `adaptive = FALSE`, most often). `bw.gwr()`
  searches adaptive bandwidths from 20 neighbours up, so below 20
  observations its choice is capped the same way, with a warning, and is
  not an optimised bandwidth.

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
AICc (`NA`, with a warning, where GWmodel's AICc is undefined: its
effective number of parameters \\tr(S)\\ is not below \\n - 2\\, the
local regressions all but interpolate the data, and the large negative
value GWmodel reports would rank the fit above any other),
`bandwidth_is_fallback` (`TRUE` when automatic selection failed and the
arbitrary fallback was used), `condition_index` (the global index),
`local_collinearity` (one row per observation: `row`, `x`, `y`,
`n_window`, `cn` and `cn_slopes`; see **Collinearity diagnostics**),
`n_local_collinear` (the locations whose slopes count as collinear),
`n_local_singular` (the locations whose local coefficients came back
non-finite), `nonfinite_coef` (a logical matrix, one row per observation
and one column per term, with `Intercept` first, `TRUE` where the local
coefficient came back non-finite, so the count in `n_local_singular` can
be placed) and `n_dropped` (the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`). The raw GWmodel result is in
`$engine`.

## Collinearity diagnostics

The function computes **scaled condition indices** of the design and
warns when one exceeds 30, the conventional threshold, which Wheeler &
Tiefelsdorf (2005) carry over to the local designs of GWR. An index is
the ratio of the largest to the smallest singular value (from an SVD)
after each column is scaled (Belsley, Kuh & Welsch 1980); an exactly
singular design gives `Inf`, which counts as the worst case, not an
exempt one. Scaling makes the index independent of the predictors'
units; [`kappa()`](https://rdrr.io/r/base/kappa.html) on the raw matrix
is not, and a threshold on it is a threshold on nothing in particular.

A **global** index is computed on the predictors centred at their means
(with the intercept, which centring makes orthogonal to them), and kept
as `info$condition_index`. It measures how nearly the predictors are
collinear with one another over the whole study area; it is 1 for a
single predictor, and a change of origin (degrees C or kelvin, a year or
years since 2000) does not move it.

**Local** indices are then computed at **every** location, on the design
the local regression there inverts: each row weighted by the square root
of its kernel weight at the bandwidth the model is fitted with (supplied
or selected), with rows of negligible weight dropped. A window left with
fewer rows than columns counts as singular. Two indices are kept for
each window, as columns of `info$local_collinearity`:

- `cn`:

  Belsley's index of the intercept plus the predictors, scaled to unit
  length but not centred. A predictor whose values in the window are far
  from 0 against their spread (a year, a temperature in kelvin) is
  collinear with the intercept and raises it: the local intercept is
  then an extrapolation to 0 and is ill-determined, but the slopes are
  not. It is what GWmodel's own solve sees.

- `cn_slopes`:

  The index for the slopes: the predictors centred at their weighted
  mean in the window and each divided by its standard deviation over the
  whole study area. It is 1 when the predictors vary as much, and as
  independently, inside the window as they do across the study area; it
  grows as a predictor becomes nearly constant inside the window (a
  regional covariate) or two predictors move together there. It does not
  depend on the predictors' origin or units.

A window's slopes count as collinear when `cn_slopes` is above 30 or
singular, or when `cn` is above 1e6, where GWmodel's uncentred solve
starts to lose precision in the slopes too. Those windows are counted in
`info$n_local_collinear`. A predictor that is constant, or nearly so,
inside a window is caught this way whether it is alone or has company,
so a single predictor is surveyed too.

A warning is issued whenever **any** location has collinear slopes; the
wording reports a percentage when more than 25% of locations are
affected and a count otherwise. Both are real R warnings, not log lines.
Coefficients at a near-singular window are unstable and can be
implausibly large. An **exactly** singular window (an indicator that is
constant inside it, or fewer observations than parameters) makes GWmodel
stop, so the fit fails with an error that says so; the window is not
returned as `NaN`. A window with only `cn` above 30 raises no warning:
`plot(fit, type = "coefficients", term = "Intercept")` masks it, and
slope maps do not. Centre such a predictor if you want an interpretable
local intercept.

After the fit, the local coefficient surfaces are scanned and a further
warning counts local regressions that came back non-finite. GWmodel
returns those where the kernel weights are undefined: with an adaptive
bandwidth of `k`, a location where `k` or more observations share the
same coordinates has a kernel of zero width, and every kernel but the
boxcar divides 0 by 0 there. The warning names that cause when it
applies. [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
all drop those rows, so when this warning fires the metrics describe
only the part of the study area that fitted.

## See also

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md),
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
  print(summary(fit))     # print(): only a block's last value shows on its own
  head(predict(fit, newdata = dat))   # newdata is re-projected if needed
}
#> Summary of <gwr_fit> fit (n = 60)
#> 
#>   Formula: price ~ elev
#> 
#>   In-sample metrics:
#>     RMSE    = 1.1889
#>     MAE     = 0.9995
#>     R^2     = 0.8583
#>     SMAPE   = 6.87%
#> [1] 16.47047 13.69566 17.02994 18.13618 11.49180 18.37175
```
