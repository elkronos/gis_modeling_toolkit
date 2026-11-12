# Extract the coefficients of a fitted model

[`coef()`](https://rdrr.io/r/stats/coef.html) on a GWR returns the local
coefficients, one row per observation and one column per term. On a
Bayesian GP fit it returns the posterior summary of the global
regression terms. A random forest has no coefficients, so
[`coef()`](https://rdrr.io/r/stats/coef.html) on an `rf_fit` always
stops with an error rather than returning `NULL`; a forest's variable
importance is in `fit$info$importance` (see
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)).
None of the three returns `NULL` (see "The coef() contract" in
[`new_spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/new_spatial_fit.md)).

## Usage

``` r
# S3 method for class 'gwr_fit'
coef(object, se = FALSE, ...)

# S3 method for class 'bayesian_fit'
coef(object, ...)

# S3 method for class 'rf_fit'
coef(object, ...)
```

## Arguments

- object:

  A fitted `gwr_fit`, `bayesian_fit` or `rf_fit`.

- se:

  GWR only. Logical; also return the standard errors and t values.
  Default `FALSE`.

- ...:

  Not used.

## Value

For a GWR with `se = FALSE`, a data.frame of local coefficient
estimates: one row per observation, one column per model term. With
`se = TRUE`, a list of three such data.frames: `estimate`, `se` and `t`.
An engine without an `SDF` component is an error.

For a Bayesian GP fit, a matrix of fixed-effect posterior summaries, as
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
returns them, on the fitted scale (per standard deviation of each
predictor under `standardize_predictors = TRUE`), in the units of the
response: where the response was handed to brms less a centre or divided
by a scale, the intercept has the centre added back and every column is
multiplied by the scale (see "The response's centre and scale" in
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) on
`object$engine` answers on the scale brms was given. A missing brms or a
failing `fixef()` call is an error.

For a random forest, nothing: the call is an error.

## Details

The GWR table shows where, and by how much, a predictor's relationship
with the response changes across the study area: read the spread of its
column, or join the table back to `object$data_sf` to map it.
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md)
draws that map.

The Bayesian GP summary gives the estimate, error and credible interval
of each predictor, as
[`brms::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
reports them. The Gaussian-process term has already taken up the
spatially structured part of the signal, so these are effects net of
location (see "Spatial confounding" in
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).

## GWR columns

Only the model terms are returned: the intercept and one column per
predictor. The engine's `SDF` data slot holds more beside them (standard
errors, t values, the observed response, the fitted values, the
residuals, `Local_R2`); read `object$engine$SDF` for those, and ask for
the standard errors and t values with `se = TRUE`.

A coefficient is `NA` at a location whose local window cannot estimate
it: a predictor that does not vary among the observations that carry
weight there has no slope, and the intercept is then not identified
either unless the predictor is 0 there. `object$info$nonfinite_coef` is
the mask, and
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
raises a warning that counts them.

If the model terms cannot be located in the `SDF` (an engine that names
its coefficient columns differently), the whole slot is returned, with a
WARN line in the log saying so (a log record, not an R warning; see
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)).
With `se = TRUE` the call stops instead. It never returns a silently
short table.

## GWR standard errors

With `se = TRUE` the local standard errors and t values come with the
estimates. A standard error is \\\hat\sigma\sqrt{(C C^\top)\_{jj}}\\,
with \\C = (X^\top W X)^{-1} X^\top W\\ the weights the local regression
puts on the responses and \\\hat\sigma^2\\ the residual sum of squares
over the fit's effective degrees of freedom: the definition
[`GWmodel::gwr.basic()`](https://rdrr.io/pkg/GWmodel/man/gwr.basic.html)
uses, computed from the QR decomposition that gave the estimate. They
treat the bandwidth as known. The t values are one test per term and
location, tens to thousands of them and strongly dependent, so a map of
\\\|t\| \> 2\\ overstates what is significant; no p-value, adjusted or
not, is returned. An adjustment for the number of tests can use the t
values here, the effective number of parameters `object$info$enp` and
the residual degrees of freedom `object$info$edf`. When observations
share a location the standard errors are too small as well (see
"Observations that share a location" in
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)).

## Bayesian GP standardised predictors

The summaries are on the scale the model was fitted on. A fit made with
`standardize_predictors = TRUE` was fitted on centred and scaled numeric
predictors, so each slope is the change in the linear predictor per
standard deviation of its predictor and the intercept is its value at
the predictor means, not the raw-unit numbers
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) reports on the same
formula. The returned matrix does not say so;
[`print()`](https://rdrr.io/r/base/print.html) on the fit does, and the
centre and scale of each predictor are in
`object$info$predictor_scaling`. To put a slope back in raw units divide
its `Estimate`, `Est.Error` and interval bounds by that predictor's
`scale`. The intercept's `Estimate` follows by linearity (subtract each
raw-unit slope times its predictor's `center`), but its `Est.Error` and
interval depend on the posterior covariance of the coefficients:
transform the draws from `brms::as_draws_df(object$engine)` for those,
or refit without standardising.

## See also

Other fit methods:
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`spatial_fit-fitted`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md),
[`spatial_fit-predict`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-predict.md),
[`spatial_fit-residuals`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 80
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             a = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
# The slope on a rises from west to east.
east <- (st_coordinates(pts)[, 1] - 5e5) / 1000
pts$z <- 10 + (1 + 2 * east) * pts$a + rnorm(n, 0, 0.3)
gw <- fit_gwr_model(pts, "z", "a", bandwidth = 30)
str(coef(gw, se = TRUE))
#> List of 3
#>  $ estimate:'data.frame':    80 obs. of  2 variables:
#>   ..$ Intercept: num [1:80] 10.09 9.9 10.06 10.21 9.99 ...
#>   ..$ a        : num [1:80] 1.53 1.62 2.17 2.74 1.45 ...
#>  $ se      :'data.frame':    80 obs. of  2 variables:
#>   ..$ Intercept: num [1:80] 0.0757 0.0736 0.079 0.0786 0.0734 ...
#>   ..$ a        : num [1:80] 0.0792 0.0852 0.0802 0.0845 0.0875 ...
#>  $ t       :'data.frame':    80 obs. of  2 variables:
#>   ..$ Intercept: num [1:80] 133 135 127 130 136 ...
#>   ..$ a        : num [1:80] 19.3 19 27 32.4 16.6 ...
summary(coef(gw)$a)     # the local slopes span the west-east change
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.325   1.596   1.967   2.075   2.603   2.780 
```
