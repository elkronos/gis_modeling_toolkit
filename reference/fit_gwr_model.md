# Fit a Geographically Weighted Regression (GWR)

Fits a Gaussian GWR on an sf dataset with an adaptive or a fixed
bandwidth. The kernels, the bandwidth conventions and the diagnostics
are those of the GWmodel package. The local regressions are solved by
this package, each as a weighted least squares by pivoted QR, so a
window whose design is rank deficient neither stops the fit nor returns
coefficients made of rounding error.

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
  once). `character(0)` fits the model with a local intercept alone: the
  fitted value at a location is the kernel-weighted mean of the response
  around it. That is the baseline a GWR with predictors has to beat,
  since any GWR beats a global mean through its local intercept, and
  [`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
  fits it as the null model of a GWR sweep.

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
  beforehand to fix the units yourself. If NULL (default), the bandwidth
  with the lowest AICc is searched for (the lowest
  leave-one-location-out error when observations share a location); see
  **Bandwidth selection**. With `adaptive = FALSE`, a bandwidth smaller
  than a ten-thousandth of the data's extent raises a warning naming the
  extent and the CRS the fit runs in: every local window then holds its
  own observation alone. With `adaptive = TRUE` the count is rounded,
  and one too small for the model is raised, with a warning, to the
  smallest that gives every local regression more points of non-zero
  weight than parameters: the number of predictors plus 3 for the
  bisquare and tricube kernels (which give the farthest neighbour in a
  window weight 0), plus 2 for the others. That floor is enough unless
  several neighbours tie at the kernel's edge (a regular grid), which
  leaves a window fewer weighted points; then use a larger bandwidth. A
  count above the number of observations is capped at it, with a warning
  (a distance meant for `adaptive = FALSE`, most often).

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
AICc (`NA`, with a warning, where AICc is undefined: its effective
number of parameters \\tr(S)\\ is not below \\n - 2\\, the local
regressions all but interpolate the data, and the large negative value
the formula gives would rank the fit above any other; also where the fit
reproduces the mean of the observations at each location, see
**Observations that share a location**), `bandwidth_is_fallback` (`TRUE`
when no bandwidth could be selected: the search failed and the arbitrary
fallback was used, or no bandwidth has a defined criterion and the
widest was used), `condition_index` (the global index),
`local_collinearity` (one row per observation: `row`, `x`, `y`,
`n_window`, `cn` and `cn_slopes`, see **Collinearity diagnostics**, then
`rank`, the rank the local solve found, and `leverage`, the
observation's diagonal of the hat matrix), `n_local_collinear` (the
locations whose slopes count as collinear), `n_local_singular` (the
locations with a local coefficient that is `NA`), `n_rank_deficient`
(the windows whose rank is below the number of parameters),
`n_saturated` (the windows whose fitted value is their own observation:
leverage 1), `nonfinite_coef` (a logical matrix, one row per observation
and one column per term, with `Intercept` first, `TRUE` where the local
coefficient is `NA`, so the count in `n_local_singular` can be placed),
`n_dropped` (the rows
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)
removed for missing or non-finite values or a bad geometry, so `$n` can
be read against `nrow(data_sf)`), `n_locations` (the distinct locations
of the rows fitted; see **Observations that share a location**), `tr_S`,
`enp`, `edf` and `sigma2` (the trace of the hat matrix, the effective
number of parameters \\2\\tr(S) - tr(S^\top S)\\, the effective residual
degrees of freedom and the residual variance, as GWmodel defines them),
`global_ols` (`RSS`, `RMSE` and `AICc` of one least-squares regression
on the same rows, the AICc by the same formula with the regression's
rank as its trace; a GWR whose AICc is no lower has not shown that its
coefficients vary, and
[`summary()`](https://rdrr.io/r/base/summary.html) prints the two side
by side), `rank_tolerance`, `max_kernel_width` (the widest kernel at any
observation, which [`predict()`](https://rdrr.io/r/stats/predict.html)
measures new locations against) and `bandwidth_search` (see **Bandwidth
selection**; `NULL` when a bandwidth was supplied).

`$engine` holds the fit in the layout of the object
[`GWmodel::gwr.basic()`](https://rdrr.io/pkg/GWmodel/man/gwr.basic.html)
returns, under the same names: `SDF` (a `SpatialPointsDataFrame` of the
local coefficients, `y`, `yhat`, `residual`, `CV_Score`,
`Stud_residual`, the standard errors `<term>_SE`, the t values
`<term>_TV` and `Local_R2`), `GW.diagnostic`, `GW.arguments` and `lm`
(the global least-squares fit). It is a plain list, not a `gwrm` object.

## The local solve

At every observation the model is a weighted least squares on the
observations whose kernel weight is above zero. It is solved by a
Householder QR with column pivoting of the window's weighted design,
never by inverting \\X^\top W X\\. Before the factorisation each row is
multiplied by the square root of its weight over the window's largest
weight, and each predictor is divided by its standard deviation over the
whole data and centred in the window, so nothing below depends on a
predictor's origin or units: an indicator coded 0/1 and one coded
288.7/289.3 give the same fit.

One tolerance, `1e-7`, decides what a window can estimate:

- **Rank.** A predictor is aliased in a window when the part of it that
  the terms entered before it do not explain varies there by less than
  `1e-7` of its standard deviation over the study area (the length of
  that part over the weighted rows). A predictor that is constant among
  the weighted observations, two predictors that are exactly collinear
  there, and a window with fewer weighted observations than parameters
  are the cases. No weight is cut off: under a Gaussian or exponential
  kernel every observation enters, and distant ones stop counting only
  through this rule.

- **Coefficients.** A coefficient is returned when it has one value over
  all least-squares solutions of the window, and is `NA` otherwise. The
  slope of an aliased predictor is `NA`. So is the intercept when a
  predictor is constant in the window at a value other than 0, because
  it would absorb that value times the slope; when the constant is 0 the
  intercept is identified and returned. The slopes of the other
  predictors are returned.

- **Fitted values.** The fitted value at an observation is always
  estimable (its own row is in its window) and is returned, also where a
  coefficient is `NA`.

A window of deficient rank does not stop the fit. A warning of class
`"spatialkit_gwr_rank_deficient"` says how many windows are affected and
how many coefficients of each term are `NA`; `info$n_rank_deficient`,
`info$n_local_singular`, `info$nonfinite_coef` and the `rank` column of
`info$local_collinearity` locate them. When many windows hold a
predictor constant, a local coefficient is the wrong model for it: its
effect is identified only as one global coefficient (a mixed GWR, as
[`GWmodel::gwr.mixed()`](https://rdrr.io/pkg/GWmodel/man/gwr.mixed.html)
fits).

The leverage of each observation, the trace of the hat matrix, the
residual sum of squares, AICc and the coefficients' standard errors come
from the same decomposition. On a window of full rank every one of them
equals what
[`GWmodel::gwr.basic()`](https://rdrr.io/pkg/GWmodel/man/gwr.basic.html)
returns, to rounding error: the largest relative difference in any cell
of the six fits the package's tests compare is `3e-12`. They are not
bit-identical, since GWmodel inverts the normal equations.

## Bandwidth selection

With `bandwidth = NULL` the AICc is scored on 16 bandwidths spaced
evenly in the logarithm over the whole range, and then refined by a
golden-section search between the grid points either side of the best
(and, for a count of neighbours, at the counts next to the one it ends
on). That is about 20 to 35 fits of the model in all. The range runs,
for an adaptive bandwidth, from the smallest count the model's
parameters allow to every observation; for a fixed one, from the
distance at which half of the windows hold one observation more than the
model has parameters (a third of it for the Gaussian and exponential
kernels) to the largest distance between two observations. The bandwidth
returned is the best of all those scored.

A golden-section search over the whole range, which
[`GWmodel::bw.gwr()`](https://rdrr.io/pkg/GWmodel/man/bw.gwr.html) runs,
assumes the criterion has one minimum. It has several when the
coefficients vary at a range that only a few dozen observations resolve,
and that search then returns a bandwidth near the top of the range, a
nearly global regression, as the selected one.

A candidate bandwidth with windows of deficient rank is scored on its
fitted values, which those windows do estimate. One at which AICc is
undefined (\\tr(S)\\ not below \\n - 2\\) cannot be selected, and
neither can one at which a term of the model is estimable at no
location: the fit would have no coefficient for it anywhere.

The criterion at a candidate comes from the local solve. Under a
Gaussian or an exponential kernel, which weight every observation at
every location, the fitted value and leverage of a window that is far
from rank deficient are computed from its centred cross-products, many
windows to one matrix product, and the windows that are not are left to
the QR; the two routes agree to rounding error (`1e-9` relative in the
package's tests) and take the same rank decisions. A default-bandwidth
fit with one of those kernels and two predictors then takes 1.1 to 1.6
times as long as
[`GWmodel::bw.gwr()`](https://rdrr.io/pkg/GWmodel/man/bw.gwr.html)
followed by `gwr.basic()` (2.0 s against 1.6 s at 1,000 observations and
7.4 s against 5.7 s at 2,000, fixed Gaussian, on one machine), and 3 to
5 times as long by the QR alone. With several strongly correlated
predictors more windows are left to the QR and the fit is slower than
that. The fit at the selected bandwidth is always the QR's.

Everything scored is kept in `info$bandwidth_search`: `approach` (the
criterion minimised, `"AICc"` or `"CV"`), `by_location` and
`n_locations` (see **Observations that share a location**), `curve`
(bandwidth, criterion, trace of the hat matrix, windows of deficient
rank, and whether the point belongs to the grid or the refinement),
`undefined` (how many of the bandwidths scored have no criterion, by
reason), `n_local_minima` and `multimodal` (local minima on the grid),
`within` and `flat` (the range of bandwidths whose AICc is within 2 of
the minimum, and whether half the grid is), `near_global` (the widest
bandwidth of the range is within 2 of the minimum), `at_upper` and
`at_lower` (the minimum is at an end of the range). When any of these
holds, one logged line says so (on the console, unless
[`spatialkit_quiet()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit_quiet.md)
is on; in the trace only for the folds of a cross-validation); it is a
message, not a warning, because a minimum at or near the top of the
range is the ordinary outcome for a process whose coefficients do not
vary. Read `near_global` as "AICc does not tell this fit from one global
regression".

AICc is an in-sample criterion and rewards interpolation between
neighbours. The bandwidth it selects can score worse than a wider one on
held-out spatial blocks;
[`cv_gwr`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
with a supplied `bandwidth` is the check.

## Observations that share a location

Rows with the same coordinates (several years at one station, replicate
samples of one plot, records geocoded to one address) get one kernel
weight at every regression point. `info$n_locations` is the number of
distinct locations, and when it is below `n` a warning of class
`"spatialkit_gwr_repeated_locations"` gives the counts. Three things
follow.

**The bandwidth is not selected by AICc.** AICc counts the rows at one
location as independent observations. When a location has an effect of
its own, which is the usual case, a window that holds little but the
location reproduces its mean, the residual sum of squares falls to the
spread within locations, and AICc is lowest where the windows
interpolate: on 100 simulated sites with 3 observations each it is 58 at
5 to 10 neighbours, 90 at 13 and 221 at 58, while the error of
prediction at new sites is 46 at 10 neighbours, 1.35 at 13 and 0.15 at
150. A floor on the bandwidth does not repair that, because the
criterion keeps falling towards the floor. With `bandwidth = NULL` the
search therefore minimises the **leave-one-location-out** sum of
squares: each location's observations are predicted by the local
regression there without any of them, the kernel keeping its width, and
the squared errors are summed. It is computed from the hat matrix of the
fit at each bandwidth, without refitting. A bandwidth at which some
location cannot be predicted from the others cannot be selected.
`info$bandwidth_search$approach` is then `"CV"` and `by_location` is
`TRUE`. On a layer without repeat locations nothing changes: the
criterion is AICc.

**The range of an adaptive search counts locations.** It starts at the
smallest count that gives every window one location more than the model
has parameters, and for a fixed bandwidth at the distance at which half
the windows hold that many locations.

**AICc and the standard errors still assume independence.** `info$AICc`,
the comparison with the global regression in
[`summary()`](https://rdrr.io/r/base/summary.html) and the standard
errors of `coef(se = TRUE)` and `predict(se.fit = TRUE)` treat every row
as an independent observation. With an effect shared by the rows of a
location the standard errors are too small, and AICc compared across
bandwidths favours the small ones. `info$AICc` is `NA`, with a warning,
when the fit (nearly) reproduces the mean of the observations at each
location: \\tr(S)\\ is then not below the number of rows that differ in
location or in a predictor, minus 2, and the value the formula gives
measures only the spread among observations that share a location. If
the repeats are replicates of one measurement, fit their means, one row
per location.

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

**Local** indices are then computed at **every** location, on the
weighted design of the local regression there: each row weighted by the
square root of its kernel weight at the bandwidth the model is fitted
with (supplied or selected), with rows of negligible weight (below 1e-8)
dropped. A window left with fewer rows than columns counts as singular.
Two indices are kept for each window, as columns of
`info$local_collinearity`:

- `cn`:

  Belsley's index of the intercept plus the predictors, scaled to unit
  length but not centred. A predictor whose values in the window are far
  from 0 against their spread (a year, a temperature in kelvin) is
  collinear with the intercept and raises it: the local intercept is
  then an extrapolation to 0 and is ill-determined, but the slopes are
  not.

- `cn_slopes`:

  The index for the slopes: the predictors centred at their weighted
  mean in the window and each divided by its standard deviation over the
  whole study area. It is 1 when the predictors vary as much, and as
  independently, inside the window as they do across the study area; it
  grows as a predictor becomes nearly constant inside the window (a
  regional covariate) or two predictors move together there. It does not
  depend on the predictors' origin or units.

A window's slopes count as collinear when `cn_slopes` is above 30 or
singular. Those windows are counted in `info$n_local_collinear`. A
predictor that is constant, or nearly so, inside a window is caught this
way whether it is alone or has company, so a single predictor is
surveyed too. `cn` does not enter: the local solve centres each
predictor in its window, so a predictor far from zero (a year, a
temperature in kelvin, a coordinate) has the slopes it would have at any
other origin, and is not flagged, however large `cn` is. The one
exception is a `cn` that is singular to machine precision, where the
predictor's spread is lost in the rounding of its own values.

A warning is issued whenever **any** location has collinear slopes; the
wording reports a percentage when more than 25% of locations are
affected and a count otherwise. Both are real R warnings, not log lines.
Coefficients at a near-singular window are least-squares estimates with
little to rest on, and can be implausibly large. An **exactly** singular
window (an indicator that is constant inside it, or fewer observations
than parameters) is of deficient rank: the coefficients it cannot
estimate are `NA` and the fit is returned (see **The local solve**). A
window with only `cn` above 30 raises no warning:
`plot(fit, type = "coefficients", term = "Intercept")` masks it, and
slope maps do not. Centre such a predictor if you want an interpretable
local intercept. The survey is a caution about how well a coefficient is
determined, on rows whose weight is not negligible; the solve's rank is
the statement about whether it is determined at all.

A further warning counts local regressions that have no fitted value.
That happens where the kernel weights are undefined: with an adaptive
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
if (requireNamespace("sp", quietly = TRUE)) {
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
  head(coef(fit, se = TRUE)$se)       # local standard errors
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
#> 
#>   Effective parameters: 13.2 (residual degrees of freedom 46.8)
#>   Global least squares on the same rows: AICc 300.24, RMSE 2.7997
#>   This fit: AICc 218.78 (below the global fit's by 81.46)
#>   Intercept      elev
#> 1 0.3261779 0.3387700
#> 2 0.3248259 0.4382539
#> 3 0.3286578 0.3423202
#> 4 0.3203002 0.3282897
#> 5 0.3141750 0.3071722
#> 6 0.3241959 0.3390836
```
