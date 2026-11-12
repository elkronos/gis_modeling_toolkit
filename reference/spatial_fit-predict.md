# Predict from a fitted model

[`predict()`](https://rdrr.io/r/stats/predict.html) has a method for
each backend: a random forest (`rf_fit`), a GWR (`gwr_fit`) and a
Bayesian GP (`bayesian_fit`) fit. Every method reads `newdata` the same
way: a layer without a CRS is read the way the training data were (see
"Layers without a CRS" in
[spatialkit-geometry](https://elkronos.github.io/gis_modeling_toolkit/reference/spatialkit-geometry.md)),
the layer is transformed to the CRS used at fitting time, a feature that
is not a point is reduced to one, and a row with a missing or non-finite
predictor or no location gets `NA`. A layer that mixes XY with XYZ or
XYM geometries is read in x and y. The response need not be present.
With `newdata = NULL` each method returns the fit's own fitted values:
out-of-bag predictions for a random forest, in-sample values for the
others.

## Usage

``` r
# S3 method for class 'gwr_fit'
predict(object, newdata = NULL, se.fit = FALSE, ...)

# S3 method for class 'bayesian_fit'
predict(
  object,
  newdata = NULL,
  summary = c("mean", "median"),
  type = c("epred", "predict"),
  draws = FALSE,
  ...
)

# S3 method for class 'rf_fit'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `rf_fit`, `gwr_fit` or `bayesian_fit`.

- newdata:

  An sf object carrying the predictors, or `NULL` for the fitted values.

- se.fit:

  GWR only. Logical; also return standard errors. Default `FALSE`.

- ...:

  For a random forest, passed to ranger's predict method; any other name
  is an error. Arguments that make ranger return a matrix
  (`predict.all = TRUE`, `type = "quantiles"`, `type = "se"` with
  `predict.all`) are refused, because the method returns one number per
  row; call `predict(fit$engine, data = ...)` for those. What ranger
  itself refuses, such as `type = "quantiles"` on a forest grown without
  `quantreg = TRUE`, is an error that gives ranger's reason. `seed`
  defaults to a constant, so the number of
  [`predict()`](https://rdrr.io/r/stats/predict.html) calls a script
  makes (through
  [`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)'s
  `chunk_size`, say) does not shift later random draws; it does not
  change a regression forest's predictions. `num.threads` defaults to
  `getOption("mc.cores", 1L)`. For a GWR or a Bayesian GP fit, not used:
  any argument here is an error, so a misspelt `newdata` cannot slip
  through.

- summary:

  Bayesian GP only. `"mean"` (default) or `"median"` over the posterior
  draws.

- type:

  Bayesian GP only. `"epred"` (default) for the expected value, without
  observation noise, or `"predict"` for posterior predictive draws, with
  it.

- draws:

  Bayesian GP only. `TRUE` returns the full matrix of posterior draws
  instead of a summary. Default `FALSE`.

## Value

For a random forest, a numeric vector aligned to `nrow(newdata)`, with
`NA` for rows dropped as incomplete (all `NA`, with a WARN line in the
log, when every row is). The forest in `fit$engine` predicts
`(y - fit$info$response_centre) / fit$info$response_scale` (see the
Value section of
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md));
this method puts the centre and scale back. A failure inside ranger's
predict method is an error.

For a GWR with `se.fit = FALSE`, a numeric vector aligned to
`nrow(newdata)`, with `NA` for rows dropped as missing or non-finite and
for locations whose local regression cannot be estimated (see "When a
GWR prediction is NA"). If the design matrix for `newdata` cannot be
built, every value is `NA` and a warning says why. When the training
data were read as longitude and latitude, a CRS-less `newdata` with an x
outside -180 to 360 or a y beyond 90 north or south stops the call. With
`se.fit = TRUE`, a list: `fit` (that vector), `se.fit` (the standard
error of the fitted value, \\\hat\sigma\sqrt{s^\top s}\\ with \\s\\ the
weights the local regression puts on the training responses), `se.pred`
(the standard error of a new observation there, \\\hat\sigma\sqrt{1 +
s^\top s}\\, the square root of
[`GWmodel::gwr.predict()`](https://rdrr.io/pkg/GWmodel/man/gwr.predict.html)'s
prediction variance), `residual.scale` (\\\hat\sigma\\, from the fit's
residual sum of squares over its effective degrees of freedom) and `df`
(those degrees of freedom). The standard errors treat the bandwidth as
known and the local model as true at each location; they do not measure
the bias of a local fit and are `NA` where the prediction is.

For a Bayesian GP fit, a numeric vector of length `nrow(newdata)`, or an
`n_draws` by `nrow(newdata)` matrix with `draws = TRUE`, in the units of
the response (see "The response's centre and scale" in
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)).
Rows with a missing or non-finite predictor or coordinate, and rows
carrying a level the model was not fitted with, are `NA` (an `NA` column
with `draws = TRUE`). If the posterior draw fails the result is all `NA`
(a one-row matrix for `draws = TRUE`) with a warning of class
`spatialkit_predict_failed` that gives the cause. An ordinal or
categorical family is an error under `type = "epred"`: its expected
value is a probability per response category. Use
`type = "predict", draws = TRUE` for posterior draws of the predicted
category, as category indices; the share of draws in each category
estimates its probability. Without `draws = TRUE`, `type = "predict"`
returns the mean (or median) category index, an expected rank for an
ordinal family and an error for
[`brms::categorical()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
whose categories have no order. With `newdata = NULL` the cached
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) values are
returned only for the default `summary = "mean"`, `type = "epred"`,
`draws = FALSE`; any other combination is recomputed against the
training data, because the cache holds the expected values only.

## Details

A random forest predicts with ranger's predict method. A GWR solves one
local regression at each new location, with the fit's kernel and
bandwidth (an adaptive bandwidth counts neighbours among the training
points), and returns \\x^\top\hat\beta(u)\\. The local regressions are
weighted least squares solved by the pivoted QR
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
uses: on a window of full rank the values equal those of
[`GWmodel::gwr.predict()`](https://rdrr.io/pkg/GWmodel/man/gwr.predict.html)
to rounding error, and at a training location with that row's own
predictors they equal
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html). Each location
stands alone, so one that cannot be estimated does not affect the
others. A Bayesian GP fit applies the coordinate scaling and predictor
standardisation stored at fit time and then calls
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
or
[`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html).

## Random forest predictors

A categorical predictor must not carry a level the forest was grown
without, meaning a level with no training rows, not only one absent from
[`levels()`](https://rdrr.io/r/base/levels.html): an ordinary subset, or
a spatial fold that holds out a whole class, keeps the unused level
while the forest has no split for it. An unseen level is an error. A
predictor that was numeric or logical when the forest was grown must
arrive numeric or logical: ranger would otherwise code a character
column as a factor and apply the numeric split thresholds to the codes,
so a character column is refused. (ranger reads a logical as 0/1, so
either form is accepted for one.)

## When a GWR prediction is NA

A prediction is returned when it is estimable from the training points
that carry weight at the new location, and is `NA` otherwise:

- no training point carries weight there (a location beyond a fixed
  bisquare, tricube or boxcar bandwidth);

- the training points in the window cannot place the new row: a
  predictor is constant among them at a value other than the new row's
  (a regional indicator, at a location whose window lies wholly in
  another region), or fewer of them carry weight than the new row needs.

Where a predictor is constant in the window at the new row's own value,
its coefficient is not needed and the prediction is returned. The
decision uses the rule and the tolerance of the fit (see "The local
solve" in
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)).
A warning of class `"spatialkit_gwr_not_estimable"` counts the `NA`
predictions; it carries `n_na`, `n` and `n_empty`.

## GWR locations beyond the reach of the fit

An adaptive bandwidth always finds its neighbours, and a Gaussian or an
exponential kernel weights every training point, so a location far from
the data still gets a number: the local regression of training points
that are all far away, extrapolated. When a new location's nearest
training point is farther than the widest kernel used at any training
location reaches, the prediction is returned and a warning of class
`"spatialkit_gwr_extrapolation"` counts those locations and gives the
range of their distances. The widest kernel is the largest distance from
a training point to its `bandwidth`-th nearest neighbour
(`object$info$max_kernel_width`), or the bandwidth itself when it is
fixed. A bisquare, tricube or boxcar kernel reaches as far as its width;
a Gaussian or exponential kernel has no edge, and its reach is taken as
three widths. The check is on distance alone and is lenient: a location
just outside the data is not flagged. The warning carries `n_far`, `n`,
`reach` and `nearest` (the smallest and largest distance to a nearest
training point).

Each call raises the warning once, with its own counts.
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md),
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
and
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
predict in pieces and raise one warning per call with the totals. A
cross-validation gives the smallest and the largest reach among its
folds, and a sweep counts a location once for each fit that put it out
of reach. A loop of your own over
[`predict()`](https://rdrr.io/r/stats/predict.html) gets one warning per
call; the example below totals them with a handler. To drop the warning
where extrapolation is intended, use
`suppressWarnings(predict(fit, newdata), classes = "spatialkit_gwr_extrapolation")`.

## Bayesian GP boundary

brms 2.17 to 2.22 do not store the Hilbert-space boundary \\L\\ in a
fitted GP basis, and recompute it from whatever rows
[`predict()`](https://rdrr.io/r/stats/predict.html) is handed. This
method therefore appends two rows at the training coordinate extrema
before the posterior draw and drops them from the result, and when
`newdata` reaches past the training range it scales down the `c` of the
`gp()` term by as much as the range grew, so brms rebuilds the boundary
the model was fitted with. A prediction does not depend on which other
rows share the call: chunked, fold-wise and single-call predictions
agree, and
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
does not depend on `chunk_size`. brms 2.23.0 and later store \\L\\ and
reuse it; there `c` is left alone and the two extra rows change nothing.

Predictions outside the training coordinate envelope are extrapolation,
and a line in the log says so. A row further than \\L\\ from the centre
of the training coordinates on either axis is past the edge of the
basis, where the approximate GP is an odd reflection of the fitted
surface, so it is returned as `NA` (a column of `NA` with
`draws = TRUE`) with a warning. The default boundary factor puts that
edge well outside the training data.

## Bayesian GP categorical predictors

A factor or character predictor is handed to brms with the levels the
model was fitted with: the values that occur in the training rows
(`object$info$predictor_levels`). A `newdata` factor that only lists an
extra level, or a model fitted on a subset whose factor kept an unused
level, is predicted as usual. A row whose value is not among the fitted
levels has no prediction: it is returned as `NA`, the other rows are
predicted, and one warning of class `spatialkit_new_levels` names the
predictor, the levels and the number of rows. The result does not depend
on which other rows share the call.

## See also

[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md)
to predict onto a grid, and
[`area_of_applicability()`](https://elkronos.github.io/gis_modeling_toolkit/reference/area_of_applicability.md)
to mark where the predictions apply.

Other fit methods:
[`plot.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.spatial_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`spatial_fit-coef`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md),
[`spatial_fit-fitted`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-fitted.md),
[`spatial_fit-residuals`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-residuals.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)

## Examples

``` r
library(sf)
set.seed(1)
pts <- st_as_sf(
  data.frame(x = 5e5 + runif(60, 0, 1000), y = 5e6 + runif(60, 0, 1000),
             a = rnorm(60)),
  coords = c("x", "y"), crs = 32632
)
pts$z <- 10 + 2 * pts$a + rnorm(60, 0, 0.3)
new <- pts[1:5, ]
new$z <- NULL                          # the response is not needed

rf <- fit_rf_model(pts, "z", "a", num_trees = 50, seed = 1)
predict(rf, newdata = new)
#> [1] 14.185278 10.109833 11.273408  9.993768  8.506111

gw <- fit_gwr_model(pts, "z", "a", bandwidth = 30)
str(predict(gw, newdata = new, se.fit = TRUE))
#> List of 5
#>  $ fit           : num [1:5] 14.62 9.88 11.35 10.07 8.56
#>  $ se.fit        : num [1:5] 0.1884 0.0784 0.0903 0.0768 0.1026
#>  $ se.pred       : num [1:5] 0.374 0.332 0.335 0.332 0.339
#>  $ residual.scale: num 0.323
#>  $ df            : num 46.8

# A point 4 km east of the data is beyond the reach of every kernel.  A
# handler totals the extrapolation warnings of several calls; the log
# lines still print (spatialkit_quiet() turns them off).
far <- st_as_sf(data.frame(x = 5.05e5, y = 5e6 + 500, a = 0),
                coords = c("x", "y"), crs = 32632)
total <- 0
withCallingHandlers(
  for (i in 1:2) predict(gw, newdata = far),
  spatialkit_gwr_extrapolation = function(w) {
    total <<- total + w$n_far
    invokeRestart("muffleWarning")
  })
total
#> [1] 2
```
