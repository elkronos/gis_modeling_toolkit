# Predict from a GWR spatial model

When `newdata` is NULL, returns the in-sample fitted values. Otherwise
solves one local regression at each new location, with the fit's kernel
and bandwidth (an adaptive bandwidth counts neighbours among the
training points), and returns \\x^\top\hat\beta(u)\\. The local
regressions are weighted least squares solved by a pivoted QR, the solve
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
uses; on a window of full rank the values equal those of
[`GWmodel::gwr.predict()`](https://rdrr.io/pkg/GWmodel/man/gwr.predict.html)
to rounding error, and at a training location with that row's own
predictors they equal
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html). Each location
stands alone: one that cannot be estimated does not affect the others.
`newdata` is first transformed to the CRS used during fitting (via
[`ensure_projected()`](https://elkronos.github.io/gis_modeling_toolkit/reference/ensure_projected.md)),
so predictions are computed in a single coordinate system regardless of
the CRS newdata arrives in.

## Usage

``` r
# S3 method for class 'gwr_fit'
predict(object, newdata = NULL, se.fit = FALSE, ...)
```

## Arguments

- object:

  A `gwr_fit` object.

- newdata:

  An sf object with the same predictors. The response variable need not
  be present (true out-of-sample prediction is supported). NULL = fitted
  values.

- se.fit:

  Logical; also return standard errors. Default `FALSE`.

- ...:

  Ignored.

## Value

With `se.fit = FALSE`, a numeric vector aligned to `nrow(newdata)`, with
`NA` for rows dropped as missing or non-finite, and for locations whose
local regression cannot be estimated (see **When a prediction is NA**);
a warning counts those. If the design matrix for `newdata` cannot be
built, every value is `NA` and a warning says why. CRS-less `newdata`
first receives the interpretation the training data got, so the same
rows land where they did at fit time.

With `se.fit = TRUE`, a list: `fit` (that vector), `se.fit` (the
standard error of the fitted value, \\\hat\sigma\sqrt{s^\top s}\\ with
\\s\\ the weights the local regression puts on the training responses),
`se.pred` (the standard error of a new observation there,
\\\hat\sigma\sqrt{1 + s^\top s}\\, the square root of
[`GWmodel::gwr.predict()`](https://rdrr.io/pkg/GWmodel/man/gwr.predict.html)'s
prediction variance), `residual.scale` (\\\hat\sigma\\, from the fit's
residual sum of squares over its effective degrees of freedom) and `df`
(those degrees of freedom). The standard errors treat the bandwidth as
known and the local model as true at each location; they do not measure
the bias of a local fit and are `NA` where the prediction is.

## When a prediction is NA

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
decision uses the rule and the tolerance of the fit (see
[`fit_gwr_model`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
**The local solve**). A warning of class
`"spatialkit_gwr_not_estimable"` counts the `NA` predictions.

## Locations beyond the reach of the fit

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
just outside the data is not flagged.

The warning is raised once per call, with that call's counts.
[`predict_surface()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict_surface.md),
which predicts a grid in chunks, and
[`cv_gwr()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_gwr.md)
and
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
which predict fold by fold, collect the warnings of their pieces and
raise one under their own name, with the totals and the number of calls
they are over (`calls`). The folds of a cross-validation are separate
fits, so a cross-validation gives the smallest and the largest reach
among them.
[`cv_block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_block_size_sweep.md)
and
[`select_features_forward()`](https://elkronos.github.io/gis_modeling_toolkit/reference/select_features_forward.md)
run many cross-validations of one learner and raise one warning for the
sweep. A sweep predicts the same rows from several fits, so its totals
count a location once for each fit that put it out of reach, and the
message says so. `model_metrics(fit, newdata)` is one call and raises
the warning once. A loop of your own over
[`predict()`](https://rdrr.io/r/stats/predict.html) gets one warning per
call; the condition carries its numbers as fields (`n_far`, `n`,
`reach`, and `nearest`, the smallest and largest distance to a nearest
training point), so a handler can total them:
`withCallingHandlers(..., spatialkit_gwr_extrapolation = function(w) { total <<- total + w$n_far; invokeRestart("muffleWarning") })`.
To drop it where extrapolation is intended, use
`suppressWarnings(predict(fit, newdata), classes = "spatialkit_gwr_extrapolation")`.
The warning that counts `NA` predictions (class
`"spatialkit_gwr_not_estimable"`) carries `n_na`, `n` and `n_empty` and
is collected the same way.

## See also

Other methods on a fitted model:
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md),
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md),
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md),
[`fitted.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.bayesian_fit.md),
[`fitted.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.gwr_fit.md),
[`fitted.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fitted.rf_fit.md),
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md),
[`predict.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.bayesian_fit.md),
[`predict.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/predict.rf_fit.md),
[`print.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.rf_fit.md),
[`print.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/print.spatial_fit.md),
[`residuals.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.bayesian_fit.md),
[`residuals.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.gwr_fit.md),
[`residuals.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/residuals.rf_fit.md),
[`summary.spatial_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/summary.spatial_fit.md)
