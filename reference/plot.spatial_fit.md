# Plot a fitted spatial model

Diagnostic plots for a `spatial_fit`. The package previously shipped
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods but no
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), so the checks
most likely to reveal a problem (is there structure left in the
residuals, and where is it) had to be written by hand each time.

## Usage

``` r
# S3 method for class 'spatial_fit'
plot(
  x,
  type = c("residuals", "observed_predicted", "variogram", "coefficients"),
  response = TRUE,
  term = NULL,
  mask = TRUE,
  ...
)
```

## Arguments

- x:

  A `spatial_fit`.

- type:

  One of:

  `"residuals"`

  :   Residuals mapped at the training locations. Spatial structure here
      is the signal that the model has not captured the autocorrelation.

  `"observed_predicted"`

  :   Observed against fitted, with a 1:1 reference line.

  `"variogram"`

  :   Empirical variogram of the residuals with the fitted model
      overlaid, so the fit can be judged rather than trusted, and
      (unless `response = FALSE`) the variogram of the response itself
      on the same points and lags, drawn hollow with a dashed fit. The
      gap between the two curves is the spatial structure the model
      absorbed: a residual sill well below the response sill means most
      of it, two curves that coincide mean none. When both models were
      fitted the caption gives the residual sill as a share of the
      response sill, and the two effective ranges; the residual range is
      expected to come out shorter and the residual sill lower even when
      the model is right, because residuals of a fitted trend understate
      the variogram (see
      [`estimate_sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/estimate_sac_range.md),
      "Detrending and the residual-variogram bias"). The distance axis
      is labelled in the units of the CRS the variogram was actually
      fitted in, which is not necessarily the fit's own CRS (lon/lat
      data are projected first). A single-direction fit names its
      azimuth in the subtitle; a fit that did not converge says so in
      the caption, since the overlaid model line is then not a fit to
      believe. Requires 'gstat'.

  `"coefficients"`

  :   For a GWR fit only: the local coefficient of one `term` mapped at
      the training locations, which is the reason to fit GWR at all.
      Locations where the local design is collinear (the kernel-weighted
      window's scaled condition index is above 30, or the window is
      singular) are drawn hollow and grey (`mask = TRUE`), because the
      smooth surface a naive map draws over them is the picture of an
      unstable estimate, not of a relationship; the subtitle counts
      them. The condition indices are the fit's
      `info$local_collinearity`, computed for every location when the
      model was fitted. A diverging scale centred on zero is used when
      the coefficient changes sign, otherwise a sequential one.

- response:

  Logical, default `TRUE`: for `type = "variogram"`, overlay the
  response's own variogram. Ignored by the other types.

- term:

  For `type = "coefficients"`: which local coefficient to map, one of
  the names `coef(x)` returns. Default `NULL`: the first predictor.
  Ignored by the other types.

- mask:

  For `type = "coefficients"`: whether to draw locations whose local
  design is collinear (scaled condition index of the kernel-weighted
  window above 30, or singular) as hollow grey points instead of
  colouring them by a coefficient that is not to be believed there.
  Default `TRUE`. Locations whose coefficient is non-finite are masked
  either way.

- ...:

  Ignored.

## Value

A `ggplot` object.

## See also

[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md)
for the coefficients themselves.

Other plotting:
[`plot.aoa()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.aoa.md),
[`plot.block_size_sweep()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.block_size_sweep.md),
[`plot.feature_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.feature_selection.md),
[`plot.gwr_model_selection()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.gwr_model_selection.md),
[`plot.resolution_profile()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.resolution_profile.md),
[`plot.sac_range()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot.sac_range.md),
[`plot_calibration()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_calibration.md),
[`plot_cv_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_cv_metrics.md),
[`plot_folds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/plot_folds.md)

## Examples

``` r
# Works on any spatial_fit; a forest keeps the example free of the optional
# GWR/Stan backends.
if (requireNamespace("ranger", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(sf)
  set.seed(1)
  n <- 120
  pts <- st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  pts$price <- 10 + 0.01 * st_coordinates(pts)[, 1] + 2 * pts$elev + rnorm(n)
  fit <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
  plot(fit, type = "residuals")
  plot(fit, type = "observed_predicted")
  if (requireNamespace("gstat", quietly = TRUE))
    plot(fit, type = "variogram")
}
```
