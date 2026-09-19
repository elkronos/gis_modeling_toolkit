# Build a spatial_fit S3 object

The constructor for the `spatial_fit` class, and the public entry point
for plugging your own model backend into this package. The three
built-in fitters
([`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md)
and
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md))
all end by calling it, and so should a custom `fit_fn` written for
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md):
wrapping your model in a `spatial_fit` is what lets it use the package's
folds, metrics, comparison and area-of-applicability machinery
unchanged.

## Usage

``` r
new_spatial_fit(
  subclass,
  engine,
  formula,
  response_var,
  predictor_vars,
  data_sf,
  info = list()
)
```

## Arguments

- subclass:

  Character scalar naming the class to stamp on the object: one of the
  built-ins `"gwr_fit"`, `"bayesian_fit"` or `"rf_fit"`, or any name of
  your own for a custom backend (say `"lm_fit"`). It is the S3 dispatch
  key: [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`coef()`](https://rdrr.io/r/stats/coef.html) on the result all
  dispatch on it, so a custom `subclass` **requires** a matching
  `predict.<subclass>()` method to be usable with
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).

- engine:

  The raw model object your backend produced (an `lm`, a `ranger`
  object, a `brmsfit`, ...). Nothing inspects it except your own
  methods.

- formula:

  A formula.

- response_var:

  Character(1).

- predictor_vars:

  Character vector.

- data_sf:

  An sf object used for fitting.

- info:

  Named list of model-specific extras. Set `fitted_are_oob = TRUE` when
  your [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) values
  are held out instead of in-sample, so
  [`summary()`](https://rdrr.io/r/base/summary.html) labels them
  correctly.

## Value

An object of class `c(subclass, "spatial_fit")`.

## Details

There are two obligations. Return an object built here from your
`fit_fn`, and define a
[`predict()`](https://rdrr.io/r/stats/predict.html) method for the
`subclass` you chose.
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
scores folds by calling the
[`predict()`](https://rdrr.io/r/stats/predict.html) generic on the fit,
so without a matching `predict.<subclass>()` every fold fails. A
`fitted.<subclass>()` method returning one value per row of the fit's
`data_sf`, in the same order, is **required** by
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md):
both error, naming the method to define, without one.
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html) methods are optional.

## The coef() contract

[`coef()`](https://rdrr.io/r/stats/coef.html) on one of the three
built-in backends either returns the coefficients or signals an error.
It never returns `NULL`. A custom subclass inherits
[`stats::coef.default()`](https://rdrr.io/r/stats/coef.html), which
returns `NULL`, so define a `coef.<subclass>()` that errors when your
backend has no coefficients; otherwise the hazard described below
applies to your own fits.
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.rf_fit.md)
always errors, because a forest has no coefficients;
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.gwr_fit.md)
and
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/coef.bayesian_fit.md)
error when the backend cannot supply them (a missing package, an engine
without the expected component). A `NULL` return would be
indistinguishable from "this model has no fixed effects", so
`lapply(fits, coef)` would quietly produce a shorter answer than the
caller expected. Wrap in [`try()`](https://rdrr.io/r/base/try.html) or
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) when sweeping
over a heterogeneous list of fits.

## See also

[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
which consumes a custom `fit_fn`;
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
for a worked built-in fitter.

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`prep_model_data()`](https://elkronos.github.io/gis_modeling_toolkit/reference/prep_model_data.md)

## Examples

``` r
library(sf)
set.seed(1)
n <- 80
site <- st_as_sf(
  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
             elev = rnorm(n)),
  coords = c("x", "y"), crs = 32632
)
site$price <- 10 + 0.01 * (st_coordinates(site)[, 1] - 5e5) +
  2 * site$elev + rnorm(n)

# A custom backend: an ordinary linear model behind the spatial_fit interface.
lm_fit <- function(train_sf) {
  new_spatial_fit(
    subclass       = "lm_fit",
    engine         = lm(price ~ elev, st_drop_geometry(train_sf)),
    formula        = price ~ elev,
    response_var   = "price",
    predictor_vars = "elev",
    data_sf        = train_sf
  )
}

# Required: cv_spatial() scores each fold through the predict() generic,
# which dispatches on the subclass named above.
predict.lm_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data_sf
  as.numeric(stats::predict(object$engine, st_drop_geometry(newdata)))
}
registerS3method("predict", "lm_fit", predict.lm_fit)

cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1)
#> cv_spatial(): no folds supplied -- using spatial block k-fold CV (k=3).
cv$overall
#>       RMSE      MAE     MAPE    SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 3.472328 3.009688 21.29953 20.27029 0.2341257     NA     80     80      80
# Always check these two agree before trusting the metrics above.
c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
#> attempted succeeded 
#>         3         3 
```
