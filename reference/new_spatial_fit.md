# Build a spatial_fit object for a backend of your own

The constructor for the `spatial_fit` class, and the way to plug a model
backend of your own into this package. The three built-in fitters
([`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md)
and
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md))
all end by calling it, and so should a custom `fit_fn` written for
[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md).
A model wrapped in a `spatial_fit` works with the package's folds,
metrics, comparisons and area of applicability unchanged.

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
  dispatch on it, so a custom `subclass` needs a matching
  `predict.<subclass>()` method before
  [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
  can use it.

- engine:

  The raw model object your backend produced (an `lm`, a `ranger`
  object, a `brmsfit`, ...). Nothing inspects it except your own
  methods.

- formula:

  A formula: the model, as
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) show it. Its left
  side, when it has one, holds `response_var` (`log(price) ~ elev` for
  `"price"`). Its right side uses the columns `predictor_vars` names and
  no other attribute column of `data_sf` than these and the response,
  read by [`all.vars()`](https://rdrr.io/r/base/allnames.html), so
  `log(elev)` and `elev * slope` count as their columns, and so do the
  columns of an [`offset()`](https://rdrr.io/r/stats/offset.html) or a
  grouping term; names that are not columns of `data_sf` and the
  reserved names `..x` and `..y` are not checked, and a right side with
  `.` is not read. A column a `fit_fn` derives itself and keeps out of
  `predictor_vars`, such as coordinates `x` and `y` for a trend surface,
  needs a name that is not a column of `data_sf`: on a layer that keeps
  its coordinates as columns `x` and `y`, derive them as `cx` and `cy`,
  or drop those columns first.

- response_var:

  Character(1): the name of the response column of `data_sf`, an
  attribute column (not the geometry).

- predictor_vars:

  Character vector of attribute column names of `data_sf`, not including
  `response_var`: the columns the right side of `formula` uses.
  `character(0)` for a model without predictors.

- data_sf:

  An sf object with its geometry column: the rows the model was fitted
  on, with their geometry, in the order
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns its
  values.

- info:

  Named list of model-specific extras
  ([`list()`](https://rdrr.io/r/base/list.html) for none). Set
  `fitted_are_oob = TRUE` when your
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) values are
  held out instead of in-sample, so
  [`summary()`](https://rdrr.io/r/base/summary.html) labels them
  correctly.

## Value

An object of class `c(subclass, "spatial_fit")`: a list with `engine`,
`formula`, `response_var`, `predictor_vars` and `data_sf` as given, `n`
(the number of rows of `data_sf`) and `info` (the list given, with a
private cache environment added as `info$.cache`, which methods should
not modify).

## Details

A custom backend needs three things:

1.  a `fit_fn` that returns an object built by `new_spatial_fit()`;

2.  a `predict.<subclass>()` method for the `subclass` you chose.
    [`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md)
    scores each fold by calling the
    [`predict()`](https://rdrr.io/r/stats/predict.html) generic on the
    fit, so without one every fold fails;

3.  a `fitted.<subclass>()` method that returns one value per row of the
    fit's `data_sf`, in the same order.
    [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
    need it, and without one they stop with an error that names the
    method to define.

[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html) methods are optional.

The constructor checks the fields the package's functions read from a
fit, and stops with an error that names the argument: `subclass` must be
one class name, `data_sf` an sf object with its geometry column,
`formula` a formula whose left side, when it has one, holds the
response, `response_var` the name of one attribute column of `data_sf`
(a column other than its geometry), `predictor_vars` a character vector
of its attribute columns other than the response, the columns the right
side of `formula` uses, and `info` a list. It does not check the
response's type (a Bayesian fit can have a factor or ordinal response);
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`model_metrics()`](https://elkronos.github.io/gis_modeling_toolkit/reference/model_metrics.md)
refuse a response that is not numeric or logical.

## The coef() contract

[`coef()`](https://rdrr.io/r/stats/coef.html) on one of the three
built-in backends either returns the coefficients or signals an error.
It never returns `NULL`. A custom subclass inherits
[`stats::coef.default()`](https://rdrr.io/r/stats/coef.html), which
returns `NULL`, so define a `coef.<subclass>()` that errors when your
backend has no coefficients.
[`coef.rf_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)
always errors, because a forest has no coefficients;
[`coef.gwr_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)
and
[`coef.bayesian_fit()`](https://elkronos.github.io/gis_modeling_toolkit/reference/spatial_fit-coef.md)
error when the backend cannot supply them (a missing package, an engine
without the expected component). A `NULL` return could not be told from
"this model has no fixed effects", and `lapply(fits, coef)` would return
a shorter answer than the caller expected. Wrap the call in
[`try()`](https://rdrr.io/r/base/try.html) or
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) when sweeping
over a list of fits from different backends.

## See also

[`cv_spatial()`](https://elkronos.github.io/gis_modeling_toolkit/reference/cv_spatial.md),
which takes a custom `fit_fn`;
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md)
for a built-in fitter.

Other model fitting:
[`fit_bayesian_spatial_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_bayesian_spatial_model.md),
[`fit_gwr_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_gwr_model.md),
[`fit_rf_model()`](https://elkronos.github.io/gis_modeling_toolkit/reference/fit_rf_model.md),
[`gp_lengthscale_bounds()`](https://elkronos.github.io/gis_modeling_toolkit/reference/gp_lengthscale_bounds.md)

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

# Required by summary() and model_metrics(): one value per training row.
fitted.lm_fit <- function(object, ...) as.numeric(stats::fitted(object$engine))
registerS3method("fitted", "lm_fit", fitted.lm_fit)

summary(lm_fit(site))
#> Summary of <lm_fit> fit (n = 80)
#> 
#>   Formula: price ~ elev
#> 
#>   In-sample metrics:
#>     RMSE    = 2.8624
#>     MAE     = 2.4374
#>     R^2     = 0.2622
#>     SMAPE   = 16.40%
cv <- cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1)
#> cv_spatial(): no folds supplied -- using spatial block k-fold CV (k=3).
cv$overall
#>       RMSE      MAE     MAPE    SMAPE        R2 Adj_R2 n_pred n_MAPE n_SMAPE
#> 1 3.472328 3.009688 21.29953 20.27029 0.2341257     NA     80     80      80
# Check that these two agree before reading the metrics above.
c(attempted = cv$n_folds_attempted, succeeded = cv$n_folds_succeeded)
#> attempted succeeded 
#>         3         3 
```
