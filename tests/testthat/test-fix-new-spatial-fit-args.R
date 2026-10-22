# tests/testthat/test-fix-new-spatial-fit-args.R
# ---------------------------------------------------------------------------
# new_spatial_fit() refuses a fit that no method can use, and names the
# argument.
#
# Every method on a spatial_fit reads `data_sf`, `formula`, `response_var`
# and `predictor_vars`.  The constructor used to build a fit from any of
# them, and the error, when there was one, came later from another function
# and named something else: a list as `data_sf` gave n = NULL and "missing
# value where TRUE/FALSE needed" from summary(); a misspelt response gave
# "response 'Price' is NULL, not numeric"; `predictor_vars = NULL` gave
# area_of_applicability()'s "supply `model`, or both `train_sf` and
# `predictor_vars`" to a caller who had supplied `model`.  A column number as
# `response_var` gave no error at all: model_metrics() scored the fitted
# values against the first column of `data_sf`, and a factor naming the
# response did the same through its code.  The geometry column's name passed
# as a response or a predictor, and failed later: the methods read both from
# st_drop_geometry(data_sf).
#
# What the constructor still accepts is pinned too: no predictors, a factor
# or character response (a Bayesian fit can be ordinal; summary() refuses a
# non-numeric response itself), any engine, `info = NULL`, a zero-row layer.
# ---------------------------------------------------------------------------

.nsf_site <- function(n = 40) {
  set.seed(1)
  s <- sf::st_as_sf(data.frame(x = stats::runif(n, 0, 1000),
                               y = stats::runif(n, 0, 1000),
                               elev = stats::rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  s$price <- 10 + 2 * s$elev + stats::rnorm(n)
  s
}

# The ?new_spatial_fit example's arguments, one replaced at a time.  A whole
# argument is replaced: modifyList() would merge a data.frame into the sf.
.nsf_build <- function(site, ...) {
  a <- list(subclass = "nsf_fit",
            engine = stats::lm(price ~ elev, sf::st_drop_geometry(site)),
            formula = price ~ elev, response_var = "price",
            predictor_vars = "elev", data_sf = site)
  x <- list(...)
  for (k in names(x)) a[k] <- list(x[[k]])
  do.call(new_spatial_fit, a)
}

test_that("new_spatial_fit() refuses a data_sf that is not an sf object", {
  site <- .nsf_site()
  expect_error(.nsf_build(site, data_sf = sf::st_drop_geometry(site)),
               paste0("new_spatial_fit(): `data_sf` must be an sf object (the ",
                      "training rows with their geometry); got data.frame."),
               fixed = TRUE)
  # A list has no nrow(): the fit's n was NULL.
  expect_error(.nsf_build(site, data_sf = as.list(sf::st_drop_geometry(site))),
               "`data_sf` must be an sf object (the training rows with their geometry); got list.",
               fixed = TRUE)
  expect_error(.nsf_build(site, data_sf = NULL), "got NULL.", fixed = TRUE)
})

test_that("new_spatial_fit() refuses a formula that is not a formula", {
  site <- .nsf_site()
  expect_error(.nsf_build(site, formula = "price ~ elev"),
               "new_spatial_fit(): `formula` must be a formula; got character.",
               fixed = TRUE)
  expect_error(.nsf_build(site, formula = NULL),
               "`formula` must be a formula; got NULL.", fixed = TRUE)
})

test_that("new_spatial_fit() refuses a response_var that is not one column of data_sf", {
  site <- .nsf_site()
  # names(site) is elev, geometry, price: a column number was read as
  # data_sf[[1]], and the fitted prices were scored against elev.
  expect_error(.nsf_build(site, response_var = 1L),
               paste0("new_spatial_fit(): `response_var` must name one attribute ",
                      "column of `data_sf`; got integer of length 1."), fixed = TRUE)
  # A factor indexes by its code, as a number does: factor("price") is code 1,
  # so data_sf[[factor("price")]] is elev, although factor("price") %in%
  # names(data_sf) is TRUE.
  expect_identical(site[[factor("price")]], site$elev)
  expect_error(.nsf_build(site, response_var = factor("price")),
               "`response_var` must name one attribute column of `data_sf`; got factor of length 1.",
               fixed = TRUE)
  expect_error(.nsf_build(site, response_var = "Price"),
               "`response_var` must name one attribute column of `data_sf`; 'Price' is not one.",
               fixed = TRUE)
  expect_error(.nsf_build(site, response_var = c("price", "elev")),
               "`response_var` must name one attribute column of `data_sf`; got character of length 2.",
               fixed = TRUE)
  expect_error(.nsf_build(site, response_var = NA_character_),
               "'NA' is not one.", fixed = TRUE)
  expect_error(.nsf_build(site, response_var = NULL),
               "`response_var` must name one attribute column of `data_sf`; got NULL.",
               fixed = TRUE)
})

test_that("new_spatial_fit() refuses the geometry column as the response or a predictor", {
  site <- .nsf_site()
  # The methods read both from st_drop_geometry(data_sf), where the geometry
  # is not a column: model_metrics() said "response 'geometry' is NULL, not
  # numeric ... Convert the column first", and area_of_applicability() that
  # the predictor was absent from the training data.
  expect_error(.nsf_build(site, response_var = "geometry"),
               paste0("new_spatial_fit(): `response_var` must name one attribute ",
                      "column of `data_sf`; 'geometry' is its geometry column."),
               fixed = TRUE)
  expect_error(.nsf_build(site, predictor_vars = "geometry"),
               paste0("new_spatial_fit(): `predictor_vars` names 'geometry', the ",
                      "geometry column of `data_sf`; predictors are attribute columns."),
               fixed = TRUE)
  expect_error(.nsf_build(site, predictor_vars = c("elev", "geometry")),
               "`predictor_vars` names 'geometry', the geometry column", fixed = TRUE)
  # The geometry column is the one sf names, whatever it is called; a column
  # called "geometry" beside a geometry called "geom" is an attribute.
  g <- sf::st_set_geometry(site, NULL)
  g$geometry <- g$elev
  g <- sf::st_sf(g, geom = sf::st_geometry(site))
  expect_identical(attr(g, "sf_column"), "geom")
  expect_error(.nsf_build(g, response_var = "geom"),
               "'geom' is its geometry column.", fixed = TRUE)
  expect_error(.nsf_build(g, predictor_vars = "geom"),
               "`predictor_vars` names 'geom', the geometry column of `data_sf`",
               fixed = TRUE)
  f <- .nsf_build(g, formula = price ~ geometry, predictor_vars = "geometry")
  expect_identical(f$predictor_vars, "geometry")
  expect_identical(.nsf_build(g, response_var = "geometry")$response_var, "geometry")
})

test_that("new_spatial_fit() refuses predictor_vars that are not columns of data_sf", {
  site <- .nsf_site()
  expect_error(.nsf_build(site, predictor_vars = "Elev"),
               paste0("new_spatial_fit(): `predictor_vars` names a column that ",
                      "`data_sf` does not have: 'Elev'."), fixed = TRUE)
  expect_error(.nsf_build(site, predictor_vars = c("Elev", "elev", "slope")),
               "`predictor_vars` names columns that `data_sf` does not have: 'Elev', 'slope'.",
               fixed = TRUE)
  expect_error(.nsf_build(site, predictor_vars = NULL),
               paste0("new_spatial_fit(): `predictor_vars` must be a character ",
                      "vector of column names (character(0) for a model without ",
                      "predictors); got NULL."), fixed = TRUE)
  expect_error(.nsf_build(site, predictor_vars = factor("elev")),
               "must be a character vector of column names (character(0) for a model without predictors); got factor.",
               fixed = TRUE)
})

test_that("new_spatial_fit() names `subclass` when it is not one class name, and what it got", {
  site <- .nsf_site()
  got <- list("character of length 2" = c("a", "b"), "NA" = NA_character_,
              "an empty string" = "", "integer of length 1" = 1L,
              "NULL" = NULL)
  for (k in names(got))
    expect_error(.nsf_build(site, subclass = got[[k]]),
                 paste0("new_spatial_fit(): `subclass` must be one class name, a ",
                        "character string such as \"lm_fit\"; got ", k, "."),
                 fixed = TRUE)
})

test_that("new_spatial_fit() refuses an `info` that is not a list, and names it", {
  site <- .nsf_site()
  # A data.frame stopped with "replacement has 0 rows, data has 2", and a
  # character vector was coerced to a list with a warning.
  expect_error(.nsf_build(site, info = data.frame(a = 1:2)),
               paste0("new_spatial_fit(): `info` must be a list of model-specific ",
                      "extras (list() for none); got data.frame."), fixed = TRUE)
  expect_error(.nsf_build(site, info = "oob"),
               "`info` must be a list of model-specific extras (list() for none); got character.",
               fixed = TRUE)
  expect_error(.nsf_build(site, info = new.env()), "`info` must be a list", fixed = TRUE)
  # NULL is no extras, as it was.
  f <- .nsf_build(site, info = NULL)
  expect_identical(names(f$info), ".cache")
  expect_true(is.environment(f$info$.cache))
})

test_that("new_spatial_fit()'s argument errors carry no call", {
  site <- .nsf_site()
  # The errors are raised with call. = FALSE.  Through do.call(), as here and
  # in a fit_fn built that way, a call would head the message with the whole
  # engine and data_sf deparsed.
  bad <- list(list(subclass = c("a", "b")),
              list(data_sf = sf::st_drop_geometry(site)),
              list(formula = "price ~ elev"),
              list(response_var = "Price"),
              list(response_var = 1L),
              list(response_var = "geometry"),
              list(predictor_vars = NULL),
              list(predictor_vars = "Elev"),
              list(predictor_vars = "geometry"),
              list(info = "oob"))
  for (b in bad) {
    err <- tryCatch(do.call(.nsf_build, c(list(site), b)), error = identity)
    expect_s3_class(err, "error")
    expect_null(conditionCall(err))
  }
})

test_that("new_spatial_fit() still builds what its methods can use", {
  site <- .nsf_site()
  # No predictors: the intercept-only model select_features_forward() and the
  # GWR's local mean fit.
  f <- .nsf_build(site, formula = price ~ 1, predictor_vars = character(0))
  expect_identical(f$predictor_vars, character(0))
  expect_identical(f$n, 40L)
  # A factor or ordinal response (a Bayesian fit can have one), and a
  # character one, which summary() and model_metrics() refuse themselves.
  s2 <- site
  s2$ord <- factor(cut(s2$price, 3), ordered = TRUE)
  s2$price_chr <- as.character(round(s2$price, 2))
  expect_s3_class(.nsf_build(s2, response_var = "ord"), "spatial_fit")
  chr <- .nsf_build(s2, response_var = "price_chr")
  registerS3method("fitted", "nsf_fit",
                   function(object, ...) as.numeric(stats::fitted(object$engine)))
  expect_error(model_metrics(chr), "response 'price_chr' is character, not numeric")
  # Any engine, and an `info` that carries another fit's cache.
  expect_s3_class(.nsf_build(site, engine = NULL, info = list(.cache = new.env())),
                  "spatial_fit")
  # A layer with no rows: the checks are on the arguments' form and names,
  # not on the rows, and the fit's n says it holds none.
  expect_identical(.nsf_build(site, data_sf = site[0, ])$n, 0L)
})

test_that("a fit_fn whose fit names a column the training rows lack fails each fold with the reason", {
  site <- .nsf_site()
  claims_slope <- function(train_sf, ...)
    new_spatial_fit("nsf_fit", engine = NULL, formula = price ~ elev,
                    response_var = "price", predictor_vars = "slope",
                    data_sf = train_sf)
  cv <- suppressWarnings(suppressMessages(
    cv_spatial(site, "price", "elev", fit_fn = claims_slope, k = 3, seed = 1)))
  expect_equal(cv$n_folds_succeeded, 0)
  expect_true(all(grepl("`predictor_vars` names a column that `data_sf` does not have: 'slope'.",
                        cv$fold_status$message, fixed = TRUE)))
})

test_that("select_features_forward() takes the fold mean as the baseline when the fit for no predictors is refused", {
  site <- .nsf_site(60)
  site$noise <- stats::rnorm(60)
  registerS3method("predict", "nsf_sel_fit", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))
  })
  # The same learner twice: one gives character(0) for the empty set, the
  # other NULL, which the constructor refuses.  The null-model probe then
  # fails on every fold, and the sweep says so through baseline$type.
  learner <- function(empty) function(train_sf, vars) {
    f <- stats::reformulate(if (length(vars)) vars else "1", "price")
    new_spatial_fit("nsf_sel_fit", stats::lm(f, sf::st_drop_geometry(train_sf)), f,
                    "price", predictor_vars = if (length(vars)) vars else empty,
                    data_sf = train_sf)
  }
  run <- function(empty) suppressWarnings(suppressMessages(
    select_features_forward(site, "price", c("elev", "noise"), learner(empty),
                            k = 3, method = "random_kfold", seed = 1)))
  expect_identical(run(character(0))$baseline$type, "null model")
  expect_identical(run(NULL)$baseline$type, "fold mean")
})
