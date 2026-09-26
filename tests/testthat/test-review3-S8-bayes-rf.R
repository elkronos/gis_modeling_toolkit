# Regression tests for the third review's Bayesian and random-forest findings
# (slice S8).  The Bayesian tests replace brms::brm() and the posterior
# accessors with recorders, as test-review2-bayes-rf.R does; the one
# real-sampler test at the end is opt-in via SPATIALKIT_TEST_BRMS.

.r3b_rf_pts <- function(n = 80, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a - d$b + rnorm(n, 0, 0.3)
  d
}

.r3b_pts <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z    <- 2 * d$a + rnorm(n)
  d$cat3 <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
                   levels = c("lo", "mid", "hi"))
  d$ord3 <- factor(d$cat3, ordered = TRUE)
  d
}

# Every R warning an expression raises, muffled and returned with its value.
.r3b_warnings <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}


# ---------------------------------------------------------------------------
# S8-BAYES-RF-1: predict.rf_fit() on newdata with no complete row
# ---------------------------------------------------------------------------

test_that("predict.rf_fit() returns all NA when every newdata row is incomplete", {
  skip_if_not_installed("ranger")
  d   <- .r3b_rf_pts()
  fit <- fit_rf_model(d, "z", c("a", "b"), num_trees = 50)
  nd  <- d[1:4, ]
  nd$a <- NA_real_
  # It stopped with ranger's "sample_fraction too small, no observations
  # sampled", from a zero-row frame.
  expect_equal(predict(fit, newdata = nd), rep(NA_real_, 4))
  m <- model_metrics(fit, newdata = nd)
  expect_identical(m$n, 0L)
  # A failure inside ranger is still an error.
  expect_error(predict(fit, d[1:4, ], type = "se"), "keep.inbag")
})

test_that("predict_surface() on an rf_fit survives a chunk with no complete row", {
  skip_if_not_installed("ranger")
  d   <- .r3b_rf_pts()
  fit <- fit_rf_model(d, "z", c("a", "b"), num_trees = 50)
  set.seed(2)
  g <- sf::st_as_sf(
    data.frame(x = 5e5 + seq(0, 950, length.out = 20), y = 5e6 + 500,
               a = c(rep(NA, 10), rnorm(10)), b = rnorm(20)),
    coords = c("x", "y"), crs = 32632)
  out <- predict_surface(fit, grid = g, chunk_size = 10)
  expect_true(all(is.na(out$.pred[1:10])))
  expect_true(all(is.finite(out$.pred[11:20])))
  expect_equal(out$.pred[11:20], predict(fit, newdata = g[11:20, ]))
})


# ---------------------------------------------------------------------------
# S8-BAYES-RF-6: cv_rf() says the out-of-bag gap once, not once per fold
# ---------------------------------------------------------------------------

test_that("cv_rf() warns once about fold forests with no out-of-bag row", {
  skip_if_not_installed("ranger")
  d <- .r3b_rf_pts()
  r <- .r3b_warnings(suppressMessages(
    cv_rf(d, "z", c("a", "b"), k = 4, num_trees = 30, replace = FALSE,
          sample_fraction = 1)))
  # One per fold before, each ending "or score the forest with cv_rf()".
  expect_length(r$warnings, 1L)
  expect_match(r$warnings,
               paste0("cv_rf\\(\\): no training row was out of bag for any ",
                      "tree in 4 of 4 fold forest\\(s\\) \\(replace = FALSE ",
                      "with sample_fraction = 1"))
  expect_match(r$warnings, "cross-validation is unaffected")
  expect_match(r$warnings, "fitted() values or permutation importance. Use",
               fixed = TRUE)
  expect_false(grepl("score the forest with cv_rf", r$warnings, fixed = TRUE))
  # The fold counts rode in fold_metrics and are gone again.
  expect_false(any(startsWith(names(r$value$fold_metrics), "..rf_")))
  expect_identical(r$value$n_folds_succeeded, 4L)
  expect_true(all(is.finite(r$value$predictions$yhat)))

  # Impurity importance does not depend on out-of-bag rows.
  ri <- .r3b_warnings(suppressMessages(
    cv_rf(d, "z", c("a", "b"), k = 4, num_trees = 30, replace = FALSE,
          sample_fraction = 1, importance = "impurity")))
  expect_length(ri$warnings, 1L)
  expect_match(ri$warnings, "has no out-of-bag error or fitted() values. Use",
               fixed = TRUE)
})

test_that("cv_rf() counts fold forests with some rows out of no tree's bag", {
  skip_if_not_installed("ranger")
  d <- .r3b_rf_pts()
  # Three bootstrap trees leave a quarter of the rows in every tree's sample.
  r <- .r3b_warnings(suppressMessages(
    cv_rf(d, "z", c("a", "b"), k = 4, num_trees = 3)))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings,
               "cv_rf\\(\\): in 4 of 4 fold forest\\(s\\) some training rows")
  expect_false(any(startsWith(names(r$value$fold_metrics), "..rf_")))
  # A forest that covers every row says nothing.
  expect_no_warning(suppressMessages(
    cv_rf(d, "z", c("a", "b"), k = 4, num_trees = 50)))
  # fit_rf_model() called directly still warns, with its own advice.
  expect_warning(fit_rf_model(d, "z", c("a", "b"), num_trees = 30,
                              replace = FALSE, sample_fraction = 1),
                 "or score the forest with cv_rf\\(\\)")
})

test_that("cv_rf(parallel = 2) counts the folds the workers fitted", {
  skip_if_not_installed("ranger")
  skip_on_os("windows")
  skip_on_cran()
  d <- .r3b_rf_pts()
  # A warning raised in a forked worker reaches the parent once per distinct
  # text, so the count has to travel with the fold's result.
  r <- .r3b_warnings(suppressMessages(
    cv_rf(d, "z", c("a", "b"), k = 4, num_trees = 30, replace = FALSE,
          sample_fraction = 1, parallel = 2)))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "in 4 of 4 fold forest")
})


# ---------------------------------------------------------------------------
# S8-BAYES-RF-9: print() on a forest with no out-of-bag row
# ---------------------------------------------------------------------------

test_that("print.rf_fit() says the OOB error and importance are undefined", {
  skip_if_not_installed("ranger")
  d  <- .r3b_rf_pts()
  f0 <- suppressWarnings(fit_rf_model(d, "z", c("a", "b"), num_trees = 30,
                                      replace = FALSE, sample_fraction = 1))
  expect_true(is.na(f0$info$oob_rmse))
  expect_true(all(is.nan(f0$info$importance)))
  txt <- paste(utils::capture.output(print(f0)), collapse = "\n")
  # The importance line printed empty and the OOB line was missing.
  expect_match(txt, "OOB RMSE: undefined (no row is out of bag)", fixed = TRUE)
  expect_match(txt, "Importance (permutation): undefined (no row is out of bag)",
               fixed = TRUE)
  # An ordinary forest prints its numbers as before.
  f1  <- fit_rf_model(d, "z", c("a", "b"), num_trees = 50)
  txt <- paste(utils::capture.output(print(f1)), collapse = "\n")
  expect_match(txt, "OOB RMSE: [0-9.]+   OOB R\\^2")
  expect_match(txt, "Importance \\(permutation\\): a=[-0-9.e]+, b=")
})


# ---------------------------------------------------------------------------
# S8-BAYES-RF-2: the standardize_predictors slope prior keeps its dpar
# ---------------------------------------------------------------------------

.r3b_capture_fit <- function(...) {
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$args <- list(...); structure(list(), class = "r3b_stub") },
    .package = "brms")
  fit <- fit_bayesian_spatial_model(..., compute_loo = FALSE,
                                    check_convergence = FALSE)
  list(fit = fit, args = cap$args)
}

.r3b_validates <- function(args) {
  expect_no_error(suppressWarnings(brms::validate_prior(
    args$prior, formula = args$formula, data = args$data,
    family = args$family)))
}

test_that("standardize_predictors' slope prior validates for categorical and mixture fits", {
  skip_if_not_installed("brms")
  d <- .r3b_pts()
  # brms refused both before compiling: "The following priors do not
  # correspond to any model parameter: b ~ normal(0, 5)".
  rc <- suppressMessages(.r3b_capture_fit(
    d, "cat3", "a", family = brms::categorical(), gp_k = 5,
    standardize_predictors = TRUE))
  pr <- as.data.frame(rc$args$prior)
  b  <- pr[pr$class == "b", ]
  expect_setequal(b$dpar, c("mumid", "muhi"))
  expect_true(all(b$prior == "normal(0, 5)"))
  .r3b_validates(rc$args)

  rm <- suppressMessages(.r3b_capture_fit(
    d, "z", "a", family = brms::mixture(stats::gaussian(), stats::gaussian()),
    gp_k = 5, standardize_predictors = TRUE))
  pr <- as.data.frame(rm$args$prior)
  expect_setequal(pr$dpar[pr$class == "b"], c("mu1", "mu2"))
  .r3b_validates(rm$args)

  # A family with one mu keeps the single global row it always had.
  rg <- suppressMessages(.r3b_capture_fit(d, "z", "a", gp_k = 5,
                                          standardize_predictors = TRUE))
  pr <- as.data.frame(rg$args$prior)
  b  <- pr[pr$class == "b", ]
  expect_identical(nrow(b), 1L)
  expect_identical(b$dpar, "")
  expect_identical(b$coef, "")
  .r3b_validates(rg$args)
  ro <- suppressMessages(.r3b_capture_fit(d, "ord3", "a",
                                          family = brms::cumulative(), gp_k = 5,
                                          standardize_predictors = TRUE))
  .r3b_validates(ro$args)
})


# ---------------------------------------------------------------------------
# S8-BAYES-RF-3: cv_bayes() refuses a category family before any MCMC
# ---------------------------------------------------------------------------

test_that("cv_bayes() refuses a categorical or ordinal family before fitting", {
  skip_if_not_installed("brms")
  d <- .r3b_pts()
  calls <- 0L
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(...) {
      calls <<- calls + 1L
      stop("mock fit 4417")
    },
    .package = "spatialkit")
  # Every fold used to compile and sample, and then fail its scoring.
  expect_error(cv_bayes(d, "ord3", "a", k = 2,
                        fit_args = list(family = brms::cumulative())),
               "cv_bayes\\(\\): the .cumulative. family gives a probability per response category")
  expect_error(cv_bayes(d, "cat3", "a", k = 2,
                        fit_args = list(family = brms::categorical())),
               "the .categorical. family")
  expect_error(cv_bayes(d, "ord3", "a", k = 2,
                        fit_args = list(family = brms::acat)),
               "the .acat. family")
  expect_identical(calls, 0L)

  # A numeric family still reaches the fits.
  suppressWarnings(suppressMessages(
    cv_bayes(d, "z", "a", k = 2, fit_args = list(family = stats::poisson()))))
  expect_gt(calls, 0L)
})


# ---------------------------------------------------------------------------
# S8-BAYES-RF-7 and -8: what predict() says for a category family
# ---------------------------------------------------------------------------

.r3b_category_fit <- function(d, family, resp) {
  new_spatial_fit(
    "bayesian_fit",
    engine = structure(list(family = list(family = family)),
                       class = "brmsfit"),
    formula = stats::as.formula(paste(resp, "~ a")), response_var = resp,
    predictor_vars = "a", data_sf = d,
    info = list(coord_scaling = list(x_center = 5e5, x_scale = 300,
                                     y_center = 5e6, y_scale = 300),
                # Enough for .pin_gp_boundary_rows() to append its two rows.
                gp_xy_range = list(x = c(-2, 2), y = c(-2, 2))))
}

test_that("the per-category epred error counts the caller's rows", {
  skip_if_not_installed("brms")
  d   <- .r3b_pts(n = 40)
  fit <- .r3b_category_fit(d, "cumulative", "ord3")
  seen <- integer(0)
  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...) {
      seen <<- c(seen, nrow(newdata))
      array(0.3, dim = c(20L, nrow(newdata), 3L))
    },
    .package = "brms")
  err <- tryCatch(suppressMessages(predict(fit, newdata = d[1:5, ])),
                  error = conditionMessage)
  # brms was handed the two boundary rows as well ...
  expect_identical(seen, 7L)
  # ... which the message counted: "20 x 7 x 3" for five rows.
  expect_match(err, "returned a 20 x 5 x 3 array", fixed = TRUE)
  # The hint offers what works on new rows, not posterior_epred(newdata = ),
  # which the engine refuses without the scaled coordinates.
  expect_match(err, "type = \"predict\", draws = TRUE", fixed = TRUE)
  expect_match(err, "share of draws in each category", fixed = TRUE)
  expect_false(grepl("newdata = )", err, fixed = TRUE))
})

test_that("predict(type = \"predict\") without draws is refused for categorical()", {
  skip_if_not_installed("brms")
  d   <- .r3b_pts(n = 40)
  nd  <- d[1:5, ]
  local_mocked_bindings(
    posterior_predict = function(object, newdata, ...)
      matrix(rep(c(1, 3), length.out = 20L * nrow(newdata)), nrow = 20L),
    .package = "brms")
  cat_fit <- .r3b_category_fit(d, "categorical", "cat3")
  # The mean of unordered category indices came back as a number.
  expect_error(suppressMessages(predict(cat_fit, newdata = nd, type = "predict")),
               "'categorical' family's categories have no order, so the mean")
  expect_error(suppressMessages(predict(cat_fit, newdata = nd, type = "predict",
                                        summary = "median")),
               "so the median of the predicted category indices")
  dr <- suppressMessages(predict(cat_fit, newdata = nd, type = "predict",
                                 draws = TRUE))
  expect_true(is.matrix(dr))
  expect_identical(dim(dr), c(20L, 5L))
  # An ordinal family's mean index is an expected rank, and still returned.
  ord_fit <- .r3b_category_fit(d, "cumulative", "ord3")
  expect_equal(suppressMessages(predict(ord_fit, newdata = nd, type = "predict")),
               rep(2, 5))
})


# ---------------------------------------------------------------------------
# Real sampler (opt-in): a standardised categorical fit samples, and its
# predict() messages count the caller's rows
# ---------------------------------------------------------------------------

test_that("a categorical fit with standardize_predictors = TRUE samples", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  set.seed(20240817)
  n <- 40
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); z <- rnorm(n)
  resp <- 0.004 * x + 1.5 * z + rnorm(n, sd = 0.5)
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, resp = resp),
                      coords = c("x", "y"), crs = 32632)
  pts$cat3 <- cut(pts$resp, stats::quantile(pts$resp, c(0, 1/3, 2/3, 1)),
                  include.lowest = TRUE, labels = c("lo", "mid", "hi"))
  fit <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    fit <- fit_bayesian_spatial_model(pts, "cat3", "z",
                                      family = brms::categorical(),
                                      standardize_predictors = TRUE,
                                      chains = 1, iter = 300, cores = 1,
                                      compute_loo = FALSE, seed = 1234,
                                      gp_k = 6, check_convergence = FALSE))),
    type = "output")
  expect_s3_class(fit, "bayesian_fit")
  expect_false(is.null(fit$info$predictor_scaling$z))
  nd <- pts[1:5, ]
  expect_error(suppressMessages(predict(fit, newdata = nd)),
               "returned a 150 x 5 x 3 array")
  expect_error(suppressMessages(predict(fit, newdata = nd, type = "predict")),
               "categories have no order")
  dr <- suppressMessages(predict(fit, newdata = nd, type = "predict",
                                 draws = TRUE))
  expect_identical(dim(dr), c(150L, 5L))
  expect_true(all(dr %in% 1:3))
})
