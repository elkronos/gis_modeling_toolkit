# tests/testthat/test-fix-response-precision.R
# ---------------------------------------------------------------------------
# BUG-235: a response so far from zero that doubles there are spaced more
# than a millionth of its standard deviation apart (1e15 plus a spread of 1,
# stored to an eighth) was read without a word: make_folds(auto_range =
# TRUE) built other folds and select_features_forward() selected a predictor
# the same response near zero does not select.  prep_model_data(), and so
# every fit and cv_*() function, and estimate_sac_range() now warn with class
# spatialkit_response_imprecise, once per call or run.
# ---------------------------------------------------------------------------

prec_layer <- function(n = 80, seed = 3) {
  set.seed(seed)
  x <- runif(n, 0, 2000); yy <- runif(n, 0, 2000)
  a <- rnorm(n); b <- rnorm(n)
  z <- sin(2 * pi * x / 1000) + 0.5 * a + 0.3 * rnorm(n)
  d <- sf::st_as_sf(data.frame(x = x, yy = yy, z = z, a = a, b = b),
                    coords = c("x", "yy"), crs = 32617)
  d$..row_id <- seq_len(n)
  d
}

shifted <- function(d, off) { d$z <- d$z + off; d }

# The warnings of class spatialkit_response_imprecise a call raises, muffled;
# every other warning passes through.
imprecise_of <- function(expr) {
  got <- list()
  withCallingHandlers(expr, spatialkit_response_imprecise = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  got
}

test_that("prep_model_data() warns once where doubles do not hold a millionth of the spread", {
  d <- prec_layer()
  for (off in c(1e12, 1e15)) {
    w <- imprecise_of(out <- prep_model_data(shifted(d, off), "z", c("a", "b")))
    expect_length(w, 1L)
    msg <- conditionMessage(w[[1L]])
    expect_match(msg, "^prep_model_data\\(\\): `z` lies so far from zero")
    expect_match(msg, "lose their precision")
    u <- spatialkit:::.ulp_of(d$z + off)
    expect_identical(w[[1L]]$spacing, u)
    expect_match(msg, format(signif(u, 3)), fixed = TRUE)
    expect_equal(w[[1L]]$spread, sd(d$z + off))
    # The rows are kept and read as given.
    expect_identical(nrow(out), nrow(d))
    expect_identical(out$z, d$z + off)
  }
  # 0.125 at 1e15, 0.000122 at 1e12.
  expect_identical(spatialkit:::.ulp_of(1e15), 0.125)
})

test_that("a response that doubles hold to a millionth of its spread draws no such warning", {
  d <- prec_layer()
  for (off in c(0, 1e8, -1e8)) {
    w <- imprecise_of(prep_model_data(shifted(d, off), "z", c("a", "b")))
    expect_length(w, 0L)
  }
  # Units change nothing: the ratio of spacing to spread is the same.
  for (s in c(2^-530, 2^530, 1e-160, 1e160)) {
    e <- d; e$z <- e$z * s
    expect_length(imprecise_of(prep_model_data(e, "z", c("a", "b"))), 0L)
  }
  # A response with no spread is left to the tests for a constant.
  e <- d; e$z <- 1e15
  expect_length(imprecise_of(prep_model_data(e, "z", c("a", "b"))), 0L)
  # Nor does one 1e15 apart by a unit in the last place, which the package
  # reads as constant (.constant_column()).
  e$z <- 1e15 + 0.125 * (seq_len(nrow(e)) %% 2)
  expect_length(imprecise_of(prep_model_data(e, "z", c("a", "b"))), 0L)
  # Preparing newdata (no response) says nothing.
  expect_length(imprecise_of(prep_model_data(shifted(d, 1e15), "z", c("a", "b"),
                                             require_response = FALSE)), 0L)
})

test_that("the threshold is a millionth of the spread, read on the largest absolute value", {
  f <- spatialkit:::.warn_response_imprecise
  # sd 1 at 2^32 (spacing 2^-20, about 9.5e-7): no warning; at 2^33
  # (spacing 2^-19, about 1.9e-6): a warning.
  z <- c(-1, 1, -1, 1) / sqrt(4 / 3)       # sd exactly 1
  expect_false(suppressWarnings(f(2^32 + z, "f", "y")))
  expect_true(suppressWarnings(f(2^33 + z, "f", "y")))
  # A negative offset reads the same.
  expect_true(suppressWarnings(f(-2^33 + z, "f", "y")))
  # Fewer than two finite values: nothing.
  expect_false(f(c(1e15, NA, Inf), "f", "y"))
})

test_that("estimate_sac_range() warns once on a response at 1e15, and not at 1e8", {
  skip_if_not_installed("gstat")
  d <- prec_layer(n = 120)
  w <- imprecise_of(suppressMessages(estimate_sac_range(shifted(d, 1e15), "z", seed = 1)))
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "^estimate_sac_range\\(\\): `z` lies so far")
  w <- imprecise_of(suppressMessages(estimate_sac_range(shifted(d, 1e8), "z", seed = 1)))
  expect_length(w, 0L)
  # make_folds(auto_range = TRUE) reads its range there, and says it.
  w <- imprecise_of(suppressMessages(make_folds(shifted(d, 1e15), k = 4, seed = 1,
                                                method = "block_kfold",
                                                auto_range = TRUE, response_var = "z")))
  expect_length(w, 1L)
})

test_that("cv_spatial() and select_features_forward() say it once", {
  d <- prec_layer()
  fit <- function(tr, vars) lm_spatial_fit(tr, "z", vars)
  w <- imprecise_of(suppressMessages(
    cv_spatial(shifted(d, 1e15), "z", c("a", "b"),
               fit_fn = function(tr) fit(tr, c("a", "b")), k = 4, seed = 1)))
  expect_length(w, 1L)
  w <- imprecise_of(select_features_forward(shifted(d, 1e15), "z", c("a", "b"), fit,
                                            k = 4, seed = 1, quiet = TRUE))
  expect_length(w, 1L)
  w <- imprecise_of(select_features_forward(shifted(d, 1e8), "z", c("a", "b"), fit,
                                            k = 4, seed = 1, quiet = TRUE))
  expect_length(w, 0L)
})

test_that("a run says it once for the column, whichever function reads it first", {
  f <- spatialkit:::.warn_response_imprecise
  y <- 1e15 + c(-1, 1, -1, 1)
  end_run <- spatialkit:::.once_per_run(spatialkit:::.prep_warning_classes)
  w <- imprecise_of({
    f(y, "prep_model_data", "z")
    f(y, "estimate_sac_range", "z")
    f(y[-1], "prep_model_data", "z")
  })
  end_run()
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "^prep_model_data\\(\\)")
  # Outside a run each call says it.
  w <- imprecise_of({
    f(y, "prep_model_data", "z")
    f(y, "estimate_sac_range", "z")
  })
  expect_length(w, 2L)
})

test_that("kriging_adequacy() says it once, naming the response and not its own column", {
  skip_if_not_installed("gstat")
  d <- prec_layer(n = 120)
  sq <- lapply(0:3, function(i) {
    a <- (i %% 2) * 1000; b <- (i %/% 2) * 1000
    sf::st_polygon(list(rbind(c(a, b), c(a + 1000, b), c(a + 1000, b + 1000),
                              c(a, b + 1000), c(a, b))))
  })
  cells <- sf::st_sf(poly_id = sprintf("C%d", 1:4), geometry = sf::st_sfc(sq, crs = 32617))
  d$poly_id <- cells$poly_id[as.integer(sf::st_within(d, cells))]
  # The handler of imprecise_of() is the inner one, so it sees its class
  # first; the other warnings (a range refused, say) are not this test's.
  # At 1e12 the internal column is centred (.regression_centre()) and the
  # inner estimate_sac_range() call has nothing to say; at 1e10 it is not
  # centred, doubles there are 1.9e-6 apart, and that call would name the
  # internal column `..z`.
  expect_identical(spatialkit:::.regression_centre(d$z + 1e10, bounded = FALSE), 0)
  for (off in c(1e10, 1e12)) {
    w <- suppressWarnings(imprecise_of(suppressMessages(
      kriging_adequacy(shifted(d, off), "z", cells, k = 3L, nmax = 20L))))
    expect_length(w, 1L)
    expect_match(conditionMessage(w[[1L]]), "^kriging_adequacy\\(\\): `z` lies so far")
    expect_no_match(conditionMessage(w[[1L]]), "..z", fixed = TRUE)
  }
  # A supplied variogram: the response is read all the same.
  e <- shifted(d, 1e15)
  sac <- suppressWarnings(suppressMessages(estimate_sac_range(d, "z", seed = 1)))
  w <- suppressWarnings(imprecise_of(suppressMessages(
    kriging_adequacy(e, "z", cells, sac = sac, k = 3L, nmax = 20L))))
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), "^kriging_adequacy\\(\\): `z` lies so far")
})

test_that("compare_models_cv() says it once when it builds the folds itself", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("gstat")
  d <- prec_layer(n = 120)
  # folds = NULL: make_folds() estimates the range for its diagnostic with
  # every warning muffled.  That call must not count as having said it, or
  # no later preparation of the run does.
  for (auto in c(FALSE, TRUE)) {
    w <- suppressWarnings(imprecise_of(suppressMessages(
      compare_models_cv(shifted(d, 1e15), "z", c("a", "b"), models = "RF",
                        k = 4, seed = 1, auto_range = auto))))
    expect_length(w, 1L)
  }
  f <- make_folds(d, k = 4, seed = 1, method = "block_kfold")
  w <- suppressWarnings(imprecise_of(suppressMessages(
    compare_models_cv(shifted(d, 1e15), "z", c("a", "b"), models = "RF",
                      folds = f, k = 4, seed = 1))))
  expect_length(w, 1L)
})

test_that("cv_spatial() says it once when its fit_fn prepares each fold again", {
  skip_if_not_installed("ranger")
  d <- prec_layer()
  w <- suppressWarnings(imprecise_of(suppressMessages(
    cv_spatial(shifted(d, 1e15), "z", c("a", "b"),
               fit_fn = function(tr) fit_rf_model(tr, "z", c("a", "b"), num_trees = 50),
               k = 4, seed = 1))))
  expect_length(w, 1L)
  # The run closes with the call: a second call says it again.
  w <- suppressWarnings(imprecise_of(suppressMessages(
    prep_model_data(shifted(d, 1e15), "z", c("a", "b")))))
  expect_length(w, 1L)
})
