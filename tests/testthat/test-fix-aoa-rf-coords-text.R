# tests/testthat/test-fix-aoa-rf-coords-text.R
# ---------------------------------------------------------------------------
# What fit_rf_model() says about coordinates.
#
# The logged caution and the help said a forest given the coordinates will
# "fail wherever it has not been".  Across 21 datasets the forest with
# coordinates had the lower RMSE in 601 of 630 block-fold comparisons and
# was never judged worse, so that is not a general statement.  What holds is
# that the out-of-bag error and random folds overstate such a forest's
# accuracy and that the gain fades with distance, so the advice to validate
# with blocked or distance-matched folds stays.  A fit without coordinates
# prints that it uses no spatial information.
# ---------------------------------------------------------------------------

txt_pts <- function(n = 120, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000),
               a = stats::rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + stats::rnorm(n, 0, 0.3)
  d
}

test_that("the include_coords caution advises blocked folds without predicting failure", {
  skip_if_not_installed("ranger")
  rm(list = ls(spatialkit:::.warned_once), envir = spatialkit:::.warned_once)
  logged <- paste(capture_spatialkit_log(
    fit_rf_model(txt_pts(), "z", "a", num_trees = 50L, include_coords = TRUE)),
    collapse = "\n")
  expect_match(logged, "include_coords = TRUE")
  # The general claim is gone.
  expect_no_match(logged, "fail wherever it has not")
  expect_no_match(logged, "fail badly")
  # What it does say: why the easy scores are optimistic, and what to do.
  expect_match(logged, "out-of-bag error and random cross-validation overstate")
  expect_match(logged, "blocked or distance-matched folds")
  expect_match(logged, "block_kfold")
  expect_match(logged, "nndm")
  expect_match(logged, "can reverse in a region the data never visited")
  expect_match(logged, "once per session")
  # Still a logged caution about a choice, not an R warning.
  rm(list = ls(spatialkit:::.warned_once), envir = spatialkit:::.warned_once)
  expect_no_warning(spatialkit_quiet_fit <- suppressMessages(
    fit_rf_model(txt_pts(), "z", "a", num_trees = 50L, include_coords = TRUE)))
})

test_that("print() says when a forest uses no spatial information", {
  skip_if_not_installed("ranger")
  plain <- fit_rf_model(txt_pts(), "z", "a", num_trees = 50L)
  out <- utils::capture.output(print(plain))
  expect_true(any(grepl("Coords as predictors: no - no spatial information used",
                        out, fixed = TRUE)))
  with_xy <- suppressMessages(
    fit_rf_model(txt_pts(), "z", "a", num_trees = 50L, include_coords = TRUE))
  out_xy <- utils::capture.output(print(with_xy))
  expect_true(any(grepl("Coords as predictors: YES - score with blocked folds",
                        out_xy, fixed = TRUE)))
  expect_false(any(grepl("no spatial information used", out_xy, fixed = TRUE)))
})

test_that("the caution says how to compare the two forests, and the two calls pair up", {
  skip_if_not_installed("ranger")
  rm(list = ls(spatialkit:::.warned_once), envir = spatialkit:::.warned_once)
  logged <- paste(capture_spatialkit_log(
    fit_rf_model(txt_pts(), "z", "a", num_trees = 50L, include_coords = TRUE)),
    collapse = "\n")
  # The advice was "compare it on the same folds", with no word on how.
  expect_match(logged, "pass them as `folds` to two cv_rf\\(\\) calls")
  expect_match(logged, "one with include_coords = TRUE and one without")
  expect_match(logged, "\\?cv_rf shows the two calls")

  # What ?cv_rf and ?fit_rf_model show: folds built once, handed to both
  # calls.  Each fold then holds out the same rows in both, so the per-fold
  # tables pair up row for row.
  d <- txt_pts(n = 160, seed = 2)
  folds <- suppressWarnings(suppressMessages(
    make_folds(d, k = 4, method = "block_kfold", seed = 1)))
  plain <- suppressMessages(cv_rf(d, "z", "a", folds = folds, num_trees = 50L))
  with_xy <- suppressMessages(cv_rf(d, "z", "a", folds = folds, num_trees = 50L,
                                    include_coords = TRUE))
  key <- function(cv) cv$predictions[order(cv$predictions$`..row_id`),
                                     c("..row_id", "fold")]
  expect_identical(key(plain), key(with_xy))
  expect_identical(plain$fold_metrics$fold, with_xy$fold_metrics$fold)
  expect_false(is.null(plain$fold_metrics$n_test))
  expect_identical(plain$fold_metrics$n_test, with_xy$fold_metrics$n_test)
  expect_identical(plain$fold_metrics$n_train, with_xy$fold_metrics$n_train)
  expect_length(plain$fold_metrics$RMSE - with_xy$fold_metrics$RMSE, 4L)
  expect_true(all(is.finite(c(plain$overall$RMSE, with_xy$overall$RMSE))))
  # The forests do differ: one of them split on location.
  expect_false(isTRUE(all.equal(plain$overall$RMSE, with_xy$overall$RMSE)))
})
