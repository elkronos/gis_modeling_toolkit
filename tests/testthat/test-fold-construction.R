# ===========================================================================
# Basic fold construction. See test-fold-methods.R for the spatial methods
# (block, buffered LOO, leave-location-out, NNDM) and
# test-make-folds-row-ids.R for the ..row_id contract.
# ===========================================================================

.fc_points <- function(n, seed = 42) {
  set.seed(seed)
  sf::st_sf(
    id = seq_len(n),
    geometry = sf::st_sfc(
      lapply(seq_len(n), function(i) {
        sf::st_point(c(stats::runif(1, 0, 100), stats::runif(1, 0, 100)))
      }),
      crs = 32632
    )
  )
}

test_that("make_folds random_kfold produces correct fold count", {
  pts <- .fc_points(20)
  folds <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)
  expect_equal(folds$k, 5)
  expect_length(folds$folds, 5)
  # Every point should appear in exactly one test fold
  all_test <- unlist(lapply(folds$folds, `[[`, "test"))
  expect_equal(sort(all_test), 1:20)
})


test_that("random_kfold splits evenly when n divides k", {
  # Without a size assertion, an implementation that put 16 points in fold 1
  # and 1 in each of the other four passes every check above.
  pts   <- .fc_points(20)
  folds <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)

  sizes <- vapply(folds$folds, function(f) length(f$test), integer(1))
  expect_equal(sizes, rep(4L, 5L))

  # Train and test partition the data in every fold: disjoint, and together
  # they are the whole layer.
  for (f in folds$folds) {
    expect_length(intersect(f$train, f$test), 0L)
    expect_equal(sort(c(f$train, f$test)), 1:20)
    expect_equal(length(f$train), 16L)
  }
})


test_that("random_kfold spreads the remainder one row at a time", {
  # 23 into 5 must be 5,5,5,4,4 in some order -- never 7,4,4,4,4.
  pts   <- .fc_points(23, seed = 7)
  folds <- make_folds(pts, k = 5, method = "random_kfold", seed = 7)

  sizes <- vapply(folds$folds, function(f) length(f$test), integer(1))
  expect_equal(sum(sizes), 23L)
  expect_equal(max(sizes) - min(sizes), 1L)
  expect_equal(sort(sizes), c(4L, 4L, 5L, 5L, 5L))

  all_test <- unlist(lapply(folds$folds, `[[`, "test"))
  expect_equal(sort(all_test), 1:23)          # a partition, not a resample
  expect_false(anyDuplicated(all_test) > 0L)
})


test_that("random_kfold is reproducible under a fixed seed and varies without one", {
  pts <- .fc_points(30, seed = 3)

  a <- make_folds(pts, k = 4, method = "random_kfold", seed = 11)
  b <- make_folds(pts, k = 4, method = "random_kfold", seed = 11)
  expect_equal(lapply(a$folds, `[[`, "test"), lapply(b$folds, `[[`, "test"))

  c_ <- make_folds(pts, k = 4, method = "random_kfold", seed = 12)
  expect_false(isTRUE(all.equal(lapply(a$folds, `[[`, "test"),
                                lapply(c_$folds, `[[`, "test"))))
})
