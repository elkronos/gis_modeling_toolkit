# ===========================================================================
# Basic fold construction. See test-fold-methods.R for the spatial methods
# (block, buffered LOO, leave-location-out, NNDM) and
# test-make-folds-row-ids.R for the ..row_id contract.
# ===========================================================================

test_that("make_folds random_kfold produces correct fold count", {
  pts <- sf::st_sf(
    id = 1:20,
    geometry = sf::st_sfc(
      lapply(1:20, function(i) sf::st_point(c(runif(1, 0, 100), runif(1, 0, 100)))),
      crs = 32632
    )
  )
  folds <- make_folds(pts, k = 5, method = "random_kfold", seed = 42)
  expect_equal(folds$k, 5)
  expect_length(folds$folds, 5)
  # Every point should appear in exactly one test fold
  all_test <- unlist(lapply(folds$folds, `[[`, "test"))
  expect_equal(sort(all_test), 1:20)
})
