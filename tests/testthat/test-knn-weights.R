# tests/testthat/test-knn-weights.R
# ---------------------------------------------------------------------------
# .build_knn_weights() has two implementations -- an FNN kd-tree + sparse
# Matrix path and a dense O(n^2) fallback -- and until now only one of them
# could ever run.  The existing guard test skipped whenever FNN was installed,
# so on any developer machine with FNN the fallback and its size guard were
# dead code as far as the suite was concerned.
#
# The backend choice is now a parameter, so both paths can be exercised and,
# more importantly, compared against each other.
# ---------------------------------------------------------------------------

knn_fn <- function() get(".build_knn_weights", envir = asNamespace("spatialkit"))

test_that("the dense fallback refuses to allocate an n x n matrix", {
  # Previously unreachable: skipped whenever FNN was present.
  build <- knn_fn()
  set.seed(3)
  big <- matrix(runif(5001 * 2), ncol = 2)
  expect_error(build(big, k = 8L, use_fnn = FALSE),
               "requires FNN for k-NN weights")
})

test_that("the dense fallback is allowed below the size guard", {
  build <- knn_fn()
  set.seed(1)
  coords <- matrix(runif(200 * 2), ncol = 2)
  W <- build(coords, k = 5L, use_fnn = FALSE)
  expect_equal(dim(W), c(200L, 200L))
  expect_true(all(is.finite(as.matrix(W))))
})

test_that("weights are row-standardised on both paths", {
  build <- knn_fn()
  set.seed(2)
  coords <- matrix(runif(150 * 2), ncol = 2)

  W_dense <- build(coords, k = 6L, use_fnn = FALSE)
  expect_equal(unname(rowSums(as.matrix(W_dense))), rep(1, 150), tolerance = 1e-12)

  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  W_sparse <- build(coords, k = 6L)
  expect_equal(unname(rowSums(as.matrix(W_sparse))), rep(1, 150), tolerance = 1e-12)
})

test_that("the sparse and dense paths produce the same weights", {
  # The property that actually matters and has never been checked: two
  # implementations of the same thing must agree.
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  build <- knn_fn()

  set.seed(3)
  coords <- matrix(runif(300 * 2), ncol = 2)   # jittered, so no distance ties
  W_sparse <- as.matrix(build(coords, k = 8L))
  W_dense  <- as.matrix(build(coords, k = 8L, use_fnn = FALSE))

  expect_equal(dim(W_sparse), dim(W_dense))
  expect_equal(unname(W_sparse), unname(W_dense), tolerance = 1e-12)
})

test_that("k is clamped to n - 1", {
  build <- knn_fn()
  set.seed(4)
  coords <- matrix(runif(6 * 2), ncol = 2)
  W <- build(coords, k = 50L, use_fnn = FALSE)
  expect_equal(dim(W), c(6L, 6L))
  expect_equal(unname(rowSums(as.matrix(W))), rep(1, 6), tolerance = 1e-12)
  expect_true(all(diag(as.matrix(W)) == 0))   # no self-neighbours
})
