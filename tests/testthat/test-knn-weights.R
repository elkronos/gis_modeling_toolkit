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


test_that("the dense path handles k = 1 without collapsing the index matrix", {
  # At k = 1 the inner function handed to apply() returns a SCALAR, so apply()
  # simplifies to a length-n vector and t() makes it a 1 x n matrix.  The
  # `nn_idx[i, ]` lookup below then failed with "subscript out of bounds" for
  # every i > 1 -- so residual_morans_i(fit, k = 1) was broken outright on any
  # machine without FNN.
  build <- knn_fn()
  set.seed(11)
  coords <- matrix(runif(40 * 2), ncol = 2)      # jittered: no distance ties

  W <- as.matrix(build(coords, k = 1L, use_fnn = FALSE))
  expect_equal(dim(W), c(40L, 40L))
  # Exactly one neighbour per row, weight 1, never itself.
  expect_equal(unname(rowSums(W)), rep(1, 40), tolerance = 1e-12)
  expect_equal(unname(rowSums(W > 0)), rep(1L, 40))
  expect_true(all(diag(W) == 0))

  # The neighbour really is the nearest one, computed independently.
  d <- as.matrix(stats::dist(coords)); diag(d) <- Inf
  expect_equal(unname(apply(W, 1, which.max)), unname(apply(d, 1, which.min)))

  # ... and the two backends agree at k = 1 as they do at k = 8.
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  expect_equal(unname(as.matrix(build(coords, k = 1L))), unname(W),
               tolerance = 1e-12)
})

test_that("residual_morans_i(k = 1) agrees with the dense k = 1 weights", {
  # The user-visible form of the same crash.  residual_morans_i() has no
  # backend switch of its own, so on this machine it takes the FNN/Matrix path;
  # supplying the DENSE k = 1 matrix through `weights` is what puts the fixed
  # branch under the public function, and the two must give one statistic.
  build <- knn_fn()
  set.seed(12)
  n   <- 60
  xy  <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                      coords = c("x", "y"), crs = 32632)
  fit <- moran_stub_fit(pts, stats::rnorm(n))

  for (k in c(1L, 2L)) {
    lbl   <- paste("k =", k)
    dense <- build(xy, k = k, use_fnn = FALSE)
    got   <- residual_morans_i(fit, k = k)
    ref   <- residual_morans_i(fit, weights = dense)

    expect_true(is.finite(got$observed), info = lbl)
    expect_true(is.finite(got$p_value), info = lbl)
    expect_equal(got$observed, ref$observed, tolerance = 1e-12, info = lbl)
    expect_equal(got$p_value,  ref$p_value,  tolerance = 1e-12, info = lbl)
    # Cliff & Ord's null expectation, which does not depend on the weights.
    expect_equal(got$expected, -1 / (n - 1), tolerance = 1e-12, info = lbl)
  }
})


# ---------------------------------------------------------------------------
# Duplicate coordinates
#
# FNN::get.knn() answers a tied query by returning SOME of the tied points, and
# the query point's own index is eligible to be one of them.  With exact
# duplicates a k-nearest query therefore came back holding self -- putting 1/k
# on the diagonal of a matrix Moran's I is only defined for with a zero
# diagonal -- and asking for k + 1 and dropping self is NOT enough on its own,
# because the slot self occupied displaced a genuine co-located neighbour and
# left a farther point standing in for it.
#
# Repeat observations at one site are exactly what
# make_folds(method = "leave_location_out") exists for, and residual_morans_i()
# reads fit$data_sf without de-duplicating, so this is a mainstream input.
# ---------------------------------------------------------------------------

.dup_layouts <- function() list(
  "25 sites x 4 repeats" = local({
    set.seed(7)
    cbind(rep(runif(25, 0, 1000), each = 4), rep(runif(25, 0, 1000), each = 4))
  }),
  "5 sites x 20 repeats" = local({
    set.seed(11)
    cbind(rep(runif(5, 0, 100), each = 20), rep(runif(5, 0, 100), each = 20))
  }),
  "mixed dup and distinct" = local({
    set.seed(3)
    rbind(cbind(rep(runif(10, 0, 500), each = 3), rep(runif(10, 0, 500), each = 3)),
          cbind(runif(50, 0, 500), runif(50, 0, 500)))
  }),
  "all identical" = cbind(rep(5, 12), rep(9, 12))
)

test_that("duplicate coordinates leave a zero diagonal and row sums of 1", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  build <- knn_fn()
  for (nm in names(.dup_layouts())) {
    co <- .dup_layouts()[[nm]]
    for (k in c(1L, 3L, 8L)) {
      if (k >= nrow(co) - 1L) next
      W <- as.matrix(build(co, k = k, use_fnn = TRUE, use_matrix = TRUE))
      lbl <- paste(nm, "k =", k)
      expect_equal(max(abs(diag(W))), 0, info = lbl)
      expect_equal(unname(rowSums(W)), rep(1, nrow(co)), info = lbl)
      expect_equal(unique(rowSums(W != 0)), k, info = lbl)
    }
  }
})

test_that("with duplicates the kd-tree keeps the k NEAREST, not k of any", {
  # The property a k + 1 query and a self-drop does NOT give you.  Measured on
  # 25 sites x 4 repeats at k = 3, that approach retained 75 of 400 pairs at a
  # distance of 121 where a co-located neighbour at distance 0 existed.
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  build <- knn_fn()
  for (nm in names(.dup_layouts())) {
    co <- .dup_layouts()[[nm]]
    D  <- as.matrix(stats::dist(co)); diag(D) <- Inf
    for (k in c(1L, 2L, 3L, 8L)) {
      if (k >= nrow(co) - 1L) next
      W <- as.matrix(build(co, k = k, use_fnn = TRUE, use_matrix = TRUE))
      lbl <- paste(nm, "k =", k)
      for (i in seq_len(nrow(co))) {
        expect_equal(sort(unname(D[i, W[i, ] != 0])),
                     unname(sort(D[i, ])[seq_len(k)]), info = lbl)
      }
    }
  }
})

test_that("the kd-tree and the dense fallback agree on the same weights", {
  # The statistic must not depend on whether FNN happens to be installed.  On
  # distinct coordinates the two matrices are identical; with ties they may
  # break them differently, so the retained NEIGHBOUR DISTANCES are compared.
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  build <- knn_fn()
  set.seed(7)
  distinct <- cbind(runif(200, 0, 1000), runif(200, 0, 1000))
  for (k in c(1L, 3L, 8L)) {
    A <- as.matrix(build(distinct, k = k, use_fnn = TRUE,  use_matrix = TRUE))
    B <- build(distinct, k = k, use_fnn = FALSE, use_matrix = FALSE)
    expect_equal(A, B, ignore_attr = TRUE, info = paste("distinct, k =", k))
  }
  for (nm in names(.dup_layouts())) {
    co <- .dup_layouts()[[nm]]
    D  <- as.matrix(stats::dist(co)); diag(D) <- Inf
    for (k in c(1L, 3L)) {
      if (k >= nrow(co) - 1L) next
      A <- as.matrix(build(co, k = k, use_fnn = TRUE,  use_matrix = TRUE))
      B <- build(co, k = k, use_fnn = FALSE, use_matrix = FALSE)
      expect_equal(sort(D[A != 0]), sort(D[B != 0]),
                   info = paste(nm, "k =", k))
    }
  }
})

test_that("k = 1 works on the dense path", {
  # apply() simplifies a length-1 result to a VECTOR, so t() made a 1 x n
  # matrix and nn_idx[i, ] was out of bounds for every i > 1.
  build <- knn_fn()
  set.seed(5)
  co <- cbind(runif(20, 0, 100), runif(20, 0, 100))
  W  <- build(co, k = 1L, use_fnn = FALSE, use_matrix = FALSE)
  expect_equal(unname(rowSums(W)), rep(1, 20))
  expect_equal(max(abs(diag(W))), 0)
})
