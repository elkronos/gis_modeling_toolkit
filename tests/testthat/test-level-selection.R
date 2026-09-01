# ===========================================================================
# R/level-selection.R: determine_optimal_levels(), .elbow_from_wss() and
# .morans_i_for_k().
#
# Zero coverage is why a sparse-Matrix crossprod() crash in .morans_i_for_k()
# -- the same one residual_morans_i() had -- went unnoticed and made
# determine_optimal_levels(criterion != "geometric") unusable for >= 4 cells
# on any machine with FNN and Matrix installed.
# ===========================================================================

# ---------------------------------------------------------------------------
# .elbow_from_wss(): a pure numeric function with an analytic answer
# ---------------------------------------------------------------------------

# A WSS curve with a hand-placed knee at k = 4: a steep linear fall from
# k = 1 to k = 4, then an almost flat tail.  Under the perpendicular-distance
# (kneedle) rule the maximum deviation from the (k_min, k_max) chord is at the
# corner, by construction.
.knee_curve <- c(100, 70, 40, 10, 9, 8, 7, 6, 5, 4)

test_that(".elbow_from_wss finds a hand-placed knee", {
  eb <- spatialkit:::.elbow_from_wss
  out <- eb(.knee_curve)

  expect_equal(out$knee_k, 4L)
  # Neighbours of the knee, clamped into [min_k, max_k].
  expect_equal(out$candidates, c(3L, 4L, 5L))
  expect_equal(out$diagnostics$wss, .knee_curve)
  expect_equal(out$diagnostics$d1, diff(.knee_curve))
  expect_equal(out$diagnostics$d2, diff(diff(.knee_curve)))

  # return_neighbors = FALSE gives the knee alone.
  expect_equal(eb(.knee_curve, return_neighbors = FALSE)$candidates, 4L)

  # Moving the corner moves the answer -- so 4 is not a fixed point of the
  # implementation.
  later <- c(100, 90, 80, 70, 60, 50, 10, 9, 8, 7)
  expect_equal(eb(later)$knee_k, 7L)
  earlier <- c(100, 20, 19, 18, 17, 16, 15, 14, 13, 12)
  expect_equal(eb(earlier)$knee_k, 2L)
})

test_that(".elbow_from_wss is invariant to an affine rescaling of WSS", {
  # Both axes are min-max normalised before the perpendicular distance is
  # taken, so a * wss + b (a > 0) cannot move the knee.  WSS is in squared CRS
  # units, so this is what makes the answer independent of the projection.
  eb <- spatialkit:::.elbow_from_wss
  base <- eb(.knee_curve)$knee_k
  for (a in c(1e-6, 0.5, 1000)) {
    for (b in c(-50, 0, 1e5)) {
      expect_equal(eb(a * .knee_curve + b)$knee_k, base,
                   info = sprintf("a = %g, b = %g", a, b))
    }
  }
})

test_that(".elbow_from_wss honours the min_k / max_k window", {
  eb <- spatialkit:::.elbow_from_wss

  # The knee is inside the window, so it is unchanged ...
  expect_equal(eb(.knee_curve, min_k = 3L, max_k = 8L)$knee_k, 4L)
  # ... and the candidate set stays inside it.
  cand <- eb(.knee_curve, min_k = 4L, max_k = 6L)$candidates
  expect_true(all(cand >= 4L & cand <= 6L))
  # Only the windowed slice is reported back.
  expect_equal(eb(.knee_curve, min_k = 3L, max_k = 6L)$diagnostics$wss,
               .knee_curve[3:6])

  # Exactly two k values: too few for a chord, so the midpoint is returned.
  expect_equal(eb(c(5, 3))$knee_k, floor((1 + 2) / 2))
  expect_equal(eb(.knee_curve, min_k = 4L, max_k = 5L)$knee_k,
               floor((4 + 5) / 2))
})

test_that(".elbow_from_wss rejects input it cannot use", {
  eb <- spatialkit:::.elbow_from_wss
  expect_error(eb(5), "must be numeric length >= 2")
  expect_error(eb(character(2)), "must be numeric length >= 2")
  expect_error(eb(c(1, 2), min_k = 2L), "need at least two k values")
  expect_error(eb(.knee_curve, min_k = 6L, max_k = 5L),
               "need at least two k values")
})


# ---------------------------------------------------------------------------
# .morans_i_for_k()
# ---------------------------------------------------------------------------

.mifk_fixture <- function(n = 200, seed = 5) {
  set.seed(seed)
  xy   <- cbind(runif(n, 0, 100), runif(n, 0, 100))
  pred <- matrix(rnorm(n), ncol = 1)
  list(
    xy   = xy,
    pred = pred,
    # A strong east-west trend the single predictor cannot absorb, so it
    # survives into the cell-level residuals.
    spatial = 0.05 * xy[, 1] + 0.5 * pred[, 1] + rnorm(n, 0, 0.2),
    flat    = 0.5 * pred[, 1] + rnorm(n, 0, 1)
  )
}

test_that(".morans_i_for_k refuses the complete-graph value below the floor", {
  # The weight matrix uses k = min(8, n_cells - 1) neighbours, so with 9 or
  # fewer cells EVERY other cell is a neighbour and W is the complete
  # row-standardised matrix W_ij = 1/(n-1).  Then W %*% e = -e/(n-1) for any
  # mean-zero residual vector, S0 = n, and Moran's I collapses to exactly
  # -1/(n - 1): a function of the CELL COUNT ALONE that no data set can move.
  #
  # Returning that number would not merely be uninformative.  |I| = 1/(n - 1)
  # falls monotonically in n for arithmetic reasons, so criterion = "morans_i"
  # would rank the largest candidate best every single time.  The function
  # refuses it and returns NA, and determine_optimal_levels() drops those
  # candidates from the model-aware ranking.
  fx <- .mifk_fixture()
  mk <- spatialkit:::.morans_i_for_k

  for (n_cells in c(4L, 5L, 6L, 9L)) {
    clusters <- rep(seq_len(n_cells), length.out = nrow(fx$xy))
    lbl <- paste("n_cells =", n_cells)
    got  <- mk(fx$xy, fx$spatial, fx$pred, clusters)
    flat <- mk(fx$xy, fx$flat,    fx$pred, clusters)
    # The return is c(I = , z = ): BOTH must be NA, since z is the quantity
    # the ranking uses and an NA I beside a finite z would still be ranked.
    expect_identical(names(got), c("I", "z"), info = lbl)
    expect_true(all(is.na(got)),  info = lbl)
    # Refused for the flat response too -- the point is that the complete
    # graph cannot tell these two apart, not that one of them is unusable.
    expect_true(all(is.na(flat)), info = lbl)
  }
})

test_that(".morans_i_for_k reports a data-dependent value from ten cells up", {
  # Ten cells is the first count above the floor: k = min(8, 9) = 8 of the 9
  # other cells are neighbours, so W is no longer complete and I stops being a
  # function of n.  A spatially coherent 5 x 2 partition, so the cell
  # centroids are genuinely spread out rather than piled at the domain centre.
  fx <- .mifk_fixture()
  mk <- spatialkit:::.morans_i_for_k
  clusters <- as.integer(cut(fx$xy[, 1], 5)) +
    5L * (as.integer(cut(fx$xy[, 2], 2)) - 1L)
  expect_equal(length(unique(clusters)), 10L)

  # An independent reference: the dense (non-sparse) weight path plus the
  # Cliff & Ord formula written out, so it shares no code with the branch
  # under test.
  ref_I <- function(resp) {
    ids <- sort(unique(clusters)); nc <- length(ids)
    cr <- numeric(nc); cxy <- matrix(0, nc, 2); cp <- matrix(0, nc, 1)
    for (j in seq_along(ids)) {
      m <- clusters == ids[j]
      cr[j]   <- mean(resp[m])
      cxy[j, ] <- colMeans(fx$xy[m, , drop = FALSE])
      cp[j, ]  <- colMeans(fx$pred[m, , drop = FALSE])
    }
    r  <- stats::lm.fit(x = cbind(1, cp), y = cr)$residuals
    W  <- spatialkit:::.build_knn_weights(cxy, k = min(8L, nc - 1L),
                                          use_fnn = FALSE, use_matrix = FALSE)
    rc <- r - mean(r)
    (nc / sum(W)) * sum(rc * (W %*% rc)) / sum(rc^2)
  }

  I_spatial <- mk(fx$xy, fx$spatial, fx$pred, clusters)[["I"]]
  I_flat    <- mk(fx$xy, fx$flat,    fx$pred, clusters)[["I"]]

  expect_true(is.finite(I_spatial))
  expect_true(is.finite(I_flat))
  expect_equal(I_spatial, ref_I(fx$spatial), tolerance = 1e-10)
  expect_equal(I_flat,    ref_I(fx$flat),    tolerance = 1e-10)

  # The property the sub-floor case cannot have: two different residual
  # vectors on the SAME clustering give two different statistics.
  expect_false(isTRUE(all.equal(I_spatial, I_flat, tolerance = 1e-6)))
  # And neither of them is the complete-graph constant the floor rejects.
  expect_false(isTRUE(all.equal(I_spatial, -1 / 9, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(I_flat,    -1 / 9, tolerance = 1e-6)))

  # One cell fewer is below the floor and comes back NA, so ten is the
  # boundary and not an arbitrary choice.
  nine <- as.integer(cut(fx$xy[, 1], 3)) +
    3L * (as.integer(cut(fx$xy[, 2], 3)) - 1L)
  expect_equal(length(unique(nine)), 9L)
  expect_true(all(is.na(mk(fx$xy, fx$spatial, fx$pred, nine))))
})

test_that(".morans_i_for_k matches a dense hand computation for many cells", {
  # Above 9 cells the neighbour set is a genuine 8-NN graph and I depends on
  # the data.  The reference below builds the weights on the DENSE path
  # (use_fnn/use_matrix = FALSE) and applies the Cliff & Ord formula directly,
  # so it does not share the sparse code path under test.
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  fx <- .mifk_fixture()
  clusters <- as.integer(cut(fx$xy[, 1], 5)) +
    5L * (as.integer(cut(fx$xy[, 2], 4)) - 1L)

  ids <- sort(unique(clusters))
  nc  <- length(ids)
  expect_gt(nc, 9L)

  cell_resp <- numeric(nc)
  cell_xy   <- matrix(0, nc, 2)
  cell_pred <- matrix(0, nc, 1)
  for (j in seq_along(ids)) {
    m <- clusters == ids[j]
    cell_resp[j]   <- mean(fx$spatial[m])
    cell_xy[j, ]   <- colMeans(fx$xy[m, , drop = FALSE])
    cell_pred[j, ] <- colMeans(fx$pred[m, , drop = FALSE])
  }
  resid <- stats::lm.fit(x = cbind(1, cell_pred), y = cell_resp)$residuals
  W <- spatialkit:::.build_knn_weights(cell_xy, k = min(8L, nc - 1L),
                                       use_fnn = FALSE, use_matrix = FALSE)
  rc  <- resid - mean(resid)
  expected <- (nc / sum(W)) * sum(rc * (W %*% rc)) / sum(rc^2)

  got <- spatialkit:::.morans_i_for_k(fx$xy, fx$spatial, fx$pred, clusters)
  expect_equal(got[["I"]], expected, tolerance = 1e-10)

  # And the statistic discriminates: a response with no unexplained spatial
  # structure sits near zero, the trended one well above it.
  flat <- spatialkit:::.morans_i_for_k(fx$xy, fx$flat, fx$pred, clusters)
  expect_gt(expected, 0.2)
  expect_lt(abs(flat[["I"]]), 0.2)
  expect_gt(expected, flat[["I"]])

  # The z companion is what determine_optimal_levels() ranks on, and it must
  # be the Cliff & Ord standardised deviate of the SAME I -- recomputed here
  # from the dense weights so the two paths share no code.
  X   <- cbind(1, cell_pred)
  mom <- spatialkit:::.morans_residual_moments(W = W, X = X, S0 = sum(W),
                                               is_sparse = FALSE)
  expect_false(is.null(mom))
  expect_equal(got[["z"]], (expected - mom$EI) / sqrt(mom$VI), tolerance = 1e-10)
  # The trended response is the significant one; the flat one is not.
  expect_gt(got[["z"]], 2)
  expect_lt(abs(flat[["z"]]), 2)
})

test_that(".morans_i_for_k returns NA rather than a number it cannot justify", {
  fx <- .mifk_fixture()
  mk <- spatialkit:::.morans_i_for_k

  # Every refusal keeps the c(I, z) shape: a bare NA_real_ from one early exit
  # and a length-2 vector from another would make moran_z[k] pick up the I of
  # the next candidate.
  expect_na2 <- function(v) {
    expect_identical(names(v), c("I", "z"))
    expect_true(all(is.na(v)))
  }
  # Fewer than 4 cells: no usable weight graph.
  expect_na2(mk(fx$xy, fx$spatial, fx$pred, rep(1:3, length.out = nrow(fx$xy))))
  # Non-finite cell means leave fewer than 4 usable cells.
  broken <- fx$spatial; broken[] <- NA_real_
  expect_na2(mk(fx$xy, broken, fx$pred, rep(1:6, length.out = nrow(fx$xy))))
  # A perfectly explained response leaves zero residual variance.
  exact <- 3 + 2 * fx$pred[, 1]
  expect_na2(mk(fx$xy, exact, fx$pred, rep(1:6, length.out = nrow(fx$xy))))
})


test_that("the ranking statistic is flat in k where |I| is not", {
  # The reason determine_optimal_levels() ranks on z rather than |I|.  E[I] and
  # Var(I) both move with the cell count, so across datasets with NO spatial
  # structure the SAMPLING DISTRIBUTION of |I| shrinks as k grows -- an |I|
  # ranking then prefers the finest candidate for arithmetic reasons alone.
  #
  # Fresh data per replicate is load-bearing.  Holding one dataset and only
  # re-clustering it measures that realisation's own residual pattern, not the
  # sampling distribution, and shows no such trend (it can even run the other
  # way).  The claim is about the statistic, so the data must be resampled.
  skip_on_cran()
  mk <- spatialkit:::.morans_i_for_k

  grab <- function(k, reps = 25L) {
    v <- vapply(seq_len(reps), function(r) {
      set.seed(7000 + r)
      n  <- 600
      xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
      pr <- cbind(rnorm(n), rnorm(n))
      rs <- as.numeric(1 + pr %*% c(1, -1) + rnorm(n))   # no spatial structure
      cl <- stats::kmeans(xy, centers = k, iter.max = 50, nstart = 3)$cluster
      mk(xy, rs, pr, cl)
    }, c(I = 0, z = 0))
    c(absI = mean(abs(v["I", ]), na.rm = TRUE),
      absZ = mean(abs(v["z", ]), na.rm = TRUE))
  }
  lo <- grab(12L)
  hi <- grab(45L)

  # |I| falls materially across the window (measured ratio 0.60-0.65) ...
  expect_lt(hi[["absI"]], 0.8 * lo[["absI"]])
  # ... while |z| stays in the neighbourhood of E|N(0,1)| = 0.798 at both ends,
  # so the two candidates are compared on the same scale.
  for (v in c(lo[["absZ"]], hi[["absZ"]])) {
    expect_gt(v, 0.5)
    expect_lt(v, 1.3)
  }
  # And the drift in |z| is small next to the drift in |I|.
  expect_lt(abs(hi[["absZ"]] / lo[["absZ"]] - 1),
            abs(hi[["absI"]] / lo[["absI"]] - 1))
})
