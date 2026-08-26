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

test_that(".morans_i_for_k returns the analytic value when every cell is a neighbour", {
  # The weight matrix uses k = min(8, n_cells - 1) neighbours, so with 9 or
  # fewer cells EVERY other cell is a neighbour and W is the complete
  # row-standardised matrix W_ij = 1/(n-1).  Then W %*% e = -e/(n-1) for any
  # mean-zero residual vector, S0 = n, and Moran's I collapses to exactly
  # -1/(n - 1) -- independent of the data.  That is a closed-form expectation,
  # and it is the configuration the sparse crossprod() bug crashed on.
  fx <- .mifk_fixture()
  mk <- spatialkit:::.morans_i_for_k

  for (n_cells in c(4L, 5L, 6L, 9L)) {
    clusters <- rep(seq_len(n_cells), length.out = nrow(fx$xy))
    got <- mk(fx$xy, fx$spatial, fx$pred, clusters)
    expect_equal(got, -1 / (n_cells - 1), tolerance = 1e-10,
                 info = paste("n_cells =", n_cells))
    expect_true(is.finite(got))
  }
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

  expect_equal(spatialkit:::.morans_i_for_k(fx$xy, fx$spatial, fx$pred,
                                            clusters),
               expected, tolerance = 1e-10)

  # And the statistic discriminates: a response with no unexplained spatial
  # structure sits near zero, the trended one well above it.
  flat_I <- spatialkit:::.morans_i_for_k(fx$xy, fx$flat, fx$pred, clusters)
  expect_gt(expected, 0.2)
  expect_lt(abs(flat_I), 0.2)
  expect_gt(expected, flat_I)
})

test_that(".morans_i_for_k returns NA rather than a number it cannot justify", {
  fx <- .mifk_fixture()
  mk <- spatialkit:::.morans_i_for_k

  # Fewer than 4 cells: no usable weight graph.
  expect_true(is.na(mk(fx$xy, fx$spatial, fx$pred,
                       rep(1:3, length.out = nrow(fx$xy)))))
  # Non-finite cell means leave fewer than 4 usable cells.
  broken <- fx$spatial; broken[] <- NA_real_
  expect_true(is.na(mk(fx$xy, broken, fx$pred,
                       rep(1:6, length.out = nrow(fx$xy)))))
  # A perfectly explained response leaves zero residual variance.
  exact <- 3 + 2 * fx$pred[, 1]
  expect_true(is.na(mk(fx$xy, exact, fx$pred,
                       rep(1:6, length.out = nrow(fx$xy)))))
})


# ---------------------------------------------------------------------------
# determine_optimal_levels()
# ---------------------------------------------------------------------------

.dol_two_clusters <- function(seed = 1) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = c(runif(25, 0, 10), runif(25, 90, 100)),
               y = c(runif(25, 0, 10), runif(25, 90, 100))),
    coords = c("x", "y"), crs = 32632
  )
}

.dol_model_points <- function(n = 150, seed = 3) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 0.01 * sf::st_coordinates(d)[, 1] + 2 * d$a + rnorm(n)
  d
}

test_that("determine_optimal_levels finds the elbow of two separated clusters", {
  # Two tight, far-apart clusters: WSS collapses between k = 1 and k = 2 and
  # is nearly flat afterwards, so the knee is at 2 and the candidate set is
  # its immediate neighbourhood.
  out <- determine_optimal_levels(.dol_two_clusters(), max_levels = 6)

  expect_type(out, "integer")
  expect_true(2L %in% out)
  expect_true(all(out >= 1L & out <= 6L))
  # Under "geometric" the candidate set is the knee plus two neighbours, so at
  # most three values come back however large top_n is.
  expect_lte(length(out), 3L)
  expect_lte(length(determine_optimal_levels(.dol_two_clusters(),
                                             max_levels = 6, top_n = 10)), 3L)
  expect_null(attr(out, "diagnostics"))     # no model-aware run, no diagnostics
})

test_that("determine_optimal_levels runs the model-aware criteria without crashing", {
  # The regression: with FNN and Matrix installed, .morans_i_for_k() built a
  # sparse weight matrix and base::crossprod() refused it, so every call with
  # a model-aware criterion errored for >= 4 cells.
  d <- .dol_model_points()

  for (crit in c("morans_i", "combined")) {
    out <- determine_optimal_levels(d, max_levels = 8, response_var = "z",
                                    predictor_vars = "a", criterion = crit)
    expect_type(out, "integer")
    expect_gte(length(out), 1L)
    expect_true(all(out >= 1L & out <= 8L))

    diag <- attr(out, "diagnostics")
    expect_false(is.null(diag), info = crit)
    expect_true(any(is.finite(diag$moran_i)), info = crit)
    expect_true(all(is.finite(diag$wss)), info = crit)
    # Moran's I was actually evaluated over the elbow neighbourhood: finite for
    # every k that yields at least 4 cells, NA below that (the documented
    # n_cells < 4 guard in .morans_i_for_k()).
    big   <- diag$eval_ks[diag$eval_ks >= 4L]
    small <- diag$eval_ks[diag$eval_ks < 4L]
    expect_gt(length(big), 0L)
    expect_true(all(is.finite(diag$moran_i[big])), info = crit)
    if (length(small))
      expect_true(all(is.na(diag$moran_i[small])), info = crit)
    # Every k the neighbourhood did NOT evaluate stays NA.
    expect_true(all(is.na(diag$moran_i[setdiff(seq_along(diag$moran_i),
                                               diag$eval_ks)])), info = crit)
  }
})

test_that("determine_optimal_levels auto-upgrades and falls back with a log line", {
  d <- .dol_model_points()

  # Supplying model variables under the default criterion upgrades to combined.
  up <- capture_spatialkit_log(
    out <- determine_optimal_levels(d, max_levels = 8, response_var = "z",
                                    predictor_vars = "a"),
    level = logger::INFO
  )
  expect_true(log_has(up, "using combined criterion"))
  expect_equal(attr(out, "diagnostics")$criterion, "combined")

  # Asking for a model-aware criterion without the variables falls back.
  down <- capture_spatialkit_log(
    geo <- determine_optimal_levels(d, max_levels = 6, criterion = "morans_i")
  )
  expect_true(log_has(down, "requires response_var and predictor_vars"))
  expect_null(attr(geo, "diagnostics"))
  expect_equal(geo, determine_optimal_levels(d, max_levels = 6))
})

test_that("determine_optimal_levels validates input and degenerate geometry", {
  d <- .dol_model_points()

  expect_error(determine_optimal_levels(sf::st_drop_geometry(d)),
               "must be an sf object")

  d$fac <- factor(sample(letters[1:3], nrow(d), replace = TRUE))
  expect_error(
    determine_optimal_levels(d, response_var = "z", predictor_vars = "fac"),
    "`predictor_vars` must be numeric"
  )

  # Fewer than 3 rows: one level, no clustering attempted.
  expect_equal(determine_optimal_levels(d[1:2, ], max_levels = 4), 1L)

  # Only two distinct positions: k_max clamps to n_unique - 1 = 1, so again 1.
  dup <- sf::st_as_sf(
    data.frame(x = rep(c(0, 10), 10), y = rep(c(0, 10), 10)),
    coords = c("x", "y"), crs = 32632)
  expect_equal(determine_optimal_levels(dup, max_levels = 6), 1L)
})

test_that("determine_optimal_levels coerces MULTIPOINT before reading coordinates", {
  # st_coordinates() returns one row per VERTEX, so a two-vertex MULTIPOINT
  # layer produced an xy matrix twice as tall as resp_vec/pred_mat and every
  # index below read a different feature than it thought.  Coercing to
  # representative points first makes the MULTIPOINT layer and its own
  # centroids give the same answer.
  set.seed(21)
  n <- 120
  cx <- runif(n, 0, 1000); cy <- runif(n, 0, 1000)
  mp <- sf::st_sfc(lapply(seq_len(n), function(i) {
    sf::st_multipoint(rbind(c(cx[i] - 5, cy[i]), c(cx[i] + 5, cy[i])))
  }), crs = 32632)

  multi <- sf::st_sf(a = rnorm(n), geometry = mp)
  multi$z <- 0.01 * cx + 2 * multi$a + rnorm(n)

  pts <- multi
  sf::st_geometry(pts) <- sf::st_sfc(
    lapply(seq_len(n), function(i) sf::st_point(c(cx[i], cy[i]))), crs = 32632)

  # One coordinate row per feature after coercion -- the property everything
  # below depends on.
  coerced <- coerce_to_points(multi, "auto")
  expect_equal(nrow(sf::st_coordinates(coerced)), n)
  expect_equal(unname(sf::st_coordinates(coerced)[, 1]), cx, tolerance = 1e-6)

  geo_multi <- determine_optimal_levels(multi, max_levels = 6)
  geo_pts   <- determine_optimal_levels(pts,   max_levels = 6)
  expect_equal(geo_multi, geo_pts)

  mod_multi <- determine_optimal_levels(multi, max_levels = 6,
                                        response_var = "z",
                                        predictor_vars = "a",
                                        criterion = "combined")
  mod_pts   <- determine_optimal_levels(pts, max_levels = 6,
                                        response_var = "z",
                                        predictor_vars = "a",
                                        criterion = "combined")
  expect_equal(as.integer(mod_multi), as.integer(mod_pts))
  expect_equal(attr(mod_multi, "diagnostics")$moran_i,
               attr(mod_pts, "diagnostics")$moran_i)
})

test_that("determine_optimal_levels restores the caller's RNG stream", {
  d <- .dol_model_points()
  set.seed(777); expected <- runif(3)
  set.seed(777)
  invisible(determine_optimal_levels(d, max_levels = 6, set_seed = 42L))
  expect_equal(runif(3), expected)
})
