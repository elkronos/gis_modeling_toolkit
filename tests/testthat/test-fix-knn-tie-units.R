# The k-nearest-neighbour weights, and what is read from them, are the same
# in every unit.
#
# .knn_pairs() took two distances as one when they differed by less than
# 1e-9 * max(1, max|coords|) map units.  The floor of 1 tied every distance
# difference under 1e-9 map units, so on a layer whose coordinates span less
# than about 1e-5 units more and more neighbours tied as the unit shrank, up
# to every other point with equal weight: Moran's I went to -1 / (n - 1) and
# estimate_sac_range() refused the layer with "no spatial structure".  The
# tolerance is now 1e-9 of the coordinates' size, which a power-of-two
# rescale scales exactly, and on any layer whose coordinates reach 1 map unit
# it is the tolerance it was.  compare_models()'s check that two fits have
# the same rows had the same floor on the coordinates (a layer with no CRS)
# and on the response, as had the check that folds were built on the data
# (.check_fold_probe()) and the edge margin .edge_tol(), by which seeds are
# read as the centres a count record holds.

# Registered at top level: a method registered inside test_that() is
# dropped, by testthat's withr scoping, before the next test runs.
residuals.kt_lm_fit <- function(object, ...) as.numeric(stats::residuals(object$engine))
fitted.kt_lm_fit <- function(object, ...) as.numeric(stats::fitted(object$engine))
registerS3method("residuals", "kt_lm_fit", residuals.kt_lm_fit)
registerS3method("fitted", "kt_lm_fit", fitted.kt_lm_fit)

# 120 points in a 5 km square with a smooth field and a predictor.
.kt_layer <- function(n = 120, seed = 4) {
  withr::with_seed(seed, {
    xy <- cbind(runif(n, 0, 5000), runif(n, 0, 5000))
    w <- rnorm(n)
    y <- 1 + 2 * w + sin(xy[, 1] / 1500) + rnorm(n, 0, 0.5)
  })
  list(xy = xy, w = w, y = y)
}

.kt_sf <- function(xy, y, w, crs = 32615) {
  sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], y = y, p = w),
               coords = c("X", "Y"), crs = crs)
}

.kt_fit <- function(d) {
  new_spatial_fit("kt_lm_fit", stats::lm(y ~ p, data = sf::st_drop_geometry(d)),
                  y ~ p, "y", "p", d)
}

# The weights the rule gave before the change: a dense scan with the floor of
# 1 on the absolute tolerance.  On layers whose coordinates reach 1 map unit
# the rule now gives these weights bit for bit.
.kt_old_rule <- function(coords, k) {
  n <- nrow(coords)
  tol_abs <- 1e-9 * max(1, max(abs(coords)))
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    d <- sqrt((coords[, 1L] - coords[i, 1L])^2 + (coords[, 2L] - coords[i, 2L])^2)
    d[i] <- Inf
    dk <- sort(d, partial = k)[k]
    tied <- abs(d - dk) <= dk * 1e-9 + tol_abs
    close <- !tied & d < dk
    m <- sum(close); t <- sum(tied)
    W[i, close] <- 1 / k
    W[i, tied] <- (k - m) / (t * k)
  }
  W
}

.kt_powers <- c(60, 30, 0, -20, -30, -40, -45, -50, -60)

test_that(".build_knn_weights() gives the same weights at every power-of-two scale", {
  lay <- .kt_layer()
  for (k in c(4L, 8L)) {
    ref <- as.matrix(.build_knn_weights(lay$xy, k = k))
    expect_equal(max(rowSums(ref > 0)), k)
    dense_ref <- .build_knn_weights(lay$xy, k = k, use_fnn = FALSE, use_matrix = FALSE)
    expect_identical(dense_ref, ref)
    for (p in .kt_powers) {
      s <- 2^p
      expect_identical(as.matrix(.build_knn_weights(lay$xy * s, k = k)), ref,
                       info = sprintf("k = %d, coordinates * 2^%d", k, p))
      expect_identical(.build_knn_weights(lay$xy * s, k = k, use_fnn = FALSE,
                                          use_matrix = FALSE), ref,
                       info = sprintf("k = %d, dense, coordinates * 2^%d", k, p))
    }
  }
})

test_that("tied neighbours on a lattice share their slots at every scale", {
  g <- as.matrix(expand.grid(x = 0:9 * 250, y = 0:7 * 250))
  inner <- which(g[, 1] == 1000 & g[, 2] == 750)
  ref <- as.matrix(.build_knn_weights(g, k = 6))
  # Four at 250 m take 1/6 each; the four diagonals tie for the last two
  # slots and take 1/12 each.
  expect_equal(sort(ref[inner, ref[inner, ] > 0]), c(rep(1 / 12, 4), rep(1 / 6, 4)))
  for (p in .kt_powers) {
    s <- 2^p
    got <- as.matrix(.build_knn_weights(g * s, k = 6))
    expect_identical(got, ref, info = sprintf("lattice * 2^%d", p))
    expect_identical(.build_knn_weights(g * s, k = 6, use_fnn = FALSE, use_matrix = FALSE),
                     ref, info = sprintf("lattice, dense, * 2^%d", p))
  }
})

test_that("on layers in metres or feet the weights are those of the rule before", {
  lay <- .kt_layer()
  utm <- cbind(5e5 + lay$xy[, 1], 4.5e6 + lay$xy[, 2])
  feet <- utm / 0.3048
  g <- as.matrix(expand.grid(x = 0:9 * 250, y = 0:7 * 250))
  g_utm <- cbind(5e5 + g[, 1], 4.5e6 + g[, 2])
  layers <- list(xy = lay$xy, utm = utm, feet = feet, g = g, g_utm = g_utm)
  for (nm in names(layers)) {
    co <- layers[[nm]]
    for (k in c(4L, 6L, 8L)) {
      old <- .kt_old_rule(co, k)
      expect_identical(as.matrix(.build_knn_weights(co, k = k)), old,
                       info = sprintf("%s, k = %d", nm, k))
      expect_identical(.build_knn_weights(co, k = k, use_fnn = FALSE, use_matrix = FALSE),
                       old, info = sprintf("%s, dense, k = %d", nm, k))
    }
  }
})

test_that("residual_morans_i() gives one I, z and p at every scale", {
  lay <- .kt_layer()
  ref <- suppressWarnings(residual_morans_i(.kt_fit(.kt_sf(lay$xy, lay$y, lay$w)), k = 4))
  expect_true(ref$observed > 0.2)
  for (p in .kt_powers) {
    d <- .kt_sf(lay$xy * 2^p, lay$y, lay$w)
    got <- suppressWarnings(residual_morans_i(.kt_fit(d), k = 4))
    for (nm in c("observed", "expected", "sd", "z", "p_value"))
      expect_identical(got[[nm]], ref[[nm]], info = sprintf("%s at 2^%d", nm, p))
  }
})

test_that("estimate_sac_range() finds the range of 2^0 times the scale in small units", {
  lay <- .kt_layer()
  run <- function(p) {
    d <- .kt_sf(lay$xy * 2^p, lay$y, lay$w)
    suppressWarnings(suppressMessages(estimate_sac_range(d, "y", predictor_vars = "p")))
  }
  ref <- run(0)
  expect_true(is.finite(as.numeric(ref)))
  expect_null(attr(ref, "rejected_reason"))
  for (p in c(-40, -50, -60)) {
    got <- run(p)
    expect_null(attr(got, "rejected_reason"), info = sprintf("2^%d", p))
    expect_equal(as.numeric(got) / 2^p, as.numeric(ref), tolerance = 1e-10,
                 info = sprintf("2^%d", p))
  }
})

test_that("compare_models() tells different rows and responses apart in small units", {
  lay <- .kt_layer(n = 40, seed = 2)
  stub <- function(d, looic) {
    f <- .kt_fit(d)
    f$info$looic <- looic
    class(f) <- c("kt_lm_fit", "bayesian_fit", "spatial_fit")
    f
  }
  for (p in c(0, -20, -40, -60)) {
    s <- 2^p
    # No CRS: the locations are compared in their own unit.  One row moved
    # by 1000 units times the scale is a different row set.
    a <- .fit_rowset(.kt_fit(.kt_sf(lay$xy * s, lay$y, lay$w, crs = NA)))
    xy2 <- lay$xy; xy2[1, 1] <- xy2[1, 1] + 1000
    b <- .fit_rowset(.kt_fit(.kt_sf(xy2 * s, lay$y, lay$w, crs = NA)))
    expect_false(.same_rowset_xy(a, b), info = sprintf("2^%d", p))
    expect_false(.same_rowset(a, b), info = sprintf("2^%d", p))
    # The same rows in another order are the same rows.
    o <- rev(seq_len(nrow(lay$xy)))
    a2 <- .fit_rowset(.kt_fit(.kt_sf(lay$xy[o, ] * s, lay$y[o], lay$w[o], crs = NA)))
    expect_true(.same_rowset(a, a2), info = sprintf("2^%d", p))
  }
  # A response in small units and the same response times 1.5 are different
  # responses.
  for (s in c(1, 2^-40, 1e-12, 1e-200)) {
    a <- .fit_rowset(.kt_fit(.kt_sf(lay$xy, lay$y * s, lay$w)))
    b <- .fit_rowset(.kt_fit(.kt_sf(lay$xy, lay$y * s * 1.5, lay$w)))
    expect_true(.same_rowset_xy(a, b))
    expect_false(.same_rowset(a, b), info = sprintf("response * %g", s))
    expect_true(.same_rowset(a, a))
  }
  fa <- stub(.kt_sf(lay$xy, lay$y * 1e-12, lay$w), 54.8)
  fb <- stub(.kt_sf(lay$xy, lay$y * 1.5e-12, lay$w), 60.1)
  w <- character(0)
  res <- withCallingHandlers(compare_models(list(raw = fa, scaled = fb)),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
  expect_true(all(is.na(res$LOOIC)))
  expect_true(any(grepl("the values of 'y' differ between raw, scaled", w)))
})

# Sites in metres, each with an exact twin and a twin two ulps away in
# northing.  At k = 1 the near twin ties with the exact one only through the
# absolute part of the tolerance, 1e-9 of the coordinates' size.
.kt_twins <- function(n = 30) {
  lay <- .kt_layer()
  x <- 5e5 + lay$xy[seq_len(n), 1]
  y <- 4.5e6 + lay$xy[seq_len(n), 2]
  ulp <- 2^(floor(log2(y)) - 52)
  cbind(c(x, x, x), c(y, y, y + 2 * ulp))
}

test_that("near twins in metres share a slot as they did, and at every scale", {
  tw <- .kt_twins()
  n <- nrow(tw) %/% 3L
  for (k in c(1L, 2L, 4L)) {
    old <- .kt_old_rule(tw, k)
    ref <- as.matrix(.build_knn_weights(tw, k = k))
    expect_identical(ref, old, info = sprintf("twins, k = %d", k))
    expect_identical(.build_knn_weights(tw, k = k, use_fnn = FALSE, use_matrix = FALSE),
                     old, info = sprintf("twins, dense, k = %d", k))
    for (p in c(-50, -60)) {
      expect_identical(as.matrix(.build_knn_weights(tw * 2^p, k = k)), ref,
                       info = sprintf("twins * 2^%d, k = %d", p, k))
      expect_identical(.build_knn_weights(tw * 2^p, k = k, use_fnn = FALSE,
                                          use_matrix = FALSE),
                       ref, info = sprintf("twins * 2^%d, dense, k = %d", p, k))
    }
  }
  # Row 1 at k = 1: its exact twin and its near twin share the one slot.
  W1 <- as.matrix(.build_knn_weights(tw, k = 1))
  expect_identical(which(W1[1, ] > 0), c(n + 1L, 2L * n + 1L))
  expect_identical(W1[1, c(n + 1L, 2L * n + 1L)], c(0.5, 0.5))
})

test_that("rows and responses that differ in their last bit are the same in small units", {
  lay <- .kt_layer(n = 40, seed = 2)
  s <- 2^-40
  a <- .fit_rowset(.kt_fit(.kt_sf(lay$xy * s, lay$y, lay$w, crs = NA)))
  xy2 <- lay$xy * s
  xy2[1, 1] <- xy2[1, 1] + 2^(floor(log2(xy2[1, 1])) - 52)
  expect_false(identical(xy2[1, 1], lay$xy[1, 1] * s))
  b <- .fit_rowset(.kt_fit(.kt_sf(xy2, lay$y, lay$w, crs = NA)))
  expect_true(.same_rowset_xy(a, b))
  expect_true(.same_rowset(a, b))
  ya <- lay$y * 1e-12
  yb <- ya * (1 + 2e-16)
  expect_false(identical(ya, yb))
  fa <- .fit_rowset(.kt_fit(.kt_sf(lay$xy, ya, lay$w)))
  fb <- .fit_rowset(.kt_fit(.kt_sf(lay$xy, yb, lay$w)))
  expect_true(.same_rowset(fa, fb))
})

test_that("folds built on other data are refused on a layer with no CRS in small units", {
  lay <- .kt_layer()
  n <- nrow(lay$xy)
  o <- rev(seq_len(n))
  for (p in c(0, -20, -45, -60)) {
    s <- 2^p
    d <- .kt_sf(lay$xy * s, lay$y, lay$w, crs = NA)
    f <- suppressWarnings(suppressMessages(
      make_folds(d, k = 5, method = "random_kfold", seed = 1)))
    expect_false(is.null(f$params$row_probe), info = sprintf("2^%d", p))
    d$..row_id <- seq_len(n)
    # Every row at another row's location.
    d2 <- .kt_sf(lay$xy[o, ] * s, lay$y, lay$w, crs = NA)
    d2$..row_id <- seq_len(n)
    expect_error(.check_fold_probe(f, d2, "cv_probe"),
                 "built from different data", info = sprintf("2^%d", p))
    expect_error(suppressMessages(.check_fold_provenance(f, d2, "fold_separation")),
                 "built from different data", info = sprintf("2^%d", p))
    # The folds' own layer passes, and so does it with every row moved by
    # one ulp in x.
    expect_silent(.check_fold_probe(f, d, "cv_probe"))
    xy3 <- lay$xy * s
    xy3[, 1] <- xy3[, 1] + 2^(floor(log2(xy3[, 1])) - 52)
    d3 <- .kt_sf(xy3, lay$y, lay$w, crs = NA)
    d3$..row_id <- seq_len(n)
    expect_silent(.check_fold_probe(f, d3, "cv_probe"))
  }
})

test_that("the edge margin is one part in 1e10 of the coordinates in any unit", {
  # Centres 1000 units apart, seeds placed 10 units off them by a fit of
  # their own, and seeds the centres moved by rounding.
  cen <- as.matrix(expand.grid(x = 0:4 * 1000, y = 0:4 * 1000))
  off <- cen + 10
  near <- cen * (1 + 4e-16)
  for (p in .kt_powers) {
    s <- 2^p
    expect_identical(.edge_tol(c(range(cen[, 1L]), range(cen[, 2L])) * s),
                     (1e-10 * 4000) * s, info = sprintf("2^%d", p))
    for (fnn in unique(c(FALSE, .has_backend("FNN")))) {
      expect_null(.seeds_as_centres(off * s, cen * s, use_fnn = fnn),
                  info = sprintf("seeds 10 units off, * 2^%d", p))
      expect_identical(.seeds_as_centres(near * s, cen * s, use_fnn = fnn),
                       seq_len(nrow(cen)), info = sprintf("seeds at the centres, * 2^%d", p))
    }
    # Seeds 10 units off stay where they are; the cells are drawn on them.
    rec <- list(scored_cells = TRUE, scored = list(centres = cen * s),
                crs = sf::st_crs(NA), cells = nrow(cen))
    pts <- sf::st_as_sf(data.frame(x = off[, 1L] * s, y = off[, 2L] * s),
                        coords = c("x", "y"), crs = NA)
    expect_identical(unname(sf::st_coordinates(.seeds_at_centres(pts, rec))),
                     unname(off * s), info = sprintf("seeds kept, * 2^%d", p))
  }
  # Coordinates that reach 1 unit keep the margin they had.
  for (bb in list(c(4e5, 4.4e6, 4.1e5, 4.41e6), c(-1, -1, 1, 1), c(0, 0, 0.5, 2))) {
    expect_identical(.edge_tol(bb), 1e-10 * max(1, abs(bb)))
  }
  expect_identical(.edge_tol(c(0, 0, 0, 0)), 0)
  expect_identical(.edge_tol(c(NA_real_, NA_real_, NA_real_, NA_real_)), 0)
})
