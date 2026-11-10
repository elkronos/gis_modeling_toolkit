# The k-nearest-neighbour weights, and every Moran's I read from them, do not
# move with the origin of the coordinates.
#
# .knn_pairs() tied two distances when they differed by less than
# 1e-9 * dk + 1e-9 * max|coords|, and the kd-tree shortcut of
# .sac_structure_test() took the same absolute term.  That term is a part of
# the coordinates' size, not of their rounding: 4 mm at a UTM northing, 1000 m
# at an offset of 1e12 and 1e6 m at 1e15.  So the same points moved by a
# whole number of map units got 18 to 79 neighbours instead of 8 at 1e12, and
# at 1e15 every point was every other point's neighbour, Moran's I was
# exactly -1 / (n - 1) and estimate_sac_range() reported "no spatial
# structure detected".  Distances whose squares lie within .TIE_TOL of each
# other are now a tie, with no term tied to the coordinates' size.
#
# Every shift below is exact: the coordinates are multiples of 1/8 m and an
# offset of 1e15 holds them to 1/8 m, so the points, their distances and the
# response are the same at every offset, and anything that moves is the
# package's.

# 60 points with coordinates on multiples of 1/8 m in a 4000 x 3000 m box of
# UTM zone 17N, a smooth field and two predictors.
.ko_points <- function(n = 60L) {
  withr::with_seed(20261106L, {
    px <- round(8 * (600000 + stats::runif(n, 0, 4000))) / 8
    py <- round(8 * (3990000 + stats::runif(n, 0, 3000))) / 8
    ux <- (px - 600000) / 1000; uy <- (py - 3990000) / 1000
    x1 <- 0.8 * ux + stats::rnorm(n, 0, 0.3)
    x2 <- stats::runif(n, -1, 1)
    y  <- 2 + 1.5 * sin(1.2 * ux) + cos(1.4 * uy) + 0.6 * x1 - 0.4 * x2 +
      stats::rnorm(n, 0, 0.2)
  })
  stopifnot(!anyDuplicated(cbind(px, py)))
  list(xy = cbind(px, py), att = data.frame(y = y, x1 = x1, x2 = x2))
}

# An 8 x 8 lattice at 250 m: many distances tie exactly.
.ko_lattice <- function() {
  g <- expand.grid(i = 0:7, j = 0:7)
  list(xy = cbind(600000 + 250 * g$i, 3990000 + 250 * g$j),
       att = data.frame(y = 2 + sin(0.9 * g$i) + cos(0.7 * g$j),
                        x1 = g$i / 7 + 0.1 * cos(g$j),
                        x2 = ((3 * g$i + 5 * g$j) %% 7) / 7 - 0.5))
}

.ko_sf <- function(L, xy, crs = 32617) {
  sf::st_as_sf(cbind(L$att, X = xy[, 1], Y = xy[, 2]),
               coords = c("X", "Y"), crs = crs)
}

# The same points moved by `off` in x and y; the shift is checked exact.
.ko_shift <- function(xy, off) {
  s <- xy + off
  stopifnot(all(s[, 1] - off == xy[, 1]), all(s[, 2] - off == xy[, 2]))
  s
}

# The pairs in a canonical order, so two runs compare as sets.
.ko_pairs <- function(xy, k, fnn) {
  p <- .knn_pairs(xy, k, use_fnn = fnn)
  o <- order(p$i, p$j)
  list(i = p$i[o], j = p$j[o], w = p$w[o])
}

.ko_offsets <- c(1e8, 1e12, 1e15)

.ko_fit <- function(d) {
  new_spatial_fit("ko_lm_fit", stats::lm(y ~ x1 + x2, data = sf::st_drop_geometry(d)),
                  y ~ x1 + x2, "y", c("x1", "x2"), d)
}
residuals.ko_lm_fit <- function(object, ...) as.numeric(stats::residuals(object$engine))
fitted.ko_lm_fit <- function(object, ...) as.numeric(stats::fitted(object$engine))
registerS3method("residuals", "ko_lm_fit", residuals.ko_lm_fit)
registerS3method("fitted", "ko_lm_fit", fitted.ko_lm_fit)

test_that(".knn_pairs() gives the same pairs and weights at every offset", {
  for (nm in c("points", "lattice")) {
    L <- if (nm == "points") .ko_points() else .ko_lattice()
    for (fnn in c(TRUE, FALSE)) {
      if (fnn && !requireNamespace("FNN", quietly = TRUE)) next
      for (k in c(2L, 4L, 8L, 16L)) {
        ref <- .ko_pairs(L$xy, k, fnn)
        for (off in .ko_offsets) {
          got <- .ko_pairs(.ko_shift(L$xy, off), k, fnn)
          lbl <- sprintf("%s, k = %d, fnn = %s, offset %g", nm, k, fnn, off)
          expect_identical(got, ref, info = lbl)
          expect_identical(max(tabulate(got$i, nrow(L$xy))),
                           max(tabulate(ref$i, nrow(L$xy))), info = lbl)
        }
      }
    }
  }
})

test_that("on the points each row keeps its k neighbours at 1e12 and 1e15", {
  L <- .ko_points()
  for (off in c(0, .ko_offsets)) {
    p <- .knn_pairs(.ko_shift(L$xy, off), 8L, use_fnn = FALSE)
    expect_identical(range(tabulate(p$i, nrow(L$xy))), c(8L, 8L),
                     info = sprintf("offset %g", off))
  }
})

test_that("exact ties on the lattice keep their shared weights at every offset", {
  L <- .ko_lattice()
  inner <- which(L$xy[, 1] == 600000 + 750 & L$xy[, 2] == 3990000 + 750)
  for (off in c(0, .ko_offsets)) {
    W <- as.matrix(.build_knn_weights(.ko_shift(L$xy, off), k = 6))
    # Four at 250 m take 1/6 each; the four diagonals tie for the last two
    # slots and take 1/12 each.
    expect_identical(sort(W[inner, W[inner, ] > 0]), c(rep(1 / 12, 4), rep(1 / 6, 4)),
                     info = sprintf("offset %g", off))
  }
})

test_that(".knn_pairs() gives the same pairs in units of 1e-6 m and US survey feet", {
  # Rescaled coordinates round apart in their last bits; ties within
  # .TIE_TOL hold the lattice's and nothing else changes.
  for (nm in c("points", "lattice")) {
    L <- if (nm == "points") .ko_points() else .ko_lattice()
    for (fnn in c(TRUE, FALSE)) {
      if (fnn && !requireNamespace("FNN", quietly = TRUE)) next
      for (k in c(4L, 8L)) {
        ref <- .ko_pairs(L$xy, k, fnn)
        for (f in c(1e-6, 1200 / 3937)) {
          got <- .ko_pairs(L$xy / f, k, fnn)
          expect_identical(got, ref, info = sprintf("%s, k = %d, fnn = %s, unit %g",
                                                    nm, k, fnn, f))
        }
      }
    }
  }
})

test_that("the structure test's table is the same at every offset", {
  for (nm in c("points", "lattice")) {
    L <- if (nm == "points") .ko_points() else .ko_lattice()
    X <- cbind(1, L$att$x1, L$att$x2)
    ref <- .sac_structure_test(L$xy, L$att$y, X)
    for (off in .ko_offsets) {
      got <- .sac_structure_test(.ko_shift(L$xy, off), L$att$y, X)
      # The weights are identical; I, z and p differ at most by the order in
      # which the same products are summed.
      expect_equal(got$table, ref$table, tolerance = 1e-10,
                   info = sprintf("%s, offset %g", nm, off))
      expect_equal(got$p, ref$p, tolerance = 1e-10, info = sprintf("%s, offset %g", nm, off))
    }
  }
})

test_that("residual_morans_i() is the same at 0, 1e12 and 1e15", {
  L <- .ko_points()
  ref <- suppressWarnings(residual_morans_i(.ko_fit(.ko_sf(L, L$xy))))
  expect_true(ref$observed > 0.05)
  for (off in c(1e12, 1e15)) {
    got <- suppressWarnings(residual_morans_i(.ko_fit(.ko_sf(L, .ko_shift(L$xy, off)))))
    # The fit reads no coordinates, so the residuals are the same; the
    # statistic sums the same products, to rounding.
    for (nm in c("observed", "expected", "sd", "z", "p_value"))
      expect_equal(got[[nm]], ref[[nm]], tolerance = 1e-10,
                   info = sprintf("%s at offset %g", nm, off))
  }
})

test_that("estimate_sac_range() reads the same structure at 0 and 1e12", {
  L <- .ko_points()
  run <- function(xy) suppressWarnings(suppressMessages(
    estimate_sac_range(.ko_sf(L, xy), "y")))
  ref <- run(L$xy)
  got <- run(.ko_shift(L$xy, 1e12))
  expect_true(attr(ref, "structure_p") < 0.05)
  expect_equal(attr(got, "morans_i"), attr(ref, "morans_i"), tolerance = 1e-10)
  expect_equal(attr(got, "structure_p"), attr(ref, "structure_p"), tolerance = 1e-10)
  expect_false(identical(attr(got, "rejected_reason"),
                         "no spatial structure detected at the nearest-neighbour scale"))
})

test_that("determine_optimal_levels(criterion = \"combined\") and resolution_profile() agree at 0 and 1e12", {
  L <- .ko_points()
  dol <- function(xy) suppressWarnings(suppressMessages(
    determine_optimal_levels(.ko_sf(L, xy), max_levels = 59L, response_var = "y",
                             predictor_vars = c("x1", "x2"), criterion = "combined",
                             on_no_elbow = "warn")))
  ref <- dol(L$xy)
  got <- dol(.ko_shift(L$xy, 1e12))
  expect_identical(as.vector(got), as.vector(ref))
  dr <- attr(ref, "diagnostics"); dg <- attr(got, "diagnostics")
  expect_true(any(is.finite(dr$moran_i)))
  expect_equal(dg$moran_i, dr$moran_i, tolerance = 1e-10)
  expect_equal(dg$moran_z, dr$moran_z, tolerance = 1e-10)
  expect_identical(dg$combined_rank, dr$combined_rank)
  prof <- function(xy) suppressWarnings(suppressMessages(
    resolution_profile(.ko_sf(L, xy), response_var = "y", levels = c(3L, 6L, 12L),
                       estimate = "pooled")))
  pr <- prof(L$xy)
  pg <- prof(.ko_shift(L$xy, 1e12))
  expect_true(any(is.finite(pr$moran_i)))
  expect_equal(pg$moran_i, pr$moran_i, tolerance = 1e-10)
})

test_that("the structure test's Moran's I is the one of the tie-sharing weights in every unit", {
  # The kd-tree shortcut of .sac_structure_test() gives 1 / k to each of the
  # k nearest and is taken only when no row's (k+1)-th distance ties with its
  # k-th.  Where it is taken wrongly, a lattice in another unit reads other
  # graphs than .build_knn_weights().  Moran's I of the same weights sums the
  # same products in another order, hence the tolerance.
  g <- expand.grid(i = 0:7, j = 0:7)
  y <- sin(g$i) + cos(1.3 * g$j)
  e <- y - mean(y)
  for (u in c(1, 1e-6, 1200 / 3937, 0.3)) {
    xy <- cbind(600000 + 250 * g$i, 3990000 + 250 * g$j) / u
    st <- .sac_structure_test(xy, y)
    want <- vapply(c(2L, 4L, 8L, 16L), function(k) {
      W <- as.matrix(.build_knn_weights(xy, k = k))
      sum(e * (W %*% e)) / sum(e^2)
    }, numeric(1))
    expect_equal(st$table$morans_i, want, tolerance = 1e-12,
                 info = sprintf("unit %g", u))
  }
})

test_that("two distances tie when their squares are within .TIE_TOL, and only then", {
  # Row 1 at the origin with two neighbours 100 and 100 * f away and one far
  # point.  At k = 1 the slot goes to the nearer alone, or to both at 0.5
  # when the two distances tie.  (1 + 4e-8)^2 is within 1e-7 of 1;
  # (1 + 6e-8)^2 is not.
  for (fnn in c(TRUE, FALSE)) {
    if (fnn && !requireNamespace("FNN", quietly = TRUE)) next
    for (f in c(1 + 4e-8, 1 + 6e-8)) {
      xy <- rbind(c(0, 0), c(100, 0), c(0, 100 * f), c(5000, 5000))
      p <- .knn_pairs(xy, 1L, use_fnn = fnn)
      r1 <- p$i == 1L
      got <- p$w[r1][order(p$j[r1])]
      lbl <- sprintf("factor 1 + %g, fnn = %s", f - 1, fnn)
      if (f < 1 + 5e-8) {
        expect_identical(sort(p$j[r1]), c(2L, 3L), info = lbl)
        expect_identical(got, c(0.5, 0.5), info = lbl)
      } else {
        expect_identical(p$j[r1], 2L, info = lbl)
        expect_identical(got, 1, info = lbl)
      }
    }
  }
})
