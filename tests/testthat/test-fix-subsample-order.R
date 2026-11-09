# Subsamples and thinning drawn from the rows in a canonical order
# (.canonical_order()), so a result does not depend on the row order of an
# input (BUG-116 to BUG-119).  Each drew by row position from the rows as
# given: the variogram design effect of summarize_by_cell() above
# deff_max_n, the DI normaliser of area_of_applicability() above
# normalizer_max_n, the distances of gp_lengthscale_bounds() above max_n,
# and Moran's I at the GP basis floor of the pre-fit check above max_n.
# Under each cap nothing is drawn, and the results are those the package
# gave before (the figures below were taken on the code before the change).

.fso_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# A smooth field on n uniform points, with a noise column w.
.fso_layer <- function(n, side, a, seed, noise = 0.3) {
  withr::with_seed(seed, {
    x <- runif(n, 0, side); y <- runif(n, 0, side); w <- rnorm(n)
    D <- as.matrix(stats::dist(cbind(x, y)))
    f <- as.numeric(t(chol(exp(-(D / a)^2) + diag(1e-6, n))) %*% rnorm(n))
    z <- f + noise * w
  })
  data.frame(x = x, y = y, w = w, z = z)
}

# The given order, the reverse and random orders.
.fso_orders <- function(n, k = 1L) {
  c(list(seq_len(n), rev(seq_len(n))),
    lapply(seq_len(k), function(s) withr::with_seed(100 + s, sample.int(n))))
}

test_that("the variogram design effect of a cell above deff_max_n does not depend on the row order (BUG-116)", {
  d <- .fso_layer(1500, 5000, 800, seed = 7)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  pts$cell_id <- 1L
  run <- function(o, ...) sf::st_drop_geometry(.fso_quiet(summarize_by_cell(
    pts[o, ], "z", id_col = "cell_id", deff = "variogram", estimate = "mean",
    conf_level = 0.9, ...)))
  os  <- .fso_orders(1500, 1L)
  ref <- run(os[[1L]])
  cols <- c("..se_resp_z", "..neff_resp_z", "..df_resp_z", "..ci_lo_resp_z",
            "..ci_hi_resp_z", "cell_weight")
  expect_true(all(cols %in% names(ref)))
  for (o in os[-1L]) {
    r <- run(o)
    for (cc in cols) expect_equal(r[[cc]], ref[[cc]], tolerance = 1e-12, info = cc)
  }
  # The caller's stream is where it was.
  set.seed(11); a <- runif(1)
  set.seed(11); invisible(run(os[[3L]])); b <- runif(1)
  expect_identical(a, b)
})

test_that("summarize_by_cell(deff = 'variogram') gives what it gave on cells under deff_max_n (BUG-116)", {
  d <- .fso_layer(300, 3000, 600, seed = 3)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  pts$cell_id <- rep(1:2, each = 150)
  r <- sf::st_drop_geometry(.fso_quiet(summarize_by_cell(
    pts, "z", id_col = "cell_id", deff = "variogram", estimate = "mean")))
  r <- r[order(r$cell_id), ]
  # The field is drawn through chol() of a dense covariance matrix
  # (.fso_layer()), whose last bits follow the BLAS and LAPACK in use, and
  # the variogram fit carries them into the design effect: Linux x86-64 and
  # macOS arm64 differ by up to 1.4e-11 relative here.  So the figures are
  # compared within 1e-9 relative.
  expect_equal(r$..se_resp_z, c(0.28584766841440995, 0.27506304572341911),
               tolerance = 1e-9)
  expect_equal(r$cell_weight, c(8.8204054172786872, 8.8864136640242126),
               tolerance = 1e-9)
})

test_that("the DI normaliser above normalizer_max_n does not depend on the order of the training rows (BUG-117)", {
  d   <- .fso_layer(1500, 5000, 800, seed = 7)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  tr  <- pts[1:1000, ]; nd <- pts[1001:1500, ]
  run <- function(o, coords) .fso_quiet(area_of_applicability(
    nd, train_sf = tr[o, ], predictor_vars = c("w", "z"), normalizer_max_n = 300,
    coords = coords))
  os <- .fso_orders(1000, 1L)
  for (coords in c("no", "yes")) {
    ref <- run(os[[1L]], coords)
    for (o in os[-1L]) {
      r <- run(o, coords)
      expect_equal(r$normalizer, ref$normalizer, tolerance = 1e-12)
      expect_equal(r$threshold, ref$threshold, tolerance = 1e-12)
      expect_equal(r$aoa$DI, ref$aoa$DI, tolerance = 1e-12)
      expect_identical(r$aoa$AOA, ref$aoa$AOA)
      if (coords == "yes") {
        expect_equal(r$location$normalizer, ref$location$normalizer, tolerance = 1e-12)
        expect_equal(r$location$threshold, ref$location$threshold, tolerance = 1e-12)
        expect_equal(r$aoa$DI_location, ref$aoa$DI_location, tolerance = 1e-12)
      }
    }
  }
})

test_that("area_of_applicability() gives what it gave under normalizer_max_n (BUG-117)", {
  d   <- .fso_layer(300, 3000, 600, seed = 3)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  r <- .fso_quiet(area_of_applicability(pts[201:300, ], train_sf = pts[1:200, ],
                                        predictor_vars = c("w", "z"), coords = "yes"))
  expect_equal(r$normalizer, 1.759075277957445, tolerance = 1e-12)
  expect_equal(r$threshold, 0.19877958966738746, tolerance = 1e-12)
  expect_equal(r$location$normalizer, 1549.3091158202399, tolerance = 1e-12)
})

test_that("gp_lengthscale_bounds() above max_n does not depend on the row order (BUG-118)", {
  d  <- .fso_layer(1500, 5000, 800, seed = 7)
  xy <- as.matrix(d[, c("x", "y")])
  os <- .fso_orders(1500, 2L)
  ref <- gp_lengthscale_bounds(xy)
  for (o in os[-1L]) expect_identical(gp_lengthscale_bounds(xy[o, ]), ref)
  # Repeat visits are collapsed first, so they move nothing either.
  rep_xy <- rbind(xy, xy[1:300, ])
  expect_identical(gp_lengthscale_bounds(rep_xy[rev(seq_len(nrow(rep_xy))), ]), ref)
  # The caller's stream is where it was.
  set.seed(11); a <- runif(1)
  set.seed(11); invisible(gp_lengthscale_bounds(xy[os[[3L]], ])); b <- runif(1)
  expect_identical(a, b)
  # With every location the bounds are those of the subsample's population.
  all_b <- gp_lengthscale_bounds(xy, max_n = Inf)
  expect_identical(gp_lengthscale_bounds(xy[os[[2L]], ], max_n = Inf), all_b)
})

test_that("gp_lengthscale_bounds() gives what it gave under max_n (BUG-118)", {
  d  <- .fso_layer(300, 3000, 600, seed = 3)
  xy <- as.matrix(d[, c("x", "y")])
  expect_equal(gp_lengthscale_bounds(xy),
               c(lower = 398.97892886889508, upper = 1643.6464514310533),
               tolerance = 1e-14)
})

# The layer of the item: 2400 uniform points in a 6 km square, a field with
# a scale of 600 m, a nugget and noise; w is a predictor of noise alone.
# Thinned through the rows as given, its I at a floor of 1140 m fell either
# side of .gp_range_correlation (0.031 to 0.061 over 10 orders).
.fso_floor_layer <- function(n = 2400) {
  withr::with_seed(8, { x <- runif(n, 0, 6000); y <- runif(n, 0, 6000) })
  D <- as.matrix(stats::dist(cbind(x, y)))
  z <- withr::with_seed(5, {
    f <- as.numeric(t(chol(exp(-(D / 600)^2) + diag(1e-6, n))) %*% rnorm(n)) +
      rnorm(n, sd = sqrt(0.05))
    f + rnorm(n, sd = 0.3)
  })
  data.frame(x = x, y = y, w = withr::with_seed(9, rnorm(n)), z = z)
}

test_that("Moran's I at the basis floor above max_n does not depend on the row order (BUG-119)", {
  d  <- .fso_floor_layer()
  xy <- as.matrix(d[, c("x", "y")])
  X  <- cbind(w = d$w)
  os <- .fso_orders(2400, 3L)
  for (fl in c(400, 1140)) {
    for (Xo in list(NULL, X)) {
      ref <- spatialkit:::.gp_floor_moran(xy, d$z, Xo, floor = fl)
      expect_false(is.null(ref))
      for (o in os[-1L]) {
        r <- spatialkit:::.gp_floor_moran(xy[o, ], d$z[o],
                                          if (is.null(Xo)) NULL else Xo[o, , drop = FALSE],
                                          floor = fl)
        expect_identical(r[c("I", "z", "p", "n_pairs")], ref[c("I", "z", "p", "n_pairs")])
      }
    }
  }
})

test_that("the pre-fit status decided at the floor does not depend on the row order (BUG-119)", {
  d   <- .fso_floor_layer()
  lay <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  # A range below the floor, so the status is decided by the floor's test.
  testthat::local_mocked_bindings(
    estimate_sac_range = function(...) 300, .package = "spatialkit")
  run <- function(o) .fso_quiet(spatialkit:::.gp_prefit_range_check(
    lay[o, ], "z", character(0), resolution = c(x = 1140, y = 1140),
    gp_c = 1.5, S = 1, xy_scale = 1))
  # Ten orders: as given, each of the item's statuses came up for some.
  os  <- .fso_orders(2400, 8L)
  ref <- run(os[[1L]])
  expect_false(is.null(ref$moran_floor))
  for (o in os[-1L]) {
    r <- run(o)
    expect_identical(r$status, ref$status)
    expect_identical(r$moran_floor[c("I", "z", "p", "n_pairs")],
                     ref$moran_floor[c("I", "z", "p", "n_pairs")])
  }
})

test_that("Moran's I at the basis floor is what it was under max_n (BUG-119)", {
  d  <- .fso_layer(600, 3000, 300, seed = 4)
  xy <- as.matrix(d[, c("x", "y")])
  r  <- spatialkit:::.gp_floor_moran(xy, d$z, NULL, floor = 400)
  expect_equal(c(r$I, r$z, r$n_pairs),
               c(0.16499522067636652, 15.714371587441171, 8407),
               tolerance = 1e-12)
})

# Inputs on which one key does not decide the order: a regular grid (ties in
# x and in y), integer and indicator-coded predictors, and rows at one
# location.  Each order below is a random one, since reversing a grid can
# leave the order of a tie rule unchanged.
.fso_grid <- function(nx, ny, step = 100) {
  as.matrix(expand.grid(x = (seq_len(nx) - 1) * step, y = (seq_len(ny) - 1) * step))
}
.fso_field <- function(xy, a, seed) {
  withr::with_seed(seed, {
    D <- as.matrix(stats::dist(xy))
    as.numeric(t(chol(exp(-(D / a)^2) + diag(1e-6, nrow(xy)))) %*% rnorm(nrow(xy)))
  })
}
.fso_random_orders <- function(n, k = 4L, from = 200L) {
  c(list(seq_len(n)),
    lapply(seq_len(k), function(s) withr::with_seed(from + s, sample.int(n))))
}

test_that("a gridded cell above deff_max_n gives one design effect in every row order (BUG-116)", {
  g <- .fso_grid(30, 30)
  d <- data.frame(g, z = .fso_field(g, 300, 10) +
                    withr::with_seed(11, rnorm(nrow(g), sd = 0.3)), cell_id = 1L)
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  run <- function(o) sf::st_drop_geometry(.fso_quiet(summarize_by_cell(
    pts[o, ], "z", id_col = "cell_id", deff = "variogram", estimate = "mean",
    deff_max_n = 200)))
  os  <- .fso_random_orders(nrow(g))
  ref <- run(os[[1L]])
  expect_true(is.finite(ref$..se_resp_z))
  for (o in os[-1L]) {
    r <- run(o)
    expect_equal(r$..se_resp_z, ref$..se_resp_z, tolerance = 1e-12)
    expect_equal(r$cell_weight, ref$cell_weight, tolerance = 1e-12)
  }
})

test_that("the DI normaliser above normalizer_max_n is one with tied or indicator-coded predictors and on a gridded layer (BUG-117)", {
  n <- 1500
  d <- withr::with_seed(3, {
    d <- data.frame(x = runif(n, 0, 5000), y = runif(n, 0, 5000))
    d$k <- sample(1:5, n, TRUE)
    d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
    d$w <- round(rnorm(n), 1)
    d
  })
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  tr  <- pts[1:1000, ]; nd <- pts[1001:1500, ]
  os  <- .fso_random_orders(1000)
  # The first predictor column has ties, so the next one must decide.
  for (pv in list(c("k", "w"), c("k", "f"))) {
    run <- function(o) .fso_quiet(area_of_applicability(
      nd, train_sf = tr[o, ], predictor_vars = pv, normalizer_max_n = 300,
      coords = "no", categorical = "indicator"))
    ref <- run(os[[1L]])
    for (o in os[-1L]) {
      r <- run(o)
      expect_equal(r$normalizer, ref$normalizer, tolerance = 1e-12,
                   info = paste(pv, collapse = "+"))
      expect_equal(r$aoa$DI, ref$aoa$DI, tolerance = 1e-12,
                   info = paste(pv, collapse = "+"))
    }
  }
  # A gridded training layer: the location normaliser.
  g  <- data.frame(.fso_grid(40, 25))
  g$w <- withr::with_seed(4, rnorm(nrow(g)))
  gs <- sf::st_as_sf(g, coords = c("x", "y"), crs = 32617)
  os <- .fso_random_orders(nrow(gs))
  run <- function(o) .fso_quiet(area_of_applicability(
    gs[1:50, ], train_sf = gs[o, ], predictor_vars = "w",
    normalizer_max_n = 300, coords = "yes"))
  ref <- run(os[[1L]])
  for (o in os[-1L]) {
    r <- run(o)
    expect_equal(r$location$normalizer, ref$location$normalizer, tolerance = 1e-12)
    expect_equal(r$aoa$DI_location, ref$aoa$DI_location, tolerance = 1e-12)
  }
})

test_that("with a threshold supplied, the inside or outside call above normalizer_max_n does not depend on the row order (BUG-117)", {
  d <- withr::with_seed(7, {
    data.frame(x = runif(1500, 0, 5000), y = runif(1500, 0, 5000), w = rnorm(1500))
  })
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  tr  <- pts[1:1000, ]; nd <- pts[1001:1500, ]
  run <- function(o) .fso_quiet(area_of_applicability(
    nd, train_sf = tr[o, ], predictor_vars = "w", normalizer_max_n = 300,
    threshold = 0.0065))
  ref <- run(seq_len(1000))
  for (o in .fso_orders(1000, 2L)[-1L])
    expect_identical(run(o)$aoa$AOA, ref$aoa$AOA)
})

test_that("gp_lengthscale_bounds() on a grid above max_n gives one answer in every row order (BUG-118)", {
  g   <- .fso_grid(40, 40)
  os  <- .fso_random_orders(nrow(g))
  ref <- gp_lengthscale_bounds(g)
  for (o in os[-1L]) expect_identical(gp_lengthscale_bounds(g[o, ]), ref)
})

test_that("Moran's I at the basis floor above max_n is one for rows at one location and on a grid (BUG-119)", {
  s  <- withr::with_seed(5, cbind(runif(1200, 0, 6000), runif(1200, 0, 6000)))
  f  <- .fso_field(s, 600, 6)
  xy <- rbind(s, s)
  X  <- cbind(w = withr::with_seed(9, rnorm(2400)))
  os <- .fso_random_orders(2400)
  key <- c("I", "z", "p", "n_pairs")
  same_in_every_order <- function(xy, y, X) {
    ref <- spatialkit:::.gp_floor_moran(xy, y, X, floor = 1140)
    expect_false(is.null(ref))
    for (o in os[-1L]) {
      r <- spatialkit:::.gp_floor_moran(xy[o, ], y[o],
                                        if (is.null(X)) NULL else X[o, , drop = FALSE],
                                        floor = 1140)
      expect_identical(r[key], ref[key])
    }
  }
  # Each site visited twice, with a different response on each visit.
  y2 <- c(f, f) + withr::with_seed(12, rnorm(2400, sd = 0.3))
  same_in_every_order(xy, y2, X)
  # The same response on both visits, and a different predictor.
  same_in_every_order(xy, c(f, f), X)
  # A 50 x 48 grid.
  g  <- as.matrix(expand.grid(x = seq(0, 6000, length.out = 50),
                              y = seq(0, 6000, length.out = 48)))
  fg <- .fso_field(g, 600, 7) + withr::with_seed(8, rnorm(nrow(g), sd = 0.3))
  same_in_every_order(g, fg, NULL)
})
