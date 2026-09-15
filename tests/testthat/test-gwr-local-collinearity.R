# tests/testthat/test-gwr-local-collinearity.R
# ---------------------------------------------------------------------------
# The local collinearity survey fit_gwr_model() retains (10.2) and the
# coefficient map that masks on it (9.6).
#
# The survey is computed by the package itself, not by GWmodel, so it is
# tested directly against a by-hand weighted condition index; the map is
# tested on the fit the GWmodel backend (or its test stub) returns.
# ---------------------------------------------------------------------------

lc_points <- function(n = 200, seed = 3) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  a <- rnorm(n); b <- rnorm(n)
  sf::st_as_sf(data.frame(x = x, y = y, a = a, b = b,
                          z = (1 + 2 * x / 1000) * a - b + rnorm(n, 0, 0.3)),
               coords = c("x", "y"), crs = 32632)
}

# Four clusters; `soil` is a regional covariate that is all but constant inside
# each cluster (0.5 in two of them, 1.5 in the other two, plus noise of SD
# 0.01), so a window inside one cluster is near-collinear with the intercept:
# scaled condition index 87-735 at 20 neighbours, against 4-8 for a window
# spanning both levels.  It is deliberately NOT exactly constant: GWmodel
# inverts X'WX with Armadillo's inv(), which throws "matrix is singular" on an
# exactly singular window and aborts the whole fit, so an exactly-constant
# covariate cannot exercise the survey on the real backend at all.  The
# smallest rcond of X'WX here is 1e-6, comfortably invertible.
lc_clusters <- function(seed = 5) {
  set.seed(seed)
  cl <- rep(1:4, each = 50)
  d <- sf::st_as_sf(data.frame(
    x = c(runif(50, 0, 100), runif(50, 400, 500), runif(50, 0, 100), runif(50, 400, 500)),
    y = c(runif(50, 0, 100), runif(50, 0, 100), runif(50, 400, 500), runif(50, 400, 500)),
    a = rnorm(200), soil = c(0.5, 1.5, 0.5, 1.5)[cl] + rnorm(200, 0, 0.01)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + 3 * d$soil + rnorm(200, 0, 0.2)
  d
}

test_that("the kernel weights are GWmodel's definitions", {
  w <- spatialkit:::.gw_kernel_weights
  d <- c(0, 50, 100, 150)
  expect_equal(w(d, 100, "boxcar", FALSE), c(1, 1, 0, 0))
  expect_equal(w(d, 100, "bisquare", FALSE), c(1, (1 - 0.25)^2, 0, 0))
  expect_equal(w(d, 100, "tricube", FALSE), c(1, (1 - 0.125)^3, 0, 0))
  expect_equal(w(d, 100, "gaussian", FALSE), exp(-0.5 * (d / 100)^2))
  expect_equal(w(d, 100, "exponential", FALSE), exp(-d / 100))
  # Adaptive: the distance parameter is the bw-th nearest distance.
  expect_equal(w(d, 3, "bisquare", TRUE), c(1, (1 - 0.25)^2, 0, 0))
  expect_equal(w(d, 4, "boxcar", TRUE), c(1, 1, 1, 0))
})

test_that("the survey is the weighted condition index at every location", {
  pts <- lc_points()
  xy <- sf::st_coordinates(pts)
  xm <- as.matrix(sf::st_drop_geometry(pts)[, c("a", "b")])
  sv <- spatialkit:::.gwr_local_collinearity(xy, xm, adaptive = TRUE, bw = 40, kernel = "bisquare")
  expect_equal(nrow(sv), nrow(pts))
  expect_equal(names(sv), c("row", "x", "y", "n_window", "cn"))
  expect_true(all(sv$n_window == 39L))       # the 40th neighbour has weight 0
  expect_true(all(is.finite(sv$cn) & sv$cn >= 1))
  for (i in c(1L, 77L, 200L)) {
    d <- sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
    h <- sort(d)[40]; w <- ifelse(d / h < 1, (1 - (d / h)^2)^2, 0); keep <- w > 1e-8
    expect_equal(sv$cn[i], spatialkit:::.condition_index(sqrt(w[keep]) * cbind(1, xm)[keep, ]))
  }
  # A window with fewer usable rows than columns is singular, not NA.
  expect_true(all(is.infinite(spatialkit:::.gwr_local_collinearity(xy[1:3, ], xm[1:3, ], TRUE, 2, "boxcar")$cn)))
  # Nothing to survey with one predictor.
  one <- spatialkit:::.gwr_local_collinearity(xy, xm[, 1, drop = FALSE], TRUE, 40, "bisquare")
  expect_true(all(is.na(one$cn)))
})

test_that("fit_gwr_model() keeps the survey and the global index on the fit", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  pts <- lc_points()
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(pts, "z", c("a", "b"), adaptive = TRUE, bandwidth = 60)))
  lc <- fit$info$local_collinearity
  expect_s3_class(lc, "data.frame")
  expect_equal(nrow(lc), nrow(pts))
  expect_true(all(is.finite(lc$cn)))
  expect_equal(fit$info$n_local_collinear, 0L)
  expect_equal(fit$info$n_local_singular, 0L)
  expect_true(is.finite(fit$info$condition_index))
  expect_equal(fit$info$condition_index,
               spatialkit:::.condition_index(cbind(1, as.matrix(sf::st_drop_geometry(pts)[, c("a", "b")]))))
  # It is the survey at the bandwidth actually used.
  ref <- spatialkit:::.gwr_local_collinearity(sf::st_coordinates(pts),
                                              as.matrix(sf::st_drop_geometry(pts)[, c("a", "b")]),
                                              TRUE, 60, "bisquare")
  expect_equal(lc$cn, ref$cn)
  # One predictor: nothing surveyed, NA index.
  fit1 <- suppressWarnings(suppressMessages(fit_gwr_model(pts, "z", "a", adaptive = TRUE, bandwidth = 60)))
  expect_null(fit1$info$local_collinearity)
  expect_true(is.na(fit1$info$condition_index))
  expect_true(is.na(fit1$info$n_local_collinear))
})

test_that("the collinearity warning is the exact fraction, and quiet when there is nothing", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- lc_clusters()
  # Every warning raised while evaluating `expr`, muffled; the value comes
  # back beside them, so a fit that fails is a test failure rather than a
  # missing object.
  .warns <- function(expr) {
    w <- character(0)
    val <- withCallingHandlers(suppressMessages(expr),
                               warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
    list(value = val, warnings = w)
  }
  r20 <- .warns(fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE, bandwidth = 20))
  expect_s3_class(r20$value, "gwr_fit")
  expect_true(any(grepl("100% of 200 locations have a collinear local design", r20$warnings)))
  expect_equal(r20$value$info$n_local_collinear, 200L)
  r199 <- .warns(fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE, bandwidth = 199))
  expect_s3_class(r199$value, "gwr_fit")
  expect_false(any(grepl("collinear local design", r199$warnings)))
  expect_equal(r199$value$info$n_local_collinear, 0L)
})

test_that("the coefficient map draws the local coefficient and masks collinear windows", {
  skip_if_not_installed("ggplot2"); skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  pts <- lc_points()
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(pts, "z", c("a", "b"), adaptive = TRUE, bandwidth = 60)))
  p <- plot(fit, type = "coefficients")
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_equal(p$labels$title, "Local coefficient of a")
  expect_match(p$labels$subtitle, "No location masked")
  expect_match(p$labels$caption, "Adaptive bandwidth 60, bisquare kernel")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_equal(sum(geoms == "GeomSf"), 1L)                      # nothing hollow
  expect_equal(nrow(b$data[[1L]]), nrow(pts))
  # The mapped values are coef()'s, and follow the west-east gradient.
  cf <- stats::coef(fit)
  expect_gt(stats::cor(cf$a, sf::st_coordinates(pts)[, 1]), 0.8)
  pb <- plot(fit, type = "coefficients", term = "b")
  expect_equal(pb$labels$title, "Local coefficient of b")
  expect_error(plot(fit, type = "coefficients", term = "nope"), "`term` must be one of")

  # Partly collinear: masked locations are drawn hollow and counted.  At 70
  # neighbours a window reaches past its own cluster, and whether the
  # neighbouring cluster carries the other `soil` level decides its index.
  d <- lc_clusters()
  fit70 <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE, bandwidth = 70)))
  expect_s3_class(fit70, "gwr_fit")
  n_bad <- fit70$info$n_local_collinear
  expect_gt(n_bad, 0L); expect_lt(n_bad, 200L)
  pm <- plot(fit70, type = "coefficients", term = "a")
  geoms_m <- vapply(pm$layers, function(l) class(l$geom)[1L], character(1))
  expect_equal(sum(geoms_m == "GeomSf"), 2L)
  bm <- ggplot2::ggplot_build(pm)
  expect_equal(nrow(bm$data[[2L]]), n_bad)
  expect_match(pm$labels$subtitle, sprintf("^%d of 200 locations masked", n_bad))
  # mask = FALSE leaves only the non-finite ones out and says what it hides.
  pf <- plot(fit70, type = "coefficients", term = "a", mask = FALSE)
  n_nf <- fit70$info$n_local_singular
  if (n_nf > 0L) expect_match(pf$labels$subtitle, sprintf("%d with a non-finite coefficient", n_nf))
  else expect_match(pf$labels$subtitle, "mask = FALSE")
  # Everything masked: refused, not drawn.
  fit20 <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE, bandwidth = 20)))
  expect_equal(fit20$info$n_local_collinear, 200L)
  expect_error(plot(fit20, type = "coefficients", term = "a"), "every location is masked")
})

test_that("the coefficient map is refused for a fit that has no local coefficients", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points(), predictor_vars = "w")
  expect_error(plot(fit, type = "coefficients"), "GWR fits only")
})
