# Seventh audit pass.  Every block below reproduces a defect that shipped and
# is now fixed; each was found by running the code, not by reading it.

test_that(".crps_energy() survives more draws than fit in an integer square", {
  # m * m with m = nrow(draws), an integer, overflowed above 46,340 draws and
  # took every CRPS, mean_CRPS and the calibration summary to NA behind one
  # "NAs produced by integer overflow" warning.  48,000 draws is an ordinary
  # cv_bayes(fit_args = list(chains = 4, iter = 13000)) run.
  exact <- function(y) y * (2 * pnorm(y) - 1) + 2 * dnorm(y) - 1 / sqrt(pi)
  set.seed(2); y <- c(0.3, -0.7)
  for (m in c(46340L, 46341L, 50000L)) {
    dr <- matrix(rnorm(m * 2), m, 2)
    got <- spatialkit:::.crps_energy(dr, y)
    expect_true(all(is.finite(got)))
    # Against the closed-form Gaussian CRPS, to Monte-Carlo error.
    expect_equal(got, exact(y), tolerance = 0.02)
  }
})

test_that("residual_morans_i() validates k instead of silently collapsing it", {
  set.seed(3); n <- 60
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100),
                               p1 = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$value <- 2 * d$p1 + rnorm(n)
  fit <- new_spatial_fit("t7fit", engine = lm(value ~ p1, data = sf::st_drop_geometry(d)),
                         formula = value ~ p1, response_var = "value",
                         predictor_vars = "p1", data_sf = d)
  registerS3method("residuals", "t7fit", function(object, ...) residuals(object$engine))
  registerS3method("fitted", "t7fit", function(object, ...) fitted(object$engine))

  # k = c(4, 8) used to build the k = 4 matrix and return a statistic for
  # neighbours the caller never asked for, with no condition raised.
  expect_error(residual_morans_i(fit, k = c(4, 8)), "`k` must be a single number")
  expect_error(residual_morans_i(fit, k = NA), "`k` must be a single number")
  expect_error(residual_morans_i(fit, k = 0), "at least 1")

  # A weight matrix with one NA made S0 NA and the guard below it abort on
  # "missing value where TRUE/FALSE needed".
  W <- as.matrix(spatialkit:::.build_knn_weights(sf::st_coordinates(d), 8L))
  W[3, 5] <- NA
  expect_error(residual_morans_i(fit, weights = W), "non-finite value")

  # An empty POINT has a finite residual but an all-NA coordinate pair, which
  # reached FNN::get.knn() and aborted with "Data include NAs".
  d2 <- d; sf::st_geometry(d2)[7] <- sf::st_point()
  fit2 <- new_spatial_fit("t7fit", engine = lm(value ~ p1, data = sf::st_drop_geometry(d2)),
                          formula = value ~ p1, response_var = "value",
                          predictor_vars = "p1", data_sf = d2)
  # The dropped row is announced, and the warning is part of the contract.
  expect_warning(mi <- residual_morans_i(fit2), "dropping 1 row")
  expect_true(is.finite(mi$observed))
  expect_equal(mi$n, n - 1L)
})

test_that("fold labels are numbered the same way in every locale", {
  # as.factor() orders levels with LC_COLLATE, so labels differing only in
  # case landed in a different order under C than under en_US -- and the fold
  # NUMBER is the level's position, so fold_metrics$fold named different
  # groups on different machines.
  lab <- rep(c("north", "North", "south"), each = 4)
  d <- sf::st_as_sf(data.frame(x = runif(12, 0, 10), y = runif(12, 0, 10),
                               ..row_id = 1:12), coords = c("x", "y"), crs = 32632)
  f <- spatialkit:::.folds_from_labels(lab, d, "cv_spatial")
  # Radix (C) collation: "North" sorts before "north".
  expect_equal(sort(f[[1]]$test), 5:8)
  expect_equal(sort(f[[2]]$test), 1:4)
  expect_equal(sort(f[[3]]$test), 9:12)
})

test_that("a scalar argument that reaches as.integer() is validated", {
  set.seed(4); n <- 40
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100)),
                      coords = c("x", "y"), crs = 32632)
  bnd <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(pts)))

  # Each of these used to abort on "missing value where TRUE/FALSE needed" or
  # "'length = 2' in coercion to 'logical(1)'", naming nothing the caller
  # passed, because as.integer() of NA/Inf/>2^31 is NA.
  expect_error(voronoi_seeds_kmeans(pts, k = 3e9), "at most 2147483647")
  expect_error(voronoi_seeds_kmeans(pts, k = NA), "single positive number")
  expect_error(voronoi_seeds_kmeans(pts, k = c(3, 4)), "length 2")
  # k = 0 silently returned ONE seed, which "at most k" does not describe.
  expect_error(voronoi_seeds_kmeans(pts, k = 0), "at least 1")
  expect_error(voronoi_seeds_random(sf::st_as_sfc(sf::st_bbox(pts)), k = NA),
               "single positive number")
  expect_error(get_voronoi_seeds(boundary = bnd, method = "kmeans", n = 3e9,
                                 sample_points = pts), "at most 2147483647")
})

test_that("`expand` is rejected rather than silently ignored", {
  set.seed(5); n <- 40
  pts <- sf::st_as_sf(data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100)),
                      coords = c("x", "y"), crs = 32632)
  # expand = c(0.05, 0.05) returned a bbox byte-identical to expand = 0, with
  # no condition raised; a character expand did the same in the Voronoi path.
  expect_error(clip_target_for(pts, expand = c(0.05, 0.05), quiet = TRUE),
               "single non-negative number")
  expect_error(clip_target_for(pts, expand = NA_real_, quiet = TRUE),
               "single non-negative number")
  expect_error(create_voronoi_polygons(pts, expand = "5", quiet = TRUE),
               "single non-negative buffer distance")
  expect_error(create_voronoi_polygons(pts, expand = c(5, 5), quiet = TRUE),
               "single non-negative buffer distance")
  # The scalar forms still work and still differ from no expansion.
  b0 <- sf::st_bbox(clip_target_for(pts, expand = 0,    quiet = TRUE))
  b1 <- sf::st_bbox(clip_target_for(pts, expand = 0.05, quiet = TRUE))
  expect_gt(as.numeric(b0["xmin"]), as.numeric(b1["xmin"]))
})

test_that("summarize_by_cell(deff = 'variogram') works with no value column", {
  skip_if_not_installed("gstat")
  set.seed(6); n <- 250
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D <- as.matrix(dist(xy))
  xy$value <- as.numeric(t(chol(exp(-D / 60) + diag(1e-8, n))) %*% rnorm(n))
  p <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  hull <- clip_target_for(p, quiet = TRUE)
  tess <- build_tessellation(p, boundary = hull, method = "square",
                             approx_n_cells = 9, quiet = TRUE)
  asg <- assign_features_to_polygons(p, tess$cells)
  sac <- estimate_sac_range(p, "value")
  skip_if(is.na(sac), "no identified range on this draw")
  # `primary_col` is NULL when neither response_var nor predictor_vars is
  # given, and df[[NULL]] aborted with "attempt to select less than one
  # element in get1index".  The variogram design effect is a function of the
  # coordinates, so there is nothing to stop it being computed.
  out <- summarize_by_cell(asg, cells_sf = tess$cells, deff = "variogram", sac = sac)
  expect_s3_class(out, "sf")
  expect_true(all(c("n", "cell_weight") %in% names(out)))
  expect_true(any(is.finite(out$cell_weight)))
})

test_that("compare_models() says so when nothing is a spatial_fit", {
  # evaluate_insample() returns NULL when every element was skipped, and
  # `met_df$AICc <- NA_real_` turned that into a bare list, after which
  # seq_len(nrow(NULL)) aborted with "argument must be coercible to
  # non-negative integer".
  expect_error(compare_models(list(a = 1, b = 2)), "no element of `models`")
})

test_that("a fitted-value cache entry belongs to the engine that produced it", {
  set.seed(7); n <- 30
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + rnorm(n)
  f1 <- new_spatial_fit("bayesian_fit", engine = list(tag = "E1"), formula = z ~ a,
                        response_var = "z", predictor_vars = "a", data_sf = d)
  key <- spatialkit:::.fitted_cache_key(f1)
  warm <- function() assign(".fitted_values",
                            list(n = f1$n, key = key, engine = f1$engine,
                                 values = rep(10, f1$n)),
                            envir = f1$info$.cache)
  warm()
  # The cache still works for the fit that wrote it.
  expect_equal(fitted(f1)[[1]], 10)

  # The key digests n, the attribute table and the coordinates, so two fits
  # over the SAME data hash identically however different their engines are --
  # and the cache is shared by every copy of a fit.  `cp$engine <- <other>`
  # therefore used to read the first engine's fitted values back out.
  cp <- f1; cp$engine <- list(tag = "E2")
  # A cache HIT returns the stored 10 silently, so any error here proves the
  # miss.  The pattern matches the method's own message in BOTH environments:
  # with brms installed it fails inside posterior_epred(), without it at the
  # requireNamespace() guard.  Pinning the posterior_epred wording made this
  # test pass only where brms happens to be installed, which is not the CI
  # matrix.
  expect_error(fitted(cp), "fitted\\.bayesian_fit\\(\\)")
  expect_equal(fitted(f1)[[1]], 10)             # and did not evict the valid entry

  # summary() used to hand out the live cache environment, so a summary was
  # not a value snapshot and clear_fitted_cache(summary) emptied the FIT.
  s <- summary(f1)
  expect_false(".cache" %in% names(s$info))
  expect_false(any(vapply(s$info, is.environment, logical(1))))
  clear_fitted_cache(s)
  expect_true(exists(".fitted_values", envir = f1$info$.cache, inherits = FALSE))

  # An `info` list carrying another fit's cache cannot import it.
  f3 <- new_spatial_fit("bayesian_fit", engine = list(tag = "E3"), formula = z ~ a,
                        response_var = "z", predictor_vars = "a", data_sf = d,
                        info = f1$info)
  expect_false(identical(f3$info$.cache, f1$info$.cache))
})

test_that("the grid cache's order registry does not outlive its environment", {
  bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0)))), crs = 32632))
  meta <- spatialkit:::.gmt_cache_meta
  before <- length(ls(meta, all.names = TRUE))
  for (i in 1:20)
    invisible(create_grid_polygons_cached(bnd, target_cells = 9, type = "hex",
                                          cache_env = new.env()))
  gc()
  # Each entry was a permanent character vector keyed by the environment's
  # printed address; the environments are collected, so no caller could ever
  # name them to clear them.  A finalizer now removes the entry with the env.
  expect_lte(length(ls(meta, all.names = TRUE)), before + 2L)
})
