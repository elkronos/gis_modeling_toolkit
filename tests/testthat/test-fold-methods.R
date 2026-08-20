# tests/testthat/test-fold-methods.R
# ---------------------------------------------------------------------------
# leave_location_out and nndm.
#
# make_folds() previously offered random_kfold, block_kfold and buffered_loo.
# Repeated measurements at the same site were unrepresentable -- random k-fold
# splits them across folds and scores the model on sites it has already seen --
# and buffered_loo required picking a buffer with nothing to justify it.
# ---------------------------------------------------------------------------

site_points <- function(n_sites = 12, per_site = 5, extent = 1000, seed = 1) {
  set.seed(seed)
  sx <- runif(n_sites, 0, extent); sy <- runif(n_sites, 0, extent)
  i  <- rep(seq_len(n_sites), each = per_site)
  # ..row_id is set explicitly: make_folds() would assign one internally, but
  # on its own copy, so the caller's object would not carry the ids the fold
  # assignment refers to.
  out <- sf::st_as_sf(
    data.frame(x = sx[i] + rnorm(length(i), 0, 5),
               y = sy[i] + rnorm(length(i), 0, 5),
               site = paste0("s", i),
               z = rnorm(length(i))),
    coords = c("x", "y"), crs = 3857
  )
  out$..row_id <- seq_len(nrow(out))
  out
}


# ---- leave-location-out ----------------------------------------------------

test_that("no location is split across folds", {
  # The whole point: every observation from a site must share a fold, or the
  # model is scored on sites it trained on.
  pts <- site_points()
  f <- make_folds(pts, k = 4, method = "leave_location_out",
                  group_var = "site", seed = 1)

  asg <- f$assignment
  site_of <- sf::st_drop_geometry(pts)$site[match(asg$row_id, pts$..row_id)]
  per_site_folds <- tapply(asg$fold, site_of, function(z) length(unique(z)))
  expect_true(all(per_site_folds == 1L))
})

test_that("train and test never share a location", {
  pts <- site_points()
  f <- make_folds(pts, k = 4, method = "leave_location_out",
                  group_var = "site", seed = 1)
  site_by_id <- stats::setNames(sf::st_drop_geometry(pts)$site, pts$..row_id)
  expect_false(anyNA(names(site_by_id)))   # NA names would make the test vacuous

  for (fold in f$folds) {
    tr <- unique(site_by_id[as.character(fold$train)])
    te <- unique(site_by_id[as.character(fold$test)])
    expect_length(intersect(tr, te), 0L)
  }
})

test_that("every observation appears in exactly one test fold", {
  pts <- site_points()
  f <- make_folds(pts, k = 3, method = "leave_location_out",
                  group_var = "site", seed = 1)
  tested <- sort(unlist(lapply(f$folds, `[[`, "test")))
  expect_false(is.null(pts$..row_id))
  expect_equal(tested, sort(pts$..row_id))
})

test_that("k above the group count degenerates to leave-one-group-out", {
  pts <- site_points(n_sites = 5, per_site = 4)
  f <- make_folds(pts, k = 99, method = "leave_location_out",
                  group_var = "site", seed = 1)
  expect_equal(f$k, 5L)
  expect_length(f$folds, 5L)
})

test_that("leave_location_out validates its group column", {
  pts <- site_points()
  expect_error(make_folds(pts, k = 3, method = "leave_location_out"),
               "group_var")
  expect_error(make_folds(pts, k = 3, method = "leave_location_out",
                          group_var = "nope"), "group_var")

  pts$site[1] <- NA
  expect_error(make_folds(pts, k = 3, method = "leave_location_out",
                          group_var = "site"), "contains NA")
})


# ---- NNDM ------------------------------------------------------------------

nndm_setup <- function(n = 80, seed = 2) {
  set.seed(seed)
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
    coords = c("x", "y"), crs = 3857)
  pts$..row_id <- seq_len(nrow(pts))
  grid <- sf::st_as_sf(expand.grid(x = seq(50, 950, by = 100),
                                   y = seq(50, 950, by = 100)),
                       coords = c("x", "y"), crs = 3857)
  list(pts = pts, grid = grid)
}

test_that("nndm produces one leave-one-out fold per observation", {
  s <- nndm_setup()
  f <- make_folds(s$pts, method = "nndm", prediction_points = s$grid, seed = 1)
  expect_equal(f$method, "nndm")
  expect_length(f$folds, nrow(s$pts))
  expect_true(all(vapply(f$folds, function(z) length(z$test), integer(1)) == 1L))
})

test_that("nndm excludes near neighbours when predictions sit far from training", {
  # Exclusion is only warranted when the prediction locations are further from
  # the training data than training points are from each other.  Clustered
  # training against a spread grid is that case.
  set.seed(11)
  ctr <- cbind(c(200, 800, 500), c(200, 800, 800))
  i   <- rep(1:3, each = 25)
  pts <- sf::st_as_sf(
    data.frame(x = ctr[i, 1] + rnorm(75, 0, 20),
               y = ctr[i, 2] + rnorm(75, 0, 20),
               z = rnorm(75)),
    coords = c("x", "y"), crs = 3857)
  pts$..row_id <- seq_len(nrow(pts))
  grid <- sf::st_as_sf(expand.grid(x = seq(50, 950, by = 90),
                                   y = seq(50, 950, by = 90)),
                       coords = c("x", "y"), crs = 3857)

  f <- make_folds(pts, method = "nndm", prediction_points = grid, seed = 1)

  expect_gt(f$params$target_median, 0)
  expect_gt(f$params$median_excluded, 0)
  train_sizes <- vapply(f$folds, function(z) length(z$train), integer(1))
  expect_true(any(train_sizes < nrow(pts) - 1L))
})

test_that("nndm leaves LOO alone when it already matches the target", {
  # The converse, and the reason the previous expectation was wrong: when
  # prediction and training densities are comparable, plain LOO already
  # reproduces the target distribution and nothing needs removing.  Excluding
  # anyway would over-buffer and throw away usable training data.
  set.seed(12); n <- 80
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), z = rnorm(n)),
    coords = c("x", "y"), crs = 3857)
  pts$..row_id <- seq_len(n)
  grid <- sf::st_as_sf(expand.grid(x = seq(50, 950, by = 100),
                                   y = seq(50, 950, by = 100)),
                       coords = c("x", "y"), crs = 3857)

  f <- make_folds(pts, method = "nndm", prediction_points = grid, seed = 1)
  expect_lte(f$params$median_excluded, 1)
  expect_true(all(vapply(f$folds, function(z) length(z$train), integer(1)) >= 2L))
})

test_that("the realised distance distribution tracks the target", {
  # The property NNDM exists for.  Without recording both, there is no way to
  # tell whether the matching worked or the buffers were decorative.
  set.seed(13)
  ctr <- cbind(c(250, 750), c(250, 750))
  i   <- rep(1:2, each = 40)
  pts <- sf::st_as_sf(
    data.frame(x = ctr[i, 1] + rnorm(80, 0, 40),
               y = ctr[i, 2] + rnorm(80, 0, 40),
               z = rnorm(80)),
    coords = c("x", "y"), crs = 3857)
  pts$..row_id <- seq_len(80)
  grid <- sf::st_as_sf(expand.grid(x = seq(50, 950, by = 75),
                                   y = seq(50, 950, by = 75)),
                       coords = c("x", "y"), crs = 3857)

  f <- make_folds(pts, method = "nndm", prediction_points = grid, seed = 1)
  tgt <- f$params$target_median
  rea <- f$params$realised_median

  expect_true(is.finite(tgt) && is.finite(rea))

  # Realised must be closer to the target than plain LOO is.  This is the
  # regression guard for a bug where excluding everything inside a hard radius
  # removed an entire cluster, leaving the nearest survivor in the NEXT cluster
  # -- realised 602 against a target of 200, which is worse than the 12 plain
  # LOO would have given.  An implementation that overshoots is not merely
  # imprecise, it is worse than doing nothing.
  loo_median <- stats::median(vapply(seq_len(80), function(k) {
    d <- as.numeric(sf::st_distance(sf::st_geometry(pts)[k], sf::st_geometry(pts)))
    min(d[-k])
  }, numeric(1)))
  expect_lt(abs(rea - tgt), abs(loo_median - tgt))

  # And it must not overshoot the target the way the hard-radius version did.
  expect_lt(rea, 2 * tgt)
})

test_that("nndm never lands further from the target than plain LOO", {
  # The property in general, across several point configurations: matching
  # should never make the distance profile worse.
  for (sd_ in c(20, 60, 150)) {
    set.seed(20 + sd_)
    ctr <- cbind(c(250, 750), c(250, 750))
    i   <- rep(1:2, each = 40)
    pts <- sf::st_as_sf(
      data.frame(x = ctr[i, 1] + rnorm(80, 0, sd_),
                 y = ctr[i, 2] + rnorm(80, 0, sd_),
                 z = rnorm(80)),
      coords = c("x", "y"), crs = 3857)
    pts$..row_id <- seq_len(80)
    grid <- sf::st_as_sf(expand.grid(x = seq(50, 950, by = 90),
                                     y = seq(50, 950, by = 90)),
                         coords = c("x", "y"), crs = 3857)

    f <- make_folds(pts, method = "nndm", prediction_points = grid, seed = 1)
    loo <- stats::median(vapply(seq_len(80), function(k) {
      d <- as.numeric(sf::st_distance(sf::st_geometry(pts)[k], sf::st_geometry(pts)))
      min(d[-k])
    }, numeric(1)))
    expect_lte(abs(f$params$realised_median - f$params$target_median),
               abs(loo - f$params$target_median))
  }
})

test_that("a held-out point is never in its own training set", {
  s <- nndm_setup()
  f <- make_folds(s$pts, method = "nndm", prediction_points = s$grid, seed = 1)
  for (fold in f$folds) expect_false(fold$test %in% fold$train)
})

test_that("no training set is left empty", {
  # Dense prediction points imply large buffers; the fold must still be usable.
  s <- nndm_setup(n = 40)
  dense <- sf::st_as_sf(expand.grid(x = seq(0, 1000, by = 25),
                                    y = seq(0, 1000, by = 25)),
                        coords = c("x", "y"), crs = 3857)
  f <- make_folds(s$pts, method = "nndm", prediction_points = dense, seed = 1)
  expect_true(all(vapply(f$folds, function(z) length(z$train), integer(1)) >= 2L))
})

test_that("nndm requires prediction points", {
  s <- nndm_setup()
  expect_error(make_folds(s$pts, method = "nndm"), "prediction_points")
})

test_that("nndm buffers track the prediction geometry", {
  # Sparser prediction points sit further from training data, so the matched
  # buffers should be larger.
  s <- nndm_setup()
  sparse <- sf::st_as_sf(expand.grid(x = c(100, 900), y = c(100, 900)),
                         coords = c("x", "y"), crs = 3857)
  dense  <- sf::st_as_sf(expand.grid(x = seq(0, 1000, by = 50),
                                     y = seq(0, 1000, by = 50)),
                         coords = c("x", "y"), crs = 3857)

  f_sparse <- make_folds(s$pts, method = "nndm", prediction_points = sparse, seed = 1)
  f_dense  <- make_folds(s$pts, method = "nndm", prediction_points = dense,  seed = 1)
  expect_gt(f_sparse$params$median_buffer, f_dense$params$median_buffer)
})

test_that("nndm is reproducible from the seed", {
  s <- nndm_setup()
  a <- make_folds(s$pts, method = "nndm", prediction_points = s$grid, seed = 7)
  b <- make_folds(s$pts, method = "nndm", prediction_points = s$grid, seed = 7)
  expect_equal(a$params$median_buffer, b$params$median_buffer)
  expect_equal(lapply(a$folds, `[[`, "train"), lapply(b$folds, `[[`, "train"))
})

# --- regressions found in the full-package audit -----------------------------

# Plain scattered points; site_points() clusters by design, which is not what
# the RNG / geometry / block-size checks below are about.
plain_points <- function(n = 60, extent = 1000, seed = 5) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent),
               z = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
}

test_that("make_folds() restores the caller's RNG stream on every method", {
  # The nndm branch rebound `cleanup`, the same name the outer .with_seed()
  # handler used.  on.exit() stores the EXPRESSION, so both registered handlers
  # resolved to the inner closure and the outer one -- the only holder of the
  # caller's pre-call .Random.seed -- never ran, silently replacing the
  # caller's stream with the post-set.seed(seed) state.
  pts <- plain_points(60)

  # buffered_loo requires an explicit buffer, so pass the per-method arguments
  # each one actually needs.
  cases <- list(
    list(method = "random_kfold", args = list(k = 4)),
    list(method = "block_kfold",  args = list(k = 4)),
    list(method = "buffered_loo", args = list(k = 4, buffer = 50))
  )
  for (cs in cases) {
    set.seed(4242); expected <- runif(3)
    set.seed(4242)
    invisible(suppressMessages(do.call(
      make_folds,
      c(list(points_sf = pts, method = cs$method, seed = 99), cs$args))))
    expect_equal(runif(3), expected, info = cs$method)
  }

  grid <- plain_points(10, seed = 7)
  set.seed(4242); expected <- runif(3)
  set.seed(4242)
  invisible(suppressMessages(
    make_folds(pts, method = "nndm", prediction_points = grid, seed = 99)))
  expect_equal(runif(3), expected)
})

test_that("nndm survives a non-integral median exclusion count", {
  # median() of an even-length integer vector returns a double; sprintf("%d",
  # 0.5) is a hard error, so the whole call died before returning.
  pts <- sf::st_as_sf(
    data.frame(x = c(0, 3, 14, 1000, 1016, 1032), y = 0),
    coords = c("x", "y"), crs = 3857
  )
  grid <- sf::st_as_sf(data.frame(x = -15, y = 0),
                       coords = c("x", "y"), crs = 3857)
  f <- suppressMessages(make_folds(pts, method = "nndm",
                                   prediction_points = grid, seed = 1))
  expect_identical(length(f$folds), 6L)
  expect_true(is.finite(f$params$median_excluded))
})

test_that("MULTIPOINT input is coerced before folds are built", {
  # st_coordinates() returns one row per VERTEX, so a multi-vertex MULTIPOINT
  # made xy[i, ] a different feature than row i and every fold misaligned.
  g <- sf::st_sfc(
    sf::st_multipoint(rbind(c(0, 0), c(10, 10))),
    sf::st_point(c(100, 100)),
    sf::st_point(c(200, 200)),
    sf::st_point(c(300, 300)),
    crs = 32632
  )
  pts <- sf::st_sf(z = c(1, 2, 3, 4), geometry = g)
  f <- suppressMessages(make_folds(pts, k = 2, method = "random_kfold", seed = 1))
  expect_identical(sort(unlist(lapply(f$folds, `[[`, "test"))), 1:4)
})

test_that("a single block is refused rather than silently defeating blocked CV", {
  # One block means one fold with an EMPTY training set: blocked CV degenerating
  # into nothing, reported as a successful run that happens to score NA.
  pts <- plain_points(60)
  bb  <- sf::st_bbox(pts)
  huge <- 2 * max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  expect_error(
    suppressMessages(make_folds(pts, k = 4, method = "block_kfold",
                                block_size = huge, seed = 1)),
    "single block covering the whole extent"
  )
})
