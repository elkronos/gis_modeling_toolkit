# ---------------------------------------------------------------------------
# Regression tests: make_folds() must emit ..row_id values in fold splits
#
# Guards against the position-vs-ID mismatch where folds auto-generated on
# already-prepped data (with gaps in ..row_id after prep_model_data() dropped
# rows) were interpreted as row IDs by .remap_folds(), silently scrambling
# fold membership and losing observations.
# ---------------------------------------------------------------------------

test_that("make_folds(random_kfold) emits ..row_id values, not positions", {
  set.seed(1)
  n <- 20
  ids <- sort(sample(100L, n))  # non-sequential row IDs
  pts <- sf::st_sf(
    ..row_id = ids,
    geometry = sf::st_sfc(
      lapply(seq_len(n), function(i)
        sf::st_point(c(runif(1, 0, 100), runif(1, 0, 100)))),
      crs = 32632
    )
  )

  folds <- make_folds(pts, k = 4, method = "random_kfold", seed = 42)

  # Union of test sets must be exactly the row IDs (positions 1..20 would fail)
  expect_setequal(unlist(lapply(folds$folds, `[[`, "test")), ids)

  for (f in folds$folds) {
    expect_setequal(c(f$train, f$test), ids)
    expect_length(intersect(f$train, f$test), 0)
  }

  # Splits must agree with the assignment tibble
  for (j in seq_along(folds$folds)) {
    expect_setequal(folds$folds[[j]]$test,
                    folds$assignment$row_id[folds$assignment$fold == j])
  }
})


test_that("make_folds(buffered_loo) emits ..row_id values and honours buffer", {
  n <- 15
  ids <- seq(2L, by = 2L, length.out = n)  # even IDs 2, 4, ..., 30
  pts <- sf::st_sf(
    ..row_id = ids,
    geometry = sf::st_sfc(
      lapply(seq_len(n), function(i) sf::st_point(c(i * 10, 0))),
      crs = 32632
    )
  )

  folds <- make_folds(pts, k = 1, method = "buffered_loo", buffer = 15)
  expect_length(folds$folds, n)

  # Each point is a test fold exactly once, identified by its row ID
  expect_setequal(unlist(lapply(folds$folds, `[[`, "test")), ids)
  expect_equal(folds$folds[[3]]$test, ids[3])

  # Points at 10-unit spacing with buffer 15: immediate neighbours (and the
  # test point itself) must be excluded from training
  f3 <- folds$folds[[3]]
  expect_false(any(ids[c(2, 3, 4)] %in% f3$train))
  expect_true(all(ids[c(1, 5)] %in% f3$train))
})


test_that("auto-generated block folds stay aligned when prep_model_data drops rows", {
  # Replicates the cv_gwr()/cv_bayes()/cv_spatial() internal pipeline:
  # ..row_id stamped on the original data -> prep drops NA rows -> folds
  # built on the prepped subset -> .remap_folds() -> positions resolved in
  # .cv_fit_one_fold(). Before the fix, any dropped row scrambled fold
  # membership and silently excluded observations.
  set.seed(99)
  n <- 40
  df <- data.frame(
    y  = rnorm(n),
    x1 = rnorm(n),
    px = runif(n, 0, 1000),
    py = runif(n, 0, 1000)
  )
  df$x1[c(3, 17)] <- NA  # rows that prep_model_data() will drop
  data_sf <- sf::st_as_sf(df, coords = c("px", "py"), crs = 32632)
  data_sf$..row_id <- seq_len(n)

  dat_sf   <- prep_model_data(data_sf, "y", "x1")
  keep_idx <- dat_sf$..row_id
  expect_length(keep_idx, n - 2L)
  expect_false(any(c(3L, 17L) %in% keep_idx))

  folds    <- make_folds(dat_sf, k = 4, method = "block_kfold", seed = 42)
  remapped <- spatialkit:::.remap_folds(folds, keep_idx, k = 4, seed = 42)

  # No fold should be dropped, and no observation lost: the union of the
  # remapped test sets must be exactly the surviving row IDs.
  expect_length(remapped, length(folds$folds))
  expect_setequal(unlist(lapply(remapped, `[[`, "test")), keep_idx)

  for (j in seq_along(remapped)) {
    # Remapping must be the identity here: fold membership (as row IDs)
    # must survive .remap_folds() unchanged.
    expect_setequal(remapped[[j]]$test,  folds$folds[[j]]$test)
    expect_setequal(remapped[[j]]$train, folds$folds[[j]]$train)

    # And must match the assignment tibble's block/fold membership
    expect_setequal(remapped[[j]]$test,
                    folds$assignment$row_id[folds$assignment$fold == j])

    # Positional resolution in .cv_fit_one_fold() must recover the same rows
    te_pos <- match(remapped[[j]]$test, keep_idx)
    expect_false(anyNA(te_pos))
    expect_setequal(dat_sf$..row_id[te_pos], remapped[[j]]$test)
  }
})


# ---------------------------------------------------------------------------
# Fold provenance
#
# Fold splits are lists of ..row_id values, and row IDs are just
# seq_len(nrow()) unless the caller supplied them.  A `folds` object built from
# a DIFFERENT dataset of the same size therefore applied cleanly: every ID
# matched, every fold was populated, and the model was scored on splits that
# describe other observations.  Nothing in the result said so.
# ---------------------------------------------------------------------------

.probe_pts <- function(seed, n = 60) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 0.01 * sf::st_coordinates(d)[, 1] + 2 * d$a + rnorm(n)
  d
}

test_that("make_folds records a row probe in params", {
  f <- make_folds(.probe_pts(1), k = 3, method = "random_kfold", seed = 1)
  p <- f$params$row_probe
  expect_false(is.null(p))
  expect_type(p$row_id, "integer")
  expect_type(p$x, "double")
  expect_type(p$y, "double")
  expect_true(isTRUE(p$lonlat))                 # the input carried a CRS
  expect_equal(length(p$row_id), length(p$x))
  expect_true(all(p$row_id %in% seq_len(60)))
})

test_that("folds built on another dataset of the same size are refused", {
  skip_if_not_installed("ranger")
  a <- .probe_pts(1)
  b <- .probe_pts(2)                      # same n, same columns, other points
  expect_equal(nrow(a), nrow(b))

  folds_a <- make_folds(a, k = 3, method = "random_kfold", seed = 1)
  # Same data: runs.
  expect_no_error(cv_rf(a, "z", "a", folds = folds_a, num_trees = 30, seed = 1))
  # Other data: refused, and the message says why rather than reporting
  # plausible-looking metrics computed on the wrong splits.
  expect_error(
    cv_rf(b, "z", "a", folds = folds_a, num_trees = 30, seed = 1),
    "built from different data"
  )
})

test_that("the probe survives dropping rows and reprojection", {
  skip_if_not_installed("ranger")
  a <- .probe_pts(3)
  folds_a <- make_folds(a, k = 3, method = "random_kfold", seed = 1)

  # prep_model_data() drops incomplete cases, so the probe must tolerate
  # missing probe rows rather than treating them as a mismatch.
  a2 <- a; a2$a[c(2, 5, 11)] <- NA
  expect_no_error(cv_rf(a2, "z", "a", folds = folds_a, num_trees = 30, seed = 1))

  # The probe is taken in EPSG:4326, so folds built on geographic input still
  # match after make_folds()/prep_model_data() project to their own CRS.
  geo <- sf::st_transform(a, 4326)
  folds_geo <- make_folds(geo, k = 3, method = "random_kfold", seed = 1)
  expect_no_error(cv_rf(geo, "z", "a", folds = folds_geo, num_trees = 30, seed = 1))
  expect_no_error(cv_rf(a,   "z", "a", folds = folds_geo, num_trees = 30, seed = 1))
})

test_that("a folds object with no probe is passed through unchecked", {
  # Backwards compatibility: an object built by an earlier version.
  skip_if_not_installed("ranger")
  a <- .probe_pts(4)
  f <- make_folds(a, k = 3, method = "random_kfold", seed = 1)
  f$params$row_probe <- NULL
  expect_no_error(cv_rf(a, "z", "a", folds = f, num_trees = 30, seed = 1))
})


# ---------------------------------------------------------------------------
# Fourth pass: the probe must never refuse the caller's own data
# ---------------------------------------------------------------------------

test_that("the probe keeps a character ..row_id in its own type", {
  # as.integer() on a character ID gave all-NA, which matched row 1 for every
  # probe row and refused folds built on this very data.
  skip_if_not_installed("ranger")
  a <- .probe_pts(11)
  a$..row_id <- sprintf("site-%03d", seq_len(nrow(a)))
  f <- make_folds(a, k = 3, method = "random_kfold", seed = 1)
  expect_type(f$params$row_probe$row_id, "character")
  expect_no_error(cv_rf(a, "z", "a", folds = f, num_trees = 30, seed = 1))
})

test_that("folds built on polygons apply whatever pointize the cv call uses", {
  # make_folds() reduced polygons with pointize 'auto' and fingerprinted THOSE
  # points; cv_gwr()/cv_bayes()/cv_spatial() honour their own `pointize`, so
  # 'centroid' or 'bbox_center' read as different data.  Both sides now probe
  # the geometry as supplied.
  skip_if_not_installed("ranger")
  set.seed(12)
  n <- 40
  ctr <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  polys <- sf::st_sfc(lapply(seq_len(n), function(i) {
    w <- runif(1, 5, 30); h <- runif(1, 5, 30)
    sf::st_polygon(list(rbind(ctr[i, ] + c(-w, -h), ctr[i, ] + c(w, -h),
                              ctr[i, ] + c(w, h), ctr[i, ] + c(-w, 2 * h),
                              ctr[i, ] + c(-w, -h))))
  }), crs = 32632)
  d <- sf::st_sf(a = rnorm(n), geometry = polys)
  d$z <- 0.01 * ctr[, 1] + 2 * d$a + rnorm(n)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  for (pz in c("auto", "centroid", "bbox_center")) {
    expect_no_error(cv_rf(d, "z", "a", folds = f, num_trees = 30, seed = 1,
                          pointize = pz), message = pz)
  }
})

test_that("a reprojection between make_folds() and the cv call is tolerated", {
  # Comparing sprintf('%.7g') strings flipped on coordinates a round trip
  # moved by ~5e-9 degrees; the comparison is numeric with a 1e-6 degree
  # tolerance now.  A `boundary` in another projected CRS is the ordinary way
  # a cv call ends up reprojecting relative to make_folds().
  skip_if_not_installed("ranger")
  set.seed(4)
  m <- 400
  p <- sf::st_as_sf(data.frame(x = runif(m, 300000, 700000),
                               y = runif(m, 5000000, 5900000), a = rnorm(m)),
                    coords = c("x", "y"), crs = 32632)
  p$z <- 2 * p$a + rnorm(m)
  bnd <- sf::st_transform(sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(p))), 3035)
  f <- make_folds(p, k = 3, method = "block_kfold", seed = 1)
  expect_no_error(cv_rf(p, "z", "a", folds = f, num_trees = 30, seed = 1,
                        boundary = bnd))
  # And the probe itself is invariant to the route taken to EPSG:4326.
  p$..row_id <- seq_len(m)
  a <- spatialkit:::.fold_row_probe(p, max_probe = m)
  b <- spatialkit:::.fold_row_probe(sf::st_transform(p, 3035), max_probe = m)
  expect_lt(max(abs(a$x - b$x), abs(a$y - b$y)), 1e-6)
})

test_that("a duplicated ..row_id is refused up front", {
  # match() resolves a duplicated ID to its first row: the others are never
  # scored and can sit in train and test at once, with nothing said.
  skip_if_not_installed("ranger")
  a <- .probe_pts(13)
  a$..row_id <- rep(1:30, length.out = nrow(a))
  expect_error(make_folds(a, k = 3, method = "random_kfold", seed = 1), "duplicated")
  expect_error(cv_rf(a, "z", "a", k = 3, num_trees = 30, seed = 1), "duplicated")
})

test_that("n_folds_attempted counts the folds supplied, and a dropped fold warns", {
  # A fold whose test rows were all removed as incomplete vanished silently,
  # so a five-fold set reported attempted = succeeded = 4.
  skip_if_not_installed("ranger")
  a <- .probe_pts(14, n = 60)
  f <- make_folds(a, k = 5, method = "random_kfold", seed = 1)
  # Blank the predictor on every row of fold 5's test set.
  kill <- f$folds[[5]]$test
  a$a[a$..row_id %in% kill] <- NA
  if (!("..row_id" %in% names(a))) a$..row_id <- seq_len(nrow(a))
  a$a[kill] <- NA
  expect_warning(cv <- cv_rf(a, "z", "a", folds = f, num_trees = 30, seed = 1),
                 "dropped before fitting")
  expect_equal(cv$n_folds_attempted, 5L)
  expect_equal(cv$n_folds_succeeded, 4L)
})
