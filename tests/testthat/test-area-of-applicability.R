# =============================================================================
# Area of applicability
#
# The method is arithmetic, so almost everything here is checked against a
# hand-computed value rather than against the implementation's own output.
# The one test that cannot be hand-computed -- FNN's exact kd-tree distances
# versus the dense |a|^2+|b|^2-2a.b path -- is the reason the FNN backend is
# injectable.
# =============================================================================

mk_aoa_pts <- function(n = 100, seed = 3, a = NULL, b = NULL) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = if (is.null(a)) rnorm(n) else a,
               b = if (is.null(b)) rnorm(n) else b),
    coords = c("x", "y"), crs = 32632
  )
}

mk_aoa_new <- function(a, b) {
  sf::st_as_sf(
    data.frame(x = seq(10, by = 10, length.out = length(a)),
               y = seq(10, by = 10, length.out = length(a)),
               a = a, b = b),
    coords = c("x", "y"), crs = 32632
  )
}


# --- .aoa_matrix -------------------------------------------------------------

test_that(".aoa_matrix extracts a numeric matrix and drops geometry", {
  pts <- mk_aoa_pts(5)
  m <- .aoa_matrix(pts, c("a", "b"))
  expect_true(is.matrix(m))
  expect_identical(dim(m), c(5L, 2L))
  expect_identical(colnames(m), c("a", "b"))
  expect_identical(typeof(m), "double")
  expect_identical(unname(m[, "a"]), sf::st_drop_geometry(pts)$a)
})

test_that(".aoa_matrix names the offending columns", {
  pts <- mk_aoa_pts(5)
  expect_error(.aoa_matrix(pts, c("a", "nope"), "the training data"),
               "absent from the training data.*nope")
  pts$grp <- factor(letters[1:5])
  expect_error(.aoa_matrix(pts, c("a", "grp"), "`newdata`"),
               "not numeric")
  expect_error(.aoa_matrix(pts, c("a", "grp"), "`newdata`"),
               "categorical")
})


# --- scaling -----------------------------------------------------------------

test_that(".aoa_apply_scaling standardises the training data", {
  X <- cbind(a = c(1, 2, 3, 4), b = c(10, 20, 30, 40))
  sc <- .aoa_scaling(X)
  Z <- .aoa_apply_scaling(X, sc)
  expect_equal(unname(colMeans(Z)), c(0, 0))
  expect_equal(unname(apply(Z, 2, stats::sd)), c(1, 1))
})

test_that(".aoa_apply_scaling uses TRAINING statistics for new data", {
  # The classic bug is scaling newdata by its own mean/sd, which would centre
  # a far-away block at zero and make it look perfectly familiar.
  Xtr <- cbind(a = c(1, 2, 3, 4))
  Xnw <- cbind(a = c(101, 102, 103, 104))
  sc  <- .aoa_scaling(Xtr)
  Znw <- .aoa_apply_scaling(Xnw, sc)
  expect_false(isTRUE(all.equal(mean(Znw[, 1]), 0)))
  # (101 - 2.5) / sd(1:4).  Read with [[ ]]: Z[1, 1] is legitimately named
  # here, because when a fully-dropped array result has only one dimension
  # carrying dimnames, R attaches them.
  expect_equal(Znw[[1, 1]], (101 - 2.5) / stats::sd(1:4))
})

test_that(".aoa_apply_scaling returns a clean matrix", {
  # sweep() keeps names on its STATS vector, so an unguarded implementation can
  # leave a `names` attribute grafted on alongside the dimnames.
  X  <- cbind(a = c(1, 2, 3, 4), b = c(10, 20, 30, 40))
  Z  <- .aoa_apply_scaling(X, .aoa_scaling(X))
  expect_null(names(Z))
  expect_identical(dimnames(Z), list(NULL, c("a", "b")))
  expect_equal(Z[[1, 1]], (1 - 2.5) / stats::sd(1:4))
})

test_that("the returned index vectors carry no names", {
  tr  <- mk_aoa_pts(60)
  nd  <- mk_aoa_new(a = c(0, 15), b = c(0, 15))
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  expect_null(names(res$aoa$DI))
  expect_null(names(res$train_DI))
  expect_null(names(res$aoa$AOA))
  expect_identical(names(res$weights), c("a", "b"))   # these are meant to be named
})

test_that(".aoa_scaling flags zero-variance predictors and .aoa_apply_scaling drops them", {
  X <- cbind(a = c(1, 2, 3), const = c(7, 7, 7))
  sc <- .aoa_scaling(X)
  expect_identical(unname(sc$keep), c(TRUE, FALSE))
  Z <- .aoa_apply_scaling(X, sc)
  expect_identical(colnames(Z), "a")
})

test_that(".aoa_apply_scaling errors when nothing has variance", {
  X  <- cbind(a = c(1, 1, 1), b = c(2, 2, 2))
  sc <- .aoa_scaling(X)
  expect_error(.aoa_apply_scaling(X, sc), "no predictor has usable variance")
})


# --- .aoa_weight_vector ------------------------------------------------------

test_that(".aoa_weight_vector defaults to equal weights", {
  w <- .aoa_weight_vector(NULL, c("a", "b"))
  expect_identical(w, c(a = 1, b = 1))
})

test_that(".aoa_weight_vector reorders named weights to match the predictors", {
  w <- .aoa_weight_vector(c(b = 3, a = 1), c("a", "b"))
  expect_identical(names(w), c("a", "b"))
  expect_equal(unname(w), c(0.5, 1.5))          # rescaled to mean 1
})

test_that(".aoa_weight_vector accepts an unnamed vector positionally", {
  w <- .aoa_weight_vector(c(1, 3), c("a", "b"))
  expect_identical(names(w), c("a", "b"))
  expect_equal(unname(w), c(0.5, 1.5))
})

test_that(".aoa_weight_vector rejects unusable weights", {
  v <- c("a", "b")
  expect_error(.aoa_weight_vector(c(1, 2, 3), v), "one value per predictor")
  expect_error(.aoa_weight_vector(c(a = 1), v), "no entry for")
  expect_error(.aoa_weight_vector(c(a = 1, b = -1), v), "non-negative")
  expect_error(.aoa_weight_vector(c(a = 0, b = 0), v), "all .weights. are zero")
  expect_error(.aoa_weight_vector(c(a = "x", b = "y"), v), "must be numeric")
  expect_error(.aoa_weight_vector(c(a = 1, b = NA_real_), v), "finite")
})


# --- .aoa_min_dist -----------------------------------------------------------

AOA_D <- rbind(c(0, 0), c(3, 4), c(10, 0))
AOA_Q <- rbind(c(0, 1), c(3, 0))
# q1: 1, sqrt(18), sqrt(101) -> 1     q2: 3, 4, 7 -> 3

test_that(".aoa_min_dist matches hand-computed nearest distances (dense)", {
  expect_equal(.aoa_min_dist(AOA_Q, AOA_D, use_fnn = FALSE), c(1, 3))
})

test_that(".aoa_min_dist agrees between the FNN and dense paths", {
  skip_if_not_installed("FNN")
  set.seed(4)
  D <- matrix(rnorm(300 * 4), ncol = 4)
  Q <- matrix(rnorm(120 * 4), ncol = 4)
  expect_equal(.aoa_min_dist(Q, D, use_fnn = TRUE),
               .aoa_min_dist(Q, D, use_fnn = FALSE), tolerance = 1e-8)
})

test_that(".aoa_min_dist is unaffected by chunking", {
  set.seed(5)
  D <- matrix(rnorm(200 * 3), ncol = 3)
  Q <- matrix(rnorm(77 * 3), ncol = 3)
  full <- .aoa_min_dist(Q, D, use_fnn = FALSE)
  expect_equal(.aoa_min_dist(Q, D, use_fnn = FALSE, chunk_size = 1L), full)
  expect_equal(.aoa_min_dist(Q, D, use_fnn = FALSE, chunk_size = 13L), full)
})

test_that(".aoa_min_dist honours `exclude`", {
  # Excluding the identical row must give the distance to the NEXT nearest.
  d <- .aoa_min_dist(AOA_D, AOA_D, use_fnn = FALSE, exclude = seq_len(3))
  expect_equal(d, c(5, 5, sqrt(49 + 16)))
  # Without exclusion every point matches itself.
  expect_equal(.aoa_min_dist(AOA_D, AOA_D, use_fnn = FALSE), c(0, 0, 0))
})

test_that(".aoa_min_dist handles degenerate shapes", {
  expect_length(.aoa_min_dist(AOA_Q[0, , drop = FALSE], AOA_D,
                              use_fnn = FALSE), 0L)
  expect_error(.aoa_min_dist(AOA_Q, AOA_D[0, , drop = FALSE], use_fnn = FALSE),
               "no reference rows")
})


# --- .aoa_normalizer ---------------------------------------------------------

test_that(".aoa_normalizer is the mean pairwise distance", {
  X <- rbind(c(0, 0), c(3, 4), c(0, 4))        # pairwise 5, 4, 3 -> mean 4
  n <- .aoa_normalizer(X)
  expect_equal(n$value, 4)
  expect_false(n$subsampled)
  expect_identical(n$n_used, 3L)
})

test_that(".aoa_normalizer subsamples and records it", {
  set.seed(6)
  X <- matrix(rnorm(300 * 2), ncol = 2)
  n <- .aoa_normalizer(X, max_n = 50L, seed = 1L)
  expect_true(n$subsampled)
  expect_identical(n$n_used, 50L)
  # deterministic given the seed
  expect_equal(n$value, .aoa_normalizer(X, max_n = 50L, seed = 1L)$value)
})

test_that(".aoa_normalizer refuses data with no spread", {
  X <- matrix(rep(1, 10), ncol = 2)
  expect_error(.aoa_normalizer(X), "no spread")
  expect_error(.aoa_normalizer(matrix(1:2, ncol = 2)),
               "at least two training rows")
})


# --- .aoa_fold_splits --------------------------------------------------------

test_that(".aoa_fold_splits builds splits from a label vector", {
  sp <- .aoa_fold_splits(c(1, 1, 2, 2, 2), 5L)
  expect_length(sp, 2L)
  expect_identical(sp[[1]]$test, c(1L, 2L))
  expect_identical(sp[[1]]$train, c(3L, 4L, 5L))
  expect_identical(sp[[2]]$test, c(3L, 4L, 5L))
})

test_that(".aoa_fold_splits reads a make_folds()-shaped object", {
  folds <- list(method = "block_kfold", k = 2L,
                folds = list(list(train = 3:5, test = 1:2),
                             list(train = 1:2, test = 3:5)))
  sp <- .aoa_fold_splits(folds, 5L)
  expect_length(sp, 2L)
  expect_identical(sp[[1]]$test, 1:2)
})

test_that(".aoa_fold_splits keeps the supplied train set rather than the complement", {
  # A buffered fold deliberately withholds rows 2 and 3 from fold 1's training
  # set.  Reconstructing "everything not in test" would silently undo that.
  sp <- .aoa_fold_splits(list(list(test = 1L,   train = c(4L, 5L)),
                              list(test = 2:5,  train = 1L)), 5L)
  expect_identical(sp[[1]]$train, c(4L, 5L))
  expect_identical(sp[[2]]$train, 1L)
})

test_that(".aoa_fold_splits rejects incoherent fold specifications", {
  expect_error(.aoa_fold_splits(c(1, 2), 5L), "2 labels but the training data has 5")
  expect_error(.aoa_fold_splits(rep(1, 5), 5L), "at least two")
  expect_error(.aoa_fold_splits(list(list(test = 1L, train = 9L)), 5L),
               "outside 1:5")
  expect_error(.aoa_fold_splits(list(list(test = 1L, train = integer(0)),
                                     list(test = 2:5, train = 1L)), 5L),
               "empty training set")
  expect_error(.aoa_fold_splits(list(list(test = 1:3, train = 4:5),
                                     list(test = 3:5, train = 1:2)), 5L),
               "appears in the test set of both fold")
  expect_error(.aoa_fold_splits(list(list(test = 1:2, train = 3:5)), 5L),
               "in no test fold")
  expect_error(.aoa_fold_splits(list(list(a = 1)), 5L), "not a recognised|must be a make_folds")
})

test_that(".aoa_fold_splits passes NULL through", {
  expect_null(.aoa_fold_splits(NULL, 5L))
})


# --- .aoa_train_dist ---------------------------------------------------------

AOA_X3 <- rbind(c(0, 0), c(3, 4), c(0, 4))     # pairwise 5, 4, 3

test_that(".aoa_train_dist computes leave-one-out nearest neighbours", {
  expect_equal(.aoa_train_dist(AOA_X3, use_fnn = FALSE), c(4, 3, 3))
})

test_that(".aoa_train_dist agrees between the FNN and dense paths", {
  skip_if_not_installed("FNN")
  set.seed(7)
  X <- matrix(rnorm(250 * 3), ncol = 3)
  expect_equal(.aoa_train_dist(X, use_fnn = TRUE),
               .aoa_train_dist(X, use_fnn = FALSE), tolerance = 1e-8)
})

test_that(".aoa_train_dist uses each fold's own training set", {
  X <- matrix(c(0, 1, 2, 10, 11), ncol = 1)
  sp <- list(list(test = 1L,  train = c(4L, 5L)),   # buffered: 2 and 3 withheld
             list(test = 2:5, train = 1L))
  # Row 1's nearest available training point is 10 away, not 1 away.
  expect_equal(.aoa_train_dist(X, splits = sp, use_fnn = FALSE),
               c(10, 1, 2, 10, 11))
})


# --- .aoa_threshold ----------------------------------------------------------

test_that(".aoa_threshold is the outlier-removed maximum", {
  # Q3 = 1, IQR = 0, fence = 1; the 10 is an outlier and is discarded.
  expect_equal(.aoa_threshold(c(rep(1, 9), 10)), 1)
  # Q3 = 4, IQR = 2, fence = 7; nothing is an outlier, so the max stands.
  expect_equal(.aoa_threshold(c(1, 2, 3, 4, 5)), 5)
})

test_that(".aoa_threshold ignores non-finite values", {
  expect_equal(.aoa_threshold(c(1, 2, 3, 4, 5, NA, Inf)), 5)
  expect_error(.aoa_threshold(c(NA_real_, NaN)), "no finite training DI")
})


# --- area_of_applicability ---------------------------------------------------

test_that("area_of_applicability separates familiar from unfamiliar points", {
  tr <- mk_aoa_pts(100)
  nd <- mk_aoa_new(a = c(0, 20), b = c(0, 20))
  res <- area_of_applicability(nd, train_sf = tr,
                               predictor_vars = c("a", "b"))
  expect_s3_class(res, "aoa")
  expect_true(inherits(res$aoa, "sf"))
  expect_true(all(c("DI", "AOA") %in% names(res$aoa)))
  expect_type(res$aoa$AOA, "logical")
  expect_true(res$aoa$AOA[1])          # centre of the training distribution
  expect_false(res$aoa$AOA[2])         # 20 sd out
  expect_lt(res$aoa$DI[1], res$aoa$DI[2])
  expect_identical(res$n_train, 100L)
  expect_identical(res$n_new, 2L)
  expect_identical(res$n_inside + res$n_outside + res$n_na, 2L)
})

test_that("weights actually change the dissimilarity index", {
  # The previous version of this test compared weights c(1,3) against
  # c(100,300).  .aoa_weight_vector() rescales to mean 1, so BOTH normalise to
  # c(0.5, 1.5) and the two results were identical by construction -- it could
  # not have failed however the weights were used, or whether they were used at
  # all.  Weighting had no real coverage anywhere in this file.
  tr <- mk_aoa_pts(80)
  nd <- mk_aoa_new(a = c(3, 0), b = c(0, 3))

  even  <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  heavy_a <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                   weights = c(a = 10, b = 0.1))
  expect_false(isTRUE(all.equal(even$aoa$DI, heavy_a$aoa$DI)))

  # Point 1 is displaced in `a` only, point 2 in `b` only.  Up-weighting `a`
  # must push point 1 further out relative to point 2 than equal weights do.
  expect_gt(heavy_a$aoa$DI[1] / heavy_a$aoa$DI[2],
            even$aoa$DI[1] / even$aoa$DI[2])
})

test_that("weights are applied to the training data and to newdata alike", {
  # Sweeping the weights over Z_tr but not Z_nw (or vice versa) is the classic
  # asymmetry bug.  It survives any test that only compares two weighted runs
  # to each other, so compare against a hand-built equivalent instead: scaling
  # a predictor by k in BOTH data sets is undone by the training-derived
  # standardisation, so weighting `a` by k must equal that.
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(2, 0, 1), b = c(0, 2, 1))

  w <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                             weights = c(a = 4, b = 1))
  # Equivalent by hand: standardise, multiply column a by its normalised
  # weight, and compute the index from first principles.
  Xtr <- as.matrix(sf::st_drop_geometry(tr)[, c("a", "b")])
  Xnw <- as.matrix(sf::st_drop_geometry(nd)[, c("a", "b")])
  ctr <- colMeans(Xtr); scl <- apply(Xtr, 2, stats::sd)
  Ztr <- sweep(sweep(Xtr, 2, ctr, "-"), 2, scl, "/")
  Znw <- sweep(sweep(Xnw, 2, ctr, "-"), 2, scl, "/")
  wv  <- c(a = 4, b = 1); wv <- wv * (2 / sum(wv))       # rescale to mean 1
  Ztr <- sweep(Ztr, 2, wv, "*"); Znw <- sweep(Znw, 2, wv, "*")
  norm <- mean(as.numeric(stats::dist(Ztr)))
  di <- apply(Znw, 1, function(q)
    min(sqrt(colSums((t(Ztr) - q)^2)))) / norm

  expect_equal(unname(w$aoa$DI), unname(di), tolerance = 1e-8)
})

test_that("the dissimilarity index is invariant to rescaling a predictor", {
  tr <- mk_aoa_pts(80)
  nd <- mk_aoa_new(a = c(0, 1, 5), b = c(0, -1, 5))
  base <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))

  tr2 <- tr; tr2$a <- tr$a * 1000 + 50
  nd2 <- nd; nd2$a <- nd$a * 1000 + 50
  scaled <- area_of_applicability(nd2, train_sf = tr2,
                                  predictor_vars = c("a", "b"))
  expect_equal(base$aoa$DI, scaled$aoa$DI)
  expect_equal(base$threshold, scaled$threshold)
})

test_that("a shifted block of new data is not silently re-centred", {
  tr <- mk_aoa_pts(100)
  # Same internal spread as the training data, 10 sd away from it.
  nd <- mk_aoa_new(a = c(9.5, 10, 10.5), b = c(-0.5, 0, 0.5))
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  expect_false(any(res$aoa$AOA))
  expect_true(all(res$aoa$DI > res$threshold))
})

test_that("blocked folds widen the area of applicability", {
  # Two tight, well-separated clusters in predictor space.  Leave-one-out sees
  # a within-cluster neighbour a hair away; a fold that holds out a whole
  # cluster sees the other cluster.  The threshold must reflect that.
  set.seed(8)
  tr <- mk_aoa_pts(100, a = c(rnorm(50, 0, 0.1), rnorm(50, 10, 0.1)),
                   b = c(rnorm(50, 0, 0.1), rnorm(50, 10, 0.1)))
  nd <- mk_aoa_new(a = c(0, 5, 10), b = c(0, 5, 10))
  cl <- rep(1:2, each = 50)

  loo   <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  bloc  <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                 folds = cl)
  expect_gt(bloc$threshold, loo$threshold)
  expect_gte(bloc$n_inside, loo$n_inside)
  expect_true(bloc$params$folds_supplied)
  expect_identical(bloc$params$n_folds, 2L)
  expect_false(loo$params$folds_supplied)
})

test_that("area_of_applicability agrees between the FNN and dense paths", {
  skip_if_not_installed("FNN")
  tr <- mk_aoa_pts(150)
  nd <- mk_aoa_new(a = seq(-3, 6, length.out = 40),
                   b = seq(4, -3, length.out = 40))
  fast <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                use_fnn = TRUE)
  slow <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                use_fnn = FALSE)
  expect_equal(fast$aoa$DI, slow$aoa$DI, tolerance = 1e-8)
  expect_equal(fast$train_DI, slow$train_DI, tolerance = 1e-8)
  expect_equal(fast$threshold, slow$threshold, tolerance = 1e-8)
  expect_identical(fast$aoa$AOA, slow$aoa$AOA)
})

test_that("area_of_applicability takes training data from a fitted model", {
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(0, 15), b = c(0, 15))
  fit <- new_spatial_fit(subclass = "stub_fit", engine = NULL,
                         formula = stats::as.formula("z ~ a + b"),
                         response_var = "z",
                         predictor_vars = c("a", "b"), data_sf = tr)
  res <- area_of_applicability(nd, model = fit)
  expect_identical(res$predictor_vars, c("a", "b"))
  expect_identical(res$n_train, 60L)
  expect_true(res$aoa$AOA[1])
  expect_false(res$aoa$AOA[2])
  expect_error(area_of_applicability(nd, model = list(a = 1)),
               "must be a spatial_fit")
})

test_that("missing predictors in newdata give NA rather than a number", {
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(0, NA, 15), b = c(0, 1, 15))
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  expect_true(is.na(res$aoa$DI[2]))
  expect_true(is.na(res$aoa$AOA[2]))
  expect_false(is.na(res$aoa$DI[1]))
  expect_identical(res$n_na, 1L)
})

test_that("zero-variance predictors are dropped and reported", {
  tr <- mk_aoa_pts(60); tr$const <- 5
  nd <- mk_aoa_new(a = c(0, 15), b = c(0, 15)); nd$const <- 5
  res <- area_of_applicability(nd, train_sf = tr,
                               predictor_vars = c("a", "b", "const"))
  expect_identical(res$dropped_vars, "const")
  expect_identical(res$predictor_vars, c("a", "b"))
})

test_that("area_of_applicability accepts an explicit threshold", {
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(0, 15), b = c(0, 15))
  auto  <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  fixed <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                 threshold = 100)
  expect_equal(fixed$threshold, 100)
  expect_true(all(fixed$aoa$AOA))
  expect_true(fixed$params$threshold_supplied)
  expect_false(auto$params$threshold_supplied)
  expect_error(area_of_applicability(nd, train_sf = tr,
                                     predictor_vars = c("a", "b"),
                                     threshold = -1),
               "single positive number")
})

test_that("area_of_applicability validates its inputs", {
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(0, 15), b = c(0, 15))
  expect_error(area_of_applicability(nd), "supply .model., or both")
  expect_error(area_of_applicability(nd[0, ], train_sf = tr,
                                     predictor_vars = c("a", "b")),
               "has no rows")
  expect_error(area_of_applicability(nd, train_sf = tr[1, ],
                                     predictor_vars = c("a", "b")),
               "at least two training rows")
  expect_error(area_of_applicability(nd, train_sf = tr,
                                     predictor_vars = c("a", "b"),
                                     folds = c(1, 2)),
               "2 labels but the training data has 60")
})

test_that("print.aoa reports the reference scheme and the extrapolation caveat", {
  tr <- mk_aoa_pts(80)
  nd <- mk_aoa_new(a = c(0, 20), b = c(0, 20))
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  txt <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(txt, "Area of applicability")
  expect_match(txt, "no folds supplied")
  expect_match(txt, "outlier-removed max")
  expect_match(txt, "inside the AOA")
  expect_match(txt, "extrapolations")

  res2 <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                folds = rep(1:4, length.out = 80))
  txt2 <- paste(utils::capture.output(print(res2)), collapse = "\n")
  expect_match(txt2, "outside each of 4 CV folds")

  invisible(utils::capture.output(vis <- withVisible(print(res))))
  expect_false(vis$visible)
})


# --- regressions found in the full-package audit -----------------------------

test_that("a non-finite training value drops the ROW, not the predictor", {
  # complete.cases() passes Inf through; it then poisons that predictor's mean
  # and sd, .aoa_scaling() sees a non-finite scale and drops the whole
  # predictor, and the warning claims it has "no variance" -- a false diagnosis
  # that silently removes a dimension the index is supposed to measure.
  tr <- mk_aoa_pts(60)
  tr$a[1] <- Inf
  nd <- mk_aoa_new(a = c(0, 1e6), b = c(0, 0))
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))

  expect_length(res$dropped_vars, 0L)          # `a` must survive
  expect_identical(res$predictor_vars, c("a", "b"))
  expect_identical(res$n_train, 59L)           # the Inf row is gone
  expect_false(res$aoa$AOA[2])                 # a = 1e6 is nowhere near training
})

test_that("a fold with a row in both train and test is rejected", {
  # Such a row is its own nearest neighbour at distance zero, so every training
  # DI collapses to 0, the threshold collapses to 0, and nearly everything is
  # reported outside the AOA. Silent, and in the wrong direction.
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(0, 15), b = c(0, 15))
  expect_error(
    area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                          folds = list(list(test = 1:60, train = 1:60))),
    "both its train and test slots"
  )
})

test_that("unused factor levels do not become empty folds", {
  # A factor subset from a larger data set keeps its levels; each unused one
  # would otherwise be counted as a fold and reach .aoa_min_dist() with no
  # test rows.
  tr <- mk_aoa_pts(60)
  nd <- mk_aoa_new(a = c(0, 15), b = c(0, 15))
  f  <- factor(rep(1:2, each = 30), levels = 1:3)
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                               folds = f)
  expect_identical(res$params$n_folds, 2L)
  txt <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(txt, "outside each of 2 CV folds")

  # A factor with only one non-empty level is not two folds, whatever its
  # levels attribute says.
  f1 <- factor(rep("a", 60), levels = c("a", "b"))
  expect_error(
    area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                          folds = f1),
    "at least two"
  )
})

test_that("the duplicate-test-fold error names the fold the row was first seen in", {
  expect_error(
    .aoa_fold_splits(list(list(test = 1:3, train = 4:5),
                          list(test = 3:5, train = 1:2)), 5L),
    "both fold 1 and fold 2"
  )
})

test_that("a coordinate-using model is measured in coordinate space too", {
  skip_if_not_installed("ranger")
  # A forest fitted with include_coords = TRUE splits on location, so the DI
  # has to include location. Otherwise a point far outside the training extent
  # with ordinary covariate values reads as INSIDE the AOA -- the exact
  # extrapolation the index exists to catch.
  set.seed(12)
  n <- 150
  tr <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  tr$z <- 2 * tr$a + rnorm(n, 0, 0.3)
  fit <- fit_rf_model(tr, "z", "a", num_trees = 100L, include_coords = TRUE)

  # Same covariate value; one point inside the extent, one 500 km away.
  nd <- sf::st_as_sf(
    data.frame(x = c(500, 500000), y = c(500, 500000), a = c(0, 0)),
    coords = c("x", "y"), crs = 32632
  )
  res <- area_of_applicability(nd, model = fit)

  expect_true(all(c("..x", "..y") %in% res$predictor_vars))
  expect_true(res$aoa$AOA[1])
  expect_false(res$aoa$AOA[2])
  expect_gt(res$aoa$DI[2], res$aoa$DI[1])
})

test_that("a coordinate-using model reprojects newdata before measuring it", {
  skip_if_not_installed("ranger")
  set.seed(13)
  n <- 120
  tr <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), a = rnorm(n)),
    coords = c("x", "y"), crs = 32632
  )
  tr$z <- 2 * tr$a + rnorm(n, 0, 0.3)
  fit <- fit_rf_model(tr, "z", "a", num_trees = 100L, include_coords = TRUE)

  nd_proj <- sf::st_as_sf(data.frame(x = c(400, 600), y = c(400, 600),
                                     a = c(0, 0.5)),
                          coords = c("x", "y"), crs = 32632)
  nd_ll <- sf::st_transform(nd_proj, 4326)

  # Degrees compared against metres would put both points absurdly far out.
  expect_equal(area_of_applicability(nd_ll, model = fit)$aoa$DI,
               area_of_applicability(nd_proj, model = fit)$aoa$DI,
               tolerance = 1e-6)
})
