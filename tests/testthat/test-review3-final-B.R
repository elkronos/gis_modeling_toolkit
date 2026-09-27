# tests/testthat/test-review3-final-B.R
# ---------------------------------------------------------------------------
# Third review round, final pass B: fold provenance in fold_separation() and
# kriging_adequacy(), non-finite rows and the reasons for an early NA in
# estimate_sac_range(), fold_separation()'s print(), and the NaN-importance
# note in fit_rf_model().  Each test fails on the code before the fix,
# except the one that pins the polygon workflow the new provenance check must
# leave working.  Helpers are in helper-review2-folds.R.
# ---------------------------------------------------------------------------

rfb_collect_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

# A small exponential field with a smooth covariate, so both the OLS and the
# REML detrending have something to remove.
rfb_field <- function(n = 80, seed = 3) {
  set.seed(seed)
  x <- 5e5 + runif(n, 0, 1000); y <- 5e6 + runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  cv <- as.numeric(t(chol(exp(-d / 300) + diag(1e-6, n))) %*% rnorm(n))
  e  <- as.numeric(t(chol(exp(-d / 100) + diag(0.1, n))) %*% rnorm(n))
  r2_pts(x, y, cv = cv, z = 1 + 2 * cv + e)
}


# ---- S11-contracts-2: folds from another layer are refused -----------------

test_that("fold_separation() refuses folds whose rows sit elsewhere in data_sf", {
  set.seed(1)
  p <- r2_pts(runif(60, 0, 1000), runif(60, 0, 1000))
  f <- make_folds(p, k = 3, method = "random_kfold", seed = 1)
  # One row fewer: by position every later row is its neighbour's.
  expect_error(fold_separation(f, p[-5, ], sac = 200),
               "fold_separation\\(\\): the supplied `folds` were built from different data")
  # The layer the folds were built on, in another CRS, is the same layer.
  expect_s3_class(fold_separation(f, sf::st_transform(p, 3857)), "fold_separation")
  # A bare split list carries no probe and is measured as before.
  expect_s3_class(fold_separation(f$folds, p[-5, ]), "fold_separation")
})

test_that("fold_separation() skips the location check across POINT and non-POINT", {
  # Triangles: the point coerce_to_points() gives one is not its centroid,
  # which is what the probe records, so without the skip the pointized copy
  # of the layer the folds were built on was refused.
  tri <- sf::st_sf(geometry = sf::st_sfc(unlist(lapply(0:5, function(i) lapply(0:5, function(j) {
    x <- 20 * i; y <- 20 * j
    sf::st_polygon(list(rbind(c(x, y), c(x + 10, y), c(x, y + 10), c(x, y))))
  })), recursive = FALSE), crs = 32632))
  f <- make_folds(tri, k = 3, method = "random_kfold", seed = 1)
  pt <- coerce_to_points(tri, "auto")
  s_pt <- fold_separation(f, pt)
  expect_equal(s_pt$min_dist, fold_separation(f, tri)$min_dist)
})

test_that("kriging_adequacy() refuses folds built before assignment dropped points", {
  skip_if_not_installed("gstat")
  set.seed(2)
  n <- 80
  pts <- r2_pts(5e5 + runif(n, 0, 1000), 5e6 + runif(n, 0, 1000), z = rnorm(n))
  bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(
    c(xmin = 5e5 + 100, ymin = 5e6 + 100, xmax = 5e5 + 1000, ymax = 5e6 + 1000),
    crs = 32632)))
  cells <- create_grid_polygons(bnd, target_cells = 9, type = "square")
  asg <- suppressWarnings(assign_features_to_polygons(pts, cells))
  expect_lt(nrow(asg), nrow(pts))
  f <- make_folds(pts, k = 3, method = "random_kfold", seed = 1)
  expect_error(kriging_adequacy(asg, "z", cells, folds = f),
               "kriging_adequacy\\(\\): the supplied `folds` were built from different data")
})


# ---- S11-contracts-4: non-finite rows are left out of the detrending -------

test_that("one Inf predictor or response leaves the row out instead of the detrending", {
  skip_if_not_installed("gstat")
  p <- rfb_field()
  p_na <- p;  p_na$cv[5] <- NA
  p_inf <- p; p_inf$cv[5] <- Inf
  p_yinf <- p; p_yinf$z[5] <- Inf
  r_na <- estimate_sac_range(p_na, "z", predictor_vars = "cv")
  expect_true(is.finite(r_na))
  expect_true(attr(r_na, "detrended"))
  expect_no_warning(r_inf <- estimate_sac_range(p_inf, "z", predictor_vars = "cv"))
  expect_true(attr(r_inf, "detrended"))
  expect_identical(attr(r_inf, "detrend_method"), "ols")
  expect_equal(as.numeric(r_inf), as.numeric(r_na))
  expect_no_warning(r_yinf <- estimate_sac_range(p_yinf, "z", predictor_vars = "cv"))
  expect_true(attr(r_yinf, "detrended"))
  expect_equal(as.numeric(r_yinf), as.numeric(r_na))
})

test_that("a -Inf predictor under detrend = 'reml' is left out, not 'did not converge'", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  p <- rfb_field()
  p_na <- p;  p_na$cv[5] <- NA
  p_inf <- p; p_inf$cv[5] <- -Inf
  r_na <- estimate_sac_range(p_na, "z", predictor_vars = "cv", detrend = "reml")
  expect_identical(attr(r_na, "detrend_method"), "reml")
  expect_no_warning(r_inf <- estimate_sac_range(p_inf, "z", predictor_vars = "cv",
                                                detrend = "reml"))
  expect_identical(attr(r_inf, "detrend_method"), "reml")
  expect_equal(as.numeric(r_inf), as.numeric(r_na))
})


# ---- S5-FOLDS-9: an early NA says why ---------------------------------------

test_that("estimate_sac_range()'s early NA returns carry a rejected_reason and nothing else", {
  skip_if_not_installed("gstat")
  reason_of <- function(r) {
    r <- r2_quiet(r)
    expect_true(is.na(r))
    expect_false(inherits(r, "sac_range"))
    expect_identical(names(attributes(r)), "rejected_reason")
    attr(r, "rejected_reason")
  }
  set.seed(1)
  p20 <- r2_pts(5e5 + runif(20, 0, 1000), 5e6 + runif(20, 0, 1000), z = rnorm(20))
  expect_identical(reason_of(estimate_sac_range(p20, "z")),
                   "20 points, fewer than the 30 a variogram range is estimated from")

  p <- rfb_field(n = 60)
  pc <- p; pc$z <- 3
  expect_identical(reason_of(estimate_sac_range(pc, "z")), "the response is constant")
  pe <- p; pe$z <- 2 * pe$cv + 1                 # the predictor explains it exactly
  expect_match(reason_of(estimate_sac_range(pe, "z", predictor_vars = "cv")),
               "^the residuals on predictor_vars are constant")
  pf <- p; pf$z[1:40] <- NA
  expect_identical(reason_of(estimate_sac_range(pf, "z")),
                   paste0("20 point(s) with a finite value to model, fewer than ",
                          "the 30 a variogram range is estimated from"))
  ps <- r2_pts(rep(5e5, 40), rep(5e6, 40), z = rnorm(40))
  expect_match(reason_of(estimate_sac_range(ps, "z")), "^the points have no extent")

  # make_folds(auto_range = TRUE) names it in its warning.
  res <- rfb_collect_warnings(r2_quiet(make_folds(
    pf, k = 3, method = "block_kfold", auto_range = TRUE, response_var = "z", seed = 1)))
  expect_true(any(grepl("no autocorrelation range was identified (20 point(s) with a finite value",
                        res$warnings, fixed = TRUE)))
})


# ---- S11-contracts-7: the verdict fits the fold scheme ----------------------

test_that("fold_separation()'s verdict names a remedy the fold scheme has", {
  adv <- spatialkit:::.fold_separation_advice
  expect_match(adv("block_kfold", 0.9), "optimistic: widen the blocks\\.$")
  expect_match(adv("block_kfold", 0.3), "contiguous blocks at the edges")
  expect_match(adv("buffered_loo", 0.9), "widen the buffer")
  for (m in c("random_kfold", "leave_location_out", "supplied splits", NA)) {
    expect_no_match(adv(m, 0.9), "blocks")
    expect_match(adv(m, 0.9), "use blocked or buffered folds")
    expect_no_match(adv(m, 0.3), "contiguous blocks")
  }
  expect_no_match(adv("nndm", 0.9), "optimistic")
  expect_match(adv("nndm", 0.9), "prediction points")
  expect_match(adv("nndm", 0.05), "Little of the hold-out")

  set.seed(1)
  p <- r2_pts(runif(60, 0, 1000), runif(60, 0, 1000))
  flatten <- function(x) gsub("\\s+", " ", paste(x, collapse = " "))
  rnd <- make_folds(p, k = 3, method = "random_kfold", seed = 1)
  out <- flatten(utils::capture.output(print(fold_separation(rnd, p, sac = 500))))
  expect_match(out, "closer to a training point")
  expect_no_match(out, "widen the blocks")
  grid <- sf::st_as_sf(sf::st_make_grid(p, n = c(8, 8), what = "centers"))
  nn <- r2_quiet(suppressWarnings(make_folds(p, k = 5, method = "nndm",
                                             prediction_points = grid)))
  out_nn <- flatten(utils::capture.output(print(fold_separation(nn, p, sac = 500))))
  expect_match(out_nn, "closer to a training point")
  expect_no_match(out_nn, "optimistic")
  expect_no_match(out_nn, "widen the blocks")
})


# ---- fit_rf_model(): the all-NaN warning names the AOA weights ------------

test_that("a forest with no out-of-bag rows says its importance cannot weight the AOA", {
  skip_if_not_installed("ranger")
  set.seed(1)
  d <- r2_pts(runif(40, 0, 1000), runif(40, 0, 1000), a = rnorm(40), b = rnorm(40))
  d$z <- d$a + rnorm(40)
  res <- rfb_collect_warnings(r2_quiet(suppressMessages(
    fit_rf_model(d, "z", c("a", "b"), num_trees = 20, seed = 1,
                 replace = FALSE, sample_fraction = 1))))
  expect_true(any(grepl(paste0("no row is out of bag.*area_of_applicability\\(\\) ",
                               "cannot be weighted by that importance either: ",
                               "pass weights = NULL"), res$warnings)))
})
