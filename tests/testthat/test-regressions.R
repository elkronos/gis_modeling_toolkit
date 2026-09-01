# ===========================================================================
# Regression guards for fixes whose failure mode is a WRONG NUMBER rather than
# an error, and which have no natural home in a single-function test file.
#
# Every one of these was silent: the call returned, the shape was right, and
# only the values were off.
# ===========================================================================

.rg_points <- function(n = 40, seed = 5) {
  set.seed(seed)
  out <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               z = rnorm(n), w = rnorm(n)),
    coords = c("x", "y"), crs = 3857)
  out$..row_id <- seq_len(n)
  out
}


# ---------------------------------------------------------------------------
# make_folds(nndm, seed = NULL) must not re-seed from the clock
# ---------------------------------------------------------------------------

test_that("nndm with the default seed leaves the caller's RNG deterministic", {
  # `set.seed(NULL)` does NOT mean "leave the RNG alone": it RE-INITIALISES the
  # generator from the clock and process ID.  The ordinary call
  # make_folds(pts, method = "nndm", prediction_points = grid) therefore used
  # to destroy the caller's stream, so nothing downstream of it was
  # reproducible -- and the folds themselves changed run to run.
  set.seed(3)
  pts  <- sf::st_as_sf(data.frame(x = runif(30, 0, 1000), y = runif(30, 0, 1000)),
                       coords = c("x", "y"), crs = 3857)
  grid <- sf::st_as_sf(data.frame(x = runif(10, 0, 1000), y = runif(10, 0, 1000)),
                       coords = c("x", "y"), crs = 3857)

  # Unseeded means unseeded: the stream ADVANCES, by exactly the n uniforms the
  # radius draw consumes, and by nothing else.
  set.seed(999); invisible(stats::runif(nrow(pts))); expected <- stats::runif(3)
  set.seed(999)
  f1 <- suppressMessages(make_folds(pts, method = "nndm",
                                    prediction_points = grid))
  after1 <- stats::runif(3)
  expect_equal(after1, expected)

  # Same caller seed, same folds and same continuation -- which clock seeding
  # cannot deliver.
  set.seed(999)
  f2 <- suppressMessages(make_folds(pts, method = "nndm",
                                    prediction_points = grid))
  after2 <- stats::runif(3)
  expect_equal(f1$folds, f2$folds)
  expect_equal(after1, after2)

  # A different caller seed really does give different folds, so the equality
  # above is not satisfiable by a deterministic implementation.
  set.seed(1000)
  f3 <- suppressMessages(make_folds(pts, method = "nndm",
                                    prediction_points = grid))
  expect_false(isTRUE(all.equal(f1$folds, f3$folds)))

  # An explicit seed is independent of the caller's stream, as before.
  set.seed(1); a <- suppressMessages(make_folds(pts, method = "nndm",
                                                prediction_points = grid,
                                                seed = 42))
  set.seed(2); b <- suppressMessages(make_folds(pts, method = "nndm",
                                                prediction_points = grid,
                                                seed = 42))
  expect_equal(a$folds, b$folds)
})


# ---------------------------------------------------------------------------
# .remap_folds() must not renumber the folds that survive
# ---------------------------------------------------------------------------

test_that(".remap_folds keeps the original fold labels when one is dropped", {
  # Dropping an unusable fold removes an element from the list.  Renumbering
  # the survivors made fold_metrics$fold and predictions$fold disagree with
  # make_folds()$assignment$fold, so "fold 2 scored badly" pointed at the
  # wrong block on the map.
  pts   <- .rg_points(40)
  folds <- make_folds(pts, k = 4, method = "random_kfold", seed = 7)
  expect_length(folds$folds, 4L)

  dropped  <- 2L
  keep_idx <- setdiff(seq_len(40), folds$folds[[dropped]]$test)

  lines <- capture_spatialkit_log(
    remapped <- spatialkit:::.remap_folds(folds, keep_idx, k = 4, seed = 7))
  expect_true(log_has(lines, "empty test sets after remapping"))

  expect_length(remapped, 3L)
  # The surviving labels are 1, 3, 4 -- NOT 1, 2, 3.
  expect_equal(vapply(remapped, function(z) z$fold_id, integer(1)),
               c(1L, 3L, 4L))
  # And each carries the rows the ORIGINAL fold of that number carried.
  for (r in remapped)
    expect_setequal(r$test, folds$folds[[r$fold_id]]$test)
})


test_that(".remap_folds drops folds left with fewer than two training rows", {
  # An empty (or one-row) TRAINING set is as fatal as an empty test set: there
  # is nothing to fit on.  Such folds used to reach .cv_fit_one_fold(), fail
  # there one at a time, and surface only as a generic "all folds failed".
  remap <- spatialkit:::.remap_folds

  # keep_idx admits only rows 1:6, so fold 2's training set shrinks to a single
  # surviving row and fold 3's to none, while fold 1 keeps four.
  folds <- list(
    list(train = 3:10, test = 1:2),     # -> train 3:6  (kept)
    list(train = c(1L, 20L), test = 5L),# -> train 1    (dropped: < 2)
    list(train = 30:40,  test = 6L),    # -> train none (dropped)
    list(train = 1:5,    test = 6L)     # -> train 1:5  (kept)
  )
  keep <- 1:6

  lines <- capture_spatialkit_log(out <- remap(folds, keep))
  expect_true(log_has(lines, "fewer than 2 training rows"))
  expect_true(log_has(lines, "2 fold\\(s\\)"))

  expect_length(out, 2L)
  # The survivors keep their ORIGINAL labels, so the drop is not a renumbering.
  expect_equal(vapply(out, function(z) z$fold_id, integer(1)), c(1L, 4L))
  expect_equal(out[[1]]$train, 3:6)
  expect_equal(out[[2]]$train, 1:5)
  # Every fold that came back really is fittable.
  expect_true(all(vapply(out, function(z) length(z$train) >= 2L, logical(1))))

  # A fold with exactly two training rows is on the right side of the line.
  edge <- remap(list(list(train = c(1L, 2L, 90L), test = 3L)), keep)
  expect_length(edge, 1L)
  expect_equal(edge[[1]]$train, 1:2)
})


# ---------------------------------------------------------------------------
# Pooled CV R-squared must use each observation's TRAINING-fold mean
# ---------------------------------------------------------------------------

test_that("cv R-squared is scored against the training-fold mean, not the pooled one", {
  # Out-of-sample R^2 measured against the test data's OWN mean is the classic
  # flattering number: it credits the model for knowing where the held-out
  # block sits, which at prediction time it does not.  cv_spatial() therefore
  # carries a per-observation y_train_mean through to .compute_reg_metrics().
  #
  # Three explicit west-to-east strips over a steep x-gradient, so the two
  # baselines are genuinely different: each strip's training mean is the mean
  # of the OTHER two strips, which sits far from both the strip's own values
  # and the pooled mean of everything held out.  The only predictor is noise,
  # so the fitted model predicts close to its training mean and the choice of
  # baseline is the whole story.
  set.seed(11)
  n   <- 90
  x   <- runif(n, 0, 900); y <- runif(n, 0, 900)
  pts <- sf::st_as_sf(
    data.frame(x = x, y = y, w = rnorm(n), z = 0.1 * x + rnorm(n, 0, 2)),
    coords = c("x", "y"), crs = 3857)

  strip <- cut(x, breaks = c(-Inf, 300, 600, Inf), labels = FALSE)
  folds <- lapply(1:3, function(j) list(train = which(strip != j),
                                        test  = which(strip == j)))
  expect_true(all(vapply(folds, function(f) length(f$test) >= 10L, logical(1))))

  cv <- cv_spatial(pts, "z", "w",
                   fit_fn = function(tr) lm_spatial_fit(tr, "z", "w"),
                   folds = folds, seed = 1)

  p <- cv$predictions
  expect_true("y_train_mean" %in% names(p))
  expect_equal(cv$n_folds_succeeded, 3L)

  # Every prediction carries the mean of ITS OWN fold's training rows.
  for (j in seq_along(folds)) {
    in_fold <- p$fold == j
    expect_equal(sum(in_fold), length(folds[[j]]$test), info = paste("fold", j))
    expect_equal(unique(p$y_train_mean[in_fold]),
                 mean(pts$z[folds[[j]]$train]), tolerance = 1e-10,
                 info = paste("fold", j))
  }
  # The three baselines really are three different numbers.
  expect_equal(length(unique(p$y_train_mean)), 3L)

  ok  <- is.finite(p$y) & is.finite(p$yhat)
  rss <- sum((p$y[ok] - p$yhat[ok])^2)
  r2_train  <- 1 - rss / sum((p$y[ok] - p$y_train_mean[ok])^2)
  r2_pooled <- 1 - rss / sum((p$y[ok] - mean(p$y[ok]))^2)

  # The two baselines are a full R^2 unit apart, so the assertion below picks
  # one of them rather than tolerating either.
  expect_gt(abs(r2_train - r2_pooled), 0.5)
  expect_equal(cv$overall$R2, r2_train, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(cv$overall$R2, r2_pooled, tolerance = 1e-3)))

  # RMSE has no baseline, so it must be unaffected either way -- a check that
  # the fixture is not simply broken.
  expect_equal(cv$overall$RMSE, sqrt(rss / sum(ok)), tolerance = 1e-10)
})


test_that("cv_spatial reports fold labels that match make_folds()", {
  # The user-visible form of the same fix.
  pts   <- .rg_points(40)
  folds <- make_folds(pts, k = 4, method = "random_kfold", seed = 7)

  # Make every row of fold 2 unusable, so prep_model_data() drops them and
  # fold 2 has nothing left to test on.
  holed <- pts
  holed$w[folds$folds[[2]]$test] <- NA

  res <- suppressWarnings(suppressMessages(
    cv_spatial(holed, "z", "w",
               fit_fn = function(tr) lm_spatial_fit(tr, "z", "w"),
               folds = folds, seed = 1)))

  expect_equal(sort(unique(res$fold_metrics$fold)), c(1L, 3L, 4L))
  expect_equal(sort(unique(res$predictions$fold)), c(1L, 3L, 4L))
  expect_equal(res$n_folds_attempted, 3L)
  expect_equal(res$n_folds_succeeded, 3L)
  # The labels are a subset of the fold numbers make_folds() assigned.
  expect_true(all(res$fold_metrics$fold %in% unique(folds$assignment$fold)))
  # Each reported fold scored the rows make_folds() put in it.
  for (j in unique(res$fold_metrics$fold)) {
    reported <- sort(res$predictions$`..row_id`[res$predictions$fold == j])
    expect_setequal(reported, intersect(folds$folds[[j]]$test,
                                        which(!is.na(holed$w))))
  }
})


# ---------------------------------------------------------------------------
# area_of_applicability() resolves ..row_id values, it does not index with them
# ---------------------------------------------------------------------------

test_that("area_of_applicability is unchanged when ..row_id is permuted", {
  # Every make_folds() branch emits ..row_id VALUES in its train/test slots.
  # Indexing the training matrix with them computed the reference distances for
  # the WRONG rows whenever the IDs happened to land inside 1:n -- which is
  # exactly what a permutation of 1:n does -- so the threshold shifted
  # silently.  Resolving via match() makes the answer depend on fold
  # MEMBERSHIP, not on the numbers used to name the rows.
  set.seed(11)
  n  <- 60
  tr <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n)),
    coords = c("x", "y"), crs = 3857)
  nd <- sf::st_as_sf(
    data.frame(x = runif(40, 0, 1000), y = runif(40, 0, 1000),
               a = rnorm(40), b = rnorm(40)),
    coords = c("x", "y"), crs = 3857)

  # Fold membership defined by POSITION.
  by_position <- lapply(1:3, function(i) {
    te <- which((seq_len(n) - 1L) %% 3L + 1L == i)
    list(test = te, train = setdiff(seq_len(n), te))
  })

  base <- area_of_applicability(nd, train_sf = tr,
                                predictor_vars = c("a", "b"),
                                folds = by_position)
  expect_true(is.finite(base$threshold))

  # The SAME membership, expressed through a permuted ..row_id.
  set.seed(6)
  perm     <- sample(seq_len(n))
  tr_perm  <- tr
  tr_perm$..row_id <- perm
  by_id <- lapply(by_position, function(f)
    list(test = perm[f$test], train = perm[f$train]))

  permuted <- area_of_applicability(nd, train_sf = tr_perm,
                                    predictor_vars = c("a", "b"),
                                    folds = by_id)

  expect_equal(permuted$threshold, base$threshold)
  expect_equal(permuted$train_DI, base$train_DI)
  expect_equal(permuted$aoa$DI, base$aoa$DI)
  expect_equal(permuted$n_inside, base$n_inside)

  # The discriminating half: treating those id VALUES as POSITIONS -- the old
  # behaviour -- assigns different rows to each fold and gives a different
  # threshold, so the equalities above are not vacuous.
  as_positions <- area_of_applicability(nd, train_sf = tr,
                                        predictor_vars = c("a", "b"),
                                        folds = by_id)
  expect_false(isTRUE(all.equal(as_positions$threshold, base$threshold)))

  # An id that names no training row is an error, not a silent NA fold.
  bad <- lapply(by_id, function(f) f)
  bad[[1]]$test <- c(bad[[1]]$test, 999999L)
  expect_error(area_of_applicability(nd, train_sf = tr_perm,
                                     predictor_vars = c("a", "b"), folds = bad),
               "not in the training data")
})


# ---------------------------------------------------------------------------
# select_features_forward() scores every candidate on the same rows
# ---------------------------------------------------------------------------

test_that("select_features_forward drops non-finite rows, not just incomplete ones", {
  # complete.cases() alone lets Inf through, so a candidate carrying a single
  # Inf reached cv_spatial() with a SMALLER surviving row set than its rivals
  # -- and could then be preferred for having an easier subset rather than for
  # predicting better.  The pre-filter now matches prep_model_data() exactly.
  set.seed(1)
  n <- 120
  dat <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
               a = rnorm(n), b = rnorm(n), c1 = rnorm(n)),
    coords = c("x", "y"), crs = 3857)
  dat$z <- 3 * dat$a + 2 * dat$b + rnorm(n, 0, 0.3)

  holed <- dat
  holed$c1[c(4L, 11L)] <- Inf     # non-finite: complete.cases() says complete
  holed$b[7L]          <- NA      # incomplete: complete.cases() catches it

  fit_fn <- function(tr, vars) lm_spatial_fit(tr, "z", vars)

  lines <- capture_spatialkit_log(
    sel <- select_features_forward(holed, "z", c("a", "b", "c1"), fit_fn,
                                   k = 4, seed = 1, quiet = TRUE))

  # All THREE rows are dropped -- 2 non-finite plus 1 missing -- not just the
  # one complete.cases() sees.
  expect_true(log_has(lines, "dropping 3 row\\(s\\) incomplete or non-finite"))
  expect_true(log_has(lines, "scored on the same 117 observations"))

  # And the sweep is exactly the one you get by removing those rows yourself,
  # which is what "every candidate set is scored on the same observations"
  # means.
  manual <- select_features_forward(holed[-c(4L, 7L, 11L), ], "z",
                                    c("a", "b", "c1"), fit_fn, k = 4, seed = 1,
                                    quiet = TRUE)
  expect_equal(sel$selected, manual$selected)
  expect_equal(sel$score, manual$score)
  expect_equal(sel$history, manual$history)

  # Clean data logs nothing and keeps every row.
  clean_lines <- capture_spatialkit_log(
    clean <- select_features_forward(dat, "z", c("a", "b", "c1"), fit_fn,
                                     k = 4, seed = 1, quiet = TRUE))
  expect_false(log_has(clean_lines, "incomplete or non-finite"))
  expect_true(all(c("a", "b") %in% clean$selected))

  # Too few usable rows is an error naming the cause, not a crash inside CV.
  wrecked <- dat
  wrecked$c1[2:n] <- Inf          # leaves exactly one usable row
  expect_error(select_features_forward(wrecked, "z", c("a", "c1"), fit_fn,
                                       k = 4, seed = 1, quiet = TRUE),
               "fewer than two rows are complete and finite")
})


# ---------------------------------------------------------------------------
# ensure_projected()'s target_crs validation must not break CRS-less workflows
#
# ensure_projected() gained an error for a `target_crs` that does not resolve
# to a usable CRS, so that a typo cannot silently turn it into a no-op.  But
# internal callers derive the target from another object -- st_crs(training) --
# and that object is legitimately allowed to carry no CRS.  Passing NA_crs_
# straight through turned every CRS-less workflow into a hard error, and in
# cv_rf()'s case into a silent empty result (the per-fold predict() threw, so
# every fold "failed" and overall came back n_pred = 0, RMSE = NA).
# ---------------------------------------------------------------------------

.rg_nocrs <- function(n = 60, seed = 11) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100),
               z = rnorm(n), a = rnorm(n)),
    coords = c("x", "y"))          # deliberately no crs=
}

test_that("CRS-less data still flows through predict, folds and surfaces", {
  skip_if_not_installed("ranger")
  d <- .rg_nocrs()
  expect_true(is.na(sf::st_crs(d)))

  fit <- suppressMessages(fit_rf_model(d, "a", "z"))
  expect_length(suppressMessages(predict(fit, newdata = d)), nrow(d))

  expect_no_error(suppressMessages(
    make_folds(d, k = 3, method = "nndm", prediction_points = d)))
  expect_no_error(suppressMessages(
    make_folds(d, k = 3, method = "block_kfold",
               boundary = sf::st_as_sfc(sf::st_bbox(d)))))
  expect_no_error(suppressMessages(
    prep_model_data(d, "a", "z", boundary = sf::st_as_sfc(sf::st_bbox(d)))))
  expect_no_error(suppressMessages(
    predict_surface(fit, n_cells = 50, covariates = d)))
})

test_that("cv_rf on CRS-less data predicts every row instead of failing silently", {
  skip_if_not_installed("ranger")
  d <- .rg_nocrs()
  cv <- suppressMessages(cv_rf(d, "a", "z", k = 3))
  # The regression's signature was n_pred = 0 with RMSE = NA and a
  # "all folds failed" warning -- shape intact, content empty.
  expect_identical(cv$n_folds_succeeded, 3L)
  expect_identical(cv$overall$n_pred, nrow(d))
  expect_true(is.finite(cv$overall$RMSE))
})

test_that("a genuinely unusable target_crs is still rejected", {
  d <- .rg_nocrs(n = 5)
  expect_error(ensure_projected(d, target_crs = list(nonsense = TRUE)),
               "does not resolve to a usable CRS")
})


# ---------------------------------------------------------------------------
# build_tessellation(crs = ) on CRS-less input
#
# `if (!is.null(crs)) st_transform(x, crs)` fails with "cannot transform sfc
# object with missing crs" for exactly the users most likely to pass `crs =`:
# those whose data carries none.  Reprojection is impossible there, but
# assumption is not -- stamp the target and warn, as ensure_projected() does.
# ---------------------------------------------------------------------------

test_that("build_tessellation(crs=) stamps a CRS-less input rather than erroring", {
  d <- .rg_nocrs(n = 40)
  bb <- sf::st_as_sfc(sf::st_bbox(d))

  for (m in c("voronoi", "triangles")) {
    tess <- suppressMessages(build_tessellation(d, method = m, crs = 32632))
    expect_identical(sf::st_crs(tess$cells)$epsg, 32632L, label = m)
    expect_true("cell_id" %in% names(tess$cells), label = m)
  }
  for (m in c("hex", "square")) {
    tess <- suppressMessages(
      build_tessellation(d, boundary = bb, method = m, crs = 32632,
                         approx_n_cells = 16))
    expect_identical(sf::st_crs(tess$cells)$epsg, 32632L, label = m)
    expect_true("cell_id" %in% names(tess$cells), label = m)
  }
})

test_that("build_tessellation(crs=) still REPROJECTS input that has a CRS", {
  # The stamp path must not swallow the transform path: coordinates in a
  # different CRS have to actually move, not be relabelled.
  d <- .rg_nocrs(n = 40)
  sf::st_crs(d) <- 32632
  tess <- suppressMessages(build_tessellation(d, method = "voronoi", crs = 4326))
  expect_identical(sf::st_crs(tess$cells)$epsg, 4326L)
  # 32632 metres near the origin land at recognisable lon/lat, not at 0-100.
  cc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(tess$cells)))
  expect_true(all(abs(cc[, 1]) <= 180) && all(abs(cc[, 2]) <= 90))
  expect_false(isTRUE(all.equal(range(cc[, 1]), c(0, 100), tolerance = 1e-3)))
})


# ---------------------------------------------------------------------------
# compare_models_cv(): extras must REPLACE base arguments, not duplicate them
#
# base was built with c(list(...), rf_args), so any rf_args name colliding with
# a base name produced two entries of that name and do.call() died with
# "formal argument 'seed' matched by multiple actual arguments" -- reachable
# straight from the documented usage, since cv_rf() has both `k` and `seed`.
# ---------------------------------------------------------------------------

test_that("compare_models_cv accepts rf_args that collide with base arguments", {
  skip_if_not_installed("ranger")
  d <- .rg_points(n = 60)
  d$a <- d$z + 0.01 * sf::st_coordinates(d)[, 1]
  fo <- make_folds(d, k = 3, method = "random_kfold", seed = 1)

  for (extra in list(list(seed = 3), list(k = 3), list(num.trees = 50))) {
    res <- suppressWarnings(suppressMessages(
      compare_models_cv(d, "a", "z", models = "RF", folds = fo, rf_args = extra)))
    expect_true(is.data.frame(res$overall))
    expect_identical(res$overall$model, "RF")
  }
})

test_that("compare_models_cv refuses extras that would break comparability", {
  skip_if_not_installed("ranger")
  d <- .rg_points(n = 60)
  d$a <- d$z + 0.01 * sf::st_coordinates(d)[, 1]
  fo <- make_folds(d, k = 3, method = "random_kfold", seed = 1)

  expect_warning(
    res <- suppressMessages(
      compare_models_cv(d, "a", "z", models = "RF", folds = fo,
                        rf_args = list(folds = NULL))),
    "same data and the same folds")
  expect_identical(res$overall$model, "RF")   # dropped, not fatal
})


# ---------------------------------------------------------------------------
# compare_models() on a bare spatial_fit
#
# A spatial_fit *is* a list, so it passed the is.list() check and every loop
# then iterated the fit's own components as if they were models: all-NA Moran's
# I columns plus one "'fit' is not a spatial_fit object" log line per component.
# evaluate_insample() already wrapped a bare fit; compare_models() did not.
# ---------------------------------------------------------------------------

test_that("compare_models wraps a bare fit and reports real Moran's I", {
  skip_if_not_installed("ranger")
  d <- .rg_points(n = 80)
  d$a <- d$z + 0.02 * sf::st_coordinates(d)[, 1]
  fit <- suppressMessages(fit_rf_model(d, "a", "z"))

  bare <- suppressWarnings(suppressMessages(compare_models(fit)))
  expect_identical(nrow(bare), 1L)
  expect_identical(bare$model, "rf_fit")
  expect_true(is.finite(bare$resid_morans_I))
  expect_false(is.na(bare$resid_morans_p))

  # Identical to the explicitly-wrapped form apart from the row label.
  wrapped <- suppressWarnings(suppressMessages(compare_models(list(rf_fit = fit))))
  expect_equal(bare$resid_morans_I, wrapped$resid_morans_I)
  expect_equal(bare$RMSE, wrapped$RMSE)
})


# ---------------------------------------------------------------------------
# All four CV entry points must return the same shape
#
# The fold counts exist so that a run where EVERY fold failed is visible in
# the return value and not only in a warning -- `overall` is a well-formed
# all-NA row either way.  cv_spatial() and cv_rf() reported them; cv_gwr() and
# cv_bayes() built their return list by hand and omitted them, which is
# backwards: those are the two backends most likely to fail silently (Stan
# sampling, GWmodel bandwidth selection).
# ---------------------------------------------------------------------------

test_that("every cv_* wrapper returns the documented common fields", {
  common <- c("overall", "fold_metrics", "predictions", "folds",
              "n_folds_attempted", "n_folds_succeeded")

  # Checked statically for the backends that may not be installed, so this
  # test guards the contract on any machine rather than only where GWmodel
  # and brms happen to be present.
  for (fn in c("cv_gwr", "cv_bayes")) {
    src <- paste(deparse(body(get(fn, asNamespace("spatialkit")))),
                 collapse = " ")
    expect_true(grepl("n_folds_attempted", src, fixed = TRUE),
                label = paste(fn, "returns n_folds_attempted"))
    expect_true(grepl("n_folds_succeeded", src, fixed = TRUE),
                label = paste(fn, "returns n_folds_succeeded"))
  }

  skip_if_not_installed("ranger")
  d <- .rg_points(n = 60)
  d$a <- d$z + 0.01 * sf::st_coordinates(d)[, 1]
  fo <- make_folds(d, k = 3, method = "random_kfold", seed = 1)

  rf <- suppressMessages(cv_rf(d, "a", "z", folds = fo))
  sp <- suppressMessages(
    cv_spatial(d, "a", "z", folds = fo,
               fit_fn = function(tr, ...) fit_rf_model(tr, "a", "z")))

  for (nm in list(rf = rf, spatial = sp)) {
    expect_true(all(common %in% names(nm)))
    expect_identical(nm$n_folds_attempted, 3L)
    expect_identical(nm$n_folds_succeeded, 3L)
  }
})
