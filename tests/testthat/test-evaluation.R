# ===========================================================================
# R/evaluation.R: evaluate_insample(), compare_models(), compare_models_cv()
# and residual_morans_i().
#
# All four are exported and sold in DESCRIPTION as the "model comparison"
# half of the package; none had a test.  Everything except compare_models_cv()
# is reachable with helper-lmfit.R's lm_spatial_fit(), so no optional backend
# is needed for the bulk of it.
# ===========================================================================

# A stand-in gwr_fit.  Class order matters: "lmsurf_fit" first so the
# fitted()/residuals() methods registered in helper-lmfit.R still dispatch,
# "gwr_fit" second so compare_models()'s inherits() branch fires.  Constructing
# a real one needs GWmodel and a bandwidth search that succeeds by design --
# which is precisely the branch that cannot be reached that way.
.gwr_stub <- function(fit, bandwidth = 42, is_fallback = TRUE, AICc = 123.5) {
  fit$info$bandwidth             <- bandwidth
  fit$info$bandwidth_is_fallback <- is_fallback
  fit$info$AICc                  <- AICc
  class(fit) <- c("lmsurf_fit", "gwr_fit", "spatial_fit")
  fit
}


# ---------------------------------------------------------------------------
# evaluate_insample()
# ---------------------------------------------------------------------------

test_that("evaluate_insample returns one row per named fit", {
  pts <- surf_test_points(80)
  f_w <- lm_spatial_fit(pts, predictor_vars = "w")
  f_0 <- lm_spatial_fit(pts, predictor_vars = character(0))

  out <- evaluate_insample(list(with_w = f_w, intercept = f_0))

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_equal(out$model, c("with_w", "intercept"))
  expect_true(all(c("n", "RMSE", "MAE", "R2") %in% names(out)))
  expect_equal(out$n, c(80L, 80L))

  # The numbers are model_metrics() on each fit, not a shared placeholder.
  for (i in seq_len(2L)) {
    nm  <- out$model[i]
    ref <- model_metrics(list(with_w = f_w, intercept = f_0)[[nm]])
    expect_equal(out$RMSE[i], ref$RMSE)
    expect_equal(out$R2[i], ref$R2)
  }
  # An intercept-only fit explains nothing, by construction.
  expect_equal(out$R2[out$model == "intercept"], 0)
  expect_gte(out$R2[out$model == "with_w"], out$R2[out$model == "intercept"])
})


test_that("evaluate_insample accepts a bare spatial_fit and labels it by class", {
  fit <- lm_spatial_fit(surf_test_points(60), predictor_vars = "w")
  out <- evaluate_insample(fit)
  expect_equal(nrow(out), 1L)
  expect_equal(out$model, "lmsurf_fit")
  expect_equal(out$RMSE, model_metrics(fit)$RMSE)
})


test_that("evaluate_insample rejects an UNNAMED list instead of returning NULL", {
  # An unnamed list used to make lapply(names(fits), ...) iterate zero times,
  # so this returned NULL silently; compare_models() then died in
  # seq_len(nrow(NULL)) with a message naming nothing relevant.
  pts <- surf_test_points(40)
  fits <- list(lm_spatial_fit(pts, predictor_vars = "w"),
               lm_spatial_fit(pts, predictor_vars = character(0)))

  expect_error(evaluate_insample(fits), "must be a NAMED list")
  # Partially named is just as unusable.
  expect_error(evaluate_insample(stats::setNames(fits, c("a", ""))),
               "must be a NAMED list")
  expect_error(evaluate_insample(stats::setNames(fits, c("a", NA))),
               "must be a NAMED list")

  # ... and the earlier structural checks still fire first.
  expect_error(evaluate_insample(list()), "must be a spatial_fit or a named list")
  expect_error(evaluate_insample("not a fit"), "must be a spatial_fit or a named list")
})


test_that("evaluate_insample skips non-fits with a log line rather than erroring", {
  pts <- surf_test_points(50)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  lines <- capture_spatialkit_log(
    out <- evaluate_insample(list(good = fit, bad = list(a = 1)))
  )
  expect_true(log_has(lines, "'bad' is not a spatial_fit"))
  expect_equal(nrow(out), 1L)
  expect_equal(out$model, "good")
})


test_that("evaluate_insample scores newdata out of sample", {
  train <- surf_test_points(90, seed = 1)
  test  <- surf_test_points(40, seed = 77)
  fit   <- lm_spatial_fit(train, predictor_vars = "w")

  ins <- evaluate_insample(list(m = fit))
  oos <- evaluate_insample(list(m = fit), newdata = test)

  expect_equal(ins$n, 90L)
  expect_equal(oos$n, 40L)
  expect_equal(oos$RMSE, model_metrics(fit, newdata = test)$RMSE)
  # It really used newdata: the two metric sets differ.
  expect_false(isTRUE(all.equal(ins$RMSE, oos$RMSE)))
})


# ---------------------------------------------------------------------------
# compare_models()
# ---------------------------------------------------------------------------

test_that("compare_models adds information criteria and Moran columns", {
  pts <- surf_test_points(80)
  fits <- list(A = lm_spatial_fit(pts, predictor_vars = "w"),
               B = lm_spatial_fit(pts, predictor_vars = character(0)))

  cmp <- suppressWarnings(compare_models(fits))

  expect_equal(sort(cmp$model), c("A", "B"))
  expect_true(all(c("AICc", "LOOIC", "bandwidth_is_fallback",
                    "resid_morans_I", "resid_morans_z", "resid_morans_p") %in%
                    names(cmp)))
  # Neither backend supplies AICc or LOOIC, so both stay NA rather than 0.
  expect_true(all(is.na(cmp$AICc)))
  expect_true(all(is.na(cmp$LOOIC)))
  expect_true(all(is.na(cmp$bandwidth_is_fallback)))

  # The Moran columns are the ones residual_morans_i() reports for each fit.
  for (nm in cmp$model) {
    mi <- residual_morans_i(fits[[nm]])
    expect_equal(cmp$resid_morans_I[cmp$model == nm], mi$observed)
    expect_equal(cmp$resid_morans_z[cmp$model == nm], mi$z)
  }

  expect_error(compare_models(list()), "non-empty named list")
})


test_that("compare_models warns about a GWR fitted on a fallback bandwidth", {
  pts <- surf_test_points(60)
  ok   <- .gwr_stub(lm_spatial_fit(pts, predictor_vars = "w"),
                    bandwidth = 17, is_fallback = FALSE, AICc = 99)
  bad  <- .gwr_stub(lm_spatial_fit(pts, predictor_vars = "w"),
                    bandwidth = 42, is_fallback = TRUE, AICc = 123.5)

  expect_warning(compare_models(list(G = bad)), "used a fallback bandwidth")
  expect_warning(compare_models(list(G = bad)), "42")

  cmp <- suppressWarnings(compare_models(list(G = bad)))
  expect_true(cmp$bandwidth_is_fallback)
  expect_equal(cmp$AICc, 123.5)

  # A selected bandwidth must NOT warn, or the warning carries no information.
  expect_no_warning(cmp_ok <- compare_models(list(G = ok)))
  expect_false(cmp_ok$bandwidth_is_fallback)
  expect_equal(cmp_ok$AICc, 99)
})


test_that("compare_models flags residual spatial autocorrelation in the log", {
  # surf_test_points() carries a strong linear x/y trend; an intercept-only
  # fit leaves all of it in the residuals, so Moran's I is large and p tiny.
  pts <- surf_test_points(100)
  lines <- capture_spatialkit_log(
    cmp <- compare_models(list(flat = lm_spatial_fit(pts,
                                                     predictor_vars = character(0))))
  )
  expect_true(log_has(lines, "significant.*spatial autocorrelation"))
  expect_lt(cmp$resid_morans_p, 0.05)
  expect_gt(cmp$resid_morans_I, 0)
})


# ---------------------------------------------------------------------------
# residual_morans_i(): the sparse path
# ---------------------------------------------------------------------------

test_that("residual_morans_i agrees on the sparse and dense weight paths", {
  # The confirmed live bug: with FNN and Matrix installed the DEFAULT weight
  # matrix is a dgCMatrix, and base::crossprod() -- which the statistic and the
  # Cliff & Ord variance both used -- does not S4-dispatch on it, so the
  # default call raised "requires numeric/complex matrix/vector arguments".
  # The existing tests missed it: one passed a dense base matrix and the other
  # was gated on FNN being absent.
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")

  set.seed(19)
  n <- 120
  pts <- sf::st_as_sf(
    data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), w = rnorm(n)),
    coords = c("x", "y"), crs = 3857
  )
  pts$z <- 0.01 * sf::st_coordinates(pts)[, 1] + 0.5 * pts$w + rnorm(n, 0, 0.5)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  # Reproduce the weights residual_morans_i() builds for itself.
  coords <- sf::st_coordinates(
    ensure_projected(coerce_to_points(fit$data_sf, "auto")))[, 1:2, drop = FALSE]
  W_sparse <- spatialkit:::.build_knn_weights(coords, k = 8L,
                                              use_fnn = TRUE, use_matrix = TRUE)
  expect_true(inherits(W_sparse, "Matrix"))
  W_dense <- as.matrix(W_sparse)
  expect_false(is.matrix(W_sparse))            # the class that broke dispatch

  default  <- residual_morans_i(fit)                       # sparse internally
  explicit <- residual_morans_i(fit, weights = W_sparse)   # sparse, user-supplied
  dense    <- residual_morans_i(fit, weights = W_dense)

  expect_false(is.null(default))
  for (comp in c("observed", "expected", "sd", "z", "p_value")) {
    expect_equal(default[[comp]], dense[[comp]], tolerance = 1e-10)
    expect_equal(explicit[[comp]], dense[[comp]], tolerance = 1e-10)
  }
  expect_equal(default$n, n)

  # The statistic is not degenerate: a strong x-trend left in the residuals of
  # an intercept-only fit gives clear positive autocorrelation, and the fitted
  # model removes most of it.  Without this the equalities above would also
  # hold if every path returned NA.
  expect_true(is.finite(default$observed) && is.finite(default$z))
  flat <- residual_morans_i(lm_spatial_fit(pts, predictor_vars = character(0)))
  expect_gt(flat$observed, default$observed)
})


test_that("residual_morans_i refuses inputs it cannot score", {
  pts <- surf_test_points(40)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  # "Returns NULL with a warning" -- a real one, catchable by tryCatch().
  expect_warning(lines <- capture_spatialkit_log(out <- residual_morans_i(list(a = 1))),
                 "not a spatial_fit")
  expect_null(out)
  expect_true(log_has(lines, "not a spatial_fit"))

  # Constant residuals: centred they are all zero, so I is undefined.  The
  # guard has to test the CENTRED sum of squares -- a non-zero constant passes
  # sum(resid^2) and then makes the variance NaN.
  const <- fit
  registerS3method("residuals", "constresid_fit",
                   function(object, ...) rep(3, nrow(object$data_sf)))
  class(const) <- c("constresid_fit", class(fit))
  expect_warning(lines_c <- capture_spatialkit_log(out_c <- residual_morans_i(const)),
                 "zero residual variance")
  expect_null(out_c)
  expect_true(log_has(lines_c, "zero residual variance"))

  # Too few residuals to say anything.
  tiny <- lm_spatial_fit(surf_test_points(3), predictor_vars = character(0))
  expect_warning(lines_t <- capture_spatialkit_log(out_t <- residual_morans_i(tiny)),
                 "n < 4")
  expect_null(out_t)
  expect_true(log_has(lines_t, "n < 4"))
})


test_that("residual_morans_i refuses a wrongly shaped weights argument", {
  # It used to log a line and fall back to the default kNN(8) matrix, so a
  # user who supplied their own weights with one row too few -- or as a
  # data.frame -- got the statistic for a weight matrix they never chose,
  # with no R condition raised (measured: I = 0.874 for four malformed shapes
  # against 0.805 for the weights actually supplied).  An error names the
  # shape it saw.
  set.seed(23)
  pts <- surf_test_points(40)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  n   <- fit$n

  expect_error(residual_morans_i(fit, weights = matrix(1, 5, 5)), "n x n")
  expect_error(residual_morans_i(fit, weights = matrix(1, n - 1, n - 1)), "n x n")
  expect_error(residual_morans_i(fit, weights = rep(1, n)), "n x n")
  expect_error(residual_morans_i(fit, weights = matrix(1, n, n + 1)), "n x n")

  # A data.frame of the right shape is accepted as a matrix ...
  W  <- spatialkit:::.build_knn_weights(sf::st_coordinates(fit$data_sf), k = 4L,
                                        use_fnn = FALSE, use_matrix = FALSE)
  a  <- residual_morans_i(fit, weights = W)
  b  <- residual_morans_i(fit, weights = as.data.frame(W))
  expect_equal(a$observed, b$observed)
  # ... and gives the user's answer, not the default's.
  expect_false(isTRUE(all.equal(a$observed, residual_morans_i(fit)$observed)))
})


# ---------------------------------------------------------------------------
# compare_models_cv()
# ---------------------------------------------------------------------------

test_that("compare_models_cv cross-validates the RF branch", {
  skip_if_not_installed("ranger")
  pts <- surf_test_points(90, seed = 4)

  res <- compare_models_cv(pts, "z", "w", models = "RF", k = 3,
                           rf_args = list(num_trees = 60), quiet = TRUE)

  expect_named(res, c("overall", "by_fold", "rf_cv"))
  expect_equal(res$overall$model, "RF")
  expect_equal(nrow(res$overall), 1L)
  expect_true(is.finite(res$overall$RMSE))
  expect_equal(res$overall$n_pred, 90L)

  # by_fold carries one row per fold, all labelled RF.
  expect_equal(nrow(res$by_fold), 3L)
  expect_true(all(res$by_fold$model == "RF"))

  # The RF branch really ran cv_rf(), so the raw result is reachable.
  expect_true(is.list(res$rf_cv))
  expect_equal(res$overall$RMSE, res$rf_cv$overall$RMSE)
})


test_that("compare_models_cv forwards rf_args to the forest", {
  # A deterministic proof that the list reaches ranger: mtry cannot exceed the
  # predictor count, so an impossible value must surface rather than be
  # silently dropped the way a filtered argument list would drop it.
  skip_if_not_installed("ranger")
  pts <- surf_test_points(60, seed = 5)

  expect_warning(
    res <- compare_models_cv(pts, "z", "w", models = "RF", k = 3,
                             rf_args = list(mtry = 99), quiet = TRUE),
    "all folds failed"
  )
  expect_true(is.na(res$overall$RMSE))

  # ... and the same call without the bad argument succeeds, so the failure is
  # attributable to rf_args and not to the fixture.
  ok <- compare_models_cv(pts, "z", "w", models = "RF", k = 3,
                          rf_args = list(num_trees = 40), quiet = TRUE)
  expect_true(is.finite(ok$overall$RMSE))
})


test_that("compare_models_cv warns on unknown model names and errors when none remain", {
  skip_if_not_installed("ranger")
  pts <- surf_test_points(60, seed = 6)

  # Unknown names used to be dropped by a bare intersect(), so models = "RF"
  # silently ran GWR -- a wrong answer presented as the right one.
  expect_warning(
    res <- compare_models_cv(pts, "z", "w", models = c("RF", "XGBoost"),
                             k = 3, rf_args = list(num_trees = 40),
                             quiet = TRUE),
    "ignoring unknown model\\(s\\) 'XGBoost'"
  )
  expect_equal(res$overall$model, "RF")     # the recognised one still ran

  # Nothing recognised left: an error, not a silent fallback to the default.
  expect_error(
    suppressWarnings(compare_models_cv(pts, "z", "w", models = "XGBoost",
                                       k = 3, quiet = TRUE)),
    "no recognised model requested"
  )

  expect_error(compare_models_cv(sf::st_drop_geometry(pts), "z", "w"),
               "must be an sf object")
})


test_that("compare_models_cv drops a recognised model whose backend is absent", {
  skip_if_not_installed("ranger")
  absent <- c("GWR", "Bayesian")[
    !c(requireNamespace("GWmodel", quietly = TRUE) &&
         requireNamespace("sp", quietly = TRUE),
       requireNamespace("brms", quietly = TRUE))]
  skip_if(length(absent) == 0L, "both optional backends are installed")

  pts <- surf_test_points(60, seed = 8)

  # Asking for RF plus an unavailable backend runs RF alone, with a message.
  expect_message(
    res <- compare_models_cv(pts, "z", "w", models = c("RF", absent[1L]),
                             k = 3, rf_args = list(num_trees = 40)),
    "package/function unavailable"
  )
  expect_equal(res$overall$model, "RF")

  # Asking ONLY for the unavailable backend leaves nothing viable.
  expect_error(
    suppressMessages(compare_models_cv(pts, "z", "w", models = absent[1L],
                                       k = 3)),
    "no viable models"
  )
})


test_that("compare_models_cv restores the caller's RNG stream", {
  skip_if_not_installed("ranger")
  pts <- surf_test_points(60, seed = 9)

  set.seed(3); before <- runif(2)
  set.seed(3)
  invisible(compare_models_cv(pts, "z", "w", models = "RF", k = 3,
                              rf_args = list(num_trees = 30), quiet = TRUE))
  after <- runif(2)
  expect_equal(before, after)
})


# ---------------------------------------------------------------------------
# .merge_args(): the four arguments that define WHAT is being compared
# ---------------------------------------------------------------------------

test_that(".merge_args refuses per-model overrides of the comparison itself", {
  # compare_models_cv() exists to score several backends on ONE dataset with
  # ONE set of folds.  A `gwr_args = list(folds = ...)` or
  # `rf_args = list(data_sf = ...)` would quietly make the models
  # incomparable while the result still looked like a comparison table, so the
  # four defining arguments are stripped from the per-model lists.
  merge_args <- spatialkit:::.merge_args
  base <- list(data_sf = "THE DATA", response_var = "z",
               predictor_vars = "w", folds = "THE FOLDS",
               k = 5L, num_trees = 500L)

  # A NULL is the sharpest case: modifyList() DELETES an element it is given as
  # NULL, so an unstripped `folds = NULL` does not merely fail to override the
  # folds -- it removes them from the argument list entirely, and the CV
  # wrapper then silently falls back to random folds.
  # `.` stands in for the quote glyph throughout: sQuote() emits curly
  # quotes under a UTF-8 locale and straight ones otherwise.
  expect_warning(out <- merge_args(base, list(folds = NULL), "rf_args"),
                 "ignoring .folds.")
  expect_identical(out, base)
  expect_identical(out$folds, "THE FOLDS")

  # A non-NULL override of each protected name is refused the same way, while
  # everything else in the same list still takes effect.
  for (nm in c("data_sf", "response_var", "predictor_vars", "folds")) {
    extra <- stats::setNames(list("HIJACKED", 40L), c(nm, "num_trees"))
    expect_warning(o <- merge_args(base, extra, "gwr_args"),
                   sprintf("ignoring .%s.", nm))
    expect_identical(o[[nm]], base[[nm]], info = nm)
    expect_identical(o$num_trees, 40L, info = nm)
  }

  # Several at once are named together in one warning.
  expect_warning(
    both <- merge_args(base, list(data_sf = 1, folds = 2), "rf_args"),
    ".data_sf., .folds.")
  expect_identical(both, base)

  # The message says WHY, not just that something was dropped.
  expect_warning(merge_args(base, list(folds = 1), "rf_args"),
                 "same data and the same folds")
  expect_warning(merge_args(base, list(folds = 1), "bayes_args"),
                 "`bayes_args`")

  # Unprotected arguments merge by name, replacing rather than duplicating --
  # and silently, because there is nothing to complain about.
  expect_silent(ok <- merge_args(base, list(num_trees = 40L, mtry = 2L),
                                 "rf_args"))
  expect_identical(ok$num_trees, 40L)
  expect_identical(ok$mtry, 2L)
  expect_identical(ok$data_sf, "THE DATA")
  expect_equal(anyDuplicated(names(ok)), 0L)

  # Nothing supplied, nothing changed; an unnamed list is an error.
  expect_identical(merge_args(base, NULL, "rf_args"), base)
  expect_identical(merge_args(base, list(), "rf_args"), base)
  expect_error(merge_args(base, list(40L), "rf_args"),
               "every element of `rf_args` must be named")
})


# ---------------------------------------------------------------------------
# "all folds failed" must name the cause
# ---------------------------------------------------------------------------

test_that("an all-folds-failed CV warning names the first underlying error", {
  # "all 5 folds failed" on its own is not a diagnosis.  Overwhelmingly the
  # reason is a missing optional backend, and the fitter says so plainly when
  # called directly -- so the CV wrapper carries the first fold's error text
  # into both the warning and the log line.
  skip_if(requireNamespace("brms", quietly = TRUE),
          "brms is installed, so every fold would succeed")
  pts <- surf_test_points(50, seed = 4)

  lines <- capture_spatialkit_log(
    expect_warning(
      res <- suppressMessages(cv_bayes(pts, "z", "w", k = 2, seed = 1)),
      "all folds failed"
    )
  )
  # The condition the user sees carries the cause, not just the count.
  expect_warning(suppressMessages(cv_bayes(pts, "z", "w", k = 2, seed = 1)),
                 "First error: .*package 'brms' is required")
  # So does the log line, which also reports how many folds were attempted.
  expect_true(log_has(lines, "all 2 folds failed"))
  expect_true(log_has(lines, "First error: .*package 'brms' is required"))
  # And each fold said why as it failed.
  expect_true(log_has(lines, "fold 1 fit failed; skipping. Cause: .*'brms'"))

  # The result is still a well-formed, empty CV object.
  expect_equal(res$n_folds_attempted, 2L)
  expect_equal(res$n_folds_succeeded, 0L)
})


# ---------------------------------------------------------------------------
# compare_models_cv(): the shared folds can be sized from the range, and the
# overall table carries the Bayesian backend's calibration.
# ---------------------------------------------------------------------------

.fold_test_sets <- function(folds) lapply(folds, function(f) sort(f$test))

test_that("compare_models_cv forwards block_size to the shared folds", {
  skip_if_not_installed("ranger")
  pts <- surf_test_points(90, seed = 4)

  base <- compare_models_cv(pts, "z", "w", models = "RF", k = 3,
                            rf_args = list(num_trees = 40), quiet = TRUE)
  # Default: the same geometric folds as before the arguments existed, so the
  # comparison a script already made does not move.
  ref <- make_folds(pts, k = 3, method = "block_kfold", seed = 123)
  expect_identical(.fold_test_sets(base$rf_cv$folds), .fold_test_sets(ref$folds))

  # A block_size reaches make_folds(): 400-unit blocks on a 1000-unit extent
  # give a 2 x 2 grid and a different split from the default 3 x 3 one.
  big <- suppressWarnings(
    compare_models_cv(pts, "z", "w", models = "RF", k = 3, block_size = 400,
                      rf_args = list(num_trees = 40), quiet = TRUE))
  expect_false(identical(.fold_test_sets(big$rf_cv$folds), .fold_test_sets(ref$folds)))
})

test_that("compare_models_cv forwards auto_range to the shared folds", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("gstat")
  # On a response with a real range (exponential parameter 60, effective
  # range ~220 here -- short enough that the blocks it sizes still make a
  # split) make_folds() announces the range it used as the minimum block
  # size, and the shared folds differ from the geometric ones.  The message
  # alone would not prove it: compare_models_cv() catches a fold-construction
  # error and falls back to per-backend folds, so the fold sets are compared
  # as well.
  set.seed(3)
  n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  sp <- as.numeric(t(chol(exp(-d / 60) + diag(1e-4, n))) %*% rnorm(n))
  spat <- sf::st_as_sf(data.frame(x = x, y = y, w = rnorm(n)),
                       coords = c("x", "y"), crs = 32632)
  spat$z <- spat$w + sp
  expect_message(
    auto <- suppressWarnings(
      compare_models_cv(spat, "z", "w", models = "RF", k = 3, auto_range = TRUE,
                        rf_args = list(num_trees = 40), quiet = TRUE)),
    "using as minimum block size")
  geo <- make_folds(spat, k = 3, method = "block_kfold", seed = 123)
  ranged <- suppressWarnings(suppressMessages(
    make_folds(spat, k = 3, method = "block_kfold", seed = 123,
               auto_range = TRUE, response_var = "z", predictor_vars = "w")))
  expect_false(identical(.fold_test_sets(auto$rf_cv$folds), .fold_test_sets(geo$folds)))
  expect_identical(.fold_test_sets(auto$rf_cv$folds), .fold_test_sets(ranged$folds))
})

test_that(".overall_with_coverage appends a predictive_coverage summary and nothing else", {
  f  <- spatialkit:::.overall_with_coverage
  ov <- data.frame(RMSE = 1, MAE = 0.8, n_pred = 10L)
  expect_identical(f(ov, NULL), ov)
  expect_identical(f(ov, list()), ov)
  expect_identical(f(ov, list(1, 2)), ov)          # unnamed: nothing to name a column by
  pc <- list(coverage_50 = 0.48, coverage_80 = 0.79, coverage_95 = 0.96,
             mean_CRPS = 0.31)
  out <- f(ov, pc)
  expect_named(out, c("RMSE", "MAE", "n_pred", "coverage_50", "coverage_80",
                      "coverage_95", "mean_CRPS"))
  expect_identical(out$coverage_95, 0.96)
  expect_identical(out$mean_CRPS, 0.31)
  # A malformed entry becomes NA rather than an error or a list column.
  odd <- f(ov, list(coverage_95 = c(0.9, 0.95), mean_CRPS = "x"))
  expect_true(is.na(odd$coverage_95))
  expect_true(is.na(odd$mean_CRPS))
  expect_type(odd$mean_CRPS, "double")
})

test_that("compare_models_cv keeps `model` as the last overall column", {
  skip_if_not_installed("ranger")
  pts <- surf_test_points(60, seed = 4)
  res <- compare_models_cv(pts, "z", "w", models = "RF", k = 3,
                           rf_args = list(num_trees = 40), quiet = TRUE)
  expect_identical(names(res$overall)[ncol(res$overall)], "model")
  expect_false(any(c("coverage_95", "mean_CRPS") %in% names(res$overall)))
})

test_that("compare_models_cv carries coverage and CRPS into overall when a Bayesian model ran", {
  # Full Stan runs are too slow for the suite; a stand-in cv_bayes() that
  # returns the documented shape proves the plumbing from predictive_coverage
  # to the overall table, and that the point-prediction rows get NA.
  skip_if_not_installed("ranger")
  pts <- surf_test_points(60, seed = 4)
  fake_bayes <- function(data_sf, response_var, predictor_vars, folds = NULL, ...) {
    fm <- data.frame(fold = 1:2, n_train = c(30L, 30L), n_test = c(30L, 30L),
                     n_pred = c(30L, 30L), RMSE = c(1, 1.2), MAE = c(0.8, 0.9),
                     MAPE = NA_real_, SMAPE = NA_real_, R2 = c(0.5, 0.4),
                     Adj_R2 = NA_real_, n_MAPE = 0L, n_SMAPE = 0L,
                     CRPS = c(0.3, 0.35), coverage_50 = c(0.5, 0.46),
                     coverage_80 = c(0.8, 0.78), coverage_95 = c(0.97, 0.93))
    list(overall = data.frame(RMSE = 1.1, MAE = 0.85, MAPE = NA_real_,
                              SMAPE = NA_real_, R2 = 0.45, Adj_R2 = NA_real_,
                              n_pred = 60L, n_MAPE = 0L, n_SMAPE = 0L),
         fold_metrics = fm, predictions = NULL, folds = folds,
         n_folds_attempted = 2L, n_folds_succeeded = 2L,
         formula = "z ~ w",
         predictive_coverage = list(coverage_50 = 0.48, coverage_80 = 0.79,
                                    coverage_95 = 0.95, mean_CRPS = 0.325))
  }
  local_mocked_bindings(cv_bayes = fake_bayes,
                        .model_available = function(model_name) TRUE,
                        .package = "spatialkit")
  res <- compare_models_cv(pts, "z", "w", models = c("RF", "Bayesian"), k = 3,
                           rf_args = list(num_trees = 40), quiet = TRUE)
  ov <- res$overall
  # Backends run in a fixed order (GWR, Bayesian, RF) whatever `models` says.
  expect_identical(ov$model, c("Bayesian", "RF"))
  expect_true(all(c("coverage_50", "coverage_80", "coverage_95", "mean_CRPS") %in% names(ov)))
  expect_identical(names(ov)[ncol(ov)], "model")
  expect_true(is.na(ov$coverage_95[ov$model == "RF"]))
  expect_identical(ov$coverage_95[ov$model == "Bayesian"], 0.95)
  expect_identical(ov$mean_CRPS[ov$model == "Bayesian"], 0.325)
  # The point metrics on the RF row are untouched by the extra columns.
  expect_true(is.finite(ov$RMSE[ov$model == "RF"]))
})
