# tests/testthat/test-fix-integration-repeats.R
# ---------------------------------------------------------------------------
# One warning for one fact, when a function runs the same step many times.
#
# 1. predict() on a GWR counts, per call, the locations it cannot estimate
#    and those beyond the reach of the fit.  cv_spatial() with a GWR learner
#    raises each once with totals, as cv_gwr() does, and
#    cv_block_size_sweep() and select_features_forward(), which run many
#    cross-validations, raise each once for the sweep.
# 2. A parallel run hands the package's own classed warnings back with their
#    class and fields.
# 3. prep_model_data()'s warnings about a layer are raised once by a function
#    that prepares it several times (compare_models_cv(),
#    cv_block_size_sweep()).
# 4. cv_bayes(), and cv_spatial() with a Bayesian learner, raise the fit's
#    classed warnings once per class, naming the folds.
# 5. resolution_profile() does not call a pure-nugget variogram unusable.
# ---------------------------------------------------------------------------

fir_conditions <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = val, conds = w,
       messages = vapply(w, conditionMessage, character(1)),
       classes = vapply(w, function(c) class(c)[1L], character(1)))
}
fir_one <- function(r, class) r$conds[[which(r$classes == class)]]

fir_gwr_layer <- function(n = 160, seed = 9) {
  set.seed(seed)
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n, 0, 0.3)
  d
}

# The lm stand-in with the classes of a learner whose predict() raises the
# GWR's count of locations it cannot estimate: one location per call, the
# same text in every fold of equal size.
fir_counting_fit <- function(train_sf) {
  f <- lm_spatial_fit(train_sf, "z", "w")
  class(f) <- c("fir_counting_fit", class(f))
  f
}
predict.fir_counting_fit <- function(object, newdata = NULL, ...) {
  out <- NextMethod()
  n <- length(out)
  .gwr_warn("spatialkit_gwr_not_estimable", "%s",
            .gwr_not_estimable_msg("predict.gwr_fit", 1L, n, 1L, "fixed bandwidth (5)"),
            data = list(n_na = 1L, n = n, n_empty = 1L,
                        bw_text = "fixed bandwidth (5)", calls = 1L),
            tallied = TRUE)
  out
}
registerS3method("predict", "fir_counting_fit", predict.fir_counting_fit)


# --- 1. the GWR's counted warnings through the generic loops -------------------

test_that("cv_spatial() raises a GWR learner's counted warnings once, with cv_gwr()'s totals", {
  skip_if_not_installed("sp")
  d <- fir_gwr_layer()
  # Four strips as folds and a fixed bandwidth that does not reach across a
  # strip: several folds hold locations no training point carries weight for.
  strips <- cut(sf::st_coordinates(d)[, 1], 4, labels = FALSE)
  fixed <- function(tr) fit_gwr_model(tr, "z", "a", adaptive = FALSE, bandwidth = 120)
  lines <- capture_spatialkit_log(
    cs <- fir_conditions(cv_spatial(d, "z", "a", fit_fn = fixed, folds = strips)),
    level = logger::WARN)
  cg <- fir_conditions(cv_gwr(d, "z", "a", folds = strips, adaptive = FALSE, bandwidth = 120))
  expect_identical(sum(cs$classes == "spatialkit_gwr_not_estimable"), 1L)
  expect_identical(sum(cs$classes == "spatialkit_gwr_extrapolation"), 0L)
  ws <- fir_one(cs, "spatialkit_gwr_not_estimable")
  wg <- fir_one(cg, "spatialkit_gwr_not_estimable")
  expect_match(conditionMessage(ws), "^cv_spatial\\(\\): ")
  expect_gte(ws$calls, 2L)                         # more than one fold had any
  expect_identical(ws$n_na, wg$n_na)
  expect_identical(ws$n, wg$n)
  expect_identical(ws$calls, wg$calls)
  expect_identical(sum(is.na(cs$value$predictions$yhat)), ws$n_na)
  expect_equal(cs$value$overall$RMSE, cg$value$overall$RMSE)
  # Logged once, as the total, and not once per fold.
  expect_identical(sum(grepl("no estimable local", lines, fixed = TRUE)), 1L)
  expect_true(log_has(lines, "cv_spatial\\(\\): [0-9]+ of [0-9]+ location"))

  # Adaptive folds reach as far as their kernels: a strip moved away is
  # beyond every training point's reach.
  far <- d
  sf::st_geometry(far)[strips == 4] <- sf::st_geometry(d)[strips == 4] + c(6000, 0)
  sf::st_crs(far) <- 32632
  strips2 <- cut(sf::st_coordinates(far)[, 1], c(-Inf, 5e5 + 250, 5e5 + 500, 5e5 + 800, Inf),
                 labels = FALSE)
  adaptive <- function(tr) fit_gwr_model(tr, "z", "a", bandwidth = 20)
  cs2 <- fir_conditions(cv_spatial(far, "z", "a", fit_fn = adaptive, folds = strips2))
  cg2 <- fir_conditions(cv_gwr(far, "z", "a", folds = strips2, bandwidth = 20))
  expect_identical(sum(cs2$classes == "spatialkit_gwr_extrapolation"), 1L)
  w2 <- fir_one(cs2, "spatialkit_gwr_extrapolation")
  expect_match(conditionMessage(w2), "^cv_spatial\\(\\): ")
  expect_identical(w2$n_far, fir_one(cg2, "spatialkit_gwr_extrapolation")$n_far)
  expect_gt(w2$n_far, 0L)

  # A learner that raises neither class is not touched.
  expect_length(fir_conditions(cv_spatial(d, "z", "a", folds = strips,
                                          fit_fn = function(tr) lm_spatial_fit(tr, "z", "a")))$conds,
                0L)
})

test_that("cv_block_size_sweep() and select_features_forward() raise them once for the sweep", {
  skip_if_not_installed("sp")
  d <- fir_gwr_layer(n = 200, seed = 21)
  fixed <- function(tr) fit_gwr_model(tr, "z", "a", adaptive = FALSE, bandwidth = 90)
  # Each cross-validation of the sweep on its own: one warning per size.
  lines <- capture_spatialkit_log(
    sw <- fir_conditions(cv_block_size_sweep(d, "z", "a", fit_fn = fixed, k = 3,
                                             block_sizes = c(250, 330), sac = 100,
                                             quiet = TRUE)),
    level = logger::WARN)
  na_by_cv <- vapply(attr(sw$value, "results"),
                     function(r) sum(is.na(r$predictions$yhat)), integer(1))
  skip_if(sum(na_by_cv > 0L) < 2L, "fewer than two of the cross-validations had an NA here")
  expect_identical(sum(sw$classes == "spatialkit_gwr_not_estimable"), 1L)
  w <- fir_one(sw, "spatialkit_gwr_not_estimable")
  expect_match(conditionMessage(w), "^cv_block_size_sweep\\(\\): ")
  expect_identical(w$n_na, sum(na_by_cv))
  expect_match(conditionMessage(w),
               "a location predicted by several fits is counted once for each", fixed = TRUE)
  expect_identical(sum(grepl("no estimable local", lines, fixed = TRUE)), 1L)

  # The selection sweep: every candidate set is cross-validated on the same
  # folds, and each cross-validation had its own warning.
  sel_fn <- function(tr, v) fit_gwr_model(tr, "z", v, adaptive = FALSE, bandwidth = 90)
  lines2 <- capture_spatialkit_log(
    fs <- fir_conditions(select_features_forward(d, "z", c("a", "b"), sel_fn, k = 3,
                                                 seed = 7, quiet = TRUE)),
    level = logger::WARN)
  expect_identical(sum(fs$classes == "spatialkit_gwr_not_estimable"), 1L)
  wf <- fir_one(fs, "spatialkit_gwr_not_estimable")
  expect_match(conditionMessage(wf), "^select_features_forward\\(\\): ")
  expect_match(conditionMessage(wf),
               "a location predicted by several fits is counted once for each", fixed = TRUE)
  expect_gt(wf$calls, 1L)
  expect_identical(sum(grepl("no estimable local", lines2, fixed = TRUE)), 1L)
  expect_s3_class(fs$value, "feature_selection")
  # The collector is closed again: a later predict() warns under its own name.
  expect_identical(.gwr_predict_state$tally %||% 0L, 0L)
})

test_that("a total over fits with different bandwidths names no single bandwidth", {
  raise <- function(bw) .gwr_warn(
    "spatialkit_gwr_not_estimable", "%s",
    .gwr_not_estimable_msg("predict.gwr_fit", 2L, 10L, 2L, bw),
    data = list(n_na = 2L, n = 10L, n_empty = 2L, bw_text = bw, calls = 1L),
    tallied = TRUE)
  same <- fir_conditions(.gwr_tally_predict_warnings(
    { raise("fixed bandwidth (50)"); raise("fixed bandwidth (50)") }, caller = "f"))
  expect_length(same$conds, 1L)
  expect_identical(same$conds[[1]]$n_na, 4L)
  expect_match(same$messages, "within the fixed bandwidth (50).", fixed = TRUE)
  expect_false(grepl("counted once for each", same$messages, fixed = TRUE))
  differ <- fir_conditions(.gwr_tally_predict_warnings(
    { raise("fixed bandwidth (50)"); raise("fixed bandwidth (80)") }, caller = "f"))
  expect_match(differ$messages,
               "within the bandwidth of the fit that predicted it (the fits' bandwidths differ).",
               fixed = TRUE)
  # Nested: the inner total goes to the outer collector alone.
  nested <- fir_conditions(.gwr_tally_predict_warnings({
    .gwr_tally_predict_warnings(raise("fixed bandwidth (50)"), caller = "inner")
    .gwr_tally_predict_warnings({ raise("fixed bandwidth (50)"); raise("fixed bandwidth (50)") },
                                caller = "inner")
  }, caller = "outer", refits = TRUE))
  expect_length(nested$conds, 1L)
  expect_match(nested$messages, "^outer\\(\\): 6 of 30 location")
  expect_identical(nested$conds[[1]]$calls, 3L)
  expect_match(nested$messages,
               "(totals over the 3 prediction calls that found any; a location predicted by several fits is counted once for each)",
               fixed = TRUE)
})


# --- 2. classed warnings across a fork -----------------------------------------

test_that("a parallel run hands the package's classed warnings back as they were raised", {
  d <- surf_test_points(80, seed = 5)
  f <- make_folds(d, k = 4, method = "random_kfold", seed = 1)
  # Sequential: four folds, one location each, the same text in every fold.
  sq <- fir_conditions(cv_spatial(d, "z", "w", fit_fn = fir_counting_fit, folds = f, seed = 1))
  expect_identical(sq$classes, "spatialkit_gwr_not_estimable")
  expect_identical(sq$conds[[1]]$n_na, 4L)
  expect_identical(sq$conds[[1]]$calls, 4L)

  skip_on_cran()
  skip_on_os("windows")
  skip_if(.resolve_n_cores(2L) < 2L, "single-core machine")
  noting_fit <- function(tr) {
    .warn_classed("spatialkit_fir_note", "the same note from every fold", data = list(k = 7L))
    fir_counting_fit(tr)
  }
  lines <- capture_spatialkit_log(
    pl <- fir_conditions(cv_spatial(d, "z", "w", fit_fn = noting_fit, folds = f, seed = 1,
                                    parallel = 2L)),
    level = logger::WARN)
  # The counted warning: every fold's count arrives, though the four texts
  # are identical, and the total is the sequential one.
  expect_identical(sum(pl$classes == "spatialkit_gwr_not_estimable"), 1L)
  wp <- fir_one(pl, "spatialkit_gwr_not_estimable")
  expect_identical(wp$n_na, 4L)
  expect_identical(wp$calls, 4L)
  expect_identical(conditionMessage(wp), conditionMessage(sq$conds[[1]]))
  expect_identical(sum(grepl("no estimable local", lines, fixed = TRUE)), 1L)
  # Another classed warning: once per distinct text, class and fields kept,
  # and marked with the fold that raised it first.
  expect_identical(sum(pl$classes == "spatialkit_fir_note"), 1L)
  wn <- fir_one(pl, "spatialkit_fir_note")
  expect_identical(wn$k, 7L)
  expect_true(wn$fold %in% 1:4)
  expect_equal(pl$value$overall$RMSE, sq$value$overall$RMSE)
})


# --- 3. the layer's warnings, once per run -------------------------------------

fir_holed_layer <- function(n = 150, seed = 31) {
  set.seed(seed)
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n), c0 = 3),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n, 0, 0.3)
  d$z[c(5, 60)] <- NA
  d$a[c(17, 100, 101)] <- NA
  d
}

test_that("cv_block_size_sweep() raises the layer's warnings once, not once per size", {
  d <- fir_holed_layer()
  lm_fn <- function(tr) lm_spatial_fit(tr, "z", "a")
  # One cross-validation raises each once.
  one <- fir_conditions(cv_spatial(d, "z", c("a", "c0"), fit_fn = lm_fn, k = 3, seed = 1))
  expect_identical(sum(one$classes == "spatialkit_rows_dropped"), 1L)
  expect_identical(sum(one$classes == "spatialkit_constant_predictor"), 1L)
  # The sweep prepares the layer for its ladder and once per size, and the
  # random-fold reference: four preparations here.
  lines <- capture_spatialkit_log(
    sw <- fir_conditions(cv_block_size_sweep(d, "z", c("a", "c0"), fit_fn = lm_fn, k = 3,
                                             block_sizes = c(200, 300), sac = 100,
                                             quiet = TRUE)),
    level = logger::WARN)
  expect_identical(nrow(sw$value), 3L)
  expect_identical(sum(sw$classes == "spatialkit_rows_dropped"), 1L)
  expect_identical(sum(sw$classes == "spatialkit_constant_predictor"), 1L)
  expect_identical(fir_one(sw, "spatialkit_rows_dropped")$message,
                   fir_one(one, "spatialkit_rows_dropped")$message)
  expect_match(fir_one(sw, "spatialkit_rows_dropped")$message, "dropped 5 of 150 row(s)",
               fixed = TRUE)
  expect_identical(sum(grepl("dropped 5 of 150", lines, fixed = TRUE)), 1L)
  expect_identical(sum(grepl("take a single value", lines, fixed = TRUE)), 1L)
  # Every cross-validation still records the rows it dropped.
  expect_true(all(vapply(attr(sw$value, "results"), function(r) r$n_dropped, integer(1)) == 5L))

  # The run is closed when the function returns, and when it stops: the next
  # call warns again.
  expect_null(.once_state$run)
  expect_error(suppressWarnings(
    cv_block_size_sweep(d, "z", c("a", "c0"), fit_fn = lm_fn, k = 3,
                        block_sizes = c(200, 300), max_fits = 1, quiet = TRUE)),
    "above the budget")
  expect_null(.once_state$run)
  again <- fir_conditions(cv_spatial(d, "z", c("a", "c0"), fit_fn = lm_fn, k = 3, seed = 1))
  expect_identical(sum(again$classes == "spatialkit_rows_dropped"), 1L)
})

test_that("compare_models_cv() raises the layer's warnings once, not once per model", {
  skip_if_not_installed("sp")
  skip_if_not_installed("ranger")
  d <- fir_holed_layer()
  lines <- capture_spatialkit_log(
    cmp <- fir_conditions(compare_models_cv(d, "z", "a", models = c("RF", "GWR"), k = 3,
                                            quiet = TRUE,
                                            rf_args = list(num_trees = 30, num_threads = 1),
                                            gwr_args = list(bandwidth = 30))),
    level = logger::WARN)
  expect_identical(sum(cmp$classes == "spatialkit_rows_dropped"), 1L)
  expect_identical(sum(grepl("dropped 5 of 150", lines, fixed = TRUE)), 1L)
  expect_identical(sort(cmp$value$overall$model), c("GWR", "RF"))
  expect_identical(cmp$value$gwr_cv$n_dropped, 5L)
  expect_identical(cmp$value$rf_cv$n_dropped, 5L)
  expect_null(.once_state$run)
})

test_that("select_features_forward() reports a constant candidate once, not once per set", {
  d <- fir_holed_layer()
  d$b <- rnorm(nrow(d))
  # `c0` is in one set at every step: {c0}, then {a, c0}, then {a, b, c0}.
  fit_fn <- function(tr, v) lm_spatial_fit(tr, "z", setdiff(v, "c0"))
  lines <- capture_spatialkit_log(
    fs <- fir_conditions(select_features_forward(d, "z", c("a", "b", "c0"), fit_fn, k = 3,
                                                 seed = 3, quiet = TRUE)),
    level = logger::WARN)
  expect_gte(sum(grepl("c0", fs$value$history$variable, fixed = TRUE)), 2L)
  expect_identical(sum(fs$classes == "spatialkit_constant_predictor"), 1L)
  expect_identical(sum(grepl("take a single value", lines, fixed = TRUE)), 1L)
  expect_null(.once_state$run)
})

test_that("a run raises a text once, another text as usual, and nests", {
  said <- function() .warn_classed("spatialkit_fir_once", "said %d", 1L)
  expect_length(fir_conditions({ said(); said() })$conds, 2L)       # no run open
  r <- fir_conditions({
    end <- .once_per_run("spatialkit_fir_once")
    said(); said()
    .warn_classed("spatialkit_fir_once", "said %d", 2L)              # other text
    .warn_classed("spatialkit_fir_other", "said %d", 1L)             # other class,
    .warn_classed("spatialkit_fir_other", "said %d", 1L)             # not covered
    inner_end <- .once_per_run("spatialkit_fir_other")               # a run inside
    said()
    .warn_classed("spatialkit_fir_other", "said %d", 3L)
    .warn_classed("spatialkit_fir_other", "said %d", 3L)
    inner_end()
    said()                                                           # still said
    .warn_classed("spatialkit_fir_other", "said %d", 3L)             # outer run again
    end()
  })
  expect_identical(r$messages, c("said 1", "said 2", "said 1", "said 1", "said 3", "said 3"))
  expect_null(.once_state$run)
  expect_length(fir_conditions(said())$conds, 1L)
})


# --- 4. a Bayesian fit's warnings, once per class and by fold -------------------

fir_bayes_pts <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 5 + 2 * d$w + rnorm(n, 0, 0.5)
  d
}

# The lm stand-in for fit_bayesian_spatial_model(), raising what the real
# fitter raises: fold 1 and 3 a range below the basis, fold 2 a sampler that
# did not converge, fold 3 a basis too coarse and too large for its rows.
# Which fold is running is read as the collector reads it.
fir_fake_bayes <- function(data_sf, response_var, predictor_vars, ..., seed = 123) {
  i <- as.integer(.cv_fold_state$fold)
  f <- lm_spatial_fit(data_sf, "z", "w")
  f$info$convergence_ok <- i != 2L
  f$info$gp_adequacy <- list(status = if (i == 3L) "below_resolution" else "ok")
  if (i == 3L)
    .warn_classed("spatialkit_gp_basis_size",
                  "fit_bayesian_spatial_model(): 40 row(s) against 900 basis functions.")
  if (i %in% c(1L, 3L))
    .warn_classed("spatialkit_gp_resolution",
                  "fit_bayesian_spatial_model(): a variogram fitted before sampling puts the range at %d metre.",
                  100L + i)
  if (i == 2L)
    .warn_classed("spatialkit_convergence",
                  "fit_bayesian_spatial_model(): the sampler's diagnostics fail: 3 divergent transition(s).",
                  data = list(failed = c("divergent", "ess")))
  if (i == 3L)
    .warn_classed("spatialkit_gp_adequacy",
                  "fit_bayesian_spatial_model(): 23%% of the length-scale draws fall below the floor.")
  f
}

fir_expect_fold_warnings <- function(r, caller) {
  classes <- c("spatialkit_convergence", "spatialkit_gp_basis_size",
               "spatialkit_gp_resolution", "spatialkit_gp_adequacy")
  # One per class, and nothing else.
  expect_identical(sort(r$classes), sort(classes))
  expect_true(all(startsWith(r$messages, paste0(caller, "(): "))))
  conv <- fir_one(r, "spatialkit_convergence")
  expect_match(conditionMessage(conv), "did not converge in 1 of 3 fold(s) (fold 2)", fixed = TRUE)
  expect_match(conditionMessage(conv),
               "The checks that failed: fold 2: divergent transitions, effective sample size.",
               fixed = TRUE)
  expect_match(conditionMessage(conv), "adapt_delta", fixed = TRUE)
  expect_false(grepl("E-BFMI", conditionMessage(conv), fixed = TRUE))
  expect_identical(conv$folds, 2L)
  expect_identical(conv$failed, list(`2` = c("divergent", "ess")))
  res <- fir_one(r, "spatialkit_gp_resolution")
  expect_match(conditionMessage(res), "in 2 of 3 fold(s) (fold 1, 3)", fixed = TRUE)
  expect_match(conditionMessage(res),
               "Fold 1: a variogram fitted before sampling puts the range at 101 metre.",
               fixed = TRUE)
  expect_identical(res$folds, c(1L, 3L))
  expect_identical(names(res$fold_messages), c("1", "3"))
  expect_match(res$fold_messages[["3"]], "at 103 metre", fixed = TRUE)
  adq <- fir_one(r, "spatialkit_gp_adequacy")
  expect_match(conditionMessage(adq), "in 1 of 3 fold(s) (fold 3)", fixed = TRUE)
  expect_match(conditionMessage(adq), "Fold 3: 23% of the length-scale draws", fixed = TRUE)
  big <- fir_one(r, "spatialkit_gp_basis_size")
  expect_match(conditionMessage(big), "in 1 of 3 fold(s) (fold 3)", fixed = TRUE)
  expect_identical(big$folds, 3L)
}

test_that("cv_bayes() raises a fit's classed warnings once per class, naming the folds", {
  d <- fir_bayes_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  local_mocked_bindings(fit_bayesian_spatial_model = fir_fake_bayes, .package = "spatialkit")
  lines <- capture_spatialkit_log(
    r <- fir_conditions(cv_bayes(d, "z", "w", folds = f, seed = 1)),
    level = logger::WARN)
  fir_expect_fold_warnings(r, "cv_bayes")
  expect_match(conditionMessage(fir_one(r, "spatialkit_convergence")),
               "fold_metrics$convergence_ok marks them", fixed = TRUE)
  expect_match(conditionMessage(fir_one(r, "spatialkit_gp_adequacy")),
               "fold_metrics$gp_adequacy has each fold's status", fixed = TRUE)
  expect_identical(r$value$fold_metrics$convergence_ok, c(TRUE, FALSE, TRUE))
  # Each fold's own line is still in the log, beside the four totals.
  expect_identical(sum(grepl("^.*fit_bayesian_spatial_model\\(\\): a variogram fitted", lines)), 2L)
  expect_identical(sum(grepl("fit_bayesian_spatial_model(): the sampler's diagnostics fail",
                             lines, fixed = TRUE)), 1L)
  expect_identical(sum(grepl("cv_bayes(): ", lines, fixed = TRUE)), 4L)
  # It can still be caught by class.
  expect_identical(
    tryCatch(suppressMessages(cv_bayes(d, "z", "w", folds = f, seed = 1)),
             spatialkit_convergence = function(c) "caught"),
    "caught")
})

test_that("cv_spatial() does the same for a learner that fits the Bayesian model", {
  d <- fir_bayes_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  local_mocked_bindings(fit_bayesian_spatial_model = fir_fake_bayes, .package = "spatialkit")
  learner <- function(tr) fit_bayesian_spatial_model(tr, "z", "w")
  r <- fir_conditions(cv_spatial(d, "z", "w", fit_fn = learner, folds = f, seed = 1))
  fir_expect_fold_warnings(r, "cv_spatial")
  # No fold_metrics column to point at here.
  expect_false(any(grepl("fold_metrics$", r$messages, fixed = TRUE)))
  expect_match(conditionMessage(fir_one(r, "spatialkit_convergence")),
               "in the learner's call", fixed = TRUE)
  # A fold that did not converge and then failed to predict was not scored,
  # so the convergence warning, which is about scores, leaves it out.
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, ...) {
      ft <- fir_fake_bayes(data_sf)
      if (identical(as.integer(.cv_fold_state$fold), 2L)) class(ft) <- c("fir_no_predict", "spatial_fit")
      ft
    }, .package = "spatialkit")
  r2 <- fir_conditions(cv_spatial(d, "z", "w", fit_fn = learner, folds = f, seed = 1))
  expect_false("spatialkit_convergence" %in% r2$classes)
  expect_identical(r2$value$fold_status$status, c("ok", "error", "ok"))
  expect_true(any(grepl("1 of 3 fold(s) failed", r2$messages, fixed = TRUE)))

  # Forked folds: the same four, with the folds named.
  skip_on_cran()
  skip_on_os("windows")
  skip_if(.resolve_n_cores(2L) < 2L, "single-core machine")
  local_mocked_bindings(fit_bayesian_spatial_model = fir_fake_bayes, .package = "spatialkit")
  rp <- fir_conditions(cv_spatial(d, "z", "w", fit_fn = learner, folds = f, seed = 1,
                                  parallel = 2L))
  fir_expect_fold_warnings(rp, "cv_spatial")
})

test_that("a convergence summary without the folds' own warnings keeps the general advice", {
  d <- fir_bayes_pts()
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, ...) {
      ft <- lm_spatial_fit(data_sf, "z", "w")
      ft$info$convergence_ok <- !identical(as.integer(.cv_fold_state$fold), 3L)
      ft
    }, .package = "spatialkit")
  r <- fir_conditions(cv_bayes(d, "z", "w", folds = f, seed = 1))
  expect_identical(r$classes, "spatialkit_convergence")
  expect_match(r$messages, "did not converge in 1 of 3 fold(s) (fold 3)", fixed = TRUE)
  expect_match(r$messages, "The log says which checks failed", fixed = TRUE)
  expect_identical(r$conds[[1]]$folds, 3L)
})


# --- 5. a pure-nugget variogram is not "no usable model" -----------------------

test_that("resolution_profile() says what a no-structure variogram lacks", {
  skip_if_not_installed("gstat")
  set.seed(4401); n <- 240
  pts <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                                 noise = rnorm(n)), coords = c("x", "y"), crs = 32632)
  no <- suppressWarnings(estimate_sac_range(pts, "noise"))
  skip_if_not(.sac_no_structure(no), "the structure test did not refuse this draw")
  expect_true(.vgm_is_pure_nugget(attr(no, "variogram_model")))
  for (sac in list(no, NULL)) {
    lines <- capture_spatialkit_log(
      r <- fir_conditions(resolution_profile(pts, "noise", sac = sac, n_levels = 4, nstart = 2)))
    expect_identical(attr(r$value, "cp_noise")$source, "finest-level residual mean square")
    hit <- grepl("so `cp` takes its noise variance from the residual mean square", r$messages,
                 fixed = TRUE)
    expect_identical(sum(hit), 1L)
    expect_match(r$messages[hit],
                 "found no spatial structure and gives no nugget that can be told from the sill, so `cp`",
                 fixed = TRUE)
    # The model is the pure nugget the log line names, and it is used
    # elsewhere: nothing calls it unusable.
    expect_false(any(grepl("no usable variogram model", r$messages, fixed = TRUE)))
    expect_true(log_has(lines, "found no spatial structure .*its model is a pure nugget"))
  }
  expect_match(fir_conditions(resolution_profile(pts, "noise", sac = no, n_levels = 4,
                                                 nstart = 2))$messages[1],
               "^resolution_profile\\(\\): `sac` found no spatial structure")
})
