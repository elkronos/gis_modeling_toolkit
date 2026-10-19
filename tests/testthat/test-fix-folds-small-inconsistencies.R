# tests/testthat/test-fix-folds-small-inconsistencies.R
# ---------------------------------------------------------------------------
# Three labelling and silence inconsistencies; none changes a score.
#  1. fold_separation()'s n_unknown_ids counted fold entries, not rows.
#  2. The cv_*() functions and compare_models_cv() dropped block_size and
#     auto_range without a word when folds were supplied.
#  3. compare_models(newdata = ) labelled its rows "newdata" beside Moran's I
#     columns computed on the training residuals, with nothing saying so.
# Each test fails on the code before the fix.
# ---------------------------------------------------------------------------

si_pts <- function(n = 200, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1e4), y = 5e6 + runif(n, 0, 1e4),
                               a = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + rnorm(n)
  d
}
si_q <- function(expr) r2_quiet(suppressMessages(expr))
# The conditions of one class a call raises; other warnings are muffled.
si_classed <- function(expr, class) {
  hit <- list()
  withCallingHandlers(expr, warning = function(c) {
    if (inherits(c, class)) hit[[length(hit) + 1L]] <<- c
    invokeRestart("muffleWarning")
  })
  hit
}

# ---- 1. n_unknown_ids counts rows -------------------------------------------

test_that("fold_separation() counts each unknown row once, whatever the fold count", {
  u <- si_pts()
  u$..row_id <- seq_len(nrow(u))
  u195 <- u[-(1:5), ]
  for (k in c(3L, 5L, 10L)) {
    f <- si_q(make_folds(u, k = k, method = "random_kfold", seed = 1))
    # Folds built on 200 rows, measured on 195 of them.  Five rows are
    # missing; they used to read as 15, 25 and 50.
    s <- si_q(fold_separation(f, u195))
    expect_identical(attr(s, "n_unknown_ids"), 5L, info = paste("k =", k))
    expect_identical(sum(s$n_test), 195L)
    # ... which is what cv_*() report for the same folds on the same layer.
    cv <- si_q(suppressWarnings(cv_spatial(u195, "z", "a", fit_fn = r2_fit, folds = f)))
    expect_identical(cv$n_unknown_ids, 5L)
  }
  # 200 buffered leave-one-out folds read as 996.
  b <- si_q(suppressWarnings(make_folds(u, method = "buffered_loo", buffer = 300)))
  expect_identical(attr(si_q(fold_separation(b, u195)), "n_unknown_ids"), 5L)
  # A bare split list is counted the same way.
  expect_identical(attr(si_q(fold_separation(b$folds, u195)), "n_unknown_ids"), 5L)
  # Nothing missing, nothing counted.
  f <- si_q(make_folds(u, k = 5, method = "random_kfold", seed = 1))
  expect_identical(attr(fold_separation(f, u), "n_unknown_ids"), 0L)
  # The log line counts rows too.
  lines <- capture_spatialkit_log(fold_separation(f, u195), level = logger::INFO)
  expect_true(log_has(lines, "the folds name 5 row\\(s\\) `data_sf` does not have"))
})

# ---- 2. block_size and auto_range beside supplied folds ---------------------

test_that("cv_spatial() warns when block_size or auto_range arrive with folds", {
  u <- si_pts()
  f5 <- si_q(make_folds(u, k = 5, method = "random_kfold", seed = 1))
  base <- si_q(cv_spatial(u, "z", "a", fit_fn = r2_fit, folds = f5))
  for (args in list(list(auto_range = TRUE), list(block_size = 3000),
                    list(auto_range = TRUE, block_size = 3000))) {
    w <- list()
    cv <- withCallingHandlers(
      si_q(do.call(cv_spatial, c(list(u, "z", "a", fit_fn = r2_fit, folds = f5), args))),
      spatialkit_fold_args_ignored = function(c) {
        w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
      })
    expect_length(w, 1L)
    msg <- conditionMessage(w[[1L]])
    expect_match(msg, "^cv_spatial\\(\\): ")
    expect_match(msg, "no effect because `folds` is supplied")
    if (!is.null(args$block_size)) expect_match(msg, "`block_size`", fixed = TRUE)
    if (isTRUE(args$auto_range)) expect_match(msg, "`auto_range = TRUE`", fixed = TRUE)
    # The run itself is the one on the supplied folds, as before.
    expect_identical(cv$predictions, base$predictions)
  }
})

test_that("every cv_*() and compare_models_cv() raise the same classed warning", {
  u <- si_pts()
  f5 <- si_q(make_folds(u, k = 5, method = "random_kfold", seed = 1))
  # tryCatch() by class leaves the call at the warning, before any model is
  # fitted, so the backends' packages are not needed to see it raised.
  first <- function(expr) tryCatch({ expr; NULL }, spatialkit_fold_args_ignored = function(w) w)
  hits <- list(
    cv_gwr   = if (requireNamespace("GWmodel", quietly = TRUE) &&
                   requireNamespace("sp", quietly = TRUE))
      first(si_q(cv_gwr(u, "z", "a", folds = f5, bandwidth = 40, block_size = 3000))),
    cv_bayes = first(si_q(cv_bayes(u, "z", "a", folds = f5, auto_range = TRUE))),
    cv_rf    = if (requireNamespace("ranger", quietly = TRUE))
      first(si_q(cv_rf(u, "z", "a", folds = f5, num_trees = 20, block_size = 3000))),
    compare_models_cv = if (requireNamespace("ranger", quietly = TRUE))
      first(si_q(compare_models_cv(u, "z", "a", models = "RF", folds = f5,
                                   quiet = TRUE, auto_range = TRUE))))
  hits <- Filter(Negate(is.null), hits)
  # cv_bayes() needs no backend to get as far as the warning.
  expect_true("cv_bayes" %in% names(hits))
  for (nm in names(hits)) {
    expect_s3_class(hits[[nm]], "spatialkit_fold_args_ignored")
    expect_s3_class(hits[[nm]], "warning")
    # Named for the function the user called.
    expect_match(conditionMessage(hits[[nm]]), sprintf("^%s\\(\\): ", nm), info = nm)
  }
})

test_that("supplied folds without those arguments, and k, raise nothing", {
  skip_if_not_installed("ranger")
  u <- si_pts()
  f5 <- si_q(make_folds(u, k = 5, method = "random_kfold", seed = 1))
  cls <- "spatialkit_fold_args_ignored"
  # k is not checked: the package's own callers always pass it with folds.
  expect_length(si_classed(si_q(cv_rf(u, "z", "a", folds = f5, num_trees = 20, k = 10)), cls), 0L)
  expect_length(si_classed(si_q(cv_spatial(u, "z", "a", fit_fn = r2_fit, folds = f5, k = 3)), cls), 0L)
  # compare_models_cv() builds one fold set and hands it to each backend
  # with k: no warning from inside the package.
  expect_length(si_classed(si_q(compare_models_cv(u, "z", "a", models = "RF", k = 4,
                                                  rf_args = list(num_trees = 20),
                                                  quiet = TRUE)), cls), 0L)
  # block_size without folds is the documented use.
  expect_length(si_classed(si_q(compare_models_cv(u, "z", "a", models = "RF", k = 4,
                                                  block_size = 2500,
                                                  rf_args = list(num_trees = 20),
                                                  quiet = TRUE)), cls), 0L)
  # The sweep hands cv_spatial() its own folds at every size.
  expect_length(si_classed(si_q(cv_block_size_sweep(u, "z", "a", fit_fn = r2_fit, k = 3,
                                                    n_sizes = 2, include_random = FALSE,
                                                    sac = 2000, quiet = TRUE)), cls), 0L)
})

# ---- 3. compare_models(newdata = ): which residuals Moran's I is about ------

test_that("compare_models() says which residuals its Moran columns describe", {
  skip_if_not_installed("ranger")
  u <- si_pts(160, seed = 2)
  tr <- u[1:100, ]; te <- u[101:160, ]
  fits <- list(RF = fit_rf_model(tr, "z", "a", num_trees = 100, seed = 1),
               LM = lm_spatial_fit(tr, "z", "a"))
  if (requireNamespace("GWmodel", quietly = TRUE) && requireNamespace("sp", quietly = TRUE))
    fits$GWR <- si_q(suppressWarnings(fit_gwr_model(tr, "z", "a", adaptive = TRUE, bandwidth = 30)))
  with_new <- si_q(suppressWarnings(compare_models(fits, newdata = te)))
  without  <- si_q(suppressWarnings(compare_models(fits)))

  expect_true("resid_morans_basis" %in% names(with_new))
  basis <- stats::setNames(with_new$resid_morans_basis, with_new$model)
  # A forest's residuals() are out-of-bag, every other backend's in-sample.
  expect_identical(basis[["RF"]], "out-of-bag")
  expect_identical(basis[["LM"]], "in-sample")
  if ("GWR" %in% names(fits)) expect_identical(basis[["GWR"]], "in-sample")
  # The metrics are on newdata; the Moran columns are not, and never say so.
  expect_true(all(with_new$metric_basis == "newdata"))
  expect_false(any(with_new$resid_morans_basis == "newdata", na.rm = TRUE))
  # They are the same numbers with and without newdata, from the same rows.
  key <- c("resid_morans_I", "resid_morans_z", "resid_morans_p", "resid_morans_basis")
  expect_equal(with_new[match(names(fits), with_new$model), key],
               without[match(names(fits), without$model), key], ignore_attr = TRUE)
  for (nm in names(fits))
    expect_equal(with_new$resid_morans_I[with_new$model == nm],
                 suppressWarnings(residual_morans_i(fits[[nm]]))$observed)
})
