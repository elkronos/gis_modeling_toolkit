# tests/testthat/test-review2-cv.R
# ---------------------------------------------------------------------------
# Regressions from the second review of the CV runners, cv_bayes()'s
# coverage levels and plot_calibration().  Each test names the finding it
# closes.  The Bayesian backend is mocked with helper-lmfit.R's lm fit, whose
# predict(draws = TRUE) returns a draw matrix, so no Stan toolchain is
# needed.
# ---------------------------------------------------------------------------

.r2_pts <- function(n = 90, seed = 1, extent = 1000) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = runif(n, 0, extent), y = runif(n, 0, extent), w = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 5 + 2 * d$w + rnorm(n, 0, 0.5)
  d
}

.r2_lm <- function(tr) lm_spatial_fit(tr, "z", "w")

.r2_mock_bayes <- function(env = parent.frame()) {
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123)
      lm_spatial_fit(data_sf, response_var, predictor_vars),
    .package = "spatialkit", .env = env)
}


# ---------------------------------------------------------------------------
# FP9 / plotting L10: coverage levels named at full precision and carried
# explicitly; coverage_levels validated
# ---------------------------------------------------------------------------

test_that("cv_bayes() keeps close coverage levels apart and plots them where they belong", {
  skip_if_not_installed("ggplot2")
  .r2_mock_bayes()
  d <- .r2_pts(60, seed = 2)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  lv <- c(0.5, 0.975, 0.985, 0.995)
  cv <- suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = 1,
                                  coverage_levels = lv))
  # 0.975 and 0.985 were both "coverage_98", one overwriting the other, and
  # 0.995 became "coverage_100".
  cols <- c("coverage_50", "coverage_97.5", "coverage_98.5", "coverage_99.5")
  expect_true(all(cols %in% names(cv$fold_metrics)))
  expect_equal(cv$coverage_levels, stats::setNames(lv, cols))
  expect_setequal(setdiff(names(cv$predictive_coverage), "mean_CRPS"), cols)
  # The plot reads the nominal level from coverage_levels, not a rounded name.
  p <- plot_calibration(cv)
  pooled <- p$layers[[4]]$data
  expect_equal(sort(pooled$nominal), lv)
  # The default levels keep their familiar names.
  cv0 <- suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = 1))
  expect_named(cv0$coverage_levels, c("coverage_50", "coverage_80", "coverage_95"))
})

test_that("cv_bayes() refuses coverage levels outside (0, 1) or given twice", {
  # c(50, 80, 95) made alpha negative, quantile() threw inside the fold
  # extras, and every fold silently lost CRPS, n_draws and coverage (MED-2).
  # Mocked, so a regression here costs a quick lm CV rather than Stan fits.
  .r2_mock_bayes()
  d <- .r2_pts(30, seed = 3)
  expect_error(cv_bayes(d, "z", "w", coverage_levels = c(50, 80, 95)),
               "strictly between 0 and 1.*divide by 100: c\\(0.5, 0.8, 0.95\\)")
  expect_error(cv_bayes(d, "z", "w", coverage_levels = c(0.5, 1.2)),
               "strictly between 0 and 1")
  expect_error(cv_bayes(d, "z", "w", coverage_levels = c(0.5, NA)),
               "strictly between 0 and 1")
  expect_error(cv_bayes(d, "z", "w", coverage_levels = c(0.9, 0.5, 0.9)),
               "gives the level 0.9 more than once")
})


# ---------------------------------------------------------------------------
# plotting L5: the all-NA coverage message names both causes
# ---------------------------------------------------------------------------

test_that("plot_calibration() names compute_pred_intervals = FALSE as a cause", {
  skip_if_not_installed("ggplot2")
  .r2_mock_bayes()
  d <- .r2_pts(45, seed = 4)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  cv <- suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = 1,
                                  compute_pred_intervals = FALSE))
  expect_error(plot_calibration(cv), "compute_pred_intervals = FALSE")
})


# ---------------------------------------------------------------------------
# cv-runners MED-2: a fold_info_fn that throws is logged and recorded
# ---------------------------------------------------------------------------

test_that("a fold_info_fn that throws is logged and named in fold_status", {
  d <- .r2_pts(60, seed = 5)
  lab <- rep(1:3, each = 20)
  info <- function(fit, test_sf, y, yhat) {
    if (21L %in% test_sf$..row_id) stop("boom on the middle fold")
    list(n_coef = length(stats::coef(fit$engine)))
  }
  lines <- capture_spatialkit_log(
    cv <- cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab,
                     fold_info_fn = info))
  expect_equal(cv$n_folds_succeeded, 3L)
  expect_equal(cv$fold_status$status, rep("ok", 3))
  expect_match(cv$fold_status$message[2], "fold_info_fn failed: boom on the middle fold")
  expect_identical(cv$fold_status$message[c(1, 3)], c("", ""))
  expect_equal(cv$fold_metrics$n_coef, c(2, NA, 2))
  expect_true(log_has(lines, "fold 2: fold_info_fn failed: boom"))
})

test_that("cv_bayes() keeps gp_k and n_draws when coverage cannot be computed", {
  # quantile() refuses a draw matrix with an NA; that used to throw away
  # every extra of the fold, n_draws included, with nothing logged.
  registerS3method("predict", "r2_nadraws",
                   function(object, newdata = NULL, draws = FALSE, ...) {
                     mu <- predict.lmsurf_fit(object, newdata = newdata)
                     if (!isTRUE(draws)) return(mu)
                     m <- rbind(mu - 1, mu, mu + 1)
                     m[1, 1] <- NA
                     m
                   })
  local_mocked_bindings(
    fit_bayesian_spatial_model = function(data_sf, response_var, predictor_vars,
                                          ..., seed = 123) {
      fit <- lm_spatial_fit(data_sf, response_var, predictor_vars)
      fit$info$gp_k <- 7L
      class(fit) <- c("r2_nadraws", class(fit))
      fit
    },
    .package = "spatialkit")
  d <- .r2_pts(45, seed = 6)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  lines <- capture_spatialkit_log(
    cv <- suppressWarnings(cv_bayes(d, "z", "w", folds = f, seed = 1)))
  expect_equal(cv$fold_metrics$n_draws, rep(3L, 3))
  expect_equal(cv$fold_metrics$gp_k, rep(7L, 3))
  expect_true(all(is.na(cv$fold_metrics$coverage_95)))
  expect_true(log_has(lines, "coverage and CRPS could not be computed"))
})


# ---------------------------------------------------------------------------
# cv-runners MED-3: ..per_row on some folds only
# ---------------------------------------------------------------------------

test_that("cv_spatial() stacks folds whose fold_info_fn gave ..per_row on some folds only", {
  d <- .r2_pts(90, seed = 7)
  lab <- rep(1:3, each = 30)
  info <- function(fit, test_sf, y, yhat) {
    if (31L %in% test_sf$..row_id) return(list())          # fold 2: none
    list(..per_row = data.frame(abs_err = abs(y - yhat)))
  }
  cv <- cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab,
                   fold_info_fn = info)
  expect_equal(nrow(cv$predictions), 90L)
  p <- cv$predictions
  expect_true(all(is.na(p$abs_err[p$fold == 2L])))
  expect_equal(p$abs_err[p$fold != 2L], abs(p$y - p$yhat)[p$fold != 2L])

  # One of the wrong length is dropped, and the log says so.
  bad <- function(fit, test_sf, y, yhat) list(..per_row = data.frame(e = 1:2))
  lines <- capture_spatialkit_log(
    cv2 <- cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab,
                      fold_info_fn = bad))
  expect_false("e" %in% names(cv2$predictions))
  expect_true(log_has(lines, "`..per_row` has 2 rows for 30 test rows"))
})


# ---------------------------------------------------------------------------
# cv-runners MED-4: numeric fold labels numbered in numeric order
# ---------------------------------------------------------------------------

test_that("numeric fold labels keep their numeric order, a factor its own", {
  d <- .r2_pts(240, seed = 8)
  lab <- as.integer(cut(sf::st_coordinates(d)[, 1], 12))
  cv <- cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab)
  p <- cv$predictions
  # Output fold i is the user's label i; with string order fold 2 was label 10.
  expect_identical(as.integer(p$fold), lab[p$..row_id])
  expect_identical(cv$fold_metrics$n_test,
                   as.integer(table(factor(lab, levels = 1:12))))
  # A factor's own level order is the numbering.
  rev_lab <- factor(lab, levels = 12:1)
  cv_f <- cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = rev_lab)
  pf <- cv_f$predictions
  expect_identical(as.integer(pf$fold), 13L - lab[pf$..row_id])
})


# ---------------------------------------------------------------------------
# cv-runners MED-5: parallel runs match sequential ones on failures
# ---------------------------------------------------------------------------

test_that("an error that stops a sequential run stops a parallel one too", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "fewer than two cores: nothing forks")
  d <- .r2_pts(80, seed = 9)
  lab <- rep(1:4, each = 20)
  # A non-scalar extra on fold 2.  Sequentially this died with R's
  # "replacement has 0 rows"; in parallel it took fold 4 down with it and
  # the run carried on.
  info <- function(fit, test_sf, y, yhat)
    list(k = if (21L %in% test_sf$..row_id) numeric(0) else 1)
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab,
                          fold_info_fn = info),
               "element 'k' is not a single value")
  expect_error(suppressMessages(
    cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab,
               fold_info_fn = info, parallel = 2)),
    "element 'k' is not a single value.*fold 2, raised in a parallel worker")
  # The documented shape error of a user metric, which parallel runs dropped.
  dup <- function(y, yhat) c(a = 1, a = 2)
  expect_error(suppressMessages(
    cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = lab, metrics = dup,
               parallel = 2)),
    "duplicated names: a")
})

test_that("a killed parallel worker costs its own fold only, as a worker_error", {
  # Prescheduled mclapply() gave each core a chunk of folds and returned NULL
  # for every fold of a core that died: folds 2 AND 4 came back "skipped".
  # Run in a child R under `timeout`, since the fit kills its own process.
  skip_on_cran()
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "fewer than two cores: nothing forks")
  skip_if(!nzchar(Sys.which("timeout")), "coreutils `timeout` not found")
  ns_path   <- getNamespaceInfo(asNamespace("spatialkit"), "path")
  installed <- file.exists(file.path(ns_path, "Meta", "package.rds"))
  if (!installed) skip_if_not_installed("pkgload")
  load_line <- if (installed)
    sprintf("suppressMessages(library(spatialkit, lib.loc = %s))",
            deparse(dirname(ns_path)))
  else
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))",
            deparse(ns_path))
  helper <- normalizePath(test_path("helper-lmfit.R"))
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    load_line,
    sprintf("source(%s)", deparse(helper)),
    "logger::log_threshold(logger::FATAL, namespace = 'spatialkit', index = 2)",
    "set.seed(9); n <- 80",
    "d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),",
    "  w = rnorm(n)), coords = c('x', 'y'), crs = 32632)",
    "d$z <- 2 * d$w + rnorm(n)",
    "parent <- Sys.getpid()",
    "fit_fn <- function(tr) {",
    "  if (!(21L %in% tr$..row_id) && Sys.getpid() != parent)",
    "    tools::pskill(Sys.getpid(), tools::SIGKILL)",
    "  lm_spatial_fit(tr, 'z', 'w')",
    "}",
    "cv <- suppressWarnings(suppressMessages(cv_spatial(d, 'z', 'w',",
    "  fit_fn = fit_fn, folds = rep(1:4, each = 20), parallel = 2)))",
    "cat('RESULT', cv$fold_status$status, '\\n')"
  ), script)
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(
    "timeout", c("-k", "10", "120", shQuote(rscript), shQuote(script)),
    stdout = TRUE, stderr = TRUE, env = "R_TESTS="))
  res <- grep("^RESULT", out, value = TRUE)
  expect_identical(trimws(res), "RESULT ok worker_error ok ok",
                   label = paste(utils::tail(out, 15L), collapse = "\n"))
})


# ---------------------------------------------------------------------------
# cv-runners LOW-7: a user metric or a fold_info_fn may not overwrite a column
# ---------------------------------------------------------------------------

test_that("user metrics named like the backend's extras are refused", {
  .r2_mock_bayes()
  d <- .r2_pts(45, seed = 10)
  f <- make_folds(d, k = 3, method = "random_kfold", seed = 1)
  # predictive_coverage used to report the user's 999 as the mean CRPS.
  expect_error(suppressWarnings(
    cv_bayes(d, "z", "w", folds = f, metrics = function(y, yhat) c(CRPS = 999))),
    "already columns of the metrics frames: CRPS")
  expect_error(suppressWarnings(
    cv_bayes(d, "z", "w", folds = f,
             metrics = function(y, yhat) c(coverage_95 = 0.1))),
    "already columns of the metrics frames: coverage_95")
  # mean_CRPS is written by compare_models_cv() over a user column.
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = f,
                          metrics = function(y, yhat) c(mean_CRPS = 1)),
               "already columns of the metrics frames: mean_CRPS")
  # A user metric may not overwrite a fold_info_fn extra either ...
  info <- function(fit, test_sf, y, yhat) list(tuned = 1)
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = f,
                          fold_info_fn = info,
                          metrics = function(y, yhat) c(tuned = 2)),
               "already columns of the metrics frames: tuned")
  # ... nor a fold_info_fn a built-in column.
  expect_error(cv_spatial(d, "z", "w", fit_fn = .r2_lm, folds = f,
                          fold_info_fn = function(fit, test_sf, y, yhat)
                            list(RMSE = -5)),
               "`fold_info_fn` returned names that are already columns of fold_metrics: RMSE")
})


# ---------------------------------------------------------------------------
# gaps G2.4: saved folds on lon/lat polygons survive an sf_use_s2() toggle
# ---------------------------------------------------------------------------

test_that("folds on lon/lat polygons are accepted after sf_use_s2() is toggled", {
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc$z <- nc$BIR74 / 1000; nc$a <- nc$NWBIR74 / 1000
  old <- sf::sf_use_s2()
  withr::defer(suppressMessages(sf::sf_use_s2(old)))
  fit_fn <- function(tr) lm_spatial_fit(tr, "z", "a")
  run <- function(f) suppressWarnings(suppressMessages(
    cv_spatial(nc, "z", "a", fit_fn = fit_fn, folds = f)))

  # suppressWarnings(): sf warns that point-on-surface is approximate on
  # lon/lat data with s2 off, which is beside the point here.
  suppressMessages(sf::sf_use_s2(TRUE))
  f_on <- suppressWarnings(make_folds(nc, k = 3, method = "random_kfold", seed = 1))
  suppressMessages(sf::sf_use_s2(FALSE))
  f_off <- suppressWarnings(make_folds(nc, k = 3, method = "random_kfold", seed = 1))
  expect_equal(run(f_on)$n_folds_succeeded, 3L)       # built on, used off
  suppressMessages(sf::sf_use_s2(TRUE))
  expect_equal(run(f_off)$n_folds_succeeded, 3L)      # built off, used on
  expect_identical(sf::sf_use_s2(), TRUE)             # the probe restored it

  # Different data is still refused.
  shuffled <- nc[c(2:100, 1), ]
  expect_error(suppressWarnings(suppressMessages(
    cv_spatial(shuffled, "z", "a", fit_fn = fit_fn, folds = f_on))),
    "built from different data")

  # A probe saved before `kind` existed is checked the old way, not refused.
  legacy <- f_on
  legacy$params$row_probe$kind <- NULL
  sf_ids <- nc; sf_ids$..row_id <- seq_len(nrow(nc))
  old_probe <- spatialkit:::.fold_row_probe(sf_ids, legacy = TRUE)
  legacy$params$row_probe$x <- old_probe$x
  legacy$params$row_probe$y <- old_probe$y
  expect_equal(run(legacy)$n_folds_succeeded, 3L)
})


# ---------------------------------------------------------------------------
# first-pass FP15: cv_gwr() no longer calls the dead .validate_kernel()
# ---------------------------------------------------------------------------

test_that("cv_gwr() validates its kernel with match.arg() alone", {
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  src <- paste(deparse(body(cv_gwr)), collapse = " ")
  expect_false(grepl(".validate_kernel", src, fixed = TRUE))
  expect_error(cv_gwr(.r2_pts(20), "z", "w", kernel = "Gaussian"),
               "should be one of")
})
