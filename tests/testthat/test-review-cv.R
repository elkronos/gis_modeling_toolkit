# tests/testthat/test-review-cv.R
# ---------------------------------------------------------------------------
# Regressions from the review of the CV runners, forward selection and
# compare_models_cv(): a parallel cv_gwr() that hung, partial fold failures
# that were only logged, and scores compared across different row sets.
# ---------------------------------------------------------------------------

.rv_gwr_pts <- function(n = 90, seed = 21) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               elev = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$price <- 10 + 2 * d$elev + rnorm(n, 0, 0.5)
  d
}


# ---------------------------------------------------------------------------
# cv_gwr(parallel = n) after a GWR fit
# ---------------------------------------------------------------------------

test_that("cv_gwr(parallel = 2) returns after a GWR has been fitted in the session", {
  # GWmodel is built with OpenMP.  Once fit_gwr_model() had run, GNU
  # libgomp's thread pool lived in the parent, every mclapply() child blocked
  # on a futex at its first OpenMP region, and cv_gwr(parallel = 2) never
  # returned.  Run in a child R process under `timeout` so a regression fails
  # here (status 124) instead of hanging the whole suite.
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  skip_if(parallel::detectCores() < 2L, "fewer than two cores: nothing forks")
  skip_if(!nzchar(Sys.which("timeout")), "coreutils `timeout` not found")

  # Load the same spatialkit this session has: the installed copy, or the
  # source tree pkgload::load_all() put in place.
  ns_path   <- getNamespaceInfo(asNamespace("spatialkit"), "path")
  installed <- file.exists(file.path(ns_path, "Meta", "package.rds"))
  if (!installed) skip_if_not_installed("pkgload")
  load_line <- if (installed)
    sprintf("suppressMessages(library(spatialkit, lib.loc = %s))",
            deparse(dirname(ns_path)))
  else
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))",
            deparse(ns_path))

  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    load_line,
    "logger::log_threshold(logger::FATAL, namespace = 'spatialkit', index = 2)",
    "set.seed(21); n <- 90",
    "d <- sf::st_as_sf(",
    "  data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),",
    "             elev = rnorm(n)),",
    "  coords = c('x', 'y'), crs = 32632)",
    "d$price <- 10 + 2 * d$elev + rnorm(n, 0, 0.5)",
    "invisible(suppressWarnings(suppressMessages(fit_gwr_model(d, 'price', 'elev'))))",
    "seq_res <- suppressWarnings(suppressMessages(",
    "  cv_gwr(d, 'price', 'elev', k = 3, seed = 7)))",
    "w <- character(0)",
    "par_res <- withCallingHandlers(suppressMessages(",
    "  cv_gwr(d, 'price', 'elev', k = 3, seed = 7, parallel = 2L)),",
    "  warning = function(x) {",
    "    w <<- c(w, conditionMessage(x)); invokeRestart('muffleWarning') })",
    "cat('RESULT', par_res$n_folds_succeeded,",
    "    isTRUE(all.equal(par_res$overall, seq_res$overall)),",
    "    any(grepl('OpenMP', w)), '\\n')"
  ), script)

  rscript <- file.path(R.home("bin"), "Rscript")
  # `timeout` makes itself a process-group leader and signals the group, so
  # any forked worker stuck in libgomp dies with the R process.  R_TESTS is
  # cleared because R CMD check points it at a startup file relative to the
  # tests directory, which the child would fail to source.
  out <- suppressWarnings(system2(
    "timeout", c("-k", "10", "150", shQuote(rscript), shQuote(script)),
    stdout = TRUE, stderr = TRUE, env = "R_TESTS="))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L,
                   label = paste0("child R exit status (124 = hung and killed); output:\n",
                                  paste(utils::tail(out, 15L), collapse = "\n")))
  # All three folds scored, identical to the sequential run, and a warning
  # that says why the folds did not fork.
  expect_true(any(grepl("^RESULT 3 TRUE TRUE", out)),
              label = paste(utils::tail(out, 15L), collapse = "\n"))
})


test_that("compare_models_cv(gwr_args = list(parallel = 2)) never forks GWR folds", {
  # The same hang reached through compare_models_cv().  mclapply() is
  # replaced by a stub that errors, so the old code fails here at once
  # instead of hanging on the GWR fits earlier tests have left in the session.
  skip_on_os("windows")
  skip_if_not_installed("GWmodel")
  skip_if_not_installed("sp")
  skip_if(parallel::detectCores() < 2L, "fewer than two cores: nothing forks")

  d <- .rv_gwr_pts()
  invisible(suppressWarnings(suppressMessages(fit_gwr_model(d, "price", "elev"))))
  ref <- suppressWarnings(suppressMessages(
    compare_models_cv(d, "price", "elev", models = "GWR", k = 3, seed = 7,
                      quiet = TRUE)))
  local_mocked_bindings(
    mclapply = function(...) stop("GWR folds were forked"),
    .package = "parallel")
  expect_warning(
    res <- suppressMessages(
      compare_models_cv(d, "price", "elev", models = "GWR", k = 3, seed = 7,
                        quiet = TRUE, gwr_args = list(parallel = 2L))),
    "cv_gwr\\(\\): `parallel` is ignored.*OpenMP")
  expect_identical(res$gwr_cv$n_folds_succeeded, 3L)
  expect_equal(res$overall, ref$overall)
})


# ---------------------------------------------------------------------------
# Some folds fail: a real R warning, not only a log line
# ---------------------------------------------------------------------------

# 100 points; zone "core" exists only east of x = 800, and those rows are fold
# 1 of a label vector, so no training set holds the level and every learner's
# predict() fails on that fold.
.rv_zone_pts <- function(n = 100, seed = 5) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- sf::st_as_sf(data.frame(x = x, y = y, w = rnorm(n)),
                    coords = c("x", "y"), crs = 3857)
  d$zone <- factor(ifelse(x > 800, "core", sample(c("a", "b"), n, TRUE)))
  d$z <- d$w + rnorm(n)
  d$lab <- ifelse(x > 800, 1L, sample(2:4, n, TRUE))
  d
}

.rv_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

test_that("cv_spatial() warns when some folds fail, naming them and the rows scored", {
  # overall was pooled over the surviving folds -- here without the one block
  # holding the level -- with only a logger line to say so.
  d <- .rv_zone_pts()
  n_core <- sum(d$lab == 1L)
  fitf <- function(tr) lm_spatial_fit(tr, "z", c("w", "zone"))
  expect_warning(
    cv <- cv_spatial(d, "z", c("w", "zone"), fit_fn = fitf, folds = d$lab),
    sprintf(paste0("^cv_spatial\\(\\): 1 of 4 fold\\(s\\) failed \\(fold 1: ",
                   "error\\), so `overall` pools the other folds only and ",
                   "covers %d of the 100 rows"), 100L - n_core))
  expect_identical(cv$overall$n_pred, 100L - n_core)
  expect_identical(cv$fold_status$status, c("error", "ok", "ok", "ok"))
})

test_that("cv_rf() raises the same warning under its own name", {
  skip_if_not_installed("ranger")
  d <- .rv_zone_pts()
  expect_warning(
    cv <- cv_rf(d, "z", c("w", "zone"), folds = d$lab, num_trees = 30, seed = 1),
    "^cv_rf\\(\\): 1 of 4 fold\\(s\\) failed \\(fold 1: error\\)")
  expect_identical(cv$n_folds_succeeded, 3L)
})

test_that("a fold dropped before fitting is not warned about twice", {
  # .remap_folds() already raises a warning for a fold that never reaches the
  # fitter; the new one names only the folds that failed there.
  d <- .rv_zone_pts()
  d$lab[d$lab == 4L] <- 5L                   # labels 1, 2, 3, 5 ...
  d$lab[which(d$lab == 3L)[1:3]] <- 4L       # ... and a small fold 4
  d$w[d$lab == 4L] <- NA                     # whose rows prep drops
  fitf <- function(tr) lm_spatial_fit(tr, "z", c("w", "zone"))
  r <- .rv_warnings(suppressMessages(
    cv_spatial(d, "z", c("w", "zone"), fit_fn = fitf, folds = d$lab)))
  expect_length(r$warnings, 2L)
  expect_match(r$warnings[1], "1 of 5 fold\\(s\\) dropped before fitting")
  expect_match(r$warnings[2], paste0("^cv_spatial\\(\\): 1 of 5 fold\\(s\\) ",
                                     "failed \\(fold 1: error\\), so"))
  expect_identical(r$value$fold_status$status,
                   c("error", "ok", "ok", "dropped", "ok"))
  # Only a fold dropped before fitting: the one warning it always had.
  d2 <- .rv_zone_pts()
  d2$w[d2$lab == 1L] <- NA
  r2 <- .rv_warnings(suppressMessages(
    cv_spatial(d2, "z", "w", fit_fn = function(tr) lm_spatial_fit(tr, "z", "w"),
               folds = d2$lab)))
  expect_length(r2$warnings, 1L)
  expect_match(r2$warnings, "dropped before fitting")
})


# ---------------------------------------------------------------------------
# select_features_forward(): every candidate set scored on the same rows
# ---------------------------------------------------------------------------

# The response is driven by `a`; land cover `lc` is noise whose level "C" sits
# only in the north-east corner, which is also where the response is noisiest.
# Under block CV the fold holding the corner cannot be predicted by any set
# containing `lc`, so such a set used to be scored on the 192 easier rows.
.rv_corner_pts <- function(n = 250, seed = 10) {
  set.seed(seed)
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  d <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], a = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  corner <- xy[, 1] > 800 & xy[, 2] > 800
  d$lc <- factor(ifelse(corner, "C", sample(c("A", "B"), n, TRUE)))
  d$z <- 2 * d$a + rnorm(n, 0, ifelse(corner, 12, 1))
  d
}

test_that("select_features_forward() does not let a set win by losing a fold", {
  # With lm and a null model: `a` was selected, then `lc` was accepted at step
  # 2 because {a, lc} scored RMSE 0.94 on 192 rows against 2.33 for {a} on 250.
  d <- .rv_corner_pts()
  fitf <- function(tr, vars) lm_spatial_fit(tr, "z", vars)
  r <- .rv_warnings(select_features_forward(d, "z", c("a", "lc"), fitf,
                                            k = 5, seed = 1, quiet = TRUE))
  sel <- r$value
  expect_identical(sel$selected, "a")
  expect_identical(sel$params$n_scored, 250L)
  h <- sel$history
  expect_true("n_pred" %in% names(h))
  expect_identical(h$n_pred[h$variable == "<none>"], 250L)
  expect_identical(h$n_pred[h$step == 1L & h$variable == "a"], 250L)
  # Every set holding `lc` predicted 192 rows and is NA, not a score.
  lc_rows <- h[h$variable == "lc", ]
  expect_true(all(lc_rows$n_pred == 192L))
  expect_true(all(is.na(lc_rows$score)))
  # One warning per step from the sweep, naming the set; cv_spatial()'s own
  # partial-failure warning is not repeated for every candidate.
  expect_length(r$warnings, 2L)
  expect_match(r$warnings[1], paste0("^select_features_forward\\(\\): step 1: ",
                                     "\\{lc\\} \\(192 of them predicted\\)"))
  expect_match(r$warnings[2], "step 2: \\{a, lc\\} \\(192 of them predicted\\)")
})

test_that("with no null model the rows every step-1 set predicted are the reference", {
  # RF refuses an empty predictor set, so the reference is the union of the
  # step-1 candidates' rows: `lc` is still NA and `a` is chosen.  It used to
  # be the other way round (RMSE 2.35 on 192 rows beat 2.63 on 250).
  skip_if_not_installed("ranger")
  d <- .rv_corner_pts()
  fitf <- function(tr, vars) fit_rf_model(tr, "z", vars, num_trees = 100, seed = 1)
  expect_warning(
    sel <- select_features_forward(d, "z", c("a", "lc"), fitf, k = 5, seed = 1,
                                   quiet = TRUE, max_vars = 1),
    "step 1: \\{lc\\}")
  expect_identical(sel$selected, "a")
  expect_identical(sel$params$n_scored, 250L)
  expect_identical(sel$history$n_pred, c(250L, 192L))
})
