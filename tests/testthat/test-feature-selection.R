# tests/testthat/test-feature-selection.R
# ---------------------------------------------------------------------------
# select_features_forward().
#
# The statistical point is the inner loop being spatially blocked.  Random
# inner folds inside blocked outer folds select variables that look predictive
# only because nearby points leak between train and test; the outer loop then
# reports honest-looking numbers for a dishonestly chosen feature set.
# ---------------------------------------------------------------------------

# Two of five candidates carry real signal; three are noise.
fs_data <- function(n = 150, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  a <- rnorm(n); b <- rnorm(n); c1 <- rnorm(n); d <- rnorm(n); e <- rnorm(n)
  z <- 3 * a + 2 * b + rnorm(n, 0, 0.3)          # only a and b matter
  out <- sf::st_as_sf(
    data.frame(x = x, y = y, z = z, a = a, b = b, c1 = c1, d = d, e = e),
    coords = c("x", "y"), crs = 3857)
  out$..row_id <- seq_len(nrow(out))
  out
}

fs_fit <- function(tr, vars) lm_spatial_fit(tr, "z", vars)


test_that("the informative variables are selected first", {
  dat <- fs_data()
  sel <- select_features_forward(dat, "z", c("a", "b", "c1", "d", "e"),
                                 fs_fit, k = 4, seed = 1, quiet = TRUE)
  expect_true(all(c("a", "b") %in% sel$selected))
  expect_equal(sel$selected[1:2], c("a", "b"))   # a has the larger coefficient
})

test_that("noise variables are not accumulated", {
  dat <- fs_data()
  sel <- select_features_forward(dat, "z", c("a", "b", "c1", "d", "e"),
                                 fs_fit, k = 4, seed = 1, quiet = TRUE)
  # Stopping is what keeps this from being a full sweep dressed up as selection.
  expect_lt(length(sel$selected), 5L)
})

test_that("tol controls how readily variables are accepted", {
  dat <- fs_data()
  loose <- select_features_forward(dat, "z", c("a", "b", "c1", "d", "e"),
                                   fs_fit, k = 4, tol = 0, seed = 1, quiet = TRUE)
  tight <- select_features_forward(dat, "z", c("a", "b", "c1", "d", "e"),
                                   fs_fit, k = 4, tol = 1e6, seed = 1, quiet = TRUE)
  expect_lte(length(tight$selected), length(loose$selected))
  # An impossible tol now selects NOTHING, not one variable.  `best` starts at
  # the null (intercept-only) model's score rather than at +/-Inf, so the first
  # step is tested like every other one -- previously it was accepted whatever
  # its score, which is how a pure-noise response could still come back with a
  # "predictive" variable.
  expect_length(tight$selected, 0L)
})

test_that("max_vars caps the selection", {
  dat <- fs_data()
  sel <- select_features_forward(dat, "z", c("a", "b", "c1", "d", "e"),
                                 fs_fit, k = 4, max_vars = 1, seed = 1, quiet = TRUE)
  expect_length(sel$selected, 1L)
})

test_that("R2 is maximised while RMSE is minimised", {
  dat <- fs_data()
  by_rmse <- select_features_forward(dat, "z", c("a", "b", "c1"), fs_fit,
                                     k = 4, metric = "RMSE", seed = 1, quiet = TRUE)
  by_r2   <- select_features_forward(dat, "z", c("a", "b", "c1"), fs_fit,
                                     k = 4, metric = "R2", seed = 1, quiet = TRUE)
  expect_equal(by_rmse$selected[1], by_r2$selected[1])
  expect_true(by_r2$score > 0 && by_r2$score <= 1)
  expect_true(by_rmse$score > 0)
})

test_that("inner folds are spatially blocked by default", {
  dat <- fs_data()
  sel <- select_features_forward(dat, "z", c("a", "b"), fs_fit, k = 4,
                                 seed = 1, quiet = TRUE)
  expect_equal(sel$params$method, "block_kfold")
})

test_that("random inner folds are permitted but warned about", {
  # Loudly, because the failure mode is silent: an outer blocked loop will
  # report clean numbers for a feature set chosen through leakage.
  dat <- fs_data()
  lines <- capture_spatialkit_log(
    select_features_forward(dat, "z", c("a", "b"), fs_fit, k = 4,
                            method = "random_kfold", seed = 1, quiet = TRUE)
  )
  expect_true(log_has(lines, "leak spatial autocorrelation"))
})

test_that("the fit budget stops a runaway sweep", {
  dat <- fs_data()
  expect_error(
    select_features_forward(dat, "z", c("a", "b", "c1", "d", "e"), fs_fit,
                            k = 10, max_fits = 5, seed = 1, quiet = TRUE),
    "above max_fits"
  )
})

test_that("history records every candidate at every step", {
  dat <- fs_data()
  sel <- select_features_forward(dat, "z", c("a", "b", "c1"), fs_fit,
                                 k = 4, seed = 1, quiet = TRUE)
  h <- sel$history
  expect_true(all(c("step", "variable", "score") %in% names(h)))
  expect_equal(sum(h$step == 1L), 3L)            # all three tried at step 1
  expect_true(all(h$variable %in% c("<none>", "a", "b", "c1")))
  # Step 0 is the null model, the baseline the first step has to beat.
  expect_equal(sum(h$step == 0L), 1L)
  expect_equal(h$variable[h$step == 0L], "<none>")
  expect_true(is.finite(h$score[h$step == 0L]))
})

test_that("bad input is rejected", {
  dat <- fs_data()
  expect_error(select_features_forward(sf::st_drop_geometry(dat), "z", "a", fs_fit),
               "must be an sf object")
  expect_error(select_features_forward(dat, "z", "nope", fs_fit),
               "absent from")
  expect_error(select_features_forward(dat, "z", character(0), fs_fit),
               "no candidate")
  expect_error(select_features_forward(dat, "z", "a", fit_fn = "not a function"),
               "must be a function")
})

test_that("selection is reproducible from the seed", {
  dat <- fs_data()
  a <- select_features_forward(dat, "z", c("a", "b", "c1"), fs_fit, k = 4,
                               seed = 42, quiet = TRUE)
  b <- select_features_forward(dat, "z", c("a", "b", "c1"), fs_fit, k = 4,
                               seed = 42, quiet = TRUE)
  expect_equal(a$selected, b$selected)
  expect_equal(a$score, b$score)
})


test_that("auto_range is off by default and recorded either way", {
  dat <- fs_data()
  sel <- select_features_forward(dat, "z", c("a", "b"), fs_fit, k = 4,
                                 seed = 1, quiet = TRUE)
  expect_false(sel$params$auto_range)
  # With it on, the argument is recorded whatever the estimate did.  Without
  # gstat no range can be estimated and make_folds() keeps the geometric
  # blocks with a warning, so this part needs no optional package; the
  # positive proof that the blocks get sized from the range is the next test.
  sel_on <- suppressMessages(suppressWarnings(
    select_features_forward(dat, "z", c("a", "b"), fs_fit, k = 4,
                            seed = 1, quiet = TRUE, auto_range = TRUE)))
  expect_true(sel_on$params$auto_range)
})

test_that("auto_range sizes the inner blocks from the estimated range", {
  skip_if_not_installed("gstat")
  # A response with a real range: make_folds() announces the range it used as
  # the minimum block size.  A short range (exponential parameter 60,
  # effective range ~220 on this draw): a longer one estimates past half the
  # extent, and one block is not a split (make_folds() refuses it, correctly).
  set.seed(3)
  n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  sp <- as.numeric(t(chol(exp(-d / 60) + diag(1e-4, n))) %*% rnorm(n))
  spat <- sf::st_as_sf(data.frame(x = x, y = y, a = rnorm(n), b = rnorm(n)),
                       coords = c("x", "y"), crs = 3857)
  spat$z <- spat$a + sp
  spat$..row_id <- seq_len(n)
  msgs <- character(0)
  sel <- withCallingHandlers(
    suppressWarnings(
      select_features_forward(spat, "z", c("a", "b"), fs_fit, k = 4,
                              seed = 1, quiet = TRUE, auto_range = TRUE)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("using as minimum block size", msgs)))
  expect_true(sel$params$auto_range)
})


test_that("a fit_fn that ignores its second argument is called out, not returned as an empty selection", {
  # A learner written for cv_spatial() takes one argument; `function(train_sf,
  # ...)` swallows `vars` and fits the same model every time.  The training
  # layer still carries every column, so nothing errors: every candidate
  # scores exactly what the null model scored, nothing improves on it, and the
  # result was an empty `selected` and an NA `score` with no word about why.
  pts <- fs_data(200)
  one_arg <- function(train_sf, ...) lm_spatial_fit(train_sf, "z", "a")
  expect_warning(
    res <- select_features_forward(pts, "z", c("a", "b", "c1"), fit_fn = one_arg,
                                   k = 3, method = "random_kfold", seed = 1,
                                   quiet = TRUE),
    "ignoring its second argument")
  # The result is still returned -- the warning is the diagnosis, not a stop.
  expect_length(res$selected, 0L)
  expect_true(is.na(res$score))
  # Every step-1 score equals the null score to the last digit: the signature.
  h <- res$history
  expect_equal(unique(h$score[h$step == 1L]), h$score[h$step == 0L])

  # A learner that honours `vars` produces different scores per candidate and
  # must not trip it.
  two_arg <- function(tr, vars) lm_spatial_fit(tr, "z", vars)
  expect_no_warning(
    res2 <- select_features_forward(pts, "z", c("a", "b", "c1"), fit_fn = two_arg,
                                    k = 3, method = "random_kfold", seed = 1,
                                    quiet = TRUE))
  expect_true("a" %in% res2$selected)
})
