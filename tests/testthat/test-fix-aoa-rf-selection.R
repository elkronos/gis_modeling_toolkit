# tests/testthat/test-fix-aoa-rf-selection.R
# ---------------------------------------------------------------------------
# select_features_forward(): a baseline for every learner, and a step rule
# that does not accept fold noise.
#
# (a) fit_rf_model() refuses an empty predictor set, and fit_gwr_model() did
#     when this was written (it now fits the local mean for character(0)), so
#     with either as the learner the sweep had no step-0 score and the first
#     variable was accepted unconditionally: on a response unrelated to
#     every candidate a forest sweep selected 2.6 to 3.2 variables in 10
#     draws of 10.  The training-fold mean is now the baseline for such a
#     learner (it needs no model), the first variable is provisional, and
#     the finished set is kept only if it beats that baseline.
#
# (b) With tol = 0 any gain was accepted, however far inside the
#     fold-to-fold noise.  The default rule now also asks the mean per-fold
#     gain to exceed one standard error of itself; rule = "tol" is the old
#     behaviour.
#
# (c) Under the fold-mean baseline the first variable is added without a
#     test, and the test of the finished set does not show that it belongs
#     there.  A kept path of two or more is scored once more without it,
#     and it is dropped when the set is no better with it by the rule.
#
# (d) A GWR sweep on the fold-mean baseline is pointed to a `null_fn` that
#     can be written with exported functions: the kernel-weighted local
#     mean of the help page's example.
# ---------------------------------------------------------------------------

sel_noise <- function(n = 200, seed = 1) {
  set.seed(1000 + seed)
  d <- data.frame(x = stats::runif(n, 0, 2000), y = stats::runif(n, 0, 2000),
                  c1 = stats::rnorm(n), c2 = stats::rnorm(n),
                  c3 = stats::rnorm(n), c4 = stats::runif(n))
  d$z <- stats::rnorm(n)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
sel_signal <- function(n = 250, seed = 1, b = 0.5) {
  set.seed(2000 + seed)
  d <- data.frame(x = stats::runif(n, 0, 2000), y = stats::runif(n, 0, 2000),
                  a1 = stats::rnorm(n), a2 = stats::rnorm(n),
                  n1 = stats::rnorm(n), n2 = stats::rnorm(n), n3 = stats::rnorm(n))
  d$z <- b * d$a1 + b * d$a2 + stats::rnorm(n)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
sel_rf <- function(tr, v)
  fit_rf_model(tr, "z", v, num_trees = 100, importance = "none")
sel_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# The fold-mean RMSE on the folds the sweep builds, by hand.
sel_fold_mean_rmse <- function(d, cands, k, seed) {
  f <- sel_quiet(make_folds(d, k = k, method = "block_kfold", seed = seed,
                            response_var = "z", predictor_vars = cands))
  e <- unlist(lapply(f$folds, function(s) d$z[s$test] - mean(d$z[s$train])))
  sqrt(mean(e^2))
}

# A learner that, like the forest backend, refuses an empty set.
sel_lm_refusing <- function(tr, v) {
  if (!length(v)) stop("this backend needs at least one predictor", call. = FALSE)
  lm_spatial_fit(tr, "z", v)
}
sel_lm <- function(tr, v) lm_spatial_fit(tr, "z", v)


# --- (a) the baseline --------------------------------------------------------

test_that("a forest sweep has a baseline, and noise is not selected", {
  skip_if_not_installed("ranger")
  cands <- paste0("c", 1:4)
  for (s in 1:3) {
    d   <- sel_noise(seed = s)
    sel <- sel_quiet(select_features_forward(d, "z", cands, sel_rf, k = 4,
                                             seed = s, quiet = TRUE))
    # It used to return one to four of these every time.
    expect_identical(sel$selected, character(0))
    expect_true(is.na(sel$score))
    h <- sel$history
    expect_identical(h$step[1], 0L)
    expect_identical(h$variable[1], "<none>")
    expect_identical(sel$baseline$type, "fold mean")
    expect_identical(sel$params$baseline, "fold mean")
    expect_equal(sel$baseline$score, sel_fold_mean_rmse(d, cands, 4, s),
                 tolerance = 1e-10)
    expect_identical(h$score[1], sel$baseline$score)
    # The path the sweep reached is kept for inspection, with its verdict.
    expect_gt(length(sel$baseline$path), 0L)
    expect_false(sel$baseline$kept)
    expect_true(all(sel$baseline$path %in% cands))
  }
})

test_that("a forest sweep still finds real predictors", {
  skip_if_not_installed("ranger")
  # z = 0.5 a1 + 0.5 a2 + N(0, 1): a true R2 of 0.33.  No single variable
  # beats the fold mean here, which is why the first step is provisional.
  d   <- sel_signal(seed = 1)
  sel <- sel_quiet(select_features_forward(d, "z", c("a1", "a2", "n1", "n2", "n3"),
                                           sel_rf, k = 4, seed = 1, quiet = TRUE))
  expect_true(all(c("a1", "a2") %in% sel$selected))
  expect_identical(sel$baseline$type, "fold mean")
  expect_true(sel$baseline$kept)
  expect_gt(sel$baseline$gain, 0)
  expect_gt(sel$baseline$fold_gain, sel$baseline$fold_gain_se)
  expect_identical(sel$baseline$path, sel$selected)
  expect_lt(sel$score, sel$baseline$score)
})

test_that("the fold-mean baseline and its verdict are announced", {
  d <- sel_noise(n = 150, seed = 4)
  msgs <- character(0)
  sel <- withCallingHandlers(
    suppressWarnings(select_features_forward(d, "z", c("c1", "c2"),
                                             sel_lm_refusing, k = 3, seed = 4)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("baseline is the training-fold mean", msgs)))
  expect_true(any(grepl("The first variable is provisional", msgs)))
  expect_true(any(grepl("\\[provisional\\]", msgs)))
  expect_false(any(grepl("accepted unconditionally", msgs)))
  if (!length(sel$selected))
    expect_true(any(grepl("does not beat the training-fold mean by the rule", msgs)))
  # The probe of the empty set stays out of the console log.
  lines <- capture_spatialkit_log(
    sel_quiet(select_features_forward(d, "z", c("c1", "c2"), sel_lm_refusing,
                                      k = 3, seed = 4, quiet = TRUE)),
    level = logger::INFO)
  expect_false(log_has(lines, "fit failed; skipping"))
})

test_that("a learner that fits its own null model is tested against it at step 1", {
  d <- sel_signal(n = 200, seed = 2, b = 1)
  cands <- c("a1", "a2", "n1")
  own <- sel_quiet(select_features_forward(d, "z", cands, sel_lm, k = 4,
                                           seed = 2, quiet = TRUE))
  expect_identical(own$baseline$type, "null model")
  expect_true(all(c("a1", "a2") %in% own$selected))
  # The intercept-only lm IS the training-fold mean, so a refusing learner
  # falls back to the same step-0 score.
  ref <- sel_quiet(select_features_forward(d, "z", cands, sel_lm_refusing, k = 4,
                                           seed = 2, quiet = TRUE))
  expect_identical(ref$baseline$type, "fold mean")
  expect_equal(ref$baseline$score, own$baseline$score, tolerance = 1e-10)
  expect_identical(ref$selected, own$selected)
  # And `null_fn` supplies a no-predictor model for a learner without one.
  nf <- sel_quiet(select_features_forward(
    d, "z", cands, sel_lm_refusing, k = 4, seed = 2, quiet = TRUE,
    null_fn = function(tr) lm_spatial_fit(tr, "z", character(0))))
  expect_identical(nf$baseline$type, "null_fn")
  expect_equal(nf$history, own$history)
  expect_identical(nf$selected, own$selected)
})

test_that("bad rule arguments are refused by name", {
  d <- sel_noise(n = 60, seed = 5)
  expect_error(select_features_forward(d, "z", "c1", sel_lm, se_mult = -1),
               "`se_mult` must be")
  expect_error(select_features_forward(d, "z", "c1", sel_lm, se_mult = c(1, 2)),
               "`se_mult` must be")
  expect_error(select_features_forward(d, "z", "c1", sel_lm, rule = "aic"),
               "should be one of")
  expect_error(select_features_forward(d, "z", "c1", sel_lm, null_fn = "mean"),
               "`null_fn` must be NULL or a function")
  expect_error(
    sel_quiet(select_features_forward(d, "z", "c1", sel_lm, k = 3,
                                      null_fn = function(tr) stop("no"))),
    "`null_fn` gave no cross-validated prediction")
})


# --- (b) the step rule -------------------------------------------------------

test_that("history carries what the rule read, and every decision follows from it", {
  d <- sel_signal(n = 200, seed = 3, b = 1)
  cands <- c("a1", "a2", "n1", "n2", "n3")
  sel <- sel_quiet(select_features_forward(d, "z", cands, sel_lm, k = 4,
                                           seed = 3, quiet = TRUE))
  h <- sel$history
  expect_true(all(c("step", "variable", "score", "n_pred", "gain", "fold_gain",
                    "fold_gain_se") %in% names(h)))
  expect_identical(sel$params$rule, "se")
  expect_identical(sel$params$se_mult, 1)
  # gain is the pooled improvement over the set the candidate was added to.
  s1 <- h[h$step == 1L, ]
  expect_equal(s1$gain, h$score[h$step == 0L] - s1$score, tolerance = 1e-12)
  # Per-fold scores: one row per set and fold, and they reproduce fold_gain.
  fs <- sel$fold_scores
  expect_named(fs, c("step", "variable", "fold", "score"))
  expect_identical(nrow(fs), 4L * nrow(h))
  f0 <- fs[fs$step == 0L, ]
  fa <- fs[fs$step == 1L & fs$variable == sel$selected[1], ]
  d1 <- f0$score[match(fa$fold, f0$fold)] - fa$score
  row1 <- s1[s1$variable == sel$selected[1], ]
  expect_equal(row1$fold_gain, mean(d1), tolerance = 1e-12)
  expect_equal(row1$fold_gain_se, stats::sd(d1) / 2, tolerance = 1e-12)
  # Every accepted step cleared both parts of the rule ...
  for (s in seq_along(sel$selected)) {
    r <- h[h$step == s & h$variable == sel$selected[s], ]
    expect_gt(r$gain, 0)
    expect_gt(r$fold_gain, r$fold_gain_se)
  }
  # ... and the step it stopped at did not.
  nxt <- h[h$step == length(sel$selected) + 1L, ]
  if (nrow(nxt)) {
    best <- nxt[which.min(nxt$score), ]
    expect_false(isTRUE(best$gain > 0 && best$fold_gain > best$fold_gain_se))
  }
})

test_that("the default rule refuses gains the folds do not agree on", {
  # One real predictor, one weaker one, four noise candidates; ten draws.
  count <- function(rule) {
    out <- lapply(1:10, function(s) {
      set.seed(3000 + s)
      n <- 160
      d <- data.frame(x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000),
                      a = stats::rnorm(n), b = stats::rnorm(n))
      for (j in 1:4) d[[paste0("n", j)]] <- stats::rnorm(n)
      d$z <- d$a + 0.4 * d$b + stats::rnorm(n, 0, 0.7)
      d <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
      sel_quiet(select_features_forward(d, "z", c("a", "b", paste0("n", 1:4)),
                                        sel_lm, k = 5, seed = s, quiet = TRUE,
                                        rule = rule))$selected
    })
    c(noise = sum(vapply(out, function(v) sum(grepl("^n", v)), numeric(1))),
      a = sum(vapply(out, function(v) "a" %in% v, logical(1))),
      b = sum(vapply(out, function(v) "b" %in% v, logical(1))))
  }
  se  <- count("se")
  tol <- count("tol")
  # The strong predictor is never lost, and the weak one mostly kept: the
  # rule is not simply refusing everything.  (Here: 10 and 10 of 10.)
  expect_identical(unname(se[["a"]]), 10)
  expect_gte(se[["b"]], 8)
  # Fewer noise variables than when any gain is accepted.  (Here: 1 against
  # 7.  Over 40 draws with a weaker second predictor, coefficient 0.3, the
  # counts were 7 against 36, with that predictor kept in 34 and 40.)
  expect_lt(se[["noise"]], tol[["noise"]])
  expect_lte(se[["noise"]], 3)
  expect_gte(tol[["noise"]], 5)
})

test_that("rule = \"tol\" accepts any gain above tol, as before", {
  d <- sel_signal(n = 200, seed = 6, b = 1)
  cands <- c("a1", "a2", "n1", "n2", "n3")
  sel <- sel_quiet(select_features_forward(d, "z", cands, sel_lm, k = 4,
                                           seed = 6, quiet = TRUE, rule = "tol"))
  h <- sel$history
  expect_identical(sel$params$rule, "tol")
  for (s in seq_along(sel$selected)) {
    r <- h[h$step == s & h$variable == sel$selected[s], ]
    expect_gt(r$gain, 0)
  }
  nxt <- h[h$step == length(sel$selected) + 1L, ]
  if (nrow(nxt)) expect_false(isTRUE(max(nxt$gain, na.rm = TRUE) > 0))
  # A fold-mean baseline is applied under this rule too: an impossible tol
  # leaves nothing selected, where the first variable used to be accepted
  # whatever it scored.
  none <- sel_quiet(select_features_forward(d, "z", cands, sel_lm_refusing, k = 4,
                                            seed = 6, quiet = TRUE, rule = "tol",
                                            tol = 1e6))
  expect_identical(none$selected, character(0))
  expect_false(none$baseline$kept)
  expect_length(none$baseline$path, 1L)
  expect_true(none$baseline$path %in% c("a1", "a2"))
})

test_that("tol accepts Inf and -Inf, as it did before it was checked", {
  d <- sel_signal(n = 160, seed = 6, b = 1)
  cands <- c("a1", "a2", "n1")
  # Inf: nothing can clear it, so nothing is selected, under either rule and
  # with either kind of baseline.  (It stopped with "`tol` must be a single
  # finite number; got Inf".)
  for (rule in c("se", "tol")) {
    own <- sel_quiet(select_features_forward(d, "z", cands, sel_lm, k = 4,
                                             seed = 6, quiet = TRUE, tol = Inf,
                                             rule = rule))
    expect_identical(own$selected, character(0))
    expect_identical(own$baseline$type, "null model")
    expect_identical(own$params$tol, Inf)
    ref <- sel_quiet(select_features_forward(d, "z", cands, sel_lm_refusing, k = 4,
                                             seed = 6, quiet = TRUE, tol = Inf,
                                             rule = rule))
    expect_identical(ref$selected, character(0))
    expect_false(ref$baseline$kept)
  }
  # -Inf puts no floor under a step: with rule = "tol" every candidate is
  # added, whatever it scores.
  all_in <- sel_quiet(select_features_forward(d, "z", cands, sel_lm, k = 4,
                                              seed = 6, quiet = TRUE, tol = -Inf,
                                              rule = "tol"))
  expect_setequal(all_in$selected, cands)
  # What is not a single number is still refused, by name.
  expect_error(select_features_forward(d, "z", cands, sel_lm, tol = NA_real_),
               "`tol` must be a single number")
  expect_error(select_features_forward(d, "z", cands, sel_lm, tol = c(0, 1)),
               "`tol` must be a single number.*length 2")
  expect_error(select_features_forward(d, "z", cands, sel_lm, tol = "0"),
               "`tol` must be a single number.*character")
})

test_that("a GWR sweep on the fold-mean baseline is told what that baseline shows", {
  skip_if_not_installed("sp")
  d <- sel_signal(n = 120, seed = 7, b = 1)
  # Refusing the empty set here, as fit_gwr_model() does, so that the sweep
  # falls back to the fold mean whatever fit_gwr_model() does with it.
  gwr_fn <- function(tr, v) {
    if (!length(v)) stop("a GWR needs at least one predictor", call. = FALSE)
    fit_gwr_model(tr, "z", v, bandwidth = 30)
  }
  lines <- capture_spatialkit_log(
    sel <- sel_quiet(select_features_forward(d, "z", c("a1", "n1"), gwr_fn,
                                             k = 3, seed = 7, quiet = TRUE)))
  expect_identical(sel$baseline$type, "fold mean")
  expect_identical(sel$history$variable[1], "<none>")
  expect_true("a1" %in% sel$selected)
  expect_true(log_has(lines, "the learner is a GWR and the baseline is the training-fold mean"))
  expect_true(log_has(lines, "null_fn"))
  # It says where a working `null_fn` is.
  expect_true(log_has(lines, "The example in \\?select_features_forward builds one"))
  # Not said for a learner whose models are not GWRs.
  lines_lm <- capture_spatialkit_log(
    sel_quiet(select_features_forward(d, "z", c("a1", "n1"), sel_lm_refusing,
                                      k = 3, seed = 7, quiet = TRUE)))
  expect_false(log_has(lines_lm, "the learner is a GWR"))
})


# --- (c) the first variable under the fold-mean baseline ----------------------

# z = a1 + a2 + noise, and `proxy` = a1 + a2 + more noise: the best single
# predictor, and of no use once a1 and a2 are both in.
sel_proxy <- function(n = 240, seed = 1) {
  set.seed(4000 + seed)
  d <- data.frame(x = stats::runif(n, 0, 2000), y = stats::runif(n, 0, 2000),
                  a1 = stats::rnorm(n), a2 = stats::rnorm(n), n1 = stats::rnorm(n))
  d$proxy <- d$a1 + d$a2 + stats::rnorm(n)
  d$z     <- d$a1 + d$a2 + stats::rnorm(n, 0, 0.5)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

test_that("a first variable the finished set does not need is dropped", {
  d <- sel_proxy()
  cands <- c("proxy", "a1", "a2", "n1")
  msgs <- character(0)
  sel <- withCallingHandlers(
    suppressWarnings(select_features_forward(d, "z", cands, sel_lm_refusing,
                                             k = 4, seed = 1)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  expect_identical(sel$baseline$type, "fold mean")
  # The sweep walked proxy, then the two real predictors ...
  expect_identical(sel$baseline$path[1], "proxy")
  expect_true(all(c("a1", "a2") %in% sel$baseline$path))
  expect_true(sel$baseline$kept)
  # ... and the selection is that path without its first variable.  (It was
  # the whole path: `proxy` stayed, never having been tested.)
  f <- sel$baseline$first
  expect_true(is.list(f))
  expect_identical(f$variable, "proxy")
  expect_true(f$dropped)
  expect_identical(sel$selected, sel$baseline$path[-1])
  expect_false("proxy" %in% sel$selected)
  expect_match(f$note, "^dropped")
  # The numbers the rule read: the path against the path without `proxy`,
  # which did not clear it, and the rest against the fold mean, which did.
  expect_identical(f$rest$vars, sel$selected)
  expect_identical(sel$score, f$rest$score)
  last <- sel$history[sel$history$step == length(sel$baseline$path) &
                        sel$history$variable == utils::tail(sel$baseline$path, 1), ]
  expect_equal(f$gain, f$rest$score - last$score, tolerance = 1e-12)
  expect_false(isTRUE(f$gain > 0 && f$fold_gain > f$fold_gain_se))
  expect_equal(f$rest$gain, sel$baseline$score - f$rest$score, tolerance = 1e-12)
  expect_gt(f$rest$fold_gain, f$rest$fold_gain_se)
  # By hand: the selected set's RMSE on the sweep's own folds.
  fo <- sel_quiet(make_folds(d, k = 4, method = "block_kfold", seed = 1,
                             response_var = "z", predictor_vars = cands))
  e <- unlist(lapply(fo$folds, function(sp) {
    m <- stats::lm(stats::reformulate(sel$selected, "z"),
                   sf::st_drop_geometry(d)[sp$train, ])
    d$z[sp$test] - stats::predict(m, sf::st_drop_geometry(d)[sp$test, ])
  }))
  expect_equal(sel$score, sqrt(mean(e^2)), tolerance = 1e-10)
  # The sweep's record is the sweep as it ran: the extra set is not a step.
  expect_identical(max(sel$history$step), length(sel$baseline$path) +
                     as.integer(length(sel$baseline$path) < length(cands)))
  expect_true(any(grepl("proxy was added first, without a test.*so it is dropped",
                        msgs)))
  expect_true(any(grepl("its first variable only if the set is better with it",
                        msgs)))

  # A learner with a null model of its own tested `proxy` at step 1, where
  # it does beat the intercept-only model, and has no such check to make:
  # it walks the same path and keeps all of it.
  own <- sel_quiet(select_features_forward(d, "z", cands, sel_lm, k = 4,
                                           seed = 1, quiet = TRUE))
  expect_null(own$baseline$first)
  expect_identical(own$selected, sel$baseline$path)
})

test_that("a first variable the set needs is kept, and a path of one has no second test", {
  d <- sel_signal(n = 200, seed = 2, b = 1)
  msgs <- character(0)
  sel <- withCallingHandlers(
    suppressWarnings(select_features_forward(d, "z", c("a1", "a2", "n1"),
                                             sel_lm_refusing, k = 4, seed = 2)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  f <- sel$baseline$first
  expect_true(is.list(f))
  expect_false(f$dropped)
  expect_match(f$note, "^kept: the set is better with it")
  expect_identical(sel$selected, sel$baseline$path)
  expect_identical(f$variable, sel$selected[1])
  expect_gt(f$gain, 0)
  expect_gt(f$fold_gain, f$fold_gain_se)
  expect_true(any(grepl("added first without a test, is kept", msgs)))

  # One real predictor: the path has one variable, and the test against the
  # fold mean was the test of that variable.
  d1 <- d; d1$z <- d1$a1 + stats::rnorm(nrow(d1), 0, 0.5)
  one <- sel_quiet(select_features_forward(d1, "z", c("a1", "n1"),
                                           sel_lm_refusing, k = 4, seed = 2,
                                           quiet = TRUE, max_vars = 1))
  expect_identical(one$selected, "a1")
  expect_true(one$baseline$kept)
  expect_null(one$baseline$first)
  # Nothing kept, nothing to test.
  none <- sel_quiet(select_features_forward(sel_noise(seed = 2), "z",
                                            paste0("c", 1:4), sel_lm_refusing,
                                            k = 4, seed = 2, quiet = TRUE,
                                            tol = Inf))
  expect_identical(none$selected, character(0))
  expect_null(none$baseline$first)
})

test_that("a first variable that cannot be tested stays, and the warning says so", {
  d <- sel_proxy(seed = 2)
  # Fits anything but the one set the end-of-path test needs.
  picky <- function(tr, v) {
    if (!length(v) || setequal(v, c("a1", "a2")))
      stop("this backend cannot fit that set", call. = FALSE)
    lm_spatial_fit(tr, "z", v)
  }
  w <- character(0)
  sel <- withCallingHandlers(
    suppressMessages(select_features_forward(d, "z", c("proxy", "a1", "a2"),
                                             picky, k = 4, seed = 2,
                                             quiet = TRUE)),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    })
  expect_identical(sel$baseline$path, sel$selected)
  expect_identical(sel$selected[1], "proxy")
  f <- sel$baseline$first
  expect_false(f$dropped)
  expect_match(f$note, "^not tested")
  expect_true(is.na(f$rest$score))
  expect_true(any(grepl(paste0("proxy was added first, without a test, and the ",
                               "set without it, \\{a1, a2\\}, could not be scored"),
                        w)))
})

test_that("a first variable that cannot be told apart from the rest is kept and called unproven", {
  # A learner whose error in each fold is read off a table, so that the three
  # comparisons come out as wanted: the path {v1, v2} beats the fold mean in
  # every fold; the path without v1 is no worse by the rule (better in two
  # folds, worse in two); and {v2} alone does not beat the fold mean by the
  # rule.  Nothing separates v1 from v2 then, and dropping v1 would leave a
  # set the rule does not support either.
  set.seed(42)
  n <- 200
  d <- sf::st_as_sf(
    data.frame(x = stats::runif(n, 0, 2000), y = stats::runif(n, 0, 2000),
               v1 = stats::rnorm(n), v2 = stats::rnorm(n), v3 = stats::rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  fo <- sel_quiet(make_folds(d, k = 4, method = "block_kfold", seed = 1))
  fold <- integer(n)
  for (j in seq_along(fo$folds)) fold[fo$folds[[j]]$test] <- j
  # z alternates +1 / -1 inside each fold: every training-fold mean is about
  # 0, and the fold mean scores an RMSE of 1 in every fold.
  d$z <- NA_real_
  for (j in 1:4) {
    i <- which(fold == j)
    d$z[i] <- rep(c(1, -1), length.out = length(i))
  }
  key <- function(x) {
    xy <- sf::st_coordinates(x)
    paste(round(xy[, 1], 4), round(xy[, 2], 4))
  }
  look <- data.frame(key = key(d), z = d$z, fold = fold)
  err <- list("v1" = rep(0.95, 4), "v2" = c(1.25, 0.70, 1.25, 0.70),
              "v3" = rep(1.30, 4), "v1+v2" = c(0.90, 0.91, 0.90, 0.91),
              "v1+v3" = rep(0.97, 4), "v1+v2+v3" = rep(0.96, 4))
  oracle <- function(tr, v) {
    if (!length(v)) stop("this backend needs at least one predictor", call. = FALSE)
    new_spatial_fit("sel_oracle_fit",
                    engine = list(set = paste(sort(v), collapse = "+")),
                    formula = z ~ 1, response_var = "z", predictor_vars = v,
                    data_sf = tr)
  }
  registerS3method("predict", "sel_oracle_fit", function(object, newdata = NULL, ...) {
    r <- look[match(key(newdata), look$key), ]
    r$z + err[[object$engine$set]][r$fold]
  })
  lines <- capture_spatialkit_log(
    sel <- sel_quiet(select_features_forward(d, "z", c("v1", "v2", "v3"), oracle,
                                             k = 4, seed = 1, quiet = TRUE)))
  # The sweep's folds are the ones the table was laid out on.
  expect_equal(sel$baseline$score, 1, tolerance = 1e-3)
  expect_identical(sel$baseline$path, c("v1", "v2"))
  expect_true(sel$baseline$kept)
  f <- sel$baseline$first
  expect_identical(f$variable, "v1")
  expect_false(f$dropped)
  expect_identical(sel$selected, c("v1", "v2"))
  expect_match(f$note, "^kept, unproven")
  # v1's gain over {v2}: within its standard error.  {v2} over the fold
  # mean: negative on the pooled score.
  expect_lt(f$fold_gain, f$fold_gain_se)
  expect_lt(f$rest$gain, 0)
  expect_equal(f$rest$score, sqrt(sum(table(fold) * err$v2^2) / n), tolerance = 1e-10)
  expect_true(log_has(lines, "v1 was added first, without a test, and this selection does not show that it predicts"))
  expect_true(log_has(lines, "Read v1 as unproven"))
})


# --- (d) a `null_fn` for a GWR, from exported functions ------------------------

# The model the example of ?select_features_forward builds: the
# kernel-weighted local mean under fit_gwr_model()'s default kernel
# (adaptive bisquare), which is an intercept-only GWR.
sel_local_mean <- function(tr, bandwidth = 30) {
  new_spatial_fit("sel_local_mean_fit",
                  engine = list(xy = sf::st_coordinates(tr), y = tr$z,
                                bandwidth = bandwidth),
                  formula = z ~ 1, response_var = "z",
                  predictor_vars = character(0), data_sf = tr)
}
predict.sel_local_mean_fit <- function(object, newdata = NULL, ...) {
  e  <- object$engine
  at <- if (is.null(newdata)) e$xy else sf::st_coordinates(newdata)
  vapply(seq_len(nrow(at)), function(i) {
    d <- sqrt((e$xy[, 1] - at[i, 1])^2 + (e$xy[, 2] - at[i, 2])^2)
    h <- sort(d)[min(e$bandwidth, length(d))]
    w <- ifelse(d > h, 0, (1 - d^2 / h^2)^2)
    sum(w * e$y) / sum(w)
  }, numeric(1))
}
registerS3method("predict", "sel_local_mean_fit", predict.sel_local_mean_fit)

test_that("the example's local mean is the package's own kernel-weighted mean", {
  d  <- sel_signal(n = 150, seed = 8, b = 1)
  tr <- d[1:100, ]; te <- d[101:150, ]
  got <- predict(sel_local_mean(tr, 30), te)
  xy  <- sf::st_coordinates(tr); at <- sf::st_coordinates(te)
  ref <- vapply(seq_len(nrow(at)), function(i) {
    dd <- sqrt((xy[, 1] - at[i, 1])^2 + (xy[, 2] - at[i, 2])^2)
    w  <- spatialkit:::.gw_kernel_weights(dd, 30, "bisquare", adaptive = TRUE)
    sum(w * tr$z) / sum(w)
  }, numeric(1))
  expect_equal(got, ref, tolerance = 1e-12)
  expect_true(all(is.finite(got)))
})

test_that("with the local mean as `null_fn` a GWR's first variable is tested at step 1", {
  skip_if_not_installed("sp")
  d <- sel_signal(n = 120, seed = 7, b = 1)
  gwr_fn <- function(tr, v) {
    if (!length(v)) stop("a GWR needs at least one predictor", call. = FALSE)
    fit_gwr_model(tr, "z", v, bandwidth = 30)
  }
  lines <- capture_spatialkit_log(
    sel <- sel_quiet(select_features_forward(d, "z", c("a1", "n1"), gwr_fn,
                                             k = 3, seed = 7, quiet = TRUE,
                                             null_fn = sel_local_mean)))
  expect_identical(sel$baseline$type, "null_fn")
  expect_null(sel$baseline$first)
  expect_true("a1" %in% sel$selected)
  # Step 0 is the local mean on the sweep's folds, and step 1 was held to it.
  h <- sel$history
  expect_identical(h$variable[1], "<none>")
  r1 <- h[h$step == 1L & h$variable == sel$selected[1], ]
  expect_gt(r1$gain, 0)
  expect_gt(r1$fold_gain, r1$fold_gain_se)
  # No caution: the baseline is the one it asks for.
  expect_false(log_has(lines, "the learner is a GWR and the baseline is the training-fold mean"))
  # Step 0 is the local mean scored on the sweep's own folds.
  fo <- sel_quiet(make_folds(d, k = 3, method = "block_kfold", seed = 7,
                             response_var = "z", predictor_vars = c("a1", "n1")))
  e <- unlist(lapply(fo$folds, function(sp)
    d$z[sp$test] - predict(sel_local_mean(d[sp$train, ]), d[sp$test, ])))
  expect_equal(sel$baseline$score, sqrt(mean(e^2)), tolerance = 1e-8)
})
