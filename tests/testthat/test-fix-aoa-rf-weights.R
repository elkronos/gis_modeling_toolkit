# tests/testthat/test-fix-aoa-rf-weights.R
# ---------------------------------------------------------------------------
# area_of_applicability(weights = "auto").
#
# The default weighted every predictor equally even when the model carried
# variable importances, although the method the function names weights by
# importance.  With 2 predictors a forest uses and 8 it does not, equal
# weights put about 70% of rows shifted 3 SD along a used predictor (4 times
# the in-sample RMSE) inside the AOA and 43% of rows shifted 4 SD along an
# unused one; pmax(importance, 0) puts 26% and 93% inside.
#
# The default is now "auto": the model's importance when it is finite, has a
# positive value, names every predictor in use and comes from a model with
# out-of-bag skill; equal weights, with a logged note, otherwise.
# weights = NULL stays the explicit request for equal weights.
# ---------------------------------------------------------------------------

wt_pv <- c("a", "b", paste0("u", 1:8))
wt_pts <- function(n, shift_real = 0, shift_noise = 0) {
  d <- data.frame(x = stats::runif(n, 0, 1e4), y = stats::runif(n, 0, 1e4),
                  a = stats::rnorm(n) + shift_real, b = stats::rnorm(n))
  for (j in 1:8) d[[paste0("u", j)]] <- stats::rnorm(n) + shift_noise * (j == 1)
  d$z <- 2 * d$a + sin(2 * d$b) * 2 + stats::rnorm(n, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}
wt_fit <- function(seed = 11, ...) {
  set.seed(seed)
  fit_rf_model(wt_pts(600), "z", wt_pv, num_trees = 300, seed = seed, ...)
}


test_that("the default weights a forest's predictors by its importance", {
  skip_if_not_installed("ranger")
  fit <- wt_fit()
  nd  <- rbind(wt_pts(300), wt_pts(300, shift_real = 3),
               wt_pts(300, shift_noise = 4))
  grp <- rep(c("same", "used", "unused"), each = 300)

  auto <- area_of_applicability(nd, model = fit)
  hand <- area_of_applicability(nd, model = fit,
                                weights = pmax(fit$info$importance, 0))
  none <- area_of_applicability(nd, model = fit, weights = NULL)

  # The default is the importance-weighted index, number for number.
  expect_identical(auto$aoa$DI, hand$aoa$DI)
  expect_identical(auto$aoa$AOA, hand$aoa$AOA)
  expect_identical(auto$threshold, hand$threshold)
  expect_identical(auto$weights, hand$weights)
  expect_identical(auto$params$weights_source, "model importance")
  expect_identical(auto$params$weights_note, "permutation importance")
  expect_false(auto$params$weights_supplied)
  expect_identical(hand$params$weights_source, "supplied")
  expect_true(hand$params$weights_supplied)
  # weights = NULL is still equal weights.
  expect_identical(unname(none$weights), rep(1, 10))
  expect_identical(none$params$weights_source, "equal")
  expect_false(isTRUE(all.equal(auto$aoa$DI, none$aoa$DI)))

  # What the weights are for.  A 3 SD shift along a predictor the forest
  # uses is mostly outside; a 4 SD shift along one it ignores is mostly in.
  expect_lt(mean(auto$aoa$AOA[grp == "used"]), 0.40)
  expect_gt(mean(auto$aoa$AOA[grp == "unused"]), 0.85)
  # With equal weights it was the other way round.
  expect_gt(mean(none$aoa$AOA[grp == "used"]), 0.55)
  expect_lt(mean(none$aoa$AOA[grp == "unused"]), 0.65)

  # print() says where the weights came from.
  p_auto <- utils::capture.output(print(auto))
  p_hand <- utils::capture.output(print(hand))
  p_none <- utils::capture.output(print(none))
  expect_true(any(grepl("weighted    : yes, by the model's permutation importance",
                        p_auto, fixed = TRUE)))
  expect_true(any(grepl("weighted    : yes, by the weights supplied", p_hand,
                        fixed = TRUE)))
  expect_true(any(grepl("weighted    : no (all predictors count equally: weights = NULL)",
                        p_none, fixed = TRUE)))
})

test_that("\"auto\" is equal weights, without an error, when the importance means nothing", {
  skip_if_not_installed("ranger")
  set.seed(21)
  trn <- wt_pts(200); nd <- wt_pts(60)
  equal_and_says <- function(model, pattern, ...) {
    lines <- capture_spatialkit_log(
      res <- area_of_applicability(nd, model = model, ...), level = logger::INFO)
    expect_identical(unname(res$weights), rep(1, length(res$weights)))
    expect_identical(res$params$weights_source, "equal")
    expect_match(res$params$weights_note, pattern)
    expect_true(log_has(lines, paste0("weights = \"auto\": .*", pattern)))
    out <- utils::capture.output(print(res))
    expect_true(any(grepl("no (all predictors count equally: ", out, fixed = TRUE)))
    invisible(res)
  }

  # importance = "none": nothing to read.
  f_none <- fit_rf_model(trn, "z", wt_pv, num_trees = 60, importance = "none")
  equal_and_says(f_none, "carries no variable importance")

  # No out-of-bag rows: the importance is NaN, which used to be an error
  # for anyone who passed it.
  f_nan <- suppressWarnings(
    fit_rf_model(trn, "z", wt_pv, num_trees = 60, replace = FALSE,
                 sample_fraction = 1))
  expect_true(all(is.nan(f_nan$info$importance)))
  equal_and_says(f_nan, "not finite")
  expect_error(area_of_applicability(nd, model = f_nan,
                                     weights = pmax(f_nan$info$importance, 0)),
               "NaN")

  # A forest with no out-of-bag skill: its positive importances are noise.
  noise <- trn; noise$z <- stats::rnorm(nrow(noise))
  f_noise <- fit_rf_model(noise, "z", wt_pv, num_trees = 100)
  f_noise$info$oob_r_squared <- min(f_noise$info$oob_r_squared, -0.02)
  expect_true(any(f_noise$info$importance > 0))
  equal_and_says(f_noise, "no out-of-bag skill")

  # predictor_vars overridden with a column the model never saw.
  f_ab <- fit_rf_model(trn, "z", c("a", "b"), num_trees = 60)
  equal_and_says(f_ab, "no entry for .u1.", predictor_vars = c("a", "b", "u1"))

  # A model with no importance at all (a GWR, a Bayesian fit, a learner of
  # one's own).
  bare <- structure(list(data_sf = trn, predictor_vars = c("a", "b"),
                         info = list()), class = c("lmsurf_fit", "spatial_fit"))
  equal_and_says(bare, "carries no variable importance")

  # No model: equal weights, and nothing to log about it.
  lines <- capture_spatialkit_log(
    res <- area_of_applicability(nd, train_sf = trn, predictor_vars = c("a", "b")),
    level = logger::INFO)
  expect_identical(unname(res$weights), c(1, 1))
  expect_identical(res$params$weights_source, "equal")
  expect_false(log_has(lines, "weights = \"auto\""))
})

test_that("\"auto\" covers the coordinate columns of a coordinate-using forest", {
  skip_if_not_installed("ranger")
  set.seed(31)
  trn <- wt_pts(250); nd <- wt_pts(60)
  f_xy <- suppressWarnings(suppressMessages(
    fit_rf_model(trn, "z", c("a", "b"), num_trees = 100, include_coords = TRUE)))
  auto <- area_of_applicability(nd, model = f_xy)
  hand <- area_of_applicability(nd, model = f_xy,
                                weights = pmax(f_xy$info$importance, 0))
  expect_identical(names(auto$weights), c("a", "b", "..x", "..y"))
  expect_identical(auto$weights, hand$weights)
  expect_identical(auto$aoa$DI, hand$aoa$DI)
})

test_that("a weights value that is neither \"auto\", NULL nor numeric is refused", {
  trn <- wt_pts(60); nd <- wt_pts(20)
  expect_error(area_of_applicability(nd, train_sf = trn, predictor_vars = c("a", "b"),
                                     weights = "importance"),
               "`weights` must be numeric, or NULL for equal weights, or \"auto\"")
})
