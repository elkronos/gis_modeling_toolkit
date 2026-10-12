# tests/testthat/test-fix-aoa-rf-factors.R
# ---------------------------------------------------------------------------
# Factor predictors, in the forest and in its area of applicability.
#
# (a) fit_rf_model() called ranger() without respect.unordered.factors,
#     whose default for a regression forest is "ignore": a nominal factor is
#     split on the order its levels are stored in.  With a 40-level factor
#     carrying random level effects the test RMSE was 1.05 that way against
#     0.57 with "order", the handling ranger's help recommends.  "order" is
#     now what an unordered factor or character predictor gets unless the
#     caller passes the argument through `...`.
#
# (b) area_of_applicability() refused a forest fitted with a factor
#     predictor outright.  It still refuses by default, naming two ways
#     forward: categorical = "indicator" and categorical = "drop".
# ---------------------------------------------------------------------------

fac_pts <- function(n = 1500, L = 40, seed = 1) {
  set.seed(seed)
  g   <- factor(sample(sprintf("c%02d", seq_len(L)), n, TRUE))
  eff <- stats::rnorm(L, 0, 2); names(eff) <- levels(g)
  a   <- stats::rnorm(n)
  d <- sf::st_as_sf(
    data.frame(x = stats::runif(n, 0, 1e4), y = stats::runif(n, 0, 1e4), a, g),
    coords = c("x", "y"), crs = 32632)
  d$z <- unname(eff[as.character(g)]) + a + stats::rnorm(n, 0, 0.5)
  d
}
fac_rmse <- function(fit, te) sqrt(mean((predict(fit, te) - te$z)^2))


# --- (a) the forest ----------------------------------------------------------

test_that("a nominal factor is split by response order by default", {
  skip_if_not_installed("ranger")
  d  <- fac_pts()
  tr <- d[1:1000, ]; te <- d[1001:1500, ]
  f_def <- fit_rf_model(tr, "z", c("a", "g"), importance = "none", num_trees = 200)
  f_ord <- fit_rf_model(tr, "z", c("a", "g"), importance = "none", num_trees = 200,
                        respect.unordered.factors = "order")
  f_ign <- fit_rf_model(tr, "z", c("a", "g"), importance = "none", num_trees = 200,
                        respect.unordered.factors = "ignore")
  r_def <- fac_rmse(f_def, te); r_ord <- fac_rmse(f_ord, te); r_ign <- fac_rmse(f_ign, te)

  # The default is the "order" forest, not the "ignore" one it used to be.
  expect_equal(r_def, r_ord, tolerance = 0.05)
  expect_identical(predict(f_def, te), predict(f_ord, te))
  expect_lt(r_def, 0.75 * r_ign)

  expect_identical(f_def$info$respect_unordered_factors, "order")
  expect_identical(f_def$info$categorical_predictors, "g")
  # What the caller passes through `...` is what ranger gets, and is recorded.
  expect_identical(f_ign$info$respect_unordered_factors, "ignore")
  expect_false(identical(predict(f_ign, te), predict(f_def, te)))

  out <- utils::capture.output(print(f_def))
  expect_true(any(grepl("Categorical: g split by levels ordered by mean response",
                        out, fixed = TRUE)))
  expect_true(any(grepl("stored level order",
                        utils::capture.output(print(f_ign)), fixed = TRUE)))
})

test_that("predictions do not depend on the level order newdata carries", {
  skip_if_not_installed("ranger")
  d  <- fac_pts(n = 600, L = 12, seed = 2)
  tr <- d[1:400, ]; te <- d[401:600, ]
  fit <- fit_rf_model(tr, "z", c("a", "g"), importance = "none", num_trees = 100)
  p0 <- predict(fit, te)
  te_rev <- te
  te_rev$g <- factor(as.character(te$g), levels = rev(levels(te$g)))
  expect_identical(predict(fit, te_rev), p0)
  # Fewer levels than the forest was grown with, and a character column.
  te_few <- te[as.character(te$g) %in% levels(te$g)[1:5], ]
  keep <- which(as.character(te$g) %in% levels(te$g)[1:5])
  te_few$g <- droplevels(te_few$g)
  expect_identical(predict(fit, te_few), p0[keep])
  te_chr <- te; te_chr$g <- as.character(te$g)
  expect_identical(predict(fit, te_chr), p0)

  # A level the forest never saw is still an error naming it.
  te_new <- te[1:5, ]
  te_new$g <- factor(c("zz", as.character(te_new$g[-1])))
  expect_error(predict(fit, te_new), "level\\(s\\) the forest was never grown with.*zz")
})

test_that("a character predictor is treated as a nominal factor", {
  skip_if_not_installed("ranger")
  d  <- fac_pts(n = 900, L = 20, seed = 3)
  d$g <- as.character(d$g)
  tr <- d[1:600, ]; te <- d[601:900, ]
  f_def <- fit_rf_model(tr, "z", c("a", "g"), importance = "none", num_trees = 150)
  f_ign <- fit_rf_model(tr, "z", c("a", "g"), importance = "none", num_trees = 150,
                        respect.unordered.factors = "ignore")
  expect_identical(f_def$info$respect_unordered_factors, "order")
  expect_lt(fac_rmse(f_def, te), 0.85 * fac_rmse(f_ign, te))
})

test_that("ordered factors and numeric-only forests are left as they were", {
  skip_if_not_installed("ranger")
  d <- fac_pts(n = 300, L = 5, seed = 4)
  d$o <- factor(as.character(d$g), levels = levels(d$g), ordered = TRUE)
  f_o <- fit_rf_model(d, "z", c("a", "o"), importance = "none", num_trees = 50)
  expect_true(is.na(f_o$info$respect_unordered_factors))
  expect_identical(f_o$info$categorical_predictors, character(0))
  expect_false(any(grepl("Categorical", utils::capture.output(print(f_o)))))
  # ranger's own call, with nothing added: the same forest.
  r_o <- ranger::ranger(x = sf::st_drop_geometry(d)[, c("a", "o")], y = d$z,
                        num.trees = 50L, seed = 123L, importance = "none",
                        num.threads = 1L)
  expect_equal(f_o$engine$predictions, r_o$predictions)

  # splitrule = "extratrees" has a default of its own in ranger ("partition"),
  # which is left to it.
  f_x <- fit_rf_model(d, "z", c("a", "g"), importance = "none", num_trees = 50,
                      splitrule = "extratrees")
  expect_identical(f_x$info$respect_unordered_factors, "partition")
  r_x <- ranger::ranger(x = sf::st_drop_geometry(d)[, c("a", "g")], y = d$z,
                        num.trees = 50L, seed = 123L, importance = "none",
                        num.threads = 1L, splitrule = "extratrees")
  expect_equal(f_x$engine$predictions, r_x$predictions)

  f_n <- fit_rf_model(d, "z", "a", importance = "none", num_trees = 50)
  expect_true(is.na(f_n$info$respect_unordered_factors))
  r_n <- ranger::ranger(x = sf::st_drop_geometry(d)[, "a", drop = FALSE], y = d$z,
                        num.trees = 50L, seed = 123L, importance = "none",
                        num.threads = 1L)
  expect_equal(f_n$engine$predictions, r_n$predictions)
})


# --- (b) the area of applicability -------------------------------------------

fac_aoa_pts <- function(n = 200, seed = 5, levels = c("urban", "rural", "forest")) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = stats::runif(n, 0, 1e3), y = stats::runif(n, 0, 1e3),
               a = stats::rnorm(n), g = factor(sample(levels, n, TRUE))),
    coords = c("x", "y"), crs = 32632)
  d$z <- d$a + (d$g == "urban") + stats::rnorm(n, 0, 0.2)
  d
}

test_that("a forest with a factor predictor is refused with the two ways forward", {
  skip_if_not_installed("ranger")
  d <- fac_aoa_pts()
  fit <- fit_rf_model(d, "z", c("a", "g"), num_trees = 100)
  expect_error(area_of_applicability(d, model = fit),
               "not numeric or logical: .g.*categorical = \"indicator\".*categorical = \"drop\"")
})

test_that("categorical = \"indicator\" gives the forest an area of applicability", {
  skip_if_not_installed("ranger")
  d <- fac_aoa_pts()
  fit <- fit_rf_model(d, "z", c("a", "g"), num_trees = 100)
  nd  <- fac_aoa_pts(120, seed = 6)
  res <- area_of_applicability(nd, model = fit, categorical = "indicator")
  expect_s3_class(res, "aoa")
  expect_identical(res$predictor_vars, c("a", "g"))
  expect_identical(res$categorical_vars, "g")
  expect_identical(res$params$categorical, "indicator")
  expect_identical(names(res$weights), c("a", "g"))
  expect_identical(names(res$scaling$scale),
                   c("a", "g=forest", "g=rural", "g=urban"))
  # Rows on seen levels with ordinary numeric values are inside.
  expect_gt(mean(res$aoa$AOA), 0.85)
  expect_true(all(is.finite(res$aoa$DI)))

  # The scale the three indicator columns share: sqrt(1 - sum(p^2)).
  p <- as.numeric(table(d$g)) / nrow(d)
  expect_equal(unname(res$scaling$scale[-1]), rep(sqrt(1 - sum(p^2)), 3),
               tolerance = 1e-12)

  # A level no training row has is outside, with DI = Inf and a count.
  nd2 <- nd
  nd2$g <- factor(c(rep("desert", 7), as.character(nd$g[-(1:7)])))
  expect_warning(
    res2 <- area_of_applicability(nd2, model = fit, categorical = "indicator"),
    "7 of 120 prediction row\\(s\\) carry a level of .g. that no training row has")
  expect_true(all(is.infinite(res2$aoa$DI[1:7])))
  expect_false(any(res2$aoa$AOA[1:7]))
  expect_identical(res2$aoa$DI[-(1:7)], res$aoa$DI[-(1:7)])
  expect_true(any(grepl("7 outside on a level no training row has",
                        utils::capture.output(print(res2)), fixed = TRUE)))

  # A character column in newdata matches by value; a missing level is NA.
  nd3 <- nd; nd3$g <- as.character(nd$g); nd3$g[3] <- NA
  res3 <- area_of_applicability(nd3, model = fit, categorical = "indicator")
  expect_true(is.na(res3$aoa$DI[3]))
  expect_identical(res3$aoa$DI[-3], res$aoa$DI[-3])
})

test_that("the indicator distance is the same between any two levels", {
  # Three levels with very different shares and no other predictor: the
  # distance from a row to the nearest row on another level must not depend
  # on which levels they are.
  tr <- data.frame(g = rep(c("a", "b", "c"), times = c(60, 30, 10)),
                   stringsAsFactors = FALSE)
  nd <- data.frame(g = c("a", "b", "c"), stringsAsFactors = FALSE)
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = "g",
                               categorical = "indicator")
  expect_identical(res$aoa$DI, c(0, 0, 0))
  # Leave-one-out on the three distinct rows: each is sqrt(2 / (1 - sum(p^2)))
  # from the others, in units of the mean pairwise distance.
  p <- c(0.6, 0.3, 0.1)
  step <- sqrt(2 / (1 - sum(p^2)))
  expect_equal(unique(res$train_DI * res$normalizer), step, tolerance = 1e-12)
})

test_that("categorical = \"drop\" computes the index on the numeric predictors", {
  skip_if_not_installed("ranger")
  d <- fac_aoa_pts()
  fit <- fit_rf_model(d, "z", c("a", "g"), num_trees = 100)
  nd  <- fac_aoa_pts(80, seed = 7)
  lines <- capture_spatialkit_log(
    res <- area_of_applicability(nd, model = fit, categorical = "drop",
                                 weights = NULL))
  expect_true(log_has(lines, "categorical = \"drop\": the index is computed on a alone"))
  expect_identical(res$predictor_vars, "a")
  expect_identical(res$categorical_vars, "g")
  ref <- area_of_applicability(nd, train_sf = d, predictor_vars = "a")
  expect_identical(res$aoa$DI, ref$aoa$DI)
  expect_true(any(grepl("categorical : g left out of the index",
                        utils::capture.output(print(res)), fixed = TRUE)))
  # Nothing left to measure is an error, not an empty index.
  expect_error(area_of_applicability(nd, train_sf = d, predictor_vars = "g",
                                     categorical = "drop"),
               "leaves no predictor")
})

test_that("a model without categorical predictors is unaffected by `categorical`", {
  set.seed(8)
  tr <- data.frame(a = stats::rnorm(100), b = stats::rnorm(100))
  nd <- data.frame(a = stats::rnorm(40), b = stats::rnorm(40))
  ref <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  for (how in c("indicator", "drop")) {
    res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"),
                                 categorical = how)
    expect_identical(res$aoa$DI, ref$aoa$DI)
    expect_identical(res$threshold, ref$threshold)
    expect_identical(res$categorical_vars, character(0))
  }
  # Numeric when the index is defined, a factor in newdata: no setting can
  # reconcile them, and the error says that instead of offering one.
  nd_bad <- nd; nd_bad$b <- factor(round(nd$b))
  expect_error(area_of_applicability(nd_bad, train_sf = tr,
                                     predictor_vars = c("a", "b"),
                                     categorical = "indicator"),
               "numeric or logical in the training data but categorical in `newdata`")
})

test_that("dropped_vars is a character vector, by predictor or by level column", {
  # With nothing dropped it is character(0), as it always was.  The naming
  # of dropped indicator columns went through ifelse(), which returns
  # logical(0) for an empty test, so every result with no constant predictor
  # carried a logical there.
  set.seed(9)
  tr <- data.frame(a = stats::rnorm(60), b = stats::rnorm(60), k = 3)
  nd <- data.frame(a = stats::rnorm(20), b = stats::rnorm(20), k = 3)
  res <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "b"))
  expect_identical(res$dropped_vars, character(0))
  # A constant numeric predictor is named as itself.
  res_k <- area_of_applicability(nd, train_sf = tr,
                                 predictor_vars = c("a", "b", "k"))
  expect_identical(res_k$dropped_vars, "k")

  # Indicator coding.  A factor with one level among the training rows is
  # dropped whole and named as the predictor; a single level column of a
  # factor that is kept (its only row removed for a missing value elsewhere)
  # is named as the column.
  tr$g   <- rep(c("u", "f"), 30)
  tr$one <- "only"
  nd$g   <- rep(c("u", "f"), 10)
  nd$one <- "only"
  res_i <- area_of_applicability(nd, train_sf = tr,
                                 predictor_vars = c("a", "g", "one"),
                                 categorical = "indicator")
  expect_identical(res_i$dropped_vars, "one")
  expect_identical(res_i$predictor_vars, c("a", "g"))
  tr2 <- tr
  tr2$g[1] <- "r"; tr2$a[1] <- NA
  res_c <- area_of_applicability(nd, train_sf = tr2, predictor_vars = c("a", "g"),
                                 categorical = "indicator")
  expect_identical(res_c$dropped_vars, "g=r")
  expect_identical(res_c$predictor_vars, c("a", "g"))
  # And with nothing dropped under indicator coding, still character(0).
  res_n <- area_of_applicability(nd, train_sf = tr, predictor_vars = c("a", "g"),
                                 categorical = "indicator")
  expect_identical(res_n$dropped_vars, character(0))
})
