# tests/testthat/test-fix-exact-fit-rounding.R
# ---------------------------------------------------------------------------
# The residuals of a response the predictors fit exactly are rounding, and
# neither the test for `cp` in resolution_profile() nor estimate_sac_range()
# scores them.  Both read .residuals_rounding(), whose bound on the
# regression's rounding is taken on the size of the fit's terms and is read
# on the residuals taken again from the coefficients, so it holds whatever
# the number of predictors and rows.  Residuals above it are scored as
# before (BUG-61).
# ---------------------------------------------------------------------------

.ef_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# BUG-36's layer: 800 points, an exponential field (range 600 m) plus 0.5 w
# plus noise, a second predictor w2, and here a third, w3.
.ef_layer <- function() {
  .ef_seed(5)
  n <- 800
  xy <- data.frame(x = stats::runif(n, 0, 4000), y = stats::runif(n, 0, 4000))
  D <- as.matrix(stats::dist(xy))
  f <- as.numeric(t(chol(exp(-D / 600) + diag(1e-8, n))) %*% stats::rnorm(n))
  xy$w  <- stats::rnorm(n)
  xy$z  <- f + 0.5 * xy$w + stats::rnorm(n, sd = 0.3)
  xy$w2 <- stats::rnorm(n)
  xy$f  <- f
  .ef_seed(61)
  xy$w3 <- stats::rnorm(n)
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32615)
}

.ef_catch <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w, classes = vapply(w, function(c) class(c)[1L], ""))
}

.ef_exact_reason <- paste0("the residuals on predictor_vars are constant: the predictors ",
                           "explain the response exactly")

test_that("the exact fit 7 w - 2 w2 is constant to the test for cp and to the estimate", {
  L <- .ef_layer()
  d <- L; d$v <- 7 * L$w - 2 * L$w2
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = c("w", "w2")))
  expect_true(is.na(s))
  expect_identical(attr(s, "rejected_reason"), .ef_exact_reason)
  r <- .ef_catch(resolution_profile(d, "v", predictor_vars = c("w", "w2"), n_levels = 6,
                                    nstart = 2, seed = 1))
  expect_identical(attr(r$value, "variable"), "residuals")
  expect_true(attr(r$value, "one_cell")$constant)
  expect_true("spatialkit_constant_response" %in% r$classes)
  msg <- conditionMessage(r$warnings[[match("spatialkit_constant_response", r$classes)]])
  expect_match(msg, "the OLS residual of the response on `predictor_vars` is constant",
               fixed = TRUE)
  expect_true(all(is.na(r$value$cp)))
  expect_true(all(is.na(r$value$moran_i)))
  expect_true(all(is.na(r$value$moran_z)))
  expect_error(select_resolution(r$value, criterion = "cp"),
               "because the OLS residual of the response on the predictors is constant",
               fixed = TRUE)
})

test_that("exact fits on one to three predictors are constant in both functions", {
  L <- .ef_layer()
  preds <- list("w", c("w", "w2"), c("w", "w2", "w3"))
  .ef_seed(611)
  for (i in seq_len(100)) {
    pv <- preds[[(i - 1L) %% 3L + 1L]]
    cf <- round(stats::rnorm(length(pv) + 1L, sd = 5), 2)
    trend <- cf[1L] + as.numeric(as.matrix(sf::st_drop_geometry(L)[, pv, drop = FALSE]) %*%
                                   cf[-1L])
    for (off in c(0, 1e6)) {
      info <- sprintf("fit %d, %s, offset %g", i, paste(pv, collapse = " + "), off)
      d <- L; d$v <- off + trend
      s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = pv))
      expect_identical(attr(s, "rejected_reason"), .ef_exact_reason, info = info)
      r <- suppressWarnings(resolution_profile(d, "v", predictor_vars = pv, n_levels = 2,
                                               nstart = 1, seed = 1))
      expect_true(attr(r, "one_cell")$constant, info = info)
      expect_true(all(is.na(r$cp)), info = info)
    }
  }
})

test_that("residuals above the regression's rounding are scored as before", {
  skip_if_not_installed("gstat")
  L <- .ef_layer()
  # 1e3 w + 1e-11 f on w: residuals of 1e-11 of a field, about 45
  # .Machine$double.eps times the size of the fit's terms.
  d <- L; d$v <- 1e3 * L$w + 1e-11 * L$f
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "w"))
  expect_true(is.finite(s))
  expect_null(attr(s, "rejected_reason"))
  r <- .ef_catch(resolution_profile(d, "v", predictor_vars = "w", n_levels = 4, nstart = 2,
                                    seed = 1))
  expect_false(attr(r$value, "one_cell")$constant)
  expect_false("spatialkit_constant_response" %in% r$classes)
  expect_true(all(is.finite(r$value$cp)))
})

test_that(".residuals_rounding() holds the rounding of exact fits whatever the rows", {
  .ef_seed(17)
  # An integer predictor on 30,000 rows: lm.fit() leaves residuals of the
  # exact fit 2.5 x at hundreds of .Machine$double.eps times the size of the
  # terms, and the bound of .no_spread() on the response calls them spread.
  n <- 30000
  x <- sample(0:5, n, TRUE) + 0
  X <- cbind(1, x)
  y <- 2.5 * x
  fit <- stats::lm.fit(X, y)
  e <- fit$residuals
  expect_false(.no_spread(sum((e - mean(e))^2), y))
  expect_true(.residuals_rounding(e, y, y, X, fit$coefficients))
  # The same with a spread of 1e-11 of a normal column: spread.
  y2 <- y + 1e-11 * stats::rnorm(n)
  fit2 <- stats::lm.fit(X, y2)
  expect_false(.residuals_rounding(fit2$residuals, y2, y2, X, fit2$coefficients))
  # Exact fits on five predictors at an offset, raw and about the centre.
  P <- matrix(stats::rnorm(5 * 2000), 2000, 5)
  X5 <- cbind(1, P)
  y5 <- 1e6 + drop(P %*% c(3, -2, 0.5, 7, -1e3))
  for (c0 in c(0, .residual_centre(y5))) {
    f5 <- stats::lm.fit(X5, y5 - c0)
    expect_true(.residuals_rounding(f5$residuals, y5, y5 - c0, X5, f5$coefficients),
                info = format(c0))
  }
  # An aliased predictor: its coefficient is NA and counts as 0.
  Xa <- cbind(1, P[, 1], 2 * P[, 1])
  ya <- 4 + 3 * P[, 1]
  fa <- stats::lm.fit(Xa, ya)
  expect_true(anyNA(fa$coefficients))
  expect_true(.residuals_rounding(fa$residuals, ya, ya, Xa, fa$coefficients))
  # Residuals of a fit that leaves real spread are not rounding.
  yb <- ya + stats::rnorm(2000)
  fb <- stats::lm.fit(Xa, yb)
  expect_false(.residuals_rounding(fb$residuals, yb, yb, Xa, fb$coefficients))
})

test_that("Moran's I of the cell means reads the same bound on an exact fit", {
  # .morans_i_for_k() regresses the cell means again about the response's
  # centre where their residuals are small, and refuses those residuals
  # where .residuals_rounding() holds them.  An exact fit on an integer
  # predictor over 300 cells left them above 100 eps^2 times the cell
  # means' sum of squares, and I was computed on rounding.
  .ef_seed(5)
  n <- 6000
  xy <- cbind(stats::runif(n, 0, 4000), stats::runif(n, 0, 4000))
  P <- cbind(w = stats::rnorm(n), w2 = stats::rnorm(n), x = sample(0:5, n, TRUE) + 0)
  cl <- stats::kmeans(xy, 300, nstart = 1, iter.max = 50)$cluster
  cases <- list("-3.19 - 5.6 x" = list(cols = "x", y = -3.19 - 5.6 * P[, "x"]),
                "7 w - 2 w2"    = list(cols = c("w", "w2"), y = 7 * P[, "w"] - 2 * P[, "w2"]),
                "6.55 + 3.54 w - 1.09 x" = list(cols = c("w", "x"),
                                                 y = 6.55 + 3.54 * P[, "w"] - 1.09 * P[, "x"]))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    for (off in c(0, 1e6)) {
      m <- .morans_i_for_k(xy, off + cs$y, P[, cs$cols, drop = FALSE], cl, grouped = TRUE)
      expect_true(is.na(m[["I"]]), info = paste(nm, off))
      expect_true(is.na(m[["z"]]), info = paste(nm, off))
    }
  }
  # A spread on the same trend is scored.
  y <- 1e6 - 3.19 - 5.6 * P[, "x"] + 1e-3 * xy[, 1] / 4000 + stats::rnorm(n, sd = 1e-4)
  m <- .morans_i_for_k(xy, y, P[, "x", drop = FALSE], cl, grouped = TRUE)
  expect_true(is.finite(m[["I"]]))
})

test_that("exact fits on a predictor far from 0 are constant in both functions", {
  # A northing (4.5e6 plus the y coordinate) as a predictor: lm() leaves the
  # residuals of an exact fit at about 2 .Machine$double.eps times the size
  # of the fit's terms, thousands of times the response, so a floor of 1e-12
  # of the response's largest value finds spread in them.  estimate_sac_range()
  # read the shared bound only where that floor found none, and ran its test
  # for spatial structure on the rounding while the profile refused the layer.
  L <- .ef_layer()
  xy <- sf::st_coordinates(L)
  L$E <- 5e5 + xy[, 1]
  L$N <- 4.5e6 + xy[, 2]
  D <- sf::st_drop_geometry(L)
  .ef_seed(612)
  sets <- list("N", c("E", "N"), c("w", "N"))
  for (i in seq_len(15)) {
    pv <- sets[[(i - 1L) %% 3L + 1L]]
    cf <- round(stats::rnorm(length(pv) + 1L, sd = 5), 3)
    X <- as.matrix(D[, pv, drop = FALSE])
    # Near 0 at the centre of the predictors: the intercept takes off the offset.
    d <- L
    d$v <- cf[1L] + drop(X %*% cf[-1L]) - sum(colMeans(X) * cf[-1L])
    info <- sprintf("fit %d, %s", i, paste(pv, collapse = " + "))
    s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = pv))
    expect_identical(attr(s, "rejected_reason"), .ef_exact_reason, info = info)
    r <- suppressWarnings(resolution_profile(d, "v", predictor_vars = pv, n_levels = 2,
                                             nstart = 1, seed = 1))
    expect_true(attr(r, "one_cell")$constant, info = info)
    expect_true(all(is.na(r$cp)), info = info)
  }
})

test_that("the bound is taken on the absolute terms, not on their net sum", {
  # 1e-3 (N - 4.5e6) on N, plus 1e-12 of a field: the residuals' spread of
  # about 8e-13 lies below the rounding lm() leaves on the fit (4.3e-12) and
  # below 10 .Machine$double.eps times the absolute terms (2e-11), and far
  # above that bound on the terms' net sum (2.6e-15), which would score it.
  L <- .ef_layer()
  L$N <- 4.5e6 + sf::st_coordinates(L)[, 2]
  d <- L; d$v <- 1e-3 * (L$N - 4.5e6) + 1e-12 * L$f
  X <- cbind(1, d$N)
  fit <- stats::lm.fit(X, d$v)
  expect_true(.residuals_rounding(fit$residuals, d$v, d$v, X, fit$coefficients))
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "N"))
  expect_identical(attr(s, "rejected_reason"), .ef_exact_reason)
  r <- .ef_catch(resolution_profile(d, "v", predictor_vars = "N", n_levels = 3, nstart = 1,
                                    seed = 1))
  expect_true(attr(r$value, "one_cell")$constant)
  expect_true("spatialkit_constant_response" %in% r$classes)
  expect_true(all(is.na(r$value$cp)))
  # A spread of 1e-10 of the field, above the bound, is scored by both.
  skip_if_not_installed("gstat")
  d$v <- 1e-3 * (L$N - 4.5e6) + 1e-10 * L$f
  s <- suppressWarnings(estimate_sac_range(d, "v", predictor_vars = "N"))
  expect_true(is.finite(s))
  r <- suppressWarnings(resolution_profile(d, "v", predictor_vars = "N", n_levels = 3,
                                           nstart = 1, seed = 1))
  expect_false(attr(r, "one_cell")$constant)
  expect_true(all(is.finite(r$cp)))
})

test_that(".residuals_rounding() reads fits at the ends of the floating-point range", {
  # The sums are taken on the residuals and terms divided by a power of two,
  # so that at 1e-170 they do not underflow to 0 and at 1e200 they do not
  # overflow to Inf; either way a spread of 0.1 of a field read as constant.
  L <- .ef_layer()
  for (sc in c(1e-170, 1e200)) {
    X <- cbind(1, sc * L$w, L$w2)
    for (kind in c("exact", "spread")) {
      v <- sc * (7 * L$w - 2 * L$w2 + if (kind == "spread") 0.1 * L$f else 0)
      fit <- stats::lm.fit(X, v)
      expect_identical(.residuals_rounding(fit$residuals, v, v, X, fit$coefficients),
                       kind == "exact", info = paste(format(sc), kind))
    }
    d <- L; d$w <- sc * L$w
    d$v <- sc * (7 * L$w - 2 * L$w2 + 0.1 * L$f)
    r <- suppressWarnings(resolution_profile(d, "v", predictor_vars = c("w", "w2"),
                                             n_levels = 3, nstart = 1, seed = 1))
    expect_false(attr(r, "one_cell")$constant, info = format(sc))
  }
})
