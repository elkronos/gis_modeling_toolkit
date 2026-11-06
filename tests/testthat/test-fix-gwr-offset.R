# fit_gwr_model() at a large offset of a predictor or of the response.
# The global least-squares benchmark was lm() on the raw columns, which drops
# a predictor at 1e8 plus a spread of 1 as collinear with the intercept; the
# local solve divided a predictor by its standard deviation before it
# centred it, so from about 1e10 the slopes, fitted values and AICc moved;
# and a response at 1e15 plus a spread of 1 was regressed on its rounding.
# A GWR's fitted values, residuals and slopes do not change under a shift,
# so only rounding may move.  The note on a predictor with near-zero
# variance read an absolute floor, so it depended on the predictor's unit.

skip_if_not_installed("sp")

# 300 points, a predictor w and a response with a spatial field; both as
# stored at 1e15, so that every offset used holds the same spread.
.go_layer <- function(n = 300, seed = 3) {
  withr::with_seed(seed, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4)
    w <- rnorm(n)
    z <- 2 * w + sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5)
    q <- rnorm(n)
  })
  list(x = x, y = y, wq = (1e15 + w) - 1e15, zq = (1e15 + z) - 1e15, q = q,
       mk = function(p, z, ...) sf::st_as_sf(
         data.frame(x = x, y = y, z = z, p = p, ...),
         coords = c("x", "y"), crs = 32615))
}

.go_fit <- function(d, ...) suppressWarnings(fit_gwr_model(d, "z", "p", ...))

.go_rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300), na.rm = TRUE)

test_that("a predictor at a large offset gives the benchmark, slopes, fitted values and AICc of the predictor at 0", {
  L <- .go_layer()
  ref <- .go_fit(L$mk(L$wq, L$zq), bandwidth = 60)
  ref_sum <- grep("least squares|global fit", capture.output(summary(ref)),
                  value = TRUE)
  expect_identical(ref$engine$lm$rank, 2L)
  for (o in c(1e8, 1e10, 1e12, 1e15)) {
    f <- .go_fit(L$mk(o + L$wq, L$zq), bandwidth = 60)
    expect_identical(f$engine$lm$rank, 2L)
    for (k in c("RSS", "RMSE", "AICc"))
      expect_lt(.go_rel(f$info$global_ols[[k]], ref$info$global_ols[[k]]), 1e-8)
    expect_identical(grep("least squares|global fit",
                          capture.output(summary(f)), value = TRUE), ref_sum)
    expect_lt(.go_rel(f$engine$SDF$p, ref$engine$SDF$p), 1e-8)
    expect_lt(max(abs(fitted(f) - fitted(ref))), 1e-8)
    expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-8)
    expect_lt(.go_rel(coef(f$engine$lm)[["p"]], coef(ref$engine$lm)[["p"]]), 1e-8)
    # The local and global intercepts are those at the predictor's 0, in
    # the data's units: intercept plus slope times the predictor gives the
    # fitted value, and the global fit's coefficients its fitted values.
    sdf <- f$engine$SDF@data
    expect_lt(max(abs(sdf$Intercept + sdf$p * (o + L$wq) - fitted(f))),
              1e-12 * o)
    expect_lt(max(abs(drop(cbind(1, o + L$wq) %*% coef(f$engine$lm)) -
                        fitted(f$engine$lm))), 1e-12 * o)
    expect_lt(max(abs(fitted(f$engine$lm) - fitted(ref$engine$lm))), 1e-8)
    expect_identical(unname(f$engine$lm$model$p), o + L$wq)
  }
})

test_that("predict() at new rows applies the training centres of a predictor at an offset", {
  L <- .go_layer()
  nd_xy <- withr::with_seed(9, cbind(runif(20, 1e3, 9e3), runif(20, 1e3, 9e3)))
  nd_w  <- withr::with_seed(10, (1e15 + rnorm(20)) - 1e15)
  newd <- function(o) sf::st_as_sf(data.frame(x = nd_xy[, 1], y = nd_xy[, 2],
                                              p = o + nd_w),
                                   coords = c("x", "y"), crs = 32615)
  ref <- .go_fit(L$mk(L$wq, L$zq), bandwidth = 60)
  pr0 <- suppressWarnings(predict(ref, newd(0), se.fit = TRUE))
  for (o in c(1e8, 1e15)) {
    f  <- .go_fit(L$mk(o + L$wq, L$zq), bandwidth = 60)
    pr <- suppressWarnings(predict(f, newd(o), se.fit = TRUE))
    expect_lt(max(abs(pr$fit - pr0$fit)), 1e-8)
    expect_lt(.go_rel(pr$se.fit, pr0$se.fit), 1e-8)
    # At the training rows a prediction is the fitted value, to the bit.
    tr <- L$mk(o + L$wq, L$zq)
    expect_identical(suppressWarnings(predict(f, tr)), fitted(f))
  }
})

test_that("the bandwidth searched for a predictor at a large offset is the one at 0", {
  L <- .go_layer()
  ref <- .go_fit(L$mk(L$wq, L$zq))
  f   <- .go_fit(L$mk(1e15 + L$wq, L$zq))
  expect_identical(f$info$bandwidth, ref$info$bandwidth)
  expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-8)
})

test_that("a predictor that is constant or collinear at an offset is still dropped from the benchmark", {
  L <- .go_layer()
  d <- L$mk(1e8 + L$wq, L$zq, k = rep(1e8, 300), p2 = 2 * (1e8 + L$wq))
  f <- suppressWarnings(fit_gwr_model(d, "z", c("p", "k"), bandwidth = 60))
  expect_identical(f$engine$lm$rank, 2L)
  expect_true(is.na(coef(f$engine$lm)[["k"]]))
  f2 <- suppressWarnings(fit_gwr_model(d, "z", c("p", "p2"), bandwidth = 60))
  expect_identical(f2$engine$lm$rank, 2L)
  expect_true(is.na(coef(f2$engine$lm)[["p2"]]))
  # The aliased column counts 0 in the intercept at the data's origin, so
  # the global fit is the one on p alone, standard errors included.
  f1 <- .go_fit(d, bandwidth = 60)
  expect_true(is.finite(coef(f2$engine$lm)[[1L]]))
  expect_lt(.go_rel(unname(coef(f2$engine$lm)[1:2]),
                    unname(coef(f1$engine$lm))), 1e-8)
  expect_lt(.go_rel(unname(summary(f2$engine$lm)$coefficients),
                    unname(summary(f1$engine$lm)$coefficients)), 1e-8)
  expect_lt(max(abs(fitted(f2$engine$lm) - fitted(f1$engine$lm))), 1e-8)
})

test_that("a response at 1e15 gives the residuals, RSS, AICc, slopes and benchmark of the values at 0", {
  L <- .go_layer()
  ref <- .go_fit(L$mk(L$wq, L$zq), bandwidth = 60)
  f   <- .go_fit(L$mk(L$wq, 1e15 + L$zq), bandwidth = 60)
  expect_lt(max(abs(residuals(f) - residuals(ref))), 1e-8)
  expect_lt(abs(f$engine$GW.diagnostic$RSS.gw - ref$engine$GW.diagnostic$RSS.gw), 1e-8)
  expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-8)
  expect_lt(.go_rel(f$engine$SDF$p, ref$engine$SDF$p), 1e-8)
  expect_lt(.go_rel(f$engine$SDF$p_SE, ref$engine$SDF$p_SE), 1e-8)
  for (k in c("RSS", "RMSE", "AICc"))
    expect_lt(.go_rel(f$info$global_ols[[k]], ref$info$global_ols[[k]]), 1e-8)
  expect_lt(max(abs(residuals(f$engine$lm) - residuals(ref$engine$lm))), 1e-8)
  # Fitted values, local and global intercepts and predictions in the
  # response's units, to the eighth a double holds there.
  expect_lte(max(abs(fitted(f) - (1e15 + fitted(ref)))), 0.125)
  expect_lte(max(abs(f$engine$SDF$Intercept - (1e15 + ref$engine$SDF$Intercept))), 0.125)
  expect_lte(abs(coef(f$engine$lm)[[1L]] - (1e15 + coef(ref$engine$lm)[[1L]])), 0.125)
  expect_lte(max(abs(fitted(f$engine$lm) - (1e15 + fitted(ref$engine$lm)))), 0.125)
  # The global fit's standard errors and the effects of its terms are those
  # of the values at 0; its first effect is that of the response's size.
  expect_lt(.go_rel(summary(f$engine$lm)$coefficients[, 2L],
                    summary(ref$engine$lm)$coefficients[, 2L]), 1e-8)
  expect_lt(max(abs(f$engine$lm$effects[-1L] - ref$engine$lm$effects[-1L])),
            1e-8)
  expect_lt(.go_rel(f$engine$lm$effects[[1L]], unname(qr.qty(
    f$engine$lm$qr, 1e15 + L$zq)[1L])), 1e-12)
  expect_identical(f$engine$SDF$y, 1e15 + L$zq)
  tr <- L$mk(L$wq, 1e15 + L$zq)
  expect_identical(suppressWarnings(predict(f, tr)), fitted(f))
  expect_true(f$engine$centres$response != 0)
})

test_that("the bandwidth searched for a response at 1e15 is the one at 0, for both score paths", {
  L <- .go_layer()
  for (k in c("bisquare", "gaussian")) {
    ref <- .go_fit(L$mk(L$wq, L$zq), kernel = k)
    f   <- .go_fit(L$mk(L$wq, 1e15 + L$zq), kernel = k)
    expect_identical(f$info$bandwidth, ref$info$bandwidth)
    expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-8)
    expect_lt(max(abs(residuals(f) - residuals(ref))), 1e-8)
  }
})

test_that("gwr_model_selection() at a large offset of the response or a candidate gives the table at 0", {
  L <- .go_layer()
  sel <- function(p, z) {
    d <- L$mk(p, z, q = L$q)
    suppressWarnings(gwr_model_selection(d, "z", c("p", "q")))
  }
  ref <- sel(L$wq, L$zq)
  for (s in list(sel(L$wq, 1e15 + L$zq), sel(1e8 + L$wq, L$zq),
                 sel(1e15 + L$wq, 1e15 + L$zq))) {
    expect_identical(s$table$variables, ref$table$variables)
    expect_lt(max(abs(s$table$criterion - ref$table$criterion)), 1e-8)
    expect_lt(max(abs(s$table$delta - ref$table$delta)), 1e-8)
    expect_identical(s$bandwidth, ref$bandwidth)
  }
})

test_that("a response at an offset with a trend the predictor explains is fitted about its centre", {
  # The raw values show their spread through the trend (1e5 times the
  # predictor), so .regression_centre() is 0; the raw OLS residuals show
  # none, and the response is regressed about .residual_centre().
  L <- .go_layer()
  y0 <- (1e15 + (1e5 * L$wq + L$zq)) - 1e15
  y1 <- 1e15 + y0
  expect_identical(.regression_centre(y1, bounded = FALSE), 0)
  ref <- .go_fit(L$mk(L$wq, y0), bandwidth = 60)
  f   <- .go_fit(L$mk(L$wq, y1), bandwidth = 60)
  expect_true(f$engine$centres$response != 0)
  expect_lt(max(abs(residuals(f) - residuals(ref))), 1e-6)
  expect_lt(abs(f$info$AICc - ref$info$AICc), 1e-6)
  expect_lt(abs(f$info$global_ols$RSS - ref$info$global_ols$RSS), 1e-6)
})

test_that("residual_morans_i() reads the residuals of a response at 1e15 as those at 0", {
  L <- .go_layer()
  ref <- .go_fit(L$mk(L$wq, L$zq), bandwidth = 60)
  f   <- .go_fit(L$mk(L$wq, 1e15 + L$zq), bandwidth = 60)
  m0 <- suppressWarnings(residual_morans_i(ref))
  m1 <- suppressWarnings(residual_morans_i(f))
  expect_false(is.null(m0))
  expect_lt(abs(m1$observed - m0$observed), 1e-8)
})

test_that("a design and a response that show their spread are fitted as they are given", {
  L <- .go_layer()
  for (d in list(L$mk(L$wq, L$zq), L$mk(100 + L$wq, 1e5 + L$zq))) {
    df <- sf::st_drop_geometry(d)
    dsn <- .gwr_design(z ~ p, df)
    cen <- .gwr_centre(dsn$X, dsn$y)
    expect_identical(cen$X, dsn$X)
    expect_identical(cen$y, dsn$y)
    expect_null(cen$x_centres)
    expect_identical(cen$y_centre, 0)
    f <- .go_fit(d, bandwidth = 60)
    sv <- .gwr_local_solve(dsn$X, dsn$y, sf::st_coordinates(d), 60L,
                           "bisquare", TRUE, what = "fit")
    expect_identical(fitted(f), sv$fitted)
    expect_identical(unname(as.matrix(f$engine$SDF@data[, c("Intercept", "p")])),
                     unname(sv$coef))
    expect_identical(residuals(f), df$z - sv$fitted)
    expect_identical(coef(f$engine$lm), coef(stats::lm(z ~ p, data = df)))
    expect_null(f$engine$centres$predictors)
    expect_identical(f$engine$centres$response, 0)
  }
})

test_that("the near-zero-variance note reads the predictor's spread, not its unit", {
  withr::local_seed(2)
  n <- 150
  xy <- cbind(5e5 + runif(n, 0, 5000), 4.5e6 + runif(n, 0, 5000))
  w <- rnorm(n); y <- 1 + 2 * w + rnorm(n, 0, 0.5)
  note <- function(p) {
    d <- sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], y = y, p = p),
                      coords = c("X", "Y"), crs = 32615)
    lines <- capture_spatialkit_log(
      suppressWarnings(fit_gwr_model(d, "y", "p", bandwidth = 60)),
      level = logger::WARN)
    log_has(lines, "near-zero variance")
  }
  for (k in c(0, -30, -46, -50, -60)) expect_false(note(w * 2^k))
  expect_false(note(1e8 + w))
  # A column without spread beyond its rounding is still noted.
  expect_true(note(1 + (seq_len(n) %% 2) * 2^-52))
})

test_that("engine$lm at a predictor offset is the least-squares fit at the data's origin, standard errors included", {
  # The global fit is solved on the centred columns; its intercept is the
  # one at the data's origin, and so are its qr, summary(), vcov(),
  # confint() and predict().  At 1e4 lm() on the raw columns is accurate,
  # and is the reference.
  L <- .go_layer()
  nd <- c(-1, 0, 1.5)
  d <- L$mk(1e4 + L$wq, L$zq)
  f <- .go_fit(d, bandwidth = 60)
  raw <- stats::lm(z ~ p, data = sf::st_drop_geometry(d))
  expect_lt(.go_rel(summary(f$engine$lm)$coefficients,
                    summary(raw)$coefficients), 1e-6)
  expect_lt(.go_rel(confint(f$engine$lm), confint(raw)), 1e-6)
  pf <- predict(f$engine$lm, data.frame(p = 1e4 + nd), se.fit = TRUE)
  pr <- predict(raw, data.frame(p = 1e4 + nd), se.fit = TRUE)
  expect_lt(.go_rel(pf$fit, pr$fit), 1e-8)
  expect_lt(.go_rel(pf$se.fit, pr$se.fit), 1e-6)
  # Further out the raw columns lose the predictor; the reference is the
  # fit at 0, moved to the data's origin.
  ref <- .go_fit(L$mk(L$wq, L$zq), bandwidth = 60)
  V0 <- stats::vcov(ref$engine$lm)
  p0 <- predict(ref$engine$lm, data.frame(p = nd), se.fit = TRUE)
  f6 <- .go_fit(L$mk(1e6 + L$wq, L$zq), bandwidth = 60)
  pf <- predict(f6$engine$lm, data.frame(p = 1e6 + nd), se.fit = TRUE)
  expect_lt(.go_rel(pf$se.fit, p0$se.fit), 1e-8)
  for (o in c(1e8, 1e12)) {
    f <- .go_fit(L$mk(o + L$wq, L$zq), bandwidth = 60)
    M <- rbind(c(1, -o), c(0, 1))
    expect_lt(.go_rel(stats::vcov(f$engine$lm), M %*% V0 %*% t(M)), 1e-6)
    expect_lt(.go_rel(summary(f$engine$lm)$coefficients[2L, ],
                      summary(ref$engine$lm)$coefficients[2L, ]), 1e-8)
    expect_lt(.go_rel(summary(f$engine$lm)$r.squared,
                      summary(ref$engine$lm)$r.squared), 1e-8)
    pf <- predict(f$engine$lm, data.frame(p = o + nd))
    expect_lt(max(abs(pf - p0$fit)), 1e-12 * o)
    expect_lt(.go_rel(stats::anova(f$engine$lm)[["Sum Sq"]],
                      stats::anova(ref$engine$lm)[["Sum Sq"]]), 1e-8)
  }
})

test_that("the collinearity survey of a predictor at a large offset is the one at 0", {
  L <- .go_layer()
  fit2 <- function(o) suppressWarnings(fit_gwr_model(
    L$mk(o + L$wq, L$zq, q = L$q), "z", c("p", "q"), bandwidth = 60))
  ref <- fit2(0)
  for (o in c(1e8, 1e15)) {
    f <- fit2(o)
    expect_lt(.go_rel(f$info$condition_index, ref$info$condition_index), 1e-8)
    expect_lt(.go_rel(f$info$local_collinearity$cn_slopes,
                      ref$info$local_collinearity$cn_slopes), 1e-8)
    expect_identical(f$info$n_local_collinear, ref$info$n_local_collinear)
  }
})
