# tests/testthat/test-fix-gwr-se-extrapolation.R
# ---------------------------------------------------------------------------
# What the local solve gives the methods of a gwr_fit beside the estimates:
# local standard errors through coef(se = TRUE), the standard error of a
# prediction through predict(se.fit = TRUE), and a count of the predictions
# made farther from the training data than any window of the fit reached.
# ---------------------------------------------------------------------------

skip_if_not_installed("sp")

.se_catch <- function(expr) {
  w <- character(0); cls <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); cls <<- c(cls, class(cnd)[1L])
      invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w, classes = cls)
}

# The slope of `a` rises eastward; training extent 0-1000 in x and y.
.se_layer <- function(n = 250, seed = 3) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- sf::st_as_sf(data.frame(x = x, y = y, a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + (1 + x / 500) * d$a - 0.5 * d$b + rnorm(n, 0, 0.4)
  d
}
.se_new <- function(m, xlim, seed = 8) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = runif(m, xlim[1], xlim[2]), y = runif(m, 0, 1000),
                          a = rnorm(m), b = rnorm(m)),
               coords = c("x", "y"), crs = 32632)
}


# ---------------------------------------------------------------------------
# coef(se = TRUE)
# ---------------------------------------------------------------------------

test_that("coef(se = TRUE) returns the estimates, GWmodel's standard errors and the t values", {
  skip_if_not_installed("GWmodel")
  d <- .se_layer()
  for (cs in list(list(bw = 40L, adaptive = TRUE, kernel = "bisquare"),
                  list(bw = 180, adaptive = FALSE, kernel = "gaussian"))) {
    fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), adaptive = cs$adaptive,
                                          bandwidth = cs$bw, kernel = cs$kernel))
    g <- GWmodel::gwr.basic(z ~ a + b, data = .to_sp(d, c("z", "a", "b")), bw = cs$bw,
                            kernel = cs$kernel, adaptive = cs$adaptive)$SDF@data
    out <- coef(fit, se = TRUE)
    expect_named(out, c("estimate", "se", "t"))
    for (el in out) {
      expect_s3_class(el, "data.frame")
      expect_named(el, c("Intercept", "a", "b"))
      expect_identical(nrow(el), nrow(d))
    }
    # The default is unchanged: the estimates alone, as a data.frame.
    expect_identical(coef(fit), out$estimate)
    expect_identical(coef(fit, se = FALSE), out$estimate)
    expect_equal(unname(as.matrix(out$se)),
                 unname(as.matrix(g[c("Intercept_SE", "a_SE", "b_SE")])), tolerance = 1e-6)
    expect_equal(unname(as.matrix(out$t)),
                 unname(as.matrix(g[c("Intercept_TV", "a_TV", "b_TV")])), tolerance = 1e-6)
    expect_equal(as.matrix(out$t), as.matrix(out$estimate) / as.matrix(out$se))
    # From the fit's own numbers: sigma^2 = RSS / (n - 2 tr S + tr S'S).
    expect_equal(fit$info$sigma2, sum(residuals(fit)^2) / fit$info$edf)
    expect_equal(fit$info$enp + fit$info$edf, nrow(d))
  }
  expect_error(coef(fit, se = NA), "`se` must be TRUE or FALSE")
  expect_error(coef(fit, se = "yes"), "`se` must be TRUE or FALSE")
})

test_that("coef(se = TRUE) takes the standard errors by position, and refuses a table without them", {
  # A predictor called `a_SE` beside `a`: the coefficient of `a_SE` is one
  # column, the standard error of `a` another of the same name.
  set.seed(4); n <- 120
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n), a_SE = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + d$a + 2 * d$a_SE + rnorm(n, 0, 0.3)
  fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "a_SE"), bandwidth = 40))
  out <- coef(fit, se = TRUE)
  expect_named(out$estimate, c("Intercept", "a", "a_SE"))
  expect_lt(abs(mean(out$estimate$a_SE) - 2), 0.1)
  expect_true(all(out$se$a > 0 & out$se$a < 0.3))
  expect_equal(out$t$a, out$estimate$a / out$se$a)
  # An engine that carries the coefficients alone.
  bare <- fit
  bare$engine <- list(SDF = sf::st_sf(Intercept = 1, a = 2, a_SE = 3,
                                      geometry = sf::st_geometry(d)))
  expect_identical(nrow(coef(bare)), 120L)
  expect_error(coef(bare, se = TRUE), "does not carry standard errors and t values")
})


# ---------------------------------------------------------------------------
# predict(se.fit = TRUE)
# ---------------------------------------------------------------------------

test_that("predict(se.fit = TRUE) returns the fit's values with their standard errors", {
  d <- .se_layer()
  fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 40))
  nd <- .se_new(60, c(100, 900))
  p <- predict(fit, newdata = nd, se.fit = TRUE)
  expect_named(p, c("fit", "se.fit", "se.pred", "residual.scale", "df"))
  expect_identical(p$fit, predict(fit, newdata = nd))
  expect_true(all(p$se.fit > 0 & p$se.fit < p$se.pred))
  expect_equal(p$se.pred^2, p$se.fit^2 + fit$info$sigma2)
  expect_equal(p$residual.scale, sqrt(fit$info$sigma2))
  expect_equal(p$df, fit$info$edf)
  # At the training rows the fitted value's variance is sigma^2 times the sum
  # of squares of its hat row, and newdata = NULL means those rows.
  p0 <- predict(fit, se.fit = TRUE)
  expect_identical(p0$fit, unname(fitted(fit)))
  expect_identical(p0, predict(fit, newdata = d, se.fit = TRUE))
  sv <- .gwr_local_solve(cbind(Intercept = 1, a = d$a, b = d$b), d$z,
                         sf::st_coordinates(d), 40L, "bisquare", TRUE, what = "fit")
  expect_equal(p0$se.fit, sqrt(fit$info$sigma2 * sv$var_fit))
  expect_equal(sum(sv$var_fit), 2 * fit$info$tr_S - fit$info$enp)
  # Rows dropped for a missing predictor keep their place, as NA.
  nd2 <- nd[1:8, ]; nd2$a[c(2, 5)] <- NA
  p2 <- predict(fit, newdata = nd2, se.fit = TRUE)
  expect_identical(which(is.na(p2$fit)), c(2L, 5L))
  expect_identical(which(is.na(p2$se.pred)), c(2L, 5L))
  expect_identical(p2$fit[-c(2, 5)], p$fit[c(1, 3, 4, 6:8)])
  # The default return is the plain vector, and a bad `se.fit` is refused.
  expect_type(predict(fit, newdata = nd), "double")
  expect_error(predict(fit, newdata = nd, se.fit = NA), "`se.fit` must be TRUE or FALSE")
  expect_error(predict(fit, newdata = nd, interval = "prediction"), "unused argument")
})

test_that("the prediction standard error covers held-out observations about as often as it says", {
  # A stationary layer, so the local model is true everywhere: 90% intervals
  # from se.pred should hold about 90% of new observations.
  cover <- vapply(1:6, function(seed) {
    set.seed(seed); n <- 300
    x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n)
    d <- sf::st_as_sf(data.frame(x = x, y = y, a = a, z = 1 + 2 * a + rnorm(n, 0, 0.5)),
                      coords = c("x", "y"), crs = 32632)
    tr <- d[1:200, ]; te <- d[201:300, ]
    fit <- suppressMessages(fit_gwr_model(tr, "z", "a", bandwidth = 60))
    p <- suppressWarnings(predict(fit, newdata = te, se.fit = TRUE))
    mean(abs(te$z - p$fit) <= stats::qnorm(0.95) * p$se.pred)
  }, numeric(1))
  expect_gt(mean(cover), 0.85)
  expect_lt(mean(cover), 0.95)
})

test_that("a fit without a residual variance returns NA standard errors and says why", {
  d <- .se_layer(120)
  nd <- .se_new(5, c(100, 900))
  # Put together by hand: the data and the settings, nothing from a solve.
  bare <- new_spatial_fit("gwr_fit", engine = list(), formula = z ~ a + b,
                          response_var = "z", predictor_vars = c("a", "b"), data_sf = d,
                          info = list(bandwidth = 40, adaptive = TRUE, kernel = "bisquare"))
  r <- .se_catch(predict(bare, newdata = nd, se.fit = TRUE))
  expect_false(anyNA(r$value$fit))
  expect_true(all(is.na(r$value$se.fit)) && all(is.na(r$value$se.pred)))
  expect_true(is.na(r$value$residual.scale))
  expect_match(r$warnings, "the fit carries no residual variance")
  fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 40))
  expect_identical(r$value$fit, predict(fit, newdata = nd))
})


# ---------------------------------------------------------------------------
# summary(): effective parameters and the global benchmark
# ---------------------------------------------------------------------------

test_that("the fit carries a global least-squares benchmark and summary() prints it", {
  d <- .se_layer()
  fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 40))
  ols <- stats::lm(z ~ a + b, data = sf::st_drop_geometry(d))
  rss <- sum(stats::residuals(ols)^2); n <- nrow(d)
  g <- fit$info$global_ols
  expect_named(g, c("RSS", "RMSE", "AICc"))
  expect_equal(g$RSS, rss)
  expect_equal(g$RMSE, sqrt(rss / n))
  # The fit's own AICc formula, at a hat-matrix trace of 3.
  expect_equal(g$AICc, n * log(rss / n) + n * log(2 * pi) + n * (n + 3) / (n - 2 - 3))
  # The slope of `a` doubles across this layer: the GWR is far below.
  expect_lt(fit$info$AICc, g$AICc - 50)
  txt <- paste(utils::capture.output(print(summary(fit))), collapse = "\n")
  expect_match(txt, sprintf("Effective parameters: %.1f \\(residual degrees of freedom %.1f\\)",
                            fit$info$enp, fit$info$edf))
  expect_match(txt, sprintf("Global least squares on the same rows: AICc %.2f, RMSE %.4f",
                            g$AICc, g$RMSE))
  expect_match(txt, sprintf("This fit: AICc %.2f \\(below the global fit's by %.2f\\)",
                            fit$info$AICc, g$AICc - fit$info$AICc))
  expect_no_match(txt, "Rank-deficient windows")
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"),
               sprintf("Eff. parameters: %.1f of 250 observations", fit$info$enp))
  # On a process whose coefficients do not vary, the selected GWR does not
  # beat the global regression by AICc in most draws.
  beats <- vapply(1:6, function(seed) {
    set.seed(300 + seed); m <- 150
    s <- sf::st_as_sf(data.frame(x = runif(m, 0, 1000), y = runif(m, 0, 1000), a = rnorm(m)),
                      coords = c("x", "y"), crs = 32632)
    s$z <- 1 + 2 * s$a + rnorm(m, 0, 0.5)
    f <- suppressMessages(fit_gwr_model(s, "z", "a"))
    f$info$AICc < f$info$global_ols$AICc
  }, logical(1))
  expect_lte(sum(beats), 2L)
  # A fit with windows of deficient rank says so in both prints.
  set.seed(42); m <- 300
  x <- runif(m, 0, 1000); y <- runif(m, 0, 1000); a <- rnorm(m)
  r <- sf::st_as_sf(data.frame(x = x, y = y, a = a, region = ifelse(x > 500, 1.3, 0.7),
                               z = 1 + 2 * a + 3 * (x > 500) + rnorm(m, 0, 0.3)),
                    coords = c("x", "y"), crs = 32632)
  fr <- suppressWarnings(suppressMessages(fit_gwr_model(r, "z", c("a", "region"), bandwidth = 40)))
  expect_match(paste(utils::capture.output(print(fr)), collapse = "\n"),
               "Rank-deficient windows: 185 of 300 \\(their unestimable coefficients are NA; see \\$info\\$nonfinite_coef\\)")
  expect_match(paste(utils::capture.output(print(summary(fr))), collapse = "\n"),
               "Rank-deficient windows: 185 of 300 \\(their unestimable coefficients are NA\\)")
})


# ---------------------------------------------------------------------------
# Predictions beyond the reach of the fit
# ---------------------------------------------------------------------------

test_that("predict() counts the locations farther from the data than the widest kernel reaches", {
  # Adaptive bisquare, 40 neighbours.  Before: new locations three extents
  # east of the data (x from 4000 to 5000) came back with no NA and no
  # warning, at 35 times the in-extent error.
  d <- .se_layer()
  fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 40))
  xy <- sf::st_coordinates(d)
  dm <- as.matrix(stats::dist(xy))
  width <- max(apply(dm, 1L, function(v) sort(v)[40L]))
  expect_equal(fit$info$max_kernel_width, width)
  nearest <- function(nd) {
    q <- sf::st_coordinates(nd)
    vapply(seq_len(nrow(q)), function(i)
      min(sqrt((xy[, 1] - q[i, 1])^2 + (xy[, 2] - q[i, 2])^2)), numeric(1))
  }
  # Inside the extent: nothing to say.
  inside <- .se_new(200, c(0, 1000))
  expect_no_warning(p_in <- predict(fit, newdata = inside))
  expect_false(anyNA(p_in))
  # A band that straddles the reach: some locations beyond it, some not.
  band <- .se_new(200, c(1000, 1000 + 2 * width))
  far_b <- nearest(band) > width
  expect_gt(sum(far_b), 20L); expect_gt(sum(!far_b), 20L)
  r <- .se_catch(predict(fit, newdata = band))
  expect_false(anyNA(r$value))
  expect_identical(r$classes, "spatialkit_gwr_extrapolation")
  expect_match(r$warnings, sprintf(
    "^predict.gwr_fit\\(\\): %d of 200 location\\(s\\) are farther from every training point than the widest kernel of the fit reaches \\(%s in the units of the fit's CRS\\)",
    sum(far_b), format(signif(width, 4))))
  expect_match(r$warnings, sprintf("the nearest training point is %s to %s away",
                                   format(signif(min(nearest(band)[far_b]), 4)),
                                   format(signif(max(nearest(band)[far_b]), 4))))
  expect_match(r$warnings, "the prediction\\(s\\) there are extrapolations in space")
  # Three extents out: every one.
  beyond <- .se_new(150, c(4000, 5000))
  rb <- .se_catch(predict(fit, newdata = beyond))
  expect_false(anyNA(rb$value))
  expect_match(rb$warnings, "^predict.gwr_fit\\(\\): 150 of 150 location\\(s\\) are farther")
  # The class is there to be caught.
  n_caught <- 0L
  withCallingHandlers(
    p_far <- predict(fit, newdata = beyond),
    spatialkit_gwr_extrapolation = function(cnd) {
      n_caught <<- n_caught + 1L; invokeRestart("muffleWarning")
    })
  expect_identical(n_caught, 1L)
  expect_identical(p_far, rb$value)
  # The values are returned as computed, and they are what the warning is
  # about: the slope keeps rising in the truth and the fit cannot know.
  truth <- function(nd) 1 + (1 + sf::st_coordinates(nd)[, 1] / 500) * nd$a - 0.5 * nd$b
  rmse <- function(e) sqrt(mean(e^2))
  expect_gt(rmse(rb$value - truth(beyond)), 5 * rmse(p_in - truth(inside)))
})

test_that("the reach is the bandwidth for a fixed kernel with an edge, and three bandwidths for one without", {
  d <- .se_layer()
  out <- .se_new(60, c(1000, 2500), seed = 2)
  xy <- sf::st_coordinates(d); q <- sf::st_coordinates(out)
  dq <- vapply(seq_len(nrow(q)), function(i)
    sqrt((xy[, 1] - q[i, 1])^2 + (xy[, 2] - q[i, 2])^2), numeric(nrow(xy)))
  near <- apply(dq, 2L, min)
  # Bisquare, fixed 300: a location with no training point within 300 m has
  # an empty window, and one with fewer than three has too few for three
  # parameters.  Both are NA, which the other warning counts by cause, and
  # nothing finite is beyond the reach.
  fb <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), adaptive = FALSE, bandwidth = 300))
  within <- colSums(dq < 300)
  expect_gt(sum(within == 0L), 5L); expect_gt(sum(within %in% 1:2), 0L)
  rb <- .se_catch(predict(fb, newdata = out))
  expect_identical(is.na(rb$value), within < 3L)
  expect_identical(rb$classes, "spatialkit_gwr_not_estimable")
  expect_match(rb$warnings, sprintf("At %d of them no training point carries weight within the fixed bandwidth \\(300\\)",
                                    sum(within == 0L)))
  expect_match(rb$warnings, sprintf("At %d the training points in the window cannot place the new row",
                                    sum(within %in% 1:2)))
  # Gaussian, fixed 150: every location is predicted, and those more than 450
  # m from the data are counted.
  fg <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), adaptive = FALSE, bandwidth = 150,
                                       kernel = "gaussian"))
  rg <- .se_catch(predict(fg, newdata = out))
  expect_false(anyNA(rg$value))
  expect_identical(rg$classes, "spatialkit_gwr_extrapolation")
  expect_match(rg$warnings, sprintf("^predict.gwr_fit\\(\\): %d of 60 location", sum(near > 450)))
  expect_match(rg$warnings, "reaches \\(450 in the units of the fit's CRS, three times the widest gaussian bandwidth of 150\\)")
  expect_identical(.gwr_kernel_reach(fg, xy, 150, FALSE, "gaussian"), list(width = 150, reach = 450))
  expect_identical(.gwr_kernel_reach(fb, xy, 300, FALSE, "bisquare"), list(width = 300, reach = 300))
})

test_that("a fit put together by hand has its widest kernel worked out at prediction", {
  d <- .se_layer(150)
  fit <- suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 30))
  bare <- new_spatial_fit("gwr_fit", engine = list(), formula = z ~ a + b,
                          response_var = "z", predictor_vars = c("a", "b"), data_sf = d,
                          info = list(bandwidth = 30, adaptive = TRUE, kernel = "bisquare"))
  xy <- sf::st_coordinates(d)
  expect_equal(.gwr_kernel_reach(bare, xy, 30L, TRUE, "bisquare")$width,
               fit$info$max_kernel_width)
  beyond <- .se_new(20, c(3000, 4000))
  r1 <- .se_catch(predict(fit, newdata = beyond))
  r2 <- .se_catch(predict(bare, newdata = beyond))
  expect_identical(r2$value, r1$value)
  expect_identical(r2$warnings, r1$warnings)
  expect_identical(r2$classes, "spatialkit_gwr_extrapolation")
})
