# Tests for no spread that read the size of the values, and scores of a
# response at a large offset.
#
# .compute_reg_metrics() turned R2 NA wherever the total sum of squares had
# an RMS deviation under 100 eps of the RMS of y, which a response of 1e15
# plus a spread of 2 has (BUG-87).  The Kish ICC of summarize_by_cell() and
# .standardise_by() read a column again for spread only where its standard
# deviation was below eps, so a column of one value up to rounding was scored
# or scaled to a standard deviation of 1 once its size passed about 4
# (BUG-89, BUG-184).  prep_model_data() took the shift of the dropped rows
# as a difference of raw means, which at 1e15 reads their rounding (BUG-99).
# The scores of a GWR whose response fit_gwr_model() took about a centre
# read its fitted values and predictions in the response's units, held to
# the offset's rounding (BUG-188).

# A response of 200 rows and its least-squares fit, both on the grid a
# double holds at 1e15, so o + ys and o + yh are exact at every offset
# used here and keep the same spread.
.fss_metric_data <- function() {
  withr::with_seed(5, {
    n <- 200; w <- rnorm(n); y <- 2 * w + rnorm(n)
  })
  ys <- (1e15 + y) - 1e15
  yh <- as.numeric(stats::fitted(stats::lm(ys ~ w)))
  yh <- (1e15 + yh) - 1e15
  list(ys = ys, yh = yh, w = w)
}

# R2 by the test the package had: NA where the raw sum of squares shows no
# spread on the values' scale.
.fss_old_r2 <- function(y, yhat, baseline = mean(y)) {
  tol <- 100 * .Machine$double.eps
  tss <- sum((y - baseline)^2)
  if (tss > tol^2 * sum(y^2)) 1 - sum((y - yhat)^2) / tss else NA_real_
}

test_that("R2 of a response at a large offset with spread is computed (BUG-87)", {
  D <- .fss_metric_data()
  ref <- spatialkit:::.compute_reg_metrics(D$ys, D$yh, p = 1L)
  expect_true(is.finite(ref$R2))
  for (o in c(1e15, -1e15, 7e14)) {
    m <- spatialkit:::.compute_reg_metrics(o + D$ys, o + D$yh, p = 1L)
    # The raw test finds no spread here; the column read again does.
    expect_true(is.na(.fss_old_r2(o + D$ys, o + D$yh)))
    expect_false(spatialkit:::.constant_column(o + D$ys))
    # TSS is taken about a mean a double holds to an eighth at 1e15.
    expect_equal(m$R2, ref$R2, tolerance = 2e-3)
    expect_equal(m$Adj_R2, 1 - (1 - m$R2) * (200 - 1) / (200 - 1 - 1))
    expect_equal(m$RMSE, ref$RMSE, tolerance = 1e-12)
    # Out of sample: the baseline is the training mean, the test y's spread.
    mo <- spatialkit:::.compute_reg_metrics(o + D$ys, o + D$yh,
                                            y_train_mean = o + 0.25)
    expect_true(is.finite(mo$R2))
  }
})

test_that("R2 stays NA for a response of one value up to rounding at any offset (BUG-87)", {
  withr::with_seed(2, j <- sample(c(-1, 0, 1), 50, replace = TRUE))
  for (a in c(1.2345678, 4.5e6, 1e15)) {
    u <- 2^(floor(log2(a)) - 52)
    y <- a + j * u
    expect_true(spatialkit:::.constant_column(y))
    m <- spatialkit:::.compute_reg_metrics(y, rep(a, 50), p = 1L)
    expect_true(is.na(m$R2))
    expect_true(is.na(m$Adj_R2))
    expect_true(is.na(spatialkit:::.compute_reg_metrics(rep(a, 50), y)$R2))
  }
})

test_that("R2 where the raw test finds spread is computed as before, to the bit (BUG-87)", {
  D <- .fss_metric_data()
  for (o in c(0, 1e8, 1e10, 1e13, -1e12)) {
    m <- spatialkit:::.compute_reg_metrics(o + D$ys, o + D$yh, p = 1L)
    expect_identical(m$R2, .fss_old_r2(o + D$ys, o + D$yh))
  }
  for (k in c(-60, -9, 20)) {
    m <- spatialkit:::.compute_reg_metrics(D$ys * 2^k, D$yh * 2^k,
                                           y_train_mean = 0.1 * 2^k)
    expect_identical(m$R2, .fss_old_r2(D$ys * 2^k, D$yh * 2^k, 0.1 * 2^k))
  }
})

test_that("model_metrics() of a custom fit at 1e15 gives the R2 of the same values at 0 (BUG-87)", {
  D <- .fss_metric_data()
  withr::with_seed(1, xy <- cbind(5e5 + runif(200, 0, 5000), 4.5e6 + runif(200, 0, 5000)))
  e0 <- D$ys - D$yh
  registerS3method("fitted", "fss_fit", function(object, ...)
    sf::st_drop_geometry(object$data_sf)[[object$response_var]] - object$info$res)
  mm <- function(o) {
    d <- sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], v = o + D$ys, w = D$w),
                      coords = c("X", "Y"), crs = 32615)
    f <- spatialkit::new_spatial_fit("fss_fit", list(), v ~ w, "v", "w", d,
                                     info = list(res = e0))
    suppressWarnings(model_metrics(f))
  }
  r0 <- mm(0); r15 <- mm(1e15)
  expect_true(is.finite(r15$R2))
  expect_equal(r15$R2, r0$R2, tolerance = 2e-3)
})

test_that("the Kish ICC reads a column of one value up to rounding as constant at any size (BUG-89)", {
  withr::with_seed(3, {
    n <- 600
    x <- runif(n, 0, 10000); yy <- runif(n, 0, 10000)
    cells <- floor(x / 2000) * 5 + floor(yy / 2000) + 1
    up <- runif(n) < ifelse(cells %% 2 == 0, 0.8, 0.2)
    r0 <- rnorm(n) + 0.5 * (cells %% 3)
  })
  a0 <- 1.2345678
  v0 <- ifelse(up, a0 * (1 + .Machine$double.eps), a0)   # one unit in the last place apart
  icc_of <- function(col, s, extra = NULL) {
    d <- sf::st_as_sf(data.frame(x = x, y = yy, v = v0 * s, r = r0 * s, poly_id = cells),
                      coords = c("x", "y"), crs = 32615)
    out <- suppressWarnings(suppressMessages(
      summarize_by_cell(d, response_var = col, deff = "kish")))
    list(icc = attr(out, "icc"), deff = attr(out, "deff"), out = out)
  }
  r_ref <- icc_of("r", 1)$icc
  for (p in c(0, 4, 10, 20)) {
    s <- 2^p
    expect_true(spatialkit:::.constant_column(v0 * s))
    k <- icc_of("v", s)
    expect_true(is.na(as.numeric(k$icc$resp)))
    expect_identical(unname(k$icc$reason[["v"]]), spatialkit:::.icc_reason_constant)
    # No design effect: the standard errors are the uncorrected ones.
    se <- k$out[["..se_resp_v"]]; sdv <- k$out[["..sd_resp_v"]]
    ok <- is.finite(se) & sdv > 0
    expect_equal(se[ok], sdv[ok] / sqrt(k$out$n[ok]), tolerance = 1e-12)
    # A column with spread gets the same ICC under an exact power-of-two scale.
    expect_identical(as.numeric(icc_of("r", s)$icc$resp), as.numeric(r_ref$resp))
  }
})

test_that("the Kish ICC of a column with spread is as before, to the bit (BUG-89)", {
  icc <- function(v) {
    withr::with_seed(11, {
      n <- 300; x <- runif(n, 0, 6000); y <- runif(n, 0, 6000)
      cells <- floor(x / 2000) * 3 + floor(y / 2000) + 1
      z <- rnorm(n) + 0.7 * (cells %% 2)
    })
    d <- sf::st_as_sf(data.frame(x = x, y = y, z = v(z), poly_id = cells),
                      coords = c("x", "y"), crs = 32615)
    attr(suppressWarnings(suppressMessages(
      summarize_by_cell(d, response_var = "z", deff = "kish"))), "icc")
  }
  # The ICC as .estimate_icc() computes it on the raw values, by hand.
  hand <- function(v) {
    withr::with_seed(11, {
      n <- 300; x <- runif(n, 0, 6000); y <- runif(n, 0, 6000)
      cells <- floor(x / 2000) * 3 + floor(y / 2000) + 1
      z <- rnorm(n) + 0.7 * (cells %% 2)
    })
    raw <- v(z); vals <- (raw - mean(raw)) / sd(raw)
    g <- factor(paste("z", cells, sep = "\r"))
    k <- nlevels(g); N <- length(vals); ni <- tabulate(g)
    gm <- tapply(vals, g, mean)
    MSB <- sum(ni * (gm - mean(vals))^2) / (k - 1)
    MSW <- sum((vals - gm[g])^2) / (N - k)
    n0 <- (N - sum(ni^2) / N) / (k - 1)
    max(0, min(1, (MSB - MSW) / (MSB + (n0 - 1) * MSW)))
  }
  for (v in list(function(z) z, function(z) 1e6 + z, function(z) z * 1e-9,
                 function(z) 4.5e6 + 1e-6 * z))
    expect_identical(as.numeric(icc(v)$resp), hand(v))
})

# 40 rows; the six with p missing have their response shifted by `shift`.
.fss_prep <- function(o, shift) {
  withr::with_seed(2, {
    n <- 40; x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4)
    z0 <- rnorm(n); p <- rnorm(n); p[1:6] <- NA
  })
  z <- z0; z[1:6] <- z[1:6] + shift
  zq <- (1e15 + z) - 1e15
  d <- sf::st_as_sf(data.frame(x = x, y = y, z = o + zq, p = p),
                    coords = c("x", "y"), crs = 32615)
  msg <- NULL
  r <- withCallingHandlers(prep_model_data(d, "z", "p"),
    spatialkit_rows_dropped = function(w) {
      msg <<- conditionMessage(w); invokeRestart("muffleWarning")
    })
  list(sd = attr(r, "dropped")$response_std_diff, msg = msg, zq = zq)
}

test_that("prep_model_data()'s shift of the dropped rows does not move with the response's origin (BUG-99)", {
  verdict <- "more than two standard errors, so the rows are not missing at random"
  for (shift in c(0.6, 0.75, 0.9)) {
    r0  <- .fss_prep(0, shift)
    r15 <- .fss_prep(1e15, shift)
    expect_equal(r15$sd, r0$sd, tolerance = 1e-12)
    expect_identical(grepl(verdict, r15$msg, fixed = TRUE),
                     grepl(verdict, r0$msg, fixed = TRUE))
    expect_identical(sub(".*a difference of ([-0-9.]+) SD.*", "\\1", r15$msg),
                     sub(".*a difference of ([-0-9.]+) SD.*", "\\1", r0$msg))
    # The means in the message stay in the response's units.
    expect_match(r15$msg, "averages 1e\\+15 in the dropped rows")
  }
  # At 0.75 the verdict read the rounding: 0.93 SD at 1e15, 0.87 at 0.
  expect_match(.fss_prep(1e15, 0.75)$msg, "a difference of 0.87 SD.", fixed = TRUE)
  expect_false(grepl(verdict, .fss_prep(1e15, 0.75)$msg, fixed = TRUE))
})

test_that("prep_model_data()'s shift is as before, to the bit, for a response that keeps its spread (BUG-99)", {
  for (o in c(0, 1e5, -3e8)) {
    r <- .fss_prep(o, 0.8)
    y <- o + r$zq
    expect_identical(r$sd, (mean(y[1:6]) - mean(y[-(1:6)])) / stats::sd(c(y[1:6], y[-(1:6)])))
  }
})

test_that(".standardise_by() reads a column of one value up to rounding as without spread at any size (BUG-184)", {
  sb <- spatialkit:::.standardise_by
  withr::with_seed(1, j <- sample(c(-1, 0, 1), 50, replace = TRUE))
  for (a in c(1.2345678, 4.5e6, 2^40, 1e15)) {
    u <- 2^(floor(log2(a)) - 52)
    x <- a + j * u
    s <- sb(x, mean(x), sd(x), none = 1)
    expect_identical(s$scale, 1)
    expect_identical(s$center, mean(x))
    expect_lte(max(abs((x - s$center) / s$scale)), 4 * u)
    expect_identical(sb(x, mean(x), sd(x), none = NA_real_)$scale, NA_real_)
  }
  # Columns with spread keep the mean and standard deviation given.
  withr::with_seed(2, w <- rnorm(80))
  for (x in list(w, 4.5e6 + 1e-6 * w, 1e15 + w, 1e8 + w, w * 2^-70)) {
    s <- sb(x, mean(x), sd(x))
    if (sd(x) >= .Machine$double.eps) {
      expect_identical(s, list(center = mean(x), scale = sd(x)))
    } else {
      expect_equal(s$scale, sd(w) * 2^-70, tolerance = 1e-12)
    }
  }
})

test_that("a predictor or coordinate axis of one value up to rounding is fitted as one without spread (BUG-184)", {
  skip_if_not_installed("brms")
  capture <- function(...) {
    cap <- new.env()
    local_mocked_bindings(
      brm = function(...) { cap$args <- list(...); structure(list(), class = "fss_stub") },
      .package = "brms")
    fit <- suppressWarnings(suppressMessages(fit_bayesian_spatial_model(
      ..., compute_loo = FALSE, check_convergence = FALSE)))
    list(fit = fit, args = cap$args)
  }
  withr::with_seed(4, {
    n <- 120
    xy <- cbind(5e5 + runif(n, 0, 5000), 4.5e6 + runif(n, 0, 5000))
    w <- rnorm(n); y <- 1 + 2 * w + sin(xy[, 1] / 1500) + rnorm(n, 0, 0.5)
    j <- sample(c(-1, 0, 1), n, replace = TRUE)
  })
  for (a in c(1, 4.5e6, 2^40, 1e15)) {
    u <- 2^(floor(log2(a)) - 52)
    r <- capture(
      sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], y = y, p = w, c1 = a + j * u),
                   coords = c("X", "Y"), crs = 32615),
      "y", c("p", "c1"), standardize_predictors = TRUE,
      chains = 1, iter = 200, warmup = 100)
    ps <- r$fit$info$predictor_scaling
    expect_identical(ps$c1$scale, 1)
    expect_lte(max(abs(r$args$data$c1)), 4 * u)
    expect_identical(ps$p$scale, sd(w, na.rm = TRUE))
  }
  # A predictor at 4.5e6 with a spread of 1e-6 keeps its own scale.
  r <- capture(
    sf::st_as_sf(data.frame(X = xy[, 1], Y = xy[, 2], y = y, p = 4.5e6 + 1e-6 * w),
                 coords = c("X", "Y"), crs = 32615),
    "y", "p", standardize_predictors = TRUE, chains = 1, iter = 200, warmup = 100)
  pp <- 4.5e6 + 1e-6 * w
  expect_identical(r$fit$info$predictor_scaling$p,
                   list(center = mean(pp, na.rm = TRUE), scale = sd(pp, na.rm = TRUE)))

  # A transect at x = 4.5e6 whose x carries last-bit jitter reads as the
  # exactly constant one: x takes the scale of y.
  withr::with_seed(8, {
    n <- 150
    yy <- sort(runif(n, 0, 3))
    z <- sin(yy * 4) + rnorm(n, sd = 0.3)
    jx <- sample(c(-1, 0, 1), n, replace = TRUE)
  })
  run <- function(xx) capture(
    sf::st_as_sf(data.frame(X = xx, Y = yy, z = z), coords = c("X", "Y"), crs = 32615),
    "z", character(0), chains = 1, iter = 200, warmup = 100)
  exact  <- run(rep(4.5e6, n))
  jitter <- run(4.5e6 + jx * 2^-30)
  cs <- jitter$fit$info$coord_scaling
  expect_identical(cs$x_scale, cs$y_scale)
  expect_identical(cs$y_scale, exact$fit$info$coord_scaling$y_scale)
  expect_identical(jitter$fit$info$gp_resolution, exact$fit$info$gp_resolution)
  expect_identical(jitter$fit$info$gp_prefit$floor, exact$fit$info$gp_prefit$floor)
})

# The layer of the GWR probe: 300 points, a predictor and a response on the
# grid a double holds at 1e15.
.fss_gwr_layer <- function(o) {
  withr::with_seed(3, {
    n <- 300
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); w <- rnorm(n)
    z <- 2 * w + sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5)
  })
  wq <- (1e15 + w) - 1e15; zq <- (1e15 + z) - 1e15
  sf::st_as_sf(data.frame(x = x, y = y, z = o + zq, p = wq),
               coords = c("x", "y"), crs = 32615)
}

.fss_gwr <- function(d) suppressWarnings(suppressMessages(
  fit_gwr_model(d, "z", "p", bandwidth = 60)))

.fss_cols <- c("RMSE", "MAE", "R2")

test_that("a GWR of a response at a large offset scores as at 0 in sample (BUG-188)", {
  skip_if_not_installed("GWmodel")
  f0 <- .fss_gwr(.fss_gwr_layer(0))
  m0 <- model_metrics(f0)
  for (o in c(1e12, 1e15)) {
    f <- .fss_gwr(.fss_gwr_layer(o))
    yc <- f$engine$centres$response
    expect_true(yc != 0 && abs(yc - o) < 1)
    m <- model_metrics(f)
    expect_equal(unlist(m[.fss_cols]), unlist(m0[.fss_cols]), tolerance = 1e-8)
    expect_equal(unlist(summary(f)$in_sample[.fss_cols]), unlist(m0[.fss_cols]),
                 tolerance = 1e-8)
    ei <- suppressWarnings(evaluate_insample(list(g = f)))
    cm <- suppressWarnings(compare_models(list(g = f)))
    for (tab in list(ei, cm))
      expect_equal(unlist(tab[1, .fss_cols]), unlist(m0[.fss_cols]), tolerance = 1e-8)
    # MAPE and SMAPE: the centred errors over the values in their units.
    y <- sf::st_drop_geometry(f$data_sf)$z
    e <- residuals(f)
    yh_c <- (y - o) - e
    expect_equal(m$MAPE, mean(abs(e / y)) * 100, tolerance = 1e-12)
    expect_equal(m$SMAPE, mean(2 * abs(e) / (abs(y) + abs(yh_c + o))) * 100,
                 tolerance = 1e-12)
    # fitted() is still in the response's units, to the rounding there.
    half_ulp <- 2^(floor(log2(o)) - 53)
    expect_lte(max(abs((fitted(f) - o) - fitted(f0))), half_ulp + 1e-6)
  }
})

test_that("a GWR of a response at a large offset scores as at 0 on new rows and in cross-validation (BUG-188)", {
  skip_if_not_installed("GWmodel")
  d0 <- .fss_gwr_layer(0)
  tr <- 1:220; te <- 221:300
  f0 <- .fss_gwr(d0[tr, ])
  n0 <- suppressWarnings(model_metrics(f0, newdata = d0[te, ]))
  cv0 <- suppressWarnings(suppressMessages(
    cv_gwr(d0, "z", "p", bandwidth = 60, k = 5, seed = 1)))
  folds <- cv0$folds
  gfit <- function(dd, v) fit_gwr_model(dd, "z", v, bandwidth = 60)
  cs0 <- suppressWarnings(suppressMessages(
    cv_spatial(d0, "z", "p", fit_fn = function(dd) gfit(dd, "p"), folds = folds)))
  for (o in c(1e12, 1e15)) {
    d <- .fss_gwr_layer(o)
    f <- .fss_gwr(d[tr, ])
    n <- suppressWarnings(model_metrics(f, newdata = d[te, ]))
    expect_equal(unlist(n[.fss_cols]), unlist(n0[.fss_cols]), tolerance = 1e-8)
    # predict() is still in the response's units, to the rounding there.
    half_ulp <- 2^(floor(log2(o)) - 53)
    expect_lte(max(abs((suppressWarnings(predict(f, d[te, ])) - o) -
                         suppressWarnings(predict(f0, d0[te, ])))), half_ulp + 1e-6)
    cv <- suppressWarnings(suppressMessages(
      cv_gwr(d, "z", "p", bandwidth = 60, k = 5, seed = 1)))
    expect_equal(unlist(cv$overall[.fss_cols]), unlist(cv0$overall[.fss_cols]),
                 tolerance = 1e-8)
    expect_equal(as.matrix(cv$fold_metrics[.fss_cols]),
                 as.matrix(cv0$fold_metrics[.fss_cols]), tolerance = 1e-8)
    expect_lte(max(abs((cv$predictions$yhat - o) - cv0$predictions$yhat)),
               half_ulp + 1e-6)
    expect_true(all(cv$predictions$y_centre != 0 & abs(cv$predictions$y_centre - o) < 1))
    cs <- suppressWarnings(suppressMessages(
      cv_spatial(d, "z", "p", fit_fn = function(dd) gfit(dd, "p"), folds = folds)))
    expect_equal(unlist(cs$overall[.fss_cols]), unlist(cs0$overall[.fss_cols]),
                 tolerance = 1e-8)
    expect_equal(as.matrix(cs$fold_metrics[.fss_cols]),
                 as.matrix(cs0$fold_metrics[.fss_cols]), tolerance = 1e-8)
    # compare_models_cv() pools the same rows.
    expect_equal(unlist(spatialkit:::.cv_overall_metrics(cs$predictions)[.fss_cols]),
                 unlist(spatialkit:::.cv_overall_metrics(cs0$predictions)[.fss_cols]),
                 tolerance = 1e-8)
  }
})

test_that("select_features_forward() scores a GWR of a response at 1e15 as at 0 (BUG-188)", {
  skip_if_not_installed("GWmodel")
  ff <- function(dd, v) fit_gwr_model(dd, "z", v, bandwidth = 60)
  sel <- function(o) suppressWarnings(suppressMessages(select_features_forward(
    .fss_gwr_layer(o), "z", "p", fit_fn = ff, k = 3, select_on = "split",
    quiet = TRUE)))
  s0 <- sel(0); s15 <- sel(1e15)
  expect_identical(s15$selected, s0$selected)
  expect_true(is.finite(s0$score_holdout))
  expect_equal(s15$score_holdout, s0$score_holdout, tolerance = 1e-8)
  expect_equal(s15$score, s0$score, tolerance = 1e-8)
})

test_that("a GWR fitted on its raw response, and other fits, score as before (BUG-188)", {
  skip_if_not_installed("GWmodel")
  d <- .fss_gwr_layer(0)
  f <- .fss_gwr(d)
  expect_identical(f$engine$centres$response, 0)
  y <- sf::st_drop_geometry(d)$z
  expect_identical(model_metrics(f), spatialkit:::.compute_reg_metrics(y, fitted(f)))
  expect_identical(summary(f)$in_sample,
                   spatialkit:::.compute_reg_metrics(y, fitted(f)))
  nd <- suppressWarnings(model_metrics(f, newdata = d[1:50, ]))
  expect_identical(nd, spatialkit:::.compute_reg_metrics(
    y[1:50], suppressWarnings(predict(f, d[1:50, ])), y_train_mean = mean(y)))
  cv <- suppressWarnings(suppressMessages(cv_gwr(d, "z", "p", bandwidth = 60, k = 5, seed = 1)))
  expect_identical(names(cv$predictions), spatialkit:::.cv_prediction_cols)
  pr <- cv$predictions
  met <- spatialkit:::.compute_reg_metrics(pr$y, pr$yhat, y_train_mean = pr$y_train_mean)
  expect_identical(cv$overall$RMSE, met$RMSE)
  expect_identical(cv$overall$R2, met$R2)
})

test_that(".compute_reg_metrics() with a centre reads the centred values and the values' size (BUG-188)", {
  D <- .fss_metric_data()
  o <- 1e15
  m <- spatialkit:::.compute_reg_metrics(D$ys, D$yh, p = 1L, y_train_mean = 0.3,
                                         centre = o)
  ref <- spatialkit:::.compute_reg_metrics(D$ys, D$yh, p = 1L, y_train_mean = 0.3)
  expect_identical(m[c("n", "RMSE", "MAE", "R2", "Adj_R2")],
                   ref[c("n", "RMSE", "MAE", "R2", "Adj_R2")])
  e <- D$ys - D$yh
  expect_equal(m$MAPE, mean(abs(e / (D$ys + o))) * 100)
  expect_equal(m$SMAPE, mean(2 * abs(e) / (abs(D$ys + o) + abs(D$yh + o))) * 100)
  # One centre per row, filtered with the rows.
  yy <- c(D$ys, NA); hh <- c(D$yh, 1)
  cc <- c(rep(o, 100), rep(2 * o, 100), 5)
  mv <- spatialkit:::.compute_reg_metrics(yy, hh, centre = cc)
  expect_identical(mv$n, 200L)
  expect_equal(mv$MAPE, mean(abs(e / (D$ys + cc[1:200]))) * 100)
  expect_error(spatialkit:::.compute_reg_metrics(D$ys, D$yh, centre = c(1, 2)),
               "`centre` must be length 1")
  # The default keeps the arithmetic to the bit.
  expect_identical(spatialkit:::.compute_reg_metrics(D$ys, D$yh, p = 1L, centre = 0),
                   spatialkit:::.compute_reg_metrics(D$ys, D$yh, p = 1L))
})

# The review of the batch: a shift read on a response of one value up to
# rounding, the Kish ICC on the rows with a cell id, cross-validations that
# mix centred and uncentred folds, R2 in select_features_forward(), and
# predict arguments that change nothing.

test_that("prep_model_data() gives no shift for a response of one value up to rounding (BUG-99)", {
  withr::with_seed(2, {
    n <- 40; x <- runif(n, 0, 1e4); yy <- runif(n, 0, 1e4)
    p <- rnorm(n); p[1:6] <- NA
    j <- c(rep(1, 6), sample(c(0, 1), 34, replace = TRUE, prob = c(0.8, 0.2)))
  })
  run <- function(z) {
    d <- sf::st_as_sf(data.frame(x = x, y = yy, z = z, p = p),
                      coords = c("x", "y"), crs = 32615)
    msg <- NULL
    r <- withCallingHandlers(prep_model_data(d, "z", "p"),
      spatialkit_rows_dropped = function(w) {
        msg <<- conditionMessage(w); invokeRestart("muffleWarning")
      })
    list(sd = attr(r, "dropped")$response_std_diff, msg = msg)
  }
  for (a in c(1.2345678, 4.5e6, 1e15)) {
    u <- 2^(floor(log2(a)) - 52)
    z <- a + j * u
    expect_true(spatialkit:::.constant_column(z))
    # The dropped rows sit one unit in the last place above most kept rows.
    expect_gt(stats::sd(z), 0)
    r <- run(z)
    expect_identical(r$sd, NA_real_)
    expect_false(grepl("a difference of", r$msg, fixed = TRUE))
    expect_false(grepl("not missing at random", r$msg, fixed = TRUE))
    # As for a response of one value.
    r1 <- run(rep(a, n))
    expect_identical(r1$sd, NA_real_)
    expect_identical(r$msg, r1$msg)
  }
})

test_that("the Kish ICC reads the rows with a cell id again for spread (BUG-89)", {
  withr::with_seed(3, {
    n <- 600
    x <- runif(n, 0, 10000); yy <- runif(n, 0, 10000)
    cells <- floor(x / 2000) * 5 + floor(yy / 2000) + 1
    r0 <- rnorm(n) + 0.5 * (cells %% 3)
    j <- sample(c(-1, 0, 1), n, replace = TRUE)
    off <- rnorm(30)
  })
  rq <- (1e15 + r0) - 1e15
  icc_of <- function(v, id = cells) {
    d <- sf::st_as_sf(data.frame(x = x, y = yy, v = v, poly_id = id),
                      coords = c("x", "y"), crs = 32615)
    out <- suppressWarnings(suppressMessages(
      summarize_by_cell(d, response_var = "v", deff = "kish")))
    list(icc = attr(out, "icc"), out = out)
  }
  # One value up to rounding on the rows with a cell id, spread only on the
  # rows with none: nothing to estimate the ICC from.
  id2 <- cells; id2[1:30] <- NA
  v2 <- 4.5e6 + j * 2^-30; v2[1:30] <- 4.5e6 + off
  expect_false(spatialkit:::.constant_column(v2))
  expect_true(spatialkit:::.constant_column(v2[-(1:30)]))
  k <- icc_of(v2, id2)
  expect_true(is.na(as.numeric(k$icc$resp)))
  expect_identical(unname(k$icc$reason[["v"]]),
                   "its values differ only on rows that fall in no cell")
  # As for the column of one value on those rows.
  v1 <- rep(4.5e6, n); v1[1:30] <- 4.5e6 + off
  k1 <- icc_of(v1, id2)
  expect_identical(k$icc$reason, k1$icc$reason)
  expect_identical(is.na(k$out[["..se_resp_v"]]), is.na(k1$out[["..se_resp_v"]]))
  # A column with spread at 1e15 gets the ICC of the same values at 0.
  i0 <- as.numeric(icc_of(rq)$icc$resp)
  expect_true(is.finite(i0) && i0 > 0)
  for (o in c(1e15, -1e15))
    expect_equal(as.numeric(icc_of(o + rq)$icc$resp), i0, tolerance = 1e-12)
  for (s in c(2^-30, 2^30))
    expect_equal(as.numeric(icc_of((1e15 + rq) * s)$icc$resp), i0, tolerance = 1e-12)
})

test_that("a cross-validation whose folds mix centred and uncentred GWR fits scores as at 0 (BUG-188)", {
  skip_if_not_installed("GWmodel")
  d0 <- .fss_gwr_layer(0)
  z0 <- sf::st_drop_geometry(d0)$z
  # About a thousand times the spread: some folds' fits take the response
  # about a centre and others do not.
  o <- 1.0005e3 * spatialkit:::.spread_scaled(z0)$spread
  cvr <- function(d) suppressWarnings(suppressMessages(
    cv_gwr(d, "z", "p", bandwidth = 60, k = 5, seed = 1)))
  cv0 <- cvr(d0)
  cv <- cvr(.fss_gwr_layer(o))
  pr <- cv$predictions
  expect_true(all(spatialkit:::.cv_centred_cols %in% names(pr)))
  expect_true(any(pr$y_centre == 0) && any(pr$y_centre != 0))
  expect_equal(unlist(cv$overall[.fss_cols]), unlist(cv0$overall[.fss_cols]),
               tolerance = 1e-8)
  expect_equal(as.matrix(cv$fold_metrics[.fss_cols]),
               as.matrix(cv0$fold_metrics[.fss_cols]), tolerance = 1e-8)
  # The percentage errors divide by the values in the response's units.
  e <- pr$y - pr$yhat
  expect_equal(cv$overall$MAPE, mean(abs(e / pr$y)) * 100, tolerance = 1e-8)
  expect_equal(cv$overall$SMAPE,
               mean(2 * abs(e) / (abs(pr$y) + abs(pr$yhat))) * 100,
               tolerance = 1e-8)
  # The pooled frame alone gives the same scores.
  expect_equal(unlist(spatialkit:::.cv_overall_metrics(pr)[.fss_cols]),
               unlist(cv0$overall[.fss_cols]), tolerance = 1e-8)
})

test_that("select_features_forward() with metric R2 scores a GWR at 1e15 as at 0 (BUG-188)", {
  skip_if_not_installed("GWmodel")
  ff <- function(dd, v) fit_gwr_model(dd, "z", v, bandwidth = 60)
  sel <- function(o) suppressWarnings(suppressMessages(select_features_forward(
    .fss_gwr_layer(o), "z", "p", fit_fn = ff, k = 3, select_on = "split",
    metric = "R2", quiet = TRUE)))
  s0 <- sel(0); s15 <- sel(1e15)
  expect_identical(s15$selected, s0$selected)
  expect_true(is.finite(s0$score_holdout) && s0$score_holdout < 0.99)
  expect_equal(s15$score_holdout, s0$score_holdout, tolerance = 1e-8)
  expect_equal(s15$score, s0$score, tolerance = 1e-8)
})

test_that("a GWR at 1e15 given se.fit = FALSE scores as without it (BUG-188)", {
  skip_if_not_installed("GWmodel")
  d0 <- .fss_gwr_layer(0); d <- .fss_gwr_layer(1e15)
  tr <- 1:220; te <- 221:300
  f0 <- .fss_gwr(d0[tr, ]); f <- .fss_gwr(d[tr, ])
  n0 <- suppressWarnings(model_metrics(f0, newdata = d0[te, ]))
  n <- suppressWarnings(model_metrics(f, newdata = d[te, ], se.fit = FALSE))
  expect_equal(unlist(n[.fss_cols]), unlist(n0[.fss_cols]), tolerance = 1e-8)
  folds <- suppressWarnings(suppressMessages(
    cv_gwr(d0, "z", "p", bandwidth = 60, k = 5, seed = 1)))$folds
  fs <- function(dd) fit_gwr_model(dd, "z", "p", bandwidth = 60)
  cs <- function(dd, ...) suppressWarnings(suppressMessages(
    cv_spatial(dd, "z", "p", fit_fn = fs, folds = folds, ...)))
  cs0 <- cs(d0)
  cs15 <- cs(d, predict_args = list(se.fit = FALSE))
  expect_equal(unlist(cs15$overall[.fss_cols]), unlist(cs0$overall[.fss_cols]),
               tolerance = 1e-8)
  expect_equal(as.matrix(cs15$fold_metrics[.fss_cols]),
               as.matrix(cs0$fold_metrics[.fss_cols]), tolerance = 1e-8)
  # Any other argument still goes to predict(), which refuses it.
  expect_error(suppressWarnings(model_metrics(f, newdata = d[te, ], level = 0.9)),
               "level")
})
