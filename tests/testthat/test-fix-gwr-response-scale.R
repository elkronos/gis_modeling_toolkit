# fit_gwr_model(), gwr_model_selection() and cv_gwr() on a response whose
# squares overflow or underflow (1e160 or 1e-160 times a spread of 1).  The
# local regressions are linear in the response, so the fitted values were
# the factor times those of the unscaled response, but the residual sum of
# squares, AICc and the bandwidth search squared the raw residuals: at
# 1e-160 AICc lost its bits and the selection's deltas moved, at 1e160 AICc
# was NA, the search took the widest window and gwr_model_selection()
# refused with no cause.  The response is now divided by a power of two
# where its squares leave the range of doubles, and the factor is put back
# on what is returned.  gwr_model_selection() also refuses too few rows, and
# a constant response, in its own words before it converts the rows.

skip_if_not_installed("sp")

# 60 points over 4 km x 3 km, a response with a spatial field and two
# predictors.
.rs_layer <- function(n = 60, seed = 11) {
  withr::with_seed(seed, {
    x <- 600000 + runif(n, 0, 4000); y <- 3990000 + runif(n, 0, 3000)
    a <- rnorm(n); b <- rnorm(n)
    z <- 2 + 1.5 * sin(x / 900) + cos(y / 700) + 0.6 * a - 0.4 * b +
      rnorm(n, sd = 0.2)
  })
  list(n = n, mk = function(f) sf::st_as_sf(
    data.frame(x = x, y = y, a = a, b = b, z = z * f),
    coords = c("x", "y"), crs = 32617))
}

.rs_fit <- function(d, ...) suppressWarnings(fit_gwr_model(d, "z", c("a", "b"), ...))

# The largest relative difference.  A value that is NA where the reference
# is not fails: max() over no pairs would be -Inf and pass any bound.
.rs_rel <- function(a, b) {
  testthat::expect_identical(is.na(a), is.na(b))
  ok <- !is.na(a) & !is.na(b)
  if (!any(ok)) return(Inf)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300))
}

# What a bandwidth search scored, against the search on the unscaled
# response: AICc shifted by 2 n log(f), a cross-validation score (a sum of
# squares) times f^2, NA where a double cannot hold that, and Inf (a
# bandwidth with no criterion) kept.  The curves are matched on bandwidth.
.rs_search <- function(bs, br, f, n) {
  expect_identical(bs$approach, br$approach)
  expect_identical(bs$bandwidth, br$bandwidth)
  expect_identical(bs$within, br$within)
  cs <- bs$curve; cr <- br$curve
  expect_identical(cs$bandwidth, cr$bandwidth)
  got  <- c(bs$criterion, cs$criterion)
  want <- c(br$criterion, cr$criterion)
  expect_identical(is.infinite(got), is.infinite(want))
  fin <- is.finite(want)
  if (identical(br$approach, "CV")) {
    v <- f * (f * want[fin])
    held <- is.finite(v) & v >= .Machine$double.xmin
    expect_true(all(is.na(got[fin][!held])))
    if (any(held))
      expect_lt(max(abs(got[fin][held] - v[held]) / v[held]), 1e-10)
  } else {
    # Values near 4e4: an absolute 1e-6 covers the rounding of the log.
    expect_lt(max(abs(got[fin] - (want[fin] + 2 * n * log(f)))), 1e-6)
  }
}

# The factors: two exact powers of two, and the probe's 1e-160 and 1e160.
# Under an exact power of two only rounding of the logarithm moves AICc;
# under 1e+-160 the response itself is rounded once (relative 1.1e-16).
.rs_factors <- c(2^-530, 2^530, 1e-160, 1e160)

test_that("fit_gwr_model() on a response whose squares overflow or underflow gives the bandwidth, fitted values and AICc of the unscaled response", {
  L <- .rs_layer()
  n <- L$n
  # The searched bandwidth with the default kernel, a supplied count, the
  # other kernels, and a fixed distance.
  args <- list(list(), list(bandwidth = 20), list(kernel = "gaussian"),
               list(kernel = "boxcar"), list(kernel = "exponential"),
               list(adaptive = FALSE, bandwidth = 1500))
  for (a in args) {
    ref <- do.call(.rs_fit, c(list(L$mk(1)), a))
    expect_true(is.finite(ref$info$AICc))
    expect_identical(ref$engine$centres$scale, 1)
    for (f in .rs_factors) {
      d <- L$mk(f)
      w <- character(0)
      fit <- withCallingHandlers(
        do.call(fit_gwr_model, c(list(d, "z", c("a", "b")), a)),
        warning = function(cnd) {
          w <<- c(w, conditionMessage(cnd))
          invokeRestart("muffleWarning")
        })
      expect_identical(fit$info$bandwidth, ref$info$bandwidth)
      # Relative 1e-12: the fit is linear in the response, and the factor
      # 1e+-160 rounds each value once.
      expect_lt(.rs_rel(fitted(fit) / f, fitted(ref)), 1e-12)
      expect_lt(.rs_rel(residuals(fit) / f, residuals(ref)), 1e-10)
      cf <- coef(fit); cr <- coef(ref)
      for (k in c("Intercept", "a", "b"))
        expect_lt(.rs_rel(cf[[k]] / f, cr[[k]]), 1e-10)
      sdf <- fit$engine$SDF@data; sdr <- ref$engine$SDF@data
      expect_lt(.rs_rel(sdf$a_SE / f, sdr$a_SE), 1e-10)
      expect_lt(max(abs(sdf$a_TV - sdr$a_TV)), 1e-8)
      expect_lt(max(abs(sdf$Stud_residual - sdr$Stud_residual)), 1e-8)
      # AICc moves by 2 n log(factor) and by nothing else; an absolute
      # 1e-6 against values near 4e4 covers the rounding of the logarithm.
      expect_lt(abs(fit$info$AICc - (ref$info$AICc + 2 * n * log(f))), 1e-6)
      expect_lt(abs(fit$engine$GW.diagnostic$AICc -
                      (ref$engine$GW.diagnostic$AICc + 2 * n * log(f))), 1e-6)
      expect_equal(fit$info$tr_S, ref$info$tr_S, tolerance = 1e-10)
      expect_equal(fit$engine$GW.diagnostic$gw.R2,
                   ref$engine$GW.diagnostic$gw.R2, tolerance = 1e-10)
      # The global regression: its RMSE times the factor, its AICc shifted.
      expect_lt(.rs_rel(fit$info$global_ols$RMSE / f, ref$info$global_ols$RMSE),
                1e-10)
      expect_lt(abs(fit$info$global_ols$AICc -
                      (ref$info$global_ols$AICc + 2 * n * log(f))), 1e-6)
      expect_lt(.rs_rel(stats::coef(fit$engine$lm) / f,
                        stats::coef(ref$engine$lm)), 1e-10)
      expect_lt(.rs_rel(fitted(fit$engine$lm) / f, fitted(ref$engine$lm)), 1e-10)
      # The search's AICc, in the response's units, is the fit's AICc at
      # the bandwidth it selected.
      if (is.null(a$bandwidth)) {
        .rs_search(fit$info$bandwidth_search, ref$info$bandwidth_search, f, n)
        expect_lt(abs(fit$info$bandwidth_search$criterion - fit$info$AICc),
                  1e-6)
      }
      # anova() of the global fit: its sums of squares times f^2, where a
      # double holds them (1e-160 and 1e160 square beyond its range).
      ss_r <- stats::anova(ref$engine$lm)[["Sum Sq"]]
      ss_f <- f * (f * ss_r)
      if (all(is.finite(ss_f) & ss_f >= .Machine$double.xmin))
        expect_lt(.rs_rel(stats::anova(fit$engine$lm)[["Sum Sq"]], ss_f),
                  1e-10)
      expect_lt(.rs_rel(stats::effects(fit$engine$lm) / f,
                        stats::effects(ref$engine$lm)), 1e-10)
      # AIC and BIC move by 2 n log(factor) as AICc does.
      for (k in c("AIC", "BIC"))
        expect_lt(abs(fit$engine$GW.diagnostic[[k]] -
                        (ref$engine$GW.diagnostic[[k]] + 2 * n * log(f))),
                  1e-6)
      expect_true(fit$engine$centres$scale != 1)
      # The squared quantities cannot be held as doubles at these sizes:
      # NA, with a warning that says why.
      expect_true(is.na(fit$info$sigma2))
      expect_true(is.na(fit$engine$GW.diagnostic$RSS.gw))
      expect_true(is.na(fit$info$global_ols$RSS))
      expect_true(any(grepl("lose their precision", w, fixed = TRUE)))
      expect_false(any(grepl("AICc is undefined|no bandwidth from", w)))
    }
  }
})

test_that("predict() and summary() of a fit on a response near 1e160 are those of the unscaled response times the factor", {
  L <- .rs_layer()
  ref <- .rs_fit(L$mk(1), bandwidth = 20)
  nd <- L$mk(1)[c(3, 17, 40), ]
  pr <- predict(ref, newdata = nd, se.fit = TRUE)
  for (f in c(2^530, 1e160, 2^-530, 1e-160)) {
    fit <- .rs_fit(L$mk(f), bandwidth = 20)
    p <- suppressWarnings(predict(fit, newdata = L$mk(f)[c(3, 17, 40), ],
                                  se.fit = TRUE))
    expect_lt(.rs_rel(p$fit / f, pr$fit), 1e-12)
    # The residual variance is NA in the fit; the standard errors come from
    # the residual standard deviation, which a double holds.
    expect_false(anyNA(c(p$se.fit, p$se.pred, p$residual.scale)))
    expect_lt(.rs_rel(p$se.fit / f, pr$se.fit), 1e-10)
    expect_lt(.rs_rel(p$se.pred / f, pr$se.pred), 1e-10)
    expect_lt(.rs_rel(p$residual.scale / f, pr$residual.scale), 1e-10)
    expect_equal(p$df, pr$df, tolerance = 1e-10)
    # summary() prints the comparison with the global fit, as at 1.
    out <- capture.output(summary(fit))
    expect_true(any(grepl("Global least squares on the same rows: AICc", out,
                          fixed = TRUE)))
    expect_true(any(grepl("This fit: AICc", out, fixed = TRUE)))
  }
})

test_that("a response scaled past the bounds but whose squares a double holds reports its sums of squares and search scores in its units", {
  L <- .rs_layer()
  n <- L$n
  ref <- .rs_fit(L$mk(1))
  # Repeated locations: the search minimises the leave-one-location-out
  # score, a sum of squares.
  rep_layer <- function(f) {
    d <- L$mk(f)[1:30, ]
    d2 <- d
    d2$z <- d2$z + f * withr::with_seed(4, rnorm(30, sd = 0.3))
    rbind(d, d2)
  }
  refr <- .rs_fit(rep_layer(1))
  expect_identical(refr$info$bandwidth_search$approach, "CV")
  for (f in c(2^400, 2^-400)) {
    w <- character(0)
    fit <- withCallingHandlers(
      fit_gwr_model(L$mk(f), "z", c("a", "b")),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      })
    expect_true(fit$engine$centres$scale != 1)
    expect_false(any(grepl("lose their precision", w, fixed = TRUE)))
    expect_lt(.rs_rel(fit$info$sigma2, f * (f * ref$info$sigma2)), 1e-10)
    expect_lt(.rs_rel(fit$engine$GW.diagnostic$RSS.gw,
                      f * (f * ref$engine$GW.diagnostic$RSS.gw)), 1e-10)
    expect_lt(.rs_rel(fit$info$global_ols$RSS,
                      f * (f * ref$info$global_ols$RSS)), 1e-10)
    expect_lt(.rs_rel(stats::anova(fit$engine$lm)[["Sum Sq"]],
                      f * (f * stats::anova(ref$engine$lm)[["Sum Sq"]])),
              1e-10)
    .rs_search(fit$info$bandwidth_search, ref$info$bandwidth_search, f, n)
    fr <- .rs_fit(rep_layer(f))
    .rs_search(fr$info$bandwidth_search, refr$info$bandwidth_search, f,
               2L * 30L)
    expect_false(anyNA(fr$info$bandwidth_search$curve$criterion))
  }
  # Past what a double holds, the scores are NA and the warning names them.
  for (f in c(2^530, 2^-530)) {
    w <- character(0)
    fr <- withCallingHandlers(
      fit_gwr_model(rep_layer(f), "z", c("a", "b")),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      })
    .rs_search(fr$info$bandwidth_search, refr$info$bandwidth_search, f,
               2L * 30L)
    expect_true(is.na(fr$info$bandwidth_search$criterion))
    expect_true(any(grepl("the cross-validation scores of the bandwidth search",
                          w, fixed = TRUE)))
    expect_identical(fr$info$bandwidth, refr$info$bandwidth)
  }
  # The warning gives the order of the spread, not of the values: a
  # response near 2^610 with a spread near 2^601 says 2^601.
  y0 <- L$mk(1)$z
  d <- L$mk(1)
  d$z <- (1000 + (y0 - mean(y0))) * 2^600
  w <- character(0)
  fit <- withCallingHandlers(
    fit_gwr_model(d, "z", c("a", "b"), bandwidth = 20),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  k <- log2(fit$engine$centres$scale)
  msg <- grep("lose their precision", w, fixed = TRUE, value = TRUE)
  expect_length(msg, 1L)
  expect_match(msg, sprintf("the spread of the response about its median is of the order of 2^%d,", k), fixed = TRUE)
})

test_that("a response within the bounds is fitted as it is given", {
  L <- .rs_layer()
  X <- cbind(Intercept = 1, a = c(0.4, -1.1, 0.9, 2.2, -0.3))
  for (f in c(1, 1e-20, 1e20, 1e-90, 1e90)) {
    y <- c(0.3, -1.2, 2.5, 0.7, 1.1) * f
    cen <- .gwr_centre(X, y)
    expect_identical(cen$y_scale, 1)
    expect_identical(cen$y, y)
  }
  # Past the bounds the response is divided by an exact power of two.
  for (f in c(1e-110, 1e110, 2^-530, 2^530)) {
    y <- c(0.3, -1.2, 2.5, 0.7, 1.1) * f
    cen <- .gwr_centre(X, y)
    expect_true(cen$y_scale != 1)
    expect_identical(log2(cen$y_scale), round(log2(cen$y_scale)))
    expect_identical(cen$y * cen$y_scale, y)
  }
  # A fit at 1e-20 and 1e20 (the probe's scale_response) records no scale
  # and raises no warning of the probe's own.
  for (f in c(1e-20, 1e20)) {
    w <- character(0)
    fit <- withCallingHandlers(
      fit_gwr_model(L$mk(f), "z", c("a", "b"), bandwidth = 20),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      })
    expect_identical(fit$engine$centres$scale, 1)
    expect_false(any(grepl("lose their precision", w, fixed = TRUE)))
    expect_true(is.finite(fit$info$sigma2))
  }
  # The sums of squares put back: 0 stays 0, NA where a double cannot
  # hold them, the product where it can.
  expect_identical(.gwr_unscale_ss(c(0, 2, NA), 2^-600), c(0, NA, NA))
  expect_identical(.gwr_unscale_ss(2, 2^600), NA_real_)
  expect_identical(.gwr_unscale_ss(3, 2^-100), 3 * 2^-200)
  expect_identical(.gwr_unscale_ss(c(3, NA), 1), c(3, NA))
})

test_that("gwr_model_selection() on a response whose squares overflow or underflow ranks the models as on the unscaled response", {
  L <- .rs_layer()
  n <- L$n
  sel <- function(d, ...) suppressWarnings(
    gwr_model_selection(d, "z", c("a", "b"), ...))
  for (a in list(list(), list(bw_approach = "CV"), list(bandwidth = 20))) {
    ref <- do.call(sel, c(list(L$mk(1)), a))
    for (f in .rs_factors) {
      w <- character(0)
      s <- withCallingHandlers(
        do.call(gwr_model_selection, c(list(L$mk(f), "z", c("a", "b")), a)),
        warning = function(cnd) {
          w <<- c(w, conditionMessage(cnd))
          invokeRestart("muffleWarning")
        })
      expect_identical(s$bandwidth, ref$bandwidth)
      expect_identical(s$table$variables, ref$table$variables)
      expect_identical(s$best, ref$best)
      expect_lt(max(abs(s$table$delta - ref$table$delta)), 1e-6)
      expect_lt(max(abs(s$table$criterion -
                          (ref$table$criterion + 2 * n * log(f)))), 1e-6)
      if (is.null(a$bandwidth))
        .rs_search(s$bandwidth_search, ref$bandwidth_search, f, n)
      raw <- s$raw[[2]]; rr <- ref$raw[[2]]
      expect_lt(max(abs(raw[, "AIC"] - (rr[, "AIC"] + 2 * n * log(f)))), 1e-6)
      expect_true(all(is.na(raw[, "RSS"])))
      expect_true(any(grepl("lose their precision", w, fixed = TRUE)))
    }
  }
})

test_that("cv_gwr() on a response near 1e160 or 1e-160 gives the predictions of the unscaled response times the factor", {
  L <- .rs_layer()
  base <- L$mk(1)
  folds <- rep(1:3, length.out = L$n)
  cvf <- function(d, ...) suppressWarnings(
    cv_gwr(d, "z", c("a", "b"), folds = folds, ...))
  for (ad in c(TRUE, FALSE)) {
    a <- if (ad) list() else list(adaptive = FALSE)
    ref <- do.call(cvf, c(list(base), a))
    for (f in c(1e-160, 1e160, 2^530)) {
      cv <- do.call(cvf, c(list(L$mk(f)), a))
      pr <- cv$predictions[order(cv$predictions$..row_id), ]
      pref <- ref$predictions[order(ref$predictions$..row_id), ]
      expect_lt(.rs_rel(pr$yhat / f, pref$yhat), 1e-10)
      expect_lt(.rs_rel(cv$fold_metrics$RMSE / f, ref$fold_metrics$RMSE), 1e-10)
      expect_lt(max(abs(cv$fold_metrics$R2 - ref$fold_metrics$R2)), 1e-10)
      expect_lt(.rs_rel(cv$fold_metrics$bandwidth, ref$fold_metrics$bandwidth),
                1e-10)
    }
  }
})

test_that("gwr_model_selection() refuses too few rows, at any bandwidth, before it converts them", {
  L <- .rs_layer()
  d <- L$mk(1)
  for (nr in c(0L, 1L, 3L)) {
    for (ad in c(TRUE, FALSE)) {
      e <- tryCatch(
        suppressWarnings(gwr_model_selection(d[seq_len(nr), ], "z",
                                             c("a", "b"), adaptive = ad,
                                             bandwidth = if (ad) NULL else 1000)),
        error = function(e) conditionMessage(e))
      expect_identical(
        e, sprintf(paste0("gwr_model_selection(): %d observations cannot ",
                          "support a local fit of 2 candidate predictors ",
                          "plus an intercept."), nr))
    }
  }
  # Four rows for two candidates pass the check, and the sweep then finds
  # no model with a residual degree of freedom, and says so.
  e4 <- tryCatch(suppressWarnings(
    gwr_model_selection(d[1:4, ], "z", c("a", "b"), bandwidth = 4)),
    error = function(e) conditionMessage(e))
  expect_match(e4, "no model produced a finite criterion at bandwidth 4",
               fixed = TRUE)
  expect_match(e4, "too few residual degrees of freedom", fixed = TRUE)
})

test_that("gwr_model_selection() names a constant response, and says why no model has a criterion", {
  L <- .rs_layer()
  d <- L$mk(1)
  d$z <- 3.5
  expect_error(gwr_model_selection(d, "z", c("a", "b")),
               "response 'z' is constant \\(a single distinct finite value, 3.5\\)")
  # Five rows at a fixed bandwidth that leaves each window its own row:
  # every model interpolates the data, and the refusal says so.
  d5 <- L$mk(1)[1:5, ]
  e <- tryCatch(suppressWarnings(
    gwr_model_selection(d5, "z", c("a", "b"), adaptive = FALSE,
                        bandwidth = 0.01)),
    error = function(e) conditionMessage(e))
  expect_match(e, "no model produced a finite criterion at bandwidth 0.01",
               fixed = TRUE)
  expect_match(e, "in each model", fixed = TRUE)
  expect_match(e, "`data_sf`", fixed = TRUE)
})
