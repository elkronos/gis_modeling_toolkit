# Means and spreads of scores read at a scale where their arithmetic
# underflows or overflows.
#
# BUG-211: every out-of-sample R2 is scored against the mean of the
# training responses, taken with mean().  Where a long double is no wider
# than a double (as on macOS on Apple silicon), mean() of training responses
# near the largest double overflows: mean(c(x, x, x/2, -x)) is -Inf for x
# the largest double, whose true mean is 0.375 x.  R2 was then NA for that
# fold, for the pooled score and for model_metrics(newdata =), and
# select_features_forward()'s fold-mean baseline dropped that fold's rows.  On a platform with a wider
# long double mean() is finite here, so these tests pass there before the
# fix as well.
#
# BUG-236: select_features_forward()'s check that `fit_fn` ignores its
# second argument compared the scores with all.equal(tolerance = 1e-12),
# which compares absolute differences below 1e-12, so a working `fit_fn` on
# a response in units of 1e-20 drew the warning.
#
# BUG-237: cv_block_size_sweep()'s `fold_sd` and the standard error of the
# "se" rule were sd() of the fold scores, whose squared deviations underflow
# near 1e-160 and overflow near 1e160.
#
# Comparisons at other scales use a relative tolerance of 1e-12: the
# responses are scaled by powers of two, which is exact, and the scores at
# each scale are the same sums on the same bits, but the reading on a power
# of two (.mean_pow2(), .ss_pow2(), .cell_sd()) takes its sums in another
# order than the plain reading at the ordinary scale, which moves the last
# bits by platform.

# The triage's site: 120 points, price = (2 elev + N(0, 1)) times the scale.
.tmr_site <- function(s, n = 120) {
  withr::with_seed(1, {
    site <- sf::st_as_sf(
      data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
      coords = c("x", "y"), crs = 32632)
    site$price <- (2 * site$elev + rnorm(n)) * s
    site$noise <- rnorm(n)
  })
  site
}

# A fit that predicts the training median plus 2 s elev where it is handed
# `elev`, and the median alone otherwise: its predictions scale exactly with
# the response.  With `refuse_empty`, it refuses the empty set, so a sweep
# takes the training-fold mean as its baseline.
.tmr_fit_fn <- function(s, refuse_empty = FALSE) {
  force(s)
  function(train_sf, vars = "elev") {
    if (refuse_empty && !length(vars)) stop("no predictors")
    spatialkit::new_spatial_fit(
      subclass = "tmr_med",
      engine = list(m = stats::median(train_sf$price),
                    b = if ("elev" %in% vars) 2 * s else 0),
      formula = stats::reformulate(if (length(vars)) vars else "1",
                                   response = "price"),
      response_var = "price", predictor_vars = vars, data_sf = train_sf)
  }
}

.tmr_register <- function() {
  registerS3method("predict", "tmr_med", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    object$engine$m + object$engine$b * newdata$elev
  }, envir = asNamespace("spatialkit"))
}

.tmr_rel <- function(a, b) abs(a - b) / abs(b)

test_that(".train_mean() is mean() wherever mean() is finite (BUG-211)", {
  withr::with_seed(11, {
    for (i in 1:300) {
      n <- sample(1:40, 1)
      y <- rnorm(n) * 10^runif(1, -300, 300) + sample(c(0, 1e8, -1e15), 1)
      if (i %% 7 == 0) y[sample(n, 1)] <- NA
      if (i %% 11 == 0) y[sample(n, 1)] <- Inf
      ct <- if (i %% 3 == 0) stats::median(y, na.rm = TRUE) else 0
      if (!is.finite(ct)) ct <- 0
      ref <- if (ct != 0) mean(y[is.finite(y)] - ct) else mean(y[is.finite(y)])
      expect_identical(spatialkit:::.train_mean(y, ct), ref, info = i)
    }
  })
  # No finite value: NaN, as mean() of none.
  expect_identical(spatialkit:::.train_mean(c(NA, Inf)), NaN)
  # The cross-validation records mean() of each fold's training responses.
  .tmr_register()
  site <- .tmr_site(1)
  cv <- suppressMessages(cv_spatial(site, "price", "elev",
                                    fit_fn = .tmr_fit_fn(1), k = 3, seed = 1))
  for (f in unique(cv$predictions$fold)) {
    held <- cv$predictions$`..row_id`[cv$predictions$fold == f]
    tr <- site$price[-held]
    expect_identical(unique(cv$predictions$y_train_mean[cv$predictions$fold == f]),
                     mean(tr), info = f)
  }
})

test_that(".train_mean() takes the mean near the largest double (BUG-211)", {
  x <- .Machine$double.xmax
  expect_equal(spatialkit:::.train_mean(c(x, x, x / 2, -x)), 0.375 * x,
               tolerance = 1e-15)
  # Values less the centre that overflow: read on the halves.
  y <- c(x, -x, x / 2, -x / 2)
  expect_equal(spatialkit:::.train_mean(y, -x / 2), 0.5 * x, tolerance = 1e-15)
  expect_equal(spatialkit:::.train_mean(y / x, -0.5), 0.5, tolerance = 1e-15)
})

test_that("cv_spatial()'s R2 is the same up to the largest double (BUG-211)", {
  .tmr_register()
  cv0 <- suppressMessages(cv_spatial(.tmr_site(1), "price", "elev",
                                     fit_fn = .tmr_fit_fn(1), k = 3, seed = 1))
  # Before the fix, with an 8-byte long double: fold 1's y_train_mean
  # -Inf at k = 1020, its R2 NA and the overall R2 NA; RMSE unchanged.
  expect_equal(cv0$overall$R2, 0.7655072031, tolerance = 1e-9)
  for (k in c(1019, 1020)) {
    s <- 2^k
    cv <- suppressMessages(cv_spatial(.tmr_site(s), "price", "elev",
                                      fit_fn = .tmr_fit_fn(s), folds = cv0$folds))
    lab <- paste("k =", k)
    expect_true(all(is.finite(cv$predictions$y_train_mean)), label = lab)
    expect_lt(max(.tmr_rel(cv$predictions$y_train_mean / s,
                           cv0$predictions$y_train_mean)), 1e-12, label = lab)
    expect_lt(max(.tmr_rel(cv$fold_metrics$R2, cv0$fold_metrics$R2)), 1e-12,
              label = lab)
    expect_lt(.tmr_rel(cv$overall$R2, cv0$overall$R2), 1e-12, label = lab)
    expect_lt(.tmr_rel(cv$overall$RMSE / s, cv0$overall$RMSE), 1e-12, label = lab)
  }
})

test_that("model_metrics(newdata =) scores R2 against a training mean near the largest double (BUG-211)", {
  x <- .Machine$double.xmax
  mk <- function(v, x0) sf::st_as_sf(
    data.frame(px = x0 + seq_along(v) * 10, py = 0, price = v, elev = 0),
    coords = c("px", "py"), crs = 32632)
  registerS3method("predict", "tmr_fix", function(object, newdata = NULL, ...)
    c(0.8, -0.8) * object$engine$s, envir = asNamespace("spatialkit"))
  score <- function(s) {
    f <- spatialkit::new_spatial_fit(
      subclass = "tmr_fix", engine = list(s = s), formula = price ~ elev,
      response_var = "price", predictor_vars = "elev",
      data_sf = mk(c(1, 1, 0.5, -1) * s, 0))
    model_metrics(f, newdata = mk(c(0.9, -0.9) * s, 1000))
  }
  m1 <- score(1)
  # Before the fix, with an 8-byte long double: R2 NA at x.
  mx <- score(x)
  expect_equal(m1$R2, 0.9894806, tolerance = 1e-7)
  expect_lt(.tmr_rel(mx$R2, m1$R2), 1e-12)
  expect_lt(.tmr_rel(mx$RMSE / x, m1$RMSE), 1e-12)
})

test_that("a fit about a centre is scored against a training mean whose differences overflow (BUG-211)", {
  x <- .Machine$double.xmax
  # A fit that records a centre, as fit_gwr_model() does, scored on the
  # values less it.  The centre is -s / 2 where the training rows reach s,
  # so the training responses less it overflow; elsewhere it is 0.
  testthat::local_mocked_bindings(
    .scored_predictions = function(object, newdata, predict_args = list()) {
      s <- object$engine$s; ct <- object$engine$centre
      yc <- newdata$p - ct
      list(yhat = newdata$p, yhat_centred = yc, centre = ct)
    })
  fit_fn <- function(s) function(train_sf) {
    big <- max(abs(train_sf$y)) >= s
    spatialkit::new_spatial_fit(
      subclass = "tmr_centre", engine = list(s = s, centre = if (big) -s / 2 else 0),
      formula = y ~ p, response_var = "y", predictor_vars = "p",
      data_sf = train_sf)
  }
  mk <- function(s) sf::st_as_sf(data.frame(
    e = c(0, 100, 200, 300, 1000, 1100, 1200, 1300), n_ = 0,
    y = c(1, -1, 0.5, -0.5, 0.3, -0.9, 0.1, -0.2) * s,
    p = c(0.8, -0.9, 0.6, -0.4, 0.2, -0.7, 0.2, -0.1) * s),
    coords = c("e", "n_"), crs = 32617)
  folds <- c(1, 1, 1, 1, 2, 2, 2, 2)
  run <- function(s) {
    d <- mk(s)
    list(mm = model_metrics(fit_fn(s)(d[1:4, ]), newdata = d[5:8, ]),
         cv = suppressMessages(cv_spatial(d, "y", "p", fit_fn = fit_fn(s),
                                          folds = folds)))
  }
  r1 <- run(1)
  # Before the fix: the training mean less the centre was Inf, R2 NA.
  rx <- run(x)
  expect_true(is.finite(r1$mm$R2))
  # The baseline at scale 1 is the training mean (0) itself, read about the
  # centre: R2 = 1 - SSE / sum((y - mean(y_train))^2) on rows 5:8.
  y <- c(0.3, -0.9, 0.1, -0.2); p <- c(0.2, -0.7, 0.2, -0.1)
  ytr <- c(1, -1, 0.5, -0.5)
  expect_equal(r1$mm$R2, 1 - sum((y - p)^2) / sum((y - mean(ytr))^2),
               tolerance = 1e-12)
  expect_lt(.tmr_rel(rx$mm$R2, r1$mm$R2), 1e-12)
  expect_true(all(is.finite(r1$cv$fold_metrics$R2)))
  expect_lt(max(.tmr_rel(rx$cv$fold_metrics$R2, r1$cv$fold_metrics$R2)), 1e-12)
  expect_lt(.tmr_rel(rx$cv$overall$R2, r1$cv$overall$R2), 1e-12)
  expect_true("y_train_mean_centred" %in% names(rx$cv$predictions))
  expect_lt(max(.tmr_rel(rx$cv$predictions$y_train_mean_centred / x,
                         r1$cv$predictions$y_train_mean_centred)), 1e-12)
})

test_that("select_features_forward() scores R2 near the largest double as at scale 1 (BUG-211)", {
  skip_on_cran()
  .tmr_register()
  sweep <- function(k, refuse_empty, select_on) {
    s <- 2^k
    suppressWarnings(suppressMessages(select_features_forward(
      .tmr_site(s), "price", c("elev", "noise"),
      fit_fn = .tmr_fit_fn(s, refuse_empty), k = 3, seed = 1, metric = "R2",
      rule = "tol", select_on = select_on, quiet = TRUE)))
  }
  # The honest score on the held-out half.
  a <- sweep(0, FALSE, "split"); b <- sweep(1020, FALSE, "split")
  expect_identical(a$selected, "elev")
  expect_identical(b$selected, a$selected)
  expect_true(is.finite(a$score_holdout))
  expect_lt(.tmr_rel(b$score_holdout, a$score_holdout), 1e-12)
  expect_lt(.tmr_rel(b$score, a$score), 1e-12)
  # The training-fold-mean baseline: built from every fold's rows.
  a <- sweep(0, TRUE, "all"); b <- sweep(1020, TRUE, "all")
  expect_identical(a$baseline$type, "fold mean")
  expect_identical(b$baseline$type, "fold mean")
  expect_lt(abs(b$baseline$score - a$baseline$score), 1e-12)
  expect_identical(b$selected, a$selected)
  s <- 2^1020
  cv <- suppressMessages(cv_spatial(.tmr_site(s), "price", "elev",
                                    fit_fn = .tmr_fit_fn(s), k = 3, seed = 1))
  bp <- spatialkit:::.fs_fold_mean_predictions(list(cv$predictions))
  expect_setequal(bp$`..row_id`, cv$predictions$`..row_id`)
})

# ---------------------------------------------------------------- BUG-236

.tmr_pts <- function() {
  withr::with_seed(7, {
    n <- 60
    pts <- sf::st_as_sf(data.frame(e = runif(n, 0, 4000), n_ = runif(n, 0, 3000),
                                   x1 = rnorm(n), x2 = rnorm(n)),
                        coords = c("e", "n_"), crs = 32617)
    pts$y0 <- 1 + 0.8 * pts$x1 + rnorm(n)
  })
  pts
}

# A linear fit that divides its response by a power of two, so that its
# predictions scale exactly with the response.
.tmr_lm <- function(train_sf, vars) {
  df <- sf::st_drop_geometry(train_sf)
  s <- 2^floor(log2(max(abs(df$y))))
  df$y <- df$y / s
  fml <- stats::reformulate(if (length(vars)) vars else "1", response = "y")
  spatialkit::new_spatial_fit(subclass = "tmr_lm",
                              engine = list(lm = stats::lm(fml, df), s = s),
                              formula = fml, response_var = "y",
                              predictor_vars = vars, data_sf = train_sf)
}

.tmr_register_lm <- function() {
  registerS3method("predict", "tmr_lm", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    as.numeric(stats::predict(object$engine$lm, sf::st_drop_geometry(newdata))) *
      object$engine$s
  }, envir = asNamespace("spatialkit"))
}

.tmr_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

test_that("the check that fit_fn ignores its second argument reads the scores relative to their size (BUG-236)", {
  skip_on_cran()
  .tmr_register_lm()
  pts <- .tmr_pts()
  ignore_fit <- function(train_sf, vars) .tmr_lm(train_sf, character(0))
  ref <- NULL
  for (sc in c(1, 1e-20, 1e20, 2^-530, 2^530)) {
    d <- pts; d$y <- d$y0 * sc
    lab <- paste("scale", sc)
    r <- .tmr_warnings(suppressMessages(select_features_forward(
      d, "y", c("x1", "x2"), fit_fn = .tmr_lm, k = 3, rule = "tol", quiet = TRUE)))
    # Before the fix: the warning at 1e-20 and 2^-530.
    expect_false(any(grepl("ignoring its second argument", r$warnings)), label = lab)
    if (is.null(ref)) ref <- r$value$selected
    expect_identical(r$value$selected, ref, label = lab)
    r2 <- .tmr_warnings(suppressMessages(select_features_forward(
      d, "y", c("x1", "x2"), fit_fn = ignore_fit, k = 3, rule = "tol", quiet = TRUE)))
    expect_true(any(grepl("ignoring its second argument", r2$warnings)), label = lab)
  }
  expect_identical(ref, "x1")
})

test_that("the check that fit_fn ignores its second argument holds where the intercept-only R2 is 0 (BUG-236)", {
  skip_on_cran()
  .tmr_register_lm()
  registerS3method("predict", "tmr_mean", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    rep(object$engine$m, nrow(newdata))
  }, envir = asNamespace("spatialkit"))
  # Predicts the training mean whatever it is handed: each fold's R2 is
  # 1 - TSS / TSS, exactly 0, so the scores equal the null score of 0 and
  # only a comparison that admits equality draws the warning.
  mean_fit <- function(train_sf, vars) spatialkit::new_spatial_fit(
    subclass = "tmr_mean", engine = list(m = mean(train_sf$y)), formula = y ~ 1,
    response_var = "y", predictor_vars = character(0), data_sf = train_sf)
  pts <- .tmr_pts()
  for (sc in c(1, 1e-20, 2^530)) {
    d <- pts; d$y <- d$y0 * sc
    lab <- paste("scale", sc)
    r <- .tmr_warnings(suppressMessages(select_features_forward(
      d, "y", c("x1", "x2"), fit_fn = mean_fit, k = 3, metric = "R2",
      rule = "tol", quiet = TRUE)))
    expect_identical(r$value$baseline$score, 0, label = lab)
    expect_true(any(grepl("ignoring its second argument", r$warnings)), label = lab)
    r2 <- .tmr_warnings(suppressMessages(select_features_forward(
      d, "y", c("x1", "x2"), fit_fn = .tmr_lm, k = 3, metric = "R2",
      rule = "tol", quiet = TRUE)))
    expect_false(any(grepl("ignoring its second argument", r2$warnings)), label = lab)
    expect_identical(r2$value$selected, "x1", label = lab)
  }
})

# ---------------------------------------------------------------- BUG-237

test_that(".sd_any_size() is sd() at ordinary sizes and the spread at any size (BUG-237)", {
  withr::with_seed(3, {
    for (i in 1:200) {
      x <- rnorm(sample(2:12, 1)) * 10^runif(1, -100, 100) + runif(1, -1e6, 1e6)
      expect_identical(spatialkit:::.sd_any_size(x), stats::sd(x), info = i)
    }
  })
  v <- c(1, 2, 3, 4)
  for (s in c(2^-530, 1e-160, 1e-170, 2^-1000, 2^530, 1e160, 2^1000)) {
    expect_lt(.tmr_rel(spatialkit:::.sd_any_size(v * s), stats::sd(v) * s), 1e-12,
              label = paste("scale", s))
  }
  # Deviations whose squares are all 0 are spread; equal values are not.
  expect_identical(spatialkit:::.sd_any_size(c(2, 2, 2) * 1e-170), 0)
  expect_identical(spatialkit:::.sd_any_size(c(1, NA)), NA_real_)
})

test_that("cv_block_size_sweep()'s fold_sd scales with the response (BUG-237)", {
  skip_on_cran()
  .tmr_register_lm()
  pts <- .tmr_pts()
  lm1 <- function(train_sf) .tmr_lm(train_sf, c("x1", "x2"))
  run <- function(sc) {
    d <- pts; d$y <- d$y0 * sc
    suppressWarnings(suppressMessages(cv_block_size_sweep(
      d, "y", c("x1", "x2"), fit_fn = lm1, n_sizes = 3L, k = 3L, seed = 1L,
      quiet = TRUE)))
  }
  base <- run(1)
  expect_true(all(is.finite(base$fold_sd) & base$fold_sd > 0))
  for (sc in c(2^-530, 2^530)) {
    r <- run(sc)
    lab <- paste("scale 2 ^", log2(sc))
    expect_lt(max(.tmr_rel(r$value / sc, base$value)), 1e-12, label = lab)
    # Before the fix: 5.7e-4 relative at 2^-530, Inf at 2^530.
    expect_lt(max(.tmr_rel(r$fold_sd / sc, base$fold_sd)), 1e-12, label = lab)
  }
})

test_that("select_features_forward(rule = \"se\") judges a step the same at any scale (BUG-237)", {
  skip_on_cran()
  .tmr_register_lm()
  d0 <- .tmr_pts()
  withr::with_seed(5, { d0$w <- rnorm(nrow(d0)); d0$y0 <- rnorm(nrow(d0)) })
  run <- function(sc, rule) {
    d <- d0; d$y <- d$y0 * sc
    .tmr_warnings(suppressMessages(select_features_forward(
      d, "y", "w", fit_fn = .tmr_lm, k = 5, seed = 1, rule = rule, quiet = TRUE)))
  }
  # At scale 1 the "se" rule rejects w and `tol` alone accepts it.
  expect_identical(run(1, "se")$value$selected, character(0))
  expect_identical(run(1, "tol")$value$selected, "w")
  for (sc in c(2^-530, 2^530)) {
    lab <- paste("scale 2 ^", log2(sc))
    r <- run(sc, "se")
    # Before the fix, at 2^530: w selected, with the warning that fewer than
    # two folds could be scored.
    expect_identical(r$value$selected, character(0), label = lab)
    expect_false(any(grepl("needs the RMSE of at least", r$warnings)), label = lab)
  }
})

test_that("prep_model_data() measures the shift of the dropped rows at any scale of the response (BUG-237's cause)", {
  # The same cause as BUG-237, found in its pass: the shift is read in SDs
  # of the response, and sd() was Inf at 2^530 and 0 at 2^-560, so the
  # reading and the sentence that states it were left out in silence.
  pts <- .tmr_pts()
  run <- function(sc) {
    d <- pts; d$y <- d$y0 * sc; d$x2[1:6] <- NA
    msg <- NULL
    r <- withCallingHandlers(prep_model_data(d, "y", c("x1", "x2")),
      spatialkit_rows_dropped = function(w) {
        msg <<- conditionMessage(w); invokeRestart("muffleWarning")
      })
    list(sd = attr(r, "dropped")$response_std_diff, msg = msg)
  }
  verdict <- "more than two standard errors, so the rows are not missing at random"
  r1 <- run(1)
  expect_true(is.finite(r1$sd))
  expect_match(r1$msg, verdict, fixed = TRUE)
  for (sc in c(2^-560, 2^530)) {
    lab <- paste("scale 2 ^", log2(sc))
    r <- run(sc)
    # Before the fix: response_std_diff NA and no sentence on the shift.
    expect_true(is.finite(r$sd), label = lab)
    expect_lt(.tmr_rel(r$sd, r1$sd), 1e-12, label = lab)
    expect_match(r$msg, "a difference of 1.18 SD", fixed = TRUE, label = lab)
    expect_match(r$msg, verdict, fixed = TRUE, label = lab)
  }
})
