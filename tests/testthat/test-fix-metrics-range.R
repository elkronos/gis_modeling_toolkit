# Scores of a response whose squares underflow or overflow.
#
# .compute_reg_metrics() took RSS and TSS on the values as they stand.  For
# errors below about 1e-162 every squared error fell below the smallest
# double: RSS was 0, so RMSE was 0 and R2 1 whatever the fit, or R2 NA where
# TSS fell to 0 too.  For squared errors that sum past the largest double
# both sums were Inf: RMSE Inf and R2 NA.  MAE, MAPE and SMAPE were right at
# every size (BUG-209).

# 200 rows and their least-squares fit, as in the triage's probe.
.fmr_data <- function() {
  withr::with_seed(1, {
    n <- 200; w <- rnorm(n); y <- 2 * w + rnorm(n)
  })
  list(y = y, yh = as.numeric(stats::fitted(stats::lm(y ~ w))), w = w, n = n)
}

.fmr_k <- c(-900, -560, -540, -520, 0, 510, 530, 1000)

# The scores that scale with the response, divided by the scale, and those
# that do not.
.fmr_read <- function(m, s) {
  c(RMSE = m$RMSE / s, MAE = m$MAE / s, MAPE = m$MAPE, SMAPE = m$SMAPE,
    R2 = m$R2, Adj_R2 = m$Adj_R2)
}

test_that("RMSE, MAE and R2 do not depend on the response's scale (BUG-209)", {
  D <- .fmr_data()
  base <- c(1, -2)[rep(1:2, D$n / 2)]
  cases <- list(
    plain    = function(s) spatialkit:::.compute_reg_metrics(D$y * s, D$yh * s, p = 1L),
    scalar   = function(s) spatialkit:::.compute_reg_metrics(
      D$y * s, D$yh * s, p = 1L, y_train_mean = 0.1 * s),
    per_row  = function(s) spatialkit:::.compute_reg_metrics(
      D$y * s, D$yh * s, p = 1L, y_train_mean = 0.1 * base * s),
    centred  = function(s) spatialkit:::.compute_reg_metrics(
      D$y * s, D$yh * s, p = 1L, centre = 7 * s),
    centred_rows = function(s) spatialkit:::.compute_reg_metrics(
      D$y * s, D$yh * s, p = 1L, y_train_mean = 0.1 * base * s,
      centre = 7 * base * s)
  )
  for (nm in names(cases)) {
    ref <- .fmr_read(cases[[nm]](1), 1)
    expect_true(all(is.finite(ref)), info = nm)
    expect_true(ref[["R2"]] > 0.7 && ref[["R2"]] < 0.8, info = nm)
    for (k in .fmr_k) {
      s <- 2^k
      expect_equal(.fmr_read(cases[[nm]](s), s), ref, tolerance = 1e-12,
                   info = sprintf("%s, k = %d", nm, k))
    }
  }
})

test_that("a constant response keeps R2 NA and an exact fit keeps RMSE 0 and R2 1 at every scale (BUG-209)", {
  D <- .fmr_data()
  one   <- rep(1.3, D$n)
  round <- c(0.1, 0.3 / 3)[rep(1:2, D$n / 2)]
  for (k in .fmr_k) {
    s <- 2^k
    for (y in list(one, round)) {
      m <- spatialkit:::.compute_reg_metrics(y * s, (y + D$y - D$yh) * s, p = 1L)
      expect_true(is.na(m$R2), info = sprintf("k = %d", k))
      expect_true(is.na(m$Adj_R2), info = sprintf("k = %d", k))
      expect_equal(m$RMSE / s, sqrt(mean((D$y - D$yh)^2)), tolerance = 1e-12)
    }
    m <- spatialkit:::.compute_reg_metrics(D$y * s, D$y * s, p = 1L)
    expect_identical(m$RMSE, 0)
    expect_identical(m$MAE, 0)
    expect_identical(m$R2, 1)
    expect_identical(m$Adj_R2, 1)
  }
})

test_that("RMSE is Inf only where it exceeds the largest double (BUG-209)", {
  sg <- c(1, -1)[rep(1:2, 50)]
  # Errors of 0.8e308: each square overflows, the root mean square does not.
  m <- spatialkit:::.compute_reg_metrics(0.4e308 * sg, -0.4e308 * sg)
  expect_equal(m$RMSE, 0.8e308, tolerance = 1e-12)
  expect_equal(m$MAE, 0.8e308, tolerance = 1e-12)
  expect_equal(m$R2, -3, tolerance = 1e-12)
  # Errors of 1.9e308 are not doubles, and neither is their RMSE.
  m <- spatialkit:::.compute_reg_metrics(0.95e308 * sg, -0.95e308 * sg)
  expect_identical(m$RMSE, Inf)
})

test_that("the arithmetic is unchanged where the sums of squares are normal doubles (BUG-209)", {
  D <- .fmr_data()
  for (k in c(-400, -100, 0, 100, 400)) {
    s <- 2^k
    y <- D$y * s; yh <- D$yh * s
    for (b in list(NULL, 0.1 * s)) {
      m <- spatialkit:::.compute_reg_metrics(y, yh, p = 1L, y_train_mean = b)
      rss <- sum((y - yh)^2)
      tss <- sum((y - if (is.null(b)) mean(y) else b)^2)
      expect_identical(m$RMSE, sqrt(rss / D$n))
      expect_identical(m$MAE, mean(abs(y - yh)))
      expect_identical(m$R2, 1 - rss / tss)
      expect_identical(m$Adj_R2, 1 - (1 - (1 - rss / tss)) * (D$n - 1) / (D$n - 2))
    }
  }
})

.fmr_site <- function(s, n = 120) {
  withr::with_seed(1, {
    site <- sf::st_as_sf(
      data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
      coords = c("x", "y"), crs = 32632)
    site$price <- (10 + 0.01 * sf::st_coordinates(site)[, 1] + 2 * site$elev +
                     rnorm(n)) * s
  })
  site
}

.fmr_lm_fit <- function(train_sf) {
  new_spatial_fit(
    subclass = "fmr_lm", engine = stats::lm(price ~ elev, sf::st_drop_geometry(train_sf)),
    formula = price ~ elev, response_var = "price", predictor_vars = "elev",
    data_sf = train_sf)
}
predict.fmr_lm <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data_sf
  as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))
}
registerS3method("predict", "fmr_lm", predict.fmr_lm)

test_that("cv_spatial() scores a response at 2^-540 and 2^530 as at 2^0 (BUG-209)", {
  cv0 <- suppressMessages(cv_spatial(.fmr_site(1), "price", "elev",
                                     fit_fn = .fmr_lm_fit, k = 3, seed = 1))
  expect_true(is.finite(cv0$overall$R2))
  for (k in c(-540, 530)) {
    s <- 2^k
    cv <- suppressMessages(cv_spatial(.fmr_site(s), "price", "elev",
                                      fit_fn = .fmr_lm_fit, folds = cv0$folds))
    expect_equal(cv$overall$RMSE / s, cv0$overall$RMSE, tolerance = 1e-12)
    expect_equal(cv$overall$MAE / s, cv0$overall$MAE, tolerance = 1e-12)
    expect_equal(cv$overall$R2, cv0$overall$R2, tolerance = 1e-12)
    expect_equal(cv$fold_metrics$RMSE / s, cv0$fold_metrics$RMSE, tolerance = 1e-12)
    expect_equal(cv$fold_metrics$R2, cv0$fold_metrics$R2, tolerance = 1e-12)
    expect_equal(cv$fold_metrics$Adj_R2, cv0$fold_metrics$Adj_R2, tolerance = 1e-12)
  }
})

# Each sum is read again only where it underflows or overflows, on its own
# errors or deviations: a far prediction or baseline does not move the
# scale of the sum it is not part of.

test_that("a baseline far from the response leaves RMSE as it is (BUG-209)", {
  D <- .fmr_data()
  ref <- spatialkit:::.compute_reg_metrics(D$y, D$yh, p = 1L)
  for (b in list(1e200, c(1e200, rep(0, D$n - 1)))) {
    m <- spatialkit:::.compute_reg_metrics(D$y, D$yh, p = 1L, y_train_mean = b)
    expect_identical(m$RMSE, ref$RMSE)
    expect_identical(m$R2, 1)
  }
  # A response at 2^-540 against a baseline of 1.
  s <- 2^-540
  m <- spatialkit:::.compute_reg_metrics(D$y * s, D$yh * s, p = 1L, y_train_mean = 1)
  expect_equal(m$RMSE / s, ref$RMSE, tolerance = 1e-12)
  expect_identical(m$R2, 1)
})

test_that("one prediction far from the response gives R2 -Inf and its RMSE (BUG-209)", {
  D <- .fmr_data()
  for (v in c(1e160, 1e200, 1e300)) {
    yh <- D$yh; yh[1] <- v
    m <- spatialkit:::.compute_reg_metrics(D$y, yh, p = 1L)
    expect_identical(m$R2, -Inf)
    expect_equal(m$RMSE, sqrt(sum((D$y / v - yh / v)^2) / D$n) * v,
                 tolerance = 1e-12)
  }
})

test_that("a sum is read on its own terms, not on the largest value (BUG-209)", {
  # The error of 1e-200 beside a value of 1e300 fitted exactly: divided by
  # the power of two of 1e300 its square is 0.
  m <- spatialkit:::.compute_reg_metrics(c(1e300, 1e-200), c(1e300, 2e-200))
  expect_equal(m$RMSE * 1e200, 1 / sqrt(2), tolerance = 1e-12)
  # The same for the deviations about a baseline: TSS 1e-400, RSS a quarter.
  m <- spatialkit:::.compute_reg_metrics(c(1e300, 1e-200), c(1e300, 1.5e-200),
                                         y_train_mean = c(1e300, 2e-200))
  expect_equal(m$R2, 0.75, tolerance = 1e-12)
  expect_equal(m$RMSE * 1e200, 0.5 / sqrt(2), tolerance = 1e-12)
})

test_that("the largest doubles are read without overflow (BUG-209)", {
  x <- .Machine$double.xmax
  y <- c(x, x / 2, x / 4, -x / 8); yh <- c(x / 2, x / 2, x / 4, -x / 4)
  m <- spatialkit:::.compute_reg_metrics(y, yh, p = 1L)
  r <- spatialkit:::.compute_reg_metrics(y / x, yh / x, p = 1L)
  expect_equal(m$RMSE / x, r$RMSE, tolerance = 1e-12)
  expect_equal(m$R2, r$R2, tolerance = 1e-12)
  # A response whose mean() sums past the largest double where a long
  # double is no wider than a double (it gives -Inf for this one).
  y <- c(x, x, x / 2, -x)
  m <- spatialkit:::.compute_reg_metrics(y, 0.75 * y)
  r <- spatialkit:::.compute_reg_metrics(y / x, 0.75 * y / x)
  expect_equal(m$RMSE / x, r$RMSE, tolerance = 1e-12)
  expect_equal(m$R2, r$R2, tolerance = 1e-12)
  # Errors that overflow as differences, 0.75 x less -0.75 x, with an RMSE
  # of 0.75 x and R2 1 - 4.
  y <- c(0.75 * x, -0.75 * x, rep(0, 6))
  m <- spatialkit:::.compute_reg_metrics(y, -y)
  expect_equal(m$RMSE, 0.75 * x, tolerance = 1e-12)
  expect_equal(m$R2, -3, tolerance = 1e-12)
})

test_that("an exact fit of subnormal values keeps RMSE 0 and R2 1 (BUG-209)", {
  y <- c(1, 2, 3, 5) * 2^-1074
  m <- spatialkit:::.compute_reg_metrics(y, y)
  expect_identical(m$RMSE, 0)
  expect_identical(m$R2, 1)
})

# A fit that predicts the training median plus 2 elev: one training
# response recorded as 1e300 does not move it, but moves the training mean
# model_metrics(newdata) measures R2 against.
.fmr_med_fit <- function(train_sf) {
  new_spatial_fit(
    subclass = "fmr_med", engine = list(m = stats::median(train_sf$price)),
    formula = price ~ elev, response_var = "price", predictor_vars = "elev",
    data_sf = train_sf)
}
predict.fmr_med <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) newdata <- object$data_sf
  object$engine$m + 2 * newdata$elev
}
registerS3method("predict", "fmr_med", predict.fmr_med)

test_that("model_metrics(newdata) beside a training response of 1e300 keeps its RMSE (BUG-209)", {
  site <- .fmr_site(1)
  tr <- site[1:80, ]; te <- site[81:120, ]
  tr$price[1] <- 1e300
  f <- .fmr_med_fit(tr)
  yh <- predict(f, te)
  m <- model_metrics(f, newdata = te)
  expect_equal(m$RMSE, sqrt(mean((te$price - yh)^2)), tolerance = 1e-12)
  expect_true(m$RMSE > 0.5 && m$RMSE < 10)
  expect_identical(m$R2, 1)
})
