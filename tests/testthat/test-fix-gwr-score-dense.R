# tests/testthat/test-fix-gwr-score-dense.R
# ---------------------------------------------------------------------------
# The bandwidth criterion under a kernel that weights every row.
#
# A Gaussian or an exponential kernel gives every observation a weight at
# every regression point, so scoring a bandwidth by one QR per window is n
# factorisations of n rows.  A default-bandwidth fit took 10.0 s at n = 1000
# and 32.0 s at n = 2000 (fixed Gaussian), against 1.8 s and 6.2 s for
# GWmodel::bw.gwr().
#
# The score of a window that is far from rank deficient now comes from its
# centred cross-products, many windows to one matrix product
# (.gwr_score_dense()); the QR loop keeps the windows that are not.  These
# tests hold the two routes to the same numbers and the same ranks.
# ---------------------------------------------------------------------------

skip_if_not_installed("sp")

# The score by the QR loop alone: the batched path answers no window.
.sd_qr_only <- function(...) {
  testthat::local_mocked_bindings(
    .gwr_score_dense = function(X, ...) list(done = logical(nrow(X))))
  .gwr_local_solve(...)
}

.sd_layer <- function(n, p, seed, reps = 1L) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  id <- rep(seq_len(n), reps)
  P <- matrix(rnorm(n * max(p, 1L)), n, max(p, 1L))[id, , drop = FALSE]
  if (p >= 2L) P[, 2L] <- 0.7 * P[, 1L] + 0.5 * P[, 2L]
  # With repeat locations, one predictor differs between the rows of a site.
  if (reps > 1L && p >= 2L) P[, 2L] <- P[, 2L] + rnorm(nrow(P), 0, 0.3)
  X <- cbind(Intercept = 1, P)[, seq_len(p + 1L), drop = FALSE]
  if (p) colnames(X)[-1L] <- paste0("v", seq_len(p))
  list(X = X, xy = cbind(x, y)[id, , drop = FALSE],
       y = 2 + (1 + x[id] / 1000) * P[, 1L] + rnorm(n, 0, 0.3)[id] +
         rnorm(length(id), 0, 0.5))
}

.sd_rel <- function(a, b) {
  ok <- is.finite(a) | is.finite(b)
  if (!any(ok)) return(0)
  max(abs(a[ok] - b[ok]) / pmax(1, abs(a[ok])))
}


test_that("the batched score equals the QR loop's, window for window", {
  cases <- list(
    list(n = 250, p = 2, k = "gaussian",    ad = FALSE, bw = c(8, 15, 40, 150, 800)),
    list(n = 250, p = 2, k = "gaussian",    ad = TRUE,  bw = c(4, 10, 40, 250, 400)),
    list(n = 250, p = 4, k = "exponential", ad = TRUE,  bw = c(6, 18, 100)),
    list(n = 250, p = 1, k = "exponential", ad = FALSE, bw = c(2, 5, 60, 500)),
    list(n = 250, p = 0, k = "gaussian",    ad = FALSE, bw = c(15, 150)))
  answered <- 0L; left <- 0L
  for (cs in cases) {
    g <- .sd_layer(cs$n, cs$p, 1)
    for (bw in cs$bw) {
      lab <- sprintf("%s %s p=%d bw=%g", cs$k, if (cs$ad) "adaptive" else "fixed", cs$p, bw)
      qr <- .sd_qr_only(g$X, g$y, g$xy, bw, cs$k, cs$ad, what = "score")
      got <- .gwr_local_solve(g$X, g$y, g$xy, bw, cs$k, cs$ad, what = "score")
      done <- .gwr_score_dense(g$X, g$y, g$xy, bw, cs$k, cs$ad)$done
      answered <- answered + sum(done); left <- left + sum(!done)
      # Everything that is a count or a decision is identical.
      expect_identical(got$rank, qr$rank, label = lab)
      expect_identical(got$n_window, qr$n_window, label = lab)
      expect_identical(got$width, qr$width, label = lab)
      expect_identical(got$n_estimable, qr$n_estimable, label = lab)
      expect_identical(is.na(got$fitted), is.na(qr$fitted), label = lab)
      # The batched path answers only windows the QR finds of full rank.
      expect_true(all(qr$rank[done] == ncol(g$X)), label = lab)
      # The numbers agree to rounding error.
      expect_lt(.sd_rel(qr$fitted, got$fitted), 1e-9, label = lab)
      expect_lt(.sd_rel(qr$leverage, got$leverage), 1e-9, label = lab)
    }
  }
  # Both routes were taken: most windows batched, the nearly deficient ones
  # (the smallest fixed bandwidths) left to the QR.
  expect_gt(answered, 4000L)
  expect_gt(left, 100L)
})

test_that("the batched path leaves a window to the QR when it is close to rank deficient", {
  # A covariate that is one value on each side of x = 500: constant in the
  # windows away from the boundary under a fixed Gaussian of 25, where the
  # weights of the other side underflow against 1e-7.
  set.seed(5); n <- 300
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n)
  region <- ifelse(x > 500, 1.3, 0.7)
  z <- 1 + 2 * a + 3 * region + rnorm(n, 0, 0.3)
  X <- cbind(Intercept = 1, a = a, region = region); xy <- cbind(x, y)
  got <- .gwr_local_solve(X, z, xy, 25, "gaussian", FALSE, what = "score")
  qr  <- .sd_qr_only(X, z, xy, 25, "gaussian", FALSE, what = "score")
  done <- .gwr_score_dense(X, z, xy, 25, "gaussian", FALSE)$done
  deficient <- qr$rank < 3L
  expect_gt(sum(deficient), 50L)
  expect_false(any(done & deficient))
  expect_identical(got$rank, qr$rank)
  expect_identical(got$n_estimable, qr$n_estimable)
  expect_lt(.sd_rel(qr$fitted, got$fitted), 1e-9)
  expect_lt(.sd_rel(qr$leverage, got$leverage), 1e-9)
  # A predictor far from zero: centred on its mean before any product, so the
  # batched path loses nothing to the origin.
  Xs <- X; Xs[, "a"] <- Xs[, "a"] + 1e7
  got0 <- .gwr_local_solve(X, z, xy, 200, "gaussian", FALSE, what = "score")
  got7 <- .gwr_local_solve(Xs, z, xy, 200, "gaussian", FALSE, what = "score")
  expect_true(all(.gwr_score_dense(Xs, z, xy, 200, "gaussian", FALSE)$done))
  expect_lt(.sd_rel(got0$fitted, got7$fitted), 1e-7)
  # A zero-width adaptive kernel (the bandwidth's rows share one location)
  # has undefined weights: no window there, by either route.
  xy2 <- xy; xy2[1:6, ] <- xy2[rep(1, 6), ]
  got2 <- .gwr_local_solve(X, z, xy2, 5, "gaussian", TRUE, what = "score")
  expect_true(all(is.na(got2$fitted[1:6])))
  expect_identical(got2$n_window[1:6], rep(0L, 6))
  expect_false(any(.gwr_score_dense(X, z, xy2, 5, "gaussian", TRUE)$done[1:6]))
})

test_that("on a layer whose locations repeat the batched path gives the same leave-out residuals", {
  for (cs in list(list(p = 2, k = "gaussian",    ad = TRUE,  bw = c(7, 13, 40, 200), reps = 3L),
                  list(p = 2, k = "exponential", ad = FALSE, bw = c(10, 60, 400),     reps = 3L),
                  list(p = 0, k = "exponential", ad = FALSE, bw = c(10, 60),          reps = 2L))) {
    g <- .sd_layer(90, cs$p, 3, cs$reps)
    loc <- .gwr_locations(g$xy)
    for (bw in cs$bw) {
      lab <- sprintf("%s p=%d bw=%g", cs$k, cs$p, bw)
      qr  <- .sd_qr_only(g$X, g$y, g$xy, bw, cs$k, cs$ad, what = "score", loc = loc)
      got <- .gwr_local_solve(g$X, g$y, g$xy, bw, cs$k, cs$ad, what = "score", loc = loc)
      expect_identical(lengths(got$hat_loc), lengths(qr$hat_loc), label = lab)
      r_qr  <- .gwr_lolo_residuals(g$y - qr$fitted, qr$leverage, qr$hat_loc, loc)
      r_got <- .gwr_lolo_residuals(g$y - got$fitted, got$leverage, got$hat_loc, loc)
      expect_identical(is.na(r_got), is.na(r_qr), label = lab)
      expect_lt(.sd_rel(r_qr, r_got), 1e-7, label = lab)
    }
  }
})

test_that("the criterion and the selected bandwidth do not depend on the route", {
  g <- .sd_layer(200, 2, 7)
  for (cs in list(list(k = "gaussian", ad = FALSE), list(k = "exponential", ad = TRUE))) {
    a <- .gwr_bw_search(g$X, g$y, g$xy, cs$k, cs$ad)
    b <- local({
      testthat::local_mocked_bindings(
        .gwr_score_dense = function(X, ...) list(done = logical(nrow(X))))
      .gwr_bw_search(g$X, g$y, g$xy, cs$k, cs$ad)
    })
    expect_identical(a$bandwidth, b$bandwidth, label = cs$k)
    expect_identical(a$curve$bandwidth, b$curve$bandwidth, label = cs$k)
    expect_equal(a$curve$criterion, b$curve$criterion, tolerance = 1e-9, label = cs$k)
    expect_equal(a$curve$tr_S, b$curve$tr_S, tolerance = 1e-9, label = cs$k)
    expect_identical(a$curve$n_deficient, b$curve$n_deficient, label = cs$k)
  }
  # The fit at the selected bandwidth, which the QR makes, has the AICc the
  # search scored.
  d <- sf::st_as_sf(data.frame(x = 5e5 + g$xy[, 1], y = 5e6 + g$xy[, 2],
                               v1 = g$X[, 2], v2 = g$X[, 3], z = g$y),
                    coords = c("x", "y"), crs = 32632)
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("v1", "v2"), adaptive = FALSE, kernel = "gaussian")))
  expect_equal(fit$info$bandwidth_search$criterion, fit$info$AICc, tolerance = 1e-9)
})

test_that("a kernel with an edge is scored by the QR loop alone", {
  g <- .sd_layer(120, 2, 9)
  called <- 0L
  testthat::local_mocked_bindings(
    .gwr_score_dense = function(X, ...) { called <<- called + 1L; list(done = logical(nrow(X))) })
  for (k in c("bisquare", "tricube", "boxcar"))
    .gwr_local_solve(g$X, g$y, g$xy, 30, k, TRUE, what = "score")
  expect_identical(called, 0L)
  .gwr_local_solve(g$X, g$y, g$xy, 30, "gaussian", TRUE, what = "score")
  expect_identical(called, 1L)
  # The fit and the predictions are always the QR's.
  .gwr_local_solve(g$X, g$y, g$xy, 30, "gaussian", TRUE, what = "fit")
  .gwr_local_solve(g$X, g$y, g$xy, 30, "gaussian", TRUE, rp = g$xy[1:5, ],
                   X_rp = g$X[1:5, ], what = "predict")
  expect_identical(called, 1L)
})

test_that("scoring a bandwidth under a Gaussian kernel takes a fraction of the QR loop's time", {
  skip_on_cran()
  g <- .sd_layer(1500, 2, 11)
  best <- function(f) min(vapply(1:3, function(i) system.time(f())[["elapsed"]], numeric(1)))
  t_new <- best(function() .gwr_local_solve(g$X, g$y, g$xy, 80, "gaussian", FALSE, what = "score"))
  t_qr  <- best(function() .sd_qr_only(g$X, g$y, g$xy, 80, "gaussian", FALSE, what = "score"))
  # Measured 0.07 to 0.21 s against 0.33 to 0.50 s over three runs each; the
  # best of three is compared, loosely, for a loaded machine.
  expect_lt(t_new, 0.8 * t_qr)
})
