# tests/testthat/test-fix-variogram-structure.R
# ---------------------------------------------------------------------------
# estimate_sac_range() tests for spatial structure before it returns a range.
#
# Nothing did.  A variogram model can be fitted to anything, and on white
# noise (300 points on a 1000 m square) a finite range came back for 24 of 60
# layers, 7 of them with a nugget under half the sill, which removes the
# penalty of a resolution criterion that reads it.
# ---------------------------------------------------------------------------

fvs_white <- function(i, n = 300L) {
  set.seed(70000 + n + i)
  sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                          z = rnorm(n)), coords = c("x", "y"), crs = 32632)
}
# An exponential field: effective range R, nugget share nug.
fvs_field <- function(seed, n = 300L, R = 300, nug = 0.2) {
  set.seed(seed)
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  z <- drop(crossprod(chol((1 - nug) * exp(-3 * as.matrix(stats::dist(xy)) / R) +
                             diag(nug, n)), rnorm(n)))
  sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2], z = z),
               coords = c("x", "y"), crs = 32632)
}
NO_STRUCTURE <- "no spatial structure detected at the nearest-neighbour scale"

test_that("white noise gets no range: at most 2 of 20 layers come back finite", {
  skip_if_not_installed("gstat")
  res <- lapply(1:20, function(i) suppressWarnings(estimate_sac_range(fvs_white(i), "z")))
  finite <- vapply(res, function(r) is.finite(r), logical(1))
  expect_lte(sum(finite), 2L)
  refused <- vapply(res, function(r) identical(attr(r, "rejected_reason"), NO_STRUCTURE),
                    logical(1))
  # Most are refused by the test itself, with its evidence attached.
  expect_gte(sum(refused), 17L)
  for (r in res[refused]) {
    expect_true(is.na(r))
    expect_s3_class(r, "sac_range")
    expect_gte(attr(r, "structure_p"), 0.05)
    expect_true(is.na(attr(r, "rejected_range")))
    expect_s3_class(attr(r, "variogram"), "data.frame")
  }
  # A layer the test lets through is one it found structure in.
  for (r in res[!refused]) expect_lt(attr(r, "structure_p"), 0.05)
})

test_that("the refusal carries a pure-nugget model, not the model that was fitted", {
  skip_if_not_installed("gstat")
  pts <- fvs_white(1)
  lines <- capture_spatialkit_log(r <- estimate_sac_range(pts, "z"))
  expect_identical(attr(r, "rejected_reason"), NO_STRUCTURE)
  vm <- attr(r, "variogram_model")
  expect_identical(as.character(vm$model), "Nug")
  expect_equal(vm$psill, stats::var(pts$z))
  expect_true(spatialkit:::.vgm_is_pure_nugget(vm))
  expect_null(spatialkit:::.vgm_correlation_fn(vm))
  # The nugget is the variance of the response, and it is not marked as a
  # bound: it is the nugget of the model the test supports.
  expect_equal(sac_nugget(r), stats::var(pts$z))
  expect_false(attr(r, "nugget_at_bound"))
  expect_identical(attr(r, "structures"), 0L)
  # The log line says what was measured and how to override it.
  expect_true(log_has(lines, "no spatial structure detected at the nearest-neighbour scale"))
  expect_true(log_has(lines, "Moran's I of the response on the [0-9]+ nearest neighbours is"))
  expect_true(log_has(lines, "not below structure_alpha = 0.05"))
  expect_true(log_has(lines, "set structure_alpha = 1 to fit a range regardless"))
  # It prints as NA and draws, with the reason as the subtitle.
  expect_output(print(r), "^NA")
  skip_if_not_installed("ggplot2")
  p <- plot(r)
  expect_s3_class(p, "ggplot")
  expect_match(gsub("\n", " ", p$labels$subtitle),
               "No effective range: no spatial structure detected at the nearest-neighbour scale")
})

test_that("structure_alpha = 1 switches the refusal off and keeps the test attached", {
  skip_if_not_installed("gstat")
  pts <- fvs_white(1)
  r0 <- suppressWarnings(estimate_sac_range(pts, "z"))
  r1 <- suppressWarnings(estimate_sac_range(pts, "z", structure_alpha = 1))
  expect_identical(attr(r0, "rejected_reason"), NO_STRUCTURE)
  expect_false(identical(attr(r1, "rejected_reason"), NO_STRUCTURE))
  expect_identical(attr(r1, "structure_p"), attr(r0, "structure_p"))
  expect_identical(attr(r1, "structure_test"), attr(r0, "structure_test"))
  expect_identical(attr(r1, "morans_i"), attr(r0, "morans_i"))
  # A stricter or laxer level moves the line and nothing else.
  p <- attr(r0, "structure_p")
  lax <- suppressWarnings(estimate_sac_range(pts, "z", structure_alpha = min(1, p * 1.01)))
  expect_false(identical(attr(lax, "rejected_reason"), NO_STRUCTURE))
})

test_that("a field with structure passes, and the test is attached to the range", {
  skip_if_not_installed("gstat")
  for (s in 1:8) {
    r <- suppressWarnings(estimate_sac_range(fvs_field(4100 + s), "z"))
    expect_lt(attr(r, "structure_p"), 0.05)
    expect_false(identical(attr(r, "rejected_reason"), NO_STRUCTURE))
    expect_gt(attr(r, "morans_i"), 0)
  }
  r <- estimate_sac_range(fvs_field(4101), "z")
  tab <- attr(r, "structure_test")
  expect_s3_class(tab, "data.frame")
  expect_identical(names(tab), c("k", "morans_i", "z", "p"))
  expect_identical(tab$k, c(2L, 4L, 8L, 16L))
  # The adjusted p-value is the smallest of the four, times four.
  expect_equal(attr(r, "structure_p"), min(1, 4 * min(tab$p)))
  expect_equal(attr(r, "morans_i"), tab$morans_i[tab$k == 8L])
  expect_equal(tab$p, stats::pnorm(tab$z, lower.tail = FALSE))
})

test_that("the statistic is Moran's I on k nearest neighbours with its exact null moments", {
  skip_if_not_installed("gstat")
  set.seed(31); n <- 150
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000)); y <- rnorm(n)
  st <- spatialkit:::.sac_structure_test(xy, y)
  # Against the definition, from the full distance matrix.
  D <- as.matrix(stats::dist(xy)); diag(D) <- Inf
  e <- y - mean(y)
  for (j in seq_along(st$table$k)) {
    k <- st$table$k[j]
    lagged <- vapply(seq_len(n), function(i) mean(e[order(D[i, ])[seq_len(k)]]), numeric(1))
    expect_equal(st$table$morans_i[j], sum(e * lagged) / sum(e^2), tolerance = 1e-10)
  }
  # With an intercept alone the null mean is the classical -1 / (n - 1).
  W <- spatialkit:::.build_knn_weights(xy, k = 8L)
  mom <- spatialkit:::.morans_residual_moments(W, matrix(1, n, 1), sum(W), inherits(W, "Matrix"))
  expect_equal(mom$EI, -1 / (n - 1), tolerance = 1e-10)
  expect_equal(st$table$z[st$table$k == 8L],
               (st$table$morans_i[st$table$k == 8L] - mom$EI) / sqrt(mom$VI), tolerance = 1e-10)
  # No random numbers: the caller's stream is where it was, and two calls agree.
  set.seed(5); before <- runif(2)
  set.seed(5); st2 <- spatialkit:::.sac_structure_test(xy, y); after <- runif(2)
  expect_identical(after, before)
  expect_identical(st2, st)
  # With a trend's design it is the statistic of the OLS residuals, under
  # the moments that fit them.
  X <- cbind(1, xy[, 1] / 1000, (xy[, 2] / 1000)^2)
  stx <- spatialkit:::.sac_structure_test(xy, y, X)
  ex <- stats::residuals(stats::lm(y ~ X - 1))
  expect_equal(stx$table$morans_i[stx$table$k == 8L],
               (n / sum(W)) * sum(ex * as.numeric(W %*% ex)) / sum(ex^2), tolerance = 1e-10)
  momx <- spatialkit:::.morans_residual_moments(W, X, sum(W), inherits(W, "Matrix"))
  expect_lt(momx$EI, mom$EI)          # residuals are pushed below -1/(n-1)
  expect_equal(stx$table$z[stx$table$k == 8L],
               (stx$table$morans_i[stx$table$k == 8L] - momx$EI) / sqrt(momx$VI),
               tolerance = 1e-10)
})

test_that("repeat visits and a regular grid are tested without error, in any row order", {
  skip_if_not_installed("gstat")
  set.seed(32); n <- 120
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  xy <- rbind(xy, xy[1:30, ]); y <- rnorm(nrow(xy))
  a <- spatialkit:::.sac_structure_test(xy, y)
  expect_true(is.finite(a$p) && a$p >= 0 && a$p <= 1)
  g <- as.matrix(expand.grid(x = 1:12, y = 1:12)); storage.mode(g) <- "double"
  set.seed(33); yg <- rnorm(nrow(g))
  b <- spatialkit:::.sac_structure_test(g, yg)
  expect_true(is.finite(b$p))
  # Ties at the edge of a neighbourhood are shared out, so the order of the
  # rows does not decide which tied neighbour counts.
  perm <- sample(nrow(g))
  b2 <- spatialkit:::.sac_structure_test(g[perm, ], yg[perm])
  expect_equal(b2$table$morans_i, b$table$morans_i, tolerance = 1e-10)
  expect_equal(b2$p, b$p, tolerance = 1e-8)
  # Too few points for any neighbourhood: no test, not an error.
  expect_null(spatialkit:::.sac_structure_test(xy[1:3, ], y[1:3]))
})

test_that("residuals on predictors are judged as residuals", {
  skip_if_not_installed("gstat")
  # White noise regressed on smooth predictors: no structure is left, and the
  # refusal says it is about the residuals.
  set.seed(34); n <- 200
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  d <- data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2], a = xy[, 1] / 1000,
                  b = sin(xy[, 2] / 300), z = rnorm(n))
  pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
  lines <- capture_spatialkit_log(r <- estimate_sac_range(pts, "z", c("a", "b")))
  expect_identical(attr(r, "rejected_reason"), NO_STRUCTURE)
  expect_true(attr(r, "detrended"))
  expect_true(log_has(lines, "Moran's I of the residuals on predictor_vars"))
  # The statistic attached is the one for the OLS residuals on that design.
  X <- stats::model.matrix(z ~ a + b, d)
  st <- spatialkit:::.sac_structure_test(xy + matrix(c(5e5, 5e6), n, 2, byrow = TRUE), d$z, X)
  expect_equal(attr(r, "structure_p"), st$p, tolerance = 1e-8)
  # A response whose structure the predictors do not explain keeps its range.
  f <- fvs_field(4102)
  set.seed(35); f$w <- rnorm(nrow(f))
  rf <- suppressWarnings(estimate_sac_range(f, "z", "w"))
  expect_lt(attr(rf, "structure_p"), 0.05)
  expect_true(is.finite(rf))
})

test_that("the size of the test is near its level on white noise", {
  skip_if_not_installed("gstat")
  # 300 layers of 100 points: the share with an adjusted p below 0.05.  Over
  # 1000 draws it was 0.041 at n = 300 and 0.053 at n = 100.
  rej <- vapply(1:300, function(i) {
    set.seed(880000 + i)
    xy <- cbind(runif(100, 0, 1000), runif(100, 0, 1000))
    spatialkit:::.sac_structure_test(xy, rnorm(100))$p < 0.05
  }, logical(1))
  expect_gte(mean(rej), 0.02)
  expect_lte(mean(rej), 0.10)
})
