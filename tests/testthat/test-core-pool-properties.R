# tests/testthat/test-core-pool-properties.R
# ---------------------------------------------------------------------------
# What the engine of R/pool.R promises on any input, checked without the
# reference functions: the leave-one-out residuals against a loop that takes
# each row out, the behaviour under a change of scale, of row order and of
# region numbering, the outputs that are set when a column has no spread,
# the order of the degenerate cases, the neighbour searches on tied and
# coincident points, the rule for a response without structure, the bound on
# exact searches, and what `want`, `geo` and `cache` leave uncomputed.
#
# None of it needs FNN.  Each test runs on the exact searches of base R, and
# on FNN's kd-tree as well where FNN is installed; the tests that hold the
# two paths to identical() are the ones that skip without it.
#
# The package's internal functions are called by name: testthat runs a
# package's tests in an environment that inherits from its namespace.
# ---------------------------------------------------------------------------

# ---- Layouts -----------------------------------------------------------------

.pt_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# Regions of a k-means fit at `K` centres: each row in the region of its
# nearest centre, the regions that hold a row numbered 1 to K.
.pt_regions <- function(xy, K, seed = 1L) {
  .pt_seed(seed)
  km <- suppressWarnings(stats::kmeans(xy, K, nstart = 3, iter.max = 100))
  g  <- .nearest_centre(xy[, 1], xy[, 2], km$centers[, 1], km$centers[, 2])
  match(g, sort(unique(g)))
}

.pt_store <- new.env(parent = emptyenv())
.pt_keep <- function(key, expr) {
  if (is.null(.pt_store[[key]])) .pt_store[[key]] <- expr
  .pt_store[[key]]
}

.pt_layout <- function(i) .pt_keep(paste0("layout", i), do.call(pool_layout, POOL_LAYOUTS[[i]]))

# The lattice: 60 x 50 locations at a spacing of 0.7 from (0.1, 0.3).  None
# of its coordinates is exact in binary, so sums of them taken in different
# orders differ in their last bits, and every location has neighbours at
# exactly equal distances.  `centres` regions: 667 gives about 4.5 rows a
# region and none of one row; 1,500 gives 2 rows a region and many of one.
.pt_lattice <- function(centres = 667L) .pt_keep(paste0("lattice", centres), {
  xy <- as.matrix(expand.grid(x = seq(0.1, by = 0.7, length.out = 60),
                              y = seq(0.3, by = 0.7, length.out = 50)))
  dimnames(xy) <- NULL
  y <- local({
    .pt_seed(31)
    sin(xy[, 1] / 9) + cos(xy[, 2] / 7) + stats::rnorm(nrow(xy), 0, 0.5)
  })
  gi <- .pt_regions(xy, centres, 3L)
  list(xy = xy, y = y, gi = gi, K = max(gi))
})

# 400 rows in 90 regions, for the columns without spread.
.pt_small <- function() .pt_keep("small", {
  xy <- local({ .pt_seed(7); cbind(stats::runif(400, 0, 1e4), stats::runif(400, 0, 1e4)) })
  gi <- .pt_regions(xy, 90L, 1L)
  list(xy = xy, gi = gi, K = max(gi))
})

# `K` regions of 500 rows each, every row within 900 of its region's centre.
# The sum of 500 values that are not exact in binary is off by a rounding
# error, which is what the regions of the layouts above are too small to show.
.pt_large <- function(K) .pt_keep(paste0("large", K), {
  gi  <- rep(seq_len(K), each = 500L)
  cen <- local({ .pt_seed(111); cbind(stats::runif(K, 0, 1e5), stats::runif(K, 0, 1e5)) })
  xy  <- local({
    .pt_seed(112)
    cen[gi, ] + cbind(stats::runif(500L * K, -900, 900), stats::runif(500L * K, -900, 900))
  })
  list(xy = xy, gi = gi, K = K)
})

# The test for no spread inside the regions on the sum of squares the fit
# uses: about the regions' means, each mean its sum over its count.
.pt_written_rule <- function(y, gi) {
  m <- as.numeric(rowsum(y, gi)) / tabulate(gi)
  .no_spread(sum((y - m[gi])^2), y)
}

# The searches a test runs on: the exact one, and FNN's where it is installed.
.pt_paths <- function() c(FALSE, if (requireNamespace("FNN", quietly = TRUE)) TRUE)
.pt_path_tag <- function(f) if (f) "with FNN" else "exact search"

# The largest absolute difference between two vectors of the same length.
.pt_gap <- function(a, b) {
  stopifnot(length(a) == length(b))
  if (!length(a)) return(0)
  max(abs(a - b))
}


# ---- Definitions written again, row by row -----------------------------------

# Each region's neighbours by the rule alone: every other region by squared
# distance between centroids, a tie to the lower first coordinate, then the
# lower second, then the lower number.
.pt_neighbours <- function(cen, k = 8L) {
  K <- nrow(cen); k <- min(as.integer(k), K - 1L)
  out <- matrix(0L, K, k)
  for (j in seq_len(K)) {
    d2 <- (cen[, 1] - cen[j, 1])^2 + (cen[, 2] - cen[j, 2])^2
    o  <- order(d2, cen[, 1], cen[, 2], seq_len(K))
    out[j, ] <- o[o != j][seq_len(k)]
  }
  out
}

# A nearest other row by the rule alone.  A row that shares its location
# takes the first other row there, rows ordered by region and then as given;
# a row alone at its location takes the first row of the nearest other
# location, a tie to the lower first coordinate, then the lower second.
.pt_nearest_other <- function(xy, gi, rows = seq_len(nrow(xy))) {
  n <- nrow(xy)
  vapply(rows, function(i) {
    here <- which(xy[, 1] == xy[i, 1] & xy[, 2] == xy[i, 2])
    here <- here[order(gi[here], here)]
    if (length(here) > 1L) return(here[here != i][1L])
    d2 <- (xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2
    d2[i] <- Inf
    near <- which(d2 == min(d2))
    best <- near[order(xy[near, 1], xy[near, 2])][1L]
    at <- which(xy[, 1] == xy[best, 1] & xy[, 2] == xy[best, 2])
    at[order(gi[at], at)][1L]
  }, integer(1))
}

# Take each row out, compute its region's count and mean again, hold sigma2,
# tau2 and the targets, and predict it.
.pt_leave_out <- function(y, gi, xy, fit) {
  N <- length(y)
  rows_of <- split(seq_len(N), gi)
  eL <- eP <- numeric(N)
  for (i in seq_len(N)) {
    j <- gi[i]
    rest <- setdiff(rows_of[[j]], i)
    if (length(rest)) {
      mean_rest <- sum(y[rest]) / length(rest)
      v <- fit$sigma2 / length(rest)
      w <- if (fit$tau2 + v > 0) fit$tau2 / (fit$tau2 + v) else 1
      eL[i] <- y[i] - (w * mean_rest + (1 - w) * fit$target[j])
      eP[i] <- y[i] - mean_rest
    } else {
      eL[i] <- y[i] - fit$target[j]
      d2 <- (xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2
      d2[i] <- Inf
      eP[i] <- y[i] - fit$m[gi[which.min(d2)]]
    }
  }
  list(eL = eL, eP = eP)
}


# ---- 7. Leave-one-out by brute force ----------------------------------------

test_that("7. taking each row out and predicting it gives the two residual vectors", {
  for (i in c(2L, 3L, 7L)) for (f in .pt_paths()) {
    L <- .pt_layout(i); tag <- sprintf("layout %d, %s", i, .pt_path_tag(f))
    y <- L$y; gi <- L$gi; K <- L$K; xy <- L$xy; N <- length(y)
    A <- .pool_fit(y, gi, K, xy, use_fnn = f)

    # What is held while a row is out, each from its definition.
    n <- as.integer(table(gi))
    m <- vapply(seq_len(K), function(j) mean(y[gi == j]), numeric(1))
    expect_identical(A$n, n, info = tag)
    expect_lt(.pt_gap(A$m, m), 1e-12)
    expect_equal(A$sigma2, sum((y - m[gi])^2) / (N - K), tolerance = 1e-12, info = tag)
    cen <- t(vapply(seq_len(K), function(j) colMeans(xy[gi == j, , drop = FALSE]), numeric(2)))
    expect_equal(unname(A$cen), cen, tolerance = 1e-12, info = tag)
    expect_identical(A$nb, .pt_neighbours(unname(A$cen)), info = tag)
    expect_false(any(A$nb == seq_len(K)), info = tag)
    target <- vapply(seq_len(K), function(j) mean(y[gi %in% A$nb[j, ]]), numeric(1))
    expect_lt(.pt_gap(A$target, target), 1e-12)
    Q <- sum(n * (m - target)^2)
    expect_equal(A$Q, Q, tolerance = 1e-12, info = tag)
    expect_equal(A$tau2, max(0, (Q - K * A$sigma2) / N), tolerance = 1e-10, info = tag)
    w <- A$tau2 / (A$tau2 + A$sigma2 / n)
    expect_lt(.pt_gap(A$w, w), 1e-12)
    expect_lt(.pt_gap(A$value, w * m + (1 - w) * target), 1e-12)

    # The residuals.
    B <- .pt_leave_out(y, gi, xy, A)
    expect_lt(.pt_gap(A$resid_pooled, B$eL), 1e-12)
    expect_lt(.pt_gap(A$resid_mean, B$eP), 1e-12)
    expect_equal(A$loo_pooled, mean(B$eL^2), tolerance = 1e-12, info = tag)
    expect_equal(A$loo_mean, mean(B$eP^2), tolerance = 1e-12, info = tag)
    expect_equal(A$loo_pooled_se, stats::sd(B$eL^2) / sqrt(N), tolerance = 1e-12, info = tag)
    expect_equal(A$loo_mean_se, stats::sd(B$eP^2) / sqrt(N), tolerance = 1e-12, info = tag)

    # One mean for all rows: each row predicted by the mean of the others.
    e1 <- vapply(seq_len(N), function(r) y[r] - mean(y[-r]), numeric(1))
    expect_equal(A$loo_one, mean(e1^2), tolerance = 1e-12, info = tag)
    expect_equal(A$loo_one_se, stats::sd(e1^2) / sqrt(N), tolerance = 1e-12, info = tag)
    expect_identical(A$loo_one, mean(((y - mean(y)) * N / (N - 1))^2), info = tag)

    # The F test of the region means against one mean, from the sum of
    # squares between the regions.
    between <- sum(n * (m - mean(y))^2)
    p <- stats::pf((between / (K - 1)) / (sum((y - m[gi])^2) / (N - K)), K - 1, N - K,
                   lower.tail = FALSE)
    expect_gt(p, 0)
    expect_lt(abs(A$f_p / p - 1), 1e-6)
  }
  # Layouts 3 and 7 hold regions of one row, the residuals with their own rule.
  expect_gt(sum(tabulate(.pt_layout(3L)$gi) == 1L), 50L)
  expect_gt(sum(tabulate(.pt_layout(7L)$gi) == 1L), 500L)
})

test_that("7. the F test's p-value is the one-way analysis of variance's", {
  for (i in 8:13) {
    L <- .pt_layout(i)
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = FALSE)
    aov1 <- stats::anova(stats::lm(L$y ~ factor(L$gi)))
    expect_equal(A$f_p, aov1[["Pr(>F)"]][1L], tolerance = 1e-9, info = sprintf("layout %d", i))
    expect_true(A$f_p > 0 && A$f_p < 1, info = sprintf("layout %d", i))
  }
  expect_gt(.pool_fit(.pt_layout(9L)$y, .pt_layout(9L)$gi, 3L, .pt_layout(9L)$xy,
                      use_fnn = FALSE)$f_p, 0.05)
})


# ---- 8. A change of scale ----------------------------------------------------

test_that("8. y -> a + b y moves the values the same way and leaves the weights", {
  for (i in c(2L, 4L)) for (f in .pt_paths()) {
    L <- .pt_layout(i); tag <- sprintf("layout %d, %s", i, .pt_path_tag(f))
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    for (ab in list(c(37.5, -2.5), c(-4, 1000), c(0.3, 0.01))) {
      a <- ab[1L]; b <- ab[2L]
      B <- .pool_fit(a + b * L$y, L$gi, L$K, L$xy, use_fnn = f)
      rel <- function(got, want) max(abs(got - want)) / max(abs(want))
      expect_lt(rel(B$value, a + b * A$value), 1e-10)
      expect_lt(rel(B$m, a + b * A$m), 1e-10)
      expect_lt(rel(B$target, a + b * A$target), 1e-10)
      expect_lt(.pt_gap(B$w, A$w), 1e-10)
      expect_lt(rel(B$resid_pooled, b * A$resid_pooled), 1e-10)
      expect_lt(rel(B$resid_mean, b * A$resid_mean), 1e-10)
      expect_equal(B$sigma2, b^2 * A$sigma2, tolerance = 1e-10, info = tag)
      expect_equal(B$tau2, b^2 * A$tau2, tolerance = 1e-10, info = tag)
      expect_equal(B$loo_pooled, b^2 * A$loo_pooled, tolerance = 1e-10, info = tag)
      expect_equal(B$loo_mean, b^2 * A$loo_mean, tolerance = 1e-10, info = tag)
      expect_equal(B$loo_one, b^2 * A$loo_one, tolerance = 1e-10, info = tag)
      expect_equal(B$f_p, A$f_p, tolerance = 1e-8, info = tag)
      expect_identical(B$nb, A$nb, info = tag)
      expect_identical(B$cen, A$cen, info = tag)
      expect_identical(B$status, A$status, info = tag)
    }
  }
})


# ---- 9. Row order and numbering ----------------------------------------------

test_that("9. the fit does not depend on the order of the rows or the numbering of the regions", {
  cases <- list(`3,000 rows at 4.5 a region` = .pt_layout(4L),
                `the lattice in 667 regions` = .pt_lattice(667L),
                `the lattice in 1,500 regions` = .pt_lattice(1500L))
  # The lattice at 1,500 centres is there for the rows alone in their region:
  # at 667 it has none, and their neighbour is where a tie could go by row.
  expect_identical(sum(tabulate(cases[[2L]]$gi) == 1L), 0L)
  expect_gt(sum(tabulate(cases[[3L]]$gi) == 1L), 100L)
  for (nm in names(cases)) for (f in .pt_paths()) {
    L <- cases[[nm]]; tag <- sprintf("%s, %s", nm, .pt_path_tag(f))
    y <- L$y; gi <- L$gi; K <- L$K; xy <- L$xy; N <- length(y)
    A <- .pool_fit(y, gi, K, xy, use_fnn = f)
    # The centroids are the sums of the rows taken by first coordinate, then
    # by second: that expression and no other order that would do as well.
    o <- order(xy[, 1], xy[, 2])
    expect_identical(unname(A$cen), unname(rowsum(xy[o, , drop = FALSE], gi[o]) / A$n), info = tag)
    worst <- 0
    for (s in 1:20) {
      .pt_seed(100 + s)
      p   <- sample(N)                    # row i of the new order is row p[i]
      lab <- sample(K)                    # region j is now numbered lab[j]
      B <- .pool_fit(y[p], lab[gi[p]], K, xy[p, , drop = FALSE], use_fnn = f)
      # Every region has the same regions as neighbours, nearest first.
      expect_identical(B$nb[lab, , drop = FALSE], matrix(lab[A$nb], K), info = tag)
      # The centroids are the same numbers: they are summed in coordinate order.
      expect_identical(unname(B$cen[lab, , drop = FALSE]), unname(A$cen), info = tag)
      expect_identical(B$n[lab], A$n, info = tag)
      worst <- max(worst,
                   .pt_gap(B$S[lab], A$S), .pt_gap(B$m[lab], A$m),
                   .pt_gap(B$target[lab], A$target), .pt_gap(B$w[lab], A$w),
                   .pt_gap(B$value[lab], A$value),
                   .pt_gap(B$resid_pooled, A$resid_pooled[p]),
                   .pt_gap(B$resid_mean, A$resid_mean[p]))
      for (v in c("sigma2", "tau2", "Q", "loo_pooled", "loo_pooled_se", "loo_mean",
                  "loo_mean_se", "loo_one", "loo_one_se"))
        expect_equal(B[[v]], A[[v]], tolerance = 1e-12, info = paste(tag, v))
      expect_equal(B$f_p, A$f_p, tolerance = 1e-9, info = tag)
      expect_identical(B$status, A$status, info = tag)
    }
    expect_lt(worst, 1e-12)
  }
})

test_that("9. on the lattice, centroids summed in row order would move with the rows", {
  # What the test above guards against, shown on its own data: the centroid
  # of a region summed in the order of the rows is not the same number from
  # one row order to the next.
  L <- .pt_lattice(667L)
  n <- tabulate(L$gi, L$K)
  moved <- vapply(1:5, function(s) {
    .pt_seed(100 + s)
    p <- sample(length(L$y))
    !identical(unname(rowsum(L$xy[p, ], L$gi[p]) / n), unname(rowsum(L$xy, L$gi) / n))
  }, NA)
  expect_true(all(moved))
})


# ---- 10. Weights and values --------------------------------------------------

test_that("10. weights lie in [0, 1], rise with the count, and values lie between mean and target", {
  for (i in 1:7) for (f in .pt_paths()) {
    L <- .pt_layout(i); tag <- sprintf("layout %d, %s", i, .pt_path_tag(f))
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    expect_true(all(A$w >= 0 & A$w <= 1), info = tag)
    lo <- pmin(A$m, A$target); hi <- pmax(A$m, A$target)
    slack <- 1e-12 * max(abs(c(A$m, A$target)))
    expect_true(all(A$value >= lo - slack & A$value <= hi + slack), info = tag)
    # Within one fit a region with more rows has the larger weight, and
    # regions of one count have one weight.
    expect_gt(A$tau2, 0)
    by_n <- split(A$w, A$n)
    expect_true(all(vapply(by_n, function(v) length(unique(v)) == 1L, NA)), info = tag)
    one <- vapply(by_n, function(v) v[1L], numeric(1))
    expect_true(all(diff(one[order(as.integer(names(one)))]) > 0), info = tag)
    expect_equal(A$w, A$tau2 / (A$tau2 + A$sigma2 / A$n), tolerance = 1e-12, info = tag)
  }
})


# ---- 11. No spread -----------------------------------------------------------

test_that("11. a constant column gets outputs that are set, whatever its value", {
  L <- .pt_small()
  expect_identical(length(L$gi), 400L)
  expect_gte(L$K, 85L)
  for (f in .pt_paths()) for (const in c(0.1, 1234.56, 1 / 3, 52000.37, -7e-9, 0)) {
    tag <- sprintf("constant %s, %s", format(const, digits = 8), .pt_path_tag(f))
    A <- .pool_fit(rep(const, 400), L$gi, L$K, L$xy, use_fnn = f)
    expect_identical(A$sigma2, 0, info = tag)
    expect_identical(A$tau2, 0, info = tag)
    expect_identical(A$w, rep(1, L$K), info = tag)
    expect_identical(A$value, A$m, info = tag)
    expect_identical(A$resid_pooled, numeric(400), info = tag)
    expect_identical(A$resid_mean, numeric(400), info = tag)
    for (v in c("loo_pooled", "loo_pooled_se", "loo_mean", "loo_mean_se", "loo_one",
                "loo_one_se", "f_p"))
      expect_identical(A[[v]], NA_real_, info = paste(tag, v))
    expect_true(A$constant, info = tag)
    expect_true(A$flat, info = tag)
    expect_identical(A$status, "no spread inside any cell", info = tag)
    expect_lt(max(abs(A$value - const)), 1e-9 * max(abs(const), 1e-300))
  }
  # The arithmetic these outputs replace: for a constant that is not exact in
  # binary the sums of squares are rounding error and not zero, so the rule
  # has something to catch.  3.5 is exact and would not test it.
  y <- rep(0.1, 400)
  m <- as.numeric(rowsum(y, L$gi)) / tabulate(L$gi)
  expect_gt(sum((y - m[L$gi])^2), 0)
  expect_identical(sum((rep(3.5, 400) - 3.5)^2), 0)
})

test_that("11. a column constant inside every region and not overall is flat, not constant", {
  L <- .pt_small(); gi <- L$gi; K <- L$K; n <- tabulate(gi, K)
  alone <- which(n[gi] == 1L)
  shared <- n[gi] > 1L
  expect_gt(length(alone), 0L)
  for (f in .pt_paths()) {
    for (scaling in list(function(g) 0.1 * g, function(g) 1234.56 + g / 7,
                         function(g) as.numeric(g))) {
      y <- scaling(gi)
      tag <- sprintf("first value %s, %s", format(y[1L]), .pt_path_tag(f))
      A <- .pool_fit(y, gi, K, L$xy, use_fnn = f)
      expect_true(A$flat, info = tag)
      expect_false(A$constant, info = tag)
      expect_identical(A$sigma2, 0, info = tag)
      expect_identical(A$w, rep(1, K), info = tag)
      expect_identical(A$value, A$m, info = tag)
      expect_identical(A$status, "no spread inside any cell", info = tag)
      expect_identical(A$f_p, NA_real_, info = tag)
      # Residuals: exactly 0 for a row that shares its region; for a row
      # alone in its region, its distance from the target (pooled) and from
      # the mean of its nearest other row's region (plain).
      expect_true(all(A$resid_pooled[shared] == 0), info = tag)
      expect_true(all(A$resid_mean[shared] == 0), info = tag)
      expect_identical(A$resid_pooled[alone], y[alone] - A$target[gi[alone]], info = tag)
      near <- .pt_nearest_other(L$xy, gi, alone)
      expect_identical(A$resid_mean[alone], y[alone] - A$m[gi[near]], info = tag)
      expect_true(all(A$resid_mean[alone] != 0), info = tag)
      for (v in c("loo_pooled", "loo_pooled_se", "loo_mean", "loo_mean_se", "loo_one",
                  "loo_one_se", "tau2", "Q"))
        expect_true(is.finite(A[[v]]), info = paste(tag, v))
      expect_gt(A$tau2, 0)
    }
  }
})

test_that("11. the rule for no spread is relative to the size of the values", {
  L <- .pt_small()
  noise <- local({ .pt_seed(11); stats::rnorm(400) })
  A <- .pool_fit(1e6 + 1e-9 * noise, L$gi, L$K, L$xy, use_fnn = FALSE)
  expect_true(A$constant)
  expect_identical(A$status, "no spread inside any cell")
  expect_identical(A$w, rep(1, L$K))
  B <- .pool_fit(1e6 + 1e-3 * noise, L$gi, L$K, L$xy, use_fnn = FALSE)
  expect_false(B$constant)
  expect_false(B$flat)
  expect_gt(B$sigma2, 1e-7)
  expect_true(B$status %in% c("ok", "no signal between cells"))
  expect_true(is.finite(B$loo_pooled))
  # The rule itself.
  y <- c(1, 2, 3)
  lim <- 100 * .Machine$double.eps^2 * sum(y^2)
  expect_true(.no_spread(lim, y))
  expect_false(.no_spread(lim * 1.01, y))
  expect_true(.no_spread(0, c(0, 0)))
})

test_that("11. a constant column in regions of 500 rows gets the same outputs", {
  L <- .pt_large(4L); gi <- L$gi; K <- L$K; N <- length(gi)
  expect_identical(tabulate(gi), rep(500L, 4))
  for (f in .pt_paths()) for (const in c(0.1, 1234.56, 52000.37, 0.7)) {
    tag <- sprintf("constant %s, %s", format(const, digits = 8), .pt_path_tag(f))
    y <- rep(const, N)
    # Why these regions: the sum of 500 such values is off by a rounding
    # error, so the sum of squares inside the regions, as the fit computes
    # it, does not pass for no spread.  Only the test of the whole column, or
    # the one about corrected means, sees that this column has none.
    expect_false(.pt_written_rule(y, gi), info = tag)
    A <- .pool_fit(y, gi, K, L$xy, use_fnn = f)
    expect_true(A$constant, info = tag)
    expect_true(A$flat, info = tag)
    expect_identical(A$status, "no spread inside any cell", info = tag)
    expect_identical(A$sigma2, 0, info = tag)
    expect_identical(A$tau2, 0, info = tag)
    expect_identical(A$w, rep(1, K), info = tag)
    expect_identical(A$value, A$m, info = tag)
    expect_identical(A$resid_pooled, numeric(N), info = tag)
    expect_identical(A$resid_mean, numeric(N), info = tag)
    for (v in c("loo_pooled", "loo_pooled_se", "loo_mean", "loo_mean_se", "loo_one",
                "loo_one_se", "f_p"))
      expect_identical(A[[v]], NA_real_, info = paste(tag, v))
  }
})

test_that("11. a column constant inside regions of 500 rows is flat", {
  L <- .pt_large(12L); gi <- L$gi; K <- L$K; N <- length(gi); n <- tabulate(gi)
  expect_identical(n, rep(500L, 12))
  for (f in .pt_paths())
    for (scaling in list(function(g) 0.1 * g, function(g) 1234.56 + g / 7, function(g) g / 3)) {
      y <- scaling(gi)
      tag <- sprintf("first value %s, %s", format(y[1L], digits = 8), .pt_path_tag(f))
      # Not flat by the sum of squares the fit computes: the means carry the
      # rounding of their sums.
      expect_false(.pt_written_rule(y, gi), info = tag)
      A <- .pool_fit(y, gi, K, L$xy, use_fnn = f)
      expect_true(A$flat, info = tag)
      expect_false(A$constant, info = tag)
      expect_identical(A$status, "no spread inside any cell", info = tag)
      expect_identical(A$sigma2, 0, info = tag)
      expect_identical(A$w, rep(1, K), info = tag)
      expect_identical(A$value, A$m, info = tag)
      expect_identical(A$f_p, NA_real_, info = tag)
      # Every row shares its region, so every residual is exactly 0.
      expect_identical(A$resid_pooled, numeric(N), info = tag)
      expect_identical(A$resid_mean, numeric(N), info = tag)
      expect_identical(A$loo_pooled, 0, info = tag)
      expect_identical(A$loo_mean, 0, info = tag)
      expect_true(is.finite(A$loo_one) && A$loo_one > 0, info = tag)
      expect_gt(A$tau2, 0)
    }

  # In the column's record the status is the same, and `structured` is not
  # read for a column without spread.
  P <- .pool_column(0.1 * gi, gi, L$xy, obs_level = 0.9, want = "pooled", use_fnn = FALSE)
  expect_identical(P$record$status, "no spread inside any cell")
  expect_false(P$record$constant)
  expect_identical(P$record$sigma2, 0)
  expect_identical(P$record$structured, NA)
  expect_identical(P$record$offsets, c(0, 0))

  # Values that differ inside a region by a last bit or two have no spread
  # either: the rule is the bound of .no_spread(), not equality.
  y <- (0.1 * gi) * (1 + .Machine$double.eps * rep_len(c(-2, 1, 0, 2, -1), N))
  expect_gt(length(unique(y[gi == 1L])), 2L)
  expect_false(.pt_written_rule(y, gi))
  A <- .pool_fit(y, gi, K, L$xy, use_fnn = FALSE)
  expect_true(A$flat)
  expect_false(A$constant)
  expect_identical(A$status, "no spread inside any cell")
  expect_identical(A$w, rep(1, K))
  # A spread of one part in a million million inside the regions is spread.
  y <- (0.1 * gi) * (1 + 1e-12 * local({ .pt_seed(113); stats::rnorm(N) }))
  B <- .pool_fit(y, gi, K, L$xy, use_fnn = FALSE)
  expect_false(B$flat)
  expect_gt(B$sigma2, 0)
  expect_true(B$status %in% c("ok", "no signal between cells"))

  # The corrected sum of squares itself.  With no spread it is far below the
  # sum about the means as summed; with spread the two are one number to
  # rounding, so a column that has spread is never called flat by it.
  y <- 1234.56 + gi / 7
  res <- y - (as.numeric(rowsum(y, gi)) / n)[gi]
  expect_gt(sum(res^2), 100 * .Machine$double.eps^2 * sum(y^2))
  expect_lt(.ss_within(res, gi, n), 1e-6 * sum(res^2))
  y <- .pt_layout(4L)
  res <- y$y - (as.numeric(rowsum(y$y, y$gi)) / tabulate(y$gi))[y$gi]
  expect_equal(.ss_within(res, y$gi, tabulate(y$gi)), sum(res^2), tolerance = 1e-12)
  expect_identical(.ss_within(c(1, -1, 2, 4), c(1L, 1L, 2L, 2L), c(2L, 2L)), 4)
})

test_that("the own weight is tau2 over tau2 plus the mean's variance, and 1 where both are 0", {
  expect_identical(.pool_weight(1, 1), 0.5)
  expect_identical(.pool_weight(3, 1), 0.75)
  expect_identical(.pool_weight(0, 2), 0)
  expect_identical(.pool_weight(2, 0), 1)
  expect_identical(.pool_weight(0, 0), 1)
  expect_identical(.pool_weight(0.5, c(0.5, 1.5, 0)), c(0.5, 0.25, 1))
  expect_identical(.pool_weight(0, c(0, 4)), c(1, 0))
  # The reference's two forms of it.
  expect_identical(.pool_weight(0.37, 1.2 / 1:9), weight(0.37, 1.2, 1:9))
  expect_identical(.pool_weight(c(0, 0.37, 0), c(0, 0.4, 1)), wfun(c(0, 0.37, 0), c(0, 0.4, 1)))
})


# ---- 12. The order of the degenerate cases -----------------------------------

test_that("12. one region with several rows: the plain mean and the plain residuals", {
  y  <- c(2.5, 3.25, 1, 4.75, 3)
  xy <- cbind(c(0, 10, 20, 30, 45), c(5, 0, 12, 7, 3))
  for (f in .pt_paths()) {
    A <- .pool_fit(y, rep(1L, 5), 1L, xy, use_fnn = f)
    expect_identical(A$status, "one cell")
    expect_identical(A$m, 14.5 / 5)
    expect_identical(A$value, A$m)
    expect_identical(A$w, 1)
    expect_identical(A$target, NA_real_)
    expect_identical(A$tau2, NA_real_)
    expect_identical(A$Q, NA_real_)
    expect_identical(A$f_p, NA_real_)
    expect_identical(dim(A$nb), c(1L, 0L))
    expect_true(is.integer(A$nb))
    expect_equal(A$sigma2, sum((y - mean(y))^2) / 4, tolerance = 1e-14)
    expect_equal(A$sigma2, stats::var(y), tolerance = 1e-14)
    e <- vapply(1:5, function(i) y[i] - mean(y[-i]), numeric(1))
    expect_lt(.pt_gap(A$resid_pooled, e), 1e-14)
    expect_identical(A$resid_mean, A$resid_pooled)
    expect_identical(A$resid_pooled, (y - A$m) * 5 / 4)
    expect_equal(A$loo_pooled, mean(e^2), tolerance = 1e-14)
    expect_identical(A$loo_mean, A$loo_pooled)
    expect_identical(A$loo_one, A$loo_pooled)
    expect_equal(A$loo_pooled_se, stats::sd(e^2) / sqrt(5), tolerance = 1e-14)
    expect_identical(A$loo_one_se, A$loo_pooled_se)
    expect_equal(unname(A$cen), matrix(colMeans(xy), 1L), tolerance = 1e-14)
    expect_false(A$constant)
    expect_false(A$flat)
  }
})

test_that("12. one region and a constant column, and one region of one row", {
  xy <- cbind(c(0, 10, 20), c(5, 0, 12))
  y  <- rep(0.1, 3)
  # Three times 0.1 over 3 is not 0.1 in floating point, so the sum of
  # squares about the mean is rounding error and not zero.
  expect_gt(sum((y - sum(y) / 3)^2), 0)
  A <- .pool_fit(y, rep(1L, 3), 1L, xy, use_fnn = FALSE)
  expect_identical(A$status, "one cell")           # not "no spread inside any cell"
  expect_identical(A$sigma2, 0)
  expect_identical(A$resid_pooled, numeric(3))
  expect_identical(A$resid_mean, numeric(3))
  for (v in c("loo_pooled", "loo_pooled_se", "loo_mean", "loo_mean_se", "loo_one",
              "loo_one_se", "tau2", "target", "f_p"))
    expect_identical(A[[v]], NA_real_, info = v)
  expect_identical(A$w, 1)
  expect_identical(A$value, A$m)
  expect_true(A$constant)
  expect_true(A$flat)

  B <- .pool_fit(7.25, 1L, 1L, cbind(3, 4), use_fnn = FALSE)
  expect_identical(B$status, "one cell")
  expect_identical(B$value, 7.25)
  expect_identical(B$w, 1)
  expect_identical(B$N, 1L)
  expect_identical(B$K, 1L)
  expect_identical(B$sigma2, NA_real_)
  expect_identical(B$resid_pooled, NA_real_)
  expect_identical(B$resid_mean, NA_real_)
  for (v in c("loo_pooled", "loo_pooled_se", "loo_mean", "loo_mean_se", "loo_one",
              "loo_one_se", "tau2", "target", "Q", "f_p"))
    expect_identical(B[[v]], NA_real_, info = v)
  expect_identical(dim(B$nb), c(1L, 0L))
  # The interval read from these residuals has no ends.
  expect_identical(.obs_offsets(B$resid_pooled, 0.9), c(NA_real_, NA_real_))
})

test_that("12. no region holds two rows: a status and no fit", {
  y  <- c(4, 1, 3, 9, 2, 6, 5, 8, 7, 10) / 7
  xy <- cbind(seq(0, 90, by = 10), c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3))
  for (case in list(y, rep(0.1, 10))) {       # with spread, and constant
    A <- .pool_fit(case, 1:10, 10L, xy, use_fnn = FALSE)
    expect_identical(A$status, "no cell holds two values")
    expect_identical(A$n, rep(1L, 10))
    expect_identical(A$S, case)
    expect_identical(A$m, case)
    expect_identical(A$N, 10L)
    expect_identical(A$K, 10L)
    expect_null(A$resid_pooled)
    expect_null(A$resid_mean)
    for (v in c("sigma2", "cen", "nb", "target", "Q", "tau2", "w", "value", "loo_pooled",
                "loo_pooled_se", "loo_mean", "loo_mean_se", "loo_one", "loo_one_se", "f_p"))
      expect_true(all(is.na(A[[v]])), info = v)
    expect_true(A$flat)
  }
  expect_false(.pool_fit(y, 1:10, 10L, xy, use_fnn = FALSE)$constant)
  expect_true(.pool_fit(rep(0.1, 10), 1:10, 10L, xy, use_fnn = FALSE)$constant)
  # No search is made for a fit that is not made.
  local_mocked_bindings(.region_neighbours = function(...) stop("searched"),
                        .nearest_other = function(...) stop("searched"))
  expect_identical(.pool_fit(y, 1:10, 10L, xy, use_fnn = FALSE)$status,
                   "no cell holds two values")
})

test_that("12. region means that are all equal: no signal between cells", {
  # Thirty regions of four rows whose values sum to exactly zero in each, so
  # every plain mean and every target is 0 by construction.
  cen <- local({ .pt_seed(12); cbind(stats::runif(30, 0, 1e5), stats::runif(30, 0, 1e5)) })
  gi  <- rep(1:30, each = 4L)
  xy  <- cen[gi, ] + cbind(rep(c(-40, 40, -40, 40), 30), rep(c(-40, -40, 40, 40), 30))
  # Each value is followed by its negative, so the running sum of a region
  # returns to exactly 0 whatever the scale.
  y   <- rep(c(-1.5, 1.5, -0.5, 0.5), 30) * rep(seq(1, 2, length.out = 30), each = 4L)
  for (f in .pt_paths()) {
    A <- .pool_fit(y, gi, 30L, xy, use_fnn = f)
    expect_identical(A$m, numeric(30))
    expect_identical(A$status, "no signal between cells")
    expect_identical(A$tau2, 0)
    expect_gt(A$sigma2, 0)
    expect_identical(A$w, numeric(30))
    expect_identical(A$value, A$target)
    expect_identical(A$Q, 0)
    expect_false(A$flat)
    expect_false(A$constant)
    expect_true(is.finite(A$loo_pooled) && A$loo_pooled > 0)
    expect_identical(A$f_p, 1)
    # With a weight of 0 a left-out row is predicted by its region's target.
    expect_identical(A$resid_pooled, y - A$target[gi])
  }
  # A column with no spread at all has no signal either; the earlier case names it.
  expect_identical(.pool_fit(rep(0.1, 120), gi, 30L, xy, use_fnn = FALSE)$status,
                   "no spread inside any cell")
})

test_that("12. one region: `want` decides which residuals are returned", {
  y  <- c(2.5, 3.25, 1, 4.75, 3)
  xy <- cbind(c(0, 10, 20, 30, 45), c(5, 0, 12, 7, 3))
  both   <- .pool_fit(y, rep(1L, 5), 1L, xy, use_fnn = FALSE)
  pooled <- .pool_fit(y, rep(1L, 5), 1L, xy, want = "pooled", use_fnn = FALSE)
  plain  <- .pool_fit(y, rep(1L, 5), 1L, xy, want = "mean", use_fnn = FALSE)
  shared <- c("n", "S", "m", "N", "K", "sigma2", "cen", "nb", "target", "Q", "tau2", "w",
              "value", "loo_one", "loo_one_se", "f_p", "constant", "flat", "status")
  expect_identical(pooled[shared], both[shared])
  expect_identical(plain[shared], both[shared])
  expect_identical(names(pooled), names(both))
  expect_identical(names(plain), names(both))
  expect_true(is.finite(both$loo_pooled) && is.finite(both$loo_mean))
  # The one asked for is the full fit's; the other has no vector and NA for
  # its two numbers.
  expect_identical(pooled[c("resid_pooled", "loo_pooled", "loo_pooled_se")],
                   both[c("resid_pooled", "loo_pooled", "loo_pooled_se")])
  expect_null(pooled$resid_mean)
  expect_identical(pooled$loo_mean, NA_real_)
  expect_identical(pooled$loo_mean_se, NA_real_)
  expect_identical(plain[c("resid_mean", "loo_mean", "loo_mean_se")],
                   both[c("resid_mean", "loo_mean", "loo_mean_se")])
  expect_null(plain$resid_pooled)
  expect_identical(plain$loo_pooled, NA_real_)
  expect_identical(plain$loo_pooled_se, NA_real_)
})

test_that("12. the fit refuses regions that are not numbered 1 to K with a row in each", {
  xy <- cbind(1:6, c(2, 4, 1, 3, 6, 5))
  y  <- c(1, 2, 4, 3, 6, 5)
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 4L, 4L), 4L, xy, use_fnn = FALSE), "numbered 1 to K")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L, 3L), 2L, xy, use_fnn = FALSE), "numbered 1 to K")
  # Three numbers in use and K = 3, but one of them is not among 1 to K.
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 4L, 4L), 3L, xy, use_fnn = FALSE), "numbered 1 to K")
  # Every region of 1 to K holds a row, and one row is in none of them.
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L, 4L), 3L, xy, use_fnn = FALSE), "numbered 1 to K")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L, 0L), 3L, xy, use_fnn = FALSE), "numbered 1 to K")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, NA, 3L), 3L, xy, use_fnn = FALSE), "numbered 1 to K")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L, 3L), 3L, xy[1:5, ], use_fnn = FALSE),
               "one row per value")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L, 3L), 3L, xy, want = character(0)), "`want`")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L, 3L), 3L, xy, want = "median"), "`want`")
  # No value at all, and a region for some of the values only.
  expect_error(.pool_fit(numeric(0), integer(0), 0L, matrix(0, 0L, 2L), use_fnn = FALSE),
               "at least one value")
  expect_error(.pool_fit(y, c(1L, 1L, 2L, 2L, 3L), 3L, xy, use_fnn = FALSE), "the region of each")
})

test_that("12. a coordinate matrix that is not two columns of numbers is refused", {
  # With a third column the kd-tree would search in three dimensions and the
  # exact search in two, and the two paths would give different neighbours.
  L <- .pt_layout(2L); n <- tabulate(L$gi)
  xyz <- cbind(L$xy, 100 * seq_along(L$y))
  for (f in .pt_paths()) {
    tag <- .pt_path_tag(f)
    expect_error(.pool_fit(L$y, L$gi, L$K, xyz, use_fnn = f), "two columns", info = tag)
    expect_error(.pool_fit(L$y, L$gi, L$K, xyz, want = "pooled", use_fnn = f,
                           geo = .pool_geometry(L$gi, L$K, L$xy, n, f)),
                 "two columns", info = tag)
    expect_error(.pool_geometry(L$gi, L$K, xyz, n, f), "two columns", info = tag)
    expect_error(.nearest_other(xyz, L$gi, 1:5, f), "two columns", info = tag)
    expect_error(.region_neighbours(xyz[1:30, ], 8L, f), "two columns", info = tag)
    # A column's rows are handed over as their first two coordinates, so a
    # layer with heights gets the answer of the layer without, on both paths.
    expect_identical(.pool_column(L$y, L$gi, xyz, 0.9, c("pooled", "mean"), f),
                     .pool_column(L$y, L$gi, L$xy, 0.9, c("pooled", "mean"), f), info = tag)
  }
  expect_error(.pool_fit(L$y, L$gi, L$K, as.data.frame(L$xy), use_fnn = FALSE), "numeric matrix")
  expect_error(.pool_fit(L$y, L$gi, L$K, L$xy[, 1], use_fnn = FALSE), "numeric matrix")
  expect_error(.pool_fit(L$y, L$gi, L$K, matrix("a", nrow(L$xy), 2L), use_fnn = FALSE),
               "numeric matrix")
  # Two columns of numbers that are no matrix: an array of three dimensions.
  expect_error(.pool_fit(L$y, L$gi, L$K, array(L$xy, c(nrow(L$xy), 2L, 1L)), use_fnn = FALSE),
               "numeric matrix")
  expect_error(.pool_geometry(L$gi, L$K, L$xy[-1, ], n, FALSE), "one row per value")
  expect_error(.nearest_other(L$xy[-1, ], L$gi, 1:5, FALSE), "one row per value")
})

test_that("12. whole numbers are fitted as doubles, values and coordinates alike", {
  # rowsum() of integers returns NA, without a warning, once a sum passes
  # .Machine$integer.max.  Here every region's sum of the values does, and
  # every region's sum of the second coordinate.
  L <- .pt_large(4L); gi <- L$gi; K <- L$K; N <- length(gi)
  y  <- as.integer(round(5e6 + 1e5 * local({ .pt_seed(121); stats::rnorm(N) })))
  xi <- round(cbind(4e5, 5.1e6)[rep(1L, N), ] + L$xy)
  xd <- xi
  storage.mode(xi) <- "integer"
  expect_true(is.integer(y) && is.integer(xi) && is.double(xd))
  expect_true(all(is.na(suppressWarnings(rowsum(y, gi)))))
  expect_true(all(is.na(suppressWarnings(rowsum(xi[, 2], gi)))))
  for (f in .pt_paths()) {
    tag <- .pt_path_tag(f)
    A <- .pool_fit(as.numeric(y), gi, K, xd, use_fnn = f)
    expect_identical(A$status, "ok", info = tag)
    expect_false(anyNA(A$value) || anyNA(A$cen), info = tag)
    expect_identical(.pool_fit(y, gi, K, xd, use_fnn = f), A, info = tag)
    expect_identical(.pool_fit(as.numeric(y), gi, K, xi, use_fnn = f), A, info = tag)
    expect_identical(.pool_fit(y, gi, K, xi, use_fnn = f), A, info = tag)
    expect_identical(.pool_geometry(gi, K, xi, A$n, f), A[c("cen", "nb")], info = tag)
    expect_identical(.nearest_other(xi, gi, seq_len(N), f), .nearest_other(xd, gi, seq_len(N), f),
                     info = tag)
    # The same through a column of an assigned layer.
    P <- .pool_column(as.numeric(y), gi, xd, 0.9, c("pooled", "mean"), f)
    expect_identical(P$fit, A, info = tag)
    expect_identical(.pool_column(y, gi, xi, 0.9, c("pooled", "mean"), f), P, info = tag)
  }
  # One region, the case in which an integer sum past the limit gave NA and
  # the status "one cell".
  one <- .pool_fit(y, rep(1L, N), 1L, xd, use_fnn = FALSE)
  expect_identical(one$value, sum(as.numeric(y)) / N)
  expect_identical(one, .pool_fit(as.numeric(y), rep(1L, N), 1L, xd, use_fnn = FALSE))
})


# ---- 13. Neighbours of the regions -------------------------------------------

# Centroids on which a search by distance alone has no single answer.
.pt_centroid_sets <- function() {
  .pt_seed(13)
  base <- cbind(stats::runif(40, 0, 1e5), stats::runif(40, 0, 1e5))
  coincident <- base
  coincident[12, ] <- coincident[5, ]                        # a pair
  coincident[33, ] <- coincident[20, ]                       # another
  coincident[8, ] <- coincident[30, ] <- coincident[7, ]     # a triple
  list(`two pairs and a triple of coincident centroids` = coincident,
       `a lattice of centroids` =
         unname(as.matrix(expand.grid(seq(0.1, by = 0.7, length.out = 30),
                                      seq(0.3, by = 0.7, length.out = 30)))),
       `centroids rounded to the kilometre` =
         round(cbind(stats::runif(600, 0, 3e4), stats::runif(600, 0, 3e4)), -3),
       `two regions` = base[1:2, ],
       `two regions on one centroid` = base[c(3, 3), ],
       `three regions` = base[1:3, ],
       `three regions on a line` = cbind(c(0, 1, 2), 0),
       `nine regions` = base[1:9, ],
       `ten regions` = base[1:10, ],
       `no ties` = base)
}

test_that("13. no region is its own neighbour, and ties go by coordinate, then by number", {
  sets <- .pt_centroid_sets()
  expect_gt(anyDuplicated(sets[[1L]]), 0L)
  expect_gt(anyDuplicated(sets[[3L]]), 0L)
  for (nm in names(sets)) for (f in .pt_paths()) {
    cen <- sets[[nm]]; K <- nrow(cen); tag <- sprintf("%s, %s", nm, .pt_path_tag(f))
    nb <- .region_neighbours(cen, 8L, f)
    expect_true(is.integer(nb), info = tag)
    expect_identical(dim(nb), c(K, min(8L, K - 1L)), info = tag)
    expect_false(any(nb == seq_len(K)), info = tag)
    expect_true(all(nb >= 1L & nb <= K), info = tag)
    expect_true(all(apply(nb, 1L, anyDuplicated) == 0L), info = tag)
    expect_identical(nb, .pt_neighbours(cen, 8L), info = tag)
    # Nearest first.
    d <- matrix(sqrt((cen[nb, 1] - cen[, 1])^2 + (cen[nb, 2] - cen[, 2])^2), K)
    if (ncol(d) > 1L) expect_true(all(d[, -1L] >= d[, -ncol(d)]), info = tag)
    # Another number of neighbours.
    expect_identical(.region_neighbours(cen, 3L, f), .pt_neighbours(cen, 3L), info = tag)
    # Renumbered, the regions keep their neighbours, but for regions on one
    # centroid, which only their number tells apart.
    if (!anyDuplicated(cen)) {
      lab <- local({ .pt_seed(131); sample(K) })
      cen2 <- cen; cen2[lab, ] <- cen
      expect_identical(.region_neighbours(cen2, 8L, f)[lab, , drop = FALSE],
                       matrix(lab[nb], K), info = tag)
    }
  }
  # One region has no neighbour, and none is searched for.
  expect_identical(.region_neighbours(cbind(1, 2), 8L, FALSE), matrix(0L, 1L, 0L))
})

test_that("13. FNN's path and the exact path return the same neighbours", {
  skip_if_not_installed("FNN")
  sets <- .pt_centroid_sets()
  for (nm in names(sets)) {
    cen <- sets[[nm]]
    for (k in c(8L, 3L, 1L))
      expect_identical(.region_neighbours(cen, k, TRUE), .region_neighbours(cen, k, FALSE),
                       info = sprintf("%s, %d neighbours", nm, k))
  }
  # On coincident centroids the kd-tree alone returns a region as its own
  # neighbour, which is why such rows are searched again.
  co <- sets[[1L]]
  raw <- FNN::get.knn(co, k = 8L)$nn.index
  expect_true(any(raw == seq_len(nrow(co))))
  # Where nothing ties, the neighbours are the kd-tree's own.
  expect_identical(.region_neighbours(sets[["no ties"]], 8L, TRUE),
                   FNN::get.knn(sets[["no ties"]], k = 8L)$nn.index)
  # A whole fit on coincident centroids: the two paths agree, and no region's
  # own sum is in its target.
  L <- .pt_layout(2L)
  xy <- L$xy
  twins <- which(L$gi %in% c(4L, 9L))
  # Regions 4 and 9 on one point, in whole metres so that both centroids
  # are that point exactly.
  xy[twins, ] <- round(xy[rep(twins[1L], length(twins)), ])
  A <- .pool_fit(L$y, L$gi, L$K, xy, use_fnn = TRUE)
  B <- .pool_fit(L$y, L$gi, L$K, xy, use_fnn = FALSE)
  expect_identical(unname(A$cen[4, ]), unname(A$cen[9, ]))
  expect_identical(A, B)
  expect_false(any(A$nb == seq_len(L$K)))
})

# Stand-ins for FNN's two searches that return what a kd-tree with other
# rounding than R's could return.  Such a tree ranks candidates at nearly one
# distance in another order than the distances R computes, so it can hand
# back the further members of such a group and leave the nearest out.  These
# do that and nothing else: where the candidates around the end of the list
# are within one part in 1e9 of one squared distance, as the package computes
# it, the nearest of them is dropped and the next one taken in.  `count$rows`
# counts the rows it happened to.
.pt_by_distance <- function(nn, d2) {
  o  <- t(apply(d2, 1L, order))
  at <- cbind(as.vector(row(o)), as.vector(o))
  list(nn = matrix(nn[at], nrow(nn)), d2 = matrix(d2[at], nrow(nn)))
}
.pt_push_knn <- function(real, count) {
  force(real); force(count)             # FNN's own function, read before it is replaced
  function(data, k = 10, ...) {
    n <- nrow(data)
    if (k < 2L || k + 1L > n - 1L) return(real(data, k = k, ...))
    nn <- real(data, k = k + 1L, ...)$nn.index
    d2 <- matrix((data[nn, 1] - data[, 1])^2 + (data[nn, 2] - data[, 2])^2, n)
    s  <- .pt_by_distance(nn, d2)
    # Places k - 1, k and k + 1 at nearly one distance: k - 1 is left out.
    group <- s$d2[, k + 1L] <= s$d2[, k - 1L] * (1 + 1e-9)
    out <- s$nn[, seq_len(k), drop = FALSE]
    out[group, ] <- s$nn[group, c(seq_len(k - 2L), k, k + 1L), drop = FALSE]
    count$rows <- count$rows + sum(group)
    list(nn.index = out)
  }
}
.pt_push_knnx <- function(real, count) {
  force(real); force(count)
  function(data, query, k = 10, ...) {
    if (k != 3L || nrow(data) < 5L) return(real(data, query, k = k, ...))
    nn <- real(data, query, k = 4L, ...)$nn.index
    d2 <- matrix((data[nn, 1] - query[, 1])^2 + (data[nn, 2] - query[, 2])^2, nrow(query))
    s  <- .pt_by_distance(nn, d2)
    # The location asked about stands first, at distance 0.  Its three
    # nearest others at nearly one distance: the nearest of them is left out.
    group <- s$d2[, 4L] <= s$d2[, 2L] * (1 + 1e-9)
    out <- s$nn[, 1:3, drop = FALSE]
    out[group, ] <- s$nn[group, c(1L, 3L, 4L), drop = FALSE]
    count$rows <- count$rows + sum(group)
    list(nn.index = out)
  }
}

test_that("13. neighbours a kd-tree returns out of order within rounding are searched again", {
  skip_if_not_installed("FNN")
  # On this platform FNN's distances are R's to the last bit, so the lattices
  # above give one answer whether or not near ties are searched again.  The
  # margin of one part in 1e9 is for a kd-tree built with other arithmetic;
  # the stand-in is one.
  square <- .pt_lattice(667L)$xy
  a <- expand.grid(i = 0:29, j = 0:24)
  hexagonal <- unname(cbind(a$i + 0.5 * (a$j %% 2), a$j * sqrt(3) / 2)) * 100 + 0.1
  for (case in list(list("square lattice, 2 neighbours", square, 2L),
                    list("square lattice, 5 neighbours", square, 5L),
                    list("hexagonal lattice, 8 neighbours", hexagonal, 8L))) {
    u <- case[[2L]]; k <- case[[3L]]
    want  <- .region_neighbours(u, k, FALSE)
    count <- new.env(); count$rows <- 0L
    got <- local({
      local_mocked_bindings(get.knn = .pt_push_knn(FNN::get.knn, count), .package = "FNN")
      .region_neighbours(u, k, TRUE)
    })
    expect_gt(count$rows, nrow(u) %/% 2L)        # most rows of a lattice
    expect_identical(got, want, info = case[[1L]])
  }
})

# 1,316 centroids of which 811 share their point with eight others or more:
# 601 on one point, 201 on another and 9 on a third.  8 on a fourth and a
# pair on a fifth are one short of that and are searched for like any other.
.pt_point_key <- function(xy) sprintf("%.17g %.17g", xy[, 1], xy[, 2])
.pt_crowded <- function() .pt_keep("crowded", {
  .pt_seed(132)
  base <- cbind(stats::runif(500, 0, 1e5), stats::runif(500, 0, 1e5))
  cen  <- rbind(base, base[rep(1L, 600), ], base[rep(2L, 200), ], base[rep(3L, 8), ],
                base[rep(4L, 7), ], base[rep(5L, 1), ])
  cen  <- cen[sample(nrow(cen)), ]
  key  <- .pt_point_key(cen)
  at   <- match(key, key)
  list(cen = cen, key = key, crowd = tabulate(at, nrow(cen))[at] > 8L)
})

test_that("13. centroids that share a point with eight others or more are not searched for", {
  S <- .pt_crowded(); cen <- S$cen; K <- nrow(cen); crowd <- S$crowd
  expect_identical(K, 1316L)
  expect_identical(sum(crowd), 601L + 201L + 9L)
  rule <- .pt_neighbours(cen, 8L)
  expect_false(any(rule == seq_len(K)))

  # The exact path searches for the others only, and only they count as work.
  searched <- NULL
  exact <- .knn_exact_rows
  local({
    local_mocked_bindings(.knn_exact_rows = function(x, y, rows, k) {
      searched <<- c(searched, rows)
      exact(x, y, rows, k)
    })
    expect_identical(.region_neighbours(cen, 8L, FALSE), rule)
  })
  expect_length(searched, sum(!crowd))
  pairs <- as.numeric(sum(!crowd)) * K
  withr::with_options(list(spatialkit.exact_pairs_max = pairs - 1), {
    e <- tryCatch(.region_neighbours(cen, 8L, FALSE), error = function(e) e)
    expect_s3_class(e, "spatialkit_work_limit")
    expect_identical(e$work, pairs)
  })
  withr::with_options(list(spatialkit.exact_pairs_max = pairs),
                      expect_identical(.region_neighbours(cen, 8L, FALSE), rule))
  # Fewer neighbours asked for: more points are crowded, and the rule holds.
  expect_identical(.region_neighbours(cen, 1L, FALSE), .pt_neighbours(cen, 1L))
  expect_identical(.region_neighbours(cen, 7L, FALSE), .pt_neighbours(cen, 7L))

  # Every centroid on one point: its neighbours are the first eight others
  # by number, and nothing is searched, whatever the limit.
  one <- cen[rep(which(crowd)[1L], 60L), ]
  withr::with_options(list(spatialkit.exact_pairs_max = 0), {
    nb <- .region_neighbours(one, 8L, FALSE)
    expect_identical(nb, .pt_neighbours(one, 8L))
    expect_identical(nb[1L, ], 2:9)
    expect_identical(nb[5L, ], c(1:4, 6:9))
    expect_identical(nb[9L, ], 1:8)
    expect_identical(nb[60L, ], 1:8)
  })
  # A whole fit of 40 regions of two rows, every row on one point in whole
  # metres: no search, and no region's own sum in its target.
  gi <- rep(1:40, each = 2L)
  xy <- matrix(c(431250, 5123456), 80L, 2L, byrow = TRUE)
  y  <- local({ .pt_seed(133); stats::rnorm(80) })
  A <- withr::with_options(list(spatialkit.exact_pairs_max = 0),
                           .pool_fit(y, gi, 40L, xy, use_fnn = FALSE))
  expect_identical(A$nb, .pt_neighbours(unname(A$cen), 8L))
  expect_identical(A$nb[1L, ], 2:9)
  expect_equal(A$target[1L], mean(y[gi %in% 2:9]), tolerance = 1e-14)

  skip_if_not_installed("FNN")
  expect_identical(.pool_fit(y, gi, 40L, xy, use_fnn = TRUE), A)
  # With FNN: the same neighbours; no centroid of a crowded point is handed
  # to the search that settles ties, which would list every other centroid
  # on its point for each of them; and what the kd-tree is asked for stays
  # in proportion to the centroids, not to the square of a crowd.
  tied  <- NULL
  asked <- 0
  settle <- .knn_tied_rows
  knnx   <- FNN::get.knnx
  got <- local({
    local_mocked_bindings(.knn_tied_rows = function(u, rows, k) {
      tied <<- rbind(tied, u[rows, , drop = FALSE])
      settle(u, rows, k)
    })
    local_mocked_bindings(get.knnx = function(data, query, k = 10, ...) {
      asked <<- asked + as.numeric(nrow(query)) * k
      knnx(data, query, k = k, ...)
    }, .package = "FNN")
    .region_neighbours(cen, 8L, TRUE)
  })
  expect_identical(got, rule)
  expect_gt(nrow(tied), 0L)
  expect_false(any(.pt_point_key(tied) %in% S$key[crowd]))
  expect_true(any(.pt_point_key(tied) %in% S$key[!crowd]))
  # One run: 27,362 entries, 21 for each centroid.  Searching for the crowded
  # ones as for the others asked for 1,518,789.
  expect_lt(asked, 60 * K)
  for (k in c(1L, 7L))
    expect_identical(.region_neighbours(cen, k, TRUE), .region_neighbours(cen, k, FALSE),
                     info = sprintf("%d neighbours", k))
  # Every centroid on one point: the kd-tree is not asked at all.
  local({
    local_mocked_bindings(get.knn = function(...) stop("the kd-tree was asked"),
                          get.knnx = function(...) stop("the kd-tree was asked"),
                          .package = "FNN")
    expect_identical(.region_neighbours(one, 8L, TRUE), .pt_neighbours(one, 8L))
  })
})


# ---- 14. The nearest other row -----------------------------------------------

.pt_row_sets <- function() {
  .pt_seed(5)
  xy  <- cbind(stats::runif(3000, 0, 1e5), stats::runif(3000, 0, 1e5))
  lat <- unname(as.matrix(expand.grid(1:40, 1:40))) * 0.7
  rep_xy <- xy[c(1:2000, sample(2000, 1000, TRUE)), ]
  rep_gi <- sample(700L, 3000, TRUE)
  list(random = list(xy = xy, gi = .pt_regions(xy, 667L, 2L)),
       lattice = list(xy = lat, gi = .pt_regions(lat, 500L, 4L)),
       `repeated locations, regions at random` = list(xy = rep_xy, gi = rep_gi),
       `twenty rows on two locations` = list(xy = cbind(rep(c(0, 1), 10), 0), gi = rep(1:5, 4)),
       `three rows on a line` = list(xy = cbind(c(0, 1, 2), 0), gi = 1:3))
}
.pt_rows <- function() .pt_keep("rows", {
  sets <- .pt_row_sets()
  lapply(sets, function(s) c(s, list(rule = .pt_nearest_other(s$xy, s$gi))))
})

test_that("14. the nearest other row follows the rule, whatever the order of the rows", {
  sets <- .pt_rows()
  expect_identical(anyDuplicated(sets$random$xy), 0L)
  expect_gt(sum(duplicated(sets[[3L]]$xy)), 500L)
  for (nm in names(sets)) for (f in .pt_paths()) {
    xy <- sets[[nm]]$xy; gi <- sets[[nm]]$gi; n <- nrow(xy); rows <- seq_len(n)
    tag <- sprintf("%s, %s", nm, .pt_path_tag(f))
    a <- .nearest_other(xy, gi, rows, f)
    expect_true(is.integer(a), info = tag)
    expect_length(a, n)
    expect_false(any(a == rows), info = tag)
    expect_identical(a, sets[[nm]]$rule, info = tag)
    # A subset of the rows, in any order, gets the rows' own answers.
    pick <- c(n, 1L, rows[rows %% 7L == 3L])
    expect_identical(.nearest_other(xy, gi, pick, f), a[pick], info = tag)
    expect_identical(.nearest_other(xy, gi, integer(0), f), integer(0), info = tag)
    # Under another order of the rows the neighbour is at the same location
    # and in the same region.
    for (s in 1:10) {
      .pt_seed(500 + s)
      p  <- sample(n)
      ap <- .nearest_other(xy[p, , drop = FALSE], gi[p], rows, f)
      expect_identical(gi[p][ap], gi[a[p]], info = tag)
      expect_identical(xy[p, , drop = FALSE][ap, , drop = FALSE], xy[a[p], , drop = FALSE],
                       info = tag)
    }
  }
  # A lattice in whole numbers is where the rule shows: every interior
  # location has four others at exactly one distance, and the one with the
  # lower first coordinate is taken, in any order of the rows.
  grid <- cbind(as.numeric(rep(1:7, 6)), as.numeric(rep(1:6, each = 7)))
  p <- local({ .pt_seed(14); sample(nrow(grid)) })
  for (f in .pt_paths()) {
    a <- .nearest_other(grid[p, ], rep(1:6, 7), seq_len(nrow(grid)), f)
    from <- grid[p, ]; to <- grid[p, ][a, ]
    inner <- from[, 1] > 1
    expect_identical(to[inner, 1], from[inner, 1] - 1)
    expect_identical(to[inner, 2], from[inner, 2])
    # On the left edge there is no such location; the nearest is above or below.
    expect_identical(to[!inner, 1], from[!inner, 1])
    expect_identical(abs(to[!inner, 2] - from[!inner, 2]), rep(1, sum(!inner)))
  }
  # One row has no other.
  expect_identical(.nearest_other(cbind(1, 2), 1L, 1L, FALSE), NA_integer_)
})

test_that("14. FNN's path and the exact path return the same rows, and the kd-tree's where nothing ties", {
  skip_if_not_installed("FNN")
  sets <- .pt_rows()
  for (nm in names(sets)) {
    xy <- sets[[nm]]$xy; gi <- sets[[nm]]$gi; rows <- seq_len(nrow(xy))
    expect_identical(.nearest_other(xy, gi, rows, TRUE), .nearest_other(xy, gi, rows, FALSE),
                     info = nm)
  }
  r <- sets$random
  expect_identical(.nearest_other(r$xy, r$gi, seq_len(3000), TRUE),
                   FNN::get.knn(r$xy, k = 1L)$nn.index[, 1])
  expect_identical(.nearest_other(r$xy, r$gi, seq_len(3000), FALSE),
                   FNN::get.knn(r$xy, k = 1L)$nn.index[, 1])
})

test_that("14. a nearest location a kd-tree returns out of order within rounding is searched again", {
  skip_if_not_installed("FNN")
  # The stand-in of test 13 for the other search: where a location's three
  # nearest others are at nearly one distance it returns the second and the
  # third and leaves the nearest out.
  L <- .pt_lattice(1500L); rows <- seq_along(L$y)
  want  <- .nearest_other(L$xy, L$gi, rows, FALSE)
  count <- new.env(); count$rows <- 0L
  got <- local({
    local_mocked_bindings(get.knnx = .pt_push_knnx(FNN::get.knnx, count), .package = "FNN")
    .nearest_other(L$xy, L$gi, rows, TRUE)
  })
  expect_gt(count$rows, length(rows) %/% 2L)
  expect_identical(got, want)
})

test_that("14. with a cache no location is searched twice", {
  sets <- .pt_rows()
  search <- .nearest_location_search
  for (nm in c("random", "lattice", "repeated locations, regions at random"))
    for (f in .pt_paths()) {
      xy <- sets[[nm]]$xy; gi <- sets[[nm]]$gi; n <- nrow(xy); rows <- seq_len(n)
      tag <- sprintf("%s, %s", nm, .pt_path_tag(f))
      plain <- .nearest_other(xy, gi, rows, f)
      asked <- list()
      local_mocked_bindings(.nearest_location_search = function(u, q, ...) {
        asked[[length(asked) + 1L]] <<- q
        search(u, q, ...)
      })
      cache <- new.env()
      half <- rows[rows %% 2L == 0L]
      expect_identical(.nearest_other(xy, gi, half, f, cache), plain[half], info = tag)
      expect_length(asked, 1L)
      # The rest of the rows: only locations not yet found are searched.
      expect_identical(.nearest_other(xy, gi, rows, f, cache), plain, info = tag)
      expect_length(asked, 2L)
      expect_length(intersect(asked[[1L]], asked[[2L]]), 0L)
      # Again, and with the rows in other regions: nothing is searched.
      expect_identical(.nearest_other(xy, gi, rows, f, cache), plain, info = tag)
      gi2 <- rev(gi)
      expect_identical(.nearest_other(xy, gi2, rows, f, cache),
                       .nearest_other(xy, gi2, rows, f), info = tag)
      expect_length(asked, 3L)      # the call without a cache, and no other
      # A cache filled on other locations starts again and is not read.  The
      # two coordinates change places, so the locations stand in another order
      # and what the cache holds would be wrong for them.
      other <- cbind(xy[, 2], xy[, 1])
      expect_identical(.nearest_other(other, gi, rows, f, cache),
                       .nearest_other(other, gi, rows, f), info = tag)
      expect_identical(cache$u, .partition_locations(other)$u, info = tag)
    }
})

test_that("14. a fit with a cache is the fit without, and searches once over several fits", {
  L <- .pt_layout(7L)
  A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = FALSE)
  cache <- new.env()
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = FALSE, cache = cache), A)
  expect_identical(sum(!is.na(cache$nearest)), sum(A$n == 1L))
  # The same rows in other regions: the locations already found are not
  # searched again.
  gi2 <- .pt_regions(L$xy, 1200L, 9L); K2 <- max(gi2)
  want <- .pool_fit(L$y, gi2, K2, L$xy, use_fnn = FALSE)
  search <- .nearest_location_search
  asked <- integer(0)
  local_mocked_bindings(.nearest_location_search = function(u, q, ...) {
    asked <<- c(asked, q)
    search(u, q, ...)
  })
  known <- which(!is.na(cache$nearest))
  expect_identical(.pool_fit(L$y, gi2, K2, L$xy, use_fnn = FALSE, cache = cache), want)
  expect_length(intersect(asked, known), 0L)
  expect_identical(anyDuplicated(asked), 0L)
})


# ---- 15. A response without structure ----------------------------------------

test_that("15. the rule: regions beat one mean by a margin, and by more than chance", {
  # FALSE when one mean is within `tol` of the regions, or better.
  expect_false(.no_structure(1, 1.02, 1e-6, tol = 0.02)$structured)
  expect_false(.no_structure(1, 0.9, 1e-6, tol = 0.02)$structured)
  expect_false(.no_structure(1, 1.019, 0, tol = 0.02)$structured)
  expect_true(.no_structure(1, 1.021, 1e-6, tol = 0.02)$structured)
  expect_false(.no_structure(1, 1.021, 1e-6, tol = 0.05)$structured)
  # FALSE when the F test is not below 5 percent after the adjustment.
  expect_false(.no_structure(1, 2, 0.05, tol = 0.02)$structured)
  expect_true(.no_structure(1, 2, 0.0499, tol = 0.02)$structured)
  expect_false(.no_structure(1, 2, 0.01, tol = 0.02, n_searched = 5L)$structured)
  expect_true(.no_structure(1, 2, 0.0099, tol = 0.02, n_searched = 5L)$structured)
  expect_false(.no_structure(1, 2, 0.9, tol = 0.02, n_searched = 7L)$structured)
  # What it returns.
  r <- .no_structure(1, 2, 0.004, tol = 0.02, n_searched = 5L)
  expect_identical(names(r), c("structured", "structure_p", "beats"))
  expect_identical(r$structure_p, 0.02)
  expect_true(r$beats)
  expect_identical(.no_structure(1, 2, 0.4, tol = 0.02, n_searched = 5L)$structure_p, 1)
  expect_false(.no_structure(1, 1, 0.001, tol = 0.02)$beats)
  # NA with a missing error; with a missing p-value NA when the regions beat
  # one mean, and FALSE when they do not: the first condition has failed
  # whatever the test would have said.
  expect_identical(.no_structure(NA_real_, 2, 0.001, tol = 0.02)$structured, NA)
  expect_identical(.no_structure(1, NA_real_, 0.001, tol = 0.02)$structured, NA)
  expect_identical(.no_structure(Inf, 2, 0.001, tol = 0.02)$structured, NA)
  expect_identical(.no_structure(1, 2, NA_real_, tol = 0.02)$structured, NA)
  expect_identical(.no_structure(1, 2, NA_real_, tol = 0.02)$structure_p, NA_real_)
  expect_identical(.no_structure(1, 1, NA_real_, tol = 0.02)$structured, FALSE)
  expect_identical(.no_structure(1, 2, numeric(0), tol = 0.02)$structured, NA)
  expect_identical(.STRUCTURE_TOL, 0.02)
})

test_that("15. white noise is called structureless and a field with noise is not", {
  # 200 layers of each kind: 300 locations uniform on a 1e5 square, x then y,
  # then the response, then regions of the density kind at 67 centres (4.5
  # rows a region) from 3 k-means++ starts under set.seed(123).
  one <- function(seed, field) {
    .pt_seed(seed)
    xy <- cbind(stats::runif(300, 0, 1e5), stats::runif(300, 0, 1e5))
    y  <- if (field) sin(xy[, 1] / 2e4) + cos(xy[, 2] / 1.5e4) + stats::rnorm(300, 0, 1.5)
          else stats::rnorm(300)
    set.seed(123)
    kb <- .partition_centres(xy, .count_at(300, 4.5), partition = "density",
                             init = "kmeans++", nstart = 3L, iter.max = 100L)
    g  <- .nearest_centre(xy[, 1], xy[, 2], kb$centers[, 1], kb$centers[, 2])
    gi <- match(g, sort(unique(g)))
    f  <- .pool_fit(y, gi, max(gi), xy, want = "pooled")
    .no_structure(f$loo_pooled, f$loo_one, f$f_p, tol = 0.02, n_searched = 1L)$structured
  }
  expect_identical(.count_at(300, 4.5), 67L)
  noise <- vapply(1:200, function(r) one(4000 + r, FALSE), NA)
  field <- vapply(1:200, function(r) one(5000 + r, TRUE), NA)
  expect_false(anyNA(noise))
  expect_false(anyNA(field))
  # One run of these seeds: FALSE on 195 of the noise layers, TRUE on 200 of
  # the field layers.
  expect_gte(sum(!noise), 185L)
  expect_gte(sum(field), 190L)
})


# ---- 16. The bound on exact searches -----------------------------------------

test_that("16. an exact search above the limit stops with spatialkit_work_limit", {
  cen <- local({ .pt_seed(16); cbind(stats::runif(20, 0, 100), stats::runif(20, 0, 100)) })
  withr::local_options(spatialkit.exact_pairs_max = 10)
  e <- tryCatch(.region_neighbours(cen, 8L, FALSE), error = function(e) e)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_s3_class(e, "error")
  expect_identical(e$work, 400)
  expect_identical(e$limit, 10)
  expect_identical(e$what, "cell neighbours")
  expect_null(conditionCall(e))
  expect_identical(
    conditionMessage(e),
    paste0("resolution_profile(): finding each cell's neighbours among 20 cells would compare ",
           "4e+02 pairs without the FNN package, above the limit of 10 ",
           "(options(spatialkit.exact_pairs_max = )). Install FNN: the result is the same."))
  # The caller's name is the one handed down.
  e2 <- tryCatch(.region_neighbours(cen, 8L, FALSE, caller = "cell_values_at"),
                 error = function(e) e)
  expect_match(conditionMessage(e2), "^cell_values_at\\(\\): ")
  # Nothing at the limit or below it.
  withr::local_options(spatialkit.exact_pairs_max = 400)
  expect_identical(.region_neighbours(cen, 8L, FALSE), .pt_neighbours(cen, 8L))
  withr::local_options(spatialkit.exact_pairs_max = 399)
  expect_error(.region_neighbours(cen, 8L, FALSE), class = "spatialkit_work_limit")
})

test_that("16. the guard itself, the fit and the nearest other row", {
  expect_null(.work_guard(10, 10, "a search", "f"))
  expect_null(.work_guard(0, 0, "a search", "f"))
  e <- tryCatch(.work_guard(11, 10, "spacing", "cell_values_at",
                            search = "the spacing of 4 locations",
                            data = list(locations = 4L)),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$what, "spacing")
  expect_identical(e$locations, 4L)
  expect_identical(e$work, 11)
  expect_identical(e$limit, 10)
  expect_identical(
    conditionMessage(e),
    paste0("cell_values_at(): the spacing of 4 locations would compare 11 pairs without the ",
           "FNN package, above the limit of 10 (options(spatialkit.exact_pairs_max = )). ",
           "Install FNN: the result is the same."))
  expect_match(conditionMessage(tryCatch(.work_guard(2, 1, "x", "f", remedy = "Do less."),
                                         error = function(e) e)), "Do less\\.$")

  # The default limit, and an option that is not a number of pairs.
  expect_identical(.EXACT_PAIRS_MAX, 2e9)
  expect_identical(.FIT_WORK_MAX, 1e10)
  withr::local_options(spatialkit.exact_pairs_max = NULL)
  expect_identical(.exact_pairs_max(), 2e9)
  withr::local_options(spatialkit.exact_pairs_max = 5e3L)
  expect_identical(.exact_pairs_max(), 5000)
  # A limit of 0 is a limit: no exact search of even one pair.
  withr::local_options(spatialkit.exact_pairs_max = 0)
  expect_identical(.exact_pairs_max(), 0)
  expect_error(.region_neighbours(cbind(c(0, 1, 3), 0), 8L, FALSE), class = "spatialkit_work_limit")
  for (bad in list("100", NA_real_, -1, c(1, 2))) {
    withr::local_options(spatialkit.exact_pairs_max = bad)
    expect_error(.exact_pairs_max(), "spatialkit.exact_pairs_max", fixed = TRUE)
  }

  # The fit passes the error on: its neighbours first, then, for plain means,
  # the nearest other row of the rows alone in their region.
  L <- .pt_layout(3L)                              # 600 rows, 300 regions, 107 of one row
  expect_identical(L$K, 300L)
  alone <- sum(tabulate(L$gi) == 1L)
  withr::local_options(spatialkit.exact_pairs_max = NULL)
  geo <- .pool_geometry(L$gi, L$K, L$xy, tabulate(L$gi), FALSE)
  withr::local_options(spatialkit.exact_pairs_max = 300^2 - 1)
  e <- tryCatch(.pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = FALSE), error = function(e) e)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, 300^2)
  expect_identical(e$what, "cell neighbours")
  withr::local_options(spatialkit.exact_pairs_max = 300^2)
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = "pooled", use_fnn = FALSE)$status, "ok")
  # The rows alone in their region against every location: fewer pairs than
  # the neighbours of the regions took, so the limit is lowered to show it.
  withr::local_options(spatialkit.exact_pairs_max = alone * 600 - 1)
  e <- tryCatch(.pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", geo = geo, use_fnn = FALSE),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, alone * 600)
  expect_identical(e$what, "nearest other location")
  expect_match(conditionMessage(e),
               sprintf("finding the nearest other location of %d of 600 locations", alone))
  # Pooled values need no such search, so the same call for them returns;
  # and at the limit the search for plain means runs.
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = "pooled", geo = geo,
                             use_fnn = FALSE)$status, "ok")
  withr::local_options(spatialkit.exact_pairs_max = alone * 600)
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", geo = geo,
                             use_fnn = FALSE)$status, "ok")
})

test_that("16. the pairs are counted as a double, past the largest integer", {
  # 46,341 regions is the first count whose square is no integer in R, and
  # just above what the default limit allows; a product taken in integers
  # would be NA there and the stop a plain error about a missing value.  The
  # guard stops before any search, so these take no time.
  withr::local_options(spatialkit.exact_pairs_max = NULL)
  cen <- local({ .pt_seed(161); cbind(stats::runif(46341), stats::runif(46341)) })
  e <- tryCatch(.region_neighbours(cen, 8L, FALSE), error = function(e) e)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, 46341^2)
  expect_gt(e$work, .Machine$integer.max)
  expect_identical(e$limit, 2e9)
  expect_identical(e$what, "cell neighbours")
  # The same for the nearest other row of 50,000 rows, each alone at its
  # location.
  xy <- local({ .pt_seed(162); cbind(stats::runif(50000), stats::runif(50000)) })
  e <- tryCatch(.nearest_other(xy, seq_len(50000), seq_len(50000), FALSE), error = function(e) e)
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$work, 50000^2)
  expect_identical(e$what, "nearest other location")
  # In a column of an assigned layer that stop is a status.
  y <- local({ .pt_seed(163); stats::rnorm(2L * 46341L) })
  P <- .pool_column(y, rep(seq_len(46341L), 2L), rbind(cen, cen + 1e-3), NULL, "pooled", FALSE)
  expect_identical(P$record$status, "work limit")
  expect_identical(P$work, 46341^2)
  expect_identical(P$limit, 2e9)
})

test_that("the engine's constants", {
  expect_identical(.POOL_NEIGHBOURS, 8L)
  expect_identical(.OBS_PER_CELL_DEFAULT, 4.5)
  expect_identical(.OBS_PER_CELL_LADDER, c(18, 12, 9, 6, 4.5, 3, 2))
  expect_identical(.RECORD_TOL, 1e-8)
  expect_identical(.FEW_RESIDUALS, 5)
  expect_identical(.POOL_MIN_DF, 30L)
  expect_identical(.STRUCTURE_TOL, 0.02)
  expect_identical(.EXACT_PAIRS_MAX, 2e9)
  expect_identical(.FIT_WORK_MAX, 1e10)
})

test_that("16. with FNN nothing is bounded", {
  skip_if_not_installed("FNN")
  L <- .pt_layout(3L)
  want <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = TRUE)
  withr::local_options(spatialkit.exact_pairs_max = 1)
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = TRUE), want)
})

test_that("16. .has_backend() answers for a package, and one mock switches FNN off", {
  expect_true(.has_backend("stats"))
  expect_false(.has_backend("a.package.that.does.not.exist"))
  expect_identical(.has_backend("FNN"), requireNamespace("FNN", quietly = TRUE))
  expect_identical(formals(.pool_fit)$use_fnn, quote(.has_backend("FNN")))
  # With the answer FALSE the fit takes the exact path: FNN is not reached.
  L <- .pt_layout(2L)
  A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = FALSE)
  local_mocked_bindings(.has_backend = function(pkg) FALSE,
                        .knn_tied_rows = function(...) stop("the kd-tree path was taken"))
  withr::local_options(spatialkit.exact_pairs_max = 10)
  expect_error(.pool_fit(L$y, L$gi, L$K, L$xy), class = "spatialkit_work_limit")
  withr::local_options(spatialkit.exact_pairs_max = NULL)
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy), A)
})


# ---- 17. want ----------------------------------------------------------------

test_that("17. `want` decides which residuals are computed and nothing else", {
  L <- .pt_layout(7L)                              # 3,000 rows at 2 a region
  for (f in .pt_paths()) {
    tag <- .pt_path_tag(f)
    both   <- .pool_fit(L$y, L$gi, L$K, L$xy, want = c("pooled", "mean"), use_fnn = f)
    pooled <- .pool_fit(L$y, L$gi, L$K, L$xy, want = "pooled", use_fnn = f)
    plain  <- .pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", use_fnn = f)
    shared <- c("n", "S", "m", "N", "K", "sigma2", "cen", "nb", "target", "Q", "tau2", "w",
                "value", "loo_one", "loo_one_se", "f_p", "constant", "flat", "status")
    expect_identical(pooled[shared], both[shared], info = tag)
    expect_identical(plain[shared], both[shared], info = tag)
    # The default is both.
    expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f), both, info = tag)
    expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = c("mean", "pooled"), use_fnn = f),
                     both, info = tag)
    # The one asked for is the one of the full fit.
    mine_p <- c("resid_pooled", "loo_pooled", "loo_pooled_se")
    mine_m <- c("resid_mean", "loo_mean", "loo_mean_se")
    expect_identical(pooled[mine_p], both[mine_p], info = tag)
    expect_identical(plain[mine_m], both[mine_m], info = tag)
    expect_true(is.finite(both$loo_pooled) && is.finite(both$loo_mean), info = tag)
    # The one not asked for: no vector, and NA for its two numbers.
    expect_null(pooled$resid_mean)
    expect_identical(pooled$loo_mean, NA_real_, info = tag)
    expect_identical(pooled$loo_mean_se, NA_real_, info = tag)
    expect_null(plain$resid_pooled)
    expect_identical(plain$loo_pooled, NA_real_, info = tag)
    expect_identical(plain$loo_pooled_se, NA_real_, info = tag)
    # The list has the same elements in the same order whatever was asked.
    expect_identical(names(pooled), names(both), info = tag)
    expect_identical(names(plain), names(both), info = tag)
    expect_identical(
      names(both),
      c("n", "S", "m", "N", "K", "sigma2", "cen", "nb", "target", "Q", "tau2", "w", "value",
        "resid_pooled", "resid_mean", "loo_pooled", "loo_pooled_se", "loo_mean",
        "loo_mean_se", "loo_one", "loo_one_se", "f_p", "constant", "flat", "status"),
      info = tag)
  }
})

test_that("17. pooled values alone never look for a nearest other row", {
  L <- .pt_layout(7L)
  expect_gt(sum(tabulate(L$gi) == 1L), 500L)
  want <- .pool_fit(L$y, L$gi, L$K, L$xy, want = "pooled", use_fnn = FALSE)
  local_mocked_bindings(.nearest_other = function(...) stop("the nearest other row was searched for"))
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = "pooled", use_fnn = FALSE), want)
  expect_error(.pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", use_fnn = FALSE),
               "nearest other row was searched for")
  # Nor does a fit handed every row's nearest other row, nor one of a
  # constant column, whose residuals are set.
  nn1 <- .pt_nearest_other(L$xy, L$gi)
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", nn1 = nn1,
                             use_fnn = FALSE)$status, "ok")
  expect_identical(.pool_fit(rep(0.1, 3000), L$gi, L$K, L$xy, use_fnn = FALSE)$status,
                   "no spread inside any cell")
  # One region needs none either.
  expect_identical(.pool_fit(L$y, rep(1L, 3000), 1L, L$xy, use_fnn = FALSE)$status, "one cell")
})

test_that("17. a fit handed each row's nearest other row uses it", {
  L <- .pt_layout(3L)
  A <- .pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", use_fnn = FALSE)
  alone <- which(A$n[L$gi] == 1L)
  nn1 <- .pt_nearest_other(L$xy, L$gi)
  expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", nn1 = nn1, use_fnn = FALSE), A)
  # Another row handed in gives that row's region's mean.
  other <- nn1
  other[alone] <- ifelse(nn1[alone] == 1L, 2L, 1L)
  B <- .pool_fit(L$y, L$gi, L$K, L$xy, want = "mean", nn1 = other, use_fnn = FALSE)
  expect_identical(B$resid_mean[alone], L$y[alone] - A$m[L$gi[other[alone]]])
  expect_identical(B$resid_mean[-alone], A$resid_mean[-alone])
  expect_false(identical(B$loo_mean, A$loo_mean))
})


# ---- 18. geo, and column names in words --------------------------------------

test_that("18. a fit given the centroids and neighbours of an earlier fit computes neither again", {
  for (i in c(4L, 12L)) for (f in .pt_paths()) {
    L <- .pt_layout(i); tag <- sprintf("layout %d, %s", i, .pt_path_tag(f))
    A <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    geo <- A[c("cen", "nb")]
    other <- cos(L$xy[, 1] / 3e4) + L$y^2            # another column on the same rows
    C <- .pool_fit(other, L$gi, L$K, L$xy, use_fnn = f)
    local({
      local_mocked_bindings(
        .region_neighbours = function(...) stop("the neighbours were searched for"),
        .pool_geometry = function(...) stop("the centroids were summed"))
      expect_identical(.pool_fit(L$y, L$gi, L$K, L$xy, geo = geo, use_fnn = f), A, info = tag)
      expect_identical(.pool_fit(other, L$gi, L$K, L$xy, geo = geo, use_fnn = f), C, info = tag)
      expect_error(.pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f), "centroids were summed")
    })
    # The geometry is what the fit computes for itself.
    expect_identical(.pool_geometry(L$gi, L$K, L$xy, A$n, f), geo, info = tag)
    # It is read: neighbours of the right shape that are not these regions'
    # nearest give other targets.  With three regions there are no others.
    if (L$K > 9L) {
      further <- geo
      further$nb <- .region_neighbours(geo$cen, 9L, f)[, -1L]    # the second to the ninth nearest
      expect_false(identical(.pool_fit(L$y, L$gi, L$K, L$xy, geo = further, use_fnn = f)$target,
                             A$target), info = tag)
    }
  }
  # One region: its geometry is passed on like any other.
  y <- c(1, 2, 4); xy <- cbind(c(0, 1, 5), c(0, 2, 1))
  A <- .pool_fit(y, rep(1L, 3), 1L, xy, use_fnn = FALSE)
  local_mocked_bindings(.pool_geometry = function(...) stop("the centroids were summed"))
  expect_identical(.pool_fit(y, rep(1L, 3), 1L, xy, geo = A[c("cen", "nb")], use_fnn = FALSE), A)
})

test_that("18. a fit refuses centroids and neighbours that are not of its regions' shape", {
  L <- .pt_layout(2L); K <- L$K                    # 600 rows in 133 regions
  A <- .pool_fit(L$y, L$gi, K, L$xy, use_fnn = FALSE)
  geo <- A[c("cen", "nb")]
  fit <- function(geo, ...) .pool_fit(L$y, L$gi, K, L$xy, geo = geo, use_fnn = FALSE, ...)
  expect_identical(fit(geo), A)
  # Those of other regions of the same rows, fewer and more.  With fewer,
  # S[nb] would be laid out in K rows without a complaint and every target
  # would be a sum over other regions.
  for (centres in c(60L, 250L)) {
    g2 <- .pt_regions(L$xy, centres, 5L); K2 <- max(g2)
    expect_false(K2 == K)
    geo2 <- .pool_geometry(g2, K2, L$xy, tabulate(g2), FALSE)
    expect_error(fit(geo2), "`geo` must hold the centroids and neighbours of these 133 regions")
    expect_error(fit(geo2, want = "pooled"), "`geo` must hold")
    # A column of an assigned layer passes the refusal on: it is no status.
    expect_error(.pool_column(L$y, L$gi, L$xy, NULL, "pooled", FALSE, geo = geo2), "`geo` must hold")
  }
  # Each part of the shape.
  bad <- function(change) { g <- geo; change(g) }
  expect_error(fit(bad(function(g) { g$nb[7L, 3L] <- 7L; g })), "`geo` must hold")       # itself
  expect_error(fit(bad(function(g) { g$nb[7L, 3L] <- K + 1L; g })), "`geo` must hold")   # no region
  expect_error(fit(bad(function(g) { g$nb[7L, 3L] <- 0L; g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$nb[7L, 3L] <- NA; g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$nb <- g$nb[, -1L]; g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$nb <- t(g$nb); g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$cen <- g$cen[-1L, ]; g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$cen <- cbind(g$cen, 0); g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$cen[2L, 1L] <- NA; g })), "`geo` must hold")
  expect_error(fit(bad(function(g) { g$nb <- NULL; g })), "`geo` must hold")
  expect_error(fit(geo$nb), "`geo` must hold")
  # What a fit without one returns when nothing was fitted is no geometry.
  none <- .pool_fit(L$y[1:5], 1:5, 5L, L$xy[1:5, ], use_fnn = FALSE)
  expect_identical(none$status, "no cell holds two values")
  expect_error(.pool_fit(L$y[1:5], 1:5, 5L, L$xy[1:5, ], geo = none[c("cen", "nb")],
                         use_fnn = FALSE), "`geo` must hold")
})

test_that("18. column names in the words of a warning", {
  expect_identical(.cols_text("a"), "'a'")
  expect_identical(.cols_text(c("a", "b")), "'a' and 'b'")
  expect_identical(.cols_text(c("a", "b", "c")), "'a', 'b' and 'c'")
  expect_identical(.cols_text(letters[1:6]), "'a', 'b', 'c', 'd', 'e' and 'f'")
  expect_identical(.cols_text(letters[1:7]), "'a', 'b', 'c', 'd', 'e', 'f' and 1 more")
  expect_identical(.cols_text(letters[1:10]), "'a', 'b', 'c', 'd', 'e', 'f' and 4 more")
  expect_identical(.cols_text(c("income", "share of renters")), "'income' and 'share of renters'")
  expect_identical(.cols_text(character(0)), "")
})


# ---- The interval's offsets ---------------------------------------------------

test_that("the offsets are quantiles at probabilities rounded to 12 places", {
  e <- local({ .pt_seed(19); stats::rnorm(937) })
  for (level in c(0.5, 0.8, 0.9, 0.95, 0.99)) {
    p <- round((1 - level) / 2, 12)
    expect_identical(.obs_offsets(e, level),
                     stats::quantile(e, c(p, 1 - p), names = FALSE, type = 7),
                     info = format(level))
  }
  expect_identical(.obs_offsets(e, 0.9), stats::quantile(e, c(0.05, 0.95), names = FALSE))
  expect_identical(.obs_offsets(e, 0.95), stats::quantile(e, c(0.025, 0.975), names = FALSE))
  expect_identical(.obs_offsets(e, 0.8), stats::quantile(e, c(0.1, 0.9), names = FALSE))
  # What the rounding is for: (1 - 0.9) / 2 is not 0.05 in floating point.
  expect_false((1 - 0.9) / 2 == 0.05)
  off <- .obs_offsets(e, 0.9)
  expect_lt(off[1L], 0)
  expect_gt(off[2L], 0)
  # A constant column's residuals are all 0 and the interval has no width.
  expect_identical(.obs_offsets(numeric(50), 0.9), c(0, 0))
  # Fewer than two residuals, or a missing one: no interval.
  expect_identical(.obs_offsets(1.5, 0.9), c(NA_real_, NA_real_))
  expect_identical(.obs_offsets(NULL, 0.9), c(NA_real_, NA_real_))
  expect_identical(.obs_offsets(c(1, NA, 3), 0.9), c(NA_real_, NA_real_))
  expect_identical(.obs_offsets(c(-1, 3), 0.5), c(0, 2))
})


# ---- One column of an assigned layer ------------------------------------------

# 3,000 rows at 4.5 a region, as an assigned layer would hand them over: an
# ID per row of the type asked for, in an order that is not the regions'.
.pt_ids <- function(gi, K, type) {
  code <- local({ .pt_seed(77); sample(100000L, K) })
  switch(type,
         integer   = code[gi],
         double    = (code / 8)[gi],
         # Upper and lower case, which the C locale and a language's own
         # order sort differently.
         character = paste0(c("a", "B", "c", "D")[code %% 4L + 1L], code)[gi],
         factor    = factor(sprintf("u%05d", code)[gi], levels = rev(sprintf("u%05d", sort(code)))))
}

test_that("a column's values are the fit's, by region ID of any type", {
  L <- .pt_layout(4L)
  for (f in .pt_paths()) for (type in c("integer", "double", "character", "factor")) {
    tag <- sprintf("%s IDs, %s", type, .pt_path_tag(f))
    id <- .pt_ids(L$gi, L$K, type)
    P  <- .pool_column(L$y, id, L$xy, obs_level = 0.9, want = c("pooled", "mean"), use_fnn = f)
    # The regions in the order of the IDs: radix order.
    u <- unique(id); u <- u[order(u, method = "radix")]
    expect_identical(P$ids, u, info = tag)
    if (type == "character") {
      expect_identical(P$ids, withr::with_collate("C", sort(unique(id))), info = tag)
      # The same in a session whose locale sorts strings another way.
      for (loc in c("en_US.UTF-8", "en_US", "English_United States.1252")) {
        other <- suppressWarnings(tryCatch(
          withr::with_collate(loc, list(
            sorts = sort(c("B", "a")),
            ids = .pool_column(L$y, id, L$xy, obs_level = 0.9, want = c("pooled", "mean"),
                               use_fnn = f)$ids)),
          error = function(e) NULL))
        if (is.null(other) || identical(other$sorts, c("B", "a"))) next
        expect_identical(other$ids, P$ids, info = paste(tag, loc))
        expect_false(identical(withr::with_collate(loc, sort(unique(id))), P$ids),
                     info = paste(tag, loc))
        break
      }
    }
    if (type == "factor") expect_identical(as.character(P$ids),
                                           levels(id)[levels(id) %in% id], info = tag)
    gi <- match(id, u)
    A  <- .pool_fit(L$y, gi, L$K, L$xy, use_fnn = f)
    expect_identical(P$fit, A, info = tag)
    expect_identical(P$value, A$value, info = tag)
    expect_identical(P$w, A$w, info = tag)
    expect_identical(P$target, A$target, info = tag)
    expect_identical(P$n, A$n, info = tag)
    expect_identical(P$geo, A[c("cen", "nb")], info = tag)
    expect_identical(P$mask, rep(TRUE, 3000), info = tag)
    expect_identical(P$work, NA_real_, info = tag)
    expect_identical(P$limit, NA_real_, info = tag)
    off <- .obs_offsets(A$resid_pooled, 0.9)
    expect_identical(P$obs_lo, A$value + off[1L], info = tag)
    expect_identical(P$obs_hi, A$value + off[2L], info = tag)
    # Region by region these are the values of the regions numbered as the
    # layout numbers them, to rounding.
    B <- .pool_fit(L$y, L$gi, L$K, L$xy, use_fnn = f)
    at <- match(id[match(seq_len(L$K), L$gi)], P$ids)
    expect_lt(.pt_gap(P$value[at], B$value), 1e-12)
    expect_lt(.pt_gap(P$w[at], B$w), 1e-12)
    expect_lt(.pt_gap(P$target[at], B$target), 1e-12)

    # The column's entry in the record.
    r <- P$record
    expect_identical(names(r), c("column", "role", "n", "cells", "df", "neighbours", "sigma2",
                                 "tau2", "own_weight", "share_alone", "loo", "loo_se",
                                 "loo_one_cell", "structured", "structure_p", "offsets",
                                 "repeated", "skipped", "constant", "status"), info = tag)
    expect_identical(r$column, NA_character_, info = tag)
    expect_identical(r$role, NA_character_, info = tag)
    expect_identical(r$n, 3000L, info = tag)
    expect_identical(r$cells, L$K, info = tag)
    expect_identical(r$df, 3000L - L$K, info = tag)
    expect_identical(r$neighbours, 8L, info = tag)
    expect_identical(r$sigma2, A$sigma2, info = tag)
    expect_identical(r$tau2, A$tau2, info = tag)
    expect_identical(r$own_weight, stats::median(A$w), info = tag)
    expect_identical(r$share_alone, sum(A$n == 1L) / 3000, info = tag)
    expect_identical(r$loo, A$loo_pooled, info = tag)
    expect_identical(r$loo_se, A$loo_pooled_se, info = tag)
    expect_identical(r$loo_one_cell, A$loo_one, info = tag)
    s <- .no_structure(A$loo_pooled, A$loo_one, A$f_p, tol = 0.02)
    expect_true(r$structured, info = tag)
    expect_identical(r$structured, s$structured, info = tag)
    expect_identical(r$structure_p, s$structure_p, info = tag)
    expect_identical(r$offsets, off, info = tag)
    expect_identical(r$repeated, 0L, info = tag)
    expect_identical(r$skipped, 0L, info = tag)
    expect_false(r$constant, info = tag)
    expect_identical(r$status, "ok", info = tag)
  }
  # Without a level there is no interval.
  P <- .pool_column(L$y, L$gi, L$xy, obs_level = NULL, want = "pooled", use_fnn = FALSE)
  expect_null(P$obs_lo)
  expect_null(P$obs_hi)
  expect_null(P$record$offsets)
  expect_true("offsets" %in% names(P$record))
  expect_null(P$fit$resid_mean)
  expect_identical(P$record$neighbours, 8L)
})

test_that("rows without a value, a region or a location are as if absent", {
  L <- .pt_layout(2L); n <- length(L$y)
  y <- L$y; id <- sprintf("r%03d", L$gi); xy <- L$xy
  y[c(3, 50, 51)] <- NA
  y[77] <- Inf
  y[90] <- NaN
  id[c(10, 200)] <- NA
  xy[c(20, 300), 1] <- NA
  xy[400, 2] <- Inf
  xy[c(3, 10), 2] <- NA                  # also without a value, or without a region
  gone <- c(3, 50, 51, 77, 90, 10, 200, 20, 300, 400)
  for (f in .pt_paths()) {
    P <- .pool_column(y, id, xy, obs_level = 0.9, want = c("pooled", "mean"), use_fnn = f)
    Q <- .pool_column(y[-gone], id[-gone], xy[-gone, ], obs_level = 0.9,
                      want = c("pooled", "mean"), use_fnn = f)
    expect_identical(which(!P$mask), as.integer(sort(gone)))
    expect_identical(P[c("ids", "n", "value", "w", "target", "obs_lo", "obs_hi", "geo", "fit")],
                     Q[c("ids", "n", "value", "w", "target", "obs_lo", "obs_hi", "geo", "fit")])
    expect_identical(P$record[names(P$record) != "skipped"],
                     Q$record[names(Q$record) != "skipped"])
    # Rows with a value and a region and no location: 20, 300 and 400.
    expect_identical(P$record$skipped, 3L)
    expect_identical(Q$record$skipped, 0L)
    expect_identical(P$record$n, n - 10L)
    expect_identical(Q$mask, rep(TRUE, n - 10L))
  }
})

test_that("each case in which a column has no pooled value, or no pooling, has its status", {
  L <- .pt_layout(2L); n <- length(L$y); K <- L$K
  col <- function(y = L$y, id = L$gi, xy = L$xy, obs_level = 0.9, want = "pooled",
                  use_fnn = FALSE, ...)
    .pool_column(y, id, xy, obs_level = obs_level, want = want, use_fnn = use_fnn, ...)

  # No coordinates at all.
  P <- col(xy = NULL)
  expect_identical(P$record$status, "no coordinates")
  expect_identical(P$record$n, 0L)
  expect_identical(P$record$cells, 0L)
  expect_identical(P$record$skipped, n)
  expect_length(P$ids, 0L)
  expect_identical(P$value, numeric(0))
  expect_identical(P$obs_lo, numeric(0))
  expect_identical(P$mask, logical(n))
  expect_null(P$fit); expect_null(P$geo)
  expect_identical(P$record$offsets, c(NA_real_, NA_real_))

  # Coordinates, and no row with a value, a region and a location.
  for (P in list(col(y = rep(NA_real_, n)), col(id = rep(NA_integer_, n)),
                 col(xy = matrix(NA_real_, n, 2)))) {
    expect_identical(P$record$status, "no located value")
    expect_identical(P$record$cells, 0L)
    expect_length(P$ids, 0L)
    expect_identical(P$mask, logical(n))
    expect_null(P$fit)
  }
  expect_identical(col(xy = matrix(NA_real_, n, 2))$record$skipped, n)
  expect_identical(col(y = rep(NA_real_, n))$record$skipped, 0L)

  # One region.
  P <- col(id = rep("all", n))
  expect_identical(P$record$status, "one cell")
  expect_identical(P$ids, "all")
  expect_equal(P$value, mean(L$y), tolerance = 1e-14)
  expect_identical(P$value, P$fit$m)
  expect_identical(P$w, 1)
  expect_identical(P$target, NA_real_)
  expect_identical(P$record$neighbours, 0L)
  expect_identical(P$record$tau2, NA_real_)
  expect_identical(P$record$structured, NA)
  expect_identical(P$record$own_weight, 1)
  expect_identical(P$record$offsets, .obs_offsets((L$y - P$fit$m) * n / (n - 1), 0.9))
  expect_identical(P$obs_hi - P$obs_lo, diff(P$record$offsets))
  # One row in all.
  P <- col(y = 3.5, id = 9L, xy = cbind(1, 2))
  expect_identical(P$record$status, "one cell")
  expect_identical(P$value, 3.5)
  expect_identical(P$w, 1)
  expect_identical(P$record$sigma2, NA_real_)
  expect_identical(P$record$loo, NA_real_)
  expect_identical(c(P$obs_lo, P$obs_hi), c(NA_real_, NA_real_))
  expect_identical(P$record$share_alone, 1)

  # No region holds two rows.
  P <- col(y = L$y[1:12], id = letters[12:1], xy = L$xy[1:12, ])
  expect_identical(P$record$status, "no cell holds two values")
  expect_identical(P$ids, letters[1:12])
  expect_identical(P$n, rep(1L, 12))
  expect_identical(P$value, rep(NA_real_, 12))
  expect_identical(P$w, rep(NA_real_, 12))
  expect_identical(P$obs_lo, rep(NA_real_, 12))
  expect_identical(P$record$df, 0L)
  expect_identical(P$record$share_alone, 1)
  expect_identical(P$record$sigma2, NA_real_)
  expect_identical(P$record$own_weight, NA_real_)
  expect_null(P$geo)
  expect_false(P$record$constant)

  # The work limit is a status here, with the pairs and the limit kept.
  withr::with_options(list(spatialkit.exact_pairs_max = 100), {
    P <- col()
    expect_identical(P$record$status, "work limit")
    expect_identical(P$work, as.numeric(K)^2)
    expect_identical(P$limit, 100)
    expect_identical(P$value, rep(NA_real_, K))
    expect_identical(P$w, rep(NA_real_, K))
    expect_identical(P$obs_hi, rep(NA_real_, K))
    expect_identical(P$ids, seq_len(K))
    expect_identical(P$record$cells, K)
    expect_identical(P$record$n, n)
    expect_identical(P$record$share_alone, sum(tabulate(L$gi) == 1L) / n)
    expect_null(P$fit); expect_null(P$geo)
    expect_identical(P$mask, rep(TRUE, n))
    # A flat column over the limit is "work limit" too: the search comes first.
    expect_identical(col(y = rep(0.1, n))$record$status, "work limit")
    # Fewer regions than the limit allows, and a search that is not needed.
    expect_true(col(id = (L$gi %% 9L) + 1L)$record$status %in%
                  c("ok", "no signal between cells"))
    expect_identical(col(id = rep(1L, n))$record$status, "one cell")
  })
  # With plain means asked for, the search for the nearest other row of the
  # rows alone in their region is bounded the same way.
  withr::with_options(list(spatialkit.exact_pairs_max = 100), {
    geo2 <- withr::with_options(list(spatialkit.exact_pairs_max = NULL),
                                .pool_geometry(L$gi, K, L$xy, tabulate(L$gi), FALSE))
    P <- col(want = "mean", geo = geo2)
    expect_identical(P$record$status, "work limit")
    expect_identical(P$work, sum(tabulate(L$gi) == 1L) * as.numeric(n))
    expect_identical(col(want = "pooled", geo = geo2)$record$status, "ok")
  })

  # No spread inside any region; no signal between regions.
  P <- col(y = rep(52000.37, n))
  expect_identical(P$record$status, "no spread inside any cell")
  expect_true(P$record$constant)
  expect_identical(P$w, rep(1, K))
  expect_identical(P$record$offsets, c(0, 0))
  expect_identical(P$obs_lo, P$value)
  expect_identical(P$obs_hi, P$value)
  expect_identical(P$record$loo, NA_real_)
  expect_identical(P$record$structured, NA)
  P <- col(y = 0.1 * L$gi)
  expect_identical(P$record$status, "no spread inside any cell")
  expect_false(P$record$constant)
  expect_identical(P$record$structured, NA)
  expect_true(is.finite(P$record$loo))
  # Flat and no better than one mean: `structured` is read only where there
  # is spread inside the regions, so it stays NA and is not FALSE.
  few <- local({ .pt_seed(41); list(xy = cbind(stats::runif(40), stats::runif(40)),
                                    y = stats::rnorm(40)) })
  few$y[40] <- few$y[39]                       # the one region of two rows, both the same
  P <- col(y = few$y, id = c(1:39, 39L), xy = few$xy)
  expect_identical(P$record$status, "no spread inside any cell")
  expect_false(P$record$constant)
  expect_lte(P$fit$loo_one, P$fit$loo_pooled * 1.02)
  expect_identical(P$record$structured, NA)
  expect_identical(P$record$structure_p, NA_real_)
  expect_identical(P$record$df, 1L)
  noise <- local({ .pt_seed(903); stats::rnorm(n) })
  P <- col(y = noise)
  expect_identical(P$record$status, "no signal between cells")
  expect_identical(P$record$tau2, 0)
  expect_identical(P$record$own_weight, 0)
  expect_false(P$record$structured)
  expect_true(all(is.finite(P$obs_lo)))
})

test_that("a column's record counts the rows that share a location", {
  L <- .pt_layout(2L); n <- length(L$y)
  xy <- L$xy
  xy[5:7, ] <- xy[rep(4L, 3), ]            # four rows on one location
  xy[100, ] <- xy[99, ]                    # and a pair
  P <- .pool_column(L$y, L$gi, xy, obs_level = NULL, want = "pooled", use_fnn = FALSE)
  expect_identical(P$record$repeated, 6L)
  expect_identical(.repeated_rows(xy), 6L)
  expect_identical(.repeated_rows(L$xy), 0L)
  expect_identical(.repeated_rows(xy[c(1, 1, 1), ]), 3L)
  expect_identical(.repeated_rows(xy[1, , drop = FALSE]), 0L)
  # A row without a value at a repeated location is not a scored row.
  y <- L$y; y[100] <- NA
  expect_identical(.pool_column(y, L$gi, xy, NULL, "pooled", FALSE)$record$repeated, 4L)
  # Gridded locations, where many rows share one coordinate and not the
  # other: two rows at (1, 1) with a row at (1, 2) between them in the order
  # of the first coordinate alone.
  expect_identical(.repeated_rows(rbind(c(1, 1), c(1, 2), c(1, 1))), 2L)
  expect_identical(.repeated_rows(rbind(c(1, 1), c(2, 1), c(1, 2), c(2, 2))), 0L)
  grid <- unname(as.matrix(expand.grid(seq(0.1, by = 0.7, length.out = 20),
                                       seq(0.3, by = 0.7, length.out = 15))))
  twice <- local({ .pt_seed(88); sample(nrow(grid), 100L) })
  both  <- rbind(grid, grid[twice, ])
  both  <- both[local({ .pt_seed(89); sample(nrow(both)) }), ]
  expect_identical(.repeated_rows(grid), 0L)
  expect_identical(.repeated_rows(both), 200L)
  expect_identical(.repeated_rows(rbind(both, grid[twice[1:10], ])), 210L)
})

test_that("a column's `structured` is the rule's at the margin of 2 percent, for one level", {
  # What .pool_column() hands the rule, held by its arguments: on a layer
  # with a strong field, or with none, any margin gives the same answer.
  L <- .pt_layout(4L)
  rule <- .no_structure
  seen <- NULL
  local_mocked_bindings(.no_structure = function(loo, loo_one, p, tol, n_searched = 1L) {
    seen <<- list(loo = loo, loo_one = loo_one, p = p, tol = tol, n_searched = n_searched)
    rule(loo, loo_one, p, tol, n_searched)
  })
  P <- .pool_column(L$y, L$gi, L$xy, obs_level = NULL, want = "pooled", use_fnn = FALSE)
  expect_identical(seen$tol, 0.02)
  expect_identical(seen$tol, .STRUCTURE_TOL)
  expect_equal(seen$n_searched, 1)
  expect_identical(seen$loo, P$fit$loo_pooled)
  expect_identical(seen$loo_one, P$fit$loo_one)
  expect_identical(seen$p, P$fit$f_p)
  expect_identical(P$record$structured, rule(seen$loo, seen$loo_one, seen$p, 0.02, 1L)$structured)
  expect_identical(P$record$structure_p, min(1, P$fit$f_p))
  # A layer on which the margin decides: one mean's error is 3 percent above
  # the regions', and the F test is far from chance.
  f <- P$fit
  local_mocked_bindings(.pool_fit = function(...) {
    f$loo_one <- f$loo_pooled * 1.03
    f$f_p <- 1e-6
    f
  })
  expect_true(.pool_column(L$y, L$gi, L$xy, NULL, "pooled", FALSE)$record$structured)
  local_mocked_bindings(.pool_fit = function(...) {
    f$loo_one <- f$loo_pooled * 1.01
    f$f_p <- 1e-6
    f
  })
  expect_false(.pool_column(L$y, L$gi, L$xy, NULL, "pooled", FALSE)$record$structured)
})

test_that("a later column with the same rows takes the first one's centroids and neighbours", {
  L <- .pt_layout(4L)
  id <- sprintf("c%04d", L$gi)
  first <- .pool_column(L$y, id, L$xy, obs_level = 0.9, want = "pooled", use_fnn = FALSE)
  other <- cos(L$xy[, 2] / 2e4) + 0.3 * L$y
  alone <- .pool_column(other, id, L$xy, obs_level = 0.9, want = "pooled", use_fnn = FALSE)
  local_mocked_bindings(.region_neighbours = function(...) stop("the neighbours were searched for"))
  again <- .pool_column(other, id, L$xy, obs_level = 0.9, want = "pooled", use_fnn = FALSE,
                        geo = first$geo)
  expect_identical(again, alone)
  # Without them the search is made, and here it is the mock's stop: a plain
  # error is passed on, only the work limit becomes a status.
  expect_error(.pool_column(other, id, L$xy, obs_level = 0.9, want = "pooled", use_fnn = FALSE),
               "neighbours were searched for")
})
