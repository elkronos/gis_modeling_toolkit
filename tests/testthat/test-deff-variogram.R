# tests/testthat/test-deff-variogram.R
# ---------------------------------------------------------------------------
# deff = "variogram" in summarize_by_cell().
#
# Kish's 1 + (n-1)*rho assumes every pair of points in a cell is equally
# correlated no matter how far apart they are.  That degrades as cells grow.
# The variogram form uses deff = sum(R)/n with R from the fitted correlation
# function, which is the exact effective-sample-size result and reduces to
# Kish when correlation is constant.
# ---------------------------------------------------------------------------

cor_fn  <- function(...) spatialkit:::.vgm_correlation_fn(...)
cell_de <- function(...) spatialkit:::.cell_deff_variogram(...)

# A fitted-variogram stand-in: gstat models are data frames of this shape.
vgm_df <- function(nugget = 0.2, psill = 0.8, range = 50, model = "Exp") {
  data.frame(model = c("Nug", model), psill = c(nugget, psill),
             range = c(0, range), stringsAsFactors = FALSE)
}


test_that("correlation decays from the sill ratio to zero", {
  f <- cor_fn(vgm_df(nugget = 0.2, psill = 0.8, range = 50))
  expect_false(is.null(f))

  # Correlation between two DISTINCT observations: the nugget discounts them
  # even at zero separation, because each carries independent nugget noise.
  # Self-correlation of 1 lives on the matrix diagonal, not in this function.
  expect_equal(f(0), 0.8)                     # psill / (nugget + psill)
  expect_equal(f(1e-9), 0.8, tolerance = 1e-6)
  expect_lt(f(1), 0.8)
  expect_gt(f(1), 0.7)
  expect_lt(f(500), 0.01)                     # far apart -> independent
  # monotone decreasing
  h <- seq(1, 300, length.out = 50)
  expect_true(all(diff(f(h)) < 0))
})

test_that("a pure-nugget model yields no correlation between distinct points", {
  f <- cor_fn(vgm_df(nugget = 1, psill = 1e-8, range = 50))
  expect_lt(f(0), 1e-6)      # even coincident distinct observations
  expect_lt(f(10), 1e-6)
})

test_that("unusable variogram input returns NULL rather than guessing", {
  expect_null(cor_fn(NULL))
  expect_null(cor_fn(data.frame(a = 1)))
  expect_null(cor_fn(vgm_df(psill = 0, nugget = 0)))
  expect_null(cor_fn(data.frame(model = "Nug", psill = 1, range = 0)))  # no structure
})

test_that("deff reduces to Kish when correlation is constant", {
  # The key property: with a constant off-diagonal rho, sum(R)/n is exactly
  # 1 + (n-1)*rho.  A flat correlation function makes that comparison direct.
  rho <- 0.4
  flat <- function(h) ifelse(h <= 0, 1, rho)

  set.seed(1)
  coords <- matrix(runif(30 * 2, 0, 100), ncol = 2)
  ids    <- rep(c("a", "b", "c"), each = 10)

  d <- cell_de(coords, ids, flat)
  expect_equal(unname(d[["a"]]), 1 + (10 - 1) * rho, tolerance = 1e-10)
  expect_equal(unname(d[["b"]]), 1 + (10 - 1) * rho, tolerance = 1e-10)
})

test_that("deff is bounded by 1 and n", {
  set.seed(2)
  coords <- matrix(runif(40 * 2, 0, 100), ncol = 2)
  ids    <- rep(c("a", "b"), each = 20)

  # Note the two stubs return different SHAPES: ifelse() inherits the matrix
  # dim of its input, rep() does not.  .cell_deff_variogram() must cope with
  # both -- it rebuilds the matrix rather than trusting cor_fn to keep `dim`.
  independent <- function(h) ifelse(h <= 0, 0, 0)   # nothing correlates
  redundant   <- function(h) rep(1, length(h))      # bare vector, no dim

  expect_equal(unname(cell_de(coords, ids, independent)[["a"]]), 1)
  expect_equal(unname(cell_de(coords, ids, redundant)[["a"]]), 20)
})

test_that("deff tolerates a correlation function that drops matrix dim", {
  # Regression guard: an earlier version called diag(R) <- 1 directly on
  # cor_fn()'s return value, which errored for any function returning a plain
  # vector ("only matrix diagonals can be replaced").
  set.seed(9)
  coords <- matrix(runif(12 * 2, 0, 50), ncol = 2)
  vec_fn <- function(h) rep(0.5, length(h))         # returns a bare vector
  d <- cell_de(coords, rep("a", 12), vec_fn)
  expect_true(is.finite(d[["a"]]))
  expect_equal(unname(d[["a"]]), 1 + (12 - 1) * 0.5, tolerance = 1e-10)
})

test_that("deff rises as correlation range grows", {
  # More correlation between the same points must mean less effective
  # information, hence a larger design effect.
  set.seed(3)
  coords <- matrix(runif(25 * 2, 0, 100), ncol = 2)
  ids    <- rep("a", 25)

  short <- cell_de(coords, ids, cor_fn(vgm_df(range = 2)))[["a"]]
  long  <- cell_de(coords, ids, cor_fn(vgm_df(range = 200)))[["a"]]
  expect_lt(short, long)
  expect_gte(short, 1)
  expect_lte(long, 25)
})

test_that("singleton and empty cells are handled", {
  coords <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  d <- cell_de(coords, c("a", "b"), cor_fn(vgm_df()))
  expect_equal(unname(d[["a"]]), 1)
  expect_equal(unname(d[["b"]]), 1)
})

test_that("a NULL correlation function yields all-NA rather than erroring", {
  set.seed(12)
  coords <- matrix(runif(10 * 2), ncol = 2)
  d <- cell_de(coords, rep("a", 10), NULL)
  expect_true(all(is.na(d)))
})

test_that("a subsampled cell reports ITS OWN design effect, not the subsample's", {
  # deff = sum(R)/n_i = 1 + (n_i - 1) * Rbar.  Subsampling estimates Rbar just
  # as well -- the subsample's pairwise-distance distribution is the cell's --
  # but sum(R)/n_used answers for a cell of size n_used.  The old code returned
  # the design effect of `max_n` points instead of the cell's, understating it
  # by roughly n_i/max_n: measured on 4000 points with an exponential
  # correlation of range 60, true deff 1821.8 and max_n = 500 returned 228.6.
  set.seed(4)
  n <- 300
  coords <- matrix(runif(n * 2, 0, 100), ncol = 2)
  f <- cor_fn(vgm_df(range = 30))

  full <- cell_de(coords, rep("a", n), f, max_n = n)[["a"]]
  sub  <- cell_de(coords, rep("a", n), f, max_n = 50L)[["a"]]

  expect_true(is.finite(sub))
  # Bounded by the CELL's size, which is the only meaningful ceiling.
  expect_lte(sub, n)
  # And close to the un-subsampled answer, which is the whole point.
  expect_equal(sub, full, tolerance = 0.15)
  # Regression guard: the old truncating form could not exceed max_n.
  expect_gt(sub, 50)
})

test_that("an un-subsampled cell is unchanged by the rescale", {
  # With n_used == n_i the new form is algebraically identical to sum(R)/n_i,
  # so the common case must be bit-for-bit what it always was.
  set.seed(9)
  n <- 40
  coords <- matrix(runif(n * 2, 0, 100), ncol = 2)
  f <- cor_fn(vgm_df(range = 30))
  d <- as.matrix(stats::dist(coords))
  R <- matrix(f(as.numeric(d)), n, n); diag(R) <- 1
  expect_equal(cell_de(coords, rep("a", n), f, max_n = n)[["a"]],
               min(max(sum(R) / n, 1), n), tolerance = 1e-10)
})
