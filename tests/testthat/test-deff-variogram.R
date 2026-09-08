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

test_that("the spherical and Gaussian families follow the textbook formulas", {
  # Only the exponential branch was ever evaluated; (h/r)^2 in place of the
  # spherical (h/r)^3 passed every test (mutation testing, pass 6).  A user
  # can hand summarize_by_cell(deff = "variogram") a gstat Sph or Gau fit
  # directly, and estimate_sac_range() itself falls back to Sph.
  h <- c(0, 25, 50, 75, 100, 150)
  f_sph <- cor_fn(vgm_df(nugget = 0.2, psill = 0.8, range = 100, model = "Sph"))
  u <- h / 100
  expect_equal(f_sph(h),
               ifelse(u >= 1, 0, 0.8 * (1 - 1.5 * u + 0.5 * u^3)),
               tolerance = 1e-12)
  expect_equal(f_sph(50), 0.25)                 # 0.8 * (1 - 0.75 + 0.0625)
  f_gau <- cor_fn(vgm_df(nugget = 0.2, psill = 0.8, range = 100, model = "Gau"))
  expect_equal(f_gau(h), 0.8 * exp(-(h / 100)^2), tolerance = 1e-12)
  f_exp <- cor_fn(vgm_df(nugget = 0.2, psill = 0.8, range = 100, model = "Exp"))
  expect_equal(f_exp(h), 0.8 * exp(-h / 100), tolerance = 1e-12)
  # The three families are distinct at every interior lag, so the checks
  # above cannot be satisfied by one formula standing in for another.
  inner <- h[h > 0 & h < 100]
  expect_true(all(abs(f_sph(inner) - f_exp(inner)) > 1e-3))
  expect_true(all(abs(f_gau(inner) - f_exp(inner)) > 1e-3))
  expect_true(all(abs(f_gau(inner) - f_sph(inner)) > 1e-3))
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


# ---------------------------------------------------------------------------
# Nested models and unsupported families (sixth pass, reviewer A's Low list)
# ---------------------------------------------------------------------------

test_that("every structured component counts, weighted by its partial sill", {
  # The function read the single largest component, so a user-built
  # Nug + Exp + Sph model gave 0.108 at h = 200 where the model implies 0.197.
  # Reference: rho(h) = 1 - gamma(h) / total sill, written out by hand from
  # gstat's parametrisation of each family.
  vm <- data.frame(model = c("Nug", "Exp", "Sph"),
                   psill = c(0.2, 0.5, 0.3), range = c(0, 100, 500),
                   stringsAsFactors = FALSE)
  f  <- cor_fn(vm)
  expect_false(is.null(f))
  h  <- c(0, 1, 10, 50, 100, 200, 400, 499, 500, 800)
  sph <- ifelse(h >= 500, 0, 1 - 1.5 * (h / 500) + 0.5 * (h / 500)^3)
  ref <- (0.5 * exp(-h / 100) + 0.3 * sph) / 1.0
  expect_equal(f(h), ref, tolerance = 1e-12)
  expect_equal(f(200), 0.1973, tolerance = 1e-3)
  # A single-component model is unchanged by the generalisation.
  f1 <- cor_fn(vgm_df(nugget = 0.2, psill = 0.8, range = 50))
  expect_equal(f1(c(0, 25, 100)), 0.8 * exp(-c(0, 25, 100) / 50))
})

test_that("the nested-model correlation agrees with gstat::variogramLine()", {
  skip_if_not_installed("gstat")
  vm <- gstat::vgm(psill = 0.3, model = "Sph", range = 500,
                   add.to = gstat::vgm(psill = 0.5, model = "Exp", range = 100,
                                       nugget = 0.2))
  f <- cor_fn(vm)
  h <- c(1, 5, 10, 25, 50, 75, 100, 150, 200, 300, 500, 800, 1200)
  g <- gstat::variogramLine(vm, dist_vector = h)$gamma
  expect_equal(f(h), 1 - g / sum(vm$psill), tolerance = 1e-10)
})

test_that("a variogram family the function does not implement is refused, not read as exponential", {
  vm <- data.frame(model = c("Nug", "Mat"), psill = c(0.2, 0.8),
                   range = c(0, 100), stringsAsFactors = FALSE)
  expect_null(cor_fn(vm))
  vm2 <- data.frame(model = c("Nug", "Exp", "Pow"), psill = c(0.2, 0.5, 0.3),
                    range = c(0, 100, 1), stringsAsFactors = FALSE)
  expect_null(cor_fn(vm2))
  # End to end: summarize_by_cell() says which family it cannot use and
  # falls back to deff = 1 -- the naive SE, no deff_applied attribute.
  set.seed(77)
  pts <- sf::st_as_sf(
    data.frame(x = rep(c(0, 300, 600), each = 20) + runif(60, 0, 60),
               y = runif(60, 0, 60), z = rnorm(60),
               poly_id = rep(1:3, each = 20)),
    coords = c("x", "y"), crs = 32632)
  sac <- structure(300, class = c("sac_range", "numeric"),
                   variogram_model = vm, crs = sf::st_crs(pts))
  expect_warning(
    out <- summarize_by_cell(pts, "z", deff = "variogram", sac = sac),
    "supports exponential, spherical and Gaussian variogram models.*'Mat'")
  expect_null(attr(out, "deff_applied"))
  naive <- summarize_by_cell(pts, "z", deff = 1)
  expect_equal(out[["..se_resp_z"]], naive[["..se_resp_z"]])
})
