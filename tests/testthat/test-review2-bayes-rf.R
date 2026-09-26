# Regression tests for the second review's Bayesian and random-forest
# findings.  Nothing here compiles a Stan model: brms::brm() and the
# posterior accessors are replaced by recorders, while everything in front of
# the sampler -- brms::get_prior(), validate_prior() and this package's own
# checks -- runs for real.  The real-sampler counterparts are in
# test-review2-bayes-brms.R, behind SPATIALKIT_TEST_BRMS.

.r2b_pts <- function(n = 60, seed = 1) {
  set.seed(seed)
  d <- sf::st_as_sf(
    data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
               a = rnorm(n)),
    coords = c("x", "y"), crs = 32632)
  d$z    <- 2 * d$a + rnorm(n)
  d$cat3 <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
                   levels = c("lo", "mid", "hi"))
  d$ord3 <- factor(d$cat3, ordered = TRUE)
  d$bin  <- as.integer(d$z > 0)
  d$pres <- factor(ifelse(d$bin == 1L, "present", "absent"))
  d
}

# fit_bayesian_spatial_model() with brms::brm() replaced by a recorder that
# returns `engine`.  The mock lives in this helper's frame, so it is gone by
# the time the fit is returned.
.r2b_capture_fit <- function(...,
                             engine = structure(list(), class = "r2b_stub"),
                             check_convergence = FALSE) {
  cap <- new.env()
  local_mocked_bindings(
    brm = function(...) { cap$args <- list(...); engine },
    .package = "brms")
  fit <- fit_bayesian_spatial_model(..., compute_loo = FALSE,
                                    check_convergence = check_convergence)
  list(fit = fit, args = cap$args)
}

.r2b_quiet <- function(expr) suppressMessages(expr)


# ---------------------------------------------------------------------------
# models MED-2: a user's gp_c sizes the derived gp_k
# ---------------------------------------------------------------------------

test_that(".gp_basis_spec() sizes k for the boundary factor it is given", {
  set.seed(1)
  xy <- scale(cbind(runif(200), runif(200)))
  b  <- gp_lengthscale_bounds(xy)
  def <- spatialkit:::.gp_basis_spec(xy, b)
  S    <- def$S
  r_lo <- b[["lower"]] / S
  rule <- function(cc) min(max(ceiling(1.75 * cc / r_lo), 10L), 50L)

  # The default is unchanged.
  expect_identical(def$k, as.integer(rule(def$c)))
  # A wider boundary gets the k the rule gives for IT, not the default's k.
  s3 <- spatialkit:::.gp_basis_spec(xy, b, c = 3)
  expect_identical(s3$c, 3)
  expect_identical(s3$k, as.integer(rule(3)))
  expect_gt(s3$k, def$k)
  expect_false(s3$capped)
  # ...up to the cap, which is then reported.
  s5 <- spatialkit:::.gp_basis_spec(xy, b, c = 5)
  expect_identical(s5$k, 50L)
  expect_true(s5$capped)
  # A narrower one gets fewer.
  s13 <- spatialkit:::.gp_basis_spec(xy, b, c = 1.3)
  expect_lte(s13$k, def$k)
})

test_that("fit_bayesian_spatial_model(gp_c = ) derives gp_k for that gp_c", {
  skip_if_not_installed("brms")
  d <- .r2b_pts(n = 200, seed = 3)
  base <- .r2b_quiet(.r2b_capture_fit(d, "z", "a"))$fit
  wide <- .r2b_quiet(.r2b_capture_fit(d, "z", "a", gp_c = 3))$fit
  expect_identical(wide$info$gp_c, 3)
  # The rule, from the numbers the fit itself records.
  r_lo <- wide$info$gp_lengthscale_bounds[["lower"]] / wide$info$gp_S
  expect_identical(wide$info$gp_k,
                   as.integer(min(max(ceiling(1.75 * 3 / r_lo), 10), 50)))
  expect_gt(wide$info$gp_k, base$info$gp_k)
  # So the resolvable scale stays where the default basis put it, instead of
  # coarsening in proportion to gp_c.
  expect_lt(wide$info$gp_ell_min, 1.1 * base$info$gp_ell_min)
  # An explicit gp_k still passes through untouched.
  fixed <- .r2b_quiet(.r2b_capture_fit(d, "z", "a", gp_c = 3, gp_k = 12))$fit
  expect_identical(fixed$info$gp_k, 12L)
  # And an invalid gp_c is refused before the rule reads it.
  expect_error(.r2b_capture_fit(d, "z", "a", gp_c = "3"), "`gp_c` must be")
  expect_error(.r2b_capture_fit(d, "z", "a", gp_c = 0.5), "`gp_c` must be")
})
