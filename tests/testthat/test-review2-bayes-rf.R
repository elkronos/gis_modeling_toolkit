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


# ---------------------------------------------------------------------------
# gaps G1.4: the length-scale prior keeps its dpar
# ---------------------------------------------------------------------------

test_that("a categorical fit's automatic lscale prior validates in brms", {
  skip_if_not_installed("brms")
  d <- .r2b_pts()
  r <- .r2b_quiet(.r2b_capture_fit(d, "cat3", "a", family = brms::categorical(),
                                   gp_k = 5))
  pr <- r$args$prior
  ls <- as.data.frame(pr)[pr$class == "lscale", ]
  # One row per coefficient PER category, each addressed by its dpar.  The
  # prior used to carry the coef names only, so every row landed on dpar ""
  # twice and brms refused the model before sampling.
  expect_setequal(unique(ls$dpar), c("mumid", "muhi"))
  expect_false(any(duplicated(ls[, c("coef", "dpar")])))
  expect_no_error(suppressWarnings(brms::validate_prior(
    pr, formula = r$args$formula, data = r$args$data, family = r$args$family)))
})

test_that("a mixture fit's automatic lscale prior validates in brms", {
  skip_if_not_installed("brms")
  d <- .r2b_pts()
  r <- .r2b_quiet(.r2b_capture_fit(
    d, "z", "a", family = brms::mixture(stats::gaussian(), stats::gaussian()),
    gp_k = 5))
  pr <- r$args$prior
  expect_setequal(unique(pr$dpar[pr$class == "lscale"]), c("mu1", "mu2"))
  expect_no_error(suppressWarnings(brms::validate_prior(
    pr, formula = r$args$formula, data = r$args$data, family = r$args$family)))
})

test_that("a user's dpar-level lscale prior is expanded onto its own dpar only", {
  skip_if_not_installed("brms")
  d <- .r2b_pts()
  up <- brms::set_prior("normal(0, 1)", class = "lscale", dpar = "mumid") +
    brms::set_prior("normal(0, 2)", class = "lscale", dpar = "muhi",
                    coef = "gp..x..y..y")
  r <- .r2b_quiet(.r2b_capture_fit(d, "cat3", "a", family = brms::categorical(),
                                   gp_k = 5, prior = up))
  pr <- as.data.frame(r$args$prior)
  ls <- pr[pr$class == "lscale", ]
  # mumid's two coefficients get the dpar-level prior; muhi keeps the one
  # coefficient-level row it was given and nothing else.
  expect_setequal(ls$coef[ls$dpar == "mumid"], c("gp..x..y..x", "gp..x..y..y"))
  expect_true(all(ls$prior[ls$dpar == "mumid"] == "normal(0, 1)"))
  expect_identical(ls$prior[ls$dpar == "muhi"], "normal(0, 2)")
  expect_false(any(ls$dpar == ""))
  expect_no_error(suppressWarnings(brms::validate_prior(
    r$args$prior, formula = r$args$formula, data = r$args$data,
    family = r$args$family)))

  # A global prior alongside a coefficient-level one for the same parameter is
  # not expanded on top of it (that duplicate was refused as well).
  up2 <- brms::set_prior("normal(0, 1)", class = "lscale") +
    brms::set_prior("normal(0, 3)", class = "lscale", coef = "gp..x..y..x")
  r2 <- .r2b_quiet(.r2b_capture_fit(d, "z", "a", gp_k = 5, prior = up2))
  ls2 <- as.data.frame(r2$args$prior)
  ls2 <- ls2[ls2$class == "lscale", ]
  expect_identical(ls2$prior[ls2$coef == "gp..x..y..x"], "normal(0, 3)")
  expect_identical(ls2$prior[ls2$coef == "gp..x..y..y"], "normal(0, 1)")
  expect_no_error(suppressWarnings(brms::validate_prior(
    r2$args$prior, formula = r2$args$formula, data = r2$args$data,
    family = r2$args$family)))
})


# ---------------------------------------------------------------------------
# gaps G1.5: a factor response is refused unless the family is categorical
# ---------------------------------------------------------------------------

test_that("a factor response under bernoulli is refused before anything is fitted", {
  skip_if_not_installed("brms")
  d <- .r2b_pts()
  called <- FALSE
  local_mocked_bindings(brm = function(...) { called <<- TRUE; stop("reached brm") },
                        .package = "brms")
  expect_error(
    .r2b_quiet(fit_bayesian_spatial_model(d, "pres", "a",
                                          family = brms::bernoulli())),
    "response 'pres' is factor, not numeric, and the family is bernoulli")
  expect_false(called)
  expect_error(
    .r2b_quiet(fit_bayesian_spatial_model(d, "pres", "a", family = poisson())),
    "the family is poisson")
  expect_false(called)
})

test_that("the categorical, ordinal and numeric-binary cases still reach brms", {
  skip_if_not_installed("brms")
  d <- .r2b_pts()
  d$lgl <- d$bin == 1L
  ok <- function(resp, family) {
    r <- .r2b_quiet(.r2b_capture_fit(d, resp, "a", family = family, gp_k = 5))
    expect_false(is.null(r$args), info = resp)
  }
  ok("cat3", brms::categorical())
  ok("ord3", brms::cumulative())
  ok("bin",  brms::bernoulli())
  ok("lgl",  brms::bernoulli())
})


# ---------------------------------------------------------------------------
# models MED-3 (Bayesian): a per-category epred is an error, not "draw failed"
# ---------------------------------------------------------------------------

.r2b_ordinal_fit <- function(d) {
  new_spatial_fit(
    "bayesian_fit",
    engine = structure(list(family = list(family = "cumulative")),
                       class = "brmsfit"),
    formula = ord3 ~ a, response_var = "ord3", predictor_vars = "a",
    data_sf = d,
    info = list(coord_scaling = list(x_center = 5e5, x_scale = 300,
                                     y_center = 5e6, y_scale = 300)))
}

test_that("predict() and fitted() on an ordinal fit say why they cannot answer", {
  skip_if_not_installed("brms")
  d  <- .r2b_pts(n = 40)
  nd <- d[1:5, ]
  fit <- .r2b_ordinal_fit(d)
  local_mocked_bindings(
    posterior_epred = function(object, newdata, ...)
      array(0.3, dim = c(20L, nrow(newdata), 3L)),
    posterior_predict = function(object, newdata, ...)
      matrix(2, nrow = 20L, ncol = nrow(newdata)),
    .package = "brms")

  # It used to log "posterior draw failed" and return five NAs.
  expect_error(.r2b_quiet(predict(fit, newdata = nd)),
               "cumulative.*probability per response category.*type = \"predict\"")
  expect_error(.r2b_quiet(predict(fit, newdata = nd, draws = TRUE)),
               "probability per response category")
  expect_error(fitted(fit), "probability per response category")
  expect_error(residuals(fit), "probability per response category")
  # type = "predict" is unaffected.
  p <- .r2b_quiet(predict(fit, newdata = nd, type = "predict"))
  expect_equal(p, rep(2, 5))
})
