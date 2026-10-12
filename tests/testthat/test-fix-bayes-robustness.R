# Two early exits in fit_bayesian_spatial_model().
#
# The response was checked for its type and for n >= 2 only, so a constant
# response, or a handful of rows against several hundred basis functions,
# reached the compiler.

.fbx_pts <- function(n = 120, seed = 6) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = runif(n, 0, 1e5), y = runif(n, 0, 1e5),
               z = stats::rnorm(n), k = 7, a = stats::rnorm(n),
               cnt = 4L, flag = TRUE),
    coords = c("x", "y"), crs = 32617)
}

test_that("a constant response is refused before anything is compiled", {
  skip_if_not_installed("brms")
  d <- .fbx_pts()
  calls <- 0L
  local_mocked_bindings(
    brm = function(...) { calls <<- calls + 1L; structure(list(), class = "fbx_stub") },
    .package = "brms")
  expect_error(
    suppressMessages(fit_bayesian_spatial_model(d, "k", "a", compute_loo = FALSE,
                                                check_convergence = FALSE)),
    "response 'k' takes a single value \\(7\\) in all 120 row\\(s\\)")
  # Whatever the family: a count or a logical with one value has no more to
  # explain than a Gaussian one.
  expect_error(
    suppressMessages(fit_bayesian_spatial_model(d, "cnt", "a",
                                                family = stats::poisson(),
                                                compute_loo = FALSE,
                                                check_convergence = FALSE)),
    "response 'cnt' takes a single value \\(4\\)")
  expect_error(
    suppressMessages(fit_bayesian_spatial_model(d, "flag", "a",
                                                family = brms::bernoulli(),
                                                compute_loo = FALSE,
                                                check_convergence = FALSE)),
    "response 'flag' takes a single value \\(TRUE\\)")
  # It reached brms::brm() every time.
  expect_identical(calls, 0L)
})

test_that("fewer rows than a tenth of the basis functions is a classed warning", {
  skip_if_not_installed("brms")
  d <- .fbx_pts()
  local_mocked_bindings(brm = function(...) structure(list(), class = "fbx_stub"),
                        .package = "brms")
  local_mocked_bindings(estimate_sac_range = function(...) 9e4,
                        .package = "spatialkit")
  run <- function(dat, ...) {
    w <- list()
    fit <- withCallingHandlers(
      suppressMessages(fit_bayesian_spatial_model(dat, "z", "a", compute_loo = FALSE,
                                                  check_convergence = FALSE, ...)),
      warning = function(c) { w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning") })
    list(fit = fit, w = Filter(function(c) inherits(c, "spatialkit_gp_basis_size"), w))
  }
  # 30 rows against the derived basis (several hundred functions).
  small <- run(d[1:30, ])
  nb <- small$fit$info$gp_n_basis
  expect_gt(nb, 300)
  expect_length(small$w, 1L)
  expect_match(conditionMessage(small$w[[1]]),
               sprintf("30 row(s) against %d basis functions (gp_k = %d per dimension)",
                       nb, small$fit$info$gp_k), fixed = TRUE)
  # 120 rows against 36 functions, and against the derived basis: no warning
  # at more than a tenth.
  expect_length(run(d, gp_k = 6)$w, 0L)
  full <- run(d)
  expect_lt(full$fit$info$gp_n_basis / 10, 120)
  expect_length(full$w, 0L)
})
