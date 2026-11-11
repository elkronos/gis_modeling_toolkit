# At its default control the Bayesian fit runs without divergent
# transitions on ordinary data.
#
# At adapt_delta = 0.9, the previous default, four default-length fits (four
# chains of 2,000 iterations) of 250 points over a squared-exponential field
# with a 60 km effective range had 2, 4, 4 and 12 divergent transitions, and
# one divergent transition sets convergence_ok to FALSE.  At 0.99 they had
# none.  The fixtures are the ones those counts were measured on, with rstan
# 2.32; the sampler is pinned to rstan, since another backend or Stan
# version draws other chains from the same seed.
#
# Opt-in, via SPATIALKIT_TEST_BRMS (see test-bayes-smoke.R): about 8
# minutes on one core, four compiles and four fits.

.bdc_fixture <- function(id) {
  F <- sqrt(2 * log(20)); n <- 250; ell <- 6e4 / F
  if (identical(id, "A")) {
    set.seed(1); xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5)); a <- rnorm(n)
    e <- eigen(exp(-as.matrix(dist(xy))^2 / (2 * ell^2)), symmetric = TRUE)
    set.seed(5)
    f <- as.numeric(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(n)))
    seed <- 123L
  } else {
    s <- as.integer(sub("^B", "", id))
    set.seed(100 + s); xy <- cbind(runif(n, 0, 1e5), runif(n, 0, 1e5)); a <- rnorm(n)
    e <- eigen(exp(-as.matrix(dist(xy))^2 / (2 * ell^2)), symmetric = TRUE)
    f <- as.numeric(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(n)))
    seed <- s
  }
  d <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], a = a,
                               z = 1 + 0.8 * a + f + rnorm(n, 0, sqrt(0.3))),
                    coords = c("x", "y"), crs = 32617)
  list(data = d, seed = seed)
}

test_that("default-length fits at the default control have no divergent transitions on four ordinary layers", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan tests")
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  for (id in c("A", "B1", "B2", "B3")) {
    fx <- .bdc_fixture(id)
    fit <- NULL
    utils::capture.output(suppressWarnings(suppressMessages(
      fit <- fit_bayesian_spatial_model(fx$data, "z", "a", cores = 1,
                                        seed = fx$seed, compute_loo = FALSE,
                                        backend = "rstan"))),
      type = "output")
    ctl <- tryCatch(fit$engine$fit@stan_args[[1]]$control, error = function(e) NULL)
    if (!is.null(ctl$adapt_delta)) expect_equal(ctl$adapt_delta, 0.99, info = id)
    expect_equal(fit$info$convergence_diagnostics$n_divergent, 0, info = id)
    expect_true(isTRUE(fit$info$convergence_ok), info = id)
  }
})
