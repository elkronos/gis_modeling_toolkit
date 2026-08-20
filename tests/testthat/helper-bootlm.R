# tests/testthat/helper-bootlm.R
# ---------------------------------------------------------------------------
# A deliberately stochastic learner for the CV reproducibility tests.
#
# It MUST consume RNG inside every fold.  A deterministic fit_fn passes the
# parallel/sequential equivalence test whether or not fold seeding works,
# which is exactly how the original gap survived to release: the previous
# test suite only checked that a `parallel` argument existed in formals().
# ---------------------------------------------------------------------------

boot_fit_fn <- function(train_sf) {
  d   <- sf::st_drop_geometry(train_sf)
  idx <- sample.int(nrow(d), nrow(d), replace = TRUE)   # <- consumes RNG
  new_spatial_fit(
    subclass       = "bootlm_fit",
    engine         = stats::lm(z ~ w, data = d[idx, , drop = FALSE]),
    formula        = z ~ w,
    response_var   = "z",
    predictor_vars = "w",
    data_sf        = train_sf
  )
}

predict.bootlm_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  as.numeric(stats::predict(object$engine,
                            newdata = sf::st_drop_geometry(newdata)))
}

# Required.  .cv_fit_one_fold() calls predict() from inside the spatialkit
# namespace, where a method defined only in the test environment is invisible
# to S3 dispatch.  registerS3method() puts it in the dispatch table, which is
# consulted regardless of where the generic was called from.
registerS3method("predict", "bootlm_fit", predict.bootlm_fit)

make_cv_test_points <- function(n = 120, seed = 1) {
  set.seed(seed)
  sf::st_as_sf(
    data.frame(x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000),
               z = stats::rnorm(n), w = stats::rnorm(n)),
    coords = c("x", "y"), crs = 3857
  )
}
