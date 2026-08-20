# tests/testthat/helper-lmfit.R
# ---------------------------------------------------------------------------
# A minimal deterministic spatial_fit backed by lm().
#
# Lets prediction machinery -- predict_surface(), chunking, covariate joins --
# be tested end to end without GWmodel or a Stan toolchain, which would
# otherwise make these paths untestable in the default suite.
# ---------------------------------------------------------------------------

lm_spatial_fit <- function(data_sf, response_var = "z",
                           predictor_vars = character(0)) {
  d <- sf::st_drop_geometry(data_sf)
  rhs <- if (length(predictor_vars)) paste(predictor_vars, collapse = " + ") else "1"
  fml <- stats::as.formula(paste(response_var, "~", rhs))
  new_spatial_fit(
    subclass       = "lmsurf_fit",
    engine         = stats::lm(fml, data = d),
    formula        = fml,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    data_sf        = data_sf
  )
}

predict.lmsurf_fit <- function(object, newdata = NULL, draws = FALSE, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  nd  <- sf::st_drop_geometry(newdata)
  mu  <- as.numeric(stats::predict(object$engine, newdata = nd))
  if (!isTRUE(draws)) return(mu)
  # A stand-in draw matrix (n_draws x n_newdata), matching the shape
  # predict.bayesian_fit() returns, so the se path can be exercised.
  s <- stats::sigma(object$engine)
  matrix(rep(mu, each = 50) + stats::rnorm(50 * length(mu), 0, s),
         nrow = 50, ncol = length(mu))
}

# plot.spatial_fit() goes through the residuals()/fitted() generics, so the
# stand-in backend needs them too.
residuals.lmsurf_fit <- function(object, ...) stats::residuals(object$engine)
fitted.lmsurf_fit    <- function(object, ...) stats::fitted(object$engine)

registerS3method("predict",   "lmsurf_fit", predict.lmsurf_fit)
registerS3method("residuals", "lmsurf_fit", residuals.lmsurf_fit)
registerS3method("fitted",    "lmsurf_fit", fitted.lmsurf_fit)

surf_test_points <- function(n = 120, extent = 1000, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, extent); y <- runif(n, 0, extent); w <- rnorm(n)
  sf::st_as_sf(
    data.frame(x = x, y = y, w = w,
               z = 0.01 * x + 0.02 * y + 0.5 * w + rnorm(n, 0, 0.1)),
    coords = c("x", "y"), crs = 3857
  )
}
