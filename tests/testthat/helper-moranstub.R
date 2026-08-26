# tests/testthat/helper-moranstub.R
# ---------------------------------------------------------------------------
# A stand-in spatial_fit for residual_morans_i(), which needs only
# residuals(fit) and fit$data_sf.
#
# The residuals method is registered against a TEST-ONLY subclass.  Doing it
# against "spatial_fit" itself -- the package's own base class -- puts a
# method the package does not ship into the global S3 dispatch table for the
# rest of the session.  testthat runs files alphabetically, so every file
# after test-morans-*.R then saw a residuals.spatial_fit that does not exist
# in production: it would mask a genuinely missing method, and a single-file
# run would behave differently from a suite run.
#
# Registering here rather than in a test body also means it happens once, at
# helper load, instead of on every test execution.
# ---------------------------------------------------------------------------

moran_stub_fit <- function(points_sf, resid) {
  structure(
    list(data_sf = points_sf, residuals = resid, engine = list()),
    class = c("moran_stub_fit", "spatial_fit")
  )
}

residuals.moran_stub_fit <- function(object, ...) object$residuals

registerS3method("residuals", "moran_stub_fit", residuals.moran_stub_fit)

# A row-standardised k-nearest-style weight matrix built from a fixed seed.
moran_stub_weights <- function(n, k, seed, value = 1 / k) {
  set.seed(seed)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nbrs <- sample(setdiff(seq_len(n), i), k)
    W[i, nbrs] <- value
  }
  W
}
