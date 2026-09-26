# ===========================================================================
# GWR regressions from the second adversarial review.
# ===========================================================================

skip_if_not_installed("GWmodel")
skip_if_not_installed("sp")

# Every warning raised while evaluating `expr`, muffled, beside its value (or
# the error message, as a character string, when it failed).
.r2_catch <- function(expr) {
  w <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w)
}

# Four clusters of 50; `soil` is a regional covariate all but constant inside
# each (0.5 or 1.5 plus noise of SD 0.01), so any window inside one cluster is
# near-collinear with the intercept.  With `exact = TRUE` it is exactly
# constant, a 0/1 indicator.
.r2_clusters <- function(exact = FALSE, seed = 5) {
  set.seed(seed)
  cl <- rep(1:4, each = 50)
  lev <- if (exact) c(0, 1, 0, 1) else c(0.5, 1.5, 0.5, 1.5)
  d <- sf::st_as_sf(data.frame(
    x = c(runif(50, 0, 100), runif(50, 400, 500), runif(50, 0, 100), runif(50, 400, 500)),
    y = c(runif(50, 0, 100), runif(50, 0, 100), runif(50, 400, 500), runif(50, 400, 500)),
    a = rnorm(200),
    soil = lev[cl] + if (exact) 0 else rnorm(200, 0, 0.01)),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + 3 * d$soil + rnorm(200, 0, 0.2)
  d
}


# ---------------------------------------------------------------------------
# A single predictor is collinear with the intercept inside a window too, and
# the survey used to skip it (fewer than two numeric predictors).
# ---------------------------------------------------------------------------

test_that("a one-predictor GWR is surveyed for local collinearity", {
  d <- .r2_clusters()
  r <- .r2_catch(fit_gwr_model(d, "z", "soil", adaptive = TRUE, bandwidth = 20))
  expect_s3_class(r$value, "gwr_fit")
  # The local slopes of soil run to hundreds around a true 3; that has to be
  # said, as it is when a second predictor is present.
  expect_gt(max(abs(coef(r$value)$soil)), 30)
  expect_true(any(grepl("100% of 200 locations have a collinear local design",
                        r$warnings)))
  expect_equal(r$value$info$n_local_collinear, 200L)
  expect_s3_class(r$value$info$local_collinearity, "data.frame")
  expect_equal(r$value$info$condition_index,
               .condition_index(cbind(1, d$soil)))
})


# ---------------------------------------------------------------------------
# The survey's kernel weights are GWmodel's, edge and zero width included, and
# the messages say what a singular window and a non-finite coefficient are.
# ---------------------------------------------------------------------------

test_that("the survey's kernel weights equal GWmodel::gw.weight()", {
  d <- c(0, 0, 50, 100, 100, 150)
  for (k in c("bisquare", "gaussian", "tricube", "boxcar", "exponential")) {
    # Fixed, with points exactly at the kernel's edge.
    expect_identical(.gw_kernel_weights(d, 100, k, FALSE),
                     GWmodel::gw.weight(d, 100, k, FALSE), info = k)
    # Adaptive, the 4th neighbour at 100 and a tie there.
    expect_identical(.gw_kernel_weights(d, 4, k, TRUE),
                     GWmodel::gw.weight(d, 4, k, TRUE), info = k)
    # Zero width: the 2 nearest share one location.  NaN except for the boxcar.
    expect_identical(.gw_kernel_weights(d, 2, k, TRUE),
                     GWmodel::gw.weight(d, 2, k, TRUE), info = k)
    # More neighbours than points: GWmodel widens the kernel.
    expect_equal(.gw_kernel_weights(d, 9, k, TRUE),
                 GWmodel::gw.weight(d, 9, k, TRUE), info = k)
  }
  expect_true(all(is.nan(.gw_kernel_weights(d, 2, "bisquare", TRUE)[1:2])))
})

test_that("a fixed boxcar bandwidth equal to the grid spacing is not called singular", {
  set.seed(3)
  g <- expand.grid(x = seq(0, 1500, by = 100), y = seq(0, 1500, by = 100))
  g$a <- rnorm(nrow(g)); g$b <- rnorm(nrow(g))
  g$z <- 1 + 2 * g$a - g$b + rnorm(nrow(g), 0, 0.3)
  d <- sf::st_as_sf(g, coords = c("x", "y"), crs = 32632)
  r <- .r2_catch(fit_gwr_model(d, "z", c("a", "b"), adaptive = FALSE,
                               bandwidth = 100, kernel = "boxcar"))
  expect_s3_class(r$value, "gwr_fit")
  expect_false(any(grepl("collinear local design", r$warnings)))
  # The survey's windows are GWmodel's: 3 points at a corner, 5 inside.
  W <- GWmodel::gw.weight(as.matrix(stats::dist(g[, c("x", "y")])), 100,
                          "boxcar", FALSE)
  expect_equal(r$value$info$local_collinearity$n_window,
               as.integer(colSums(W > 0)))
})

test_that("an exactly singular window stops the fit with the cause named", {
  d <- .r2_clusters(exact = TRUE)
  r <- .r2_catch(fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE,
                               bandwidth = 20))
  expect_type(r$value, "character")
  expect_match(r$value, "^fit_gwr_model\\(\\): GWR fit failed: inv\\(\\)")
  expect_match(r$value, "local window's design is singular")
  expect_match(r$value, "survey found [0-9]+ singular window")
  # The collinearity warning before it does not promise NaN coefficients.
  expect_true(any(grepl("an exactly singular window makes GWmodel stop the fit",
                        r$warnings)))
})

test_that("non-finite coefficients at co-located points are blamed on the zero-width kernel", {
  set.seed(2)
  sx <- runif(40, 0, 1000); sy <- runif(40, 0, 1000)
  dd <- data.frame(x = rep(sx, each = 4), y = rep(sy, each = 4))
  dd$a <- rnorm(160); dd$b <- rnorm(160)
  dd$z <- 1 + 2 * dd$a - dd$b + rnorm(160, 0, 0.3)
  d <- sf::st_as_sf(dd, coords = c("x", "y"), crs = 32632)
  r <- .r2_catch(fit_gwr_model(d, "z", c("a", "b"), adaptive = TRUE,
                               bandwidth = 4))
  expect_s3_class(r$value, "gwr_fit")
  expect_equal(r$value$info$n_local_singular, 160L)
  # The survey sees the same windows GWmodel does: undefined, not fine.
  expect_equal(r$value$info$n_local_collinear, 160L)
  nf <- grep("returned non-finite coefficients", r$warnings, value = TRUE)
  expect_length(nf, 1L)
  expect_match(nf, "at 160 of them 4 or more observations share one location")
  expect_no_match(nf, "windows are singular")
  # A boxcar kernel keeps the co-located points (weight 1), and fits.
  rb <- .r2_catch(fit_gwr_model(d, "z", c("a", "b"), adaptive = TRUE,
                                bandwidth = 4, kernel = "boxcar"))
  expect_equal(rb$value$info$n_local_singular, 0L)
  expect_false(any(grepl("non-finite coefficients", rb$warnings)))
})
