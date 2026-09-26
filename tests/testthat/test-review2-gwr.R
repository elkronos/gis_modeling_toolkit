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
  # Four neighbours for three parameters: the Gaussian's floor, and every
  # site's four observations are its four nearest.
  r <- .r2_catch(fit_gwr_model(d, "z", c("a", "b"), adaptive = TRUE,
                               bandwidth = 4, kernel = "gaussian"))
  expect_s3_class(r$value, "gwr_fit")
  expect_equal(r$value$info$bandwidth, 4)
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


# ---------------------------------------------------------------------------
# GWmodel's AICc is defined only for tr(S) < n - 2.  Past it the penalty
# changes sign and a (near-)interpolating fit scores a huge negative AICc,
# which compare_models() and gwr_model_selection() ranked first.  The adaptive
# floor used to land bisquare and tricube fits exactly there.
# ---------------------------------------------------------------------------

.r2_pts <- function(n, seed, vars = c("a", "b", "c")) {
  set.seed(seed)
  df <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  for (v in vars) df[[v]] <- rnorm(n)
  sf::st_as_sf(df, coords = c("x", "y"), crs = 32632)
}

test_that("a too-small adaptive bandwidth is raised past the interpolating floor, with a warning", {
  d <- .r2_pts(100, 21)
  d$v <- 1 + d$a + 0.5 * d$b + rnorm(100)
  # 4 parameters; bisquare windows of k neighbours fit on k - 1 points, so 5
  # (the old floor) interpolated every window: R2 = 1, AICc = -15033.
  r <- .r2_catch(fit_gwr_model(d, "v", c("a", "b", "c"), bandwidth = 2))
  expect_s3_class(r$value, "gwr_fit")
  expect_true(any(grepl(paste0("adaptive bandwidth of 2 neighbours is too ",
                               "small for 4 parameters with the bisquare ",
                               "kernel.*using 6"), r$warnings)))
  expect_equal(r$value$info$bandwidth, 6)
  expect_true(is.finite(r$value$info$AICc))
  expect_gt(r$value$info$AICc, 0)
  # The Gaussian weights every point, so parameters + 1 is enough there.
  rg <- .r2_catch(fit_gwr_model(d, "v", c("a", "b", "c"), bandwidth = 2,
                                kernel = "gaussian"))
  expect_equal(rg$value$info$bandwidth, 5)
  expect_true(any(grepl("using 5\\.$", rg$warnings)))
})

test_that("fit_gwr_model() reports an undefined AICc as NA, and compare_models() does not rank it", {
  d <- .r2_pts(25, 1, c("a", "b"))
  d$v <- 1 + d$a + rnorm(25)
  # Five neighbours clear the floor for 3 parameters, but tr(S) = 23.7 is past
  # n - 2 = 23, where GWmodel's AICc is -1687.
  r <- .r2_catch(fit_gwr_model(d, "v", c("a", "b"), bandwidth = 5))
  expect_s3_class(r$value, "gwr_fit")
  expect_lt(r$value$engine$GW.diagnostic$AICc, 0)
  expect_true(is.na(r$value$info$AICc))
  w <- grep("AICc is undefined", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "tr\\(S\\) = 23\\.[0-9]+, is not below n - 2 = 23")
  # The recovered trace is GWmodel's own.
  g <- r$value$engine$GW.diagnostic
  expect_equal(.gwr_trace_s(g$AIC, g$RSS.gw, 25),
               g$AIC - (25 * log(g$RSS.gw / 25) + 25 * log(2 * pi) + 25))
  wide <- suppressWarnings(fit_gwr_model(d, "v", c("a", "b"), bandwidth = 20))
  expect_true(is.finite(wide$info$AICc))
  cmp <- suppressWarnings(compare_models(list(small = r$value, wide = wide)))
  expect_true(is.na(cmp$AICc[cmp$model == "small"]))
  expect_equal(cmp$AICc[cmp$model == "wide"], wide$info$AICc)
})

test_that("gwr_model_selection() ranks a model with an undefined AICc last", {
  d <- .r2_pts(30, 3, c("a", "n1", "n2"))
  d$v <- 1 + d$a + rnorm(30, 0, 0.5)
  # At 6 neighbours the three-variable model has tr(S) = 28.5 > n - 2 = 28 and
  # a GWmodel AICc of -3871, which made it the selected model.
  r <- .r2_catch(gwr_model_selection(d, "v", c("a", "n1", "n2"), bandwidth = 6))
  sel <- r$value
  expect_s3_class(sel, "gwr_model_selection")
  expect_lt(min(sel$raw[[2]][, 3]), 0)
  expect_identical(sel$best, "a")
  expect_true(any(grepl("AICc is undefined for 1 of 6 model\\(s\\) at bandwidth 6",
                        r$warnings)))
  expect_true(is.na(sel$table$criterion[6]))
  expect_identical(sel$table$n_vars[6], 3L)
  expect_true(all(is.finite(sel$table$criterion[1:5])))
})

test_that("gwr_model_selection() warns when it raises a too-small adaptive bandwidth", {
  d <- .r2_pts(60, 5, c("a", paste0("n", 1:4)))
  d$z <- 2 * d$a + rnorm(60, 0, 0.5)
  # Six parameters in the full model: bisquare needs 8 neighbours.
  r <- .r2_catch(gwr_model_selection(d, "z", c("a", paste0("n", 1:4)),
                                     bandwidth = 4))
  expect_s3_class(r$value, "gwr_model_selection")
  expect_identical(r$value$bandwidth, 8L)
  expect_true(any(grepl(paste0("adaptive bandwidth of 4 neighbours is too ",
                               "small for the full 5-predictor model with the ",
                               "bisquare kernel; using 8"), r$warnings)))
})
