# ===========================================================================
# GWR regressions from the third review (slice S7: GWR and the fit classes).
# ===========================================================================

# Every warning raised while evaluating `expr`, muffled, beside its value (or
# the error message, as a character string, when it failed).
.r3_catch <- function(expr) {
  w <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w)
}

# A smooth temperature field, in degrees C and in kelvin: the same predictor
# with a different origin, so every local slope is the same in both.
.r3_temperature <- function() {
  set.seed(1); n <- 200
  x <- runif(n, 0, 10000); y <- runif(n, 0, 10000)
  tC <- 15 + 4 * sin(x / 3000) + 3 * cos(y / 2500) + rnorm(n, 0, 1)
  d <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, tC = tC,
                               tK = tC + 273.15, nz = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$yield <- 50 - 0.8 * d$tC + rnorm(n, 0, 1)
  d
}


# ---------------------------------------------------------------------------
# S7-GWR-CLASSES-1: the collinearity verdict on the slopes must not depend on
# a predictor's origin (degrees C vs kelvin), while a regional covariate that
# is nearly constant inside a window is still caught.
# ---------------------------------------------------------------------------

test_that("the slope index is origin- and unit-free and matches its definition", {
  set.seed(11); n <- 150
  xy <- cbind(runif(n, 0, 1000), runif(n, 0, 1000))
  xm <- cbind(a = rnorm(n, 5, 2), b = rnorm(n))
  s0 <- .gwr_local_collinearity(xy, xm, TRUE, 30, "bisquare")
  expect_named(s0, c("row", "x", "y", "n_window", "cn", "cn_slopes"))
  expect_true(all(is.finite(s0$cn_slopes) & s0$cn_slopes >= 1))
  # A shift of origin and a change of units move `cn` but not `cn_slopes`.
  xs <- cbind(a = xm[, "a"] * 1000 + 1e5, b = xm[, "b"] - 40)
  s1 <- .gwr_local_collinearity(xy, xs, TRUE, 30, "bisquare")
  expect_equal(s1$cn_slopes, s0$cn_slopes, tolerance = 1e-8)
  expect_gt(median(s1$cn), 10 * median(s0$cn))
  # By hand at one location: locally centred, globally scaled, weighted.
  i <- 17L
  d <- sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
  w <- .gw_kernel_weights(d, 30, "bisquare", TRUE); k <- which(w > 1e-8)
  ww <- w[k] / sum(w[k])
  z <- sqrt(ww) * sweep(sweep(xm[k, ], 2, colSums(ww * xm[k, ])), 2,
                        apply(xm, 2, sd), "/")
  sv <- svd(z)$d
  expect_equal(s0$cn_slopes[i], max(1, max(sv)) / min(1, min(sv)))
  # A predictor exactly constant inside the window is singular outright.
  expect_identical(.gwr_slope_index(c(1, 1, 1, 1), cbind(c(3, 3, 3, 3), 1:4),
                                    c(1, 1)), Inf)
})

test_that("a regional covariate is still flagged, including a cluster at the global mean", {
  # Three clusters at 0.5 / 1.0 / 1.5 (plus noise of SD 0.01).  Centring at
  # the GLOBAL mean would score the middle cluster near 1, although its local
  # slope is undetermined; the slope index centres in the window.
  set.seed(6); cl <- rep(1:3, each = 60)
  xy <- cbind(c(0, 400, 800)[cl] + runif(180, 0, 100), runif(180, 0, 100))
  soil <- c(0.5, 1.0, 1.5)[cl] + rnorm(180, 0, 0.01)
  s <- .gwr_local_collinearity(xy, cbind(soil = soil), TRUE, 20, "bisquare")
  flagged <- tapply(.gwr_slopes_collinear(s), cl, sum)
  expect_true(all(flagged >= 55))
  expect_gte(flagged[[2L]], 55)
})

test_that("a far-from-zero origin flags the slopes only where GWmodel's solve degrades", {
  set.seed(1); n <- 200
  x <- runif(n, 0, 10000); y <- runif(n, 0, 10000)
  el <- 500 + 150 * sin(x / 3000) + 100 * cos(y / 2500) + rnorm(n, 0, 5)
  # Elevation in metres: the old uncentred rule flagged most windows.
  s <- .gwr_local_collinearity(cbind(x, y), cbind(el), TRUE, 30, "bisquare")
  expect_gt(sum(!is.finite(s$cn) | s$cn > 30), 100)
  expect_equal(sum(.gwr_slopes_collinear(s)), 0L)
  # Shifted until the uncentred index passes 1e6: flagged, on `cn` alone.
  s7 <- .gwr_local_collinearity(cbind(x, y), cbind(el + 1e7), TRUE, 30, "bisquare")
  expect_equal(s7$cn_slopes, s$cn_slopes, tolerance = 1e-6)
  expect_gt(sum(.gwr_slopes_collinear(s7)), 0L)
  expect_true(all(s7$cn[.gwr_slopes_collinear(s7)] > 1e6))
  # A survey made before cn_slopes existed falls back to cn > 30.
  old <- s[, c("row", "x", "y", "n_window", "cn")]
  expect_identical(.gwr_slopes_collinear(old), !is.finite(old$cn) | old$cn > 30)
})

test_that("the same field in degrees C and in kelvin gets the same verdict and a drawable map", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .r3_temperature()
  rC <- .r3_catch(fit_gwr_model(d, "yield", "tC", bandwidth = 30))
  rK <- .r3_catch(fit_gwr_model(d, "yield", "tK", bandwidth = 30))
  expect_s3_class(rC$value, "gwr_fit"); expect_s3_class(rK$value, "gwr_fit")
  # The slopes are the same surface ...
  expect_lt(max(abs(coef(rC$value)$tC - coef(rK$value)$tK)), 1e-6)
  # ... so neither fit warns, globally or locally.
  expect_length(rC$warnings, 0L)
  expect_length(rK$warnings, 0L)
  expect_identical(rC$value$info$n_local_collinear, rK$value$info$n_local_collinear)
  expect_identical(rK$value$info$n_local_collinear, 0L)
  expect_equal(rK$value$info$condition_index, 1)
  expect_equal(rK$value$info$local_collinearity$cn_slopes,
               rC$value$info$local_collinearity$cn_slopes, tolerance = 1e-8)
  # The intercept-inclusive index is kept, and it does see the kelvin origin.
  expect_true(all(rK$value$info$local_collinearity$cn > 30))
  # The kelvin slope map is drawn, with nothing masked (it used to be refused).
  skip_if_not_installed("ggplot2")
  pK <- plot(rK$value, type = "coefficients", term = "tK")
  expect_match(pK$labels$subtitle, "No location masked")
  # The Intercept map is masked by the intercept-inclusive index.
  lcC <- rC$value$info$local_collinearity
  n_int <- sum(!is.finite(lcC$cn) | lcC$cn > 30)
  expect_gt(n_int, 0L)
  pI <- plot(rC$value, type = "coefficients", term = "Intercept")
  expect_match(pI$labels$subtitle,
               sprintf("^%d of 200 locations masked.*%d with a collinear local design \\(condition index with the intercept > 30\\)",
                       n_int, n_int))
  pC <- plot(rC$value, type = "coefficients", term = "tC")
  expect_match(pC$labels$subtitle, "No location masked")
  # Every intercept masked: refused, and the message says what to do.
  expect_error(plot(rK$value, type = "coefficients", term = "Intercept"),
               "centre the predictors to map it")
})

test_that("a far-from-zero predictor beside a second one no longer draws the global warning", {
  # v2.0.0 already warned on this for two or more predictors.
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  d <- .r3_temperature()
  rK <- .r3_catch(fit_gwr_model(d, "yield", c("tK", "nz"), bandwidth = 30))
  rC <- .r3_catch(fit_gwr_model(d, "yield", c("tC", "nz"), bandwidth = 30))
  expect_false(any(grepl("global predictors|collinear local design", rK$warnings)))
  expect_equal(rK$value$info$condition_index, rC$value$info$condition_index)
  expect_lt(rK$value$info$condition_index, 30)
  expect_identical(rK$value$info$n_local_collinear, rC$value$info$n_local_collinear)
  # Two predictors that really move together still warn globally.
  d$tC2 <- d$tC + rnorm(nrow(d), 0, 0.01)
  r2 <- .r3_catch(fit_gwr_model(d, "yield", c("tC", "tC2"), bandwidth = 60))
  expect_true(any(grepl("global predictors \\(centred\\) have scaled condition index",
                        r2$warnings)))
})


# ---------------------------------------------------------------------------
# S7-GWR-CLASSES-3: the adaptive floor is not enough when neighbours tie at
# the kernel's edge (a regular grid); the warning and the fit error say so.
# ---------------------------------------------------------------------------

# A 10 x 10 grid of spacing 100 with predictor `a`.  Away from the four
# corners a point has three or four neighbours at 100 and none nearer, so at
# the bisquare floor of 4 neighbours, the point itself counted, the kernel's
# edge is at 100: they all get weight 0, and 96 local regressions hold their
# own point alone.  X'WX is then [1, a; a, a^2], of rank 1 for two parameters.
.r3_tied_grid <- function(a) {
  g <- expand.grid(x = 5e5 + (0:9) * 100, y = 5e6 + (0:9) * 100)
  g$a <- a; g$z <- 1 + 2 * g$a + rnorm(100)
  sf::st_as_sf(g, coords = c("x", "y"), crs = 32632)
}

test_that("the floor warning and the singular-fit error name ties at the kernel's edge", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  # Whether GWmodel calls a one-point window singular is a matter of rounding.
  # Eliminating [1, a; a, a^2] leaves a^2 - a * a as the last pivot: 0 where
  # the product is rounded before the subtraction, and the rounding error of
  # a * a where the two are fused into one operation (as on arm64), unless
  # a * a is exact.  The predictor takes multiples of 1/8 in [-1, 1]: their
  # squares are exact and none outranks the 1 as first pivot, so every step is
  # exact, the pivot is 0 and the fit stops on any platform.
  set.seed(2)
  d <- .r3_tied_grid(sample(seq(-1, 1, by = 0.125), 100, replace = TRUE))
  r <- .r3_catch(fit_gwr_model(d, "z", "a", bandwidth = 2))
  expect_true(any(grepl("using 4\\. That is enough unless several neighbours tie at the kernel's edge",
                        r$warnings)))
  expect_type(r$value, "character")
  expect_match(r$value, "survey found 96 singular window")
  expect_match(r$value, "neighbours tied at the kernel's edge \\(a regular grid\\)")
  # Kernels that keep the edge point are not told about ties.
  rg <- .r3_catch(fit_gwr_model(d, "z", "a", bandwidth = 2, kernel = "gaussian"))
  expect_false(any(grepl("tie at the kernel's edge", rg$warnings)))
})

test_that("a continuous predictor on the tied grid is flagged whether or not GWmodel stops", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  # The squares of normal draws are rounded, so GWmodel stops on the 96
  # one-point windows or returns coefficients for them according to the
  # platform's arithmetic.  The survey does not depend on that: it counts the
  # points that carry weight, and warns before the fit.
  set.seed(2)
  d <- .r3_tied_grid(rnorm(100))
  r <- .r3_catch(fit_gwr_model(d, "z", "a", bandwidth = 2))
  expect_true(any(grepl("using 4\\. That is enough unless several neighbours tie at the kernel's edge",
                        r$warnings)))
  expect_true(any(grepl("local collinearity: 96% of 100 locations have a collinear local design",
                        r$warnings)))
  if (is.character(r$value)) {
    expect_match(r$value, "survey found 96 singular window")
    expect_match(r$value, "neighbours tied at the kernel's edge \\(a regular grid\\)")
  } else {
    # The fit that comes back carries the survey: the 96 windows are there to
    # be read, each with one point and an infinite condition index.
    expect_s3_class(r$value, "gwr_fit")
    lc <- r$value$info$local_collinearity
    expect_identical(sum(lc$n_window == 1L & is.infinite(lc$cn)), 96L)
    expect_identical(r$value$info$n_local_collinear, 96L)
    skip(paste0("GWmodel returned coefficients for the 96 one-point windows ",
                "instead of stopping, so the singular-fit error is not reached ",
                "with a continuous predictor on this platform"))
  }
})


# ---------------------------------------------------------------------------
# S7-GWR-CLASSES-4: one range check for an adaptive count in all three GWR
# entry points, before anything is fitted.
# ---------------------------------------------------------------------------

test_that("gwr_model_selection() refuses an adaptive count below 1 or above R's integers", {
  set.seed(1); n <- 60
  dat <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                                 a = rnorm(n), b = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  dat$z <- 1 + 2 * dat$a + rnorm(n)
  never <- function(...) stop("engine reached")
  expect_error(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 3e9, .engine = never),
               "`bandwidth` must be a single number of nearest neighbours when adaptive = TRUE and at most 2147483647")
  expect_error(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 0.5, .engine = never),
               "`bandwidth` must be a single number of nearest neighbours when adaptive = TRUE and at least 1; got 0.5")
  # A fixed distance below 1 is a distance, not a count, and passes.
  expect_error(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 0.5,
                                   adaptive = FALSE, .engine = never),
               "engine reached")
})

test_that("cv_gwr() refuses an adaptive count below 1 up front, not in every fold", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  set.seed(1); n <- 60
  dat <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                                 a = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  dat$z <- 1 + 2 * dat$a + rnorm(n)
  expect_error(cv_gwr(dat, "z", "a", bandwidth = 0.5, k = 3),
               "^cv_gwr\\(\\): `bandwidth` must be a single number of nearest neighbours when adaptive = TRUE and at least 1")
  expect_error(cv_gwr(dat, "z", "a", bandwidth = 3e9, k = 3),
               "^cv_gwr\\(\\): `bandwidth` .* at most 2147483647")
})


# ---------------------------------------------------------------------------
# S7-GWR-CLASSES-5: "use a larger bandwidth" is not advice when the adaptive
# bandwidth is already every observation.
# ---------------------------------------------------------------------------

test_that("the undefined-AICc warning gives advice that can be followed", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  set.seed(4); n <- 4
  dat <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                                 a = rnorm(n), b = rnorm(n)),
                      coords = c("x", "y"), crs = 32632)
  dat$z <- 1 + 2 * dat$a + rnorm(n)
  r <- .r3_catch(fit_gwr_model(dat, "z", c("a", "b")))
  w <- grep("AICc is undefined", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "Even the widest adaptive window \\(all 4 observations\\) leaves too few residual degrees of freedom for 3 parameters")
  expect_no_match(w, "use a larger one")
  rs <- .r3_catch(gwr_model_selection(dat, "z", c("a", "b")))
  ws <- grep("AICc is undefined", rs$warnings, value = TRUE)
  expect_length(ws, 1L)
  expect_match(ws, "Even the widest adaptive window \\(all 4 observations\\)")
  expect_no_match(ws, "Use a larger bandwidth")
  # Below n, a larger bandwidth is still the advice.
  set.seed(1); n <- 25
  d2 <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                                a = rnorm(n), b = rnorm(n)),
                     coords = c("x", "y"), crs = 32632)
  d2$v <- 1 + d2$a + rnorm(n)
  r2 <- .r3_catch(fit_gwr_model(d2, "v", c("a", "b"), bandwidth = 5))
  w2 <- grep("AICc is undefined", r2$warnings, value = TRUE)
  expect_length(w2, 1L)
  expect_match(w2, "The bandwidth \\(5 neighbours\\) is too small for 3 parameters; use a larger one\\.")
})


# ---------------------------------------------------------------------------
# S7-GWR-CLASSES-6: gwr_model_selection() says what fit_gwr_model() says
# about a fixed bandwidth in the wrong units, and about a singular window.
# ---------------------------------------------------------------------------

test_that("gwr_model_selection() warns about a tiny fixed bandwidth and explains a singular window", {
  expect_match(.gwr_singular_hint(), "^ -- at least one local window's design is singular")
  expect_no_match(.gwr_singular_hint(), "survey found")
  expect_match(.gwr_singular_hint(3L), "The collinearity survey found 3 singular window")
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  set.seed(1); n <- 60
  dat <- sf::st_as_sf(data.frame(x = -80 + runif(n, 0, 0.5), y = 35 + runif(n, 0, 0.5),
                                 a = rnorm(n), b = rnorm(n)),
                      coords = c("x", "y"), crs = 4326)
  dat$z <- 1 + 2 * dat$a + rnorm(n)
  r <- .r3_catch(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 0.2,
                                     adaptive = FALSE))
  expect_true(any(grepl("^gwr_model_selection\\(\\): a fixed bandwidth of 0.2 is less than a ten-thousandth",
                        r$warnings)))
  # What GWmodel then does with empty windows is a backend property; when it
  # is the singular-matrix error, the error explains it.
  if (is.character(r$value) && grepl("singular", r$value))
    expect_match(r$value, "local window's design is singular")
  # A bandwidth in the units the sweep runs in draws no such warning.
  ok <- .r3_catch(gwr_model_selection(dat, "z", c("a", "b"), bandwidth = 20000,
                                      adaptive = FALSE))
  expect_s3_class(ok$value, "gwr_model_selection")
  expect_false(any(grepl("ten-thousandth", ok$warnings)))
})


# ---------------------------------------------------------------------------
# S7-GWR-CLASSES-8: print() shows a fixed bandwidth in full, with its unit.
# ---------------------------------------------------------------------------

test_that("print() shows a fixed bandwidth with its unit and an adaptive one as neighbours", {
  pts <- surf_test_points(30)          # EPSG:3857, metres
  gwr <- new_spatial_fit(
    "gwr_fit", engine = list(), formula = z ~ w, response_var = "z",
    predictor_vars = "w", data_sf = pts,
    info = list(bandwidth = 122372.3, adaptive = FALSE, kernel = "bisquare",
                AICc = NA_real_, bandwidth_is_fallback = FALSE))
  txt <- paste(utils::capture.output(print(gwr)), collapse = "\n")
  expect_match(txt, "Bandwidth: 122,372 metre (fixed, bisquare kernel)", fixed = TRUE)
  expect_no_match(txt, "e+05", fixed = TRUE)
  gwr$info$bandwidth <- 42; gwr$info$adaptive <- TRUE
  txt2 <- paste(utils::capture.output(print(gwr)), collapse = "\n")
  expect_match(txt2, "Bandwidth: 42 neighbours (adaptive, bisquare kernel)", fixed = TRUE)
})


# ---------------------------------------------------------------------------
# S10-DOCS-PKG-7: a character or factor response is refused with a message
# naming it, as README and getting-started say, before any bandwidth search.
# ---------------------------------------------------------------------------

test_that("fit_gwr_model() refuses a non-numeric response up front", {
  skip_if_not_installed("GWmodel"); skip_if_not_installed("sp")
  s <- surf_test_points(60)
  s$zc <- as.character(round(s$z, 3))
  r <- .r3_catch(fit_gwr_model(s, "zc", "w"))
  expect_type(r$value, "character")
  expect_match(r$value, "^fit_gwr_model\\(\\): response 'zc' is not numeric \\(it is character\\)")
  expect_length(r$warnings, 0L)       # no failed search, no fallback bandwidth
  s$zf <- factor(sample(c("lo", "hi"), nrow(s), replace = TRUE))
  rf <- .r3_catch(fit_gwr_model(s, "zf", "w"))
  expect_match(rf$value, "response 'zf' is not numeric \\(it is factor\\)")
  expect_length(rf$warnings, 0L)
})
