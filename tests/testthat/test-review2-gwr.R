# ===========================================================================
# GWR regressions from the second adversarial review.
# ===========================================================================

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
  # The global index is on the centred predictors: 1 for a single one.
  expect_equal(r$value$info$condition_index,
               .condition_index(cbind(1, d$soil - mean(d$soil))))
  expect_equal(r$value$info$condition_index, 1)
})


# ---------------------------------------------------------------------------
# The survey's kernel weights are GWmodel's, edge and zero width included, and
# the messages say what a singular window and a non-finite coefficient are.
# ---------------------------------------------------------------------------

test_that("the survey's kernel weights equal GWmodel::gw.weight()", {
  skip_if_not_installed("GWmodel")
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
  skip_if_not_installed("GWmodel")
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

test_that("an exactly singular window returns NA for what it cannot estimate, with the cause named", {
  # `soil` is a 0/1 indicator, constant inside every 20-neighbour window.  The
  # fit used to stop here ("inv(): matrix is singular").  The local solve
  # returns the fit: the slope of `soil` is NA in every window, the intercept
  # is NA where soil is 1 (it would absorb the slope there) and an estimate
  # where soil is 0, and the slope of `a` and every fitted value are
  # estimates.
  d <- .r2_clusters(exact = TRUE)
  r <- .r2_catch(fit_gwr_model(d, "z", c("a", "soil"), adaptive = TRUE,
                               bandwidth = 20))
  fit <- r$value
  expect_s3_class(fit, "gwr_fit")
  co <- coef(fit)
  expect_true(all(is.na(co$soil)))
  expect_identical(is.na(co$Intercept), d$soil == 1)
  expect_false(anyNA(co$a))
  expect_lt(max(abs(co$a - 2)), 0.5)
  expect_false(anyNA(fitted(fit)))
  expect_lt(sqrt(mean((fitted(fit) - (2 * d$a + 3 * d$soil))^2)), 0.15)
  expect_identical(fit$info$n_local_singular, 200L)
  expect_identical(fit$info$n_rank_deficient, 200L)
  expect_true(all(fit$info$local_collinearity$rank == 2L))
  w <- grep("local regressions are rank deficient", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "^fit_gwr_model\\(\\): 200 of 200 local regressions are rank deficient")
  expect_match(w, "Intercept at 100, a at 0, soil at 200")
  expect_match(w, "a predictor that is constant inside a window")
  # The collinearity warning before it says what such a window returns.
  expect_true(any(grepl("a coefficient that an exactly singular window cannot estimate is NA",
                        r$warnings)))
  expect_false(any(grepl("GWmodel stop", r$warnings)))
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
  # Of the models with a predictor, `a` alone is best, as it was.  The model
  # without predictors is scored too and ranks above it: at 6 neighbours each
  # window fits two parameters on five points, which costs more in AICc than
  # `a` explains.
  expect_identical(sel$table$variables[sel$table$n_vars > 0L][1L], "a")
  expect_identical(sel$table$variables[1L], "<none>")
  expect_identical(sel$best, character(0))
  # Seven models: the six of the sweep and the one without predictors, which
  # the table carries as `<none>`.
  expect_true(any(grepl("AICc is undefined for 1 of 7 model\\(s\\) at bandwidth 6",
                        r$warnings)))
  expect_identical(nrow(sel$table), 7L)
  expect_true(is.na(sel$table$criterion[7]))
  expect_identical(sel$table$n_vars[7], 3L)
  expect_true(all(is.finite(sel$table$criterion[1:6])))
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


# ---------------------------------------------------------------------------
# An adaptive bandwidth above n was capped at n in silence: a distance passed
# with adaptive left at TRUE.  GWmodel::bw.gwr()'s choice below 20 points,
# whose search range [20, n] is reversed, was capped the same way; the
# package's own search covers the range the model allows, from its floor to n.
# ---------------------------------------------------------------------------

test_that("an adaptive bandwidth above n is capped with a warning that suggests adaptive = FALSE", {
  set.seed(7)
  dd <- data.frame(x = runif(200, 0, 5000), y = runif(200, 0, 5000), a = rnorm(200))
  dd$z <- (1 + 3 * dd$x / 5000) * dd$a + rnorm(200, 0, 0.3)
  d <- sf::st_as_sf(dd, coords = c("x", "y"), crs = 32632)
  r <- .r2_catch(fit_gwr_model(d, "z", "a", bandwidth = 1500))
  expect_equal(r$value$info$bandwidth, 200)
  w <- grep("exceeds the 200 observations", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "adaptive bandwidth of 1500 neighbours.*set adaptive = FALSE")
  # n itself is a valid count and says nothing.
  r200 <- .r2_catch(fit_gwr_model(d, "z", "a", bandwidth = 200))
  expect_false(any(grepl("exceeds", r200$warnings)))

  s <- .r2_pts(14, 1, c("a", "b"))
  s$z <- 1 + s$a + rnorm(14, 0, 0.3)
  rs <- .r2_catch(gwr_model_selection(s, "z", c("a", "b"), bandwidth = 50))
  expect_identical(rs$value$bandwidth, 14L)
  expect_true(any(grepl("adaptive bandwidth of 50 neighbours exceeds the 14 observations",
                        rs$warnings)))
})

test_that("below 20 points the bandwidth is searched from the model's floor to n", {
  # GWmodel::bw.gwr() searches adaptive counts from 20 up, a range that is
  # empty for 12 observations: it returned 18, which was capped at 12 with a
  # warning that said neither was an optimised bandwidth.  The package's
  # search runs from the floor for two parameters (4 neighbours with the
  # bisquare kernel) to 12, so the bandwidth is a selected one and there is
  # nothing to cap.
  set.seed(1)
  n <- 12
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n)
  d <- sf::st_as_sf(data.frame(x = x, y = y, a = a,
                               z = 1 + (x / 1000) * a + rnorm(n, 0, 0.3)),
                    coords = c("x", "y"), crs = 32617)
  r <- .r2_catch(fit_gwr_model(d, "z", "a"))
  bs <- r$value$info$bandwidth_search
  expect_identical(bs$range, c(4, 12))
  expect_gte(r$value$info$bandwidth, 4)
  expect_lte(r$value$info$bandwidth, 12)
  expect_false(r$value$info$bandwidth_is_fallback)
  expect_false(any(grepl("bw\\.gwr|range that is empty|exceeds the 12 observations",
                         r$warnings)))
  # Every count in the range was a candidate, and the one used has the lowest
  # AICc of those scored.
  expect_identical(bs$curve$bandwidth[which.min(bs$curve$criterion)],
                   r$value$info$bandwidth)
  expect_equal(min(bs$curve$criterion), r$value$info$AICc, tolerance = 1e-8)

  s <- .r2_pts(14, 1, c("a", "b"))
  s$z <- 1 + s$a + rnorm(14, 0, 0.3)
  rs <- .r2_catch(gwr_model_selection(s, "z", c("a", "b")))
  expect_gte(rs$value$bandwidth, 5L)
  expect_lte(rs$value$bandwidth, 14L)
  expect_false(any(grepl("empty for 14 observations", rs$warnings)))
})


# ---------------------------------------------------------------------------
# A predictor named twice is one term, and the collinearity checks and the
# coefficient map treated it as two perfectly collinear ones.
# ---------------------------------------------------------------------------

test_that("a predictor named twice is fitted, surveyed and mapped once", {
  set.seed(1)
  n <- 120
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- sf::st_as_sf(data.frame(x = x, y = y, a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "y"), crs = 32617)
  d$z <- 1 + (x / 1000) * d$a + 0.5 * d$b + rnorm(n, 0, 0.3)
  dup <- .r2_catch(fit_gwr_model(d, "z", c("a", "b", "a"), bandwidth = 40))
  ok  <- .r2_catch(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 40))
  expect_identical(dup$warnings, ok$warnings)
  expect_false(any(grepl("singular|collinear", dup$warnings)))
  expect_identical(dup$value$predictor_vars, c("a", "b"))
  expect_equal(dup$value$info$condition_index, ok$value$info$condition_index)
  expect_equal(dup$value$info$n_local_collinear, 0L)
  expect_identical(coef(dup$value), coef(ok$value))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(dup$value, type = "coefficients"), "ggplot")
})
