# tests/testthat/test-fix-gwr-local-solve.R
# ---------------------------------------------------------------------------
# The local regressions of a GWR are solved by the package: one weighted
# least squares per window, by pivoted QR, with one tolerance.
#
# A window of deficient rank (a predictor constant among the weighted
# observations, exactly collinear predictors, fewer weighted observations
# than parameters) used to go to an inverse of X'WX, which stopped the whole
# fit or returned coefficients made of rounding error according to how the
# predictor was coded: an indicator at 0/1 stopped, the same indicator at
# 0.7/1.3 returned slopes of -117 to 271 around a true 5.  The tests below
# take each such case in several codings and require the same outcome for
# all of them: the fit is returned, a coefficient the window cannot estimate
# is NA, and a fitted or predicted value is returned when it is estimable.
#
# On windows of full rank the results are GWmodel's, to rounding error.
# ---------------------------------------------------------------------------

skip_if_not_installed("sp")

# Every warning raised while evaluating `expr`, muffled, with its class,
# beside the value (or the error message, as a character string).
.ls_catch <- function(expr) {
  w <- character(0); cls <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); cls <<- c(cls, class(cnd)[1L])
      invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w, classes = cls)
}

.ls_rmse <- function(e) sqrt(mean(e^2))

# `n` points over a square of side `extent`, `np` standard-normal predictors
# (a, b, ...) and a response whose slope on `a` drifts from west to east.
.ls_pts <- function(n, seed, np = 1L, extent = 1000) {
  set.seed(seed)
  x <- runif(n, 0, extent); y <- runif(n, 0, extent)
  df <- data.frame(x = 5e5 + x, y = 5e6 + y)
  for (j in seq_len(np)) df[[letters[j]]] <- rnorm(n)
  df$z <- 1 + (1 + x / extent) * df$a + rnorm(n, 0, 0.3)
  sf::st_as_sf(df, coords = c("x", "y"), crs = 32632)
}

# 300 points on a 1 km square.  `region` takes `lev[1]` west of x = 500 and
# `lev[2]` east of it, so a 40-neighbour window away from that line holds it
# constant.  The response is 1 + 2 a + 3 (east of the line) + noise of SD 0.3:
# the same points, predictor, noise and response for every `lev`, so two
# codings are one data set, and the slope of `region` times `step` (the
# distance between the two levels) is the same number in both.
.ls_regional <- function(lev) {
  set.seed(42); n <- 300
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n)
  east <- x > 500; eps <- rnorm(n, 0, 0.3)
  region <- ifelse(east, lev[2], lev[1])
  mu <- 1 + 2 * a + 3 * east
  d <- sf::st_as_sf(data.frame(x = x, y = y, a = a, region = region, z = mu + eps),
                    coords = c("x", "y"), crs = 32632)
  # Windows in which `region` is constant among the rows of weight above 0.
  flat <- vapply(seq_len(n), function(i) {
    w <- .gw_kernel_weights(sqrt((x - x[i])^2 + (y - y[i])^2), 40L, "bisquare", TRUE)
    length(unique(region[w > 0])) == 1L
  }, logical(1))
  list(d = d, mu = mu, flat = flat, east = east, x = x, y = y, region = region,
       step = lev[2] - lev[1])
}


# ---------------------------------------------------------------------------
# Full rank: GWmodel's numbers
# ---------------------------------------------------------------------------

test_that("a well-posed fit agrees with GWmodel::gwr.basic() in every column and diagnostic", {
  skip_if_not_installed("GWmodel")
  cases <- list(
    list(d = .ls_pts(150, 9), p = "a", bw = 30L, adaptive = TRUE, kernel = "bisquare"),
    list(d = .ls_pts(200, 4, 2L), p = c("a", "b"), bw = 200, adaptive = FALSE, kernel = "gaussian"),
    list(d = .ls_pts(250, 5, 3L), p = c("a", "b", "c"), bw = 40L, adaptive = TRUE, kernel = "tricube"),
    list(d = .ls_pts(180, 7), p = "a", bw = 300, adaptive = FALSE, kernel = "bisquare"),
    list(d = .ls_pts(160, 8, 2L), p = c("a", "b"), bw = 25L, adaptive = TRUE, kernel = "exponential"),
    list(d = .ls_pts(170, 6, 2L), p = c("a", "b"), bw = 20L, adaptive = TRUE, kernel = "boxcar"))
  for (i in seq_along(cases)) {
    cs <- cases[[i]]; lab <- paste("case", i, cs$kernel)
    r <- .ls_catch(fit_gwr_model(cs$d, "z", cs$p, adaptive = cs$adaptive,
                                 bandwidth = cs$bw, kernel = cs$kernel))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    expect_false("spatialkit_gwr_rank_deficient" %in% r$classes, info = lab)
    g <- GWmodel::gwr.basic(
      formula = stats::reformulate(cs$p, "z"), data = .to_sp(cs$d, c("z", cs$p)),
      bw = cs$bw, kernel = cs$kernel, adaptive = cs$adaptive)
    ours <- fit$engine$SDF@data; theirs <- g$SDF@data
    # The same table, column for column: coefficients, fitted values,
    # residuals, studentised residuals, standard errors, t values, local R2.
    expect_identical(names(ours), names(theirs), info = lab)
    nm <- c("Intercept", cs$p)
    expect_equal(as.matrix(ours[nm]), as.matrix(theirs[nm]), tolerance = 1e-8, info = lab)
    expect_equal(ours$yhat, theirs$yhat, tolerance = 1e-10, info = lab)
    se <- paste0(nm, "_SE")
    expect_lt(max(abs(as.matrix(ours[se]) / as.matrix(theirs[se]) - 1)), 1e-8, label = lab)
    expect_equal(as.matrix(ours), as.matrix(theirs), tolerance = 1e-7, info = lab)
    # coef() and fitted() are that table's columns.
    expect_equal(as.matrix(coef(fit)), as.matrix(theirs[nm]), tolerance = 1e-8, info = lab)
    expect_equal(as.numeric(fitted(fit)), theirs$yhat, tolerance = 1e-10, info = lab)
    # The global numbers: trace of the hat matrix, effective parameters, RSS,
    # AIC, AICc, BIC.
    dg <- g$GW.diagnostic; n <- nrow(cs$d)
    expect_equal(fit$info$tr_S, .gwr_trace_s(dg$AIC, dg$RSS.gw, n), tolerance = 1e-6, info = lab)
    expect_equal(fit$info$enp, dg$enp, tolerance = 1e-8, info = lab)
    expect_equal(fit$info$edf, dg$edf, tolerance = 1e-8, info = lab)
    expect_equal(fit$info$AICc, dg$AICc, tolerance = 1e-9, info = lab)
    expect_equal(unlist(fit$engine$GW.diagnostic), unlist(dg), tolerance = 1e-8, info = lab)
    expect_equal(sum(fit$info$local_collinearity$leverage), fit$info$tr_S, info = lab)
    # Every window has full rank, and nothing is NA.
    expect_true(all(fit$info$local_collinearity$rank == length(nm)), info = lab)
    expect_identical(fit$info$n_rank_deficient, 0L, info = lab)
    expect_identical(fit$info$n_local_singular, 0L, info = lab)
  }
})

test_that("predictions on full-rank windows agree with GWmodel::gwr.predict(), variance included", {
  skip_if_not_installed("GWmodel")
  d <- .ls_pts(200, 4, 2L)
  set.seed(3)
  nd <- sf::st_as_sf(data.frame(x = 5e5 + runif(40, 100, 900), y = 5e6 + runif(40, 100, 900),
                                a = rnorm(40), b = rnorm(40)),
                     coords = c("x", "y"), crs = 32632)
  for (cs in list(list(bw = 35L, adaptive = TRUE, kernel = "bisquare"),
                  list(bw = 250, adaptive = FALSE, kernel = "gaussian"))) {
    fit <- suppressWarnings(suppressMessages(
      fit_gwr_model(d, "z", c("a", "b"), adaptive = cs$adaptive, bandwidth = cs$bw,
                    kernel = cs$kernel)))
    g <- .gwr_quietly(GWmodel::gwr.predict(
      z ~ a + b, data = .to_sp(d, c("z", "a", "b")),
      predictdata = .to_sp(nd, c("a", "b")), bw = cs$bw, kernel = cs$kernel,
      adaptive = cs$adaptive))
    p <- predict(fit, newdata = nd, se.fit = TRUE)
    expect_equal(p$fit, g$SDF$prediction, tolerance = 1e-9, info = cs$kernel)
    expect_identical(predict(fit, newdata = nd), p$fit)
    expect_equal(p$se.pred^2, g$SDF$prediction_var, tolerance = 1e-8, info = cs$kernel)
    expect_true(all(p$se.fit < p$se.pred))
    expect_equal(p$residual.scale^2, fit$info$sigma2)
    expect_equal(p$df, fit$info$edf)
  }
})

test_that("predict() at the training rows is fitted(), to the last bit, for every kernel", {
  d <- .ls_pts(120, 3, 2L)
  for (kernel in c("bisquare", "gaussian", "tricube", "boxcar", "exponential"))
    for (ad in c(TRUE, FALSE)) {
      fit <- suppressWarnings(suppressMessages(
        fit_gwr_model(d, "z", c("a", "b"), adaptive = ad, kernel = kernel,
                      bandwidth = if (ad) 25 else 350)))
      expect_false(anyNA(fitted(fit)))
      expect_identical(predict(fit, newdata = d), unname(fitted(fit)),
                       info = paste(kernel, ad))
    }
  # And where windows are of deficient rank.
  rg <- .ls_regional(c(0.7, 1.3))
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(rg$d, "z", c("a", "region"), bandwidth = 40)))
  expect_gt(fit$info$n_rank_deficient, 100L)
  expect_identical(predict(fit, newdata = rg$d), unname(fitted(fit)))
})

test_that("the window helper returns the rows and weights of the full kernel vector", {
  set.seed(1)
  for (rep in 1:60) {
    n <- sample(c(3L, 10L, 60L), 1L)
    d <- c(0, round(runif(n - 1L, 0, 300), sample(0:2, 1L)))
    if (rep %% 3L == 0L) d[2:3] <- 0                # points at the regression point
    if (rep %% 4L == 0L) d <- pmin(d, 100)          # ties at the kernel's edge
    for (k in c("bisquare", "gaussian", "tricube", "boxcar", "exponential"))
      for (ad in c(TRUE, FALSE)) {
        bw <- if (ad) sample(seq_len(n + 2L), 1L) else sample(c(1e-9, 50, 100, 1e4), 1L)
        w   <- .gw_kernel_weights(d, bw, k, ad)
        win <- .gw_window(d, .gw_kernel_width(d, bw, ad), k)
        expect_identical(win$rows, which(w > 0))
        expect_identical(win$w, w[which(w > 0)])
      }
  }
})


# ---------------------------------------------------------------------------
# A predictor constant inside a window, in three codings
# ---------------------------------------------------------------------------

test_that("a regional covariate gives the same fit coded 0/1, 0.7/1.3 and 288.7/289.3", {
  # Before the local solve: coded 0/1 the fit stopped ("inv(): matrix is
  # singular"); coded 0.7/1.3 it returned with fitted values 5.3 off the
  # noise-free mean in the flat windows (root mean square) and an in-sample
  # R2 of -1.9; coded 288.7/289.3 they were 13.9 off and R2 was -18.9.
  fits <- list(); sets <- list()
  for (lev in list(c(0, 1), c(0.7, 1.3), c(288.7, 289.3))) {
    rg <- .ls_regional(lev); lab <- paste(lev, collapse = "/")
    r <- .ls_catch(fit_gwr_model(rg$d, "z", c("a", "region"), bandwidth = 40))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    expect_identical(sum(rg$flat), 185L)
    co <- coef(fit); fv <- fitted(fit)
    # The slope of `region` is NA exactly in the flat windows.
    expect_identical(is.na(co$region), rg$flat, info = lab)
    # So is the intercept, which would absorb 3 * region there, unless
    # `region` is 0 in the window: then it is identified, and returned.
    expect_identical(is.na(co$Intercept), rg$flat & rg$region != 0, info = lab)
    # The slope of `a` is an estimate everywhere, near its true 2.
    expect_false(anyNA(co$a))
    expect_lt(abs(mean(co$a[rg$flat]) - 2), 0.05)
    # Every fitted value is estimable, and close to the noise-free mean.
    expect_false(anyNA(fv))
    expect_lt(.ls_rmse(fv[rg$flat] - rg$mu[rg$flat]), 0.15)
    expect_lt(.ls_rmse(fv[!rg$flat] - rg$mu[!rg$flat]), 0.15)
    expect_gt(model_metrics(fit)$R2, 0.95)
    # Counted, placed, and warned about with a class.
    expect_identical(fit$info$n_rank_deficient, 185L, info = lab)
    expect_identical(fit$info$n_local_singular, 185L, info = lab)
    expect_identical(fit$info$local_collinearity$rank, ifelse(rg$flat, 2L, 3L), info = lab)
    expect_identical(unname(fit$info$nonfinite_coef[, "region"]), rg$flat, info = lab)
    expect_identical(sum(r$classes == "spatialkit_gwr_rank_deficient"), 1L, info = lab)
    w <- r$warnings[r$classes == "spatialkit_gwr_rank_deficient"]
    expect_match(w, "^fit_gwr_model\\(\\): 185 of 300 local regressions are rank deficient at the bandwidth in use \\(40 neighbours, bisquare kernel\\)")
    expect_match(w, sprintf("Intercept at %d, a at 0, region at 185", sum(rg$flat & rg$region != 0)))
    expect_match(w, "GWmodel::gwr.mixed\\(\\)")
    expect_no_match(w, "No bandwidth was supplied")
    # Standard errors: NA with their coefficient, estimates elsewhere.
    cs <- coef(fit, se = TRUE)
    expect_identical(is.na(cs$se$region), rg$flat, info = lab)
    expect_false(anyNA(cs$se$a))
    fits[[lab]] <- fit; sets[[lab]] <- rg
  }
  # One fit, three codings: the same ranks, trace, AICc, slopes of `a`, and
  # fitted values up to the shift of the response's level.
  ref <- fits[[1L]]
  # The two numbers are this fit's own, recorded when the solve was written:
  # they hold the solve to what it returned then.  That the three codings
  # agree is tested below without them.
  expect_equal(ref$info$tr_S, 42.2988, tolerance = 1e-5)
  expect_equal(ref$info$AICc, 145.2796, tolerance = 1e-5)
  for (lab in names(fits)[-1L]) {
    f <- fits[[lab]]
    expect_identical(f$info$local_collinearity$rank, ref$info$local_collinearity$rank)
    expect_equal(f$info$tr_S, ref$info$tr_S, tolerance = 1e-9, info = lab)
    expect_equal(f$info$AICc, ref$info$AICc, tolerance = 1e-9, info = lab)
    expect_equal(coef(f)$a, coef(ref)$a, tolerance = 1e-8, info = lab)
    expect_equal(fitted(f), fitted(ref), tolerance = 1e-8, info = lab)
    keep <- !sets[[lab]]$flat
    expect_equal(coef(f)$region[keep] * sets[[lab]]$step,
                 coef(ref)$region[keep] * sets[[1L]]$step, tolerance = 1e-7, info = lab)
  }
})

test_that("the rank does not depend on a predictor's origin, a million units out included", {
  # A rank decision on the uncentred design depends on where 0 is: at
  # 1000000.7 / 1000001.3 it called 30 of the 115 windows in which `region`
  # does vary rank deficient.
  rg0 <- .ls_regional(c(0, 1))
  rg6 <- .ls_regional(c(1e6 + 0.7, 1e6 + 1.3))
  f0 <- suppressWarnings(suppressMessages(fit_gwr_model(rg0$d, "z", c("a", "region"), bandwidth = 40)))
  f6 <- suppressWarnings(suppressMessages(fit_gwr_model(rg6$d, "z", c("a", "region"), bandwidth = 40)))
  expect_identical(f6$info$local_collinearity$rank, f0$info$local_collinearity$rank)
  expect_equal(f6$info$tr_S, f0$info$tr_S, tolerance = 1e-7)
  expect_equal(fitted(f6), fitted(f0), tolerance = 1e-6)
  # A temperature in degrees C and in kelvin: the same slopes and fitted values.
  set.seed(1); n <- 200
  x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4)
  tC <- 15 + 4 * sin(x / 3000) + 3 * cos(y / 2500) + rnorm(n)
  d <- sf::st_as_sf(data.frame(x = x, y = y, tC = tC, tK = tC + 273.15),
                    coords = c("x", "y"), crs = 32632)
  d$yield <- 50 - 0.8 * d$tC + rnorm(n)
  fC <- suppressMessages(fit_gwr_model(d, "yield", "tC", bandwidth = 30))
  fK <- suppressMessages(fit_gwr_model(d, "yield", "tK", bandwidth = 30))
  expect_equal(coef(fK)$tK, coef(fC)$tC, tolerance = 1e-10)
  expect_equal(fitted(fK), fitted(fC), tolerance = 1e-10)
  expect_equal(coef(fK, se = TRUE)$se$tK, coef(fC, se = TRUE)$se$tC, tolerance = 1e-9)
})

test_that("predict() is finite where the new row's region is its window's, and NA where it is not", {
  # Before: coded 0.7/1.3, predict() returned a number for every row, 119
  # off at worst, with no warning.
  set.seed(7); m <- 400
  nx <- runif(m, 0, 1000); ny <- runif(m, 0, 1000); na_ <- rnorm(m)
  out <- list()
  for (lev in list(c(0, 1), c(0.7, 1.3), c(288.7, 289.3))) {
    rg <- .ls_regional(lev); lab <- paste(lev, collapse = "/")
    fit <- suppressWarnings(suppressMessages(
      fit_gwr_model(rg$d, "z", c("a", "region"), bandwidth = 40)))
    # The value `region` is constant at among the weighted training rows of
    # each new location's window, or NA where it varies there.
    held <- vapply(seq_len(m), function(j) {
      w <- .gw_kernel_weights(sqrt((rg$x - nx[j])^2 + (rg$y - ny[j])^2), 40L, "bisquare", TRUE)
      u <- unique(rg$region[w > 0]); if (length(u) == 1L) u else NA_real_
    }, numeric(1))
    expect_identical(sum(!is.na(held)), 234L)
    # (a) Every new row carries the region of where it lies: all estimable.
    nd <- sf::st_as_sf(data.frame(x = nx, y = ny, a = na_,
                                  region = ifelse(nx > 500, lev[2], lev[1])),
                       coords = c("x", "y"), crs = 32632)
    p <- expect_no_warning(predict(fit, newdata = nd))
    expect_false(anyNA(p))
    truth <- 1 + 2 * na_ + 3 * (nx > 500)
    expect_lt(.ls_rmse(p - truth), 0.15)
    expect_lt(max(abs(p - truth)), 1)
    # (b) The boundary of the new rows' region moved 300 m east: 53 rows now
    # carry a value their window's training rows do not have.  Those, and
    # only those, are NA.
    nd2 <- nd; nd2$region <- ifelse(nx > 800, lev[2], lev[1])
    unseen <- !is.na(held) & held != nd2$region
    expect_identical(sum(unseen), 53L)
    r <- .ls_catch(predict(fit, newdata = nd2))
    expect_identical(is.na(r$value), unseen, info = lab)
    expect_identical(r$classes, "spatialkit_gwr_not_estimable")
    expect_match(r$warnings, "^predict.gwr_fit\\(\\): 53 of 400 location\\(s\\) have no estimable local regression")
    expect_match(r$warnings, "At 0 of them no training point carries weight")
    expect_match(r$warnings, "At 53 the training points in the window cannot place the new row")
    truth2 <- 1 + 2 * na_ + 3 * (nx > 800)
    expect_lt(.ls_rmse((r$value - truth2)[!unseen]), 0.2)
    # A standard error is NA with its prediction.
    ps <- suppressWarnings(predict(fit, newdata = nd2, se.fit = TRUE))
    expect_identical(is.na(ps$se.fit), unseen)
    expect_identical(is.na(ps$se.pred), unseen)
    out[[lab]] <- p
  }
  expect_equal(out[[2L]], out[[1L]], tolerance = 1e-8)
  expect_equal(out[[3L]], out[[1L]], tolerance = 1e-8)
})


# ---------------------------------------------------------------------------
# Fewer weighted observations than parameters, and exact collinearity
# ---------------------------------------------------------------------------

# Fifty points in a 500 m square and one more 400 m beyond its corner, whose
# predictor value is `a_far`.  With a fixed bisquare bandwidth of 150 the far
# point's window holds the far point alone: one weighted point for two
# parameters.
.ls_far_point <- function(a_far) {
  set.seed(1); n <- 50
  x <- c(runif(n, 0, 500), 900); y <- c(runif(n, 0, 500), 900)
  d <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, a = c(rnorm(n), a_far)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + c(rnorm(n, 0, 0.3), 0)
  d
}

test_that("a window with one weighted point for two parameters returns its fitted value and no slope", {
  # Before: with the far point's predictor at 0.5 or 0 the fit stopped; at
  # 0.3, 0.7 and 1.3 it returned local coefficients of (4, -8), (19.6, -156)
  # and (9.0, -13.0) there, for a true (1, 2).
  res <- lapply(c(0.5, 0.3, 0.7, 1.3, 0), function(a_far) {
    d <- .ls_far_point(a_far)
    r <- .ls_catch(fit_gwr_model(d, "z", "a", adaptive = FALSE, bandwidth = 150))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    lc <- fit$info$local_collinearity
    expect_identical(lc$rank[51], 1L)
    expect_true(all(lc$rank[-51] == 2L))
    expect_equal(lc$leverage[51], 1)
    expect_equal(fitted(fit)[51], d$z[51])
    co <- coef(fit)
    expect_true(is.na(co$a[51]))
    # One point at a = a_far says what the line is at a_far and nothing about
    # where it crosses a = 0, unless a_far is 0.
    expect_identical(is.na(co$Intercept[51]), a_far != 0)
    if (a_far == 0) expect_equal(co$Intercept[51], d$z[51])
    expect_false(anyNA(co[-51, ]))
    expect_identical(fit$info$n_rank_deficient, 1L)
    expect_identical(fit$info$n_saturated, 1L)
    w <- r$warnings[r$classes == "spatialkit_gwr_rank_deficient"]
    expect_length(w, 1L)
    expect_match(w, "^fit_gwr_model\\(\\): 1 of 51 local regressions are rank deficient at the bandwidth in use \\(150, bisquare kernel\\)")
    expect_match(w, "a bandwidth that leaves a window fewer observations than parameters")
    expect_no_match(w, "No bandwidth was supplied")
    fit
  })
  # The other fifty windows do not hold the far point and are the same fit
  # (to rounding: the far point's value moves the predictor's overall scale).
  for (f in res[-1L])
    expect_equal(coef(f)[-51, ], coef(res[[1L]])[-51, ], tolerance = 1e-12)
})

test_that("exactly collinear predictors: no slope is identified, every fitted value is", {
  set.seed(3); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n)
  d <- sf::st_as_sf(data.frame(x = x, y = y, a = a, b = 2 * a + 1, k = 5),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + a + rnorm(n, 0, 0.1)
  # b is 2 a + 1 everywhere: neither slope has one value, and neither has the
  # intercept (the 1 in b can move into it).  Before, this fit stopped.
  r <- .ls_catch(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 30))
  fit <- r$value
  expect_s3_class(fit, "gwr_fit")
  expect_true(all(is.na(as.matrix(coef(fit)))))
  expect_false(anyNA(fitted(fit)))
  expect_lt(.ls_rmse(fitted(fit) - (1 + a)), 0.1)
  expect_true(all(fit$info$local_collinearity$rank == 2L))
  expect_true("spatialkit_gwr_rank_deficient" %in% r$classes)
  # A predictor constant over the whole layer: its slope and the intercept go,
  # the slope of `a` stays.  Before, this fit returned slopes of -11 to 16
  # for `k` and fitted values 5 off the mean.
  rk <- .ls_catch(fit_gwr_model(d, "z", c("a", "k"), bandwidth = 30))
  ck <- coef(rk$value)
  expect_true(all(is.na(ck$k)) && all(is.na(ck$Intercept)))
  expect_false(anyNA(ck$a))
  expect_equal(fitted(rk$value),
               fitted(suppressMessages(fit_gwr_model(d, "z", "a", bandwidth = 30))),
               tolerance = 1e-10)
})

# Four clusters of 50, 400 m apart; `soil` is exactly constant inside each.
.ls_clusters <- function(lev) {
  set.seed(5)
  cl <- rep(1:4, each = 50)
  d <- sf::st_as_sf(data.frame(
    x = c(runif(50, 0, 100), runif(50, 400, 500), runif(50, 0, 100), runif(50, 400, 500)),
    y = c(runif(50, 0, 100), runif(50, 0, 100), runif(50, 400, 500), runif(50, 400, 500)),
    a = rnorm(200), soil = rep(lev, 2L)[cl]),
    coords = c("x", "y"), crs = 32632)
  d$z <- 2 * d$a + 3 * (d$soil - lev[1]) + rnorm(200, 0, 0.2)
  d
}

test_that("a predictor constant in every window gives the same fit at 0/1 and at 0.3/1.3", {
  # Before: 0/1 stopped; 0.3/1.3 returned slopes of -461 to 186 around a true 3.
  f01 <- suppressWarnings(suppressMessages(
    fit_gwr_model(.ls_clusters(c(0, 1)), "z", c("a", "soil"), bandwidth = 20)))
  f03 <- suppressWarnings(suppressMessages(
    fit_gwr_model(.ls_clusters(c(0.3, 1.3)), "z", c("a", "soil"), bandwidth = 20)))
  for (f in list(f01, f03)) {
    expect_true(all(is.na(coef(f)$soil)))
    expect_false(anyNA(coef(f)$a))
    expect_false(anyNA(fitted(f)))
    expect_identical(f$info$n_rank_deficient, 200L)
  }
  expect_equal(coef(f03)$a, coef(f01)$a, tolerance = 1e-9)
  expect_equal(fitted(f03), fitted(f01), tolerance = 1e-9)
  expect_identical(sum(is.na(coef(f01)$Intercept)), 100L)   # the soil = 1 clusters
  expect_identical(sum(is.na(coef(f03)$Intercept)), 200L)
})

test_that("a zero-width adaptive kernel leaves a window no rows: rank 0, nothing estimable", {
  set.seed(2)
  sx <- runif(30, 0, 1000); sy <- runif(30, 0, 1000)
  dd <- data.frame(x = rep(sx, each = 4), y = rep(sy, each = 4))
  dd$a <- rnorm(120); dd$z <- 1 + 2 * dd$a + rnorm(120, 0, 0.3)
  d <- sf::st_as_sf(dd, coords = c("x", "y"), crs = 32632)
  r <- .ls_catch(fit_gwr_model(d, "z", "a", bandwidth = 4, kernel = "gaussian"))
  fit <- r$value
  expect_identical(fit$info$local_collinearity$rank, rep(0L, 120L))
  expect_identical(fit$info$n_rank_deficient, 120L)
  expect_true(all(is.na(fitted(fit))))
  expect_identical(fit$info$n_local_singular, 120L)
  expect_true(is.na(fit$info$AICc))
  expect_true(any(grepl("4 or more observations share one location", r$warnings)))
  # Not the rank-deficiency warning: these windows have no fitted value.
  expect_false("spatialkit_gwr_rank_deficient" %in% r$classes)
})


# ---------------------------------------------------------------------------
# Gaussian and exponential kernels: every point has weight, none is cut off
# ---------------------------------------------------------------------------

# Sixty points in a 500 m square and, 400 m and more from them, two points U
# and T on one line, `gap` metres apart.
.ls_pair <- function(gap) {
  set.seed(7); n <- 60
  xy <- rbind(cbind(runif(n, 0, 500), runif(n, 0, 500)), c(900 - gap, 900), c(900, 900))
  d <- sf::st_as_sf(data.frame(x = 5e5 + xy[, 1], y = 5e6 + xy[, 2], a = rnorm(n + 2L)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n + 2L, 0, 0.3)
  d
}

# Weighted least squares at row `i` by another route: LINPACK's QR
# (lm.wfit) on the predictors centred at row i's own values, so that its
# intercept is the fitted value there.
.ls_wls_at <- function(d, i, vars, w) {
  X <- as.matrix(sf::st_drop_geometry(d)[, vars, drop = FALSE])
  Xc <- cbind(1, sweep(X, 2L, X[i, ]))
  f <- stats::lm.wfit(Xc, d$z, w / max(w), tol = 1e-13)
  list(fitted = unname(f$coefficients[1L]), slopes = unname(f$coefficients[-1L]),
       rank = f$rank)
}

test_that("an isolated pair under a Gaussian kernel is fitted, with the least-squares coefficients", {
  # A count of the points whose weight is above 1e-8 of the window's largest
  # calls U's and T's windows one-point windows at a gap of 310 m (6.2
  # bandwidths) and would refuse the fit.  The windows have full rank: every
  # other point enters with its weight, however small against the pair's.
  for (gap in c(300, 310)) {
    d <- .ls_pair(gap)
    r <- .ls_catch(fit_gwr_model(d, "z", "a", adaptive = FALSE, bandwidth = 50,
                                 kernel = "gaussian"))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    expect_false("spatialkit_gwr_rank_deficient" %in% r$classes)
    expect_true(all(fit$info$local_collinearity$rank == 2L))
    xy <- sf::st_coordinates(d)
    for (i in 61:62) {
      w <- exp(-0.5 * ((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2) / 50^2)
      ref <- .ls_wls_at(d, i, "a", w)
      expect_identical(ref$rank, 2L)
      expect_equal(fitted(fit)[i], ref$fitted, tolerance = 1e-9, info = paste(gap, i))
      expect_equal(coef(fit)$a[i], ref$slopes, tolerance = 1e-7, info = paste(gap, i))
    }
  }
})

test_that("an observation many Gaussian bandwidths from the rest keeps its fitted value and loses its slope", {
  # One observation `far` metres beyond 60 others, fixed Gaussian 100.  At
  # 3 and 6.5 bandwidths its window has full rank and the coefficients are the
  # weighted least squares.  At 10 and 15 the others' weights are 1e-23 and
  # less of its own: by the tolerance they say nothing about the slope, which
  # is NA, and the fitted value is the observation.  Before, the solver
  # returned (0.88, 1.52) there at both distances, which no window of the
  # point supports, and a fitted value 0.17 from the observation.
  one <- function(far) {
    set.seed(1); m <- 60
    x <- c(runif(m, 0, 500), 250); y <- c(runif(m, 0, 500), 500 + far)
    d <- sf::st_as_sf(data.frame(x = x, y = y, a = rnorm(m + 1L)),
                      coords = c("x", "y"), crs = 32632)
    d$z <- 1 + 2 * d$a + rnorm(m + 1L, 0, 0.3)
    list(d = d, w = exp(-0.5 * ((x - 250)^2 + (y - 500 - far)^2) / 100^2))
  }
  for (far in c(310, 651)) {
    s <- one(far)
    fit <- suppressWarnings(suppressMessages(
      fit_gwr_model(s$d, "z", "a", adaptive = FALSE, bandwidth = 100, kernel = "gaussian")))
    ref <- .ls_wls_at(s$d, 61L, "a", s$w)
    expect_identical(fit$info$local_collinearity$rank[61], 2L)
    expect_equal(coef(fit)$a[61], ref$slopes, tolerance = 1e-8, info = far)
    expect_equal(fitted(fit)[61], ref$fitted, tolerance = 1e-10, info = far)
  }
  for (far in c(1026, 1523)) {
    s <- one(far)
    r <- .ls_catch(fit_gwr_model(s$d, "z", "a", adaptive = FALSE, bandwidth = 100,
                                 kernel = "gaussian"))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    expect_identical(fit$info$local_collinearity$rank[61], 1L)
    expect_true(is.na(coef(fit)$a[61]))
    expect_equal(fitted(fit)[61], s$d$z[61], tolerance = 1e-12)
    expect_equal(fit$info$local_collinearity$leverage[61], 1, tolerance = 1e-12)
    expect_false(anyNA(coef(fit)$a[-61]))
    expect_identical(sum(r$classes == "spatialkit_gwr_rank_deficient"), 1L)
  }
})

test_that("a regional covariate under a fixed Gaussian or exponential kernel is fitted, whatever its origin", {
  # Fixed Gaussian 60 and exponential 20 on the regional layer.  Before, the
  # Gaussian fit coded 0/1.3 was 0.18 off the noise-free mean and the one
  # coded 288/289.3 was 2.6 off, by the origin alone; both are 0.13 off.
  for (cs in list(list(kernel = "gaussian", bw = 60, n_def = 15L),
                  list(kernel = "exponential", bw = 20, n_def = 0L))) {
    fits <- lapply(list(c(0, 1.3), c(288, 289.3)), function(lev) {
      rg <- .ls_regional(lev)
      r <- .ls_catch(fit_gwr_model(rg$d, "z", c("a", "region"), adaptive = FALSE,
                                   bandwidth = cs$bw, kernel = cs$kernel))
      expect_s3_class(r$value, "gwr_fit")
      expect_false(anyNA(fitted(r$value)))
      expect_lt(.ls_rmse(fitted(r$value) - rg$mu), 0.25)
      list(fit = r$value, rg = rg)
    })
    f0 <- fits[[1L]]$fit; f288 <- fits[[2L]]$fit; rg <- fits[[1L]]$rg
    rk <- f0$info$local_collinearity$rank
    expect_identical(f288$info$local_collinearity$rank, rk, info = cs$kernel)
    expect_identical(sum(rk < 3L), cs$n_def, info = cs$kernel)
    expect_equal(fitted(f288), fitted(f0), tolerance = 1e-7, info = cs$kernel)
    expect_equal(coef(f288)$region, coef(f0)$region, tolerance = 1e-6, info = cs$kernel)
    # Where the window has full rank the coefficients are the weighted least
    # squares, also where the other region's rows carry little weight.
    xy <- sf::st_coordinates(rg$d)
    kern <- function(dd) if (cs$kernel == "gaussian") exp(-0.5 * dd^2 / cs$bw^2) else exp(-dd / cs$bw)
    for (i in c(which(rk == 3L)[c(1, 50, 120)], which.max(xy[, 1] * (rk == 3L)))) {
      ref <- .ls_wls_at(rg$d, i, c("a", "region"),
                        kern(sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)))
      expect_equal(fitted(f0)[i], ref$fitted, tolerance = 1e-9, info = paste(cs$kernel, i))
      expect_equal(unlist(coef(f0)[i, c("a", "region")], use.names = FALSE), ref$slopes,
                   tolerance = 1e-6, info = paste(cs$kernel, i))
    }
    # Where it has not, the other region's rows are too far to say anything
    # about `region`: its slope is NA and the fitted value is returned.
    expect_true(all(is.na(coef(f0)$region[rk < 3L])))
  }
})

test_that("a Gaussian window that holds one region predicts NA for a row of the other", {
  # Fixed Gaussian 40.  A new location more than about 330 m from the line
  # x = 500 has a window of rank 2: the other region's training rows are too
  # far, by the tolerance, to say anything about `region`.  Before, a row
  # that carries the eastern value was predicted there (west of x = 170) 2.6
  # below the model's mean on average and up to 28 off, with no warning.  It
  # is NA; a row with the window's own value is predicted.
  rg <- .ls_regional(c(0.7, 1.3))
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(rg$d, "z", c("a", "region"), adaptive = FALSE, bandwidth = 40,
                  kernel = "gaussian")))
  set.seed(7); m <- 400
  nx <- runif(m, 0, 1000); ny <- runif(m, 0, 1000); na_ <- rnorm(m)
  mk <- function(region) sf::st_as_sf(data.frame(x = nx, y = ny, a = na_, region = region),
                                      coords = c("x", "y"), crs = 32632)
  own <- ifelse(nx > 500, 1.3, 0.7)
  p_own <- expect_no_warning(predict(fit, newdata = mk(own)))
  expect_false(anyNA(p_own))
  # Every row given the eastern value: in the west that is the other region's.
  sp_train <- .to_sp(rg$d, c("z", "a", "region"))
  sv <- .gwr_predict_at(z ~ a + region, sp_train, data.frame(a = na_, region = 1.3),
                        cbind(nx, ny), 40, "gaussian", FALSE)
  west2 <- sv$rank == 2L & nx < 500
  expect_gt(sum(west2), 30L)
  r <- .ls_catch(predict(fit, newdata = mk(rep(1.3, m))))
  expect_identical(is.na(r$value), west2)
  expect_identical(r$classes, "spatialkit_gwr_not_estimable")
  # In the eastern windows of rank 2 the value is the window's own: predicted.
  expect_false(anyNA(r$value[sv$rank == 2L & nx > 500]))
  expect_identical(r$value[nx > 500], p_own[nx > 500])
})

test_that("far beyond the data a Gaussian kernel still predicts, on weights that compare", {
  # 2.5 km from the data with a bandwidth of 200 no weight reaches 1e-8, yet
  # many are within a few orders of the largest, and the solve depends only
  # on how they compare.
  d <- .ls_pts(200, 4, 2L)
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("a", "b"), adaptive = FALSE, bandwidth = 200, kernel = "gaussian")))
  nd <- sf::st_as_sf(data.frame(x = 5e5 + 3000, y = 5e6 + 500, a = 0.4, b = -0.2),
                     coords = c("x", "y"), crs = 32632)
  xy <- sf::st_coordinates(d)
  w <- exp(-0.5 * ((xy[, 1] - (5e5 + 3000))^2 + (xy[, 2] - (5e6 + 500))^2) / 200^2)
  expect_lt(max(w), 1e-8)
  r <- .ls_catch(predict(fit, newdata = nd))
  wls <- stats::lm.wfit(cbind(1, d$a, d$b), d$z, w / max(w))$coefficients
  expect_equal(r$value, sum(wls * c(1, 0.4, -0.2)), tolerance = 1e-8)
  # It is an extrapolation, and said to be one.
  expect_identical(r$classes, "spatialkit_gwr_extrapolation")
})

test_that("600 locations east of a regional layer are all predicted, as they were", {
  # Fixed Gaussian 110; every new row carries the eastern region's value.  In
  # 358 of the windows the training rows that count share one value of
  # `region`; the predictions there were sound before (within 0.43 of the
  # truth) and must stay finite and as close.
  set.seed(8); n <- 300; x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, a = rnorm(n),
                               region = 1.3 * (x > 500)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + 3 * d$region + rnorm(n, 0, 0.3)
  fit <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("a", "region"), adaptive = FALSE, bandwidth = 110,
                  kernel = "gaussian")))
  set.seed(4); m <- 600
  nd <- data.frame(x = 5e5 + runif(m, 1000, 1500), y = 5e6 + runif(m, 0, 1000),
                   a = rnorm(m), region = 1.3)
  r <- .ls_catch(predict(fit, newdata = sf::st_as_sf(nd, coords = c("x", "y"), crs = 32632)))
  expect_false(anyNA(r$value))
  expect_false("spatialkit_gwr_not_estimable" %in% r$classes)
  err <- abs(r$value - (1 + 2 * nd$a + 3 * 1.3))
  expect_lt(max(err), 0.43)
  expect_lt(stats::median(err), 0.05)
})


# ---------------------------------------------------------------------------
# Cross-validation folds are fitted and scored
# ---------------------------------------------------------------------------

test_that("cv_gwr() with strips held out scores every fold near the noise SD, whatever the coding", {
  # Four strips in x.  Before: coded 0/1 all four folds failed and no score
  # came back; coded 0.7/1.3 the fold RMSEs were 0.53, 0.29, 0.33 and 1.15
  # against a noise SD of 0.3.
  score <- lapply(list(c(0, 1), c(0.7, 1.3)), function(lev) {
    rg <- .ls_regional(lev)
    lab <- as.integer(cut(rg$x, c(-1, 250, 500, 750, 1001)))
    cv <- suppressWarnings(suppressMessages(
      cv_gwr(rg$d, "z", c("a", "region"), folds = lab, bandwidth = 40)))
    expect_identical(cv$fold_status$status, rep("ok", 4L))
    expect_identical(cv$overall$n_pred, 300L)
    expect_true(all(cv$fold_metrics$RMSE > 0.2 & cv$fold_metrics$RMSE < 0.4))
    cv$fold_metrics$RMSE
  })
  expect_equal(score[[2L]], score[[1L]], tolerance = 1e-7)
  expect_equal(score[[1L]], c(0.297, 0.291, 0.329, 0.340), tolerance = 0.01)
})

test_that("cv_gwr() scores five of five folds under Gaussian and exponential kernels", {
  # Random folds.  An isolated observation under fixed Gaussian 100, and the
  # regional covariate under fixed Gaussian 60 and exponential 20: a rule that
  # refuses a fit for one thin window loses most or all of these folds.
  # Before, all folds were scored, the Gaussian regional case at a pooled
  # RMSE of 0.88 against a noise SD of 0.3; it is 0.31.
  set.seed(11); n <- 60
  x <- c(runif(n, 0, 400), 400 + 620); y <- c(runif(n, 0, 400), 200); a <- c(rnorm(n), 0.3)
  d1 <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, a = a,
                                z = 1 + 2 * a + rnorm(n + 1L, 0, 0.3)),
                     coords = c("x", "y"), crs = 32632)
  cv1 <- suppressWarnings(suppressMessages(
    cv_gwr(d1, "z", "a", folds = make_folds(d1, k = 5, method = "random_kfold", seed = 1),
           adaptive = FALSE, bandwidth = 100, kernel = "gaussian")))
  expect_identical(sum(cv1$fold_status$status == "ok"), 5L)
  expect_lt(cv1$overall$RMSE, 0.45)
  rg <- .ls_regional(c(0, 1.3))
  fl <- make_folds(rg$d, k = 5, method = "random_kfold", seed = 1)
  for (cs in list(list(kernel = "gaussian", bw = 60), list(kernel = "exponential", bw = 20))) {
    cv <- suppressWarnings(suppressMessages(
      cv_gwr(rg$d, "z", c("a", "region"), folds = fl, adaptive = FALSE,
             bandwidth = cs$bw, kernel = cs$kernel)))
    expect_identical(sum(cv$fold_status$status == "ok"), 5L, info = cs$kernel)
    expect_lt(cv$overall$RMSE, 0.45)
    expect_gte(cv$overall$n_pred, 290L)
  }
})
