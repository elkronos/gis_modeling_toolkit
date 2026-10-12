# tests/testthat/test-fix-gwr-repeated-locations.R
# ---------------------------------------------------------------------------
# A layer with several observations per location.
#
# AICc counts the rows at one location as independent observations.  When a
# location has an effect of its own, a window that holds little but the
# location reproduces its mean, the residual sum of squares falls to the
# spread within locations, and AICc is lowest where the windows interpolate.
# On 100 sites with 3 observations each the default bandwidth search settled
# at 7 neighbours: every coefficient NA, every prediction at a new site NA.
#
# The search now minimises the leave-one-location-out sum of squares on such
# a layer, from the hat matrix; its range counts locations; a bandwidth at
# which a term is estimable nowhere cannot be selected; the fit says how many
# locations the observations share.
# ---------------------------------------------------------------------------

skip_if_not_installed("sp")

.rl_catch <- function(expr) {
  w <- character(0); cls <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); cls <<- c(cls, class(cnd)[1L])
      invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w, classes = cls)
}

# `sites` locations with `reps` observations each.  The predictors belong to
# the site, the slope of `a` drifts eastward, and each site has an effect of
# its own (SD 0.3) under measurement noise (SD 0.2).
.rl_sites <- function(sites, reps, seed) {
  set.seed(seed)
  x <- runif(sites, 0, 1000); y <- runif(sites, 0, 1000)
  a <- rnorm(sites); b <- rnorm(sites); u <- rnorm(sites, 0, 0.3)
  id <- rep(seq_len(sites), reps); n <- length(id)
  mu <- 1 + (1 + x[id] / 1000) * a[id] - 0.5 * b[id]
  sf::st_as_sf(data.frame(x = 5e5 + x[id], y = 5e6 + y[id], a = a[id], b = b[id],
                          z = mu + u[id] + rnorm(n, 0, 0.2), mu = mu, site = id),
               coords = c("x", "y"), crs = 32632)
}

# 40 new sites inside the study area, and the noise-free mean there.
.rl_new_sites <- function() {
  set.seed(99); m <- 40
  nd <- sf::st_as_sf(data.frame(x = 5e5 + runif(m, 100, 900), y = 5e6 + runif(m, 100, 900),
                                a = rnorm(m), b = rnorm(m)),
                     coords = c("x", "y"), crs = 32632)
  nd$truth <- 1 + (1 + (sf::st_coordinates(nd)[, 1] - 5e5) / 1000) * nd$a - 0.5 * nd$b
  nd
}

.rl_parts <- function(d, vars = c("a", "b")) {
  list(X = cbind(Intercept = 1, as.matrix(sf::st_drop_geometry(d)[, vars, drop = FALSE])),
       y = d$z, xy = sf::st_coordinates(d))
}

# The leave-one-location-out residuals by their definition: each location's
# rows predicted by the local regression there on the rows of the other
# locations, the kernel keeping the width it has in the fit.
.rl_brute <- function(X, y, xy, bw, kernel, adaptive) {
  loc  <- .gwr_locations(xy)
  full <- .gwr_local_solve(X, y, xy, bw, kernel, adaptive, what = "score")
  out  <- rep(NA_real_, length(y))
  for (g in seq_len(loc$n)) {
    r <- which(loc$id == g)
    if (!is.finite(full$width[r[1L]]) || full$width[r[1L]] <= 0) next
    pr <- .gwr_local_solve(X[-r, , drop = FALSE], y[-r], xy[-r, , drop = FALSE],
                           full$width[r[1L]], kernel, FALSE,
                           rp = xy[r, , drop = FALSE], X_rp = X[r, , drop = FALSE],
                           what = "predict")$fitted
    out[r] <- y[r] - pr
  }
  out
}


test_that("locations are told apart exactly, and the distinct rows are counted", {
  xy <- cbind(c(1, 2, 1, 3, 2, 1), c(5, 5, 5, 5, 5, 6))
  loc <- .gwr_locations(xy)
  expect_identical(loc$id, c(1L, 2L, 1L, 3L, 2L, 4L))
  expect_identical(loc$n, 4L)
  expect_identical(loc$size, c(2L, 2L, 1L, 1L))
  # One unit in the last place is another location; printing both to 15
  # digits gives the same text.
  near <- cbind(c(0.1 + 0.2, 0.3), c(1, 1))
  expect_identical(format(near[1, 1], digits = 15), format(near[2, 1], digits = 15))
  expect_identical(.gwr_locations(near)$n, 2L)
  expect_identical(.gwr_locations(xy[1, , drop = FALSE])$n, 1L)
  # Rows that share a location and every predictor are one distinct row.
  X <- cbind(1, c(7, 8, 7, 9, 8.5, 7))
  expect_identical(.gwr_n_distinct(X, loc$id), 5L)
  expect_identical(.gwr_n_distinct(X, seq_len(6L)), 6L)
})

test_that("the default bandwidth of a layer with 3 observations per site gives coefficients and predictions", {
  d  <- .rl_sites(100, 3, 2)
  nd <- .rl_new_sites()
  r  <- .rl_catch(fit_gwr_model(d, "z", c("a", "b")))
  fit <- r$value
  expect_s3_class(fit, "gwr_fit")
  # Before: 7 neighbours, 900 of 900 coefficients NA, 40 of 40 predictions NA.
  # GWmodel's search gave 58 neighbours and a prediction error of 0.197.
  expect_gte(fit$info$bandwidth, 58)
  expect_identical(sum(is.na(as.matrix(coef(fit)))), 0L)
  expect_identical(fit$info$n_rank_deficient, 0L)
  p <- suppressWarnings(predict(fit, newdata = nd))
  expect_false(anyNA(p))
  expect_lt(sqrt(mean((p - nd$truth)^2)), 0.197)
  expect_lt(sqrt(mean((fitted(fit) - d$mu)^2)), 0.2)
  # One warning, classed, with the counts and how the bandwidth was chosen.
  expect_identical(r$classes, "spatialkit_gwr_repeated_locations")
  expect_match(r$warnings, "^fit_gwr_model\\(\\): the 300 observations are at 100 distinct location\\(s\\) \\(100 hold more than one observation, the most 3\\)")
  expect_match(r$warnings, "No bandwidth was supplied: it was selected by leave-one-location-out cross-validation \\(each location predicted from the others\\), not by AICc")
  expect_identical(fit$info$n_locations, 100L)
  # What the search did is on the fit.
  bs <- fit$info$bandwidth_search
  expect_identical(bs$approach, "CV")
  expect_identical(bs$requested, "AICc")
  expect_true(bs$by_location)
  expect_identical(bs$n_locations, 100L)
  # Four sites of 3 rows and the row the bisquare kernel gives weight 0.
  expect_identical(bs$range, c(13, 300))
  expect_identical(bs$curve$bandwidth[which.min(bs$curve$criterion)], fit$info$bandwidth)
  # The summary says it too.
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("Locations: 100 for 300 observations", out, fixed = TRUE)))
})

test_that("AICc falls towards the bandwidths that interpolate the site means, and they cannot be selected", {
  d <- .rl_sites(100, 3, 2); p <- .rl_parts(d)
  cr <- lapply(c(5, 7, 10, 13, 58), function(k)
    .gwr_criterion(p$X, p$y, p$xy, k, "bisquare", TRUE, "AICc"))
  # The numbers the help page quotes: 58 at 5 to 10 neighbours, 90 at 13,
  # 221 at 58.
  expect_equal(vapply(cr, `[[`, numeric(1), "AICc"),
               c(58.43, 58.43, 58.43, 90.02, 220.59), tolerance = 1e-4)
  # At 5 to 10 neighbours tr(S) is the number of sites: the fit is the site
  # means.  AICc is not a criterion there, and the reason is given.
  expect_equal(vapply(cr[1:3], `[[`, numeric(1), "tr_S"), rep(100, 3), tolerance = 1e-8)
  expect_identical(vapply(cr, `[[`, numeric(1), "n_distinct"), rep(100, 5))
  expect_identical(vapply(cr, `[[`, character(1), "why"),
                   c("df_distinct", "df_distinct", "df_distinct", "", ""))
  expect_identical(vapply(cr[1:3], `[[`, numeric(1), "value"), rep(Inf, 3))
  # The leave-one-location-out score has no value where a site cannot be
  # predicted from the others, and is far higher at 16 than at 58.
  cv <- lapply(c(7, 10, 16, 58, 150), function(k)
    .gwr_criterion(p$X, p$y, p$xy, k, "bisquare", TRUE, "CV"))
  v <- vapply(cv, `[[`, numeric(1), "value")
  expect_identical(v[1:2], c(Inf, Inf))
  expect_identical(vapply(cv, `[[`, character(1), "why"), c("leave_out", "leave_out", "", "", ""))
  expect_gt(v[3], 5 * v[4])
  expect_lt(v[5], v[4])
  # At 13, the bottom of the range, one site's window holds a fourth site at
  # weight 1.4e-5: without the site itself two sites are left for three
  # parameters, to the solve's tolerance, and the score is undefined there.
  cv13 <- .gwr_criterion(p$X, p$y, p$xy, 13, "bisquare", TRUE, "CV")
  expect_identical(cv13$value, Inf)
  expect_identical(cv13$why, "leave_out")
  # The error of prediction at new sites, the other numbers of the help page:
  # 46 at 10 neighbours, 1.35 at 13, 0.15 at 150.
  set.seed(99); m <- 200
  nd <- sf::st_as_sf(data.frame(x = 5e5 + runif(m, 100, 900), y = 5e6 + runif(m, 100, 900),
                                a = rnorm(m), b = rnorm(m)), coords = c("x", "y"), crs = 32632)
  truth <- 1 + (1 + (sf::st_coordinates(nd)[, 1] - 5e5) / 1000) * nd$a - 0.5 * nd$b
  err <- vapply(c(10, 13, 150), function(k) {
    f <- suppressWarnings(suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = k)))
    sqrt(mean((suppressWarnings(predict(f, newdata = nd)) - truth)^2))
  }, numeric(1))
  expect_equal(err, c(46.1, 1.35, 0.149), tolerance = 0.01)
})

test_that("the leave-one-location-out residuals are those of refitting without the location", {
  d <- .rl_sites(60, 3, 5); p <- .rl_parts(d)
  # A second design in which one predictor differs between the rows of a site.
  set.seed(8)
  Xv <- p$X; Xv[, "b"] <- Xv[, "b"] + rnorm(nrow(Xv), 0, 0.5)
  loc <- .gwr_locations(p$xy)
  for (cs in list(list(X = p$X, bw = 13, k = "bisquare", ad = TRUE),
                  list(X = p$X, bw = 40, k = "bisquare", ad = TRUE),
                  list(X = Xv,  bw = 25, k = "tricube", ad = TRUE),
                  list(X = Xv,  bw = 30, k = "boxcar", ad = TRUE),
                  list(X = p$X, bw = 150, k = "gaussian", ad = FALSE),
                  list(X = Xv,  bw = 20, k = "exponential", ad = TRUE),
                  list(X = Xv,  bw = 260, k = "bisquare", ad = FALSE))) {
    lab <- paste(cs$k, cs$bw)
    sv  <- .gwr_local_solve(cs$X, p$y, p$xy, cs$bw, cs$k, cs$ad, what = "score", loc = loc)
    got <- .gwr_lolo_residuals(p$y - sv$fitted, sv$leverage, sv$hat_loc, loc)
    ref <- .rl_brute(cs$X, p$y, p$xy, cs$bw, cs$k, cs$ad)
    expect_false(anyNA(ref), label = lab)
    expect_equal(got, ref, tolerance = 1e-8, label = lab)
    expect_equal(.gwr_criterion(cs$X, p$y, p$xy, cs$bw, cs$k, cs$ad, "CV")$value,
                 sum(ref^2), tolerance = 1e-8, label = lab)
  }
  # Where the other sites of a window cannot place a site, both say so: at 10
  # neighbours a bisquare window holds 3 sites, and 2 are left for 3
  # parameters.
  sv  <- .gwr_local_solve(p$X, p$y, p$xy, 10, "bisquare", TRUE, what = "score", loc = loc)
  got <- .gwr_lolo_residuals(p$y - sv$fitted, sv$leverage, sv$hat_loc, loc)
  expect_true(all(is.na(got)))
  expect_true(all(is.na(.rl_brute(p$X, p$y, p$xy, 10, "bisquare", TRUE))))
  # A layer without repeat locations: the leave-one-out residual e / (1 - h).
  one <- !duplicated(loc$id)
  l1  <- .gwr_locations(p$xy[one, ])
  s1  <- .gwr_local_solve(p$X[one, ], p$y[one], p$xy[one, ], 20, "bisquare", TRUE, what = "score")
  expect_equal(.gwr_lolo_residuals(p$y[one] - s1$fitted, s1$leverage, NULL, l1),
               (p$y[one] - s1$fitted) / (1 - s1$leverage), tolerance = 1e-12)
})

test_that("the range of the search counts locations, for either kind of bandwidth", {
  d <- .rl_sites(100, 3, 2); p <- .rl_parts(d)
  # Adaptive: the rows of the four nearest sites (itself included), and one
  # more for the kernels that give the farthest row weight 0.
  expect_identical(.gwr_bw_range(p$xy, 3L, "bisquare", TRUE), c(13, 300))
  expect_identical(.gwr_bw_range(p$xy, 3L, "tricube", TRUE), c(13, 300))
  expect_identical(.gwr_bw_range(p$xy, 3L, "boxcar", TRUE), c(10, 300))
  expect_identical(.gwr_bw_range(p$xy, 3L, "gaussian", TRUE), c(10, 300))
  # At that count every bisquare window holds four sites.
  at13 <- vapply(seq_len(300L), function(i) {
    dd <- sqrt((p$xy[, 1] - p$xy[i, 1])^2 + (p$xy[, 2] - p$xy[i, 2])^2)
    length(unique(d$site[.gw_kernel_weights(dd, 13, "bisquare", TRUE) > 0]))
  }, integer(1))
  expect_true(all(at13 >= 4L))
  # Sites of unequal size: the count is the largest any window needs.
  xy <- rbind(cbind(rep(0, 5), 0), cbind(c(1, 2, 3, 10, 11, 12, 13), 0))
  expect_identical(.gwr_bw_range(xy, 2L, "bisquare", TRUE), c(8, 12))
  # Fixed: the median, over the rows, of the distance to the site that gives
  # a window four of them.
  one <- !duplicated(d$site)
  dm  <- as.matrix(stats::dist(p$xy[one, ]))
  d5  <- apply(dm, 1L, function(v) sort(v)[5L])
  expect_equal(.gwr_bw_range(p$xy, 3L, "bisquare", FALSE), c(stats::median(d5), max(dm)),
               tolerance = 1e-12)
  # Counted in rows it would be the distance to the next site only.
  expect_gt(.gwr_bw_range(p$xy, 3L, "bisquare", FALSE)[1L],
            2 * stats::median(apply(dm, 1L, function(v) sort(v)[2L])))
  # Without repeat locations nothing changes.
  expect_identical(.gwr_bw_range(p$xy[one, ], 3L, "bisquare", TRUE), c(5, 100))
  expect_identical(.gwr_bw_range(p$xy[one, ], 3L, "gaussian", TRUE), c(4, 100))
})

test_that("a fixed bandwidth, two observations per site and a kernel without an edge are selected the same way", {
  nd <- .rl_new_sites()
  # GWmodel's search before: prediction error 0.941 (fixed, 100 x 3), 0.145
  # and 0.151 (150 x 2).  The first draft of the new search: all NA, 0.234 and
  # 0.229.
  for (cs in list(list(sites = 100, reps = 3, seed = 2, ad = FALSE, k = "bisquare", below = 0.2),
                  list(sites = 150, reps = 2, seed = 3, ad = TRUE,  k = "bisquare", below = 0.145),
                  list(sites = 150, reps = 2, seed = 3, ad = FALSE, k = "bisquare", below = 0.151),
                  list(sites = 60,  reps = 2, seed = 1, ad = TRUE,  k = "gaussian", below = 0.2))) {
    lab <- sprintf("%d x %d, %s %s", cs$sites, cs$reps, if (cs$ad) "adaptive" else "fixed", cs$k)
    d <- .rl_sites(cs$sites, cs$reps, cs$seed)
    r <- .rl_catch(fit_gwr_model(d, "z", c("a", "b"), adaptive = cs$ad, kernel = cs$k))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    expect_identical(r$classes, "spatialkit_gwr_repeated_locations", label = lab)
    expect_true(fit$info$bandwidth_search$by_location, label = lab)
    expect_identical(sum(is.na(as.matrix(coef(fit)))), 0L, label = lab)
    p <- suppressWarnings(predict(fit, newdata = nd))
    expect_false(anyNA(p), label = lab)
    expect_lt(sqrt(mean((p - nd$truth)^2)), cs$below, label = lab)
  }
})

test_that("a supplied bandwidth that interpolates the site means has no AICc, and the fit says why", {
  d <- .rl_sites(100, 3, 2)
  r <- .rl_catch(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 10))
  fit <- r$value
  expect_s3_class(fit, "gwr_fit")
  expect_equal(fit$info$tr_S, 100, tolerance = 1e-8)
  # The formula gives 58.43, which would rank this fit above the one at 58
  # neighbours (220.59).
  expect_true(is.na(fit$info$AICc))
  w <- grep("AICc is reported as NA", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "The 300 observations are 100 distinct combinations of location and predictor values")
  expect_match(w, "tr\\(S\\) = 100, is not below that number minus 2 \\(98\\)")
  expect_match(w, "the value the formula gives \\(58\\.4334\\)")
  expect_match(w, "Use a larger bandwidth\\.$")
  # The warning about the shared locations does not claim a selection.
  rep_w <- r$warnings[r$classes == "spatialkit_gwr_repeated_locations"]
  expect_length(rep_w, 1L)
  expect_no_match(rep_w, "No bandwidth was supplied")
  expect_match(rep_w, "If the repeats are replicates of one measurement, fit their means, one row per location\\.$")
  # At a bandwidth that leaves the windows room, AICc is reported.
  f58 <- suppressWarnings(suppressMessages(fit_gwr_model(d, "z", c("a", "b"), bandwidth = 58)))
  expect_equal(f58$info$AICc, 220.59, tolerance = 1e-4)
  # A predictor that differs between the rows of a site: every row is its own
  # distinct row, and the rule is the usual one (n - 2).
  set.seed(4); d$w <- d$a + rnorm(300, 0, 0.5)
  p <- .rl_parts(d, c("w", "b"))
  cr <- .gwr_criterion(p$X, p$y, p$xy, 13, "bisquare", TRUE, "AICc")
  expect_identical(cr$n_distinct, 300L)
  expect_identical(cr$why, "")
})

test_that("a fold of a cross-validation logs the shared locations and does not warn", {
  d <- .rl_sites(60, 2, 1)
  dat <- prep_model_data(d, "z", c("a", "b"))
  r <- .rl_catch(fit_gwr_model(dat, "z", c("a", "b"), .already_prepped = TRUE))
  expect_s3_class(r$value, "gwr_fit")
  expect_length(r$warnings, 0L)
  expect_true(r$value$info$bandwidth_search$by_location)
})

test_that("a layer without repeat locations is searched by AICc and raises nothing", {
  d <- .rl_sites(100, 1, 2)
  r <- .rl_catch(fit_gwr_model(d, "z", c("a", "b")))
  fit <- r$value
  expect_length(r$warnings, 0L)
  expect_identical(fit$info$n_locations, 100L)
  bs <- fit$info$bandwidth_search
  expect_identical(bs$approach, "AICc")
  expect_false(bs$by_location)
  expect_identical(bs$range, c(5, 100))
  expect_false(any(grepl("Locations:", capture.output(print(summary(fit))), fixed = TRUE)))
  # One pair of rows at one location is enough to change the criterion.
  d2 <- rbind(d, d[1, ]); d2$z[101] <- d2$z[101] + 0.1
  r2 <- .rl_catch(fit_gwr_model(d2, "z", c("a", "b")))
  expect_identical(r2$classes, "spatialkit_gwr_repeated_locations")
  expect_match(r2$warnings, "the 101 observations are at 100 distinct location\\(s\\) \\(1 hold more than one observation, the most 2\\)")
  expect_identical(r2$value$info$bandwidth_search$approach, "CV")
})

test_that("a bandwidth at which a term is estimable nowhere cannot be selected", {
  # A covariate that takes one value west of x = 500 and another east of it.
  # Bisquare windows of 40 m hold one side only, wherever they are.
  set.seed(3); n <- 200
  x <- c(runif(n / 2, 0, 400), runif(n / 2, 600, 1000)); y <- runif(n, 0, 1000); a <- rnorm(n)
  region <- as.numeric(x > 500)
  z <- 1 + 2 * a + 3 * region + rnorm(n, 0, 0.3)
  X <- cbind(Intercept = 1, a = a, region = region); xy <- cbind(x, y)
  sv <- .gwr_local_solve(X, z, xy, 150, "bisquare", FALSE, what = "score")
  # The slope of `a` is estimable in every window; the covariate's and, with
  # it, the intercept's in none on the east side, the intercept's on the west
  # (where the covariate is 0).
  expect_identical(sv$n_estimable[2:3], c(200L, 0L))
  expect_identical(sv$n_estimable[1L], 100L)
  cr <- .gwr_criterion(X, z, xy, 150, "bisquare", FALSE, "AICc")
  expect_identical(cr$value, Inf)
  expect_identical(cr$why, "term")
  expect_true(is.finite(cr$AICc))
  # A bandwidth that reaches across the gap estimates it in some windows.
  wide <- .gwr_criterion(X, z, xy, 400, "bisquare", FALSE, "AICc")
  expect_identical(wide$why, "")
  expect_true(is.finite(wide$value))
  expect_gt(.gwr_local_solve(X, z, xy, 400, "bisquare", FALSE, what = "score")$n_estimable[3L], 0L)
  # The fit counts the same windows, term by term, as the score does.
  d <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, a = a, region = region, z = z),
                    coords = c("x", "y"), crs = 32632)
  f <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("a", "region"), adaptive = FALSE, bandwidth = 400)))
  expect_identical(as.integer(colSums(!f$info$nonfinite_coef)),
                   .gwr_local_solve(X, z, xy, 400, "bisquare", FALSE, what = "score")$n_estimable)
  # The search does not settle where the covariate has no coefficient at all.
  fs <- suppressWarnings(suppressMessages(fit_gwr_model(d, "z", c("a", "region"), adaptive = FALSE)))
  expect_gt(sum(!is.na(coef(fs)$region)), 0L)
  # A predictor that is one number over the whole layer has no coefficient at
  # any bandwidth: none can be selected, and the warning gives that reason.
  d$flat <- 2
  r <- .rl_catch(fit_gwr_model(d, "z", c("a", "flat")))
  expect_s3_class(r$value, "gwr_fit")
  expect_true(r$value$info$bandwidth_search$none_defined)
  expect_identical(names(r$value$info$bandwidth_search$undefined), "term")
  w <- grep("has a defined AICc", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "so none could be selected: at each of the [0-9]+ tried a term can be estimated at no location\\. The widest, 200, is used")
})

test_that("gwr_model_selection() takes its bandwidth and its warnings from the same rules", {
  d <- .rl_sites(100, 3, 2)
  r <- .rl_catch(gwr_model_selection(d, "z", c("a", "b")))
  sel <- r$value
  expect_s3_class(sel, "gwr_model_selection")
  # Before: 7 neighbours, the three models tied at 58.43343, `a` selected.
  fit <- suppressWarnings(suppressMessages(fit_gwr_model(d, "z", c("a", "b"))))
  expect_identical(as.numeric(sel$bandwidth), fit$info$bandwidth)
  expect_identical(sel$best, c("a", "b"))
  expect_identical(sel$bandwidth_source, "leave-one-location-out CV search on the full model")
  expect_true(sel$bandwidth_search$by_location)
  expect_gt(sel$table$delta[2L], 100)
  expect_identical(r$classes, "spatialkit_gwr_repeated_locations")
  expect_match(r$warnings, "^gwr_model_selection\\(\\): the 300 observations are at 100 distinct location\\(s\\) \\(100 hold more than one observation, the most 3\\)")
  expect_match(r$warnings, "No bandwidth was supplied: it was selected by leave-one-location-out cross-validation")
  # Asked for by CV, the search leaves out locations too: with a twin left
  # in, the leave-one-out score is lowest where the windows hold one site.
  cvs <- suppressWarnings(suppressMessages(gwr_model_selection(d, "z", c("a", "b"), bw_approach = "CV")))
  expect_identical(as.numeric(cvs$bandwidth), fit$info$bandwidth)
  expect_identical(cvs$bandwidth_search$requested, "CV")
  # A supplied bandwidth at which every model reproduces the site means: the
  # three are not ranked, and the warning names them.
  r7 <- .rl_catch(gwr_model_selection(d, "z", c("a", "b"), bandwidth = 7))
  t7 <- r7$value$table
  expect_true(all(is.na(t7$criterion[t7$variables != "<none>"])))
  expect_identical(t7$variables[1L], "<none>")
  w <- grep("are not ranked", r7$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "^gwr_model_selection\\(\\): 3 of 4 model\\(s\\) are not ranked at bandwidth 7; their criterion is NA and they are listed last\\.")
  expect_match(w, "3 \\(a; b; a \\+ b\\) reproduce the mean of the observations at each location, or nearly")
  expect_match(w, "Use a larger bandwidth\\.$")
  rep7 <- r7$warnings[r7$classes == "spatialkit_gwr_repeated_locations"]
  expect_length(rep7, 1L)
  expect_no_match(rep7, "No bandwidth was supplied")
})

test_that("cv_gwr(bandwidth = NULL) fits its folds at bandwidths that predict", {
  skip_on_cran()
  # Before, the folds of this layer were fitted at 7 to 10 neighbours.
  d <- .rl_sites(100, 3, 2)
  cv <- suppressWarnings(suppressMessages(
    cv_gwr(d, "z", c("a", "b"), folds = make_folds(d, k = 5, method = "random_kfold", seed = 1))))
  expect_identical(sum(cv$fold_status$status == "ok"), 5L)
  expect_true(all(cv$fold_metrics$bandwidth >= 58))
  expect_identical(cv$overall$n_pred, 300L)
  # Noise SD 0.2 within a site and 0.3 between sites: 0.36 together.
  expect_lt(cv$overall$RMSE, 0.45)
})
