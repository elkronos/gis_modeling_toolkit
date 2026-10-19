# tests/testthat/test-fix-gwr-bandwidth-search.R
# ---------------------------------------------------------------------------
# The default bandwidth of a GWR is searched for on a grid first.
#
# GWmodel::bw.gwr() runs a golden-section search over the whole range, which
# finds the minimum of a criterion that has one.  When the coefficients vary
# at a range only a few dozen observations resolve, AICc has a minimum at a
# small bandwidth and falls again towards a global fit; the search compares
# its first two probes (at 0.38 and 0.62 of the range), drops the lower part
# of the range, and returns a bandwidth of nearly every observation as the
# selected one.  It also scores any bandwidth with a window of deficient rank
# as infinite or as rounding error.
#
# The package's search scores 16 log-spaced bandwidths over the whole range,
# refines round the best of them, and keeps what it scored.
# ---------------------------------------------------------------------------

skip_if_not_installed("sp")

.bs_catch <- function(expr) {
  w <- character(0); cls <- character(0)
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); cls <<- c(cls, class(cnd)[1L])
      invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, warnings = w, classes = cls)
}

# 300 points; the slope of `a` varies with a period of about 500 m.
.bs_local_slope <- function(n, seed) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n); b <- rnorm(n)
  mu <- 2 + (1 + sin(x / 80) * cos(y / 80)) * a + 0.5 * b
  sf::st_as_sf(data.frame(x = x, y = y, a = a, b = b, z = mu + rnorm(n, 0, 0.5), mu = mu),
               coords = c("x", "y"), crs = 32632)
}

# The design, response and coordinates fit_gwr_model() works on.
.bs_parts <- function(d, vars) {
  list(X = cbind(Intercept = 1, as.matrix(sf::st_drop_geometry(d)[, vars, drop = FALSE])),
       y = d$z, xy = sf::st_coordinates(d))
}


test_that("the default bandwidth is the criterion's minimum, not a near-global one", {
  skip_if_not_installed("GWmodel")
  skip_on_cran()
  # Before: this layer got 299 of 300 neighbours, 52 AICc above the minimum
  # over the counts from 5 up, which lies at 36.
  d <- .bs_local_slope(300, 2)
  r <- .bs_catch(fit_gwr_model(d, "z", c("a", "b")))
  fit <- r$value
  expect_s3_class(fit, "gwr_fit")
  expect_length(r$warnings, 0L)
  expect_lt(fit$info$bandwidth, 0.2 * 300)
  p <- .bs_parts(d, c("a", "b"))
  # Every count to 100 and every fifth beyond.
  counts <- c(5:100, seq(105, 300, by = 5))
  full <- vapply(counts, function(k)
    .gwr_criterion(p$X, p$y, p$xy, k, "bisquare", TRUE)$value, numeric(1))
  expect_identical(counts[which.min(full)], 36)
  expect_lt(fit$info$AICc - min(full), 1)
  expect_lt(abs(fit$info$bandwidth - 36), 10)
  expect_false(fit$info$bandwidth_is_fallback)
  # GWmodel's search on the same data, scored the same way.
  gw <- .gwr_quietly(suppressWarnings(GWmodel::bw.gwr(
    z ~ a + b, data = .to_sp(d, c("z", "a", "b")), approach = "AICc",
    kernel = "bisquare", adaptive = TRUE)))
  expect_gt(gw, 0.9 * 300)
  expect_gt(.gwr_criterion(p$X, p$y, p$xy, gw, "bisquare", TRUE)$value - fit$info$AICc, 20)
  # What the search scored is kept, and the bandwidth is its minimum.
  bs <- fit$info$bandwidth_search
  expect_identical(bs$approach, "AICc")
  expect_identical(bs$range, c(5, 300))
  expect_named(bs$curve, c("bandwidth", "criterion", "tr_S", "n_deficient", "stage"))
  expect_identical(sum(bs$curve$stage == "grid"), 16L)
  expect_lte(bs$n_evaluations, 30L)
  expect_identical(nrow(bs$curve), bs$n_evaluations)
  expect_false(is.unsorted(bs$curve$bandwidth))
  expect_identical(bs$curve$bandwidth[which.min(bs$curve$criterion)], fit$info$bandwidth)
  expect_equal(bs$criterion, fit$info$AICc, tolerance = 1e-8)
  # Every scored value is the criterion at that bandwidth.
  on_grid <- bs$curve$bandwidth %in% counts
  expect_equal(bs$curve$criterion[on_grid],
               full[match(bs$curve$bandwidth[on_grid], counts)], tolerance = 1e-10)
  # Two minima on the grid: the one selected and the global fit's.
  expect_true(bs$multimodal)
  expect_gte(bs$n_local_minima, 2L)
  expect_false(bs$near_global)
  expect_false(bs$at_upper); expect_false(bs$at_lower); expect_false(bs$none_defined)
  expect_match(.gwr_bw_search_note(bs, TRUE), "the criterion has 2 local minima on the 16-point grid")
  # Another draw of the same process: 299 neighbours before, 35 now.
  f3 <- suppressMessages(fit_gwr_model(.bs_local_slope(300, 3), "z", c("a", "b")))
  expect_lt(f3$info$bandwidth, 0.2 * 300)
})

test_that("on a stationary process the search ends near the top of the range, and info says so", {
  skip_if_not_installed("GWmodel")
  set.seed(201); n <- 300
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 2 + d$a + 0.5 * d$b + rnorm(n, 0, 0.5)
  r <- .bs_catch(fit_gwr_model(d, "z", c("a", "b")))
  fit <- r$value
  # A minimum at or near the top is the ordinary outcome here: no warning.
  expect_length(r$warnings, 0L)
  bs <- fit$info$bandwidth_search
  expect_gt(fit$info$bandwidth, 0.8 * n)
  # Every observation in each window scores within 2 of the minimum: AICc
  # does not tell this fit from a global regression, and info records it.
  expect_true(bs$near_global)
  expect_identical(bs$within[2L], 300)
  expect_identical(bs$n_local_minima, 1L)
  expect_false(bs$multimodal)
  # No worse than GWmodel's search, which is right on such data.
  gw <- .gwr_quietly(suppressWarnings(GWmodel::bw.gwr(
    z ~ a + b, data = .to_sp(d, c("z", "a", "b")), approach = "AICc",
    kernel = "bisquare", adaptive = TRUE)))
  p <- .bs_parts(d, c("a", "b"))
  expect_lte(fit$info$AICc, .gwr_criterion(p$X, p$y, p$xy, gw, "bisquare", TRUE)$value + 0.2)
  # The one sentence the log line carries.
  note <- .gwr_bw_search_note(bs, TRUE)
  expect_match(note, "the widest bandwidth of the range \\(300 neighbours\\) scores within 2 of the minimum, so the criterion does not tell this fit from one global regression")
  expect_identical(.gwr_bw_search_note(NULL, TRUE), "")
})

test_that("the note names a floor minimum, several minima and a flat criterion", {
  bs <- list(none_defined = FALSE, bandwidth = 5, criterion = 100, approach = "AICc",
             range = c(5, 300), near_global = FALSE,
             at_upper = FALSE, at_lower = TRUE, multimodal = TRUE, n_local_minima = 3L,
             flat = TRUE, within = c(5, 90),
             curve = data.frame(stage = rep(c("grid", "refine"), c(16, 4))))
  note <- .gwr_bw_search_note(bs, TRUE)
  top <- utils::modifyList(bs, list(bandwidth = 300, at_upper = TRUE, at_lower = FALSE,
                                    near_global = TRUE, multimodal = FALSE, flat = FALSE))
  expect_identical(.gwr_bw_search_note(top, TRUE),
                   "the minimum is at the top of the range (300 neighbours), so the fit is close to one global regression.")
  expect_match(note, "the minimum is at the bottom of the range \\(5 neighbours\\)")
  expect_match(note, "the criterion has 3 local minima on the 16-point grid")
  expect_match(note, "the criterion is flat: bandwidths from 5 to 90 neighbours are within 2 of the minimum")
  bs$approach <- "CV"
  expect_match(.gwr_bw_search_note(bs, FALSE), "bandwidths from 5 to 90 are within 1 percent of the minimum")
  bs$none_defined <- TRUE
  expect_identical(.gwr_bw_search_note(bs, TRUE), "")
})

test_that("a fixed bandwidth is searched from where half the windows are complete, and beats bw.gwr()", {
  skip_if_not_installed("GWmodel")
  d <- .bs_local_slope(300, 2)
  p <- .bs_parts(d, c("a", "b"))
  sp_d <- .to_sp(d, c("z", "a", "b"))
  for (kernel in c("bisquare", "gaussian")) {
    fit <- suppressWarnings(suppressMessages(
      fit_gwr_model(d, "z", c("a", "b"), adaptive = FALSE, kernel = kernel)))
    bs <- fit$info$bandwidth_search
    # The range: from the median distance to the 5th-nearest point (a third of
    # it for the Gaussian) to the largest distance between two points.
    dm <- as.matrix(stats::dist(p$xy))
    d5 <- stats::median(apply(dm, 1L, function(v) sort(v)[5L]))
    expect_equal(bs$range, c(if (kernel == "gaussian") d5 / 3 else d5, max(dm)),
                 tolerance = 1e-10, info = kernel)
    expect_gt(fit$info$bandwidth, bs$range[1L])
    expect_lt(fit$info$bandwidth, bs$range[2L])
    gw <- .gwr_quietly(suppressWarnings(GWmodel::bw.gwr(
      z ~ a + b, data = sp_d, approach = "AICc", kernel = kernel, adaptive = FALSE)))
    a_gw <- .gwr_criterion(p$X, p$y, p$xy, gw, kernel, FALSE)$value
    expect_lte(fit$info$AICc, a_gw + 0.01, label = kernel)
  }
  # With the bisquare kernel GWmodel's search ends at 1048 m, 56 AICc above
  # this one's 215 m.
  f_b <- suppressWarnings(suppressMessages(
    fit_gwr_model(d, "z", c("a", "b"), adaptive = FALSE)))
  gw_b <- .gwr_quietly(suppressWarnings(GWmodel::bw.gwr(
    z ~ a + b, data = sp_d, approach = "AICc", kernel = "bisquare", adaptive = FALSE)))
  expect_gt(.gwr_criterion(p$X, p$y, p$xy, gw_b, "bisquare", FALSE)$value - f_b$info$AICc, 20)
})

test_that("a bandwidth with windows of deficient rank is scored on its fitted values, so a small one can win", {
  # A slope that varies at short range and a regional covariate.  GWmodel's
  # scorer gives every bandwidth with a deficient window Inf (coded 0/1) or
  # 950 to 3958 (coded 0.7/1.3, at 12 to 130 neighbours), so the search
  # settled at about the smallest bandwidth with none: 163 and 154 neighbours
  # of 400, fitted values 0.69 off the noise-free mean.  Scored on what the
  # windows estimate, the minimum is at 22 for both codings.
  set.seed(11); n <- 400
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000); a <- rnorm(n)
  mu <- 1 + (2 + 1.5 * sin(x / 90) * cos(y / 90)) * a + 3 * (x > 500)
  z <- mu + rnorm(n, 0, 0.3)
  got <- lapply(list(c(0, 1), c(0.7, 1.3)), function(lev) {
    d <- sf::st_as_sf(data.frame(x = x, y = y, a = a, z = z,
                                 region = ifelse(x > 500, lev[2], lev[1])),
                      coords = c("x", "y"), crs = 32632)
    r <- .bs_catch(fit_gwr_model(d, "z", c("a", "region")))
    fit <- r$value
    expect_s3_class(fit, "gwr_fit")
    expect_identical(fit$info$bandwidth, 22)
    expect_lt(sqrt(mean((fitted(fit) - mu)^2)), 0.3)
    expect_gt(fit$info$bandwidth_search$curve$n_deficient[
      fit$info$bandwidth_search$curve$bandwidth == 22], 200L)
    # The warning about the deficient windows says the bandwidth was selected.
    w <- r$warnings[r$classes == "spatialkit_gwr_rank_deficient"]
    expect_length(w, 1L)
    expect_match(w, "No bandwidth was supplied: the search selected this one \\(22\\)")
    fit
  })
  expect_equal(got[[2L]]$info$AICc, got[[1L]]$info$AICc, tolerance = 1e-8)
  expect_equal(got[[2L]]$info$bandwidth_search$curve$criterion,
               got[[1L]]$info$bandwidth_search$curve$criterion, tolerance = 1e-8)
})

test_that("the criterion at a bandwidth is GWmodel's AICc and its leave-one-out score", {
  skip_if_not_installed("GWmodel")
  set.seed(9); n <- 120
  d <- sf::st_as_sf(data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                               a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + (1 + sf::st_coordinates(d)[, 1] / 1000) * d$a - d$b + rnorm(n, 0, 0.4)
  p <- .bs_parts(d, c("a", "b"))
  dm <- GWmodel::gw.dist(p$xy)
  for (cs in list(list(bw = 30L, ad = TRUE, k = "bisquare"), list(bw = 250, ad = FALSE, k = "gaussian"),
                  list(bw = 18L, ad = TRUE, k = "exponential"))) {
    aicc <- GWmodel:::gwr.aic(cs$bw, p$X, p$y, cs$k, cs$ad, p$xy, dMat = dm, verbose = FALSE)
    cv   <- GWmodel:::gwr.cv(cs$bw, p$X, p$y, cs$k, cs$ad, p$xy, dMat = dm, verbose = FALSE)
    ours <- .gwr_criterion(p$X, p$y, p$xy, cs$bw, cs$k, cs$ad, "AICc")
    expect_equal(ours$value, aicc, tolerance = 1e-9, info = cs$k)
    expect_equal(.gwr_criterion(p$X, p$y, p$xy, cs$bw, cs$k, cs$ad, "CV")$value, cv,
                 tolerance = 1e-8, info = cs$k)
    expect_identical(ours$n_deficient, 0L)
  }
  # A bandwidth at which every window holds its own point alone has no AICc
  # and no leave-one-out prediction: Inf, so no search settles there.
  expect_identical(.gwr_criterion(p$X, p$y, p$xy, 1e-6, "bisquare", FALSE, "AICc")$value, Inf)
  expect_identical(.gwr_criterion(p$X, p$y, p$xy, 1e-6, "bisquare", FALSE, "CV")$value, Inf)
})

test_that("when no bandwidth has a defined AICc the widest is used and the fit says it was not selected", {
  set.seed(4); n <- 4
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n), b = rnorm(n)),
                    coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n)
  r <- .bs_catch(fit_gwr_model(d, "z", c("a", "b")))
  fit <- r$value
  expect_s3_class(fit, "gwr_fit")
  expect_identical(fit$info$bandwidth, 4)
  expect_true(fit$info$bandwidth_search$none_defined)
  expect_true(fit$info$bandwidth_is_fallback)
  w <- grep("has a defined AICc", r$warnings, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "^fit_gwr_model\\(\\): no bandwidth from 4 to 4 neighbours has a defined AICc, so none could be selected")
  expect_match(w, "The widest, 4, is used")
  expect_false(any(grepl("arbitrary fallback", r$warnings)))
})

test_that("cv_gwr(bandwidth = NULL) selects each fold's bandwidth by the same search", {
  skip_on_cran()
  # Before, the five folds of this layer were fitted at 222 to 238 of 240
  # neighbours, five nearly global regressions.  They are fitted at 30 to 37.
  # (On a draw whose AICc has two minima of nearly equal depth the folds can
  # still fall on either; each fit's info$bandwidth_search says when.)
  d <- .bs_local_slope(300, 1)
  cv <- suppressWarnings(suppressMessages(
    cv_gwr(d, "z", c("a", "b"), folds = make_folds(d, k = 5, method = "random_kfold", seed = 1))))
  expect_identical(sum(cv$fold_status$status == "ok"), 5L)
  bw <- cv$fold_metrics$bandwidth
  expect_lt(max(bw) / min(bw), 3)
  expect_lt(max(bw), 0.25 * 240)
})

test_that("gwr_model_selection() takes its default bandwidth from the same search, for either criterion", {
  d <- .bs_local_slope(300, 2)
  sel <- suppressWarnings(suppressMessages(gwr_model_selection(d, "z", c("a", "b"))))
  expect_identical(sel$bandwidth_source, "AICc search on the full model")
  p <- .bs_parts(d, c("a", "b"))
  ref <- .gwr_bw_search(p$X, p$y, p$xy, "bisquare", TRUE, approach = "AICc")
  expect_identical(as.numeric(sel$bandwidth), ref$bandwidth)
  expect_identical(sel$bandwidth_search$curve, ref$curve)
  expect_lt(sel$bandwidth, 100L)
  cvs <- suppressWarnings(suppressMessages(
    gwr_model_selection(d, "z", c("a", "b"), bw_approach = "CV")))
  expect_identical(cvs$bandwidth_source, "CV search on the full model")
  expect_identical(cvs$bandwidth_search$approach, "CV")
  expect_null(suppressWarnings(suppressMessages(
    gwr_model_selection(d, "z", c("a", "b"), bandwidth = 40)))$bandwidth_search)
})
