# tests/testthat/test-plotting-fits.R
# ---------------------------------------------------------------------------
# plot.spatial_fit() and plot_folds().
#
# The package shipped print() and summary() for fitted models but no plot(),
# so the diagnostics most likely to reveal a problem -- residual structure, and
# whether blocks actually separate the folds -- had to be hand-rolled.
#
# ggplot objects are built lazily, so constructing one does not prove it will
# render.  ggplot2::ggplot_build() forces the computation, which is what these
# tests assert on.
# ---------------------------------------------------------------------------

test_that("residual and observed-predicted plots build", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points(), predictor_vars = "w")

  for (ty in c("residuals", "observed_predicted")) {
    p <- plot(fit, type = ty)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))   # forces the layer computation
  }
})

test_that("the residual plot carries one point per observation", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 90)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  b <- ggplot2::ggplot_build(plot(fit, type = "residuals"))
  expect_equal(nrow(b$data[[1]]), 90)
})

test_that("observed = fitted + residual is what gets plotted", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 60)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  b <- ggplot2::ggplot_build(plot(fit, type = "observed_predicted"))
  pt_layer <- b$data[[2]]                        # layer 1 is the abline
  expect_equal(sort(pt_layer$x), sort(pts$z), tolerance = 1e-8)
})

test_that("an unknown plot type is rejected", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points())
  # A bare expect_error() passes on ANY error, including a typo in the test
  # itself.  match.arg() gives a stable message that names the valid values,
  # so pin it -- and pin that the valid values are still the documented four.
  expect_error(plot(fit, type = "nonsense"), "'arg' should be one of")
  expect_error(plot(fit, type = "nonsense"), "residuals")
  expect_error(plot(fit, type = "nonsense"), "observed_predicted")
  expect_error(plot(fit, type = "nonsense"), "variogram")
  expect_error(plot(fit, type = "nonsense"), "coefficients")
})

test_that("a fit without geometry is rejected", {
  skip_if_not_installed("ggplot2")
  fit <- lm_spatial_fit(surf_test_points())
  fit$data_sf <- sf::st_drop_geometry(fit$data_sf)
  expect_error(plot(fit, type = "residuals"), "no training geometry")
})

test_that("the residual variogram plot builds", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # A field with real structure at the lags the variogram resolves (the
  # default binning is cutoff/15 ~ 40 units here), so the range is
  # identifiable.  An earlier a = 25 field was flat by the first lag bin and
  # fitted only by luck.
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
                      coords = c("x", "y"), crs = 3857)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")

  p <- plot(fit, type = "variogram")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  # The fitted model line and the range marker are on it.
  layers <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomLine" %in% layers)
  expect_true("GeomVline" %in% layers)

  # Residuals with NO spatial structure give a flat, nugget-only variogram
  # that no model fits: every start is singular, or the optimiser never
  # converges.  That is the picture most worth seeing, and it used to be
  # refused with "could not be fitted; there may be too few finite
  # residuals": estimate_sac_range() now returns NA WITH the empirical
  # variogram attached whichever way the fit failed, and the plot draws the
  # points, with a model line only when a (non-converged) model exists.
  flat <- pts
  set.seed(6); flat$z <- rnorm(n)
  fit_flat <- lm_spatial_fit(flat, predictor_vars = "w")
  sac_flat <- suppressWarnings(estimate_sac_range(
    sf::st_sf(.resid = residuals(fit_flat), geometry = sf::st_geometry(flat)),
    ".resid"))
  expect_true(is.na(sac_flat))
  expect_s3_class(attr(sac_flat, "variogram"), "data.frame")
  # White noise is refused for want of spatial structure, by a test made
  # before any model is looked at.  It used to end as "no variogram model" or
  # "did not converge", whichever way gstat's fit happened to fail.
  expect_identical(attr(sac_flat, "rejected_reason"),
                   "no spatial structure detected at the nearest-neighbour scale")
  # response = FALSE: the layer assertions below are about the RESIDUAL
  # curve alone; the response overlay has its own test file
  # (test-plotting-diagnostics.R).
  p2 <- plot(fit_flat, type = "variogram", response = FALSE)
  expect_s3_class(p2, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p2))
  layers2 <- vapply(p2$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomPoint" %in% layers2)
  expect_false("GeomVline" %in% layers2)               # no range marker
  expect_equal("GeomLine" %in% layers2, !is.null(attr(sac_flat, "variogram_model")))
  expect_match(p2$labels$subtitle, "^No effective range")

  # The no-model outcome, made deterministic: the fit reports that no model
  # could be fitted.  The NA carries the empirical variogram and no model,
  # and the plot says why there is no range.  The stub replaces
  # .sac_fit_model(), the fit estimate_sac_range() makes (it used to replace
  # gstat::fit.variogram(), which is no longer called), and the layer is the
  # structured one, because white noise is refused before any fit is read.
  local_mocked_bindings(.sac_fit_model = function(...) NULL)
  sac_sing <- suppressWarnings(estimate_sac_range(
    sf::st_sf(.resid = residuals(fit), geometry = sf::st_geometry(pts)),
    ".resid"))
  expect_true(is.na(sac_sing))
  expect_s3_class(attr(sac_sing, "variogram"), "data.frame")
  expect_null(attr(sac_sing, "variogram_model"))
  expect_match(attr(sac_sing, "rejected_reason"), "^no variogram model could be fitted")
  p3 <- plot(fit, type = "variogram", response = FALSE)
  expect_no_error(ggplot2::ggplot_build(p3))
  layers3 <- vapply(p3$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomPoint" %in% layers3)
  expect_false("GeomLine" %in% layers3)
  expect_match(p3$labels$subtitle, "no variogram model could be fitted")
})


test_that("plot_folds maps every point to a fold and draws the block design", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 100)
  f <- make_folds(pts, k = 4, method = "block_kfold", seed = 1)

  p <- plot_folds(f, pts)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Two sf layers: the block outlines first, one row per block, then the
  # points, one row each.
  expect_length(b$data, 2L)
  expect_equal(nrow(b$data[[1]]), nrow(f$params$blocks))
  expect_equal(nrow(b$data[[1]]), f$params$blocks_used)
  expect_equal(nrow(b$data[[2]]), 100)
  # colour is discrete by fold, so no more distinct colours than folds
  expect_lte(length(unique(b$data[[2]]$colour)), length(f$folds))
  # blocks = FALSE leaves the points alone; folds without a design draw none.
  p0 <- plot_folds(f, pts, blocks = FALSE)
  expect_length(ggplot2::ggplot_build(p0)$data, 1L)
  fr <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
  expect_length(ggplot2::ggplot_build(plot_folds(fr, pts))$data, 1L)
})

test_that("plot_folds rejects mismatched inputs", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 50)
  f <- make_folds(pts, k = 3, method = "block_kfold", seed = 1)

  expect_error(plot_folds(list(), pts), "returned by make_folds")
  expect_error(plot_folds(f, sf::st_drop_geometry(pts)), "must be an sf object")

  # Folds built from a different layer share no row ids.
  other <- surf_test_points(n = 50, seed = 99)
  other$..row_id <- seq_len(50) + 10000L
  expect_error(plot_folds(f, other), "no points matched")
})

test_that("plot_folds accepts a boundary", {
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(n = 80)
  f <- make_folds(pts, k = 3, method = "block_kfold", seed = 1)
  bnd <- sf::st_as_sfc(sf::st_bbox(pts))
  p <- plot_folds(f, pts, boundary = bnd)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
})


# --- regressions found in the full-package audit -----------------------------

test_that("a sill-less variogram is drawn rather than refused", {
  # estimate_sac_range() returns NA when the fitted range exceeds the longest
  # lag, and used to discard the variogram with it -- so plot() refused in
  # exactly the case worth looking at.
  skip_if_not_installed("gstat")
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(200)
  fit <- lm_spatial_fit(pts, predictor_vars = character(0))
  # z carries a strong linear x/y trend, so an intercept-only fit leaves it
  # all in the residuals: the variogram rises monotonically and never reaches a
  # sill.  The deterministic attribute check lives in test-sac-range.R; this
  # asserts the user-visible outcome, that a plot comes back.
  # gstat cannot converge on a variogram that never reaches a sill -- that is
  # the whole point of this case, so its warning is expected here.
  p <- suppressWarnings(plot(fit, type = "variogram"))
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_gt(nrow(b$data[[1]]), 0L)          # the empirical points are there
})

test_that("all-NA residuals error instead of painting a uniformly grey map", {
  # limits = c(Inf, -Inf) builds without complaint and maps every value to
  # na.value, which looks like a result.
  skip_if_not_installed("ggplot2")
  pts <- surf_test_points(40)
  fit <- lm_spatial_fit(pts)
  broken <- fit
  registerS3method("residuals", "broken_fit",
                   function(object, ...) rep(NA_real_, nrow(object$data_sf)))
  class(broken) <- c("broken_fit", class(fit))
  for (ty in c("residuals", "observed_predicted", "variogram"))
    expect_error(plot(broken, type = ty), "no finite residuals")
})

# ---------------------------------------------------------------------------
# plot.sac_range(): the method estimate_sac_range()'s own messages used to
# point at without it existing ("plot(type = \"variogram\")" is
# plot.spatial_fit, which needs a fit).
# ---------------------------------------------------------------------------

test_that("plot.sac_range draws an identified range with model line and marker", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z), coords = c("x", "y"), crs = 3857)
  r <- estimate_sac_range(pts, "z")
  expect_true(is.finite(r))

  p <- plot(r)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_identical(p$labels$title, "Empirical variogram")
  expect_match(p$labels$subtitle, "^Effective range = ")
  layers <- vapply(p$layers, function(l) class(l$geom)[1L], character(1))
  expect_true("GeomLine" %in% layers)    # the fitted model
  expect_true("GeomVline" %in% layers)   # the range marker
})

test_that("plot.sac_range draws a rejected range and says why", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # The same white-noise residuals the residual-variogram test above pins as
  # rejected ("no variogram model" or "did not converge"): the range is NA
  # with the empirical variogram attached, and the picture is still worth
  # drawing.  Not every white-noise draw is rejected -- a lucky one fits a
  # finite range -- so the fixture is reused rather than re-rolled.
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  flat <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
                       coords = c("x", "y"), crs = 3857)
  set.seed(6); flat$z <- rnorm(n)
  fit_flat <- lm_spatial_fit(flat, predictor_vars = "w")
  r <- suppressWarnings(estimate_sac_range(
    sf::st_sf(.resid = residuals(fit_flat), geometry = sf::st_geometry(flat)),
    ".resid"))
  expect_true(is.na(r))
  expect_s3_class(attr(r, "variogram"), "data.frame")
  p <- plot(r)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_identical(p$labels$title, "Empirical variogram")
  expect_match(p$labels$subtitle, "^No effective range")
})

test_that("plot.sac_range refuses a bare NA and a non-sac_range object informatively", {
  skip_if_not_installed("ggplot2")
  bare <- structure(NA_real_, class = c("sac_range", "numeric"))
  expect_error(plot(bare), "carries no empirical variogram")
  expect_error(spatialkit:::plot.sac_range(1), "must be an object returned by")
})

test_that("plot.spatial_fit(type = 'variogram') and plot.sac_range share one drawer", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  set.seed(5); n <- 200
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  d <- as.matrix(stats::dist(cbind(x, y)))
  z <- as.numeric(t(chol(exp(-d / 50) + diag(1e-4, n))) %*% rnorm(n))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, w = rnorm(n)),
                      coords = c("x", "y"), crs = 3857)
  fit <- lm_spatial_fit(pts, predictor_vars = "w")
  p_fit <- plot(fit, type = "variogram")
  # Same residuals through the public estimator give the same subtitle.
  pts$.resid <- as.numeric(residuals(fit))
  p_sac <- plot(estimate_sac_range(pts, ".resid"))
  expect_identical(p_fit$labels$title, "Residual variogram")
  expect_identical(p_sac$labels$title, "Empirical variogram")
  expect_identical(p_fit$labels$subtitle, p_sac$labels$subtitle)
})

test_that("plot_folds' subtitle describes the scheme that was built", {
  skip_if_not_installed("ggplot2")
  # It used to be one fixed line about block size, printed under random folds
  # too, where there are no blocks.
  pts <- surf_test_points(n = 100)
  sub <- function(f) plot_folds(f, pts)$labels$subtitle
  fr <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
  expect_match(sub(fr), "Random folds")
  expect_false(grepl("[Bb]lock", sub(fr)))
  fb <- make_folds(pts, k = 4, method = "block_kfold", block_size = 300, seed = 1)
  # The blocks as built, and the block_size asked for, in the unit of the CRS
  # the folds were built in, not its identifier ("EPSG:3857 units", or a
  # whole WKT for a CRS with no EPSG code).
  expect_match(sub(fb), sprintf("^Blocks %s x %s m, block_size 300\n",
                                format(signif(fb$params$block_width, 3)),
                                format(signif(fb$params$block_height, 3))))
  expect_false(grepl(as.character(fb$params$crs), sub(fb), fixed = TRUE))
  # Two lines: one long one is clipped at the width these are drawn at, and
  # neither line is long enough to clip on its own.
  lines <- strsplit(sub(fb), "\n", fixed = TRUE)[[1]]
  expect_length(lines, 2L)
  expect_true(all(nchar(lines) <= 55L))
  fg <- make_folds(pts, k = 4, method = "block_kfold", block_nx = 3, block_ny = 3,
                   seed = 1)
  expect_match(sub(fg), "3 x 3 block grid")
  fl <- suppressWarnings(make_folds(pts, k = nrow(pts), method = "buffered_loo",
                                    buffer = 150))
  expect_match(sub(fl), "Leave-one-out with a 150 m buffer", fixed = TRUE)
  pts$g <- rep(c("a", "b", "c"), length.out = nrow(pts))
  fo <- make_folds(pts, k = 3, method = "leave_location_out", group_var = "g",
                   seed = 1)
  expect_match(sub(fo), "3 groups of `g`")
})

test_that("plot_folds' subtitle gives the blocks built, never in scientific notation", {
  skip_if_not_installed("ggplot2")
  # Blocks of 400,000 US survey feet were subtitled "Block size 4e+05 (US
  # survey foot)": the size asked for, in scientific notation, while the
  # blocks drawn were larger.
  set.seed(4)
  ft <- sf::st_as_sf(data.frame(x = runif(80, 1e6, 3e6), y = runif(80, 5e5, 1e6)),
                     coords = c("x", "y"), crs = 2264)
  fb <- suppressWarnings(make_folds(ft, k = 4, method = "block_kfold",
                                    block_size = 4e5, seed = 1))
  sub <- strsplit(plot_folds(fb, ft)$labels$subtitle, "\n", fixed = TRUE)[[1]][1L]
  expect_false(grepl("e+", sub, fixed = TRUE))
  expect_match(sub, ", block_size 400,000$")
  expect_match(sub, sprintf("^Blocks %s x %s US ft,",
                            format(signif(fb$params$block_width, 3), big.mark = ",",
                                   scientific = FALSE),
                            format(signif(fb$params$block_height, 3), big.mark = ",",
                                   scientific = FALSE)))
  expect_true(fb$params$block_width >= 4e5)
  # A grid set by its dimensions names the blocks it built as well.
  pts <- surf_test_points(n = 100)
  fg <- make_folds(pts, k = 4, method = "block_kfold", block_nx = 3, block_ny = 3,
                   seed = 1)
  expect_match(plot_folds(fg, pts)$labels$subtitle,
               "^3 x 3 block grid of [0-9.]+ x [0-9.]+ m blocks\n")
})

test_that("the variogram labels give the unit in short form, the caption's ranges too", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gstat")
  # GDAL's singular name read "Effective range = 463.9 metre" and
  # "Distance (metre)", and the caption's two ranges had no unit.
  set.seed(3); n <- 250
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  d <- as.matrix(dist(xy))
  xy$z <- as.numeric(t(chol(exp(-d / 150) + diag(0.5, n))) %*% rnorm(n))
  xy$elev <- rnorm(n); xy$price <- 2 * xy$elev + xy$z
  pts <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  r <- estimate_sac_range(pts, "z")
  skip_if_not(is.finite(r))
  p <- plot(r)
  expect_identical(p$labels$x, "Distance (m)")
  expect_match(p$labels$subtitle, "^Effective range = [0-9.]+ m$")
  fit <- fit_rf_model(pts, "price", "elev", num_trees = 50, seed = 1)
  cap <- gsub("\n", " ", plot(fit, type = "variogram")$labels$caption)
  if (grepl("effective ranges", cap, fixed = TRUE))
    expect_match(cap, "effective ranges [0-9.]+ m \\(residuals\\) and [0-9.]+ m \\(response\\)")
  # A layer with no CRS keeps "CRS units".
  expect_identical(.crs_unit_short(NA_character_), NA_character_)
  expect_identical(.crs_unit_short("EPSG:2264"), "US ft")
})

test_that("the sweep's axis names the CRS's unit and its table fits 80 columns", {
  skip_if_not_installed("ggplot2")
  # The axis read "Block edge length (EPSG:32632, log scale)", and the printed
  # table, with a method column and a long header, wrapped fold_max onto a
  # block of its own.
  sw <- data.frame(block_size = c(NA, 100, 200), method = c("random_kfold", "block_kfold",
                                                           "block_kfold"),
                   blocks_used = c(NA, 16L, 4L), k = 4L, n_folds_succeeded = 4L,
                   value = c(1.0123, 1.1234, 1.2345), fold_min = c(0.9, 1, 1.1),
                   fold_max = c(1.1333, 1.2444, 1.3555), fold_sd = 0.1,
                   block_width = c(NA, 120, 240), block_height = c(NA, 125, 250))
  class(sw) <- c("block_size_sweep", "data.frame")
  attr(sw, "metric") <- "RMSE"; attr(sw, "sac_range") <- 150
  attr(sw, "crs") <- "EPSG:32632"; attr(sw, "k") <- 4L; attr(sw, "n_fits") <- 12L
  expect_identical(plot(sw)$labels$x, "Smallest block edge asked for (m, log scale)")
  out <- utils::capture.output(print(sw))
  expect_true(all(nchar(out) <= 80L))
  head <- out[grepl("^ *block_size", out)]
  expect_length(head, 1L)
  expect_false(grepl("method", head, fixed = TRUE))
  expect_match(head, "folds_ok", fixed = TRUE)
  expect_match(head, "fold_max", fixed = TRUE)
})


test_that("print(<sac_range>) wraps the directional ranges within 80 columns", {
  # On one line, with a refused direction among them, the directional ranges
  # ran past 100 characters and scrolled sideways on a help page.
  r <- structure(463.88531, class = "sac_range",
                 directional = c(`0` = 248.16083, `45` = 413.73181,
                                 `90` = NA, `135` = 426.58801),
                 directional_status = c(`0` = "ok", `45` = "ok",
                                        `90` = "over_cutoff", `135` = "ok"),
                 directional_fitted = c(`0` = 248.16083, `45` = 413.73181,
                                        `90` = 2617.96761, `135` = 426.58801),
                 anisotropy = 1.719, crs = sf::st_crs(32632), detrended = FALSE)
  out <- utils::capture.output(print(r))
  expect_true(all(nchar(out) <= 80L))
  dir <- paste(trimws(out[seq(grep("directional:", out, fixed = TRUE),
                              grep("^  in ", out) - 1L)]), collapse = " ")
  expect_identical(dir, paste("directional: 0 deg = 248.1608, 45 deg = 413.7318,",
                              "90 deg = 2617.968 (past the fitted lags),",
                              "135 deg = 426.588  (ratio 1.72)"))
  # No anisotropy attribute: the line still prints, without a ratio.
  attr(r, "anisotropy") <- NULL
  expect_false(any(grepl("ratio", utils::capture.output(print(r)), fixed = TRUE)))
})

test_that("plot_folds draws its axes in the projected CRS of the map", {
  skip_if_not_installed("ggplot2")
  # The axes were a lon/lat graticule over UTM data, in degrees to the third
  # decimal, while the subtitle and block_size are in metres.
  set.seed(1)
  pts <- sf::st_as_sf(data.frame(x = 5e5 + runif(80, 0, 1000),
                                 y = 5e6 + runif(80, 0, 1000)),
                      coords = c("x", "y"), crs = 32632)
  f <- make_folds(pts, k = 5, method = "block_kfold", block_size = 300)
  expect_true(sf::st_crs(plot_folds(f, pts)$coordinates$datum) == sf::st_crs(32632))
  # Lon/lat points: the blocks, drawn first, are in the UTM CRS the folds
  # were built in, and so is the graticule.
  ll <- sf::st_transform(pts, 4326)
  fl <- make_folds(ll, k = 5, method = "block_kfold", block_size = 300)
  p <- plot_folds(fl, ll)
  expect_true(sf::st_crs(p$coordinates$datum) == sf::st_crs(fl$params$blocks))
  expect_silent(ggplot2::ggplot_build(p))
  # A map in lon/lat keeps the default lon/lat graticule.
  fr <- make_folds(ll, k = 5, method = "random_kfold", seed = 1)
  expect_true(sf::st_is_longlat(plot_folds(fr, ll)$coordinates$datum))
})
