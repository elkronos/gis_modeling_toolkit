# tests/testthat/test-fix-gwr-origin-tally-label.R
# ---------------------------------------------------------------------------
# Three repairs to what a GWR fit says about itself.
#
# 1. A predictor far from zero was still flagged as collinear (every window
#    counted in info$n_local_collinear, the coefficient map refused) once the
#    uncentred condition index passed 1e6, although the local solve centres
#    each window and returns the same slopes at any origin.
# 2. predict()'s two counted warnings carry their counts as fields, and a
#    caller that predicts in pieces can raise each once with totals.
# 3. print() headed the fit "<GWR (GWmodel)>", though no GWmodel code makes it.
# ---------------------------------------------------------------------------

skip_if_not_installed("sp")

.otl_catch <- function(expr) {
  conds <- list()
  val <- tryCatch(
    withCallingHandlers(suppressMessages(expr), warning = function(cnd) {
      conds[[length(conds) + 1L]] <<- cnd
      invokeRestart("muffleWarning")
    }),
    error = function(e) conditionMessage(e))
  list(value = val, conds = conds,
       warnings = vapply(conds, conditionMessage, character(1)),
       classes = vapply(conds, function(cnd) class(cnd)[1L], character(1)))
}

# 200 points, one predictor whose slope drifts eastward.
.otl_layer <- function() {
  set.seed(1); n <- 200
  x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4); t0 <- rnorm(n)
  d <- sf::st_as_sf(data.frame(x = x, y = y, t0 = t0), coords = c("x", "y"), crs = 32632)
  d$yield <- 50 - 0.8 * (1 + x / 1e4) * t0 + rnorm(n)
  d
}


test_that("a predictor shifted far from zero is not flagged, and its slope map is drawn", {
  d <- .otl_layer()
  ref <- .otl_catch(fit_gwr_model(d, "yield", "t0", bandwidth = 30))
  expect_length(ref$warnings, 0L)
  expect_identical(ref$value$info$n_local_collinear, 0L)
  for (off in c(1e5, 1e7, 1e9)) {
    lab <- sprintf("origin %g", off)
    d$tt <- d$t0 + off
    r <- .otl_catch(fit_gwr_model(d, "yield", "tt", bandwidth = 30))
    f <- r$value
    expect_s3_class(f, "gwr_fit")
    # Before, at 1e7 and 1e9: "100% of 200 locations have a collinear local
    # design", n_local_collinear 200, and the map stopped with "every
    # location is masked".
    expect_length(r$warnings, 0L)
    expect_identical(f$info$n_local_collinear, 0L, label = lab)
    # The slopes are the origin-0 fit's, to the precision the shifted column
    # has (its spacing is `off` times 2e-16).
    expect_lt(max(abs(coef(f)$tt - coef(ref$value)$t0)), off * 1e-15, label = lab)
    expect_lt(max(abs(fitted(f) - fitted(ref$value))), off * 1e-15, label = lab)
    # The index with the intercept is kept, for the Intercept map.
    expect_true(all(f$info$local_collinearity$cn > 1e4), label = lab)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- plot(f, type = "coefficients", term = "tt")
      expect_match(p$labels$subtitle, "No location masked", label = lab)
      pI <- plot(f, type = "coefficients", term = "Intercept", mask = FALSE)
      expect_s3_class(pI, "ggplot")
    }
  }
  # At 1e7 the index is past the 1e6 a normal-equation solve breaks at.
  d$tt <- d$t0 + 1e7
  f7 <- suppressMessages(fit_gwr_model(d, "yield", "tt", bandwidth = 30))
  expect_true(all(f7$info$local_collinearity$cn > 1e6))
  expect_lt(max(abs(coef(f7)$tt - coef(ref$value)$t0)), 2e-9)
  # A window in which the predictor really is constant is still counted.
  d$flat <- ifelse(sf::st_coordinates(d)[, 1] > 5e3, 1e7 + 1, 1e7)
  rf <- .otl_catch(fit_gwr_model(d, "yield", "flat", bandwidth = 30))
  expect_gt(rf$value$info$n_local_collinear, 100L)
  expect_true(any(grepl("local collinearity", rf$warnings)))
})

test_that("the extrapolation warning carries its counts, and can be dropped by class", {
  set.seed(4); n <- 150
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n, 0, 0.3)
  fit <- suppressMessages(fit_gwr_model(d, "z", "a", bandwidth = 30))
  nd <- sf::st_as_sf(data.frame(x = 5e5 + c(500, 3000, 4000, 9000), y = 5e6 + 500,
                                a = 0), coords = c("x", "y"), crs = 32632)
  r <- .otl_catch(predict(fit, newdata = nd))
  expect_identical(r$classes, "spatialkit_gwr_extrapolation")
  w <- r$conds[[1L]]
  expect_identical(w$n_far, 3L)
  expect_identical(w$n, 4L)
  expect_identical(w$calls, 1L)
  expect_equal(w$reach, fit$info$max_kernel_width)
  expect_length(w$nearest, 2L)
  expect_lt(w$nearest[1L], w$nearest[2L])
  expect_gt(w$nearest[1L], w$reach)
  expect_match(conditionMessage(w), "^predict\\.gwr_fit\\(\\): 3 of 4 location\\(s\\) are farther from every training point than the widest kernel of the fit reaches")
  # The documented way to drop it where extrapolation is intended.
  expect_silent(p <- suppressWarnings(predict(fit, newdata = nd),
                                      classes = "spatialkit_gwr_extrapolation"))
  expect_false(anyNA(p))
  # The warning that counts NA predictions carries its counts too.
  fx <- suppressMessages(fit_gwr_model(d, "z", "a", adaptive = FALSE, bandwidth = 150))
  rx <- .otl_catch(predict(fx, newdata = nd))
  expect_identical(rx$classes, "spatialkit_gwr_not_estimable")
  expect_identical(rx$conds[[1L]]$n_na, 3L)
  expect_identical(rx$conds[[1L]]$n, 4L)
  expect_identical(rx$conds[[1L]]$n_empty, 3L)
})

test_that("a caller that predicts in pieces can raise each counted warning once, with totals", {
  set.seed(4); n <- 150
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n, 0, 0.3)
  fit <- suppressMessages(fit_gwr_model(d, "z", "a", bandwidth = 30))
  set.seed(5); m <- 120
  nd <- sf::st_as_sf(data.frame(x = 5e5 + runif(m, -1500, 2500), y = 5e6 + runif(m, -1500, 2500),
                                a = rnorm(m)), coords = c("x", "y"), crs = 32632)
  whole <- .otl_catch(predict(fit, newdata = nd))
  expect_identical(whole$classes, "spatialkit_gwr_extrapolation")
  total <- whole$conds[[1L]]$n_far
  expect_gt(total, 20L)
  chunks <- split(seq_len(m), rep(1:6, each = 20))
  # In pieces: one warning per piece, each with that piece's counts.
  loose <- .otl_catch(unlist(lapply(chunks, function(i) predict(fit, newdata = nd[i, ]))))
  expect_identical(loose$classes, rep("spatialkit_gwr_extrapolation", 6L))
  expect_identical(sum(vapply(loose$conds, function(w) w$n_far, integer(1))), total)
  # Collected: one warning, the totals, the caller's name, the same values.
  tallied <- .otl_catch(.gwr_tally_predict_warnings(
    unlist(lapply(chunks, function(i) predict(fit, newdata = nd[i, ]))),
    caller = "predict_surface", n_total = m))
  expect_identical(tallied$classes, "spatialkit_gwr_extrapolation")
  w <- tallied$conds[[1L]]
  expect_identical(w$n_far, total)
  expect_identical(w$n, 120L)
  expect_identical(w$calls, 6L)
  expect_equal(w$nearest, whole$conds[[1L]]$nearest)
  expect_match(conditionMessage(w), sprintf(
    "^predict_surface\\(\\): %d of 120 location\\(s\\) are farther from every training point", total))
  expect_match(conditionMessage(w), "\\(totals over the 6 prediction calls that found any\\)")
  expect_equal(unname(tallied$value), unname(whole$value))
  # Without `n_total` the denominator is the locations of the calls that warned.
  t2 <- .otl_catch(.gwr_tally_predict_warnings(
    lapply(chunks[1:2], function(i) predict(fit, newdata = nd[i, ])), caller = "cv_gwr"))
  expect_identical(t2$conds[[1L]]$n, 40L)
  expect_match(conditionMessage(t2$conds[[1L]]), "^cv_gwr\\(\\): ")
  # Nested: the inner total is one more warning to the outer collector.
  nested <- .otl_catch(.gwr_tally_predict_warnings({
    .gwr_tally_predict_warnings(lapply(chunks[1:3], function(i) predict(fit, newdata = nd[i, ])),
                                caller = "inner")
    lapply(chunks[4:6], function(i) predict(fit, newdata = nd[i, ]))
  }, caller = "outer", n_total = m))
  expect_identical(nested$classes, "spatialkit_gwr_extrapolation")
  expect_identical(nested$conds[[1L]]$n_far, total)
  expect_identical(nested$conds[[1L]]$calls, 6L)
  # Nothing to report: the value, and no warning.  Other warnings pass through.
  inside <- nd[sf::st_coordinates(nd)[, 1] > 5e5 + 200 & sf::st_coordinates(nd)[, 1] < 5e5 + 800 &
                 sf::st_coordinates(nd)[, 2] > 5e6 + 200 & sf::st_coordinates(nd)[, 2] < 5e6 + 800, ]
  quiet <- .otl_catch(.gwr_tally_predict_warnings(predict(fit, newdata = inside), caller = "x"))
  expect_length(quiet$warnings, 0L)
  expect_length(quiet$value, nrow(inside))
  other <- .otl_catch(.gwr_tally_predict_warnings({ warning("something else"); 1 }, caller = "x"))
  expect_identical(other$warnings, "something else")
  # The collector leaves no trace behind, also after an error.
  expect_error(.gwr_tally_predict_warnings(stop("boom"), caller = "x"), "boom")
  expect_identical(.gwr_predict_state$tally, 0L)
  # The fixed-bandwidth warning about NA predictions is totalled the same way.
  fx <- suppressMessages(fit_gwr_model(d, "z", "a", adaptive = FALSE, bandwidth = 150))
  one <- .otl_catch(predict(fx, newdata = nd))
  tx <- .otl_catch(.gwr_tally_predict_warnings(
    lapply(chunks, function(i) predict(fx, newdata = nd[i, ])), caller = "predict_surface",
    n_total = m))
  na_one <- one$conds[[which(one$classes == "spatialkit_gwr_not_estimable")]]
  na_all <- tx$conds[[which(tx$classes == "spatialkit_gwr_not_estimable")]]
  expect_identical(sum(tx$classes == "spatialkit_gwr_not_estimable"), 1L)
  expect_identical(na_all$n_na, na_one$n_na)
  expect_identical(na_all$n_empty, na_one$n_empty)
  expect_match(conditionMessage(na_all), sprintf(
    "^predict_surface\\(\\): %d of 120 location\\(s\\) have no estimable local regression", na_one$n_na))
})

test_that("predict_surface() in chunks, wrapped in the collector, warns once with the whole grid's count", {
  skip_if_not_installed("ggplot2")
  set.seed(4); n <- 150
  d <- sf::st_as_sf(data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
                               a = rnorm(n)), coords = c("x", "y"), crs = 32632)
  d$z <- 1 + 2 * d$a + rnorm(n, 0, 0.3)
  fit <- suppressMessages(fit_gwr_model(d, "z", "a", bandwidth = 30))
  g <- expand.grid(x = 5e5 + seq(-1500, 2500, length.out = 20),
                   y = 5e6 + seq(-1500, 2500, length.out = 15))
  g$a <- 0
  grid <- sf::st_as_sf(g, coords = c("x", "y"), crs = 32632)
  one <- .otl_catch(predict_surface(fit, grid, chunk_size = 5000L))
  far <- one$conds[[which(one$classes == "spatialkit_gwr_extrapolation")]]
  # The change predict_surface() needs is this wrapper round its chunk loop.
  many <- .otl_catch(.gwr_tally_predict_warnings(
    predict_surface(fit, grid, chunk_size = 50L), caller = "predict_surface",
    n_total = nrow(grid)))
  expect_identical(sum(many$classes == "spatialkit_gwr_extrapolation"), 1L)
  got <- many$conds[[which(many$classes == "spatialkit_gwr_extrapolation")]]
  expect_identical(got$n_far, far$n_far)
  expect_identical(got$n, 300L)
  expect_identical(got$calls, 6L)
})

test_that("print() heads a GWR fit without naming GWmodel", {
  d <- .otl_layer()
  fit <- suppressMessages(fit_gwr_model(d, "yield", "t0", bandwidth = 30))
  out <- capture.output(print(fit))
  expect_identical(out[1L], "<GWR> spatial model fit")
  expect_false(any(grepl("GWmodel", out, fixed = TRUE)))
  expect_true(any(grepl("Eff. parameters", out, fixed = TRUE)))
})
