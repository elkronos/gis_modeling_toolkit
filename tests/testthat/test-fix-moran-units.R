# Moran's I of the cell means, and its z, do not change under a + b * y.
# .morans_i_for_k() refused a response in small units (an absolute floor on
# the residuals' sum of squares) and regressed the raw cell means of a
# response at a large offset, so resolution_profile() left moran_i and
# moran_z NA for a response in units of 1e-9, and determine_optimal_levels()
# fell back to the geometric ranking there and ranked other counts at an
# offset of 1e15.

# A smooth field and noise on 600 points; the spread as stored at 1e15.
.mu_layer <- function(n = 600, seed = 3) {
  withr::with_seed(seed, {
    x <- runif(n, 0, 1e4); y <- runif(n, 0, 1e4)
    w <- rnorm(n)
    z <- sin(x / 2000) + cos(y / 1500) + rnorm(n, sd = 0.5)
  })
  list(xy = cbind(x, y), w = w, z = (1e15 + z) - 1e15)
}

.mu_cells <- function(xy, k, seed = 1) {
  withr::with_seed(seed, stats::kmeans(xy, k, nstart = 1, iter.max = 50)$cluster)
}

test_that(".morans_i_for_k() gives one I and one z in every unit and at every offset", {
  lay <- .mu_layer()
  none <- matrix(numeric(0), nrow(lay$xy), 0L)
  for (k in c(20L, 40L)) {
    cl <- .mu_cells(lay$xy, k)
    for (grouped in c(FALSE, TRUE)) {
      ref <- .morans_i_for_k(lay$xy, lay$z, none, cl, grouped = grouped)
      expect_true(all(is.finite(ref)))
      # 1e160 and 1e-170 lie beyond the power-of-two scaling bounds and
      # within the squares' range.
      for (b in c(1e-200, 1e-170, 1e-12, 1e-9, 1e-6, 1e6, 1e160, 1e200)) {
        got <- .morans_i_for_k(lay$xy, lay$z * b, none, cl, grouped = grouped)
        expect_equal(got, ref, tolerance = 1e-8,
                     info = sprintf("k = %d, grouped = %s, z * %g", k, grouped, b))
      }
      # 1e10 to 1e12 show their spread as stored and are not centred first.
      for (a in c(1e8, 1e10, 1e11, 1e12, 1e15)) {
        got <- .morans_i_for_k(lay$xy, a + lay$z, none, cl, grouped = grouped)
        expect_equal(got, ref, tolerance = 1e-8,
                     info = sprintf("k = %d, grouped = %s, %g + z", k, grouped, a))
      }
      # Offsets beyond the pooling's bounds, the spread as stored there: lost
      # in the raw values (1e-14 of the offset) or shown in them (1e-11).
      for (a in c(1e-200, 1e200)) for (f in c(1e-14, 1e-11)) {
        v <- a + lay$z * (a * f)
        got <- .morans_i_for_k(lay$xy, v, none, cl, grouped = grouped)
        want <- .morans_i_for_k(lay$xy, (v - a) / (a * f), none, cl, grouped = grouped)
        expect_true(all(is.finite(want)))
        expect_equal(got, want, tolerance = 1e-8,
                     info = sprintf("k = %d, grouped = %s, %g + z * %g", k, grouped, a, a * f))
      }
    }
  }
})

test_that(".morans_i_for_k() on the residuals of the cell means is the same in every unit", {
  lay <- .mu_layer()
  pr <- matrix(lay$w, ncol = 1L)
  y <- 3 + 2 * lay$w + lay$z
  cl <- .mu_cells(lay$xy, 40L)
  ref <- .morans_i_for_k(lay$xy, y, pr, cl)
  expect_true(all(is.finite(ref)))
  for (b in c(1e-200, 1e-170, 1e-12, 1e-9, 1e6, 1e160, 1e200))
    expect_equal(.morans_i_for_k(lay$xy, y * b, pr, cl), ref, tolerance = 1e-8,
                 info = sprintf("y * %g", b))
})

test_that(".morans_i_for_k() reads a response with a trend at a large offset", {
  # The raw values show their spread through the trend, so the response is
  # not centred first, and the residuals of the cell means at 1e15 hold the
  # rounding of the regression there; it is regressed again about its centre.
  lay <- .mu_layer()
  pr <- matrix(lay$w, ncol = 1L)
  v15 <- 1e15 + 1e5 * lay$w + lay$z           # the spread as stored at 1e15
  v0 <- v15 - 1e15
  for (k in c(20L, 40L)) {
    cl <- .mu_cells(lay$xy, k)
    for (grouped in c(FALSE, TRUE)) {
      info <- sprintf("k = %d, grouped = %s", k, grouped)
      ref <- .morans_i_for_k(lay$xy, v0, pr, cl, grouped = grouped)
      expect_true(all(is.finite(ref)), info = info)
      for (a in c(1e6, 1e10, 1e13, 1e14, 1e15))
        expect_equal(.morans_i_for_k(lay$xy, a + v0, pr, cl, grouped = grouped), ref,
                     tolerance = 1e-8, info = sprintf("%s, offset %g", info, a))
      # The predictor fits it exactly: the residuals are the rounding of the
      # stored values, and stay refused.
      ex <- .morans_i_for_k(lay$xy, 1e15 + 1e5 * lay$w, pr, cl, grouped = grouped)
      expect_true(all(is.na(ex)), info = paste(info, "exact fit at 1e15"))
    }
  }
})

test_that(".morans_i_for_k() still refuses cell means that differ by rounding alone", {
  lay <- .mu_layer()
  none <- matrix(numeric(0), nrow(lay$xy), 0L)
  pr <- matrix(lay$w, ncol = 1L)
  cl <- .mu_cells(lay$xy, 40L)
  na2 <- function(v, info) {
    expect_identical(names(v), c("I", "z"), info = info)
    expect_true(all(is.na(v)), info = info)
  }
  n <- nrow(lay$xy)
  for (v in c(0, 5, 1e-12, 1e-200, 1e15, 1e200))
    na2(.morans_i_for_k(lay$xy, rep(v, n), none, cl), paste("constant", v))
  # 0.1 beside 0.3 / 3: values that differ only by their rounding.
  na2(.morans_i_for_k(lay$xy, rep(c(0.1, 0.3 / 3), length.out = n), none, cl),
      "0.1 and 0.3 / 3")
  # The predictors fit the response exactly, in any unit.
  for (b in c(1e-200, 1e-12, 1, 1e10, 1e200))
    na2(.morans_i_for_k(lay$xy, b * (3 + 2 * lay$w), pr, cl), paste("exact fit times", b))
})

test_that("resolution_profile() gives moran_i and moran_z in small units", {
  lay <- .mu_layer()
  pts <- sf::st_as_sf(data.frame(x = lay$xy[, 1], y = lay$xy[, 2], z = lay$z),
                      coords = c("x", "y"), crs = 32615)
  prof <- function(b) {
    d <- pts; d$v <- d$z * b
    suppressWarnings(resolution_profile(d, "v", n_levels = 4, seed = 1, sac = 2000))
  }
  ref <- prof(1)
  expect_true(sum(is.finite(ref$moran_z)) >= 3L)
  for (b in c(1e-9, 1e-12)) {
    p <- prof(b)
    expect_identical(p$levels, ref$levels)
    expect_identical(is.finite(p$moran_z), is.finite(ref$moran_z), info = paste("z *", b))
    expect_equal(p$moran_i, ref$moran_i, tolerance = 1e-8, info = paste("z *", b))
    expect_equal(p$moran_z, ref$moran_z, tolerance = 1e-8, info = paste("z *", b))
  }
})

test_that("resolution_profile() gives one moran_i and moran_z at offsets the raw values show", {
  lay <- .mu_layer()
  pts <- sf::st_as_sf(data.frame(x = lay$xy[, 1], y = lay$xy[, 2], z = lay$z),
                      coords = c("x", "y"), crs = 32615)
  prof <- function(a) {
    d <- pts; d$v <- a + d$z
    suppressWarnings(resolution_profile(d, "v", n_levels = 4, seed = 1, sac = 2000))
  }
  ref <- prof(0)
  for (a in c(1e11, 1e12)) {
    p <- prof(a)
    expect_identical(p$levels, ref$levels)
    expect_equal(p$moran_i, ref$moran_i, tolerance = 1e-8, info = paste(a, "+ z"))
    expect_equal(p$moran_z, ref$moran_z, tolerance = 1e-8, info = paste(a, "+ z"))
  }
})

test_that("determine_optimal_levels(criterion = \"morans_i\") ranks the same counts in any unit and at any offset", {
  withr::local_seed(11)
  K <- 40
  cen <- cbind(runif(K, 0, 1e5), runif(K, 0, 1e5))
  id <- sample(K, 2000, TRUE)
  cx <- cen[id, 1] + rnorm(2000, sd = 800); cy <- cen[id, 2] + rnorm(2000, sd = 800)
  cz <- sin(cx / 2e4) + cos(cy / 1.5e4) + rnorm(2000, sd = 0.5)
  cz <- (1e15 + cz) - 1e15
  lay <- sf::st_as_sf(data.frame(x = cx, y = cy, z = cz), coords = c("x", "y"), crs = 32615)
  run <- function(a, b) {
    d <- lay; d$v <- a + d$z * b
    lines <- capture_spatialkit_log(
      r <- suppressWarnings(determine_optimal_levels(d, max_levels = 80, response_var = "v",
                                                     criterion = "morans_i",
                                                     on_no_elbow = "warn")),
      level = logger::WARN)
    list(r = r, log = lines)
  }
  ref <- run(0, 1)
  dr <- attr(ref$r, "diagnostics")
  expect_false(is.null(dr))
  expect_true(any(is.finite(dr$moran_z)))
  for (cs in list(c(0, 1e-9), c(0, 1e-12), c(1e15, 1))) {
    info <- sprintf("%g + z * %g", cs[1], cs[2])
    got <- run(cs[1], cs[2])
    d <- attr(got$r, "diagnostics")
    expect_identical(as.vector(got$r), as.vector(ref$r), info = info)
    expect_false(is.null(d), info = info)
    expect_identical(d$eval_ks, dr$eval_ks, info = info)
    expect_equal(d$moran_z, dr$moran_z, tolerance = 1e-8, info = info)
    expect_equal(d$moran_i, dr$moran_i, tolerance = 1e-8, info = info)
    expect_false(log_has(got$log, "falling back to geometric"), info = info)
  }
})

test_that("a constant response falls back in any unit, and the log lines name cell means that differ by rounding alone", {
  withr::local_seed(11)
  K <- 12
  cen <- cbind(runif(K, 0, 1e5), runif(K, 0, 1e5))
  id <- sample(K, 600, TRUE)
  lay <- sf::st_as_sf(data.frame(x = cen[id, 1] + rnorm(600, sd = 800),
                                 y = cen[id, 2] + rnorm(600, sd = 800)),
                      coords = c("x", "y"), crs = 32615)
  geo <- as.vector(determine_optimal_levels(lay, max_levels = 30, on_no_elbow = "warn"))
  for (v in c(1e-9, 1, 1e15)) {
    lay$const <- v
    for (crit in c("morans_i", "combined")) {
      info <- sprintf("%s, constant %g", crit, v)
      lines <- capture_spatialkit_log(
        r <- suppressWarnings(determine_optimal_levels(lay, max_levels = 30,
                                                       response_var = "const", criterion = crit,
                                                       on_no_elbow = "warn")),
        level = logger::WARN)
      expect_true(log_has(lines, paste0("Moran's I could not be computed at any k from .*",
                                        "or the cell means differ by rounding alone\\)")),
                  info = info)
      if (crit == "morans_i") {
        # Rounding is not ranked: the geometric ranking comes back.
        expect_identical(as.vector(r), geo, info = info)
        expect_null(attr(r, "diagnostics"), info = info)
        expect_true(log_has(lines, "falling back to geometric"), info = info)
      } else {
        expect_true(all(is.na(attr(r, "diagnostics")$moran_z)), info = info)
      }
    }
  }
})

test_that("determine_optimal_levels(criterion = \"morans_i\") with a predictor ranks the same counts at a large offset", {
  withr::local_seed(11)
  K <- 40
  cen <- cbind(runif(K, 0, 1e5), runif(K, 0, 1e5))
  id <- sample(K, 2000, TRUE)
  cx <- cen[id, 1] + rnorm(2000, sd = 800); cy <- cen[id, 2] + rnorm(2000, sd = 800)
  pw <- rnorm(2000)
  cz <- sin(cx / 2e4) + cos(cy / 1.5e4) + rnorm(2000, sd = 0.5)
  v15 <- 1e15 + 1e5 * pw + cz                  # the spread as stored at 1e15
  lay <- sf::st_as_sf(data.frame(x = cx, y = cy, p = pw, v0 = v15 - 1e15, v15 = v15),
                      coords = c("x", "y"), crs = 32615)
  run <- function(v) {
    lines <- capture_spatialkit_log(
      r <- suppressWarnings(determine_optimal_levels(lay, max_levels = 80, response_var = v,
                                                     predictor_vars = "p",
                                                     criterion = "morans_i",
                                                     on_no_elbow = "warn")),
      level = logger::WARN)
    list(r = r, log = lines)
  }
  ref <- run("v0")
  got <- run("v15")
  dr <- attr(ref$r, "diagnostics"); d <- attr(got$r, "diagnostics")
  expect_false(is.null(dr))
  expect_true(any(is.finite(dr$moran_z)))
  expect_false(is.null(d))
  expect_identical(as.vector(got$r), as.vector(ref$r))
  expect_identical(d$eval_ks, dr$eval_ks)
  expect_equal(d$moran_z, dr$moran_z, tolerance = 1e-8)
  expect_equal(d$moran_i, dr$moran_i, tolerance = 1e-8)
  expect_false(log_has(got$log, "falling back to geometric"))
})
