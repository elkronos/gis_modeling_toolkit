# tests/testthat/test-gp-boundary.R
# ---------------------------------------------------------------------------
# predict.bayesian_fit() holds the HSGP boundary L at its fitted value.
#
# brms 2.17 to 2.22 rebuild L = c * max(1, pooled range of the unique newdata
# rows centred on the training cmeans) at every predict call.  Appending the
# training extrema (.pin_gp_boundary_rows()) can only widen that range, so one
# row past the training envelope widened L and moved the prediction of EVERY
# row in the call, and predict_surface() depended on chunk_size.  The repair
# hands brms c * S_fit / S_new instead of c.  The arithmetic is tested here
# without Stan; the last two tests fit a real model and are opt-in, like
# test-bayes-smoke.R.
# ---------------------------------------------------------------------------

# brms:::choose_L() applied the way brms:::.data_gp() applies it.
brms_L <- function(xy, cmeans, c) {
  Xc <- sweep(unique(as.matrix(xy)[, 1:2, drop = FALSE]), 2L, cmeans)
  c * max(1, max(Xc) - min(Xc))
}

gpb_train_xy <- function(n = 60, seed = 3) {
  set.seed(seed)
  xy <- cbind(runif(n, -2, 1.5), runif(n, -1, 2.5))
  rbind(xy, xy[1:5, ])                   # replicated locations, as brms sees them
}

gpb_fake_fit <- function(xy, c = 1.6, store_cmeans = TRUE) {
  spec <- spatialkit:::.gp_basis_spec(xy, c(lower = 0.3, upper = 1))
  info <- list(gp_c = c, gp_S = spec$S,
               gp_xy_range = list(x = range(xy[, 1]), y = range(xy[, 2])))
  if (store_cmeans) info$gp_cmeans <- spec$cmeans
  structure(list(info = info,
                 engine = list(data = data.frame(..x = xy[, 1], ..y = xy[, 2],
                                                 a = 0))),
            class = c("bayesian_fit", "spatial_fit"))
}

test_that(".gp_basis_spec returns the centre brms uses", {
  xy <- gpb_train_xy()
  spec <- spatialkit:::.gp_basis_spec(xy, c(lower = 0.3, upper = 1))
  u <- unique(xy)
  expect_equal(spec$cmeans, unname(colMeans(u)))
  expect_equal(brms_L(xy, spec$cmeans, 1), spec$S)
})

test_that(".gp_c_scale makes brms rebuild the fitted L from any rows", {
  xy <- gpb_train_xy()
  spec <- spatialkit:::.gp_basis_spec(xy, c(lower = 0.3, upper = 1))
  c_fit <- 1.6
  L_fit <- brms_L(xy, spec$cmeans, c_fit)
  expect_equal(L_fit, c_fit * spec$S)

  set.seed(9)
  cases <- list(
    wider    = rbind(xy, c(4, 0.5)),                     # one row far out
    both     = rbind(xy[1:3, ], c(-5, 6), c(-5, 6)),     # widened, duplicated
    narrow   = matrix(c(0.1, 0.2, 0.15, 0.25), 2),       # range < 1: brms's floor
    single   = matrix(c(0.3, 0.4), 1),
    training = xy
  )
  for (nm in names(cases)) {
    nd <- cases[[nm]]
    s  <- spatialkit:::.gp_c_scale(nd, spec$cmeans, spec$S)
    expect_equal(brms_L(nd, spec$cmeans, c_fit * s), L_fit, tolerance = 1e-12,
                 info = nm)
  }
  expect_identical(spatialkit:::.gp_c_scale(xy, spec$cmeans, spec$S), 1)
  expect_identical(spatialkit:::.gp_c_scale(matrix(NA_real_, 1, 2),
                                            spec$cmeans, spec$S), 1)
})

test_that(".pin_gp_boundary_rows scales c only when the rows widen the range", {
  xy  <- gpb_train_xy()
  fit <- gpb_fake_fit(xy)
  inside <- data.frame(a = 0, ..x = c(0, 0.5), ..y = c(1, 0))
  p_in <- spatialkit:::.pin_gp_boundary_rows(fit, inside)
  expect_identical(p_in$n_pad, 2L)
  expect_identical(p_in$c_scale, 1)       # the padding rows alone are exact
  expect_identical(p_in$beyond, c(FALSE, FALSE))

  out <- rbind(inside, data.frame(a = 0, ..x = 3, ..y = 3))
  p_out <- spatialkit:::.pin_gp_boundary_rows(fit, out)
  expect_lt(p_out$c_scale, 1)
  expect_equal(brms_L(cbind(p_out$df$..x, p_out$df$..y), fit$info$gp_cmeans,
                      fit$info$gp_c * p_out$c_scale),
               fit$info$gp_c * fit$info$gp_S, tolerance = 1e-12)

  # A fit saved before gp_cmeans was stored reads the centre off the brmsfit's
  # data, as brms did, and gets the same answer.
  old <- gpb_fake_fit(xy, store_cmeans = FALSE)
  expect_equal(spatialkit:::.pin_gp_boundary_rows(old, out)[c("c_scale", "beyond")],
               p_out[c("c_scale", "beyond")])
})

test_that(".pin_gp_boundary_rows flags rows beyond +/- L of the training centre", {
  xy  <- gpb_train_xy()
  fit <- gpb_fake_fit(xy)
  cm  <- fit$info$gp_cmeans
  L   <- fit$info$gp_c * fit$info$gp_S
  nd  <- data.frame(a = 0,
                    ..x = cm[1] + c(0, 0.99 * L, 1.01 * L, 0, -1.2 * L),
                    ..y = cm[2] + c(0, 0, 0, -1.01 * L, 0))
  expect_identical(spatialkit:::.pin_gp_boundary_rows(fit, nd)$beyond,
                   c(FALSE, FALSE, TRUE, TRUE, TRUE))
})

test_that(".scale_gp_c rescales c in the gp() term the package builds, and nothing else", {
  f <- stats::as.formula(paste("resp ~ z + w +",
                               spatialkit:::.gp_formula_term(12L, 1.5)))
  obj <- list(formula = list(formula = f), other = 1)
  got <- spatialkit:::.scale_gp_c(obj, 0.5)
  gp_call <- got$formula$formula[[3L]][[3L]]
  expect_identical(gp_call[[1L]], as.name("gp"))
  expect_identical(gp_call[["c"]], 0.75)
  gp_call[["c"]] <- 1.5
  expect_identical(gp_call, f[[3L]][[3L]])          # k, scale, iso untouched
  expect_identical(got$formula$formula[[2L]], f[[2L]])
  expect_identical(got$other, 1)

  expect_null(spatialkit:::.scale_gp_c(list(formula = list(formula = resp ~ z)), 0.5))
})


# ---------------------------------------------------------------------------
# Real brms.  Opt-in via SPATIALKIT_TEST_BRMS; see test-bayes-smoke.R.
# ---------------------------------------------------------------------------

gpb_points <- function(n = 50, seed = 11) {
  set.seed(seed)
  x <- runif(n, 0, 1000)
  y <- runif(n, 0, 1000)
  a <- rnorm(n)
  # A surface with real spatial structure, so the GP term carries the fit.
  resp <- 2 * sin(x / 250) + 1.5 * cos(y / 300) + 0.5 * a + rnorm(n, sd = 0.3)
  sf::st_as_sf(data.frame(x = x, y = y, a = a, resp = resp),
               coords = c("x", "y"), crs = 32632)
}

gpb_cache <- new.env(parent = emptyenv())
gpb_fit <- function() {
  if (is.null(gpb_cache$fit)) {
    fit <- NULL
    utils::capture.output(
      suppressWarnings(suppressMessages(
        fit <- fit_bayesian_spatial_model(
          gpb_points(), response_var = "resp", predictor_vars = "a",
          chains = 1, iter = 300, warmup = 150, cores = 1,
          compute_loo = FALSE, check_convergence = FALSE, seed = 4321)
      )),
      type = "output"
    )
    gpb_cache$fit <- fit
  }
  gpb_cache$fit
}

gpb_sf <- function(x, y) {
  sf::st_as_sf(data.frame(x = x, y = y, a = 0), coords = c("x", "y"), crs = 32632)
}

test_that("an out-of-envelope row does not move the other rows' predictions", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  fit <- gpb_fit()
  bb  <- sf::st_bbox(fit$data_sf)
  ix  <- c(500, 300, 700); iy <- c(500, 700, 250)
  far_x <- bb[["xmax"]] + 0.1 * (bb[["xmax"]] - bb[["xmin"]])
  far_y <- bb[["ymax"]] + 0.1 * (bb[["ymax"]] - bb[["ymin"]])

  # The far row must widen the range brms builds L from, or this proves nothing.
  cs <- fit$info$coord_scaling
  tr <- cbind(fit$engine$data$..x, fit$engine$data$..y)
  cm <- colMeans(unique(tr))
  far_c <- c((far_x - cs$x_center) / cs$x_scale,
             (far_y - cs$y_center) / cs$y_scale) - cm
  expect_gt(max(far_c), max(sweep(tr, 2L, cm)))

  p_inner <- suppressMessages(predict(fit, newdata = gpb_sf(ix, iy)))
  p_alone <- vapply(1:3, function(i)
    suppressMessages(predict(fit, newdata = gpb_sf(ix[i], iy[i]))), numeric(1))
  p_far   <- suppressMessages(predict(fit, newdata = gpb_sf(c(ix, far_x),
                                                            c(iy, far_y))))
  d_far   <- suppressMessages(predict(fit, newdata = gpb_sf(c(ix, far_x),
                                                            c(iy, far_y)),
                                      draws = TRUE))
  d_inner <- suppressMessages(predict(fit, newdata = gpb_sf(ix, iy), draws = TRUE))

  expect_true(all(is.finite(p_far)))
  expect_equal(p_alone, p_inner, tolerance = 1e-10)
  # Before the fix these moved by ~1e-2 on a response of SD ~1.5.
  expect_equal(p_far[1:3], p_inner, tolerance = 1e-10)
  expect_equal(d_far[, 1:3], d_inner, tolerance = 1e-10)
  # In-sample prediction is untouched.
  expect_equal(suppressMessages(predict(fit, newdata = fit$data_sf)), fitted(fit),
               tolerance = 1e-10)

  # A row past the edge of the basis is NA, with a warning, and still moves
  # nothing else.
  L  <- fit$info$gp_c * fit$info$gp_S
  wx <- cs$x_center + (cm[1] + 1.1 * L) * cs$x_scale
  expect_warning(
    p_wall <- suppressMessages(predict(fit, newdata = gpb_sf(c(ix, wx), c(iy, 500)))),
    "^predict\\.bayesian_fit\\(\\): 1 of 4 row\\(s\\) of `newdata` lie beyond the GP boundary")
  expect_true(is.na(p_wall[4]))
  expect_equal(p_wall[1:3], p_inner, tolerance = 1e-10)
})

test_that("predict_surface() on a grid past the training bbox does not depend on chunk_size", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SPATIALKIT_TEST_BRMS")),
          "set SPATIALKIT_TEST_BRMS=true to run the Stan smoke tests")
  skip_if_not_installed("brms")

  fit <- gpb_fit()
  bb  <- sf::st_bbox(fit$data_sf)
  w <- bb[["xmax"]] - bb[["xmin"]]; h <- bb[["ymax"]] - bb[["ymin"]]
  g <- expand.grid(x = seq(bb[["xmin"]] - 0.15 * w, bb[["xmax"]] + 0.15 * w,
                           length.out = 16),
                   y = seq(bb[["ymin"]] - 0.15 * h, bb[["ymax"]] + 0.15 * h,
                           length.out = 16))
  grid <- gpb_sf(g$x, g$y)
  inb  <- g$x >= bb[["xmin"]] & g$x <= bb[["xmax"]] &
          g$y >= bb[["ymin"]] & g$y <= bb[["ymax"]]

  one   <- suppressMessages(predict_surface(fit, grid = grid, chunk_size = 1e6))$.pred
  many  <- suppressMessages(predict_surface(fit, grid = grid, chunk_size = 23))$.pred
  in_bb <- suppressMessages(predict_surface(fit, grid = grid[inb, ],
                                            chunk_size = 1e6))$.pred

  expect_true(all(is.finite(one)))
  expect_equal(many, one, tolerance = 1e-10)
  expect_equal(one[inb], in_bb, tolerance = 1e-10)
})
