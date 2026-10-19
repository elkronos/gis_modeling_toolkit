# tests/testthat/test-fix-variogram-callers.R
# ---------------------------------------------------------------------------
# The callers of estimate_sac_range() keep working with what it returns now:
# an accepted range, the two refusals it did not have before (no spatial
# structure, a flat variogram) and a model with two structures.
# resolution_profile(), make_folds(auto_range = TRUE),
# summarize_by_cell(deff = "variogram"), kriging_adequacy() and plot() read
# the class, the value, `rejected_reason`, `variogram`, `variogram_model`,
# `crs` and `detrended`, and none of those changed its meaning.
# ---------------------------------------------------------------------------

NO_STRUCTURE_REASON <- "no spatial structure detected at the nearest-neighbour scale"

# 240 points on a 1000 m square with an exponential field (effective range
# 150, nugget 0.2 on a sill of 1), a white-noise column, and a 4 x 4 grid of
# cells they are assigned to.
fvc_setup <- local({
  val <- NULL
  function() {
    if (is.null(val)) {
      set.seed(4401); n <- 240
      x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
      D <- as.matrix(stats::dist(cbind(x, y)))
      z <- drop(crossprod(chol(0.8 * exp(-3 * D / 150) + diag(0.2, n)), rnorm(n)))
      pts <- sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z, noise = rnorm(n)),
                          coords = c("x", "y"), crs = 32632)
      bnd <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(pts)))
      cells <- create_grid_polygons(bnd, target_cells = 16, type = "square")
      val <<- list(pts = pts, cells = cells,
                   asg = suppressWarnings(assign_features_to_polygons(pts, cells)))
    }
    val
  }
})
fvc_warnings <- function(expr) {
  w <- list()
  val <- withCallingHandlers(expr, warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val, messages = vapply(w, conditionMessage, character(1)),
       classes = lapply(w, class))
}

test_that("the object keeps its class, its value and the attributes callers read", {
  skip_if_not_installed("gstat")
  s <- fvc_setup()
  ok <- estimate_sac_range(s$pts, "z")
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))
  expect_true(is.finite(ok)); expect_true(is.na(no))
  expect_identical(attr(no, "rejected_reason"), NO_STRUCTURE_REASON)
  for (r in list(ok, no)) {
    expect_s3_class(r, c("sac_range", "numeric"), exact = TRUE)
    expect_true(is.numeric(r) && length(r) == 1L)
    expect_s3_class(attr(r, "variogram"), "data.frame")
    expect_true(all(c("np", "dist", "gamma") %in% names(attr(r, "variogram"))))
    vm <- attr(r, "variogram_model")
    expect_s3_class(vm, "data.frame")
    expect_true(all(c("model", "psill", "range") %in% names(vm)))
    expect_s3_class(attr(r, "crs"), "crs")
    expect_identical(attr(r, "detrended"), FALSE)
    expect_identical(sac_nugget(r), attr(r, "nugget"))
    expect_true(attr(vm, "converged"))
  }
  expect_null(attr(ok, "rejected_reason"))
  # A nugget and one structure, or two when the equal-lag fit has a range
  # at least twice as long that the adaptive variogram supports (it does
  # here: 307 against 120).  Either way the value is the effective range of
  # the longest structure in the model attached.
  vm <- attr(ok, "variogram_model")
  expect_identical(as.character(vm$model)[1:2], c("Nug", "Exp"))
  expect_identical(attr(ok, "structures"), nrow(vm) - 1L)
  expect_equal(as.numeric(ok), 3 * max(vm$range))
  one <- estimate_sac_range(s$pts, "z", lags = "equal")
  expect_identical(as.character(attr(one, "variogram_model")$model), c("Nug", "Exp"))
  expect_equal(as.numeric(one), 3 * attr(one, "variogram_model")$range[2])
})

test_that("make_folds(auto_range = TRUE) sizes blocks from a range and says why when there is none", {
  skip_if_not_installed("gstat")
  s <- fvc_setup()
  ok <- estimate_sac_range(s$pts, "z")
  f <- suppressMessages(make_folds(s$pts, k = 3, method = "block_kfold", auto_range = TRUE,
                                   response_var = "z", seed = 1))
  expect_equal(as.numeric(f$params$sac_range), as.numeric(ok))
  expect_equal(as.numeric(f$params$block_size), as.numeric(ok))
  res <- fvc_warnings(make_folds(s$pts, k = 3, method = "block_kfold", auto_range = TRUE,
                                 response_var = "noise", seed = 1))
  hit <- grepl("no autocorrelation range was identified", res$messages)
  expect_identical(sum(hit), 1L)
  expect_match(res$messages[hit], NO_STRUCTURE_REASON, fixed = TRUE)
  expect_match(res$messages[hit], "falling back to geometric blocks", fixed = TRUE)
  expect_gt(res$value$params$grid_nx * res$value$params$grid_ny, 1L)
})

test_that("resolution_profile() takes no nugget from a layer with no structure", {
  skip_if_not_installed("gstat")
  s <- fvc_setup()
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))
  res <- fvc_warnings(resolution_profile(s$pts, "noise", sac = no, n_levels = 4, nstart = 2))
  p <- res$value
  # Cp's noise variance is Mallows' own, read at the finest level, and there
  # is no reliability curve: the pure-nugget model has no correlation to give.
  expect_identical(attr(p, "cp_noise")$source, "finest-level residual mean square")
  expect_true(all(is.finite(p$cp)))
  expect_true(all(is.na(p$reliability)))
  expect_false(any(grepl("nugget is 0", res$messages, fixed = TRUE)))
  # The same when the profile estimates the variogram itself.
  p2 <- suppressWarnings(resolution_profile(s$pts, "noise", n_levels = 4, nstart = 2))
  expect_identical(attr(attr(p2, "sac"), "rejected_reason"), NO_STRUCTURE_REASON)
  expect_identical(attr(p2, "cp_noise")$source, "finest-level residual mean square")
  # And with an accepted range the nugget is read from the model, as before.
  ok <- estimate_sac_range(s$pts, "z")
  p3 <- suppressWarnings(resolution_profile(s$pts, "z", sac = ok, n_levels = 4, nstart = 2))
  expect_identical(attr(p3, "cp_noise")$source, "variogram nugget")
  expect_equal(attr(p3, "cp_noise")$value, sac_nugget(ok))
  expect_true(any(is.finite(p3$reliability)))
})

test_that("a fit with a minimum, refused as past the lags, still gives cp its nugget", {
  skip_if_not_installed("gstat")
  # The refusal is about the range.  The nugget of a fit that has a minimum
  # was fitted at the short lags, and the model says so (converged = TRUE).
  s <- fvc_setup()
  sac <- suppressWarnings(estimate_sac_range(s$pts, "z", range_frac = 0.05))
  expect_true(is.na(sac))
  expect_identical(attr(sac, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_identical(attr(attr(sac, "variogram_model"), "converged", exact = TRUE), TRUE)
  res <- fvc_warnings(resolution_profile(s$pts, "z", sac = sac, n_levels = 4, nstart = 2))
  expect_identical(attr(res$value, "cp_noise")$source, "variogram nugget")
  expect_equal(attr(res$value, "cp_noise")$value, sac_nugget(sac))
  expect_true(all(is.na(res$value$reliability)))
  expect_true(any(grepl("`cp` still takes its nugget", res$messages, fixed = TRUE)))
  # A variogram with no sill inside four cutoffs has no minimum: the model is
  # marked, and neither its nugget nor its sill is used.
  tr <- s$pts
  tr$t <- 0.02 * sf::st_coordinates(tr)[, 1] + tr$z
  sac2 <- suppressWarnings(estimate_sac_range(tr, "t"))
  expect_identical(attr(sac2, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_identical(attr(attr(sac2, "variogram_model"), "converged", exact = TRUE), FALSE)
  expect_equal(attr(sac2, "rejected_range"), 4 * attr(sac2, "cutoff_dist"))
  res2 <- fvc_warnings(resolution_profile(tr, "t", sac = sac2, n_levels = 4, nstart = 2))
  expect_identical(attr(res2$value, "cp_noise")$source, "finest-level residual mean square")
  expect_false(any(grepl("still takes its nugget", res2$messages, fixed = TRUE)))
})

test_that("summarize_by_cell(deff = 'variogram') uses an accepted model and falls back on a refusal", {
  skip_if_not_installed("gstat")
  s <- fvc_setup()
  ok <- estimate_sac_range(s$pts, "z")
  sm <- summarize_by_cell(s$asg, "z", deff = "variogram", sac = ok, estimate = "mean")
  da <- attr(sm, "deff_applied")
  expect_identical(da$method, "variogram")
  expect_true(all(da$deff[is.finite(da$deff)] >= 1))
  expect_true(any(da$deff > 1.05, na.rm = TRUE))
  # No structure: the model that comes with the refusal is a pure nugget,
  # under which distinct observations are uncorrelated.  The design effect
  # is 1 exactly, by the variogram route, and nothing falls back.  (It used
  # to be set aside with the other refusals: one warning of class
  # spatialkit_deff_fallback, and Kish's correction.)
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))
  expect_identical(attr(no, "rejected_reason"), NO_STRUCTURE_REASON)
  res <- fvc_warnings(summarize_by_cell(s$asg, "noise", deff = "variogram", sac = no,
                                        estimate = "mean"))
  expect_length(res$messages, 0L)
  da <- attr(res$value, "deff_applied")
  expect_identical(da$method, "variogram")
  expect_true(all(da$deff[is.finite(da$deff)] == 1))
  expect_true(all(res$value$deff_method == "variogram"))
  # The standard errors are the uncorrected ones, and the same when the
  # function estimates the variogram itself.
  plain <- summarize_by_cell(s$asg, "noise", estimate = "mean")
  expect_equal(res$value$..se_resp_noise, plain$..se_resp_noise)
  own <- fvc_warnings(summarize_by_cell(s$asg, "noise", deff = "variogram", estimate = "mean"))
  expect_false(any(vapply(own$classes, function(k) "spatialkit_deff_fallback" %in% k,
                          logical(1))))
  expect_equal(own$value$..se_resp_noise, plain$..se_resp_noise)
  # Any other refusal is still set aside, with the classed fallback warning
  # when nothing can replace it.
  far <- suppressWarnings(estimate_sac_range(s$pts, "z", range_frac = 0.05))
  expect_identical(attr(far, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  res2 <- fvc_warnings(summarize_by_cell(s$asg, NULL, predictor_vars = "z",
                                         deff = "variogram", sac = far, estimate = "mean"))
  cls <- vapply(res2$classes, function(k) "spatialkit_deff_fallback" %in% k, logical(1))
  expect_identical(sum(cls), 1L)
  expect_match(res2$messages[cls], "fitted range exceeds the largest lag fitted", fixed = TRUE)
})

test_that("kriging_adequacy() kriges with the model, whichever verdict the range got", {
  skip_if_not_installed("gstat")
  s <- fvc_setup()
  ok <- estimate_sac_range(s$pts, "z")
  ka <- kriging_adequacy(s$asg, "z", s$cells, sac = ok, k = 3)
  expect_s3_class(ka, "kriging_adequacy")
  expect_true(all(is.finite(ka$kr_var[ka$n > 0])))
  # With no structure the model is a pure nugget.  It is used, with a warning
  # that the range was not identified and why.
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))
  res <- fvc_warnings(kriging_adequacy(s$asg, "noise", s$cells, sac = no, k = 3))
  expect_s3_class(res$value, "kriging_adequacy")
  expect_true(any(grepl(paste0("range was not identified (", NO_STRUCTURE_REASON, ")"),
                        res$messages, fixed = TRUE)))
})

test_that("a model with two structures is read by every caller", {
  skip_if_not_installed("gstat")
  s <- fvc_setup()
  ok <- estimate_sac_range(s$pts, "z")
  # Built by hand the way estimate_sac_range() builds one: nugget, a short
  # structure and a long one, with the long range as the value.
  vm <- gstat::vgm(psill = 0.3, model = "Exp", range = 150,
                   add.to = gstat::vgm(psill = 0.5, model = "Exp", range = 20, nugget = 0.2))
  attr(vm, "converged") <- TRUE
  two <- ok
  attr(two, "variogram_model") <- vm
  attr(two, "nugget") <- 0.2
  two[] <- 450
  expect_identical(as.character(vm$model), c("Nug", "Exp", "Exp"))
  expect_identical(sac_nugget(two), 0.2)
  expect_equal(spatialkit:::.vgm_nugget_of(vm), 0.2)
  # The correlation function is the sum of both structures over the sill.
  cf <- spatialkit:::.vgm_correlation_fn(vm)
  expect_equal(cf(c(10, 100)), (0.5 * exp(-c(10, 100) / 20) + 0.3 * exp(-c(10, 100) / 150)) / 1)
  p <- suppressWarnings(resolution_profile(s$pts, "z", sac = two, n_levels = 4, nstart = 2))
  expect_identical(attr(p, "cp_noise")$source, "variogram nugget")
  expect_equal(attr(p, "cp_noise")$value, 0.2)
  expect_equal(attr(p, "variogram")$psill, 0.8)
  sm <- summarize_by_cell(s$asg, "z", deff = "variogram", sac = two)
  expect_identical(attr(sm, "deff_applied")$method, "variogram")
  ka <- kriging_adequacy(s$asg, "z", s$cells, sac = two, k = 3)
  expect_s3_class(ka, "kriging_adequacy")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(two), "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot(two)))
})

test_that("plot() draws each verdict and print() stays one line of number", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("ggplot2")
  s <- fvc_setup()
  ok <- estimate_sac_range(s$pts, "z")
  no <- suppressWarnings(estimate_sac_range(s$pts, "noise"))
  flat <- suppressWarnings(estimate_sac_range(s$pts, "noise", structure_alpha = 1))
  for (r in list(ok, no, flat)) {
    p <- plot(r)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
    out <- utils::capture.output(print(r))
    expect_match(out[1], "^(NA|[0-9])")
    expect_false(any(grepl("psill|gamma", out)))
  }
  expect_match(plot(ok)$labels$subtitle, "^Effective range = ")
  expect_match(gsub("\n", " ", plot(no)$labels$subtitle), NO_STRUCTURE_REASON, fixed = TRUE)
})
