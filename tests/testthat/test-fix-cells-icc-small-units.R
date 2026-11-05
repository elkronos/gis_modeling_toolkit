# tests/testthat/test-fix-cells-icc-small-units.R
# ---------------------------------------------------------------------------
# summarize_by_cell(deff = "kish"): the ICC of a column does not depend on
# its unit.
#
# .estimate_icc() skipped a column as one that does not vary when its
# standard deviation was below .Machine$double.eps, on all its values and
# again on the rows with a cell id.  A column that varies, in units that
# make its values of order 1e-16 or smaller, got an ICC of NA with "the
# column does not vary", no design effect and the uncorrected standard
# errors, with no warning; where only the rows with a cell id were that
# small, the ICC was NA with "its values differ only on rows that fall in no
# cell".  A standard deviation that overflowed gave NA with "the variance
# decomposition is undefined".  Where the standard deviation is below the
# floor or not finite, the column is now read again by .constant_column(),
# and a column with spread gets its ICC from .spread_scaled() of it.
# ---------------------------------------------------------------------------

# 60 cells of 12: an unclustered response and a clustered predictor.
.iu_points <- function() {
  set.seed(3)
  k <- 60L; per <- 12L; n <- k * per; cell <- rep(seq_len(k), each = per)
  d <- data.frame(x = 4e5 + runif(n, 0, 1e5), y = 4.6e6 + runif(n, 0, 1e5), poly_id = cell)
  d$resp  <- rnorm(n)
  d$clust <- rnorm(k)[cell] + rnorm(n, 0, 0.33)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
}
.iu_run <- function(expr) {
  w <- list()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val, classes = vapply(w, function(c) class(c)[1L], character(1)))
}

test_that("a column's ICC and standard errors over its scale do not depend on its unit", {
  pts <- .iu_points()
  ref <- .iu_run(summarize_by_cell(pts, "clust", "clust", deff = "kish", estimate = "mean"))$value
  rho <- attr(ref, "icc")$resp
  expect_gt(rho, 0.5)
  for (sc in c(1e-14, 1e-17, 1e-20, 1e-100)) {
    p <- pts; p$c <- p$clust * sc
    res <- .iu_run(summarize_by_cell(p, "c", "c", deff = "kish", estimate = "mean"))
    out <- res$value
    icc <- attr(out, "icc")
    expect_equal(icc$resp, rho, tolerance = 1e-12, info = format(sc))
    expect_equal(icc$pred[["c"]], rho, tolerance = 1e-12, info = format(sc))
    expect_length(icc$reason, 0L)
    expect_false("spatialkit_icc_unavailable" %in% res$classes)
    expect_false(is.null(attr(out, "deff_applied")))
    expect_true(all(out$deff_applied))
    for (s in c("resp", "pred"))
      expect_equal(out[[paste0("..se_", s, "_c")]] / sc, ref[[paste0("..se_", s, "_clust")]],
                   tolerance = 1e-10, info = paste(format(sc), s))
    expect_equal(out$cell_weight, ref$cell_weight, tolerance = 1e-10, info = format(sc))
  }
})

test_that("rows with a cell id in small units beside rows with no cell id get the ICC", {
  pts <- .iu_points()[, c("poly_id", "clust")]
  noid <- sf::st_as_sf(data.frame(x = 4e5, y = 4.6e6, poly_id = NA_integer_,
                                  clust = rep(c(0, 1), 5)), coords = c("x", "y"), crs = 32615)
  ref <- attr(.iu_run(summarize_by_cell(pts, predictor_vars = "clust", deff = "kish",
                                        estimate = "mean"))$value, "icc")$pred[["clust"]]
  p <- rbind(pts, noid)
  id <- !is.na(p$poly_id)
  p$clust[id] <- p$clust[id] * 1e-20
  res <- .iu_run(summarize_by_cell(p, predictor_vars = "clust", deff = "kish", estimate = "mean"))
  expect_equal(attr(res$value, "icc")$pred[["clust"]], ref, tolerance = 1e-12)
  expect_false("spatialkit_icc_unavailable" %in% res$classes)
  expect_true(all(is.finite(res$value[["..se_pred_clust"]])))
})

test_that("a column whose standard deviation overflows gets its ICC", {
  pts <- .iu_points()
  rho <- attr(.iu_run(summarize_by_cell(pts, predictor_vars = "clust", deff = "kish",
                                        estimate = "mean"))$value, "icc")$pred[["clust"]]
  p <- pts; p$c <- p$clust * 1e200
  expect_identical(stats::sd(p$c), Inf)
  res <- .iu_run(summarize_by_cell(p, predictor_vars = "c", deff = "kish", estimate = "mean"))
  expect_equal(attr(res$value, "icc")$pred[["c"]], rho, tolerance = 1e-12)
  expect_length(attr(res$value, "icc")$reason, 0L)
  expect_false("spatialkit_icc_unavailable" %in% res$classes)
})

test_that("a column of one value, or of values equal to their rounding, still does not vary", {
  pts <- .iu_points()
  cols <- list(seven = rep(7, nrow(pts)), tiny = rep(1e-20, nrow(pts)),
               rounding = replace(rep(0.1, nrow(pts)), c(1, 13, 25), 0.3 / 3))
  for (nm in names(cols)) {
    p <- pts; p$c <- cols[[nm]]
    res <- .iu_run(summarize_by_cell(p, "c", deff = "kish", estimate = "mean"))
    expect_identical(attr(res$value, "icc")$reason, c(c = "the column does not vary"), info = nm)
    expect_true(is.na(attr(res$value, "icc")$resp), info = nm)
    expect_null(attr(res$value, "deff_applied"))
    expect_length(res$classes, 0L)
    expect_identical(res$value$cell_weight, rep(12, 60), info = nm)
  }
  p <- pts; p$c <- 7
  expect_identical(.iu_run(summarize_by_cell(p, "c", deff = "kish", estimate = "mean"))$value[["..se_resp_c"]],
                   rep(0, 60))
})

# The standard deviation inside a cell squared the raw deviations, which
# overflow for values near 1e200 and underflow for a spread near 1e-170:
# the standard errors were NA, or 0, in every cell, with no warning, under
# any `deff`.  Where the raw standard deviation is not finite or below
# .POOL_Y_MIN it is now taken on .spread_scaled() of the cell's values.
test_that("the standard deviation and standard error inside a cell do not depend on the unit", {
  pts <- .iu_points()
  for (deff in list(1, "kish")) {
    ref <- .iu_run(summarize_by_cell(pts, "clust", "clust", deff = deff, estimate = "mean"))$value
    for (sc in c(1e200, 1e300, 1e-170, 1e-300)) {
      p <- pts; p$c <- p$clust * sc
      res <- .iu_run(summarize_by_cell(p, "c", "c", deff = deff, estimate = "mean"))
      out <- res$value
      what <- paste(format(sc), deff)
      expect_length(res$classes, 0L)
      for (s in c("resp", "pred")) {
        expect_true(all(out[[paste0("..se_", s, "_c")]] > 0), info = paste(what, s))
        expect_equal(out[[paste0("..sd_", s, "_c")]] / sc, ref[[paste0("..sd_", s, "_clust")]],
                     tolerance = 1e-10, info = paste(what, s))
        expect_equal(out[[paste0("..se_", s, "_c")]] / sc, ref[[paste0("..se_", s, "_clust")]],
                     tolerance = 1e-10, info = paste(what, s))
      }
      expect_equal(out$cell_weight, ref$cell_weight, tolerance = 1e-10, info = what)
    }
  }
})

test_that("the interval of a cell in large or small units is the interval of the column", {
  pts <- .iu_points()
  ref <- .iu_run(summarize_by_cell(pts, "clust", deff = "kish", estimate = "mean",
                                   conf_level = 0.9))$value
  for (sc in c(1e200, 1e-300)) {
    p <- pts; p$c <- p$clust * sc
    out <- .iu_run(summarize_by_cell(p, "c", deff = "kish", estimate = "mean",
                                     conf_level = 0.9))$value
    expect_equal(out[["..ci_lo_resp_c"]] / sc, ref[["..ci_lo_resp_clust"]], tolerance = 1e-10,
                 info = format(sc))
    expect_equal(out[["..ci_hi_resp_c"]] / sc, ref[["..ci_hi_resp_clust"]], tolerance = 1e-10,
                 info = format(sc))
  }
})

test_that("a cell whose standard deviation is beyond the largest double is named in a warning", {
  pts <- .iu_points()
  p <- pts
  p$c <- p$clust
  p$c[1:12] <- rep(c(-1.79e308, 1.79e308), 6)
  res <- .iu_run(summarize_by_cell(p, predictor_vars = "c", estimate = "mean"))
  out <- res$value
  expect_identical(res$classes, "spatialkit_sd_overflow")
  expect_identical(out[["..sd_pred_c"]][1], Inf)
  expect_true(is.na(out[["..se_pred_c"]][1]))
  expect_true(all(is.finite(out[["..se_pred_c"]][-1])))
  cond <- tryCatch(summarize_by_cell(p, predictor_vars = "c", estimate = "mean"),
                   spatialkit_sd_overflow = function(w) w)
  expect_identical(cond$columns, "c")
  expect_identical(cond$cells, 1L)
  expect_match(conditionMessage(cond), "'c' (1 cell)", fixed = TRUE)
  # A column whose standard deviation is finite in every cell raises none.
  expect_length(.iu_run(summarize_by_cell(pts, predictor_vars = "clust",
                                          estimate = "mean"))$classes, 0L)
})

# The Kish ICC keeps the raw arithmetic wherever the raw standard deviation
# is at or above .Machine$double.eps (the item's Decision): a column just
# above that floor, at an offset of 1, gets the ICC of the raw one-way
# ANOVA to the last bit.
test_that("a column just above the floor keeps the ICC of its raw values", {
  pts <- .iu_points()
  p <- pts; p$c <- 1 + p$clust * 2^-45
  expect_gt(stats::sd(p$c), .Machine$double.eps)
  expect_lt(stats::sd(p$c), 1e3 * .Machine$double.eps)
  rho <- attr(.iu_run(summarize_by_cell(p, predictor_vars = "c", deff = "kish",
                                        estimate = "mean"))$value, "icc")$pred[["c"]]
  # The arithmetic of .estimate_icc() on the raw values, step by step.
  vals <- (p$c - mean(p$c)) / stats::sd(p$c)
  grps <- factor(paste("c", as.character(p$poly_id), sep = "\r"))
  k <- nlevels(grps); N <- length(vals)
  ni <- tabulate(grps)
  gm <- tapply(vals, grps, mean)
  msb <- sum(ni * (gm - mean(vals))^2) / (k - 1)
  msw <- sum((vals - gm[grps])^2) / (N - k)
  n0 <- (N - sum(ni^2) / N) / (k - 1)
  expect_identical(rho, max(0, min(1, (msb - msw) / (msb + (n0 - 1) * msw))))
})
