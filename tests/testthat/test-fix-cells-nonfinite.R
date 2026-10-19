# tests/testthat/test-fix-cells-nonfinite.R
# ---------------------------------------------------------------------------
# summarize_by_cell(): an infinite value in a summarised column.
#
# It used to pass as a valid observation everywhere the function counts or
# averages: the cell's mean came back Inf (NaN when both signs met), its SE
# NA, while ..neff_ and cell_weight still counted the row; and under
# deff = "kish" one infinite row set the ICC to 0, so every cell's standard
# error came back uncorrected (0.169 against 1.365).  No condition was raised.
# Infinite values are treated as missing now, once, before anything reads the
# columns, with a classed warning that gives the count.
# ---------------------------------------------------------------------------

# Three cells of 30 with real between-cell differences.
.fn_points <- function(seed = 1) {
  set.seed(seed)
  id <- rep(1:3, each = 30L)
  cx <- c(0, 1000, 2000)[id] + runif(90, 0, 800)
  cy <- runif(90, 0, 800)
  y  <- 10 + c(-1, 0, 1.5)[id] + rnorm(90, sd = 1)
  sf::st_as_sf(data.frame(poly_id = id, y = y, x = cx, yy = cy),
               coords = c("x", "yy"), crs = 32615)
}
.fn_sac <- function() structure(
  300, class = c("sac_range", "numeric"),
  variogram_model = data.frame(model = c("Nug", "Exp"), psill = c(0.2, 0.8),
                               range = c(0, 200)),
  crs = sf::st_crs(32615))

.fn_run <- function(expr) {
  w <- character(0); cls <- list(); conds <- list()
  val <- withCallingHandlers(
    expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); cls[[length(cls) + 1L]] <<- class(c)
      conds[[length(conds) + 1L]] <<- c
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = w, classes = cls, conditions = conds)
}
# The frame and the attributes a comparison should cover.
.fn_same <- function(a, b) {
  expect_equal(as.data.frame(a), as.data.frame(b), tolerance = 1e-12)
  expect_equal(attr(a, "icc"), attr(b, "icc"), tolerance = 1e-12)
  expect_equal(attr(a, "deff_applied"), attr(b, "deff_applied"), tolerance = 1e-12)
}

test_that("one Inf in a cell of 30 is left out of every statistic, with one warning", {
  pts0 <- .fn_points()
  ptsI <- pts0;  ptsI$y[5] <- Inf
  ptsN <- pts0;  ptsN$y[5] <- NA

  res <- .fn_run(summarize_by_cell(ptsI, response_var = "y", conf_level = 0.95, estimate = "mean"))
  out <- res$value
  expect_equal(out$resp_mean_y[1], mean(pts0$y[1:30][-5]), tolerance = 1e-12)
  expect_identical(out$n[1], 30L)
  expect_equal(out[["..neff_resp_y"]][1], 29)
  expect_equal(out$cell_weight[1], 29)
  expect_true(all(is.finite(unlist(out[1, c("..sd_resp_y", "..se_resp_y",
                                            "..ci_lo_resp_y", "..ci_hi_resp_y")]))))
  # Exactly one warning, classed, naming 1 row and 1 cell.
  expect_length(res$warnings, 1L)
  expect_true("spatialkit_nonfinite" %in% res$classes[[1L]])
  expect_match(res$warnings, "1 infinite value\\(s\\) are treated as missing: response 'y' 1 row\\(s\\) in 1 cell\\(s\\)")
  expect_identical(res$conditions[[1L]]$n, 1L)
  expect_identical(res$conditions[[1L]]$columns, "y")
  # The result is the one with that row set to NA.
  .fn_same(out, summarize_by_cell(ptsN, response_var = "y", conf_level = 0.95, estimate = "mean"))
})

test_that("an Inf row gives the row-as-NA result under every design effect", {
  pts0 <- .fn_points()
  ptsI <- pts0;  ptsI$y[5] <- Inf
  ptsN <- pts0;  ptsN$y[5] <- NA

  # "kish": the ICC is the clean one, so every cell is corrected.  It was 0,
  # with deff_applied FALSE and uncorrected SEs on all three cells.
  kI <- suppressWarnings(summarize_by_cell(ptsI, "y", deff = "kish", conf_level = 0.95,
                                           estimate = "mean"))
  kN <- summarize_by_cell(ptsN, "y", deff = "kish", conf_level = 0.95, estimate = "mean")
  .fn_same(kI, kN)
  expect_gt(attr(kI, "icc")$resp, 0.5)
  expect_true(all(kI$deff_applied))
  expect_gt(kI[["..se_resp_y"]][2], 1)

  # "variogram" with a supplied sac: cell_weight equals ..neff_.
  vI <- suppressWarnings(summarize_by_cell(ptsI, "y", deff = "variogram",
                                           sac = .fn_sac(), conf_level = 0.95, estimate = "mean"))
  vN <- summarize_by_cell(ptsN, "y", deff = "variogram", sac = .fn_sac(),
                          conf_level = 0.95, estimate = "mean")
  .fn_same(vI, vN)
  expect_equal(vI$cell_weight, vI[["..neff_resp_y"]], tolerance = 1e-12)

  # A fixed number.
  .fn_same(suppressWarnings(summarize_by_cell(ptsI, "y", deff = 2, conf_level = 0.95,
                                              estimate = "mean")),
           summarize_by_cell(ptsN, "y", deff = 2, conf_level = 0.95, estimate = "mean"))
})

test_that("user aggregation functions see the same rows as the built-in statistics", {
  ptsI <- .fn_points();  ptsI$y[5] <- Inf
  # The documented mean, passed explicitly: finite, and inside its interval.
  # It was Inf beside a finite SE when only the built-ins were patched.
  out <- suppressWarnings(summarize_by_cell(
    ptsI, "y", agg_funs = list(mean = function(x) mean(x, na.rm = TRUE)),
    conf_level = 0.95))
  expect_true(is.finite(out$resp_mean_y[1]))
  expect_gt(out$resp_mean_y[1], out[["..ci_lo_resp_y"]][1])
  expect_lt(out$resp_mean_y[1], out[["..ci_hi_resp_y"]][1])
  # max and sum never see the infinite value.
  mx <- suppressWarnings(summarize_by_cell(
    ptsI, "y", agg_funs = list(max = function(x) max(x, na.rm = TRUE),
                               sum = function(x) sum(x, na.rm = TRUE))))
  expect_equal(mx$resp_max_y[1], max(.fn_points()$y[1:30][-5]))
  expect_true(is.finite(mx$resp_sum_y[1]))
})

test_that("Inf and -Inf in one cell, an Inf in a predictor, and several columns", {
  pts0 <- .fn_points()
  # Both signs in one cell used to give a mean of NaN.
  ptsB <- pts0; ptsB$y[5] <- Inf; ptsB$y[6] <- -Inf
  rB <- .fn_run(summarize_by_cell(ptsB, "y", estimate = "mean"))
  expect_equal(rB$value$resp_mean_y[1], mean(pts0$y[1:30][-(5:6)]), tolerance = 1e-12)
  expect_match(rB$warnings, "2 infinite value\\(s\\).*response 'y' 2 row\\(s\\) in 1 cell\\(s\\)")
  expect_equal(rB$value$cell_weight[1], 28)

  # A predictor column behaves the same way, and the response is untouched.
  set.seed(2)
  ptsD <- pts0; ptsD$p <- rnorm(90); p_clean <- ptsD$p; ptsD$p[70] <- Inf
  rD <- .fn_run(summarize_by_cell(ptsD, "y", "p", conf_level = 0.9, estimate = "mean"))
  expect_equal(rD$value$pred_mean_p[3], mean(p_clean[61:90][-10]), tolerance = 1e-12)
  expect_equal(rD$value[["..neff_pred_p"]][3], 29)
  expect_true(is.finite(rD$value[["..se_pred_p"]][3]))
  expect_match(rD$warnings, "predictor 'p' 1 row\\(s\\) in 1 cell\\(s\\)")
  expect_equal(rD$value$resp_mean_y, summarize_by_cell(pts0, "y", estimate = "mean")$resp_mean_y)
  expect_equal(rD$value$cell_weight, rep(30, 3))

  # Two columns, three cells: still one warning, naming both.
  ptsE <- ptsD; ptsE$y[c(5, 40)] <- c(Inf, -Inf)
  rE <- .fn_run(summarize_by_cell(ptsE, "y", "p", estimate = "mean"))
  expect_length(rE$warnings, 1L)
  expect_match(rE$warnings, "3 infinite value\\(s\\) are treated as missing: response 'y' 2 row\\(s\\) in 2 cell\\(s\\); predictor 'p' 1 row\\(s\\) in 1 cell\\(s\\)")
  expect_identical(rE$conditions[[1L]]$columns, c("y", "p"))
  expect_identical(rE$conditions[[1L]]$n, 3L)
})

test_that("data with no infinite value raise nothing and give the same frame", {
  pts0 <- .fn_points()
  ptsN <- pts0; ptsN$y[c(5, 40)] <- c(NA, NaN)       # missing already, as before
  for (p in list(pts0, ptsN)) {
    res <- .fn_run(summarize_by_cell(p, "y", deff = "kish", conf_level = 0.95, estimate = "mean"))
    expect_length(res$warnings, 0L)
  }
  out <- summarize_by_cell(pts0, "y", conf_level = 0.95, estimate = "mean")
  d <- sf::st_drop_geometry(pts0)
  expect_equal(out$resp_mean_y, as.numeric(tapply(d$y, d$poly_id, mean)), tolerance = 1e-12)
  expect_equal(out[["..se_resp_y"]],
               as.numeric(tapply(d$y, d$poly_id, stats::sd)) / sqrt(30), tolerance = 1e-12)
  expect_equal(out$cell_weight, rep(30, 3))
})
