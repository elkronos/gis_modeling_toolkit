# tests/testthat/test-fix-cells-deff-unidentified.R
# ---------------------------------------------------------------------------
# summarize_by_cell(): two design effects the data do not identify.
#
# 1. deff = "variogram" with a variogram that finds no spatial structure.
#    Its pure nugget gives every column a design effect of 1, which is the
#    answer for the variable it was fitted to and says nothing about the
#    predictor columns summarised beside it.  They kept the uncorrected
#    standard errors in silence; a warning of class
#    spatialkit_deff_no_structure now names them.
# 2. deff = "kish" (and a "variogram" request that falls back to it) with an
#    ICC that cannot be estimated: values in one cell, or fewer than two
#    cells with two or more values.  .estimate_icc() returned 0 for that and
#    for an estimate at or below 0 alike, so both gave the uncorrected
#    standard error under one label.  The first is NA now, with the reason on
#    attr(, "icc")$reason and a warning of class spatialkit_icc_unavailable.
# ---------------------------------------------------------------------------

.du_run <- function(expr) {
  w <- list()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(c) {
    w[[length(w) + 1L]] <<- c; invokeRestart("muffleWarning")
  })
  list(value = val, conds = w,
       classes = vapply(w, function(c) class(c)[1L], character(1)))
}
.du_of <- function(res, class) Filter(function(c) inherits(c, class), res$conds)

# 700 points on a 100 km square: a white-noise response, a predictor with an
# exponential effective range of 30 km and a white-noise predictor; 30
# k-means Voronoi cells.
.du_layer <- local({
  val <- NULL
  function() {
    if (is.null(val)) {
      set.seed(12); n <- 700
      x <- runif(n, 0, 1e5); y <- runif(n, 0, 1e5)
      D <- as.matrix(stats::dist(cbind(x, y)))
      x1 <- as.numeric(t(chol(exp(-3 * D / 30000) + diag(1e-8, n))) %*% rnorm(n))
      pts <- sf::st_as_sf(data.frame(x = x, y = y, z = rnorm(n), x1 = x1, x2 = rnorm(n)),
                          coords = c("x", "y"), crs = 32617)
      seeds <- suppressWarnings(voronoi_seeds_kmeans(pts, k = 30, partition = "kmeans"))
      tess  <- suppressWarnings(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
      asg   <- suppressWarnings(assign_features_to_polygons(pts, tess$cells,
                                                            polygon_id_col = "cell_id"))
      val <<- list(pts = pts, asg = asg)
    }
    val
  }
})

test_that("predictors summarised under a no-structure verdict are named in a classed warning", {
  skip_if_not_installed("gstat")
  L <- .du_layer()
  no <- suppressWarnings(estimate_sac_range(L$pts, "z"))
  expect_true(.sac_no_structure(no))
  iid <- summarize_by_cell(L$asg, "z", c("x1", "x2"), id_col = "cell_id", estimate = "mean")
  for (sac in list(NULL, no)) {
    res <- .du_run(summarize_by_cell(L$asg, "z", c("x1", "x2"), id_col = "cell_id",
                                     deff = "variogram", sac = sac, estimate = "mean"))
    hit <- .du_of(res, "spatialkit_deff_no_structure")
    expect_length(hit, 1L)
    expect_identical(hit[[1L]]$columns, c("x1", "x2"))
    expect_identical(hit[[1L]]$reason, attr(no, "rejected_reason"))
    msg <- conditionMessage(hit[[1L]])
    expect_match(msg, "every column gets a design effect of 1, the predictor columns 'x1' and 'x2' included",
                 fixed = TRUE)
    expect_match(msg, if (is.null(sac)) "estimated from `response_var`" else "of the supplied `sac`",
                 fixed = TRUE)
    # No conf_level, so no effective n or intervals to speak of.
    expect_match(msg, "the standard errors of those columns are the uncorrected ones", fixed = TRUE)
    expect_false(grepl("effective n and intervals", msg, fixed = TRUE))
    expect_length(.du_of(res, "spatialkit_deff_fallback"), 0L)
    # The warning changes no number: the design effect is still 1.
    expect_identical(unique(res$value$deff_method), "variogram")
    for (col in c("..se_resp_z", "..se_pred_x1", "..se_pred_x2"))
      expect_equal(res$value[[col]], iid[[col]], tolerance = 1e-12, info = col)
  }
  # With predictors only and the response's verdict supplied, the same.
  res <- .du_run(summarize_by_cell(L$asg, predictor_vars = "x1", id_col = "cell_id",
                                   deff = "variogram", sac = no, conf_level = 0.95,
                                   estimate = "mean"))
  hit <- .du_of(res, "spatialkit_deff_no_structure")
  expect_length(hit, 1L)
  expect_identical(hit[[1L]]$columns, "x1")
  expect_match(conditionMessage(hit[[1L]]),
               "the predictor column 'x1' gets a design effect of 1", fixed = TRUE)
  expect_match(conditionMessage(hit[[1L]]),
               "the standard errors, effective n and intervals of that column", fixed = TRUE)
  # The response named again among the predictors is not one of them here.
  res <- .du_run(summarize_by_cell(L$asg, "z", c("z", "x1"), id_col = "cell_id",
                                   deff = "variogram", sac = no, estimate = "mean"))
  hit <- .du_of(res, "spatialkit_deff_no_structure")
  expect_length(hit, 1L)
  expect_identical(hit[[1L]]$columns, "x1")
  expect_match(conditionMessage(hit[[1L]]),
               "every column gets a design effect of 1, the predictor column 'x1' included",
               fixed = TRUE)
})

test_that("no such warning without predictors, under a structured variogram, or for the response named twice", {
  skip_if_not_installed("gstat")
  L <- .du_layer()
  no <- suppressWarnings(estimate_sac_range(L$pts, "z"))
  alone <- .du_run(summarize_by_cell(L$asg, "z", id_col = "cell_id", deff = "variogram",
                                     sac = no, estimate = "mean"))
  expect_length(alone$conds, 0L)
  twice <- .du_run(summarize_by_cell(L$asg, "z", "z", id_col = "cell_id", deff = "variogram",
                                     sac = no, estimate = "mean"))
  expect_length(.du_of(twice, "spatialkit_deff_no_structure"), 0L)
  fitted <- suppressWarnings(estimate_sac_range(L$pts, "x1"))
  expect_true(is.finite(as.numeric(fitted)))
  ok <- .du_run(summarize_by_cell(L$asg, "x1", "x2", id_col = "cell_id", deff = "variogram",
                                  sac = fitted, estimate = "mean"))
  expect_length(.du_of(ok, "spatialkit_deff_no_structure"), 0L)
  # A pure nugget the caller built, with no verdict behind it, is applied as
  # given (test-review3-S3-assignment.R), and names nothing.
  nug <- structure(1, class = c("sac_range", "numeric"),
                   variogram_model = data.frame(model = "Nug", psill = 1, range = 0),
                   crs = sf::st_crs(32617))
  own <- .du_run(summarize_by_cell(L$asg, "z", "x1", id_col = "cell_id", deff = "variogram",
                                   sac = nug, estimate = "mean"))
  expect_length(.du_of(own, "spatialkit_deff_no_structure"), 0L)
})

# Layers on which the ICC is estimated, and on which it cannot be.
.du_mk <- function(d) sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
.du_flat <- function() {
  set.seed(5); n <- 150
  .du_mk(data.frame(x = runif(n, 0, 600), y = runif(n, 0, 60),
                    z = rep(c(-1, 1), length.out = n) + rnorm(n, 0, 0.01),
                    poly_id = rep(1:6, each = 25)))
}
.du_one_cell <- function() {
  set.seed(6); n <- 150
  .du_mk(data.frame(x = runif(n, 0, 600), y = runif(n, 0, 60),
                    z = rep(rnorm(6, sd = 3), each = 25) + rnorm(n), poly_id = 1L))
}
.du_lone <- function() {
  set.seed(7)
  .du_mk(data.frame(x = runif(30, 0, 600), y = runif(30, 0, 60),
                    z = c(rnorm(25, 5), rnorm(5)), poly_id = c(rep(1L, 25), 2:6)))
}
# .du_flat() with a clustered column and one that has values in one cell.
.du_mixed <- function() {
  pts <- .du_flat()
  set.seed(9)
  pts$clust <- rep(rnorm(6, sd = 3), each = 25) + rnorm(150)
  pts$few <- NA_real_
  pts$few[pts$poly_id == 2L] <- rnorm(25)
  pts
}

test_that("an ICC that cannot be estimated gives NA standard errors, with the reason", {
  for (case in list(list(.du_one_cell(), "its values fall in one cell"),
                    list(.du_lone(), "fewer than two cells hold two or more of its values"))) {
    pts <- case[[1L]]
    res <- .du_run(summarize_by_cell(pts, "z", deff = "kish", conf_level = 0.95,
                                     estimate = "mean"))
    out <- res$value
    hit <- .du_of(res, "spatialkit_icc_unavailable")
    expect_length(hit, 1L)
    expect_identical(hit[[1L]]$columns, "z")
    expect_identical(unname(hit[[1L]]$reason), case[[2L]])
    expect_match(conditionMessage(hit[[1L]]), sprintf("'z' (%s)", case[[2L]]), fixed = TRUE)
    expect_match(conditionMessage(hit[[1L]]), "and so is cell_weight", fixed = TRUE)
    icc <- attr(out, "icc")
    expect_identical(names(icc), c("resp", "pred", "reason"))
    expect_true(is.na(icc$resp))
    expect_identical(icc$reason, c(z = case[[2L]]))
    big <- out$n >= 2L
    for (col in c("..se_resp_z", "..neff_resp_z", "..ci_lo_resp_z", "..ci_hi_resp_z"))
      expect_true(all(is.na(out[[col]])), info = col)
    expect_true(all(is.na(out$cell_weight[big])))
    expect_identical(out$cell_weight[!big], rep(1, sum(!big)))
    # The degrees of freedom of s^2 do not depend on the ICC.
    expect_identical(out[["..df_resp_z"]][big], as.numeric(out$n[big] - 1L))
    # The means are untouched.
    expect_equal(out$resp_mean_z[big],
                 as.numeric(tapply(pts$z, pts$poly_id, mean))[big], tolerance = 1e-12)
    expect_identical(unique(out$deff_method), "none")
    expect_false(any(out$deff_applied))
    expect_null(attr(out, "deff_applied"))
  }
})

test_that("an estimate at or below 0 is still 0, with the uncorrected standard errors", {
  pts <- .du_flat()
  res <- .du_run(summarize_by_cell(pts, "z", deff = "kish", estimate = "mean"))
  expect_length(.du_of(res, "spatialkit_icc_unavailable"), 0L)
  expect_identical(attr(res$value, "icc"), list(resp = 0, pred = NA_real_))
  iid <- summarize_by_cell(pts, "z", estimate = "mean")
  expect_equal(res$value[["..se_resp_z"]], iid[["..se_resp_z"]], tolerance = 1e-12)
  expect_identical(res$value$cell_weight, rep(25, 6))
})

test_that("one column's missing ICC leaves the others as they were", {
  pts <- .du_mixed()
  both <- .du_run(summarize_by_cell(pts, "clust", c("few", "z"), deff = "kish",
                                    conf_level = 0.95, estimate = "mean"))
  alone <- summarize_by_cell(pts, "clust", "z", deff = "kish", conf_level = 0.95,
                             estimate = "mean")
  hit <- .du_of(both, "spatialkit_icc_unavailable")
  expect_length(hit, 1L)
  expect_identical(hit[[1L]]$columns, "few")
  # 'few' does not set cell_weight, so the warning does not say it is NA.
  expect_false(grepl("cell_weight", conditionMessage(hit[[1L]]), fixed = TRUE))
  for (col in c("..se_resp_clust", "..se_pred_z", "..ci_lo_pred_z", "cell_weight"))
    expect_equal(both$value[[col]], alone[[col]], tolerance = 1e-12, info = col)
  expect_true(all(is.na(both$value[["..se_pred_few"]])))
  icc <- attr(both$value, "icc")
  expect_true(is.na(icc$pred[["few"]]))
  expect_identical(icc$reason, c(few = "its values fall in one cell"))
  expect_identical(unique(both$value$deff_method), "kish")
  applied <- attr(both$value, "deff_applied")
  expect_true(is.na(applied$icc_pred[["few"]]))
  expect_identical(applied$deff, attr(alone, "deff_applied")$deff)
})

test_that("the primary variable's missing ICC leaves its design effect NA beside a predictor's positive one", {
  # The response has ten values in cell 1, one in cell 2 and none elsewhere;
  # the predictor is clustered in all six cells.
  set.seed(70)
  d <- data.frame(x = runif(60, 0, 600), y = runif(60, 0, 60), poly_id = rep(1:6, each = 10))
  d$z <- NA_real_
  d$z[d$poly_id == 1L] <- rnorm(10)
  d$z[d$poly_id == 2L][1L] <- 0.3
  d$x1 <- rep(rnorm(6, sd = 3), each = 10) + rnorm(60)
  res <- .du_run(summarize_by_cell(.du_mk(d), "z", "x1", deff = "kish", estimate = "mean"))
  out <- res$value
  expect_identical(.du_of(res, "spatialkit_icc_unavailable")[[1L]]$columns, "z")
  icc <- attr(out, "icc")
  expect_true(is.na(icc$resp))
  expect_gt(icc$pred[["x1"]], 0)
  # The predictor's ICC is applied, so the attribute is there, and its
  # `deff` is the response's: unknown in the cell of ten, 1 elsewhere.
  applied <- attr(out, "deff_applied")
  expect_identical(applied$method, "kish")
  expect_true(is.na(applied$icc_resp))
  expect_identical(applied$deff, c(NA, 1, 1, 1, 1, 1))
  expect_identical(out$cell_weight, c(NA, 1, 0, 0, 0, 0))
  expect_true(all(out$deff_applied))
  expect_identical(unique(out$deff_method), "kish")
  expect_true(all(is.finite(out[["..se_pred_x1"]])))
})

test_that("rows with no cell id give no cell to the ICC, and are not read as a column that does not vary", {
  na_id <- function(z, id = c(1:6, rep(NA, 10))) {
    force(z)
    set.seed(31)
    .du_mk(data.frame(x = runif(length(id), 0, 600), y = runif(length(id), 0, 60),
                      z = z, poly_id = id))
  }
  check <- function(pts, reason, warned = TRUE) {
    res <- .du_run(summarize_by_cell(pts, "z", deff = "kish", estimate = "mean"))
    out <- res$value
    expect_identical(attr(out, "icc")$reason, c(z = reason))
    hit <- .du_of(res, "spatialkit_icc_unavailable")
    expect_length(hit, as.integer(warned))
    out
  }
  set.seed(30)
  # One value in each of six cells, and ten in rows with no cell id: the
  # only summary row of two or more values is theirs.
  out <- check(na_id(rnorm(16)), "fewer than two cells hold two or more of its values")
  noid <- is.na(out$poly_id)
  expect_identical(out$n[noid], 10L)
  expect_true(is.na(out[["..se_resp_z"]][noid]))
  expect_true(is.na(out$cell_weight[noid]))
  expect_identical(out$cell_weight[!noid], rep(1, 6))
  # Equal on the rows with a cell id, varying on the others.
  id <- c(rep(1:6, each = 5), rep(NA, 10))
  out <- check(na_id(c(rep(2, 30), rnorm(10)), id),
               "its values differ only on rows that fall in no cell")
  expect_true(all(is.na(out[["..se_resp_z"]])))
  expect_true(all(is.na(out$cell_weight)))
  # Every value in rows with no cell id.
  check(na_id(c(rep(NA, 30), rnorm(10)), id), "its values fall in no cell")
  # One value in a cell, the others in rows with no cell id.
  check(na_id(c(0.5, rep(NA, 29), rnorm(10)), id), "only one of its values falls in a cell")
  # Equal everywhere, those rows included: standard errors of 0, no warning.
  out <- check(na_id(rep(2, 40), id), "the column does not vary", warned = FALSE)
  expect_identical(out[["..se_resp_z"]], rep(0, 7))
  expect_identical(out$cell_weight, c(rep(5, 6), 10))
})

test_that("nothing is withheld where nothing could be computed, or where the column does not vary", {
  set.seed(8)
  singles <- .du_mk(data.frame(x = seq(5, 795, by = 100), y = 50, z = rnorm(8),
                               poly_id = seq_len(8)))
  res <- .du_run(summarize_by_cell(singles, "z", deff = "kish", estimate = "mean"))
  expect_length(.du_of(res, "spatialkit_icc_unavailable"), 0L)
  expect_identical(attr(res$value, "icc")$reason,
                   c(z = "fewer than two cells hold two or more of its values"))
  expect_identical(res$value$cell_weight, rep(1, 8))
  # One value in all: not a column that does not vary.
  one <- .du_flat()
  one$w <- NA_real_
  one$w[1L] <- 4
  res <- .du_run(summarize_by_cell(one, "z", "w", deff = "kish", estimate = "mean"))
  expect_length(.du_of(res, "spatialkit_icc_unavailable"), 0L)
  expect_identical(attr(res$value, "icc")$reason, c(w = "it has fewer than four values"))
  const <- .du_mk(data.frame(x = runif(150, 0, 600), y = runif(150, 0, 60), z = 3,
                             poly_id = rep(1:6, each = 25)))
  res <- .du_run(summarize_by_cell(const, "z", deff = "kish", conf_level = 0.95,
                                   estimate = "mean"))
  expect_length(.du_of(res, "spatialkit_icc_unavailable"), 0L)
  expect_identical(attr(res$value, "icc")$reason, c(z = "the column does not vary"))
  expect_true(is.na(attr(res$value, "icc")$resp))
  # A mean of equal values has no error whatever the correlation.
  expect_identical(res$value[["..se_resp_z"]], rep(0, 6))
  expect_identical(res$value[["..neff_resp_z"]], rep(25, 6))
  expect_identical(res$value$cell_weight, rep(25, 6))
})

test_that("a variogram request that falls back to Kish says which ICCs cannot be estimated", {
  # Predictors only: there is no response to fit a variogram to.
  pts <- .du_one_cell()
  res <- .du_run(summarize_by_cell(pts, predictor_vars = "z", deff = "variogram",
                                   estimate = "mean"))
  fb <- .du_of(res, "spatialkit_deff_fallback")
  expect_length(fb, 1L)
  expect_match(conditionMessage(fb[[1L]]),
               "whose ICC is 0 or cannot be estimated (predictor 'z' NA)", fixed = TRUE)
  # No column keeps an uncorrected standard error, so the text does not
  # speak of one.
  expect_match(conditionMessage(fb[[1L]]),
               "so no column carries a correction: the standard errors are NA for 'z', whose ICC cannot be estimated;",
               fixed = TRUE)
  expect_false(grepl("uncorrected", conditionMessage(fb[[1L]]), fixed = TRUE))
  expect_identical(fb[[1L]]$used, "none")
  icc_na <- .du_of(res, "spatialkit_icc_unavailable")
  expect_length(icc_na, 1L)
  # No conf_level, so no effective n or intervals to speak of.
  expect_match(conditionMessage(icc_na[[1L]]), "its standard errors are NA in every cell", fixed = TRUE)
  expect_false(grepl("effective n and intervals", conditionMessage(icc_na[[1L]]), fixed = TRUE))
  expect_true(all(is.na(res$value[["..se_pred_z"]])))
  expect_identical(unique(res$value$deff_method), "none (variogram fallback)")

  # Beside a column whose ICC is 0, which keeps its uncorrected standard
  # errors, the text names both.
  pts <- .du_mixed()
  res <- .du_run(summarize_by_cell(pts, predictor_vars = c("z", "few"), deff = "variogram",
                                   estimate = "mean"))
  fb <- .du_of(res, "spatialkit_deff_fallback")
  expect_length(fb, 1L)
  expect_identical(fb[[1L]]$used, "none")
  expect_match(conditionMessage(fb[[1L]]), paste0(
    "so no column carries a correction: the standard errors are the uncorrected ones, ",
    "and NA for 'few', whose ICC cannot be estimated;"), fixed = TRUE)
  iid <- summarize_by_cell(pts, predictor_vars = "z", estimate = "mean")
  expect_equal(res$value[["..se_pred_z"]], iid[["..se_pred_z"]], tolerance = 1e-12)
  expect_true(all(is.na(res$value[["..se_pred_few"]][res$value$poly_id == 2L])))

  # A positive ICC beside one that cannot be estimated: Kish's correction,
  # and the text names the column it could not correct.
  pts <- .du_mixed()
  res <- .du_run(summarize_by_cell(pts, predictor_vars = c("clust", "few"),
                                   deff = "variogram", estimate = "mean"))
  kish <- suppressWarnings(summarize_by_cell(pts, predictor_vars = "clust", deff = "kish",
                                             estimate = "mean"))
  fb <- .du_of(res, "spatialkit_deff_fallback")
  expect_length(fb, 1L)
  expect_identical(fb[[1L]]$used, "kish")
  expect_match(conditionMessage(fb[[1L]]), paste0(
    "the standard errors carry Kish's correction in place of the variogram's ",
    "(NA for 'few', whose ICC cannot be estimated), and `deff_method` reads"), fixed = TRUE)
  expect_identical(.du_of(res, "spatialkit_icc_unavailable")[[1L]]$columns, "few")
  expect_equal(res$value[["..se_pred_clust"]], kish[["..se_pred_clust"]], tolerance = 1e-12)
  expect_true(all(is.na(res$value[["..se_pred_few"]])))
  expect_identical(unique(res$value$deff_method), "kish (variogram fallback)")

  # A column that does not vary: its ICC is NA and its standard errors are 0,
  # so the text withholds nothing and no second warning is raised.
  const <- .du_mk(data.frame(x = runif(150, 0, 600), y = runif(150, 0, 60), z = 3,
                             poly_id = rep(1:6, each = 25)))
  res <- .du_run(summarize_by_cell(const, predictor_vars = "z", deff = "variogram",
                                   estimate = "mean"))
  fb <- .du_of(res, "spatialkit_deff_fallback")
  expect_length(fb, 1L)
  expect_match(conditionMessage(fb[[1L]]), paste0(
    "whose ICC is 0 or cannot be estimated (predictor 'z' NA), so no column carries ",
    "a correction: the standard errors are the uncorrected ones; `deff_method`"), fixed = TRUE)
  expect_identical(fb[[1L]]$used, "none")
  expect_length(.du_of(res, "spatialkit_icc_unavailable"), 0L)
  expect_identical(res$value[["..se_pred_z"]], rep(0, 6))
})

test_that("the pooled path takes an ICC that cannot be estimated, and its table says so", {
  own <- paste0("the cell's own mean, when its observations are spread through the cell; ",
                "not a range for one observation")
  pts <- .du_lone()
  res <- .du_run(summarize_by_cell(pts, "z", deff = "kish", conf_level = 0.95,
                                   estimate = "pooled"))
  expect_true("..pooled_resp_z" %in% names(res$value))
  expect_true(all(is.na(res$value[["..ci_lo_resp_z"]])))
  # Nothing was applied, and the column's intervals are NA: the row says
  # both.
  iv <- attr(res$value, "intervals")
  expect_identical(iv$method[iv$centre == "plain mean"],
                   "t, deff = 1; NA for 'z', whose ICC cannot be estimated")
  expect_identical(iv$covers[iv$centre == "plain mean"], own)
  shown <- utils::capture.output(print(attr(res$value, "pooled")))
  expect_identical(
    grep("..ci_lo_", shown, fixed = TRUE, value = TRUE),
    paste0("  ..ci_lo_, ..ci_hi_ : 95% for THE CELL'S MEAN (t, deff = 1; NA for 'z', ",
           "whose ICC cannot be estimated), centred on the plain mean;"))
  # A column that does not vary first, so nothing is applied, and one whose
  # ICC cannot be estimated second.
  pts <- .du_mixed()
  pts$k <- 3
  res <- .du_run(summarize_by_cell(pts, predictor_vars = c("k", "few"), deff = "kish",
                                   conf_level = 0.9, estimate = "pooled"))
  expect_null(attr(res$value, "deff_applied"))
  iv <- attr(res$value, "intervals")
  expect_identical(iv$method[iv$centre == "plain mean"],
                   "t, deff = 1; NA for 'few', whose ICC cannot be estimated")
  expect_identical(iv$covers[iv$centre == "plain mean"], own)
  # One ICC applied, one that cannot be estimated.
  pts <- .du_mixed()
  res <- .du_run(summarize_by_cell(pts, "clust", c("few", "z"), deff = "kish",
                                   conf_level = 0.95, estimate = "pooled"))
  iv <- attr(res$value, "intervals")
  expect_identical(iv$method[iv$centre == "plain mean"],
                   "t, deff = \"kish\"; 1 for 'z', whose ICC is 0; NA for 'few', whose ICC cannot be estimated")
})

test_that(".se_with_deff() returns NA for an unknown design effect, and .deff_label() names such columns", {
  expect_identical(.se_with_deff(2, 5, NA_real_), NA_real_)
  expect_equal(.se_with_deff(2, 5, 1), 2 / sqrt(5))
  applied <- list(method = "kish")
  expect_identical(.deff_label(applied), "\"kish\"")
  expect_identical(.deff_label(applied, kish_na = c("a", "b")),
                   "\"kish\"; NA for 'a' and 'b', whose ICC cannot be estimated")
  # Nothing applied, and a column left NA: the label still reads 1 first,
  # and the table still describes the default's interval.
  expect_identical(.deff_label(NULL, kish_na = "a"),
                   "1; NA for 'a', whose ICC cannot be estimated")
  expect_identical(.deff_label(NULL, kish_zero = "u", kish_na = "a"),
                   "1; NA for 'a', whose ICC cannot be estimated")
  own <- .intervals_table(NULL, 0.95, "1")$covers
  expect_identical(.intervals_table(NULL, 0.95, .deff_label(NULL, kish_na = "a"))$covers, own)
  expect_false(identical(.intervals_table(NULL, 0.95, "1.5")$covers, own))
  expect_false(identical(.intervals_table(NULL, 0.95, "10")$covers, own))
})

test_that("the help says where the NA falls", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "summarize_by_cell.Rd")
  # The Kish attribute's `deff` when the primary variable's ICC is unknown.
  expect_match(rd, paste0("all 1 when only a predictor's ICC was positive, except that it is ",
                          "\\code{NA} in each cell of two or more of the primary variable's ",
                          "values when that variable's ICC cannot be estimated"), fixed = TRUE)
  # `reason` is added when an ICC is NA; it is not a condition of the NA.
  expect_match(rd, paste0("\\code{NA} for a variable type with no numeric column and for a ",
                          "column whose ICC cannot be estimated; \\code{reason}, a character ",
                          "vector named by column, is then added and says why"), fixed = TRUE)
  expect_false(grepl("when \\code{reason}, a character vector named by column, says why",
                     rd, fixed = TRUE))
  # Rows with no cell id.
  expect_match(rd, paste0("give a summary row of their own but no cell to the ICC, so a column ",
                          "whose values all fall in them, or differ only there, has an ICC ",
                          "that cannot be estimated"), fixed = TRUE)
})
