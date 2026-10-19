# tests/testthat/test-core-pool-edge.R
# ---------------------------------------------------------------------------
# summarize_by_cell(estimate = "pooled", obs_level = ) on the inputs where
# pooling has little or nothing to work with: one test for each case the
# design lists for the pooled value and for the interval for one observation,
# through the exported function, with the values returned and the class and
# fields of the condition raised.
#
# The cells here are an ID column on the rows: the pooled columns read each
# row's ID and its coordinates and nothing of a polygon, so a layer can be
# given exactly the cells a case needs.  Every call names `estimate`.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.pe_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.pe_points <- function(x, y, data, crs = 32617)
  sf::st_as_sf(cbind(data.frame(x = x, y = y), data), coords = c("x", "y"), crs = crs)

.pe_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.pe_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.pe_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)
.pe_fields <- function(w) unclass(w)[setdiff(names(unclass(w)), c("message", "call"))]
.pe_new <- function(x) grep("^\\.\\.(pooled|obs)_", names(x), value = TRUE)

.pe_without_fnn <- function(env = parent.frame())
  testthat::local_mocked_bindings(
    .has_backend = function(pkg)
      if (identical(pkg, "FNN")) FALSE else requireNamespace(pkg, quietly = TRUE),
    .env = env)

.pe_pooled <- function(layer, ..., id_col = "cell_id")
  summarize_by_cell(layer, ..., id_col = id_col, estimate = "pooled")

# 720 rows of a smooth field with noise on a 10 x 6 grid of cells, the cell an
# arithmetic function of the coordinates: x, then y, then z.  Every cell
# holds several rows, and the layer raises nothing on the pooled path.
.pe_store <- new.env(parent = emptyenv())
.pe_layer <- function() {
  if (is.null(.pe_store$layer)) {
    d <- local({
      .pe_seed(3001)
      x <- stats::runif(720, 0, 1e5); y <- stats::runif(720, 0, 1e5)
      data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(720, 0, 0.5))
    })
    d$cell_id <- as.integer(6 * floor(d$x / 1e4) + floor(d$y / (1e5 / 6)) + 1)
    .pe_store$layer <- .pe_points(d$x, d$y, d[c("z", "cell_id")])
  }
  .pe_store$layer
}

# The pooled residuals of a column written out from what the frame and the
# record hold: each row predicted by its cell's value with the row taken out,
# sigma2, tau2 and the targets held.
.pe_residuals <- function(y, id, out, rec, name, id_col = "cell_id") {
  at <- match(id, out[[id_col]])
  n  <- as.integer(table(id)[as.character(id)])
  S  <- as.numeric(tapply(y, id, sum)[as.character(id)])
  t  <- out[[paste0("..pooled_target_", name)]][at]
  m1 <- (S - y) / pmax(n - 1L, 1L)
  w1 <- ifelse(n > 1L, rec$tau2 / (rec$tau2 + rec$sigma2 / pmax(n - 1L, 1L)), 0)
  y - (w1 * m1 + (1 - w1) * t)
}

# For each region the regions whose centroids are nearest, by the rule of the
# design: distance, then the lower x, then the lower y, then the lower number.
.pe_neighbours <- function(cen, k = 8L) {
  K <- nrow(cen)
  t(vapply(seq_len(K), function(j) {
    d2 <- (cen[, 1L] - cen[j, 1L])^2 + (cen[, 2L] - cen[j, 2L])^2
    o  <- order(d2, cen[, 1L], cen[, 2L], seq_len(K))
    o[o != j][seq_len(min(k, K - 1L))]
  }, integer(min(k, K - 1L))))
}

test_that("the layer of these tests has cells of several rows and raises nothing", {
  L <- .pe_layer()
  expect_identical(length(unique(L$cell_id)), 60L)
  expect_gte(min(table(L$cell_id)), 4L)
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  expect_identical(attr(r$value, "pooled")$columns$resp_z$status, "ok")
})


# ---- A cell with one row, and a cell with none -----------------------------------

test_that("a cell with one row has a value between its row and its target, and the same interval", {
  L <- .pe_layer()
  L$cell_id[7L] <- 999L                                # row 7 is now a cell of its own
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9, conf_level = 0.95))
  expect_length(r$warnings, 0L)
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_z
  i <- which(out$cell_id == 999L)
  expect_identical(out$n[i], 1L)
  expect_identical(out$resp_mean_z[i], L$z[7L])
  # The plain columns of a one-row cell are as they were: no spread, no error.
  expect_true(all(is.na(c(out$..sd_resp_z[i], out$..se_resp_z[i], out$..ci_lo_resp_z[i],
                          out$..ci_hi_resp_z[i]))))
  # Its pooled value: the weight of one row on its value, the rest on its target.
  w1 <- rec$tau2 / (rec$tau2 + rec$sigma2)
  expect_equal(out$..pooled_w_resp_z[i], w1, tolerance = 1e-12)
  expect_identical(out$..pooled_w_resp_z[i], min(out$..pooled_w_resp_z))
  expect_equal(out$..pooled_resp_z[i], w1 * L$z[7L] + (1 - w1) * out$..pooled_target_resp_z[i],
               tolerance = 1e-12)
  expect_false(is.na(out$..pooled_target_resp_z[i]))
  expect_identical(rec$share_alone, 1 / 720)
  # Left out, the row is predicted by the target of its cell.
  e <- .pe_residuals(L$z, L$cell_id, out, rec, "resp_z")
  expect_equal(e[7L], L$z[7L] - out$..pooled_target_resp_z[i], tolerance = 1e-12)
  expect_equal(rec$loo, mean(e^2), tolerance = 1e-10)
  # The interval has the offsets of every other cell.
  expect_identical(out$..obs_lo_resp_z[i], out$..pooled_resp_z[i] + rec$offsets[1L])
  expect_identical(out$..obs_hi_resp_z[i], out$..pooled_resp_z[i] + rec$offsets[2L])
  expect_equal(rec$offsets, stats::quantile(e, c(0.05, 0.95), names = FALSE, type = 7),
               tolerance = 1e-10)
  expect_equal(unique(round(out$..obs_hi_resp_z - out$..obs_lo_resp_z, 10)),
               round(diff(rec$offsets), 10))
})

test_that("a cell with no scored row of a column gets NA, and is no cell's neighbour", {
  L <- .pe_layer()
  L$p <- L$z + 1
  cell <- 23L
  L$p[L$cell_id == cell] <- NA                         # every row of the cell misses p
  r <- .pe_catch(.pe_pooled(L, "z", "p", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  out <- r$value
  i <- which(out$cell_id == cell)
  # No value is made up from the neighbours, and no interval.
  for (col in grep("pred_p$", .pe_new(out), value = TRUE)) {
    expect_identical(out[[col]][i], NA_real_, info = col)
    expect_false(anyNA(out[[col]][-i]), info = col)
  }
  expect_false(is.na(out$..pooled_resp_z[i]))
  expect_true(is.nan(out$pred_mean_p[i]) || is.na(out$pred_mean_p[i]))
  expect_identical(attr(out, "pooled")$columns$pred_p$cells, 59L)
  # The other cells' values are those of a layer without the cell.
  less <- .pe_pooled(L[L$cell_id != cell, ], predictor_vars = "p", obs_level = 0.9)
  for (col in grep("pred_p$", .pe_new(out), value = TRUE))
    expect_identical(out[[col]][-i], less[[col]], info = col)
  # A polygon no row fell in: NA in the joined frame.
  squares <- sf::st_sf(cell_id = c(cell, 5000L, 1L),
                       geometry = sf::st_sfc(lapply(1:3, function(k) sf::st_polygon(list(
                         cbind(c(0, 1, 1, 0, 0) + 2 * k, c(0, 0, 1, 1, 0))))), crs = 32617))
  joined <- suppressWarnings(.pe_pooled(L, "z", "p", obs_level = 0.9, cells_sf = squares))
  expect_identical(joined$cell_id, c(cell, 5000L, 1L))
  expect_true(all(is.na(sf::st_drop_geometry(joined)[2L, .pe_new(joined)])))
  expect_identical(joined$..pooled_resp_z[c(1L, 3L)], out$..pooled_resp_z[match(c(cell, 1L), out$cell_id)])
  expect_identical(joined$..pooled_pred_p[1L], NA_real_)
})


# ---- Rows that take no part -----------------------------------------------------

test_that("rows with a missing or infinite value are not rows of the column", {
  L <- .pe_layer()
  L$z[c(3L, 9L)] <- c(NA, NaN)
  L$z[11L] <- Inf
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_nonfinite")
  expect_identical(.pe_fields(r$warnings[[1L]])[c("n", "columns")], list(n = 1L, columns = "z"))
  rec <- attr(r$value, "pooled")$columns$resp_z
  expect_identical(rec$n, 717L)
  keep <- !(seq_len(nrow(L)) %in% c(3L, 9L, 11L))
  want <- .pe_pooled(L[keep, ], "z", obs_level = 0.9)
  for (col in .pe_new(want)) expect_identical(r$value[[col]], want[[col]], info = col)
  expect_identical(rec, attr(want, "pooled")$columns$resp_z)
})

test_that("rows without a cell take no part, and their summary row has NA in the new columns", {
  L <- .pe_layer()
  loose <- c(4L, 100L, 250L, 400L, 555L)
  K <- L
  K$cell_id[loose] <- NA
  r <- .pe_catch(.pe_pooled(K, "z", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  out <- r$value
  i <- which(is.na(out$cell_id))
  expect_length(i, 1L)
  expect_identical(out$n[i], 5L)
  expect_equal(out$resp_mean_z[i], mean(L$z[loose]))             # the plain columns are kept
  expect_false(is.na(out$..se_resp_z[i]))
  expect_true(all(is.na(sf::st_drop_geometry(out)[i, .pe_new(out)])))
  want <- .pe_pooled(L[-loose, ], "z", obs_level = 0.9)
  for (col in .pe_new(want)) expect_identical(out[[col]][-i], want[[col]], info = col)
  expect_identical(attr(out, "pooled")$columns$resp_z, attr(want, "pooled")$columns$resp_z)
  expect_identical(attr(out, "pooled")$columns$resp_z$n, 715L)
})

test_that("rows without a location count in the plain columns and not in the pooled ones", {
  L <- .pe_layer()
  L$p <- L$z * 2
  # Two rows of one cell and one of another lose their point.
  nowhere <- c(which(L$cell_id == 17L)[1:2], which(L$cell_id == 44L)[1L])
  L$p[nowhere[2L]] <- NA                               # this row has no value of p
  g <- sf::st_geometry(L)
  for (i in nowhere) g[[i]] <- sf::st_point()
  sf::st_geometry(L) <- g
  r <- .pe_catch(.pe_pooled(L, "z", "p", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_rows_skipped")
  # Rows are counted once, whichever columns they have a value of.
  # The columns are those with such a row, in the order of the frame.
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(n = 3L, cells = 2L, columns = c("z", "p")))
  only_z <- L
  only_z$p[nowhere] <- NA
  expect_identical(.pe_fields(.pe_catch(.pe_pooled(only_z, "z", "p"))$warnings[[1L]]),
                   list(n = 3L, cells = 2L, columns = "z"))
  rec <- attr(r$value, "pooled")$columns
  expect_identical(c(rec$resp_z$skipped, rec$pred_p$skipped), c(3L, 2L))
  expect_identical(c(rec$resp_z$n, rec$pred_p$n), c(717L, 717L))
  # n and the plain mean still hold the rows.
  plain <- summarize_by_cell(L, "z", "p", id_col = "cell_id", estimate = "mean")
  expect_identical(r$value$n, plain$n)
  expect_identical(sum(r$value$n), 720L)
  expect_identical(r$value$resp_mean_z, plain$resp_mean_z)
  # The pooled columns are those of the located rows.
  want <- .pe_pooled(L[-nowhere, ], "z", "p", obs_level = 0.9)
  for (col in .pe_new(want)) expect_identical(r$value[[col]], want[[col]], info = col)
  # A row without a location and without a cell is no row of any cell: nothing is said.
  M <- .pe_layer()
  g <- sf::st_geometry(M); g[[12L]] <- sf::st_point(); sf::st_geometry(M) <- g
  M$cell_id[12L] <- NA
  expect_length(.pe_catch(.pe_pooled(M, "z"))$warnings, 0L)
  # Nor is a row without a location and without a value of any pooled column.
  M <- .pe_layer()
  g <- sf::st_geometry(M); g[[12L]] <- sf::st_point(); sf::st_geometry(M) <- g
  M$z[12L] <- NA
  expect_length(.pe_catch(.pe_pooled(M, "z"))$warnings, 0L)
})


# ---- Pooling unavailable ----------------------------------------------------------

test_that("a layer without geometry gets NA pooled columns and one warning", {
  L <- .pe_layer()
  L$p <- L$z + 2
  df <- sf::st_drop_geometry(L)
  r <- .pe_catch(.pe_pooled(df, "z", "p", obs_level = 0.9, conf_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_unavailable")
  w <- r$warnings[[1L]]
  expect_identical(.pe_fields(w), list(reason = "no coordinates", columns = c("z", "p")))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): estimate = \"pooled\" needs the observations' locations, and ",
           "`assigned_points_sf` carries no geometry: the pooled columns are NA. The plain ",
           "columns are unaffected."))
  expect_null(conditionCall(w))
  out <- r$value
  expect_length(.pe_new(out), 10L)
  for (col in .pe_new(out)) expect_identical(out[[col]], rep(NA_real_, nrow(out)), info = col)
  rec <- attr(out, "pooled")
  expect_true(is.na(rec$crs))
  for (nm in c("resp_z", "pred_p")) {
    expect_identical(rec$columns[[nm]]$status, "no coordinates", info = nm)
    expect_identical(rec$columns[[nm]]$n, 0L, info = nm)
    expect_identical(rec$columns[[nm]]$offsets, c(NA_real_, NA_real_), info = nm)
  }
  # It is a warning and not an error, and the plain columns are the ones of
  # the call of plain means.
  plain <- summarize_by_cell(df, "z", "p", id_col = "cell_id", conf_level = 0.9,
                             estimate = "mean")
  expect_identical(setdiff(names(out), .pe_new(out)), names(plain))
  for (col in names(plain)) expect_identical(out[[col]], plain[[col]], info = col)
  # The table of intervals is there all the same.
  expect_identical(nrow(attr(out, "intervals")), 2L)
})

test_that("a column with no row that has a value, a cell and a location is NA", {
  L <- .pe_layer()
  L$none <- NA_real_
  r <- .pe_catch(.pe_pooled(L, "none", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_unavailable")
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(reason = "no located value", columns = "none"))
  expect_identical(
    conditionMessage(r$warnings[[1L]]),
    paste0("summarize_by_cell(): the pooled columns of 'none' are NA: no row has a value of ",
           "it, a cell and a location."))
  out <- r$value
  for (col in .pe_new(out)) expect_identical(out[[col]], rep(NA_real_, 60L), info = col)
  rec <- attr(out, "pooled")$columns$resp_none
  expect_identical(c(rec$n, rec$cells, rec$skipped), c(0L, 0L, 0L))
  expect_identical(rec$status, "no located value")
  expect_identical(c(rec$sigma2, rec$tau2, rec$loo), rep(NA_real_, 3L))
  # A layer of no rows: a frame of no rows with the columns, and the same status.
  r <- .pe_catch(.pe_pooled(L[0L, ], "z", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_unavailable")
  expect_identical(r$warnings[[1L]]$reason, "no located value")
  expect_identical(nrow(r$value), 0L)
  expect_length(.pe_new(r$value), 5L)
})

test_that("with no cell of two rows the pooled columns are NA", {
  # Twenty locations with two rows each, every row a cell of its own.
  L <- .pe_layer()[c(1:20, 1:20), ]
  L$z[21:40] <- L$z[21:40] + 0.25
  L$cell_id <- 101:140
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  # One warning: a column without an interval is not warned about its
  # repeated locations, nor about the few residuals it would rest on.
  expect_identical(.pe_classes(r), "spatialkit_pool_unavailable")
  expect_identical(attr(r$value, "pooled")$columns$resp_z$repeated, 40L)
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(reason = "no cell holds two values", columns = "z", n = c(z = 40L),
                        cells = c(z = 40L)))
  out <- r$value
  for (col in .pe_new(out)) expect_identical(out[[col]], rep(NA_real_, 40L), info = col)
  rec <- attr(out, "pooled")$columns$resp_z
  expect_identical(rec$status, "no cell holds two values")
  expect_identical(c(rec$n, rec$cells, rec$df), c(40L, 40L, 0L))
  expect_identical(rec$share_alone, 1)
  expect_identical(rec$sigma2, NA_real_)
  expect_identical(out$resp_mean_z, L$z)               # the plain mean of each row's cell
  # Five such cells are few cells, and nothing is said of their targets:
  # there are none.
  r <- .pe_catch(.pe_pooled(L[1:5, ], "z"))
  expect_identical(.pe_classes(r), "spatialkit_pool_unavailable")
  expect_identical(r$warnings[[1L]]$cells, c(z = 5L))
  # One cell of two rows is enough for a fit, and then the weights rest on one
  # value beyond one per cell.
  L$cell_id[2L] <- 101L
  r <- .pe_catch(.pe_pooled(L, "z"))
  expect_true("spatialkit_pool_few_df" %in% .pe_classes(r))
  expect_false("spatialkit_pool_unavailable" %in% .pe_classes(r))
  expect_identical(.pe_of(r, "spatialkit_pool_few_df")[[1L]]$df, c(z = 1L))
  expect_false(anyNA(r$value$..pooled_resp_z))
})

test_that("the limit on an exact search is 2e9 pairs, and a status in the summary", {
  expect_identical(.EXACT_PAIRS_MAX, 2e9)
  expect_identical(.exact_pairs_max(), 2e9)
  # The neighbours of K cells are K^2 pairs: about 44,700 cells at the limit.
  expect_identical(floor(sqrt(.EXACT_PAIRS_MAX)), 44721)
  L <- .pe_layer()
  .pe_without_fnn()
  withr::local_options(spatialkit.exact_pairs_max = 3599)      # 60 cells are 3,600 pairs
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_unavailable")
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(reason = "work limit", columns = "z", work = 3600, limit = 3599,
                        cells = c(z = 60L)))
  expect_match(conditionMessage(r$warnings[[1L]]), "Install FNN: the result is the same.",
               fixed = TRUE)
  expect_identical(attr(r$value, "pooled")$columns$resp_z$status, "work limit")
  expect_true(all(is.na(r$value$..obs_lo_resp_z)))
  withr::local_options(spatialkit.exact_pairs_max = 3600)
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  expect_identical(attr(r$value, "pooled")$columns$resp_z$status, "ok")
})


# ---- One cell ------------------------------------------------------------------------

test_that("one cell: the value is the plain mean, and the interval is read from the plain residuals", {
  L <- .pe_layer()[1:150, ]
  L$cell_id <- 5L
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  w <- r$warnings[[1L]]
  expect_identical(.pe_fields(w), list(reason = "one cell", columns = "z", n = c(z = 150L)))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): the values of 'z' lie in one cell, so there is nothing to ",
           "pool toward: the pooled value is the plain mean and its weight is 1."))
  out <- r$value
  expect_identical(nrow(out), 1L)
  expect_equal(out$..pooled_resp_z, mean(L$z), tolerance = 1e-12)
  expect_equal(out$..pooled_resp_z, out$resp_mean_z, tolerance = 1e-12)
  expect_identical(out$..pooled_w_resp_z, 1)
  expect_identical(out$..pooled_target_resp_z, NA_real_)
  rec <- attr(out, "pooled")$columns$resp_z
  expect_identical(rec$status, "one cell")
  expect_identical(rec$tau2, NA_real_)
  expect_identical(c(rec$cells, rec$neighbours, rec$df), c(1L, 0L, 149L))
  expect_equal(rec$sigma2, stats::var(L$z), tolerance = 1e-12)
  expect_identical(rec$structured, NA)
  # The residuals are the plain ones, (y - m) n / (n - 1).
  e <- (L$z - mean(L$z)) * 150 / 149
  expect_equal(rec$loo, mean(e^2), tolerance = 1e-12)
  expect_equal(rec$loo_one_cell, rec$loo, tolerance = 1e-12)
  expect_equal(rec$offsets, stats::quantile(e, c(0.05, 0.95), names = FALSE, type = 7),
               tolerance = 1e-10)
  expect_identical(c(out$..obs_lo_resp_z, out$..obs_hi_resp_z), out$..pooled_resp_z + rec$offsets)
  # Twenty rows in one cell: nothing is said of the values beyond one per
  # cell, since no weight rests on them.
  r <- .pe_catch(.pe_pooled(L[1:20, ], "z"))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  expect_identical(attr(r$value, "pooled")$columns$resp_z$df, 19L)
})

test_that("one row in one cell: a value, a weight of 1, and no spread, error or interval", {
  L <- .pe_layer()[42L, ]
  r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(reason = "one cell", columns = "z", n = c(z = 1L)))
  out <- r$value
  expect_identical(out$..pooled_resp_z, L$z)
  expect_identical(out$..pooled_w_resp_z, 1)
  expect_identical(c(out$..pooled_target_resp_z, out$..obs_lo_resp_z, out$..obs_hi_resp_z),
                   rep(NA_real_, 3L))
  rec <- attr(out, "pooled")$columns$resp_z
  expect_identical(rec$status, "one cell")
  expect_identical(c(rec$sigma2, rec$tau2, rec$loo, rec$loo_se, rec$loo_one_cell),
                   rep(NA_real_, 5L))
  expect_identical(rec$offsets, c(NA_real_, NA_real_))
  expect_identical(c(rec$n, rec$cells, rec$df), c(1L, 1L, 0L))
  expect_identical(rec$share_alone, 1)
})


# ---- Few cells, and few values beyond one per cell -------------------------------------

test_that("with 2 to 9 cells the target is the mean of all the other cells, with a warning", {
  base <- .pe_layer()
  for (K in c(2L, 3L, 9L)) {
    L <- base
    tag <- sprintf("%d cells", K)
    # K bands of equal width from west to east, so their counts differ.
    L$cell_id <- as.integer(cut(sf::st_coordinates(L)[, 1L], K, labels = FALSE))
    r <- .pe_catch(.pe_pooled(L, "z", obs_level = 0.9))
    expect_identical(.pe_classes(r), "spatialkit_pool_few_cells", info = tag)
    out <- r$value
    expect_gt(length(unique(out$n)), 1L)
    expect_gt(diff(range(out$..pooled_w_resp_z)), 0)
    w <- r$warnings[[1L]]
    expect_identical(.pe_fields(w),
                     list(columns = "z", cells = c(z = K), own_weight = c(z = min(out$..pooled_w_resp_z))),
                     info = tag)
    expect_identical(
      conditionMessage(w),
      sprintf(paste0(
        "summarize_by_cell(): 'z' has values in fewer than 10 cells (%d), so each cell's target ",
        "is the mean of all the other cells, not of cells near it. A cell with a small own ",
        "weight then gets a pooled value on the far side of the overall mean from its plain ",
        "mean, and two cells can change order. The `..pooled_w_*` columns hold the weights; ",
        "estimate = \"mean\" gives the plain means."), K))
    # Each target is the mean of every row outside the cell.
    S <- as.numeric(tapply(L$z, L$cell_id, sum))
    expect_equal(out$..pooled_target_resp_z, (sum(L$z) - S) / (720 - out$n), tolerance = 1e-12,
                 info = tag)
    rec <- attr(out, "pooled")$columns$resp_z
    expect_identical(c(rec$cells, rec$neighbours), c(K, K - 1L), info = tag)
    expect_identical(rec$status, "ok", info = tag)
    # The interval is computed.
    expect_false(anyNA(c(out$..obs_lo_resp_z, out$..obs_hi_resp_z)), info = tag)
    expect_equal(rec$offsets,
                 stats::quantile(.pe_residuals(L$z, L$cell_id, out, rec, "resp_z"),
                                 c(0.05, 0.95), names = FALSE, type = 7), tolerance = 1e-10)
  }
  # Ten cells are not few.
  L <- base
  L$cell_id <- as.integer(cut(sf::st_coordinates(L)[, 1L], 10L, labels = FALSE))
  expect_false("spatialkit_pool_few_cells" %in% .pe_classes(.pe_catch(.pe_pooled(L, "z"))))
})

test_that("two cells with no difference between them: the pooled values are the plain means swapped", {
  # 60 rows in each of two cells, one distribution for both.
  d <- local({
    .pe_seed(3011)
    data.frame(x = c(stats::runif(60, 0, 100), stats::runif(60, 200, 300)),
               y = stats::runif(120, 0, 100), v = stats::rnorm(120))
  })
  L <- .pe_points(d$x, d$y, data.frame(v = d$v, cell_id = rep(c("east", "west"), each = 60L)))
  r <- .pe_catch(.pe_pooled(L, "v", obs_level = 0.9))
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_v
  expect_identical(rec$tau2, 0)
  expect_identical(rec$status, "no signal between cells")
  # Three events for one column, in the order of the design's table.
  expect_identical(.pe_classes(r), c("spatialkit_pool_degenerate", "spatialkit_pool_few_cells",
                                     "spatialkit_no_structure"))
  expect_identical(.pe_of(r, "spatialkit_pool_degenerate")[[1L]]$reason, "no signal between cells")
  expect_identical(.pe_of(r, "spatialkit_pool_few_cells")[[1L]]$own_weight, c(v = 0))
  # The two plain means differ, and each cell's pooled value is the other's mean.
  expect_gt(abs(diff(out$resp_mean_v)), 0.01)
  expect_equal(out$..pooled_resp_v, rev(out$resp_mean_v), tolerance = 1e-12)
  expect_identical(out$..pooled_w_resp_v, c(0, 0))
  expect_identical(out$..pooled_resp_v, out$..pooled_target_resp_v)
  # So each lies on the far side of the overall mean from its own plain mean.
  overall <- mean(L$v)
  expect_true(all(sign(out$..pooled_resp_v - overall) == -sign(out$resp_mean_v - overall)))
  # The interval is computed: the residuals are real.
  expect_false(anyNA(c(out$..obs_lo_resp_v, out$..obs_hi_resp_v)))
  expect_gt(diff(rec$offsets), 1)
})

test_that("few values beyond one per cell: the values are returned with spatialkit_pool_few_df", {
  # 120 cells of which 117 hold one row and 3 hold two.
  d <- local({
    .pe_seed(3002)
    x <- stats::runif(123, 0, 1e5); y <- stats::runif(123, 0, 1e5)
    data.frame(x = x, y = y, v = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(123, 0, 0.3))
  })
  L <- .pe_points(d$x, d$y, data.frame(v = d$v, cell_id = c(1:120, 1:3)))
  r <- .pe_catch(.pe_pooled(L, "v"))
  # Everything the layer raises, in the order of the design's table.
  expect_identical(.pe_classes(r), c("spatialkit_pool_degenerate", "spatialkit_pool_few_df",
                                     "spatialkit_no_structure"))
  w <- .pe_of(r, "spatialkit_pool_few_df")
  expect_length(w, 1L)
  expect_identical(.pe_fields(w[[1L]]),
                   list(columns = "v", needed = 30L, df = c(v = 3L), n = c(v = 123L),
                        cells = c(v = 120L)))
  expect_identical(
    conditionMessage(w[[1L]]),
    paste0("summarize_by_cell(): the spread of 'v' inside a cell is estimated from fewer than ",
           "30 values beyond one per cell (3), and every own weight rests on it: the pooled ",
           "values and their weights are unstable."))
  rec <- attr(r$value, "pooled")$columns$resp_v
  expect_identical(rec$df, 3L)
  expect_equal(rec$share_alone, 117 / 123)
  # sigma2 is the spread of the three pairs, and the values are as computed.
  pairs <- L$v[c(1:3, 121:123)]
  expect_equal(rec$sigma2, sum((pairs[1:3] - pairs[4:6])^2 / 2) / 3, tolerance = 1e-12)
  expect_false(anyNA(r$value$..pooled_resp_v))
  expect_true(all(r$value$..pooled_w_resp_v >= 0 & r$value$..pooled_w_resp_v <= 1))
  # The field is not called structureless.  Its cells predict a left-out row
  # several times better than one mean; what fails is the F test, which has
  # three values beyond one per cell behind it, and the warning says that.
  expect_false(rec$structured)
  expect_gt(rec$loo_one_cell, 3 * rec$loo)
  expect_gte(rec$structure_p, 0.05)
  ns <- .pe_of(r, "spatialkit_no_structure")[[1L]]
  expect_identical(ns$beats, c(v = TRUE))
  expect_identical(
    conditionMessage(ns),
    sprintf(paste0(
      "summarize_by_cell(): the pooled values of 'v' predict a left-out observation better ",
      "than one mean for all rows, but the F test does not tell the cell means from one mean ",
      "(leave-one-out error %s against %s; F test of the cells against one mean, p = %s; ",
      "median own weight %.2f), so no spatial structure these cells can use is shown. The ",
      "values are returned."),
      format(signif(rec$loo, 4)), format(signif(rec$loo_one_cell, 4)),
      format(signif(rec$structure_p, 2)), rec$own_weight))
  expect_false(grepl("no better than one mean", conditionMessage(ns), fixed = TRUE))
  expect_false(grepl("looks smooth", conditionMessage(ns), fixed = TRUE))
  # Beside a column of white noise, for which one mean does as well as the
  # cells: one warning, a sentence for each of the two, and the numbers of
  # both behind the sentences, in the order of the frame.
  L$wn <- local({ .pe_seed(3021); stats::rnorm(123) })
  r2  <- .pe_catch(.pe_pooled(L, "v", "wn"))
  ns2 <- .pe_of(r2, "spatialkit_no_structure")
  expect_length(ns2, 1L)
  expect_identical(ns2[[1L]]$columns, c("v", "wn"))
  expect_identical(ns2[[1L]]$beats, c(v = TRUE, wn = FALSE))
  rw <- attr(r2$value, "pooled")$columns$pred_wn
  expect_lt(rw$loo_one_cell, rw$loo)
  piece <- function(r1) sprintf(
    paste0("'%s': leave-one-out error %s against %s; F test of the cells against one mean, ",
           "p = %s; median own weight %.2f."),
    r1$column, format(signif(r1$loo, 4)), format(signif(r1$loo_one_cell, 4)),
    format(signif(r1$structure_p, 2)), r1$own_weight)
  expect_identical(
    conditionMessage(ns2[[1L]]),
    paste(paste0(
      "summarize_by_cell(): the pooled values of 'wn' predict a left-out observation no better ",
      "than one mean for all rows, so it shows no spatial structure these cells can use. The ",
      "values are returned. A map of them looks smooth because each is mostly its neighbours' ",
      "mean, not because the response varies smoothly. The pooled values of 'v' predict a ",
      "left-out observation better than one mean for all rows, but the F test does not tell ",
      "the cell means from one mean, so no spatial structure these cells can use is shown for ",
      "it either."), piece(rec), piece(rw)))
  # Thirty such values are enough; twenty-nine are not.
  d <- local({
    .pe_seed(3004)
    x <- stats::runif(90, 0, 1e5); y <- stats::runif(90, 0, 1e5)
    data.frame(x = x, y = y, v = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(90, 0, 0.3))
  })
  M <- .pe_points(d$x, d$y, data.frame(v = d$v, cell_id = c(1:60, 1:30)))
  r30 <- .pe_catch(.pe_pooled(M, "v"))
  expect_identical(attr(r30$value, "pooled")$columns$resp_v$df, 30L)
  expect_false("spatialkit_pool_few_df" %in% .pe_classes(r30))
  r29 <- .pe_catch(.pe_pooled(M[-90L, ], "v"))
  expect_identical(.pe_of(r29, "spatialkit_pool_few_df")[[1L]]$df, c(v = 29L))
})


# ---- No spread, no signal, no structure --------------------------------------------------

test_that("a constant column: weights of 1, a zero-width interval, no leave-one-out error", {
  L <- .pe_layer()
  L$k <- 0.1
  r <- .pe_catch(.pe_pooled(L, "k", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(reason = "no spread inside any cell", columns = "k", constant = c(k = TRUE)))
  out <- r$value
  expect_identical(out$..pooled_w_resp_k, rep(1, 60L))
  expect_equal(out$..pooled_resp_k, rep(0.1, 60L), tolerance = 1e-14)
  # Zero width: both ends are the value itself.
  expect_identical(out$..obs_lo_resp_k, out$..pooled_resp_k)
  expect_identical(out$..obs_hi_resp_k, out$..pooled_resp_k)
  rec <- attr(out, "pooled")$columns$resp_k
  expect_identical(rec$status, "no spread inside any cell")
  expect_true(rec$constant)
  expect_identical(c(rec$sigma2, rec$tau2), c(0, 0))
  expect_identical(c(rec$loo, rec$loo_se, rec$loo_one_cell), rep(NA_real_, 3L))
  expect_identical(rec$offsets, c(0, 0))
  expect_identical(rec$structured, NA)
  expect_identical(rec$own_weight, 1)
})

test_that("a column constant inside every cell and not overall: weights of 1, residuals of 0", {
  L <- .pe_layer()
  L$step <- 0.1 * L$cell_id                             # one value a cell, 60 values in all
  r <- .pe_catch(.pe_pooled(L, "step", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(reason = "no spread inside any cell", columns = "step",
                        constant = c(step = FALSE)))
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_step
  expect_identical(rec$status, "no spread inside any cell")
  expect_false(rec$constant)
  expect_identical(rec$sigma2, 0)
  expect_gt(rec$tau2, 0)
  expect_identical(out$..pooled_w_resp_step, rep(1, 60L))
  expect_equal(out$..pooled_resp_step, 0.1 * out$cell_id, tolerance = 1e-12)
  # Every row shares its cell, so every residual is 0: an error of 0 and an
  # interval of no width.
  expect_identical(c(rec$loo, rec$loo_se), c(0, 0))
  expect_identical(rec$offsets, c(0, 0))
  expect_identical(out$..obs_lo_resp_step, out$..pooled_resp_step)
  expect_gt(rec$loo_one_cell, 0)
  # With a fifth of the rows alone in their cell, those rows are predicted by
  # their targets and the ends of the interval are theirs.
  alone <- seq(1L, nrow(L), by = 5L)
  L$cell_id[alone] <- 1000L + seq_along(alone)
  L$step <- 0.1 * (L$cell_id %% 97L)
  r <- .pe_catch(.pe_pooled(L, "step", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_step
  expect_identical(rec$status, "no spread inside any cell")
  expect_equal(rec$share_alone, 0.2)
  e <- L$step - out$..pooled_target_resp_step[match(L$cell_id, out$cell_id)]
  e[!(seq_len(nrow(L)) %in% alone)] <- 0
  expect_equal(rec$loo, mean(e^2), tolerance = 1e-12)
  expect_gt(rec$loo, 0)
  expect_equal(rec$offsets, stats::quantile(e, c(0.05, 0.95), names = FALSE, type = 7),
               tolerance = 1e-10)
  expect_lt(rec$offsets[1L], 0); expect_gt(rec$offsets[2L], 0)
  # With 3 percent of the rows alone, fewer than 5 percent of the residuals
  # lie on either side of 0, and both ends are 0.
  L <- .pe_layer()
  alone <- seq(1L, nrow(L), by = 34L)
  L$cell_id[alone] <- 1000L + seq_along(alone)
  L$step <- 0.1 * (L$cell_id %% 97L)
  rec <- attr(suppressWarnings(.pe_pooled(L, "step", obs_level = 0.9)), "pooled")$columns$resp_step
  expect_lt(rec$share_alone, 0.05)
  expect_gt(rec$loo, 0)
  expect_identical(rec$offsets, c(0, 0))
})

test_that("a column without spread is warned of that, and of nothing its weights rest on", {
  # 30 rows in 5 cells: few cells, 25 values beyond one per cell, and fewer
  # residuals than a 90 percent interval wants.  `v` varies inside the cells,
  # `k` is constant, `step` is constant inside every cell and not overall.
  d <- local({
    .pe_seed(3012)
    x <- stats::runif(30, 0, 1e5); y <- stats::runif(30, 0, 1e5)
    data.frame(x = x, y = y, noise = stats::rnorm(30, 0, 0.2))
  })
  L <- .pe_points(d$x, d$y, data.frame(cell_id = rep(1:5, 6L)))
  L$v    <- 2 * L$cell_id + d$noise
  L$k    <- 0.1
  L$step <- 0.1 * L$cell_id
  # The column with spread raises the three, so the layer is one on which
  # they are due.
  r <- .pe_catch(.pe_pooled(L, "v", obs_level = 0.9))
  expect_identical(.pe_classes(r), c("spatialkit_pool_few_cells", "spatialkit_pool_few_df",
                                     "spatialkit_few_residuals"))
  expect_identical(attr(r$value, "pooled")$columns$resp_v$status, "ok")
  # The count the interval wants is the first whole number at which five
  # residuals are expected beyond each end: 34 at a level of 0.7, where
  # 5 / 0.15 is 33.3.
  fr <- .pe_of(.pe_catch(.pe_pooled(L, "v", obs_level = 0.7)), "spatialkit_few_residuals")[[1L]]
  expect_identical(fr$needed, 34L)
  expect_match(conditionMessage(fr),
               "the 70% interval of 'v' rests on fewer than 34 leave-one-out residuals (30)",
               fixed = TRUE)
  # A constant column: its weights are set to 1, its spread to 0 and its
  # residuals to 0.  One warning says so, and no other.
  r <- .pe_catch(.pe_pooled(L, "k", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_pool_degenerate")
  rec <- attr(r$value, "pooled")$columns$resp_k
  expect_identical(c(rec$cells, rec$df, rec$n), c(5L, 25L, 30L))
  expect_identical(rec$offsets, c(0, 0))
  expect_identical(r$value$..pooled_w_resp_k, rep(1, 5L))
  # Constant inside every cell: nothing of the cells or of the values beyond
  # one per cell.  Its interval is read from residuals, so the count of
  # those is still said.
  r <- .pe_catch(.pe_pooled(L, "step", obs_level = 0.9))
  expect_identical(.pe_classes(r), c("spatialkit_pool_degenerate", "spatialkit_few_residuals"))
  expect_identical(.pe_classes(.pe_catch(.pe_pooled(L, "step"))), "spatialkit_pool_degenerate")
  # The three columns in one call: each warning names the columns it is true of.
  r <- .pe_catch(.pe_pooled(L, "v", c("k", "step"), obs_level = 0.9))
  expect_identical(.pe_classes(r), c("spatialkit_pool_degenerate", "spatialkit_pool_few_cells",
                                     "spatialkit_pool_few_df", "spatialkit_few_residuals"))
  expect_identical(lapply(r$warnings, function(w) w$columns),
                   list(c("k", "step"), "v", "v", c("v", "step")))
  # Rows that share a location: said of the column whose interval is read
  # from residuals, not of the constant one.
  twice <- L[c(seq_len(30L), 1:5), ]
  expect_identical(attr(suppressWarnings(.pe_pooled(twice, "k")), "pooled")$columns$resp_k$repeated,
                   10L)
  expect_identical(.pe_classes(.pe_catch(.pe_pooled(twice, "k", obs_level = 0.9))),
                   "spatialkit_pool_degenerate")
  expect_true("spatialkit_repeated_locations" %in%
                .pe_classes(.pe_catch(.pe_pooled(twice, "v", obs_level = 0.9))))
})

test_that("no signal between cells: every weight is 0 and each value is its target", {
  L <- .pe_layer()
  L$w <- local({ .pe_seed(3101); stats::rnorm(nrow(L)) })
  r <- .pe_catch(.pe_pooled(L, "w", obs_level = 0.9))
  expect_identical(.pe_classes(r), c("spatialkit_pool_degenerate", "spatialkit_no_structure"))
  rec <- attr(r$value, "pooled")$columns$resp_w
  expect_identical(rec$status, "no signal between cells")
  expect_identical(rec$tau2, 0)
  expect_gt(rec$sigma2, 0)
  w <- r$warnings[[1L]]
  expect_identical(.pe_fields(w), list(reason = "no signal between cells", columns = "w",
                                       sigma2 = c(w = rec$sigma2)))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): every own weight of 'w' is 0: the cell means differ from ",
           "their neighbours' means by no more than the spread inside the cells explains. Each ",
           "pooled value is then the mean of the cell's nearest cells and holds none of the ",
           "cell's own values."))
  out <- r$value
  expect_identical(out$..pooled_w_resp_w, rep(0, 60L))
  expect_identical(out$..pooled_resp_w, out$..pooled_target_resp_w)
  expect_false(isTRUE(all.equal(out$..pooled_resp_w, out$resp_mean_w)))
  # The interval is computed: the residuals are those of the targets.
  e <- L$w - out$..pooled_target_resp_w[match(L$cell_id, out$cell_id)]
  expect_equal(rec$offsets, stats::quantile(e, c(0.05, 0.95), names = FALSE, type = 7),
               tolerance = 1e-10)
  expect_equal(rec$loo, mean(e^2), tolerance = 1e-12)
  expect_false(anyNA(out$..obs_lo_resp_w))
})

test_that("a response without structure: the values and the interval are returned, with a warning", {
  L <- .pe_layer()
  L$w <- local({ .pe_seed(3102); stats::rnorm(nrow(L)) })
  r <- .pe_catch(.pe_pooled(L, "w", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_no_structure")
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_w
  expect_identical(rec$status, "ok")
  expect_false(rec$structured)
  expect_gt(rec$tau2, 0)
  # One mean for all rows does as well as the cells, to the margin of the rule.
  expect_lte(rec$loo_one_cell, rec$loo * 1.02)
  expect_false(anyNA(out$..pooled_resp_w))
  # The interval is then about as wide as the response's own 5th to 95th percentiles.
  own <- diff(stats::quantile(L$w, c(0.05, 0.95), names = FALSE))
  expect_gt(diff(rec$offsets) / own, 0.9)
  expect_lt(diff(rec$offsets) / own, 1.15)
  # Each value is mostly its neighbours' mean.
  expect_lt(rec$own_weight, 0.4)
  # The field of the same layer is structured, and says nothing.
  expect_true(attr(.pe_pooled(L, "z"), "pooled")$columns$resp_z$structured)
})


# ---- Locations that repeat or coincide ------------------------------------------------------

test_that("rows at one location each count, and are left out one at a time", {
  # 240 locations with three rows each.
  d <- local({
    .pe_seed(3005)
    x <- rep(stats::runif(240, 0, 1e5), 3L); y <- rep(stats::runif(240, 0, 1e5), 3L)
    data.frame(x = x, y = y, v = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(720, 0, 0.5))
  })
  L <- .pe_points(d$x, d$y, data.frame(v = d$v))
  L$cell_id <- as.integer(6 * floor(d$x / 1e4) + floor(d$y / (1e5 / 6)) + 1)
  r <- .pe_catch(.pe_pooled(L, "v", obs_level = 0.9))
  expect_identical(.pe_classes(r), "spatialkit_repeated_locations")
  expect_identical(.pe_fields(r$warnings[[1L]]),
                   list(columns = "v", n = c(v = 720L), total = c(v = 720L)))
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_v
  expect_identical(rec$repeated, 720L)
  expect_identical(sum(out$n), 720L)                    # each row counts in its cell
  expect_true(all(out$n %% 3L == 0L))
  # One row is left out and the others at its location stay.
  e <- .pe_residuals(L$v, L$cell_id, out, rec, "resp_v")
  expect_equal(rec$loo, mean(e^2), tolerance = 1e-10)
  expect_equal(rec$offsets, stats::quantile(e, c(0.05, 0.95), names = FALSE, type = 7),
               tolerance = 1e-10)
  # Without the level the count is in the record and nothing is said.
  quiet <- .pe_catch(.pe_pooled(L, "v"))
  expect_length(quiet$warnings, 0L)
  expect_identical(attr(quiet$value, "pooled")$columns$resp_v$repeated, 720L)
  # A layer whose locations are all different counts none.
  expect_identical(attr(.pe_pooled(.pe_layer(), "z"), "pooled")$columns$resp_z$repeated, 0L)
})

test_that("rows at one location in different cells each count in their own cell", {
  # 360 locations, two rows at each, the two in cells of two overlapping
  # sets: 20 cells from west to east and 12 from south to north.
  d <- local({
    .pe_seed(3006)
    x <- stats::runif(360, 0, 1e5); y <- stats::runif(360, 0, 1e5)
    data.frame(x = rep(x, 2L), y = rep(y, 2L),
               v = sin(rep(x, 2L) / 2e4) + cos(rep(y, 2L) / 1.5e4) + stats::rnorm(720, 0, 0.4))
  })
  cell <- c(as.integer(floor(d$x[1:360] / 5e3)) + 1L,
            100L + as.integer(floor(d$y[361:720] / (1e5 / 12))) + 1L)
  L <- .pe_points(d$x, d$y, data.frame(v = d$v, cell_id = cell))
  r <- .pe_catch(.pe_pooled(L, "v"))
  expect_length(r$warnings, 0L)
  out <- r$value
  expect_identical(nrow(out), 32L)
  expect_identical(sum(out$n[out$cell_id < 100L]), 360L)
  expect_identical(sum(out$n[out$cell_id > 100L]), 360L)
  expect_identical(attr(out, "pooled")$columns$resp_v$repeated, 720L)
  # The values are the engine's on that grouping.
  xy <- unname(sf::st_coordinates(L)[, 1:2])
  u  <- sort(unique(cell))
  f  <- .pool_fit(d$v, match(cell, u), length(u), xy, want = "pooled")
  expect_identical(out$..pooled_resp_v, f$value)
  expect_identical(out$..pooled_target_resp_v, f$target)
  # Cells of the two sets are each other's neighbours where their rows'
  # centroids are near, whatever polygons they came from.
  cen <- rowsum(xy, match(cell, u)) / tabulate(match(cell, u))
  nb  <- .pe_neighbours(cen)
  expect_true(any(u[nb[1L, ]] > 100L) && any(u[nb[1L, ]] < 100L))
})

test_that("cells that share a centroid are not their own neighbours", {
  # A 5 x 4 grid of cells of four rows each on whole metres, and two more cells
  # whose rows have one centroid: (250, 250) for both.
  g <- expand.grid(cx = seq(100, 900, by = 200), cy = seq(100, 700, by = 200))
  x <- c(rep(g$cx, each = 4L) + c(-20, 20, -20, 20), 240, 260, 250, 250)
  y <- c(rep(g$cy, each = 4L) + c(-20, -20, 20, 20), 250, 250, 230, 270)
  cell <- c(rep(seq_len(nrow(g)), each = 4L), 21L, 21L, 22L, 22L)
  v <- local({ .pe_seed(3007); sin(x / 300) + cos(y / 250) + stats::rnorm(length(x), 0, 0.3) })
  L <- .pe_points(x, y, data.frame(v = v, cell_id = cell))
  out <- suppressWarnings(.pe_pooled(L, "v"))
  cen <- rowsum(cbind(x, y), cell) / tabulate(cell)
  expect_identical(unname(cen[21L, ]), unname(cen[22L, ]))
  nb <- .pe_neighbours(unname(cen))
  expect_identical(nb[21L, 1L], 22L)                    # each is the other's nearest
  expect_identical(nb[22L, 1L], 21L)
  S <- as.numeric(rowsum(v, cell)); n <- tabulate(cell)
  target <- rowSums(matrix(S[nb], 22L)) / rowSums(matrix(n[nb], 22L))
  expect_equal(out$..pooled_target_resp_v, target, tolerance = 1e-12)
  # A cell that was its own neighbour would have its own rows in its target.
  own <- (rowSums(matrix(S[nb[, -8L]], 22L)) + S) / (rowSums(matrix(n[nb[, -8L]], 22L)) + n)
  expect_gt(abs(out$..pooled_target_resp_v[21L] - own[21L]), 1e-6)
  # The exact search gives the kd-tree's answer.  (Without FNN the call above
  # ran on the exact search already.)
  skip_if_not_installed("FNN")
  exact <- local({ .pe_without_fnn(); suppressWarnings(.pe_pooled(L, "v")) })
  expect_identical(exact, out)
})

test_that("on gridded locations a tie between neighbours goes to the lower x, then the lower y", {
  # A 30 x 20 lattice at 100 m in cells of 2 x 2 locations: 150 cells whose
  # centroids are a lattice themselves, so candidate neighbours tie.
  g <- expand.grid(i = 0:29, j = 0:19)
  x <- 1000 + 100 * g$i; y <- 2000 + 100 * g$j
  cell <- as.integer(15 * (g$j %/% 2L) + g$i %/% 2L + 1L)
  v <- local({ .pe_seed(3008); sin(x / 900) + cos(y / 700) + stats::rnorm(length(x), 0, 0.3) })
  L <- .pe_points(x, y, data.frame(v = v, cell_id = cell))
  r <- .pe_catch(.pe_pooled(L, "v", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  out <- r$value
  cen <- unname(rowsum(cbind(x, y), cell) / tabulate(cell))
  nb  <- .pe_neighbours(cen)
  # The eighth and the ninth nearest tie for the two cells beside each corner
  # of the lattice: there the rule decides which cells are in the target.
  d2 <- as.matrix(stats::dist(cen))^2
  tied <- vapply(seq_len(150L), function(j) {
    s <- sort(d2[j, -j]); s[8L] == s[9L] }, NA)
  expect_identical(sum(tied), 8L)
  S <- as.numeric(rowsum(v, cell)); n <- tabulate(cell)
  expect_equal(out$..pooled_target_resp_v,
               rowSums(matrix(S[nb], 150L)) / rowSums(matrix(n[nb], 150L)), tolerance = 1e-12)
  # Ties settled the other way give other targets in the tied cells.
  other <- t(vapply(seq_len(150L), function(j) {
    o <- order(d2[j, ], -cen[, 1L], -cen[, 2L]); o[o != j][1:8] }, integer(8)))
  alt <- rowSums(matrix(S[other], 150L)) / rowSums(matrix(n[other], 150L))
  expect_identical(which(abs(alt - out$..pooled_target_resp_v) > 1e-9), which(tied))
  # The same values in any row order, and with either search.
  for (s in 1:3) {
    p <- local({ .pe_seed(3080 + s); sample.int(nrow(L)) })
    again <- .pe_pooled(L[p, ], "v", obs_level = 0.9)
    for (col in .pe_new(out)) expect_equal(again[[col]], out[[col]], tolerance = 1e-12,
                                           info = paste(col, s))
  }
  skip_if_not_installed("FNN")
  exact <- local({ .pe_without_fnn(); .pe_pooled(L, "v", obs_level = 0.9) })
  expect_identical(exact, out)
})


# ---- Columns ---------------------------------------------------------------------------------

test_that("a column named as response and as predictor is fitted once and returned under both", {
  L <- .pe_layer()
  real <- .pool_column
  force(real)
  fits <- 0L
  testthat::local_mocked_bindings(.pool_column = function(...) {
    fits <<- fits + 1L
    real(...)
  })
  out <- .pe_pooled(L, "z", "z", obs_level = 0.9)
  expect_identical(fits, 1L)
  for (stem in c("..pooled_", "..pooled_w_", "..pooled_target_", "..obs_lo_", "..obs_hi_"))
    expect_identical(out[[paste0(stem, "pred_z")]], out[[paste0(stem, "resp_z")]], info = stem)
  expect_identical(.pe_new(out),
                   c(paste0(c("..pooled_", "..pooled_w_", "..pooled_target_", "..obs_lo_",
                              "..obs_hi_"), "resp_z"),
                     paste0(c("..pooled_", "..pooled_w_", "..pooled_target_", "..obs_lo_",
                              "..obs_hi_"), "pred_z")))
  rec <- attr(out, "pooled")$columns
  expect_identical(names(rec), c("resp_z", "pred_z"))
  expect_identical(c(rec$resp_z$role, rec$pred_z$role), c("response", "predictor"))
  expect_identical(rec$resp_z[-2L], rec$pred_z[-2L])
  # A warning names the column once.
  L$w <- local({ .pe_seed(3102); stats::rnorm(nrow(L)) })
  r <- .pe_catch(.pe_pooled(L, "w", c("w", "z", "z")))
  expect_identical(.pe_classes(r), "spatialkit_no_structure")
  expect_identical(r$warnings[[1L]]$columns, "w")
  expect_identical(names(attr(r$value, "pooled")$columns), c("resp_w", "pred_w", "pred_z"))
})

test_that("a column that is not numeric gets no pooled column", {
  L <- .pe_layer()
  L$flag <- L$z > 0
  L$kind <- ifelse(L$z > 0, "high", "low")
  L$fac  <- factor(L$kind)
  # A logical response alone: today's warning, no new column and no record.
  r <- .pe_catch(.pe_pooled(L, "flag", obs_level = 0.9))
  expect_length(r$warnings, 1L)
  expect_identical(class(r$warnings[[1L]])[1L], "simpleWarning")
  expect_match(conditionMessage(r$warnings[[1L]]),
               "response column(s) 'flag' are not numeric and are not summarised", fixed = TRUE)
  expect_length(.pe_new(r$value), 0L)
  expect_null(attr(r$value, "pooled"))
  expect_null(attr(r$value, "intervals"))
  plain <- suppressWarnings(summarize_by_cell(L, "flag", id_col = "cell_id", estimate = "mean"))
  expect_identical(r$value, plain)
  # Beside numeric columns the others are skipped and the numeric ones pooled.
  r <- .pe_catch(.pe_pooled(L, "flag", c("kind", "z", "fac")))
  expect_length(r$warnings, 2L)
  expect_identical(.pe_new(r$value), paste0(c("..pooled_", "..pooled_w_", "..pooled_target_"),
                                            "pred_z"))
  expect_identical(names(attr(r$value, "pooled")$columns), "pred_z")
  # A response that is not a column, and no column at all.
  r <- .pe_catch(.pe_pooled(L, "absent"))
  expect_length(r$warnings, 1L)
  expect_null(attr(r$value, "pooled"))
  none <- .pe_pooled(L)
  expect_identical(none, summarize_by_cell(L, id_col = "cell_id", estimate = "mean"))
  expect_null(attr(none, "pooled"))
  # An integer column is numeric: it is pooled as its doubles.
  L$count <- as.integer(round(10 * L$z))
  a <- .pe_pooled(L, "count")
  L$count_dbl <- as.numeric(L$count)
  b <- .pe_pooled(L, "count_dbl")
  expect_identical(a$..pooled_resp_count, b$..pooled_resp_count_dbl)
  expect_type(a$..pooled_resp_count, "double")
})

test_that("both families of interval are returned under their own names", {
  L <- .pe_layer()
  both <- .pe_pooled(L, "z", obs_level = 0.9, conf_level = 0.8)
  obs  <- .pe_pooled(L, "z", obs_level = 0.9)
  conf <- .pe_pooled(L, "z", conf_level = 0.8)
  expect_true(all(c("..ci_lo_resp_z", "..ci_hi_resp_z", "..obs_lo_resp_z", "..obs_hi_resp_z") %in%
                    names(both)))
  # Neither changes the other.
  expect_identical(both$..obs_lo_resp_z, obs$..obs_lo_resp_z)
  expect_identical(both$..obs_hi_resp_z, obs$..obs_hi_resp_z)
  expect_identical(both$..ci_lo_resp_z, conf$..ci_lo_resp_z)
  expect_identical(both$..ci_hi_resp_z, conf$..ci_hi_resp_z)
  # They are different intervals: one is centred on the plain mean and narrows
  # with the count, the other has one width for every cell.
  expect_equal(both$..ci_lo_resp_z + both$..ci_hi_resp_z, 2 * both$resp_mean_z, tolerance = 1e-12)
  expect_length(unique(round(both$..obs_hi_resp_z - both$..obs_lo_resp_z, 10)), 1L)
  expect_gt(length(unique(round(both$..ci_hi_resp_z - both$..ci_lo_resp_z, 10))), 50L)
  expect_identical(attr(both, "intervals")$level, c(0.9, 0.8))
})


# ---- Coordinates -----------------------------------------------------------------------------

test_that("coordinates without a CRS are used as they are, or projected when they look like lon/lat", {
  L <- .pe_layer()
  base <- .pe_pooled(L, "z", obs_level = 0.9)
  # Planar numbers without a CRS: the same values, and the record says no CRS.
  bare <- sf::st_set_crs(L, NA)
  r <- .pe_catch(.pe_pooled(bare, "z", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  for (col in .pe_new(base)) expect_identical(r$value[[col]], base[[col]], info = col)
  expect_true(is.na(attr(r$value, "pooled")$crs))
  # Numbers that look like degrees, without a CRS: projected as the lon/lat
  # layer is, with the warning ensure_projected() raises for the assumption.
  ll <- sf::st_transform(L, 4326)
  want <- .pe_pooled(ll, "z", obs_level = 0.9)
  r <- .pe_catch(.pe_pooled(sf::st_set_crs(ll, NA), "z", obs_level = 0.9))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), "CRS is missing; assuming EPSG:4326", fixed = TRUE)
  for (col in .pe_new(want)) expect_identical(r$value[[col]], want[[col]], info = col)
  expect_false(isTRUE(sf::st_is_longlat(attr(r$value, "pooled")$crs)))
  # The projected values are the planar layer's to the accuracy of the projection.
  for (col in .pe_new(want)) expect_equal(want[[col]], base[[col]], tolerance = 1e-6, info = col)
})

test_that("features that are not points are read at one point each", {
  L <- .pe_layer()[1:300, ]
  squares <- sf::st_buffer(L, 300, endCapStyle = "SQUARE")
  expect_true(all(sf::st_geometry_type(squares) == "POLYGON"))
  r <- .pe_catch(.pe_pooled(squares, "z", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  expect_identical(r$value, .pe_pooled(coerce_to_points(squares, "auto"), "z", obs_level = 0.9))
  # A square around a point has the point as its centre, to rounding.
  pts <- .pe_pooled(L, "z", obs_level = 0.9)
  for (col in .pe_new(pts)) expect_equal(r$value[[col]], pts[[col]], tolerance = 1e-9, info = col)
  # Points stored as MULTIPOINT are the points themselves.
  multi <- sf::st_cast(L, "MULTIPOINT")
  expect_identical(.pe_pooled(multi, "z", obs_level = 0.9), pts)
  # The point of a polygon is the one coerce_to_points(, "auto") gives, a point
  # on its surface, and not its centroid: for an L the two are far apart.
  ells <- sf::st_sf(
    cell_id = 1:3,
    geometry = sf::st_sfc(lapply(c(0, 5000, 10000), function(dx) sf::st_polygon(list(cbind(
      dx + c(0, 3000, 3000, 200, 200, 0, 0), c(0, 0, 200, 200, 2000, 2000, 0))))), crs = 32617))
  at <- .pooled_xy(ells)
  expect_identical(at$xy, unname(sf::st_coordinates(coerce_to_points(ells, "auto"))[, 1:2]))
  cen <- unname(sf::st_coordinates(coerce_to_points(ells, "centroid"))[, 1:2])
  expect_gt(min(sqrt(rowSums((at$xy - cen)^2))), 100)
  expect_identical(at$crs, sf::st_crs(ells))
  expect_null(.pooled_xy(sf::st_drop_geometry(ells)))
  # A layer with a third coordinate: the first two are read.
  xyz <- sf::st_zm(L, drop = FALSE, what = "Z")
  expect_identical(ncol(sf::st_coordinates(xyz)), 3L)
  z3 <- .pe_pooled(xyz, "z", obs_level = 0.9)
  for (col in .pe_new(pts)) expect_identical(z3[[col]], pts[[col]], info = col)
})

test_that("points and other features in one layer are each read at their own point", {
  L <- .pe_layer()[1:300, ]
  ell <- function(p) sf::st_polygon(list(rbind(
    p, p + c(3000, 0), p + c(3000, 200), p + c(200, 200), p + c(200, 2000), p + c(0, 2000), p)))
  xy <- unname(sf::st_coordinates(L))
  ells <- sf::st_set_geometry(L, sf::st_sfc(lapply(seq_len(300L), function(i) ell(xy[i, ])),
                                            crs = sf::st_crs(L)))
  # Whether the layer needs reducing is asked of every feature, not of the
  # first or the last: points followed by polygons, and polygons by points.
  for (first in c("points", "polygons")) {
    mixed <- if (first == "points") rbind(L[1:150, ], ells[151:300, ])
             else rbind(ells[1:150, ], L[151:300, ])
    expect_setequal(as.character(unique(sf::st_geometry_type(mixed))), c("POINT", "POLYGON"))
    at <- .pooled_xy(mixed)
    expect_identical(dim(at$xy), c(300L, 2L), info = first)
    as_points <- coerce_to_points(mixed, "auto")
    expect_identical(at$xy, unname(sf::st_coordinates(as_points)[, 1:2]), info = first)
    # The points are kept as they are, and a polygon is read away from its corner.
    kept <- if (first == "points") 1:150 else 151:300
    expect_identical(at$xy[kept, ], xy[kept, ], info = first)
    expect_gt(min(sqrt(rowSums((at$xy[-kept, ] - xy[-kept, ])^2))), 50)
    r <- .pe_catch(.pe_pooled(mixed, "z", obs_level = 0.9))
    expect_length(r$warnings, 0L)
    expect_identical(r$value, .pe_pooled(as_points, "z", obs_level = 0.9), info = first)
  }
})

test_that("lon/lat features are reduced to points first and projected then", {
  # The order resolution_profile() takes.  Projected first, a polygon's point
  # on its surface is another point.
  L <- .pe_layer()[1:300, ]
  ell <- function(p) sf::st_polygon(list(rbind(
    p, p + c(4000, 0), p + c(4000, 500), p + c(500, 500), p + c(500, 4000), p + c(0, 4000), p)))
  xy <- unname(sf::st_coordinates(L))
  ells <- sf::st_transform(
    sf::st_set_geometry(L, sf::st_sfc(lapply(seq_len(300L), function(i) ell(xy[i, ])),
                                      crs = sf::st_crs(L))), 4326)
  expect_true(sf::st_is_longlat(ells))
  r <- .pe_catch(.pooled_xy(ells))
  # sf warns that a point on the surface of a lon/lat polygon may not be
  # correct; the point is projected before any distance is taken, and the
  # warning is not passed on.  The call of plain means raises none either.
  expect_length(r$warnings, 0L)
  as_points <- ensure_projected(suppressWarnings(coerce_to_points(ells, "auto")))
  expect_identical(r$value$xy, unname(sf::st_coordinates(as_points)[, 1:2]))
  expect_identical(r$value$crs, sf::st_crs(as_points))
  expect_false(isTRUE(sf::st_is_longlat(r$value$crs)))
  other <- unname(sf::st_coordinates(
    suppressWarnings(coerce_to_points(ensure_projected(ells), "auto")))[, 1:2])
  expect_gt(max(sqrt(rowSums((r$value$xy - other)^2))), 0.1)        # metres
  # Through the exported function: the values of those points, and no warning
  # that is sf's.
  got  <- .pe_catch(.pe_pooled(ells, "z", obs_level = 0.9))
  want <- .pe_catch(.pe_pooled(as_points, "z", obs_level = 0.9))
  expect_identical(.pe_classes(got), .pe_classes(want))
  expect_false(any(grepl("st_point_on_surface", vapply(got$warnings, conditionMessage, ""),
                         fixed = TRUE)))
  for (col in .pe_new(want$value))
    expect_identical(got$value[[col]], want$value[[col]], info = col)
  expect_identical(attr(got$value, "pooled")$crs, sf::st_crs(as_points))
  expect_length(.pe_catch(summarize_by_cell(ells, "z", id_col = "cell_id",
                                            estimate = "mean"))$warnings, 0L)
})

test_that("a coordinate too large to take a distance from is a row without a location", {
  expect_identical(.POOL_XY_MAX, 1e150)
  L  <- .pe_layer()
  at <- unname(sf::st_coordinates(L)[7L, ])
  put <- function(x, y) {
    far <- L
    g <- sf::st_geometry(far)
    g[[7L]] <- sf::st_point(c(x, y))
    sf::st_geometry(far) <- g
    far
  }
  want <- .pe_pooled(L[-7L, ], "z", obs_level = 0.9)
  # Both coordinates too large, the second alone, the first alone.
  for (p in list(c(1e200, 1e200), c(at[1L], .Machine$double.xmax), c(-1e151, at[2L]))) {
    tag <- paste(format(p), collapse = ", ")
    far <- put(p[1L], p[2L])
    r <- .pe_catch(.pe_pooled(far, "z", obs_level = 0.9))
    expect_identical(.pe_classes(r), "spatialkit_pool_rows_skipped", info = tag)
    expect_identical(.pe_fields(r$warnings[[1L]]), list(n = 1L, cells = 1L, columns = "z"),
                     info = tag)
    # The pooled columns are those of the layer without the row, the plain
    # ones still hold it, and the exact search gives the same frame.
    for (col in .pe_new(want)) expect_identical(r$value[[col]], want[[col]],
                                                info = paste(tag, col))
    expect_identical(sum(r$value$n), 720L)
    exact <- local({ .pe_without_fnn(); suppressWarnings(.pe_pooled(far, "z", obs_level = 0.9)) })
    expect_identical(exact, r$value, info = tag)
  }
  # A large coordinate a distance can be taken from is a location: the row
  # takes part, with either search.
  for (p in list(c(1e30, 1e30), c(at[1L], -1e150), c(1e150, at[2L]))) {
    tag <- paste(format(p), collapse = ", ")
    far <- put(p[1L], p[2L])
    r <- .pe_catch(.pe_pooled(far, "z", obs_level = 0.9))
    expect_false("spatialkit_pool_rows_skipped" %in% .pe_classes(r), info = tag)
    expect_identical(attr(r$value, "pooled")$columns$resp_z$n, 720L, info = tag)
    expect_false(anyNA(r$value$..pooled_resp_z), info = tag)
    exact <- local({ .pe_without_fnn(); suppressWarnings(.pe_pooled(far, "z", obs_level = 0.9)) })
    expect_identical(exact, r$value, info = tag)
  }
})
