# tests/testthat/test-core-pool-summary.R
# ---------------------------------------------------------------------------
# summarize_by_cell(estimate = "pooled", obs_level = ): the pooled columns,
# their record and their warnings, through the exported function.
#
# The pooled value, its own weight, its target, the interval for one
# observation and the leave-one-out error in the record are held to the
# reference functions of helper-pool-reference.R on cells built from real
# polygons.  The columns a call shares with estimate = "mean" are held
# identical() to that call's, so nothing the function returned before moves.
#
# Every call names `estimate`.  Unless a test says otherwise the cells are of
# the density kind from 3 k-means++ starts, built as a user builds them:
# voronoi_seeds_kmeans(), build_tessellation(), assign_features_to_polygons().
# A backend is switched off by replacing .has_backend(), as the code under
# test asks for FNN through that one function.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.ps_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.ps_store <- new.env(parent = emptyenv())
.ps_keep <- function(key, expr) {
  if (is.null(.ps_store[[key]])) .ps_store[[key]] <- expr
  .ps_store[[key]]
}

.ps_points <- function(x, y, data, crs = 32617)
  sf::st_as_sf(cbind(data.frame(x = x, y = y), data), coords = c("x", "y"), crs = crs)

# Cells of the density kind at `k` centres from 3 k-means++ starts, and the
# layer assigned to them.
.ps_assign <- function(pts, k) {
  seeds <- suppressWarnings(
    voronoi_seeds_kmeans(pts, k, set_seed = 123, nstart = 3, partition = "density",
                         init = "kmeans++"))
  tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  list(asg = assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id"),
       cells = tess$cells)
}

# A call with its warnings kept and muffled.
.ps_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.ps_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.ps_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)

# What a classed condition carries besides its text.
.ps_fields <- function(w) unclass(w)[setdiff(names(unclass(w)), c("message", "call"))]

.ps_without_fnn <- function(env = parent.frame())
  testthat::local_mocked_bindings(
    .has_backend = function(pkg)
      if (identical(pkg, "FNN")) FALSE else requireNamespace(pkg, quietly = TRUE),
    .env = env)

.ps_pooled_names <- function(x) grep("^\\.\\.(pooled|obs)_", names(x), value = TRUE)

# 900 rows with a response and two predictors that all vary smoothly over the
# layer, in 150 cells: x, then y, then z, p1 and p2.  No column of it raises
# a warning of the pooled path.
.ps_layer <- function() .ps_keep("layer", {
  d <- local({
    .ps_seed(2001)
    x <- stats::runif(900, 0, 1e5); y <- stats::runif(900, 0, 1e5)
    data.frame(x = x, y = y,
               z  = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(900, 0, 0.6),
               p1 = cos(x / 1.2e4) - sin(y / 2.5e4) + stats::rnorm(900, 0, 0.5),
               p2 = (x + y) / 5e4 + stats::rnorm(900, 0, 0.4))
  })
  .ps_assign(.ps_points(d$x, d$y, d[c("z", "p1", "p2")]), 150L)
})

# The (3000, 4.5) layout of the engine's reference tests, assigned through
# real polygons: Voronoi cells on the seeds of the k-means kind from 10 starts.
.ps_reference <- function() .ps_keep("reference", {
  L   <- do.call(pool_layout, POOL_LAYOUTS[[4L]])
  pts <- .ps_points(L$xy[, 1], L$xy[, 2], data.frame(v = L$y))
  seeds <- voronoi_seeds_kmeans(pts, 667, set_seed = 456, nstart = 10, partition = "kmeans")
  tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  list(layout = L, cells = tess$cells,
       asg = assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id"))
})

# The lattice of the engine's property tests: 60 x 50 locations at a spacing
# of 0.7 from (0.1, 0.3), in the regions of a k-means fit at 667 centres.
# None of its coordinates is exact in binary and every location has
# neighbours at exactly equal distances.  The regions are an ID column: a row
# on a boundary between two polygons would be assigned by the tie rule of
# assign_features_to_polygons(), which is not what is tested here.
.ps_lattice <- function() .ps_keep("lattice", {
  xy <- as.matrix(expand.grid(x = seq(0.1, by = 0.7, length.out = 60),
                              y = seq(0.3, by = 0.7, length.out = 50)))
  dimnames(xy) <- NULL
  v <- local({
    .ps_seed(31)
    sin(xy[, 1] / 9) + cos(xy[, 2] / 7) + stats::rnorm(nrow(xy), 0, 0.5)
  })
  g <- local({
    .ps_seed(3)
    km <- suppressWarnings(stats::kmeans(xy, 667L, nstart = 3, iter.max = 100))
    .nearest_centre(xy[, 1], xy[, 2], km$centers[, 1], km$centers[, 2])
  })
  .ps_points(xy[, 1], xy[, 2], data.frame(v = v, cell_id = match(g, sort(unique(g)))))
})

# 900 rows on a 10 x 10 grid of squares, the cell an arithmetic function of
# the coordinates, so that the cells are the same on every platform.
.ps_grid_layer <- function() .ps_keep("grid", {
  d <- local({
    .ps_seed(2002)
    x <- stats::runif(900, 0, 1e5); y <- stats::runif(900, 0, 1e5)
    data.frame(x = x, y = y,
               income = 52000 + 9000 * sin(x / 2e4) + 7000 * cos(y / 1.5e4) +
                 stats::rnorm(900, 0, 6000))
  })
  d$cell_id <- as.integer(10 * floor(d$x / 1e4) + floor(d$y / 1e4) + 1)
  .ps_points(d$x, d$y, d[c("income", "cell_id")])
})

# 400 rows of a Gaussian field with an exponential covariance and a nugget,
# on which a variogram is identified, in 60 cells: x, then y, then z, p1, p2.
.ps_field_layer <- function() .ps_keep("field", {
  d <- local({
    .ps_seed(2010)
    x <- 5e5 + stats::runif(400, 0, 1000); y <- 5e6 + stats::runif(400, 0, 1000)
    D <- as.matrix(stats::dist(cbind(x, y)))
    z <- 10 + as.numeric(t(chol(exp(-D / 90) + diag(0.3, 400))) %*% stats::rnorm(400))
    data.frame(x = x, y = y, z = z, p1 = 0.5 * z + stats::rnorm(400), p2 = stats::rnorm(400))
  })
  .ps_assign(.ps_points(d$x, d$y, d[c("z", "p1", "p2")], crs = 32632), 60L)
})


# ---- The argument omitted, and the columns the two estimates share -----------

test_that("estimate = \"pooled\" named is the call without the argument", {
  L <- .ps_layer()
  calls <- list(
    function(...) summarize_by_cell(L$asg, "z", id_col = "cell_id", ...),
    function(...) summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id",
                                    deff = "kish", conf_level = 0.9, ...),
    function(...) summarize_by_cell(L$asg, "z", id_col = "cell_id", cells_sf = L$cells,
                                    area = TRUE, deff = 2, ...),
    function(...) summarize_by_cell(L$asg, "z", id_col = "cell_id", agg_funs = median, ...),
    function(...) summarize_by_cell(sf::st_drop_geometry(L$asg), "z", "p1",
                                    id_col = "cell_id", ...))
  for (f in calls) {
    # The default is the pooled value beside the plain mean, with the same
    # warnings as when it is named (a layer with no geometry has none to
    # pool by).
    plain <- .ps_catch(f())
    named <- .ps_catch(f(estimate = "pooled"))
    expect_identical(named$value, plain$value)
    expect_identical(.ps_classes(named), .ps_classes(plain))
    expect_identical(.ps_catch(f(estimate = "pooled", obs_level = NULL))$value, plain$value)
    expect_false(is.null(attr(plain$value, "pooled")))
    # No attribute of the pooled path on the frame of plain means.
    expect_null(attr(f(estimate = "mean"), "pooled"))
    expect_null(attr(f(estimate = "mean"), "intervals"))
    expect_length(.ps_pooled_names(f(estimate = "mean")), 0L)
  }
  # The two new arguments come last, so a call by position is as it was.
  expect_identical(utils::tail(names(formals(summarize_by_cell)), 2L),
                   c("estimate", "obs_level"))
  # Which of the two `estimate` defaults to is held in test-core-defaults.R;
  # the test above (the argument named against the argument omitted) is the
  # one expectation of this file that reads the default.
  expect_setequal(eval(formals(summarize_by_cell)$estimate), c("mean", "pooled"))
  expect_null(formals(summarize_by_cell)$obs_level)
  expect_identical(
    summarize_by_cell(L$asg, "z", NULL, "cell_id", list(mean = function(x) mean(x, na.rm = TRUE)),
                      NULL, 1, NULL, 500L, TRUE, 0.9, FALSE),
    summarize_by_cell(L$asg, "z", id_col = "cell_id", conf_level = 0.9))
  # A layer with structure in every column raises nothing on the pooled path,
  # and the call draws no random number.
  .ps_seed(77)
  before <- get(".Random.seed", envir = globalenv())
  r <- .ps_catch(summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id",
                                   estimate = "pooled", obs_level = 0.9))
  expect_length(r$warnings, 0L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_error(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "median"),
               "should be one of")
})

test_that("every column and attribute the two estimates share is identical", {
  settings <- list(`deff = 1` = list(), `deff = "kish"` = list(deff = "kish"),
                   `deff = 2` = list(deff = 2))
  layers <- list(`deff = 1` = .ps_layer(), `deff = "kish"` = .ps_layer(),
                 `deff = 2` = .ps_layer())
  if (requireNamespace("gstat", quietly = TRUE)) {
    # A variogram of the response, on a layer where one is identified.
    layers$`a supplied sac` <- .ps_field_layer()
    sac <- suppressWarnings(estimate_sac_range(layers$`a supplied sac`$asg, "z"))
    expect_true(is.finite(as.numeric(sac)))
    settings$`a supplied sac` <- list(deff = "variogram", sac = sac)
  }
  compared <- 0L
  for (nm in names(settings)) for (conf in list(NULL, 0.95))
    for (shape in c("frame", "cells", "area")) {
      L <- layers[[nm]]
      args <- c(list(L$asg, "z", c("p1", "p2"), id_col = "cell_id", conf_level = conf),
                settings[[nm]],
                switch(shape, frame = list(), cells = list(cells_sf = L$cells),
                       area = list(cells_sf = L$cells, area = TRUE)))
      tag <- sprintf("%s, conf_level %s, %s", nm, format(conf), shape)
      a <- suppressWarnings(do.call(summarize_by_cell, c(args, list(estimate = "mean"))))
      b <- suppressWarnings(do.call(summarize_by_cell, c(args, list(estimate = "pooled"))))
      new <- .ps_pooled_names(b)
      expect_length(new, 9L)
      # Every present column keeps its name, its values and its position; the
      # new ones stand together behind the grouped summary's.
      expect_identical(setdiff(names(b), new), names(a), info = tag)
      for (col in names(a)) expect_identical(b[[col]], a[[col]], info = paste(tag, col))
      at <- match(new, names(b))
      expect_identical(at, seq(at[1L], length.out = 9L), info = tag)
      last <- if (!identical(settings[[nm]], list())) "deff_method" else "cell_weight"
      expect_identical(names(b)[at[1L] - 1L], last, info = tag)
      if (shape != "frame")
        expect_identical(names(b)[at[9L] + 1L], attr(b, "sf_column"), info = tag)
      for (att in c("deff_applied", "icc", "area_error"))
        expect_identical(attr(b, att), attr(a, att), info = paste(tag, att))
      expect_identical(setdiff(names(attributes(b)), names(attributes(a))),
                       c("pooled", "intervals"), info = tag)
      expect_identical(class(b), class(a), info = tag)
      compared <- compared + 1L
    }
  expect_gte(compared, 18L)
  # The design effect was applied where it was asked for, so the comparison
  # above is of corrected columns and not of two frames without any.
  k <- suppressWarnings(summarize_by_cell(.ps_layer()$asg, "z", id_col = "cell_id",
                                          deff = "kish", estimate = "pooled"))
  expect_true(all(k$deff_applied))
  expect_identical(attr(k, "deff_applied")$method, "kish")
  if (requireNamespace("gstat", quietly = TRUE)) {
    v <- suppressWarnings(summarize_by_cell(.ps_field_layer()$asg, "z", id_col = "cell_id",
                                            deff = "variogram", sac = sac, estimate = "pooled"))
    expect_identical(unique(v$deff_method), "variogram")
    expect_identical(attr(v, "deff_applied")$method, "variogram")
  }
})

test_that("deff does not enter the pooled columns", {
  L <- .ps_layer()
  one <- summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id", deff = 1,
                           estimate = "pooled", obs_level = 0.9)
  for (d in list("kish", 3)) {
    other <- summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id", deff = d,
                               estimate = "pooled", obs_level = 0.9, conf_level = 0.95)
    for (col in .ps_pooled_names(one))
      expect_identical(other[[col]], one[[col]], info = paste(format(d), col))
    expect_identical(attr(other, "pooled")$columns, attr(one, "pooled")$columns)
    # The standard errors did move.
    expect_false(isTRUE(all.equal(other$..se_resp_z, one$..se_resp_z)))
  }
})


# ---- Reproduction ------------------------------------------------------------

test_that("the pooled columns and the record are the reference functions' on real polygons", {
  skip_if_not_installed("FNN")
  R <- .ps_reference()
  asg <- R$asg
  expect_identical(nrow(asg), 3000L)
  expect_true(is.integer(asg$cell_id))
  xy <- unname(sf::st_coordinates(asg)[, 1:2])
  expect_identical(xy, R$layout$xy)                 # the rows are in fixture order
  y  <- asg$v
  g  <- asg$cell_id
  u  <- sort(unique(g)); K <- length(u); gi <- match(g, u)
  expect_gt(K, 600L)
  ref <- region_estimates(y, gi, K, xy, want = c("P", "L"))
  fr  <- fit_regions(g, g, y, xy, xy, pool_nearest_other(xy))
  li  <- loo_int(unname(ref$est[, "L"]), fr$eL, 1, 1)
  N   <- length(y)

  for (fnn in c(TRUE, FALSE)) {
    tag <- if (fnn) "FNN" else "exact"
    out <- local({
      if (!fnn) .ps_without_fnn()
      r <- .ps_catch(summarize_by_cell(asg, "v", id_col = "cell_id", estimate = "pooled",
                                       obs_level = 0.90))
      expect_length(r$warnings, 0L)
      r$value
    })
    expect_identical(out$cell_id, u, info = tag)
    expect_identical(out$..pooled_resp_v, unname(ref$est[, "L"]), info = tag)
    expect_identical(out$..pooled_w_resp_v, fr$w, info = tag)
    expect_identical(out$..pooled_target_resp_v, fr$t, info = tag)
    expect_identical(out$..obs_lo_resp_v, li[, 1], info = tag)
    expect_identical(out$..obs_hi_resp_v, li[, 2], info = tag)
    rec <- attr(out, "pooled")$columns$resp_v
    expect_identical(rec$loo, fr$loo_L, info = tag)
    expect_identical(rec$loo_se, fr$loo_L_se, info = tag)
    expect_identical(rec$own_weight, stats::median(fr$w), info = tag)
    expect_identical(rec$share_alone, sum(fr$n == 1L) / N, info = tag)
    expect_identical(rec$df, N - K, info = tag)
    expect_identical(rec$loo_one_cell, mean(((y - mean(y)) * N / (N - 1))^2), info = tag)
    expect_identical(rec$sigma2, fr$sigma2, info = tag)
    expect_identical(rec$tau2, fr$tau2, info = tag)
    expect_identical(rec$n, N, info = tag)
    expect_identical(rec$cells, K, info = tag)
    expect_identical(rec$offsets,
                     stats::quantile(fr$eL, c(0.05, 0.95), names = FALSE, type = 7), info = tag)
    expect_identical(rec$status, "ok", info = tag)
  }
  # The layout has cells of one row, where the weight is smallest.
  expect_gt(sum(fr$n == 1L), 20L)

  # In another row order and with IDs that are not whole numbers the sums are
  # taken in another order, and the equality is to rounding.
  p   <- local({ .ps_seed(41); sample.int(N) })
  lab <- local({ .ps_seed(42); sample(sprintf("cell-%04d", seq_len(K))) })
  alt <- asg[p, ]
  alt$cell_id <- lab[gi[p]]
  out <- summarize_by_cell(alt, "v", id_col = "cell_id", estimate = "pooled", obs_level = 0.90)
  at  <- match(lab, out$cell_id)
  expect_false(is.unsorted(out$cell_id))
  expect_false(identical(at, seq_len(K)))
  expect_equal(out$..pooled_resp_v[at], unname(ref$est[, "L"]), tolerance = 1e-12)
  expect_equal(out$..pooled_w_resp_v[at], fr$w, tolerance = 1e-12)
  expect_equal(out$..pooled_target_resp_v[at], fr$t, tolerance = 1e-12)
  expect_equal(out$..obs_lo_resp_v[at], li[, 1], tolerance = 1e-12)
  expect_equal(out$..obs_hi_resp_v[at], li[, 2], tolerance = 1e-12)
  rec <- attr(out, "pooled")$columns$resp_v
  expect_equal(rec$loo, fr$loo_L, tolerance = 1e-12)
  expect_equal(rec$own_weight, stats::median(fr$w), tolerance = 1e-12)
  expect_identical(rec$share_alone, sum(fr$n == 1L) / N)
  expect_identical(rec$df, N - K)
})

test_that("the pooled value is its weight times the plain mean and the rest times the target", {
  L <- .ps_layer()
  out <- summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id", estimate = "pooled")
  for (v in c("resp_z", "pred_p1", "pred_p2")) {
    w <- out[[paste0("..pooled_w_", v)]]
    m <- out[[sub("_", "_mean_", v, fixed = TRUE)]]
    expect_equal(out[[paste0("..pooled_", v)]],
                 w * m + (1 - w) * out[[paste0("..pooled_target_", v)]],
                 tolerance = 1e-12, info = v)
    expect_true(all(w > 0 & w < 1), info = v)
    # The weights of one column rise with the cell's count.
    expect_false(is.unsorted(w[order(out$n)]), info = v)
  }
})

test_that("a valued row without a location is left out of the pooled mean, with a warning", {
  L <- .ps_layer()
  asg <- L$asg
  # A row of a cell of four rows or more loses its point; its value is the
  # cell's largest by far, so the two means differ visibly.
  cell <- as.integer(names(which(table(asg$cell_id) >= 4L))[1L])
  row  <- which(asg$cell_id == cell)[1L]
  asg$z[row] <- 50
  g <- sf::st_geometry(asg)
  g[[row]] <- sf::st_point()
  sf::st_geometry(asg) <- g
  r <- .ps_catch(summarize_by_cell(asg, "z", id_col = "cell_id", estimate = "pooled"))
  expect_identical(.ps_classes(r), "spatialkit_pool_rows_skipped")
  w <- r$warnings[[1L]]
  expect_identical(.ps_fields(w), list(n = 1L, cells = 1L, columns = "z"))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): 1 row(s) in 1 cell(s) have empty or non-finite coordinates. ",
           "They count in `n` and the plain columns and take no part in the pooled columns."))
  out <- r$value
  rec <- attr(out, "pooled")$columns$resp_z
  expect_identical(rec$skipped, 1L)
  expect_identical(rec$n, nrow(asg) - 1L)
  i <- which(out$cell_id == cell)
  expect_identical(out$n[i], sum(asg$cell_id == cell))          # the row still counts in n
  identity_gap <- abs(out$..pooled_resp_z -
                        (out$..pooled_w_resp_z * out$resp_mean_z +
                           (1 - out$..pooled_w_resp_z) * out$..pooled_target_resp_z))
  # In that cell the plain mean holds the row and the pooled value's own mean
  # does not, so the identity fails there.  Its neighbours' targets changed
  # with it, and everywhere else the identity holds.
  expect_gt(identity_gap[i], 1)
  expect_lt(max(identity_gap[-i]), 1e-12)
  # The pooled columns are those of the layer without the row.
  less <- summarize_by_cell(asg[-row, ], "z", id_col = "cell_id", estimate = "pooled")
  for (col in .ps_pooled_names(out)) expect_identical(out[[col]], less[[col]], info = col)
  rec$skipped <- 0L
  expect_identical(rec, attr(less, "pooled")$columns$resp_z)
})


# ---- Columns, attributes and the join -----------------------------------------

test_that("the new columns follow the present ones, the response's first", {
  L <- .ps_layer()
  out <- summarize_by_cell(L$asg, "z", c("p2", "p1"), id_col = "cell_id", deff = "kish",
                           conf_level = 0.95, estimate = "pooled", obs_level = 0.9)
  stems <- c("..pooled_", "..pooled_w_", "..pooled_target_", "..obs_lo_", "..obs_hi_")
  # Behind every present column, the design effect's two included; then the
  # response's, then each predictor's in the order of `predictor_vars`.
  expect_identical(utils::tail(names(out), 18L),
                   c("cell_weight", "deff_applied", "deff_method", paste0(stems, "resp_z"),
                     paste0(stems, "pred_p2"), paste0(stems, "pred_p1")))
  # Without a level, three columns a variable; without a design effect they
  # follow cell_weight.
  out <- summarize_by_cell(L$asg, "z", "p1", id_col = "cell_id", estimate = "pooled")
  expect_identical(utils::tail(names(out), 7L),
                   c("cell_weight", paste0(stems[1:3], "resp_z"), paste0(stems[1:3], "pred_p1")))
  # Predictors alone.
  out <- summarize_by_cell(L$asg, predictor_vars = c("p1", "p2"), id_col = "cell_id",
                           estimate = "pooled")
  expect_identical(.ps_pooled_names(out),
                   c(paste0(stems[1:3], "pred_p1"), paste0(stems[1:3], "pred_p2")))
  expect_identical(names(attr(out, "pooled")$columns), c("pred_p1", "pred_p2"))
  # With the cells joined, the geometry and the area follow.
  out <- summarize_by_cell(L$asg, "z", id_col = "cell_id", cells_sf = L$cells, area = TRUE,
                           estimate = "pooled", obs_level = 0.9)
  expect_identical(utils::tail(names(out), 9L),
                   c("cell_weight", paste0(stems, "resp_z"), "geometry", "cell_area",
                     "n_per_area"))
})

test_that("the record has the fields of the design, and each column its entry", {
  L <- .ps_layer()
  out <- summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id", conf_level = 0.95,
                           estimate = "pooled", obs_level = 0.9)
  rec <- attr(out, "pooled")
  expect_s3_class(rec, "spatialkit_pooled")
  expect_identical(names(rec), c("toward", "neighbours", "crs", "obs_level", "conf_level",
                                 "intervals", "columns"))
  expect_identical(rec$toward, "neighbours")
  expect_identical(rec$neighbours, 8L)
  expect_identical(rec$crs, sf::st_crs(L$asg))
  expect_identical(rec$obs_level, 0.9)
  expect_identical(rec$conf_level, 0.95)
  expect_identical(names(rec$columns), c("resp_z", "pred_p1", "pred_p2"))
  fields <- c("column", "role", "n", "cells", "df", "neighbours", "sigma2", "tau2",
              "own_weight", "share_alone", "loo", "loo_se", "loo_one_cell", "structured",
              "structure_p", "offsets", "repeated", "skipped", "constant", "status")
  for (nm in names(rec$columns)) {
    r <- rec$columns[[nm]]
    expect_identical(names(r), fields, info = nm)
    expect_identical(r$column, sub("^(resp|pred)_", "", nm), info = nm)
    expect_identical(r$role, if (startsWith(nm, "resp_")) "response" else "predictor", info = nm)
    expect_identical(r$n, 900L, info = nm)
    expect_identical(r$cells, length(unique(L$asg$cell_id)), info = nm)
    expect_identical(r$df, r$n - r$cells, info = nm)
    expect_identical(r$neighbours, 8L, info = nm)
    expect_identical(r$own_weight, stats::median(out[[paste0("..pooled_w_", nm)]]), info = nm)
    expect_identical(r$share_alone, sum(out$n == 1L) / 900, info = nm)
    expect_true(r$structured, info = nm)
    expect_lt(r$structure_p, 0.05)
    expect_lt(r$loo, r$loo_one_cell)
    expect_gt(r$loo_se, 0)
    expect_length(r$offsets, 2L)
    expect_identical(out[[paste0("..obs_lo_", nm)]], out[[paste0("..pooled_", nm)]] + r$offsets[1L])
    expect_identical(out[[paste0("..obs_hi_", nm)]], out[[paste0("..pooled_", nm)]] + r$offsets[2L])
    expect_identical(c(r$repeated, r$skipped), c(0L, 0L), info = nm)
    expect_false(r$constant, info = nm)
    expect_identical(r$status, "ok", info = nm)
  }
  # Without the levels the fields are there and NULL.
  bare <- attr(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled"), "pooled")
  expect_identical(names(bare), names(rec))
  expect_null(bare$obs_level)
  expect_null(bare$conf_level)
  expect_identical(names(bare$columns$resp_z), fields)
  expect_null(bare$columns$resp_z$offsets)
  # The engine's own fit of the same rows is what the record and the frame hold.
  xy <- unname(sf::st_coordinates(L$asg)[, 1:2])
  u  <- sort(unique(L$asg$cell_id))
  f  <- .pool_fit(L$asg$z, match(L$asg$cell_id, u), length(u), xy, want = "pooled")
  expect_identical(out$..pooled_resp_z, f$value)
  expect_identical(rec$columns$resp_z$loo, f$loo_pooled)
  expect_identical(rec$columns$resp_z$loo_se, f$loo_pooled_se)
  expect_identical(rec$columns$resp_z$loo_one_cell, f$loo_one)
  expect_identical(rec$columns$resp_z$structure_p, f$f_p)
})

test_that("the attributes survive the join, and a cell without a row is NA", {
  L <- .ps_layer()
  ids  <- sort(unique(L$asg$cell_id))
  gone <- ids[c(3L, 77L)]
  asg  <- L$asg[!(L$asg$cell_id %in% gone), ]
  plain  <- summarize_by_cell(asg, "z", "p1", id_col = "cell_id", estimate = "pooled",
                              obs_level = 0.9, conf_level = 0.95, deff = "kish")
  joined <- summarize_by_cell(asg, "z", "p1", id_col = "cell_id", estimate = "pooled",
                              obs_level = 0.9, conf_level = 0.95, deff = "kish",
                              cells_sf = L$cells)
  expect_s3_class(joined, "sf")
  expect_identical(nrow(joined), nrow(L$cells))
  expect_identical(attr(joined, "pooled"), attr(plain, "pooled"))
  expect_identical(attr(joined, "intervals"), attr(plain, "intervals"))
  expect_identical(attr(joined, "icc"), attr(plain, "icc"))
  new <- .ps_pooled_names(plain)
  expect_length(new, 10L)
  empty <- joined$cell_id %in% gone
  expect_identical(sum(empty), 2L)
  for (col in new) {
    expect_true(all(is.na(joined[[col]][empty])), info = col)
    at <- match(plain$cell_id, joined$cell_id)
    expect_identical(joined[[col]][at], plain[[col]], info = col)
    expect_false(anyNA(joined[[col]][!empty]), info = col)
  }
  # The empty cells are no cell's neighbour: the values are those of a layer
  # that never held them.
  expect_identical(attr(joined, "pooled")$columns$resp_z$cells, length(ids) - 2L)
  # The frame returned keeps the cells' CRS, and IDs of another class than the
  # cells' are joined as today, with the attributes kept.
  expect_identical(sf::st_crs(joined), sf::st_crs(L$cells))
  cells_chr <- L$cells
  cells_chr$cell_id <- as.character(cells_chr$cell_id)
  mixed <- summarize_by_cell(asg, "z", "p1", id_col = "cell_id", estimate = "pooled",
                             obs_level = 0.9, conf_level = 0.95, deff = "kish",
                             cells_sf = cells_chr)
  expect_type(mixed$cell_id, "character")
  expect_identical(attr(mixed, "pooled"), attr(plain, "pooled"))
  expect_identical(attr(mixed, "intervals"), attr(plain, "intervals"))
  expect_identical(mixed$..pooled_resp_z[match(as.character(plain$cell_id), mixed$cell_id)],
                   plain$..pooled_resp_z)
  # A `cells_sf` that cannot be joined leaves the frame, with its attributes.
  r <- .ps_catch(summarize_by_cell(asg, "z", id_col = "cell_id", estimate = "pooled",
                                   cells_sf = "not a layer"))
  expect_length(r$warnings, 1L)
  expect_s3_class(attr(r$value, "pooled"), "spatialkit_pooled")
  expect_identical(attr(r$value, "intervals"), attr(r$value, "pooled")$intervals)
})


# ---- Missing and infinite values ----------------------------------------------

test_that("a column is pooled on its own rows: the rows that miss it are as if absent", {
  L <- .ps_layer()
  asg <- L$asg
  miss1 <- local({ .ps_seed(51); sample.int(nrow(asg), 60L) })
  miss2 <- local({ .ps_seed(52); sample.int(nrow(asg), 45L) })
  both_miss <- length(intersect(miss1, miss2))      # they overlap in a few rows only
  expect_gt(both_miss, 0L); expect_lt(both_miss, 10L)
  asg$p1[miss1] <- NA
  asg$p2[miss2] <- NA
  # A cell whose rows all miss p1, so the column has a cell fewer than the others.
  cell <- as.integer(names(which(table(asg$cell_id) == 3L))[1L])
  asg$p1[asg$cell_id == cell] <- NA
  both <- summarize_by_cell(asg, "z", c("p1", "p2"), id_col = "cell_id", estimate = "pooled",
                            obs_level = 0.9)
  for (v in c("p1", "p2")) {
    keep <- !is.na(asg[[v]])
    own  <- summarize_by_cell(asg[keep, ], predictor_vars = v, id_col = "cell_id",
                              estimate = "pooled", obs_level = 0.9)
    at <- match(own$cell_id, both$cell_id)
    for (col in .ps_pooled_names(own))
      expect_identical(both[[col]][at], own[[col]], info = col)
    expect_identical(attr(both, "pooled")$columns[[paste0("pred_", v)]],
                     attr(own, "pooled")$columns[[paste0("pred_", v)]])
    expect_identical(attr(own, "pooled")$columns[[paste0("pred_", v)]]$n, sum(keep))
  }
  # The response misses nothing and is as it is without the predictors.
  alone <- summarize_by_cell(asg, "z", id_col = "cell_id", estimate = "pooled", obs_level = 0.9)
  for (col in .ps_pooled_names(alone)) expect_identical(both[[col]], alone[[col]], info = col)
  # The cell whose rows all miss p1 has no pooled value of it, and one of the others.
  i <- which(both$cell_id == cell)
  expect_true(all(is.na(unlist(sf::st_drop_geometry(both)[i, grep("pred_p1$", .ps_pooled_names(both), value = TRUE)]))))
  expect_false(is.na(both$..pooled_pred_p2[i]))
  expect_identical(attr(both, "pooled")$columns$pred_p1$cells,
                   attr(both, "pooled")$columns$resp_z$cells - 1L)
})

test_that("an infinite value is a missing value, and is said once", {
  L <- .ps_layer()
  asg <- L$asg
  asg$z[c(5L, 300L)] <- c(Inf, -Inf)
  asg$p1[17L] <- Inf
  na <- asg
  na$z[c(5L, 300L)] <- NA
  na$p1[17L] <- NA
  r <- .ps_catch(summarize_by_cell(asg, "z", "p1", id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.9))
  expect_identical(.ps_classes(r), "spatialkit_nonfinite")
  expect_identical(r$warnings[[1L]]$columns, c("z", "p1"))
  expect_identical(r$warnings[[1L]]$n, 3L)
  # The same condition as the call of plain means raises.
  m <- .ps_catch(summarize_by_cell(asg, "z", "p1", id_col = "cell_id", estimate = "mean"))
  expect_identical(conditionMessage(m$warnings[[1L]]), conditionMessage(r$warnings[[1L]]))
  want <- summarize_by_cell(na, "z", "p1", id_col = "cell_id", estimate = "pooled",
                            obs_level = 0.9)
  expect_identical(r$value, want)
  expect_identical(attr(r$value, "pooled")$columns$resp_z$n, 898L)
})


# ---- One warning per call and reason -------------------------------------------

test_that("three structureless columns raise one spatialkit_no_structure that names all three", {
  # White noise in the response and both predictors: x, then y, then the three
  # columns, 3,000 rows in 667 cells.
  d <- local({
    .ps_seed(6101)
    x <- stats::runif(3000, 0, 1e5); y <- stats::runif(3000, 0, 1e5)
    data.frame(x = x, y = y, a = stats::rnorm(3000), b = stats::rnorm(3000),
               c = stats::rnorm(3000))
  })
  A <- .ps_assign(.ps_points(d$x, d$y, d[c("a", "b", "c")]), 667L)
  r <- .ps_catch(summarize_by_cell(A$asg, "b", c("c", "a"), id_col = "cell_id",
                                   estimate = "pooled"))
  expect_identical(.ps_classes(r), "spatialkit_no_structure")
  w <- r$warnings[[1L]]
  rec <- attr(r$value, "pooled")$columns
  f <- .ps_fields(w)
  expect_identical(names(f), c("criterion", "n", "one_cell", "columns", "estimate", "loo",
                               "structure_p", "beats"))
  expect_identical(f$columns, c("b", "c", "a"))             # the order of the frame
  expect_identical(f$criterion, "loo")
  expect_identical(f$estimate, "pooled")
  by <- function(field) stats::setNames(
    vapply(rec, function(r) as.numeric(r[[field]]), numeric(1)), c("b", "c", "a"))
  expect_identical(f$n, stats::setNames(rep(rec$resp_b$cells, 3L), c("b", "c", "a")))
  expect_identical(f$one_cell, by("loo_one_cell"))
  expect_identical(f$loo, by("loo"))
  expect_identical(f$structure_p, by("structure_p"))
  expect_true(all(!vapply(rec, function(r) r$structured, NA)))
  # For each of the three one mean does as well as the cells, to the margin.
  expect_identical(f$beats, c(b = FALSE, c = FALSE, a = FALSE))
  for (r1 in rec) expect_lte(r1$loo_one_cell, r1$loo * 1.02)
  # Several columns: the sentences first, then the numbers of each column.
  piece <- function(nm) {
    r1 <- rec[[nm]]
    sprintf(paste0("'%s': leave-one-out error %s against %s; F test of the cells against one ",
                   "mean, p = %s; median own weight %.2f."), r1$column,
            format(signif(r1$loo, 4)), format(signif(r1$loo_one_cell, 4)),
            format(signif(r1$structure_p, 2)), r1$own_weight)
  }
  expect_identical(
    conditionMessage(w),
    paste(paste0("summarize_by_cell(): the pooled values of 'b', 'c' and 'a' predict a left-out ",
                 "observation no better than one mean for all rows, so they show no spatial ",
                 "structure these cells can use. The values are returned. A map of them looks ",
                 "smooth because each is mostly its neighbours' mean, not because the response ",
                 "varies smoothly."),
          piece("resp_b"), piece("pred_c"), piece("pred_a")))
  expect_null(conditionCall(w))
  # One column alone: the singular, and no name before its numbers.
  one <- .ps_catch(summarize_by_cell(A$asg, "a", id_col = "cell_id", estimate = "pooled"))
  expect_identical(.ps_classes(one), "spatialkit_no_structure")
  ra <- attr(one$value, "pooled")$columns$resp_a
  expect_identical(
    conditionMessage(one$warnings[[1L]]),
    sprintf(paste0(
      "summarize_by_cell(): the pooled values of 'a' predict a left-out observation no ",
      "better than one mean for all rows (leave-one-out error %s against %s; F test of the ",
      "cells against one mean, p = %s; median own weight %.2f), so it shows no spatial ",
      "structure these cells can use. The values are returned. A map of them looks smooth ",
      "because each is mostly its neighbours' mean, not because the response varies ",
      "smoothly."),
      format(signif(ra$loo, 4)), format(signif(ra$loo_one_cell, 4)),
      format(signif(ra$structure_p, 2)), ra$own_weight))
  # The values are returned.
  expect_false(anyNA(one$value$..pooled_resp_a))
  # Plain means of the same layer say nothing.
  expect_length(.ps_catch(summarize_by_cell(A$asg, "a", id_col = "cell_id",
                                            estimate = "mean"))$warnings, 0L)
})

test_that("a warning that lists many columns says what it has to say before their numbers", {
  # A response and ten predictors, all white noise, with long names: x, then
  # y, then the eleven columns; 1,200 rows in 100 cells.
  nm <- sprintf("a_long_column_name_%02d", 1:11)
  d <- local({
    .ps_seed(6102)
    x <- stats::runif(1200, 0, 1e5); y <- stats::runif(1200, 0, 1e5)
    cbind(data.frame(x = x, y = y),
          as.data.frame(stats::setNames(lapply(nm, function(i) stats::rnorm(1200)), nm)))
  })
  cell <- function(x, y) as.integer(10 * floor(x / 1e4) + floor(y / 1e4) + 1)
  pts <- .ps_points(d$x, d$y, d[nm])
  pts$cell_id <- cell(d$x, d$y)
  # R prints the first getOption("warning.length") characters of a warning.
  withr::local_options(warning.length = 1000L)
  r <- .ps_catch(summarize_by_cell(pts, nm[1L], nm[-1L], id_col = "cell_id",
                                   estimate = "pooled"))
  w <- .ps_of(r, "spatialkit_no_structure")
  expect_length(w, 1L)
  cols <- w[[1L]]$columns
  expect_gte(length(cols), 8L)
  # For each of them one mean does as well as the cells, to the margin of the
  # rule: for one at least its error is above the cells', by less than that.
  expect_true(all(!w[[1L]]$beats))
  rec <- attr(r$value, "pooled")$columns
  expect_true(any(vapply(rec, function(c1)
    isFALSE(c1$structured) && c1$loo_one_cell > c1$loo && c1$loo_one_cell <= 1.02 * c1$loo, NA)))
  msg <- conditionMessage(w[[1L]])
  # The whole text passes what R prints; the closing sentence does not.
  expect_gt(nchar(msg), 1000L)
  closing <- "not because the response varies smoothly."
  ends <- regexpr(closing, msg, fixed = TRUE) + nchar(closing) - 1L
  expect_gt(ends, 0L)
  expect_lt(ends, 600L)
  expect_true(grepl(closing, substr(msg, 1L, getOption("warning.length")), fixed = TRUE))
  expect_true(grepl("The values are returned.", substr(msg, 1L, ends), fixed = TRUE))
  # Six names and six columns' numbers, the numbers behind the sentences; the
  # condition holds every column.
  expect_match(msg, sprintf("'%s' and %d more predict", cols[6L], length(cols) - 6L), fixed = TRUE)
  at <- vapply(cols, function(v)
    regexpr(sprintf("'%s': leave-one-out error ", v), msg, fixed = TRUE)[1L], 1L)
  expect_true(all(at[1:6] > ends))
  expect_false(is.unsorted(at[1:6]))
  expect_true(all(at[-(1:6)] == -1L))
  expect_true(endsWith(msg, sprintf("And %d more, in attr(, \"pooled\").", length(cols) - 6L)))
  expect_identical(names(w[[1L]]$loo), cols)
  # With six columns every number is given, and nothing is said of more.
  six <- .ps_catch(summarize_by_cell(pts, cols[1L], cols[2:6], id_col = "cell_id",
                                     estimate = "pooled"))
  m6 <- conditionMessage(.ps_of(six, "spatialkit_no_structure")[[1L]])
  expect_identical(lengths(regmatches(m6, gregexpr("leave-one-out error ", m6, fixed = TRUE))), 6L)
  expect_false(grepl("more, in attr", m6, fixed = TRUE))
  expect_true(grepl(closing, substr(m6, 1L, getOption("warning.length")), fixed = TRUE))

  # The same for rows that share a location: every row of the layer twice.
  twice <- .ps_points(rep(d$x, 2L), rep(d$y, 2L), rbind(d[nm], d[nm]))
  twice$cell_id <- cell(rep(d$x, 2L), rep(d$y, 2L))
  r <- .ps_catch(summarize_by_cell(twice, nm[1L], nm[-1L], id_col = "cell_id",
                                   estimate = "pooled", obs_level = 0.9))
  w <- .ps_of(r, "spatialkit_repeated_locations")
  expect_length(w, 1L)
  expect_identical(w[[1L]]$columns, nm)
  expect_identical(w[[1L]]$n, stats::setNames(rep(2400L, 11L), nm))
  expect_identical(
    conditionMessage(w[[1L]]),
    paste0("summarize_by_cell(): rows with a value share a location with another row (",
           paste(sprintf("2400 of the 2400 rows of '%s'", nm[1:6]), collapse = "; "),
           "; and 5 more). The leave-one-out error leaves one row out and keeps the others at ",
           "its location, so it, and an interval built from it, describes a further observation ",
           "at a place already observed, not an observation at a new place."))
  expect_lt(nchar(conditionMessage(w[[1L]])), 1000L)
})

test_that("ten constant predictors raise one spatialkit_pool_degenerate, six of them named", {
  L <- .ps_layer()
  asg <- L$asg
  consts <- c(0.1, 1234.56, 1 / 3, 52000.37, 0, 0.7, -2.5, 1e-9, 3.5, 99)
  names(consts) <- sprintf("k%02d", 1:10)
  for (nm in names(consts)) asg[[nm]] <- consts[[nm]]
  r <- .ps_catch(summarize_by_cell(asg, "z", names(consts), id_col = "cell_id",
                                   estimate = "pooled"))
  expect_identical(.ps_classes(r), "spatialkit_pool_degenerate")
  w <- r$warnings[[1L]]
  expect_identical(
    .ps_fields(w),
    list(reason = "no spread inside any cell", columns = names(consts),
         constant = stats::setNames(rep(TRUE, 10L), names(consts))))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): the values of 'k01', 'k02', 'k03', 'k04', 'k05', 'k06' and 4 ",
           "more do not vary inside any cell, so every own weight is 1 and the pooled value is ",
           "the plain mean. attr(, \"pooled\") says which columns are constant."))
  for (nm in names(consts)) {
    expect_identical(r$value[[paste0("..pooled_w_pred_", nm)]], rep(1, nrow(r$value)), info = nm)
    expect_equal(r$value[[paste0("..pooled_pred_", nm)]], rep(consts[[nm]], nrow(r$value)),
                 tolerance = 1e-12, info = nm)
    rec <- attr(r$value, "pooled")$columns[[paste0("pred_", nm)]]
    expect_identical(rec$status, "no spread inside any cell", info = nm)
    expect_true(rec$constant, info = nm)
    expect_identical(c(rec$sigma2, rec$tau2), c(0, 0), info = nm)
    expect_identical(rec$loo, NA_real_, info = nm)
  }
})

test_that("stations visited several times raise one spatialkit_repeated_locations", {
  # 300 stations visited three times: 900 rows on 300 locations, three columns.
  d <- local({
    .ps_seed(2003)
    x <- stats::runif(300, 0, 1e5); y <- stats::runif(300, 0, 1e5)
    x <- rep(x, 3L); y <- rep(y, 3L)
    data.frame(x = x, y = y,
               a = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(900, 0, 0.5),
               b = cos(x / 1.2e4) + stats::rnorm(900, 0, 0.5),
               c = (x - y) / 4e4 + stats::rnorm(900, 0, 0.5))
  })
  A <- .ps_assign(.ps_points(d$x, d$y, d[c("a", "b", "c")]), 60L)
  # The third column misses two visits of each of two stations, so each has
  # one row left that shares its location with no other row of the column.
  A$asg$c[c(1L, 301L, 2L, 302L)] <- NA
  r <- .ps_catch(summarize_by_cell(A$asg, "a", c("b", "c"), id_col = "cell_id",
                                   estimate = "pooled", obs_level = 0.9))
  expect_identical(.ps_classes(r), "spatialkit_repeated_locations")
  w <- r$warnings[[1L]]
  rep_c <- attr(r$value, "pooled")$columns$pred_c$repeated
  expect_identical(
    .ps_fields(w),
    list(columns = c("a", "b", "c"), n = c(a = 900L, b = 900L, c = rep_c),
         total = c(a = 900L, b = 900L, c = 896L)))
  expect_identical(rep_c, 894L)
  expect_identical(
    conditionMessage(w),
    sprintf(paste0(
      "summarize_by_cell(): rows with a value share a location with another row (900 of the ",
      "900 rows of 'a'; 900 of the 900 rows of 'b'; %d of the 896 rows of 'c'). The ",
      "leave-one-out error leaves one row out and keeps the others at its location, so it, ",
      "and an interval built from it, describes a further observation at a place already ",
      "observed, not an observation at a new place."), rep_c))
  # The record counts them whether or not an interval is asked for; the warning
  # is the interval's.
  quiet <- .ps_catch(summarize_by_cell(A$asg, "a", c("b", "c"), id_col = "cell_id",
                                       estimate = "pooled"))
  expect_length(quiet$warnings, 0L)
  expect_identical(attr(quiet$value, "pooled")$columns$resp_a$repeated, 900L)
  # Each row counts in its cell, and the pooled values are the engine's.
  expect_identical(sum(r$value$n), 900L)
})

test_that("two reasons for NA columns raise two spatialkit_pool_unavailable", {
  # `a` has one value in each cell; `b` has values only on rows without a point.
  d <- local({
    .ps_seed(2004)
    x <- stats::runif(120, 0, 1e5); y <- stats::runif(120, 0, 1e5)
    data.frame(x = x, y = y, z = sin(x / 2e4) + stats::rnorm(120, 0, 0.3),
               cell_id = rep(1:30, 4L))
  })
  asg <- .ps_points(d$x, d$y, d[c("z", "cell_id")])
  asg$a <- ifelse(seq_len(120) <= 30L, asg$z, NA_real_)
  asg$b <- NA_real_
  nowhere <- c(40L, 41L, 99L)
  asg$b[nowhere] <- c(1, 2, 3)
  g <- sf::st_geometry(asg)
  for (i in nowhere) g[[i]] <- sf::st_point()
  sf::st_geometry(asg) <- g
  r <- .ps_catch(summarize_by_cell(asg, "a", "b", id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.9))
  expect_identical(.ps_classes(r), c("spatialkit_pool_unavailable", "spatialkit_pool_unavailable",
                                     "spatialkit_pool_rows_skipped"))
  # In the order of the design's table: no located value, then no cell with two.
  w1 <- r$warnings[[1L]]; w2 <- r$warnings[[2L]]
  expect_identical(.ps_fields(w1), list(reason = "no located value", columns = "b"))
  expect_identical(
    conditionMessage(w1),
    paste0("summarize_by_cell(): the pooled columns of 'b' are NA: no row has a value of it, ",
           "a cell and a location."))
  expect_identical(.ps_fields(w2),
                   list(reason = "no cell holds two values", columns = "a", n = c(a = 30L),
                        cells = c(a = 30L)))
  expect_identical(
    conditionMessage(w2),
    paste0("summarize_by_cell(): the pooled columns of 'a' are NA: no cell holds two values of ",
           "it, so the spread of values inside a cell cannot be estimated and nothing can set a ",
           "weight. The plain columns are unaffected."))
  expect_identical(.ps_fields(r$warnings[[3L]]), list(n = 3L, cells = 3L, columns = "b"))
  out <- r$value
  for (col in .ps_pooled_names(out)) expect_true(all(is.na(out[[col]])), info = col)
  expect_length(.ps_pooled_names(out), 10L)
  rec <- attr(out, "pooled")$columns
  expect_identical(rec$resp_a$status, "no cell holds two values")
  expect_identical(rec$pred_b$status, "no located value")
  expect_identical(c(rec$pred_b$n, rec$pred_b$cells, rec$pred_b$skipped), c(0L, 0L, 3L))
  # The plain columns are there.
  expect_false(anyNA(out$resp_mean_a))
  expect_identical(sum(!is.na(out$pred_mean_b)), 3L)
  # Two columns with one reason are named in one warning, in the plural.
  asg$a2 <- asg$a
  two <- .ps_catch(summarize_by_cell(asg, "a", "a2", id_col = "cell_id", estimate = "pooled"))
  expect_identical(.ps_classes(two), "spatialkit_pool_unavailable")
  expect_identical(two$warnings[[1L]]$columns, c("a", "a2"))
  expect_identical(
    conditionMessage(two$warnings[[1L]]),
    paste0("summarize_by_cell(): the pooled columns of 'a' and 'a2' are NA: for none of them ",
           "does any cell hold two values, so the spread of values inside a cell cannot be ",
           "estimated and nothing can set a weight. The plain columns are unaffected."))
  asg$b2 <- asg$b
  two <- .ps_catch(summarize_by_cell(asg, "b", "b2", id_col = "cell_id", estimate = "pooled"))
  expect_identical(
    conditionMessage(two$warnings[[1L]]),
    paste0("summarize_by_cell(): the pooled columns of 'b' and 'b2' are NA: no row has a value ",
           "of any of them, a cell and a location."))
})


test_that("the warnings of one call come in the order of the design's table", {
  # 600 rows in 100 cells, and one column for each thing that can happen to a
  # column: x, then y, then the response.  Two rows share a location and three
  # rows have none.
  d <- local({
    .ps_seed(2020)
    x <- stats::runif(600, 0, 1e5); y <- stats::runif(600, 0, 1e5)
    x[2L] <- x[1L]; y[2L] <- y[1L]
    data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(600, 0, 0.5))
  })
  d$cell_id <- as.integer(10 * floor(d$x / 1e4) + floor(d$y / 1e4) + 1)
  cells   <- sort(unique(d$cell_id))
  count   <- table(d$cell_id)
  nowhere <- c(50L, 51L, 52L)
  # Values only on the rows without a location.
  d$none <- NA_real_
  d$none[nowhere] <- c(1, 2, 3)
  # One value in each of 30 cells at most, none of them on a row without a
  # location.
  d$single <- ifelse(!duplicated(d$cell_id) & d$cell_id %in% cells[1:30], d$z, NA_real_)
  d$single[nowhere] <- NA
  # Values in one cell only.
  d$one <- ifelse(d$cell_id == cells[40L], d$z, NA_real_)
  # A constant.
  d$k <- 0.1
  # Five rows in each of five cells, the same five values in every one: the
  # cell means are equal, so nothing tells the cells apart.
  five <- as.integer(names(count)[count >= 5L & names(count) != cells[40L]])[c(1L, 5L, 9L, 13L, 17L)]
  d$same <- NA_real_
  for (cl in five) d$same[which(d$cell_id == cl)[1:5]] <- c(-2, -1, 0, 1, 2)
  asg <- .ps_points(d$x, d$y, d[setdiff(names(d), c("x", "y"))])
  g <- sf::st_geometry(asg)
  for (i in nowhere) g[[i]] <- sf::st_point()
  sf::st_geometry(asg) <- g
  preds <- c("none", "single", "one", "k", "same")
  said <- function(r) vapply(r$warnings, function(w)
    paste(c(class(w)[1L], w$reason), collapse = ": "), character(1))

  r <- .ps_catch(summarize_by_cell(asg, "z", preds, id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.9))
  status <- vapply(attr(r$value, "pooled")$columns, function(c1) c1$status, "")
  expect_identical(
    status,
    c(resp_z = "ok", pred_none = "no located value", pred_single = "no cell holds two values",
      pred_one = "one cell", pred_k = "no spread inside any cell",
      pred_same = "no signal between cells"))
  expect_identical(
    said(r),
    c("spatialkit_pool_unavailable: no located value",
      "spatialkit_pool_unavailable: no cell holds two values",
      "spatialkit_pool_degenerate: one cell",
      "spatialkit_pool_degenerate: no spread inside any cell",
      "spatialkit_pool_degenerate: no signal between cells",
      "spatialkit_pool_few_cells",
      "spatialkit_pool_few_df",
      "spatialkit_no_structure",
      "spatialkit_pool_rows_skipped",
      "spatialkit_repeated_locations",
      "spatialkit_few_residuals"))
  # Each names the columns it is true of, in the order of the frame.  The
  # rows without a location hold a value of the response, of `none` and of the
  # constant.
  expect_identical(
    lapply(r$warnings, function(w) w$columns),
    list("none", "single", "one", "k", "same", "same", "same", "same", c("z", "none", "k"), "z",
         c("one", "same")))

  # Above the limit of an exact search: that reason stands after the two that
  # need no search and before the columns whose values were computed.
  .ps_without_fnn()
  withr::local_options(spatialkit.exact_pairs_max = 10)
  r <- .ps_catch(summarize_by_cell(asg, "z", preds, id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.9))
  expect_identical(
    said(r),
    c("spatialkit_pool_unavailable: no located value",
      "spatialkit_pool_unavailable: no cell holds two values",
      "spatialkit_pool_unavailable: work limit",
      "spatialkit_pool_degenerate: one cell",
      "spatialkit_pool_rows_skipped",
      "spatialkit_few_residuals"))
  expect_identical(r$warnings[[3L]]$columns, c("z", "k", "same"))
})


# ---- The work limit is a status -------------------------------------------------

test_that("above the limit of an exact search the pooled columns are NA and nothing stops", {
  d <- local({
    .ps_seed(2005)
    x <- stats::runif(400, 0, 1e5); y <- stats::runif(400, 0, 1e5)
    data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(400, 0, 0.5),
               p = cos(x / 1.2e4) + stats::rnorm(400, 0, 0.5))
  })
  A <- .ps_assign(.ps_points(d$x, d$y, d[c("z", "p")]), 40L)
  expect_identical(length(unique(A$asg$cell_id)), 40L)
  # The predictor has no value in one cell, so it has 39 cells where the
  # response has 40: the warning speaks of the larger search.
  A$asg$p[A$asg$cell_id == A$asg$cell_id[1L]] <- NA
  mean_call <- summarize_by_cell(A$asg, "z", "p", id_col = "cell_id", cells_sf = A$cells,
                                 conf_level = 0.9, estimate = "mean")
  .ps_without_fnn()
  withr::local_options(spatialkit.exact_pairs_max = 100)
  seen <- list()
  r <- withCallingHandlers(
    .ps_catch(summarize_by_cell(A$asg, "z", "p", id_col = "cell_id", cells_sf = A$cells,
                                conf_level = 0.9, estimate = "pooled", obs_level = 0.9)),
    condition = function(c) seen[[length(seen) + 1L]] <<- class(c))
  # The guard's error is caught inside the function: no condition of its class
  # reaches the caller, as an error or otherwise.
  expect_false(any(vapply(seen, function(k) "spatialkit_work_limit" %in% k, NA)))
  expect_identical(.ps_classes(r), "spatialkit_pool_unavailable")
  w <- r$warnings[[1L]]
  expect_identical(
    .ps_fields(w),
    list(reason = "work limit", columns = c("z", "p"), work = 1600, limit = 100,
         cells = c(z = 40L, p = 39L)))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): the pooled columns of 'z' and 'p' are NA: finding each cell's ",
           "neighbours among 40 cells without the FNN package would compare 1.6e+03 pairs, ",
           "above the limit of 1e+02 (options(spatialkit.exact_pairs_max = )). Install FNN: ",
           "the result is the same. The plain columns are unaffected."))
  out <- r$value
  new <- .ps_pooled_names(out)
  expect_length(new, 10L)
  for (col in new) expect_true(all(is.na(out[[col]])), info = col)
  rec <- attr(out, "pooled")$columns
  expect_identical(rec$resp_z$status, "work limit")
  expect_identical(rec$pred_p$status, "work limit")
  expect_identical(c(rec$resp_z$n, rec$resp_z$cells, rec$resp_z$df), c(400L, 40L, 360L))
  expect_identical(rec$pred_p$cells, 39L)
  # The predictor first: the same pairs and the same count in the text.
  flipped <- .ps_catch(summarize_by_cell(A$asg, "p", "z", id_col = "cell_id",
                                         estimate = "pooled"))
  expect_identical(.ps_fields(flipped$warnings[[1L]])[c("work", "cells")],
                   list(work = 1600, cells = c(p = 39L, z = 40L)))
  expect_match(conditionMessage(flipped$warnings[[1L]]),
               "neighbours among 40 cells without the FNN package would compare 1.6e+03 pairs",
               fixed = TRUE)
  expect_identical(rec$resp_z$offsets, c(NA_real_, NA_real_))
  # Every column the call shares with the call of plain means is that call's.
  expect_identical(setdiff(names(out), new), names(mean_call))
  for (col in names(mean_call)) expect_identical(out[[col]], mean_call[[col]], info = col)
  # At the limit the search runs, and the values are those with FNN.
  withr::local_options(spatialkit.exact_pairs_max = 1600)
  at_limit <- .ps_catch(summarize_by_cell(A$asg, "z", "p", id_col = "cell_id",
                                          estimate = "pooled", obs_level = 0.9))
  expect_length(at_limit$warnings, 0L)
  expect_identical(attr(at_limit$value, "pooled")$columns$resp_z$status, "ok")
})

test_that("without FNN the pooled columns are the ones with it", {
  skip_if_not_installed("FNN")
  for (layer in list(.ps_layer()$asg, .ps_lattice())) {
    v <- if ("z" %in% names(layer)) "z" else "v"
    with <- summarize_by_cell(layer, v, id_col = "cell_id", estimate = "pooled", obs_level = 0.9)
    without <- local({
      .ps_without_fnn()
      testthat::local_mocked_bindings(
        get.knn = function(...) stop("FNN reached"), get.knnx = function(...) stop("FNN reached"),
        .package = "FNN")
      summarize_by_cell(layer, v, id_col = "cell_id", estimate = "pooled", obs_level = 0.9)
    })
    expect_identical(without, with)
  }
})


# ---- Shared neighbour searches -----------------------------------------------------

test_that("columns with the same scored rows share one neighbour search", {
  L <- .ps_layer()
  real  <- .region_neighbours
  force(real)
  calls <- 0L
  testthat::local_mocked_bindings(.region_neighbours = function(...) {
    calls <<- calls + 1L
    real(...)
  })
  full <- summarize_by_cell(L$asg, "z", c("p1", "p2"), id_col = "cell_id", estimate = "pooled")
  expect_identical(calls, 1L)
  # A missing value in one predictor: that column has other scored rows, and a
  # search of its own.
  asg <- L$asg
  asg$p1[10L] <- NA
  calls <- 0L
  part <- summarize_by_cell(asg, "z", c("p1", "p2"), id_col = "cell_id", estimate = "pooled")
  expect_identical(calls, 2L)
  # Whichever column comes first, and a third mask is a third search.
  calls <- 0L
  summarize_by_cell(asg, "p1", c("z", "p2"), id_col = "cell_id", estimate = "pooled")
  expect_identical(calls, 2L)
  asg$p2[20L] <- NA
  calls <- 0L
  summarize_by_cell(asg, "z", c("p1", "p2"), id_col = "cell_id", estimate = "pooled")
  expect_identical(calls, 3L)
  # A row without a location is scored for no column: the columns still have
  # the same scored rows, and share one search.
  nowhere <- L$asg
  g <- sf::st_geometry(nowhere)
  g[[33L]] <- sf::st_point()
  sf::st_geometry(nowhere) <- g
  calls <- 0L
  suppressWarnings(summarize_by_cell(nowhere, "z", c("p1", "p2"), id_col = "cell_id",
                                     estimate = "pooled"))
  expect_identical(calls, 1L)
  # Sharing changes no number: each column is what it is when summarised alone.
  for (v in c("p1", "p2")) {
    own <- summarize_by_cell(L$asg, predictor_vars = v, id_col = "cell_id", estimate = "pooled")
    for (col in .ps_pooled_names(own)) expect_identical(full[[col]], own[[col]], info = col)
  }
  expect_identical(part$..pooled_resp_z, full$..pooled_resp_z)
  expect_identical(part$..pooled_pred_p2, full$..pooled_pred_p2)
  expect_false(identical(part$..pooled_pred_p1, full$..pooled_pred_p1))
})


# ---- The intervals table and the print ----------------------------------------------

test_that("the intervals table has a row per family in the frame, and is in the record", {
  G <- .ps_grid_layer()
  call <- function(obs, conf, ...)
    summarize_by_cell(G, "income", id_col = "cell_id", estimate = "pooled",
                      obs_level = obs, conf_level = conf, ...)
  cases <- list(neither = call(NULL, NULL), obs = call(0.9, NULL), conf = call(NULL, 0.95),
                both = call(0.9, 0.95))
  cols <- c("columns", "covers", "centre", "level", "method")
  for (nm in names(cases)) {
    tab <- attr(cases[[nm]], "intervals")
    expect_identical(tab, attr(cases[[nm]], "pooled")$intervals, info = nm)
    expect_s3_class(tab, "data.frame")
    expect_identical(names(tab), cols, info = nm)
    expect_identical(vapply(tab, typeof, ""),
                     c(columns = "character", covers = "character", centre = "character",
                       level = "double", method = "character"), info = nm)
  }
  expect_identical(vapply(cases, function(x) nrow(attr(x, "intervals")), 1L),
                   c(neither = 0L, obs = 1L, conf = 1L, both = 2L))
  obs_row <- data.frame(
    columns = "..obs_lo_, ..obs_hi_",
    covers  = paste0("one further observation in the cell, located like the data; on average ",
                     "over observations, not cell by cell"),
    centre  = "pooled value", level = 0.9,
    method  = "quantiles of leave-one-out residuals", stringsAsFactors = FALSE)
  ci_row <- data.frame(
    columns = "..ci_lo_, ..ci_hi_",
    covers  = paste0("the cell's own mean, when its observations are spread through the cell; ",
                     "not a range for one observation"),
    centre  = "plain mean", level = 0.95, method = "t, deff = 1", stringsAsFactors = FALSE)
  expect_identical(attr(cases$obs, "intervals"), obs_row)
  expect_identical(attr(cases$conf, "intervals"), ci_row)
  expect_identical(attr(cases$both, "intervals"), rbind(obs_row, ci_row))
  # `conf_level` is in the record, as `obs_level` is.
  expect_null(attr(cases$obs, "pooled")$conf_level)
  expect_identical(attr(cases$conf, "pooled")$conf_level, 0.95)
  expect_identical(attr(cases$both, "pooled")[c("obs_level", "conf_level")],
                   list(obs_level = 0.9, conf_level = 0.95))
  # The columns the table names are in the frame exactly when it names them.
  expect_identical(any(grepl("^\\.\\.obs_", names(cases$conf))), FALSE)
  expect_identical(any(grepl("^\\.\\.ci_", names(cases$obs))), FALSE)
  expect_true(all(c("..obs_lo_resp_income", "..ci_lo_resp_income") %in% names(cases$both)))
  # The helper alone.
  expect_identical(.intervals_table(0.9, 0.95, "1"), rbind(obs_row, ci_row))
  expect_identical(.intervals_table(NULL, NULL, "1"), rbind(obs_row, ci_row)[0L, ])
})

test_that("what the t interval covers follows the design effect", {
  L <- .ps_layer()
  own <- paste0("the cell's own mean, when its observations are spread through the cell; ",
                "not a range for one observation")
  pop <- "the cell mean as an estimate of the population mean; not a range for one observation"
  row <- function(...) {
    out <- suppressWarnings(summarize_by_cell(L$asg, "z", id_col = "cell_id", conf_level = 0.95,
                                              estimate = "pooled", ...))
    tab <- attr(out, "intervals")
    expect_identical(nrow(tab), 1L)
    expect_identical(tab$columns, "..ci_lo_, ..ci_hi_")
    expect_identical(tab$centre, "plain mean")
    unlist(tab[c("covers", "method")])
  }
  expect_identical(row(deff = 1), c(covers = own, method = "t, deff = 1"))
  expect_identical(row(), c(covers = own, method = "t, deff = 1"))
  expect_identical(row(deff = "kish"), c(covers = pop, method = "t, deff = \"kish\""))
  expect_identical(row(deff = 2), c(covers = pop, method = "t, deff = 2"))
  expect_identical(row(deff = 1.5), c(covers = pop, method = "t, deff = 1.5"))
  # A design effect the function refuses is replaced by 1, and the row says 1.
  expect_identical(row(deff = 0.5), c(covers = own, method = "t, deff = 1"))
  # A number above 1 by less than seven digits show is applied, and the row
  # does not read 1.
  near <- suppressWarnings(summarize_by_cell(L$asg, "z", id_col = "cell_id", conf_level = 0.95,
                                             estimate = "pooled", deff = 1.0000001))
  expect_identical(unique(near$deff_method), "fixed")
  expect_identical(unlist(attr(near, "intervals")[c("covers", "method")]),
                   c(covers = pop, method = "t, deff = 1.0000001"))
  # Kish's correction with an ICC of 0 corrects nothing: the interval is the
  # one of deff = 1, `deff_method` says "none", and the row describes the
  # numbers the frame holds.  `u` has the same mean in every cell.
  flat <- L$asg
  flat$u <- flat$p2 - stats::ave(flat$p2, flat$cell_id)
  ci <- function(resp, pred = NULL, deff)
    suppressWarnings(summarize_by_cell(flat, resp, pred, id_col = "cell_id", conf_level = 0.95,
                                       estimate = "pooled", deff = deff))
  k0 <- ci("u", deff = "kish")
  expect_identical(attr(k0, "icc")$resp, 0)
  expect_identical(unique(k0$deff_method), "none")
  expect_identical(k0$..ci_lo_resp_u, ci("u", deff = 1)$..ci_lo_resp_u)
  expect_identical(unlist(attr(k0, "intervals")[c("covers", "method")]),
                   c(covers = own, method = "t, deff = 1"))
  expect_match(utils::capture.output(print(attr(k0, "pooled")))[4L],
               "95% for THE CELL'S MEAN (t, deff = 1), centred on the plain mean;", fixed = TRUE)
  # Beside columns Kish does correct, the row is Kish's and says which column
  # it left as it was.
  both <- ci("z", c("u", "p1"), deff = "kish")
  expect_identical(unique(both$deff_method), "kish")
  expect_gt(attr(both, "icc")$resp, 0.1)
  expect_identical(attr(both, "icc")$pred[["u"]], 0)
  expect_identical(unlist(attr(both, "intervals")[c("covers", "method")]),
                   c(covers = pop, method = "t, deff = \"kish\"; 1 for 'u', whose ICC is 0"))
  # The same when it is the response Kish leaves as it was, and a predictor
  # it corrects.
  resp0 <- ci("u", c("z", "p1"), deff = "kish")
  expect_identical(attr(resp0, "icc")$resp, 0)
  expect_identical(unique(resp0$deff_method), "kish")
  expect_identical(attr(resp0, "intervals")$method,
                   "t, deff = \"kish\"; 1 for 'u', whose ICC is 0")
  txt <- utils::capture.output(print(attr(resp0, "pooled")))
  expect_identical(
    grep("..ci_lo_", txt, fixed = TRUE, value = TRUE),
    paste0("  ..ci_lo_, ..ci_hi_ : 95% for THE POPULATION MEAN (t, deff = \"kish\"; 1 for 'u', ",
           "whose ICC is 0), centred on the plain mean;"))
  # The label alone, and the table told what was applied.
  expect_identical(.deff_label(NULL), "1")
  expect_identical(.deff_label(NULL, kish_zero = "u"), "1")
  expect_identical(.deff_label(list(method = "fixed", deff = 2.5)), "2.5")
  expect_identical(.deff_label(list(method = "variogram")), "\"variogram\"")
  expect_identical(.deff_label(list(method = "kish")), "\"kish\"")
  expect_identical(.deff_label(list(method = "kish"), c("a", "b")),
                   "\"kish\"; 1 for 'a' and 'b', whose ICC is 0")
  tiny <- .deff_label(list(method = "fixed", deff = 1 + 4 * .Machine$double.eps))
  expect_false(identical(tiny, "1"))
  expect_gt(as.numeric(tiny), 1)
  expect_identical(.intervals_table(NULL, 0.95, "1")$covers, own)
  expect_identical(.intervals_table(NULL, 0.95, tiny)$covers, pop)
  expect_identical(.intervals_table(NULL, 0.95, "\"kish\"; 1 for 'u', whose ICC is 0")$covers, pop)
  skip_if_not_installed("gstat")
  # A supplied variogram, on the layer where one is identified.
  L <- .ps_field_layer()
  sac <- suppressWarnings(estimate_sac_range(L$asg, "z"))
  expect_true(is.finite(as.numeric(sac)))
  expect_identical(row(deff = "variogram", sac = sac),
                   c(covers = pop, method = "t, deff = \"variogram\""))
  # A "variogram" request that falls back to Kish carries Kish's correction,
  # and the row names the design effect the standard errors carry.
  testthat::local_mocked_bindings(.vgm_correlation_fn = function(...) NULL,
                                  .vgm_is_pure_nugget = function(...) FALSE)
  fell <- .ps_catch(summarize_by_cell(L$asg, "z", id_col = "cell_id", conf_level = 0.95,
                                      estimate = "pooled", deff = "variogram", sac = sac))
  expect_true("spatialkit_deff_fallback" %in% .ps_classes(fell))
  expect_identical(unique(fell$value$deff_method), "kish (variogram fallback)")
  expect_identical(unlist(attr(fell$value, "intervals")[c("covers", "method")]),
                   c(covers = pop, method = "t, deff = \"kish\""))
})

test_that("the record prints what was fitted, and each interval under its own words", {
  G <- .ps_grid_layer()
  shown <- function(obs, conf, ...)
    utils::capture.output(print(attr(
      summarize_by_cell(G, "income", id_col = "cell_id", estimate = "pooled",
                        obs_level = obs, conf_level = conf, ...), "pooled")))
  block <- c(
    "Cell values: income, 900 rows in 100 cells",
    "  value       : ..pooled_resp_income, pooled toward the 8 nearest cells (centroids of the rows)",
    "  variances   : inside cells 3.80e+07 (from 800 values beyond one per cell);",
    "                between a cell and its target 3.84e+06",
    "  own weight  : median 0.48 (..pooled_w_resp_income)",
    "  alone       : 0.0% of the rows are the only one in their cell",
    "  loo         : 4.02e+07 (leave-one-out mean squared error of the pooled value);",
    "                one mean for all rows: 1.01e+08")
  interval <- c(
    "  interval    : 90% for ONE OBSERVATION (..obs_lo_resp_income, ..obs_hi_resp_income):",
    "                value - 10500 to value + 10300, quantiles of 900 leave-one-out residuals")
  also <- c(
    "Also in the frame",
    "  ..ci_lo_, ..ci_hi_ : 95% for THE CELL'S MEAN (t, deff = 1), centred on the plain mean;",
    "                       not a range for one observation")
  closing <- "Away from the observations no interval here is calibrated."
  expect_identical(shown(NULL, NULL), c(block, closing))
  expect_identical(shown(0.9, NULL), c(block, interval, closing))
  expect_identical(shown(NULL, 0.95), c(block, also, closing))
  expect_identical(shown(0.9, 0.95), c(block, interval, also, closing))
  # Under a design effect the other interval is for the population mean.
  expect_identical(
    shown(NULL, 0.9, deff = "kish")[10:11],
    c("  ..ci_lo_, ..ci_hi_ : 90% for THE POPULATION MEAN (t, deff = \"kish\"), centred on the plain mean;",
      "                       not a range for one observation"))
  # print() returns its argument, invisibly, and reads nothing but the record.
  rec <- attr(summarize_by_cell(G, "income", id_col = "cell_id", estimate = "pooled"), "pooled")
  utils::capture.output(seen <- withVisible(print(rec)))
  expect_false(seen$visible)
  expect_identical(seen$value, rec)
  expect_identical(utils::capture.output(print(structure(unclass(rec), class = "spatialkit_pooled"))),
                   c(block, closing))
})

test_that("the print says when cells do no better, when all the others are the target, and a status", {
  # A structureless predictor beside the response: its loo line says so.
  G <- .ps_grid_layer()
  G$noise <- local({ .ps_seed(2006); stats::rnorm(nrow(G)) })
  r <- .ps_catch(summarize_by_cell(G, "income", "noise", id_col = "cell_id",
                                   estimate = "pooled", obs_level = 0.9))
  expect_identical(.ps_classes(r), "spatialkit_no_structure")
  txt <- utils::capture.output(print(attr(r$value, "pooled")))
  expect_length(txt, 21L)
  expect_identical(txt[11L], "Cell values: noise, 900 rows in 100 cells")
  expect_match(txt[18L], "^                one mean for all rows: [0-9.e+-]+: these cells do no better than one mean$")
  expect_false(grepl("do no better", txt[8L], fixed = TRUE))
  expect_match(txt[19L], "(..obs_lo_pred_noise, ..obs_hi_pred_noise):", fixed = TRUE)
  expect_identical(txt[21L], "Away from the observations no interval here is calibrated.")
  # Five cells: the target is all the others.
  five <- G
  five$cell_id <- as.integer((five$cell_id - 1L) %/% 20L + 1L)
  r <- .ps_catch(summarize_by_cell(five, "income", id_col = "cell_id", estimate = "pooled"))
  expect_identical(.ps_classes(r), "spatialkit_pool_few_cells")
  txt <- utils::capture.output(print(attr(r$value, "pooled")))
  expect_identical(txt[1:2],
                   c("Cell values: income, 900 rows in 5 cells",
                     "  value       : ..pooled_resp_income, pooled toward all the other cells (4 of them)"))
  # Nine cells have eight others, all of them in each target; ten do not.
  xy <- sf::st_coordinates(G)
  value_line <- function(cell) {
    g <- G
    g$cell_id <- as.integer(cell)
    utils::capture.output(print(attr(suppressWarnings(
      summarize_by_cell(g, "income", id_col = "cell_id", estimate = "pooled")), "pooled")))[1:2]
  }
  third <- 1e5 / 3
  expect_identical(value_line(3 * floor(xy[, 1L] / third) + floor(xy[, 2L] / third) + 1),
                   c("Cell values: income, 900 rows in 9 cells",
                     "  value       : ..pooled_resp_income, pooled toward all the other cells (8 of them)"))
  expect_identical(value_line(2 * floor(xy[, 1L] / 2e4) + floor(xy[, 2L] / 5e4) + 1),
                   c("Cell values: income, 900 rows in 10 cells",
                     "  value       : ..pooled_resp_income, pooled toward the 8 nearest cells (centroids of the rows)"))
  # Rows that are the only one in their cell are counted in the `alone` line.
  lone <- .ps_grid_layer()
  lone$cell_id[c(10L, 200L, 400L, 600L, 800L)] <- 501:505
  txt <- utils::capture.output(print(attr(suppressWarnings(
    summarize_by_cell(lone, "income", id_col = "cell_id", estimate = "pooled")), "pooled")))
  expect_identical(txt[1L], "Cell values: income, 900 rows in 105 cells")
  expect_identical(txt[6L], "  alone       : 0.6% of the rows are the only one in their cell")
  # A field in 120 cells of which 117 hold one row.  The cells predict a
  # left-out row several times better than one mean, and the F test, with
  # three values beyond one per cell behind it, does not tell them from one
  # mean: the loo line says that, and does not say the cells do no better.
  d <- local({
    .ps_seed(3262)
    x <- stats::runif(123, 0, 1e5); y <- stats::runif(123, 0, 1e5)
    data.frame(x = x, y = y, v = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(123, 0, 0.3))
  })
  few <- .ps_points(d$x, d$y, data.frame(v = d$v, cell_id = c(1:120, 1:3)))
  r <- .ps_catch(summarize_by_cell(few, "v", id_col = "cell_id", estimate = "pooled"))
  expect_identical(.ps_classes(r), c("spatialkit_pool_few_df", "spatialkit_no_structure"))
  rec <- attr(r$value, "pooled")$columns$resp_v
  expect_identical(rec$status, "ok")
  expect_false(rec$structured)
  expect_gt(rec$loo_one_cell, 3 * rec$loo)
  txt <- utils::capture.output(print(attr(r$value, "pooled")))
  expect_identical(
    txt[8L],
    sprintf(paste0("                one mean for all rows: %.2e: the F test does not tell these ",
                   "cells from one mean (p = %s)"),
            rec$loo_one_cell, format(signif(rec$structure_p, 2))))
  expect_identical(txt[6L], "  alone       : 95.1% of the rows are the only one in their cell")
  expect_match(conditionMessage(r$warnings[[2L]]),
               "predict a left-out observation better than one mean for all rows, but the F test",
               fixed = TRUE)
  # A status other than "ok": the first line and the status, and the interval
  # the block has not shown is listed below it.
  G$k <- 0.1
  r <- .ps_catch(summarize_by_cell(G, "k", id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.9, conf_level = 0.95))
  expect_identical(
    utils::capture.output(print(attr(r$value, "pooled"))),
    c("Cell values: k, 900 rows in 100 cells",
      "  status      : no spread inside any cell",
      "Also in the frame",
      "  ..obs_lo_, ..obs_hi_ : 90% for ONE OBSERVATION (quantiles of leave-one-out residuals), centred on the pooled value;",
      "                         on average over observations, not cell by cell",
      "  ..ci_lo_, ..ci_hi_ : 95% for THE CELL'S MEAN (t, deff = 1), centred on the plain mean;",
      "                       not a range for one observation",
      "Away from the observations no interval here is calibrated."))
  # One row in one cell.
  r <- .ps_catch(summarize_by_cell(G[1L, ], "income", id_col = "cell_id", estimate = "pooled"))
  expect_identical(utils::capture.output(print(attr(r$value, "pooled"))),
                   c("Cell values: income, 1 row in 1 cell", "  status      : one cell",
                     "Away from the observations no interval here is calibrated."))
})


# ---- A response without structure ----------------------------------------------------

test_that("white noise is reported as structureless, and a field with noise is not", {
  skip_on_cran()
  # 40 layers of each kind: 3,000 locations uniform on a 1e5 square, x then y,
  # then the response; the points in EPSG:32617; the cells as a user builds
  # them at 667 seeds of the density kind from 3 starts.
  one <- function(seed, field) {
    pts <- local({
      .ps_seed(seed)
      x <- stats::runif(3000, 0, 1e5); y <- stats::runif(3000, 0, 1e5)
      z <- if (field) sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(3000, 0, 1.5)
           else stats::rnorm(3000)
      .ps_points(x, y, data.frame(z = z))
    })
    seeds <- suppressWarnings(
      voronoi_seeds_kmeans(pts, 667, set_seed = 123, nstart = 3, partition = "density"))
    tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id")
    r <- .ps_catch(summarize_by_cell(asg, "z", id_col = "cell_id", estimate = "pooled"))
    list(r = r, seeds = seeds, pts = pts, rows = nrow(asg),
         tau2 = attr(r$value, "pooled")$columns$resp_z$tau2)
  }
  raised <- function(x) length(.ps_of(x$r, "spatialkit_no_structure")) == 1L
  noise <- lapply(1:40, function(r) one(6000 + r, FALSE))
  field <- lapply(1:40, function(r) one(7000 + r, TRUE))
  expect_true(all(vapply(c(noise, field), function(x) x$rows, 1L) == 3000L))
  # One run of these seeds: raised on 40 of the noise layers, with tau2 exactly
  # 0 on 2 of them; raised on none of the field layers, with tau2 0 on none.
  expect_gte(sum(vapply(noise, raised, NA)), 38L)
  expect_lte(sum(vapply(noise, function(x) x$tau2 == 0, NA)), 8L)
  expect_lte(sum(vapply(field, raised, NA)), 1L)
  expect_lte(sum(vapply(field, function(x) x$tau2 == 0, NA)), 1L)
  # Never more than one for a column, and nothing else but the warning for
  # every own weight being 0, which white noise brings by chance.
  expect_true(all(unlist(lapply(c(noise, field), function(x) .ps_classes(x$r))) %in%
                    c("spatialkit_no_structure", "spatialkit_pool_degenerate")))
  expect_true(all(vapply(noise, function(x) sum(.ps_classes(x$r) == "spatialkit_no_structure"),
                         1L) <= 1L))
  # The condition's fields, on the first noise layer.
  w   <- .ps_of(noise[[1L]]$r, "spatialkit_no_structure")[[1L]]
  rec <- attr(noise[[1L]]$r$value, "pooled")$columns$resp_z
  expect_identical(
    .ps_fields(w),
    list(criterion = "loo", n = c(z = rec$cells), one_cell = c(z = rec$loo_one_cell),
         columns = "z", estimate = "pooled", loo = c(z = rec$loo),
         structure_p = c(z = rec$structure_p), beats = c(z = FALSE)))
  expect_false(rec$structured)
  # One mean does as well: its error is within 2 percent of the cells', or lower.
  expect_lte(rec$loo_one_cell, rec$loo * 1.02)
  # The seeds are the centres of one fit under set.seed(123), which is what a
  # profile keeps at 667 cells of this kind, so later stages meet these cells.
  for (x in list(noise[[1L]], field[[1L]])) {
    xy <- unname(sf::st_coordinates(x$pts)[, 1:2])
    kb <- local({
      .ps_seed(123)
      .partition_centres(xy, 667L, partition = "density", init = "kmeans++", nstart = 3L,
                         iter.max = 100L)
    })
    expect_identical(unname(sf::st_coordinates(x$seeds)[, 1:2]), unname(kb$centers))
  }
})


# ---- Any IDs, any polygons, any row order ---------------------------------------------

test_that("the numbering and the type of the IDs change no value", {
  L <- .ps_layer()
  asg <- L$asg
  u <- sort(unique(asg$cell_id)); K <- length(u)
  numbered <- asg
  numbered$cell_id <- match(asg$cell_id, u)                      # 1..K
  base <- summarize_by_cell(numbered, "z", "p1", id_col = "cell_id", estimate = "pooled",
                            obs_level = 0.9)
  new <- .ps_pooled_names(base)
  same_values <- function(out, key, tag) {
    at <- match(key, out$cell_id)
    expect_false(anyNA(at), info = tag)
    for (col in new) expect_equal(out[[col]][at], base[[col]], tolerance = 1e-12,
                                  info = paste(tag, col))
    a <- attr(out, "pooled")$columns$resp_z; b <- attr(base, "pooled")$columns$resp_z
    expect_equal(a, b, tolerance = 1e-12, info = tag)
  }
  # Character IDs in an order that is not the numbers'.
  lab <- local({ .ps_seed(61); sample(c(sprintf("z%03d", 1:75), sprintf("A-%d", 1:75))) })
  expect_lte(K, length(lab))
  lab <- lab[seq_len(K)]
  chr <- numbered
  chr$cell_id <- lab[numbered$cell_id]
  out <- summarize_by_cell(chr, "z", "p1", id_col = "cell_id", estimate = "pooled",
                           obs_level = 0.9)
  expect_false(identical(order(lab, method = "radix"), seq_len(K)))
  expect_identical(out$cell_id, sort(lab, method = "radix"))
  same_values(out, lab, "character")
  # Factor, double and negative IDs.
  fac <- numbered
  fac$cell_id <- factor(lab[numbered$cell_id], levels = rev(lab))
  same_values(summarize_by_cell(fac, "z", "p1", id_col = "cell_id", estimate = "pooled",
                                obs_level = 0.9), factor(lab, levels = rev(lab)), "factor")
  dbl <- numbered
  dbl$cell_id <- (K:1)[numbered$cell_id] * 1.5 - 100
  same_values(summarize_by_cell(dbl, "z", "p1", id_col = "cell_id", estimate = "pooled",
                                obs_level = 0.9), (K:1) * 1.5 - 100, "double")
})

test_that("polygons the user brings, overlapping ones and a grid, are read by their IDs alone", {
  L <- .ps_layer()
  pts <- L$asg[c("z", "p1")]
  bnd <- sf::st_as_sfc(sf::st_bbox(pts))
  grid <- create_grid_polygons(bnd, target_cells = 100, type = "square", quiet = TRUE)
  id <- intersect(c("poly_id", "polygon_id", "id", "cell_id", "grid_id"), names(grid))[1L]
  # Overlapping squares: each cell of the grid grown by a third of its width.
  wide <- grid
  sf::st_geometry(wide) <- sf::st_buffer(sf::st_geometry(grid),
                                         dist = sqrt(as.numeric(sf::st_area(grid))[1L]) / 3,
                                         endCapStyle = "SQUARE", joinStyle = "MITRE")
  wide$GEOID <- sprintf("U%03d", rev(seq_len(nrow(wide))))
  expect_gt(sum(lengths(sf::st_overlaps(wide)) > 0L), 50L)
  for (case in list(list("grid", grid, id), list("overlapping", wide, "GEOID"))) {
    cells <- case[[2L]]; col <- case[[3L]]
    asg <- suppressWarnings(assign_features_to_polygons(pts, cells, polygon_id_col = col))
    expect_identical(nrow(asg), nrow(pts), info = case[[1L]])
    r <- .ps_catch(summarize_by_cell(asg, "z", "p1", id_col = col, cells_sf = cells,
                                     estimate = "pooled", obs_level = 0.9))
    expect_length(r$warnings, 0L)
    out <- r$value
    u <- unique(asg[[col]])
    u <- u[order(u, method = "radix")]
    numbered <- asg
    numbered$region <- match(asg[[col]], u)
    base <- summarize_by_cell(numbered, "z", "p1", id_col = "region", estimate = "pooled",
                              obs_level = 0.9)
    at <- match(u, out[[col]])
    for (nm in .ps_pooled_names(base))
      expect_equal(out[[nm]][at], base[[nm]], tolerance = 1e-12, info = paste(case[[1L]], nm))
    expect_equal(attr(out, "pooled")$columns, attr(base, "pooled")$columns, tolerance = 1e-12)
    # A polygon no row was given is NA and is no cell's neighbour.
    expect_true(all(is.na(out$..pooled_resp_z[!(out[[col]] %in% u)])), info = case[[1L]])
    # No record of a count on cells a user brings.
    expect_null(attr(out, "count_rule"))
  }
})

test_that("the order of the rows changes no value, on scattered and on gridded locations", {
  R <- .ps_reference()
  layers <- list(`the (3000, 4.5) layout` = list(R$asg, "v"), `the lattice` = list(.ps_lattice(), "v"))
  for (nm in names(layers)) {
    layer <- layers[[nm]][[1L]]; v <- layers[[nm]][[2L]]
    base <- summarize_by_cell(layer, v, id_col = "cell_id", estimate = "pooled", obs_level = 0.9)
    for (s in 1:3) {
      p <- local({ .ps_seed(70 + s); sample.int(nrow(layer)) })
      out <- summarize_by_cell(layer[p, ], v, id_col = "cell_id", estimate = "pooled",
                               obs_level = 0.9)
      expect_identical(out$cell_id, base$cell_id, info = nm)
      for (col in .ps_pooled_names(base))
        expect_equal(out[[col]], base[[col]], tolerance = 1e-12, info = paste(nm, col, s))
      expect_equal(attr(out, "pooled")$columns, attr(base, "pooled")$columns,
                   tolerance = 1e-12, info = nm)
    }
  }
  # On the lattice the centroids of cells tie in distance, which is where a
  # result that followed the rows would move by far more than rounding.
  lat <- .ps_lattice()
  xy  <- unname(sf::st_coordinates(lat)[, 1:2])
  cen <- rowsum(xy, lat$cell_id) / tabulate(lat$cell_id)
  d   <- as.matrix(stats::dist(cen)); diag(d) <- Inf
  near <- t(apply(d, 1L, sort))[, 1:9]
  expect_gt(sum(abs(near[, 9L] - near[, 8L]) < 1e-9), 20L)
})

test_that("lon/lat rows give the values of the same rows projected by ensure_projected()", {
  L <- .ps_layer()
  ll   <- sf::st_transform(L$asg, 4326)
  proj <- ensure_projected(ll)
  expect_false(isTRUE(sf::st_is_longlat(proj)))
  a <- summarize_by_cell(ll, "z", "p1", id_col = "cell_id", estimate = "pooled", obs_level = 0.9)
  b <- summarize_by_cell(proj, "z", "p1", id_col = "cell_id", estimate = "pooled",
                         obs_level = 0.9)
  for (col in .ps_pooled_names(a))
    expect_equal(a[[col]], b[[col]], tolerance = 1e-8, info = col)
  # The distances were measured in the projected CRS, never in degrees.
  expect_identical(attr(a, "pooled")$crs, sf::st_crs(proj))
  expect_false(isTRUE(sf::st_is_longlat(attr(a, "pooled")$crs)))
  # Degrees read as distances would give other neighbours: the layer is wider
  # in degrees of longitude than in metres, relative to its height.
  deg <- sf::st_set_crs(sf::st_set_crs(ll, NA), 32617)
  c0  <- summarize_by_cell(deg, "z", id_col = "cell_id", estimate = "pooled")
  expect_false(isTRUE(all.equal(c0$..pooled_target_resp_z, a$..pooled_target_resp_z,
                                tolerance = 1e-8)))
  # With the cells joined the frame keeps the cells' CRS.
  joined <- summarize_by_cell(ll, "z", id_col = "cell_id", cells_sf = L$cells,
                              estimate = "pooled")
  expect_identical(sf::st_crs(joined), sf::st_crs(L$cells))
})

test_that("on lon/lat rows the projection is chosen from every row of the layer", {
  # Rows that hold no value of a column are as if absent for it, with one
  # stated exception: lon/lat coordinates are projected in the CRS
  # ensure_projected() picks for the whole layer, as the profile projects
  # them, and the rows without a value are rows of the layer.
  L  <- .ps_layer()
  ll <- sf::st_transform(L$asg, 4326)
  bb <- sf::st_bbox(ll)
  more <- function(lon, lat) {
    e <- sf::st_as_sf(data.frame(lon = lon, lat = lat, z = NA_real_, p1 = NA_real_,
                                 p2 = NA_real_, cell_id = 9000L + seq_along(lon)),
                      coords = c("lon", "lat"), crs = 4326)
    rbind(ll, e[names(ll)])
  }
  call <- function(x) summarize_by_cell(x, "z", "p1", id_col = "cell_id", estimate = "pooled",
                                        obs_level = 0.9)
  alone <- call(ll)
  new   <- .ps_pooled_names(alone)
  # Ten valueless rows inside the extent: the same CRS, and the same values.
  near <- call(more(seq(bb[["xmin"]], bb[["xmax"]], length.out = 12L)[2:11],
                    seq(bb[["ymin"]], bb[["ymax"]], length.out = 12L)[2:11]))
  expect_identical(attr(near, "pooled")$crs, attr(alone, "pooled")$crs)
  at <- match(alone$cell_id, near$cell_id)
  for (col in new) expect_identical(near[[col]][at], alone[[col]], info = col)
  expect_identical(attr(near, "pooled")$columns, attr(alone, "pooled")$columns)
  # Forty valueless rows across a continent: another CRS is picked for the
  # layer, the distances between centroids are measured in it, and the values
  # are those of the valued rows projected to that CRS beforehand.
  far <- call(more(seq(-124, -70, length.out = 40L), rep(c(26, 48), 20L)))
  crs_far <- attr(far, "pooled")$crs
  expect_false(isTRUE(sf::st_is_longlat(crs_far)))
  expect_false(identical(crs_far, attr(alone, "pooled")$crs))
  at <- match(alone$cell_id, far$cell_id)
  expect_false(isTRUE(all.equal(far$..pooled_target_resp_z[at], alone$..pooled_target_resp_z,
                                tolerance = 1e-6)))
  first <- call(sf::st_transform(ll, crs_far))
  for (col in new) expect_equal(far[[col]][at], first[[col]], tolerance = 1e-10, info = col)
  expect_equal(attr(far, "pooled")$columns, attr(first, "pooled")$columns, tolerance = 1e-10)
  # The cells of the valueless rows have no pooled value.
  expect_true(all(is.na(far$..pooled_resp_z[far$cell_id > 9000L])))
  # A layer projected beforehand is used in its own CRS, whatever rows it holds.
  proj <- sf::st_transform(more(seq(-124, -70, length.out = 40L), rep(c(26, 48), 20L)), 5070)
  a <- call(proj)
  b <- call(proj[!is.na(proj$z), ])
  expect_identical(attr(a, "pooled")$crs, sf::st_crs(5070))
  at <- match(b$cell_id, a$cell_id)
  for (col in new) expect_identical(a[[col]][at], b[[col]], info = col)
})


# ---- agg_funs --------------------------------------------------------------------------

test_that("the pooled value pools the mean whatever agg_funs computes", {
  L <- .ps_layer()
  base <- summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled")
  med  <- summarize_by_cell(L$asg, "z", id_col = "cell_id", agg_funs = median,
                            estimate = "pooled")
  expect_true("resp_median_z" %in% names(med))
  expect_false("resp_mean_z" %in% names(med))
  for (col in .ps_pooled_names(base)) expect_identical(med[[col]], base[[col]], info = col)
  expect_identical(attr(med, "pooled"), attr(base, "pooled"))
  # An entry named `pooled` keeps its own column beside the pooled value's.
  both <- summarize_by_cell(L$asg, "z", id_col = "cell_id",
                            agg_funs = list(pooled = function(x) max(x, na.rm = TRUE),
                                            mean = function(x) mean(x, na.rm = TRUE)),
                            estimate = "pooled")
  expect_true(all(c("resp_pooled_z", "..pooled_resp_z") %in% names(both)))
  expect_identical(both$..pooled_resp_z, base$..pooled_resp_z)
  expect_identical(both$resp_pooled_z,
                   as.numeric(tapply(L$asg$z, L$asg$cell_id, max))[match(both$cell_id, sort(unique(L$asg$cell_id)))])
  expect_false(isTRUE(all.equal(both$resp_pooled_z, both$..pooled_resp_z)))
})


# ---- obs_level ---------------------------------------------------------------------------

test_that("the interval's offsets are the stated quantiles of the leave-one-out residuals", {
  R <- .ps_reference()
  asg <- R$asg
  xy <- unname(sf::st_coordinates(asg)[, 1:2])
  u  <- sort(unique(asg$cell_id))
  e  <- .pool_fit(asg$v, match(asg$cell_id, u), length(u), xy, want = "pooled")$resid_pooled
  for (level in c(0.5, 0.9, 0.99)) {
    r <- .ps_catch(summarize_by_cell(asg, "v", id_col = "cell_id", estimate = "pooled",
                                     obs_level = level))
    expect_length(r$warnings, 0L)
    out <- r$value
    p   <- c((1 - level) / 2, (1 + level) / 2)
    off <- attr(out, "pooled")$columns$resp_v$offsets
    expect_equal(off, stats::quantile(e, p, names = FALSE, type = 7), tolerance = 1e-12,
                 info = format(level))
    expect_identical(off, stats::quantile(e, round(p, 12), names = FALSE, type = 7))
    expect_identical(out$..obs_lo_resp_v, out$..pooled_resp_v + off[1L])
    expect_identical(out$..obs_hi_resp_v, out$..pooled_resp_v + off[2L])
    # The residuals straddle zero, so the value lies inside its interval, and
    # the two offsets are the same in every cell.
    expect_lt(off[1L], 0); expect_gt(off[2L], 0)
    expect_true(all(out$..obs_lo_resp_v <= out$..pooled_resp_v &
                      out$..pooled_resp_v <= out$..obs_hi_resp_v))
    expect_identical(attr(out, "pooled")$obs_level, level)
    expect_identical(attr(out, "intervals")$level, level)
    # The share of residuals inside the offsets is the level.
    expect_equal(mean(e >= off[1L] & e <= off[2L]), level, tolerance = 2 / length(e))
  }
  # A wider level is a wider interval.
  width <- vapply(c(0.5, 0.9, 0.99), function(level)
    diff(attr(summarize_by_cell(asg, "v", id_col = "cell_id", estimate = "pooled",
                                obs_level = level), "pooled")$columns$resp_v$offsets), 1)
  expect_false(is.unsorted(width, strictly = TRUE))
})

test_that("obs_level is refused with plain means and when it is not a level", {
  L <- .ps_layer()
  expect_error(
    summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "mean", obs_level = 0.9),
    paste0("summarize_by_cell(): `obs_level` needs estimate = \"pooled\": the interval for an ",
           "observation is built from the leave-one-out residuals of the pooled value and is ",
           "centred on it. `conf_level` gives the plain mean's t interval, which is not a range ",
           "for an observation."), fixed = TRUE)
  # A confidence level beside it changes nothing.
  expect_error(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "mean",
                                 obs_level = 0.9, conf_level = 0.9),
               "`obs_level` needs estimate = \"pooled\"", fixed = TRUE)
  text <- paste0("summarize_by_cell(): `obs_level` must be a single number strictly between 0 ",
                 "and 1 (0.90 for a 90% interval), or NULL for no interval.")
  # 1 - 1e-13 and 1e-13: the interval reads a level through its tail
  # probability rounded to 12 decimals, to which these are 1 and 0.
  for (bad in list(0, 1, -0.1, 1.5, 90, c(0.8, 0.9), NA_real_, "0.9", TRUE, Inf, numeric(0),
                   1 - 1e-13, 1e-13))
    for (est in c("pooled", "mean")) {
      e <- tryCatch(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = est,
                                      obs_level = bad), error = function(e) e)
      expect_s3_class(e, "error")
      expect_identical(conditionMessage(e), text, info = paste(format(bad), collapse = ","))
      expect_null(conditionCall(e))
    }
  # Both checks run before any summary: a layer that cannot be summarised
  # (no ID column) still stops on the level.
  expect_error(summarize_by_cell(L$asg["z"], "z", estimate = "mean", obs_level = 0.9),
               "`obs_level` needs estimate", fixed = TRUE)
  expect_error(summarize_by_cell(L$asg["z"], "z", estimate = "pooled", obs_level = 2),
               "`obs_level` must be a single number", fixed = TRUE)
})

test_that("an interval on few residuals is returned with spatialkit_few_residuals", {
  L <- .ps_layer()
  # 96 rows of the layer in 12 cells: fewer than the 100 residuals at which
  # five are expected beyond each end of a 90 percent interval.
  d <- local({
    .ps_seed(2007)
    x <- stats::runif(96, 0, 1e5); y <- stats::runif(96, 0, 1e5)
    data.frame(x = x, y = y, z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(96, 0, 0.3),
               q = cos(x / 1.5e4) + stats::rnorm(96, 0, 0.3))
  })
  small <- .ps_points(d$x, d$y, d[c("z", "q")])
  small$cell_id <- as.integer(4 * floor(d$x / 2.5e4) + floor(d$y / 3.4e4) + 1)
  expect_identical(length(unique(small$cell_id)), 12L)
  r <- .ps_catch(summarize_by_cell(small, "z", "q", id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.9))
  expect_identical(.ps_classes(r), "spatialkit_few_residuals")
  w <- r$warnings[[1L]]
  expect_identical(.ps_fields(w),
                   list(columns = c("z", "q"), needed = 100L, level = 0.9,
                        n = c(z = 96L, q = 96L)))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): the 90% interval of 'z' and 'q' rests on fewer than 100 ",
           "leave-one-out residuals (96, 96), the number at which five are expected beyond ",
           "each end; its ends are unstable."))
  expect_false(anyNA(r$value$..obs_lo_resp_z))
  # Without the level nothing is said; at 0.5 twenty residuals are enough.
  expect_length(.ps_catch(summarize_by_cell(small, "z", "q", id_col = "cell_id",
                                            estimate = "pooled"))$warnings, 0L)
  expect_length(.ps_catch(summarize_by_cell(small, "z", "q", id_col = "cell_id",
                                            estimate = "pooled", obs_level = 0.5))$warnings, 0L)
  # At 0.99 a thousand are needed, and the 900 rows of the layer are too few;
  # at 0.9 they are enough.
  r <- .ps_catch(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.99))
  expect_identical(.ps_classes(r), "spatialkit_few_residuals")
  expect_identical(.ps_fields(r$warnings[[1L]]),
                   list(columns = "z", needed = 1000L, level = 0.99, n = c(z = 900L)))
  expect_match(conditionMessage(r$warnings[[1L]]), "the 99% interval of 'z' rests on fewer than 1000",
               fixed = TRUE)
  # Exactly at the number nothing is said: 100 rows at 0.9.
  d100 <- rbind(sf::st_drop_geometry(small)[c("z", "cell_id")],
                data.frame(z = c(0.1, 0.2, 0.3, 0.4), cell_id = 1:4))
  hundred <- .ps_points(c(d$x, 100, 200, 300, 400), c(d$y, 100, 200, 300, 400), d100)
  expect_identical(nrow(hundred), 100L)
  expect_false("spatialkit_few_residuals" %in%
                 .ps_classes(.ps_catch(summarize_by_cell(hundred, "z", id_col = "cell_id",
                                                         estimate = "pooled", obs_level = 0.9))))
})


test_that("a level close to 1 is written with its digits, and the count it needs as a number", {
  L <- .ps_layer()
  seen <- list()
  r <- withCallingHandlers(
    .ps_catch(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled",
                                obs_level = 0.999999999)),
    warning = function(w) seen[[length(seen) + 1L]] <<- w)
  # One warning, of the class: none of R's own about a number that does not
  # fit a whole number.
  expect_identical(.ps_classes(r), "spatialkit_few_residuals")
  w <- r$warnings[[1L]]
  expect_identical(.ps_fields(w),
                   list(columns = "z", needed = 1e10, level = 0.999999999, n = c(z = 900L)))
  expect_identical(
    conditionMessage(w),
    paste0("summarize_by_cell(): the 99.9999999% interval of 'z' rests on fewer than ",
           "10000000000 leave-one-out residuals (900), the number at which five are expected ",
           "beyond each end; its ends are unstable."))
  expect_false(grepl("NA", conditionMessage(w), fixed = TRUE))
  txt <- utils::capture.output(print(attr(r$value, "pooled")))
  expect_identical(
    txt[9L],
    "  interval    : 99.9999999% for ONE OBSERVATION (..obs_lo_resp_z, ..obs_hi_resp_z):")
  # The offsets are the quantiles at the level, not the range of the residuals.
  xy <- unname(sf::st_coordinates(L$asg)[, 1:2])
  u  <- sort(unique(L$asg$cell_id))
  e  <- .pool_fit(L$asg$z, match(L$asg$cell_id, u), length(u), xy, want = "pooled")$resid_pooled
  off <- attr(r$value, "pooled")$columns$resp_z$offsets
  expect_identical(off, stats::quantile(e, c(5e-10, 1 - 5e-10), names = FALSE, type = 7))
  expect_gt(off[1L], min(e)); expect_lt(off[2L], max(e))
  # Where the count fits a whole number it is one.
  r <- .ps_catch(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled",
                                   obs_level = 0.999999))
  expect_identical(r$warnings[[1L]]$needed, 10000000L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "the 99.9999% interval of 'z' rests on fewer than 10000000 leave-one-out",
               fixed = TRUE)
  # A level that is 1 or 0 to the twelve decimals the interval reads is
  # refused, as 1 and 0 are, before any summary.
  for (bad in c(1 - 1e-13, 1e-13))
    expect_error(summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled",
                                   obs_level = bad),
                 "`obs_level` must be a single number strictly between 0 and 1", fixed = TRUE)
  # The nearest levels the interval can tell from them are taken.
  for (ok in c(1 - 1e-11, 1e-11))
    expect_false(anyNA(suppressWarnings(
      summarize_by_cell(L$asg, "z", id_col = "cell_id", estimate = "pooled",
                        obs_level = ok))$..obs_lo_resp_z))
})


# ---- Cost ------------------------------------------------------------------------------

test_that("the pooled path adds little to a summary of 100,000 rows", {
  skip_on_cran()
  skip_if_not_installed("FNN")
  # 100,000 rows in a 149 x 148 grid of cells, about 22,000, one column.
  big <- local({
    .ps_seed(2008)
    x <- stats::runif(1e5, 0, 1e6); y <- stats::runif(1e5, 0, 1e6)
    pts <- .ps_points(x, y, data.frame(z = sin(x / 2e5) + cos(y / 1.5e5) + stats::rnorm(1e5)))
    pts$cell_id <- as.integer(148 * floor(x / (1e6 / 149)) + floor(y / (1e6 / 148)) + 1)
    pts
  })
  cells <- length(unique(big$cell_id))
  expect_gt(cells, 21500L); expect_lt(cells, 22100L)
  secs <- function(est) min(vapply(1:2, function(i) system.time(
    suppressWarnings(summarize_by_cell(big, "z", id_col = "cell_id", estimate = est)))[["elapsed"]],
    numeric(1)))
  plain  <- secs("mean")
  pooled <- secs("pooled")
  expect_lt(pooled - plain, 2)
})
