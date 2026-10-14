# tests/testthat/test-core-values-at.R
# ---------------------------------------------------------------------------
# cell_values_at(): the value and interval of a cell at other locations.
# Which cell a location lies in (its nearest seed, or the polygon that holds
# it), the value it takes, the fill for a location with no cell or a cell
# without a value, the distance from the nearest observation, the two tie
# rules, the record of the call, and the coverage of the interval at
# held-out rows of simulated layers.
#
# The reference for the values at held-out rows is predict_partition() of
# the experiments (helper-pool-reference.R).
#
# Every profile names `partition`, `init`, `nstart`, `estimate`,
# `obs_per_cell` and `sac`, every selection its `criterion`, every summary
# its `estimate`.  The regions are of the density kind from 3 k-means++
# starts under seed 123.  A simulated layer is drawn under the seed its
# test states: the locations (x, then y), then the response, then the rows
# held out.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.va_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# A call with its warnings kept and muffled.
.va_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.va_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.va_error <- function(expr) tryCatch({ expr; NULL }, error = function(e) e)
.va_print <- function(x) utils::capture.output(print(x))
# The columns of a frame and nothing else: no row names, no record.
.va_cols <- function(d) {
  d <- as.list(d)
  attributes(d) <- list(names = names(d))
  d
}
.va_no_fnn <- function(env = parent.frame())
  testthat::local_mocked_bindings(.has_backend = function(pkg) FALSE, .env = env)

.va_store <- new.env(parent = emptyenv())
.va_keep <- function(key, expr) {
  if (is.null(.va_store[[key]])) .va_store[[key]] <- expr
  .va_store[[key]]
}

# `n` locations uniform on a square of side 1e5 with its corner at (5e5,
# 5e6) in EPSG:32617; the response `z` is sin(x / 2e4) + cos(y / 1.5e4) plus
# normal noise of standard deviation 0.6, and `w` is 2 z + 1; `hold` rows are
# then held out by sample.int().
.va_draw <- function(n, seed, hold) {
  d <- local({
    .va_seed(seed)
    x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
    z <- sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 0.6)
    te <- sample.int(n, hold)
    list(d = data.frame(x = 5e5 + x, y = 5e6 + y, z = z, w = 2 * z + 1), te = te)
  })
  pts <- sf::st_as_sf(d$d, coords = c("x", "y"), crs = 32617)
  list(rows = pts[-d$te, ], held = pts[d$te, ])
}

# The calls of the design's first block, from the profile to the summary:
# the count from a size of 4.5, the cells, the rows in them, and the two
# summaries, pooled values with both intervals at 0.90 and plain means.
.va_path <- function(rows) {
  prof <- suppressWarnings(
    resolution_profile(rows, "z", partition = "density", init = "kmeans++", nstart = 3,
                       sac = "auto", seed = 123, estimate = "pooled", obs_per_cell = 4.5))
  sel   <- select_resolution(prof, criterion = "size")
  seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = rows)
  tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  asg   <- assign_features_to_polygons(rows, tess$cells, polygon_id_col = "cell_id",
                                       quiet = TRUE)
  pooled <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                              estimate = "pooled", obs_level = 0.90, conf_level = 0.90)
  # Plain means on cells sized for pooled values: the summary says so.
  plain <- .va_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                       estimate = "mean"))
  stopifnot(identical(.va_classes(plain), "spatialkit_estimate_mismatch"))
  seeds$cell_id <- tess$index
  list(sel = sel, seeds = seeds, tess = tess, asg = asg, pooled = pooled,
       plain = plain$value)
}

# 600 locations (seed 6102), 120 held out; every row has a response.
.va_full <- function() .va_keep("full", {
  L <- .va_draw(600, 6102, 120)
  c(L, .va_path(L$rows))
})
# 1,500 locations (seed 6101), 300 held out; the other 1,200 have no response
# east of x = 580,000.  The centres are fitted to all 1,200 and the size
# divides the rows with a response, so the centres in the east hold no
# scored row, and the held-out rows there lie in cells without a value.
.va_gap <- function() .va_keep("gap", {
  L <- .va_draw(1500, 6101, 300)
  L$rows$z[sf::st_coordinates(L$rows)[, 1L] > 5.8e5] <- NA
  c(L, .va_path(L$rows))
})

.va_at <- function(L, newdata, values = L$pooled, seeds = L$seeds, points = L$asg,
                   response = "z")
  cell_values_at(newdata, values, points, response, id_col = "cell_id", seeds = seeds)

# The coordinates of a layer and, for locations `q`, the nearest of the
# locations `u` by the rule of the design: the least squared distance, a tie
# to the lower x and then the lower y.
.va_xy <- function(layer) unname(sf::st_coordinates(layer)[, 1:2, drop = FALSE])
.va_nearest <- function(q, u) {
  vapply(seq_len(nrow(q)), function(i) {
    d2 <- (q[i, 1L] - u[, 1L])^2 + (q[i, 2L] - u[, 2L])^2
    at <- which(d2 == min(d2))
    at[order(u[at, 1L], u[at, 2L])][1L]
  }, integer(1))
}

# A layer made by hand: points at (5e5 + x, 5e6 + y) in EPSG:32617.
.va_points <- function(x, y, ..., crs = 32617)
  sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, ...), coords = c("x", "y"), crs = crs)
.va_means <- function(rows)
  summarize_by_cell(rows, "z", id_col = "cell_id", estimate = "mean")

# Scored rows on a 6 x 5 lattice of spacing 10, two rows at each location.
# The location that is `i`-th in (x, y) order is cell 31 - i, so neither the
# order of the IDs nor that of the rows follows the coordinates; the rows
# are shuffled under seed 6201.  The response of cell `k` is 10 k -+ 1.
.va_lattice <- function() .va_keep("lattice", {
  g <- expand.grid(y = seq(0, 40, by = 10), x = seq(0, 50, by = 10))
  cell <- 31L - seq_len(30L)
  d <- data.frame(x = rep(g$x, 2L), y = rep(g$y, 2L), cell_id = rep(cell, 2L),
                  z = 10 * rep(cell, 2L) + rep(c(-1, 1), each = 30L))
  d <- local({ .va_seed(6201); d[sample.int(60L), ] })
  rows <- .va_points(d$x, d$y, cell_id = d$cell_id, z = d$z)
  list(rows = rows, values = .va_means(rows), at = g, cell = cell)
})


# ---- The pieces ----------------------------------------------------------------

test_that("the columns of the summary that hold a column's value: the response's first, the pooled value first", {
  expect_null(.values_at_columns(c("cell_id", "n", "resp_median_z"), "z"))
  expect_identical(.values_at_columns(c("cell_id", "resp_mean_z"), "z"),
                   list(value = "resp_mean_z", estimate = "mean", lo = NULL, hi = NULL))
  expect_identical(
    .values_at_columns(c("resp_mean_z", "..pooled_resp_z", "pred_mean_z", "..pooled_pred_z"), "z"),
    list(value = "..pooled_resp_z", estimate = "pooled", lo = NULL, hi = NULL))
  expect_identical(
    .values_at_columns(c("pred_mean_z", "resp_mean_z", "..pooled_pred_z"), "z")$value,
    "resp_mean_z")
  expect_identical(
    .values_at_columns(c("pred_mean_w", "..pooled_pred_w", "..obs_lo_pred_w", "..obs_hi_pred_w",
                         "..obs_lo_resp_z", "..obs_hi_resp_z"), "w"),
    list(value = "..pooled_pred_w", estimate = "pooled", lo = "..obs_lo_pred_w",
         hi = "..obs_hi_pred_w"))
  # One end of the interval alone is no interval; an interval beside a plain mean is not read.
  expect_null(.values_at_columns(c("..pooled_resp_z", "..obs_lo_resp_z"), "z")$lo)
  expect_null(.values_at_columns(c("resp_mean_z", "..obs_lo_resp_z", "..obs_hi_resp_z"), "z")$lo)
  # The column is matched whole.
  expect_null(.values_at_columns(c("resp_mean_zz", "..pooled_resp_zz"), "z"))
})

test_that("IDs are written in the class of the summary's ID column, also those the summary has no row for", {
  expect_identical(.values_at_ids(c("3", NA, "1", "7"), c(1L, 2L, 3L)), c(3L, NA, 1L, 7L))
  # A whole number held as a double is matched by its digits: 1e5 is "100000".
  expect_identical(.values_at_ids(c("100000", "250000"), c(1e5, 2e5)), c(1e5, 250000))
  expect_identical(.values_at_ids(c("b", NA, "zz"), c("a", "b")), c("b", NA, "zz"))
  f <- .values_at_ids(c("b", NA, "zz", "a"), factor(c("b", "a"), levels = c("b", "a")))
  expect_identical(f, factor(c("b", NA, "zz", "a"), levels = c("b", "a", "zz")))
  # A missing ID in the summary (a row of unassigned points) is no cell.
  expect_identical(.values_at_ids(c(NA, "2"), c(NA, 2L)), c(NA, 2L))
  # A class an ID cannot be written in from its text: as text.
  d <- as.Date(c("2026-01-01", "2026-01-02"))
  expect_identical(.values_at_ids(c("2026-01-02", NA), d), d[c(2L, NA)])
  expect_identical(.values_at_ids(c("2026-01-02", "x"), d), c("2026-01-02", "x"))
  # An ordered factor stays one, with the new level after its own.
  f <- .values_at_ids(c("lo", "zz", NA), factor(c("lo", "hi"), levels = c("lo", "hi"),
                                                 ordered = TRUE))
  expect_identical(f, factor(c("lo", "zz", NA), levels = c("lo", "hi", "zz"), ordered = TRUE))
  expect_false(is.ordered(.values_at_ids("zz", factor(c("lo", "hi")))))
  # A number the summary has no row for is the number its text says: a fraction in a
  # double column is not cut to its whole part.
  expect_identical(.values_at_ids(c("3.5", "2", "-0.25"), c(1, 2)), c(3.5, 2, -0.25))
  # An ID the class of the column cannot hold is not turned into NA or into another
  # ID: every ID is then returned as text.
  expect_identical(.values_at_ids(c("new-7", "3", NA), 1:5), c("new-7", "3", NA))
  expect_identical(.values_at_ids(c("3.5", "3"), 1:5), c("3.5", "3"))
  expect_identical(.values_at_ids(c("new-7", "2"), c(1, 2)), c("new-7", "2"))
  expect_identical(.values_at_ids(c("1e5", "2"), c(1, 2)), c("1e5", "2"))
  expect_identical(.values_at_ids(c("3000000000", "2"), 1:5), c("3000000000", "2"))
})

test_that("the nearest centre is .nearest_centre()'s, a tie to the lower number, however many points tie", {
  skip_if_not_installed("FNN")
  # A lattice of centres, and points on it, half-way between two and between four.
  g <- expand.grid(y = seq(0, 90, by = 10), x = seq(0, 110, by = 10))
  cx <- g$x; cy <- g$y
  q <- rbind(cbind(cx, cy), cbind(cx + 5, cy), cbind(cx, cy + 5), cbind(cx + 5, cy + 5),
             cbind(cx + 3, cy + 1))
  exact <- .nearest_exact(q[, 1L], q[, 2L], cx, cy)
  expect_identical(.values_at_nearest(q[, 1L], q[, 2L], cx, cy, TRUE), exact)
  expect_identical(.values_at_nearest(q[, 1L], q[, 2L], cx, cy, FALSE), exact)
  expect_identical(.nearest_centre(q[, 1L], q[, 2L], cx, cy, use_fnn = TRUE), exact)
  # The centres in another order: the lower number is another centre.
  o <- rev(seq_along(cx))
  expect_identical(.values_at_nearest(q[, 1L], q[, 2L], cx[o], cy[o], TRUE),
                   .nearest_exact(q[, 1L], q[, 2L], cx[o], cy[o]))
  expect_false(identical(o[.values_at_nearest(q[, 1L], q[, 2L], cx[o], cy[o], TRUE)], exact))
  # Centres that share a point, and centres all on one point.
  local({
    .va_seed(6301)
    for (K in c(17L, 40L, 300L)) {
      ax <- round(stats::runif(K, 0, 12)); ay <- round(stats::runif(K, 0, 12))
      px <- round(stats::runif(400, 0, 24)) / 2; py <- round(stats::runif(400, 0, 24)) / 2
      expect_identical(.values_at_nearest(px, py, ax, ay, TRUE), .nearest_exact(px, py, ax, ay))
    }
  })
  expect_identical(.values_at_nearest(c(1, 9), c(2, 8), rep(5, 40), rep(5, 40), TRUE), c(1L, 1L))
  # A tied point is settled from a shortlist, not by a search of every centre.
  searched <- 0L
  real <- .nearest_exact
  force(real)
  testthat::local_mocked_bindings(.nearest_exact = function(x, y, cx, cy) {
    searched <<- searched + length(x)
    real(x, y, cx, cy)
  })
  expect_identical(.values_at_nearest(q[, 1L], q[, 2L], cx, cy, TRUE), exact)
  expect_identical(searched, 0L)
})

test_that("the nearest centre is the exact one where centres are nearly the same distance away", {
  skip_if_not_installed("FNN")
  # Among 40 centres, centre 3 is 10 from the query and centre 2 further by two parts in
  # 1e10: the nearer one has the higher number.
  cen <- local({
    .va_seed(6302)
    list(x = stats::runif(40, 50, 90), y = stats::runif(40, 50, 90))
  })
  cen$x[2L] <- 10 * (1 + 2e-10); cen$y[2L] <- 0
  cen$x[3L] <- 0;                cen$y[3L] <- 10
  expect_identical(.nearest_exact(0, 0, cen$x, cen$y), 3L)
  expect_identical(.values_at_nearest(0, 0, cen$x, cen$y, TRUE), 3L)
  expect_identical(.values_at_nearest(0, 0, cen$x, cen$y, FALSE), 3L)
  # A ring of 60 centres about the query at radii 100 (1 + j 1e-11), j from 0 to 59 in
  # an order drawn under seed 6303, alone and among 40 centres far away.  The nearest is
  # the one with j = 0, which is not the first.
  j <- local({ .va_seed(6303); sample.int(60L) - 1L })
  angle <- seq(0, 2 * pi, length.out = 61L)[-61L]
  rx <- 100 * (1 + j * 1e-11) * cos(angle)
  ry <- 100 * (1 + j * 1e-11) * sin(angle)
  expect_gt(which(j == 0L), 1L)
  expect_lt(j[1L], 50L)
  qx <- c(0, 3); qy <- c(0, 4)
  for (far in list(numeric(0), 1000 + cen$x)) {
    ax <- c(rx, far); ay <- c(ry, far)
    exact <- .nearest_exact(qx, qy, ax, ay)
    expect_identical(exact[1L], which(j == 0L))
    expect_identical(.values_at_nearest(qx, qy, ax, ay, TRUE), exact)
    expect_identical(.values_at_nearest(qx, qy, ax, ay, FALSE), exact)
  }
})

test_that("a kd-tree that ranks nearly tied centres in another order than R gives the same nearest centre", {
  skip_if_not_installed("FNN")
  # Three centres at 1, 1 + 2e-12 and 1 + 4e-12 from the query, among 33.  The stand-in
  # for FNN's search returns, for a shortlist of two, the second and the third nearest,
  # as a tree would whose rounding ranked the nearest of the three last.
  cen <- local({
    .va_seed(6304)
    list(x = c(1, 1 + 2e-12, 1 + 4e-12, stats::runif(30, 50, 90)),
         y = c(0, 0, 0, stats::runif(30, 50, 90)))
  })
  expect_identical(.nearest_exact(0, 0, cen$x, cen$y), 1L)
  real <- FNN::get.knnx
  force(real)
  asked <- integer(0)
  testthat::local_mocked_bindings(get.knnx = function(data, query, k, ...) {
    asked <<- c(asked, as.integer(k))
    r <- real(data, query, k = min(nrow(data), k + 1L), ...)
    keep <- if (k == 2L) 2:3 else seq_len(k)
    list(nn.index = r$nn.index[, keep, drop = FALSE], nn.dist = r$nn.dist[, keep, drop = FALSE])
  }, .package = "FNN")
  expect_identical(.values_at_nearest(0, 0, cen$x, cen$y, TRUE), 1L)
  # The two centres of the first shortlist are too close to settle the point; the next
  # shortlist, of four, holds the nearest.
  expect_identical(asked, c(2L, 4L))
})


# ---- At the training rows ------------------------------------------------------

test_that("at the training rows each row gets its own cell's value, unfilled and at distance 0", {
  L <- .va_full()
  expect_true("cell_id" %in% names(L$asg))
  own <- match(L$asg$cell_id, L$pooled$cell_id)
  for (seeds in list(NULL, L$seeds)) {
    r <- .va_catch(.va_at(L, L$asg, seeds = seeds))
    # No warning: none about the ID column `newdata` carries, and no other.
    expect_identical(.va_classes(r), character(0))
    at <- r$value
    expect_identical(nrow(at), nrow(L$asg))
    expect_identical(names(at), c("cell_id", "value", "estimate", "obs_lo", "obs_hi", "filled",
                                  "filled_from", "dist_nearest", "dist_ratio"))
    expect_identical(at$cell_id, L$asg$cell_id)
    expect_identical(at$value, L$pooled$..pooled_resp_z[own])
    expect_identical(at$obs_lo, L$pooled$..obs_lo_resp_z[own])
    expect_identical(at$obs_hi, L$pooled$..obs_hi_resp_z[own])
    expect_identical(at$estimate, rep("pooled", nrow(at)))
    expect_identical(at$filled, rep(FALSE, nrow(at)))
    expect_identical(at$filled_from, rep(NA_integer_, nrow(at)))
    expect_identical(at$dist_nearest, rep(0, nrow(at)))
    expect_identical(at$dist_ratio, rep(0, nrow(at)))
    # The plain means of the same cells.
    p <- .va_catch(.va_at(L, L$asg, values = L$plain, seeds = seeds))
    expect_identical(.va_classes(p), character(0))
    expect_identical(names(p$value), c("cell_id", "value", "estimate", "filled", "filled_from",
                                       "dist_nearest", "dist_ratio"))
    expect_identical(p$value$value, L$plain$resp_mean_z[match(L$asg$cell_id, L$plain$cell_id)])
    expect_identical(p$value$estimate, rep("mean", nrow(at)))
    expect_identical(p$value$filled, rep(FALSE, nrow(at)))
  }
  # The assignment itself does warn about that column, which is why it is not handed over.
  expect_warning(assign_features_to_polygons(L$asg[1:5, ], L$pooled, polygon_id_col = "cell_id",
                                             keep_unassigned = TRUE, quiet = TRUE),
                 "already carries column")
})

test_that("the same holds on lon/lat rows, measured in the CRS the count was scored in", {
  L <- .va_keep("lonlat", {
    F <- .va_full()
    rows <- sf::st_transform(F$rows, 4326)
    c(list(rows = rows, held = sf::st_transform(F$held, 4326)), .va_path(rows))
  })
  rec <- attr(L$pooled, "count_rule")
  expect_true(isTRUE(sf::st_is_longlat(L$asg)))
  expect_false(isTRUE(sf::st_is_longlat(rec$crs)))
  own <- match(L$asg$cell_id, L$pooled$cell_id)
  for (seeds in list(NULL, L$seeds)) {
    r <- .va_catch(.va_at(L, L$asg, seeds = seeds))
    expect_identical(.va_classes(r), character(0))
    expect_identical(r$value$cell_id, L$asg$cell_id)
    expect_identical(r$value$value, L$pooled$..pooled_resp_z[own])
    expect_identical(r$value$filled, rep(FALSE, nrow(L$asg)))
    expect_identical(r$value$dist_nearest, rep(0, nrow(L$asg)))
    expect_identical(attr(r$value, "values_at")$crs, rec$crs)
  }
  # Distances are metres of that CRS, not degrees.
  at <- .va_at(L, L$held)
  expect_gt(stats::median(at$dist_nearest), 100)
  expect_gt(attr(at, "values_at")$spacing, 100)
})


# ---- Held-out rows, against the experiments' function ---------------------------

test_that("with seeds the values at held-out rows are predict_partition()'s, the rows of empty regions included", {
  skip_if_not_installed("FNN")
  L <- .va_gap()
  rec <- attr(L$pooled, "count_rule")
  expect_identical(sf::st_crs(L$asg), rec$crs)
  expect_false(anyNA(L$tess$index))
  scored <- is.finite(L$asg$z)
  cen   <- .va_xy(L$seeds)
  xy_tr <- .va_xy(L$asg)[scored, ]
  xy_te <- .va_xy(L$held)
  g_tr  <- .nearest_centre(xy_tr[, 1L], xy_tr[, 2L], cen[, 1L], cen[, 2L], use_fnn = TRUE)
  g_te  <- .nearest_centre(xy_te[, 1L], xy_te[, 2L], cen[, 1L], cen[, 2L], use_fnn = TRUE)
  # The cells of the summary are the regions of these centres.
  expect_identical(L$asg$cell_id[scored], L$tess$index[g_tr])
  ref <- predict_partition(g_tr, g_te, L$asg$z[scored], xy_tr, xy_te, want = c("P", "L"))
  empty <- !(g_te %in% g_tr)
  expect_identical(ref$n_empty, sum(empty))
  # Centres in the east hold no scored row, and held-out rows lie in their regions.
  expect_gt(ref$n_empty, 20L)
  expect_lt(ref$n_empty, 150L)
  # The cell of the nearest scored row, by the experiments' own expression.
  from <- L$tess$index[g_tr[FNN::get.knnx(xy_tr, xy_te[empty, , drop = FALSE],
                                          k = 1L)$nn.index[, 1L]]]

  for (est in c("L", "P")) {
    values <- if (est == "L") L$pooled else L$plain
    r <- .va_catch(.va_at(L, L$held, values = values))
    expect_identical(.va_classes(r), character(0))
    at <- r$value
    expect_lt(max(abs(at$value - ref$pred[, est]) / pmax(1, abs(ref$pred[, est]))), 1e-12)
    expect_identical(sum(at$filled), ref$n_empty)
    expect_identical(at$filled, empty)
    expect_identical(at$cell_id, L$tess$index[g_te])
    expect_identical(at$filled_from[empty], from)
    expect_true(all(is.na(at$filled_from[!empty])))
    expect_identical(at$estimate, rep(if (est == "L") "pooled" else "mean", nrow(at)))
    expect_identical(attr(at, "values_at")$filled, ref$n_empty)
    # A filled row lies in one cell and takes the value of another.
    expect_true(all(at$cell_id[empty] != at$filled_from[empty]))
    expect_identical(at$value[empty],
                     values[[if (est == "L") "..pooled_resp_z" else "resp_mean_z"]][
                       match(from, values$cell_id)])
    # The answer does not depend on the FNN package.
    local({
      .va_no_fnn()
      expect_identical(.va_at(L, L$held, values = values), at)
    })
  }
})

test_that("the cell, the fill and its source are those of a search of every pair, with and without FNN", {
  # The rule written out, without the experiments' function and without FNN: the cell of
  # the nearest seed; filled when that cell holds no scored row; from the cell of the
  # nearest scored row.  No location of this layer is the same distance from two.
  L <- .va_gap()
  scored <- is.finite(L$asg$z)
  q <- .va_xy(L$held)
  cell <- L$tess$index[.va_nearest(q, .va_xy(L$seeds))]
  holds <- cell %in% L$asg$cell_id[scored]
  from  <- L$asg$cell_id[scored][.va_nearest(q, .va_xy(L$asg)[scored, ])]
  expect_gt(sum(!holds), 20L)
  for (values in list(L$pooled, L$plain)) {
    column <- if ("..pooled_resp_z" %in% names(values)) "..pooled_resp_z" else "resp_mean_z"
    for (fnn in c(TRUE, FALSE)) local({
      if (!fnn) .va_no_fnn()
      at <- .va_at(L, L$held, values = values)
      expect_identical(at$cell_id, cell)
      expect_identical(at$filled, !holds)
      expect_identical(at$filled_from, ifelse(holds, NA_integer_, from))
      expect_identical(at$value, values[[column]][match(ifelse(holds, cell, from), values$cell_id)])
      expect_false(anyNA(at$value))
    })
  }
})

test_that("the interval of a held-out row is its value plus the summary's two offsets, filled or not", {
  L <- .va_gap()
  at <- .va_at(L, L$held)
  off <- attr(L$pooled, "pooled")$columns$resp_z$offsets
  expect_true(all(is.finite(off)) && off[1L] < 0 && off[2L] > 0)
  expect_true(any(at$filled))
  expect_equal(at$obs_lo - at$value, rep(off[1L], nrow(at)), tolerance = 1e-12)
  expect_equal(at$obs_hi - at$value, rep(off[2L], nrow(at)), tolerance = 1e-12)
  src <- match(ifelse(at$filled, at$filled_from, at$cell_id), L$pooled$cell_id)
  expect_identical(at$obs_lo, L$pooled$..obs_lo_resp_z[src])
  expect_identical(at$obs_hi, L$pooled$..obs_hi_resp_z[src])
})


# ---- Without seeds: the polygon that holds the location ---------------------------

test_that("without seeds a location takes the value of the polygon that holds it, and one in no polygon is filled", {
  L <- .va_gap()
  # The held-out rows and three locations outside every polygon: west, south and far north-east.
  far <- .va_points(c(-2e4, 5e4, 3e5), c(5e4, -3e4, 3e5))
  newdata <- rbind(sf::st_sf(geometry = sf::st_geometry(L$held)),
                   sf::st_sf(geometry = sf::st_geometry(far)))
  outside <- c(rep(FALSE, nrow(L$held)), rep(TRUE, 3L))
  r <- .va_catch(.va_at(L, newdata, seeds = NULL))
  expect_identical(.va_classes(r), character(0))
  at <- r$value
  by_seeds <- .va_at(L, newdata)
  inside <- !is.na(at$cell_id)
  # Every location inside a polygon: the value, the cell and the fill of the seeds' call.
  expect_gt(sum(inside), 290L)
  expect_identical(.va_cols(at[inside, ]), .va_cols(by_seeds[inside, ]))
  expect_identical(at$dist_nearest, by_seeds$dist_nearest)
  # A location in no polygon: no cell, filled, from the cell of its nearest scored row.
  expect_identical(inside[outside], rep(FALSE, 3L))
  expect_true(all(at$filled[!inside]))
  scored <- is.finite(L$asg$z)
  u <- .va_xy(L$asg)[scored, ]
  nearest <- .va_nearest(.va_xy(newdata)[!inside, , drop = FALSE], u)
  expect_identical(at$filled_from[!inside], L$asg$cell_id[scored][nearest])
  expect_identical(at$value[!inside],
                   L$pooled$..pooled_resp_z[match(at$filled_from[!inside], L$pooled$cell_id)])
  q <- .va_xy(newdata)[!inside, , drop = FALSE]
  expect_identical(at$dist_nearest[!inside],
                   sqrt((q[, 1L] - u[nearest, 1L])^2 + (q[, 2L] - u[nearest, 2L])^2))
  # A location in a cell with no scored row keeps that cell's ID and has another in filled_from.
  kept <- inside & at$filled
  expect_gt(sum(kept), 20L)
  expect_false(anyNA(at$cell_id[kept]))
  expect_true(all(at$cell_id[kept] != at$filled_from[kept]))
  expect_true(all(is.na(L$pooled$..pooled_resp_z[match(at$cell_id[kept], L$pooled$cell_id)])))
  expect_false(any(at$cell_id[kept] %in% L$asg$cell_id[scored]))
  expect_true(all(at$filled_from[kept] %in% L$asg$cell_id[scored]))
  expect_false(anyNA(at$value[kept]))
  expect_identical(attr(at, "values_at")$assignment, "polygons")
  expect_identical(attr(at, "values_at")$filled, sum(at$filled))
})

test_that("with every location outside every polygon the call returns, all filled, and no warning leaves it", {
  L <- .va_full()
  far <- .va_points(c(-2e4, -3e4, 3e5, 5e4, 5e4), c(5e4, -3e4, 3e5, -1e4, 2e5))
  # The assignment warns of exactly this.
  expect_warning(
    assign_features_to_polygons(far, L$pooled, polygon_id_col = "cell_id",
                                keep_unassigned = TRUE, quiet = TRUE),
    class = "spatialkit_unassigned")
  r <- .va_catch(.va_at(L, far, seeds = NULL))
  expect_identical(.va_classes(r), character(0))
  at <- r$value
  expect_identical(nrow(at), 5L)
  expect_identical(at$cell_id, rep(NA_integer_, 5L))
  expect_identical(at$filled, rep(TRUE, 5L))
  u <- .va_xy(L$asg)
  nearest <- .va_nearest(.va_xy(far), u)
  expect_identical(at$filled_from, L$asg$cell_id[nearest])
  expect_identical(at$value, L$pooled$..pooled_resp_z[match(at$filled_from, L$pooled$cell_id)])
  expect_false(anyNA(at$value))
  q <- .va_xy(far)
  expect_identical(at$dist_nearest,
                   sqrt((q[, 1L] - u[nearest, 1L])^2 + (q[, 2L] - u[nearest, 2L])^2))
  expect_true(all(at$dist_ratio > 1))
  expect_identical(attr(at, "values_at")$filled, 5L)
})


# ---- IDs of two classes -------------------------------------------------------------

test_that("IDs are matched as text: a summary with character IDs beside rows and seeds with integer ones", {
  L <- .va_gap()
  cells <- L$tess$cells
  cells$cell_id <- as.character(cells$cell_id)
  values <- summarize_by_cell(L$asg, "z", id_col = "cell_id", cells_sf = cells,
                              estimate = "pooled", obs_level = 0.90)
  expect_type(values$cell_id, "character")
  expect_type(L$asg$cell_id, "integer")
  expect_type(L$seeds$cell_id, "integer")
  same <- .va_at(L, L$held)
  expect_type(same$cell_id, "integer")
  expect_true(any(same$filled))
  for (seeds in list(L$seeds, NULL)) {
    want <- .va_at(L, L$held, seeds = seeds)
    at   <- .va_at(L, L$held, values = values, seeds = seeds)
    expect_type(at$cell_id, "character")
    expect_type(at$filled_from, "character")
    expect_identical(at$cell_id, as.character(want$cell_id))
    expect_identical(at$filled_from, as.character(want$filled_from))
    expect_equal(at$value, want$value, tolerance = 1e-12)
    expect_false(anyNA(at$value))
    expect_identical(at[c("filled", "dist_nearest", "dist_ratio")],
                     want[c("filled", "dist_nearest", "dist_ratio")])
  }
  # The other way round: character seeds beside a summary and rows with integer IDs.
  chr <- L$seeds
  chr$cell_id <- as.character(chr$cell_id)
  expect_identical(.va_at(L, L$held, seeds = chr), same)
})

test_that("a whole number held as a double is matched by its digits, whichever layer holds it so", {
  L <- .va_gap()
  same <- .va_at(L, L$held)
  # A cell that holds held-out rows and that filled rows take their value from gets the
  # ID 100000, which as.character() writes in another way for a double than for an integer.
  expect_identical(as.character(1e5), "1e+05")
  expect_identical(as.character(100000L), "100000")
  k0 <- intersect(same$filled_from[same$filled], same$cell_id[!same$filled])[1L]
  expect_false(is.na(k0))
  shift <- function(x) x + (100000L - k0)
  int <- function(layer) { layer$cell_id <- shift(layer$cell_id); layer }
  dbl <- function(layer) { layer$cell_id <- as.numeric(shift(layer$cell_id)); layer }
  cases <- list(values = list(dbl, int, int), points = list(int, dbl, int),
                seeds = list(int, int, dbl), all = list(dbl, dbl, dbl))
  for (nm in names(cases)) {
    f <- cases[[nm]]
    at <- .va_at(L, L$held, values = f[[1L]](L$pooled), points = f[[2L]](L$asg),
                 seeds = f[[3L]](L$seeds))
    expect_type(at$cell_id, if (nm %in% c("values", "all")) "double" else "integer")
    expect_type(at$filled_from, typeof(at$cell_id))
    expect_identical(as.integer(at$cell_id), shift(same$cell_id))
    expect_identical(as.integer(at$filled_from), shift(same$filled_from))
    expect_true(any(at$cell_id == 1e5) && any(at$filled_from == 1e5, na.rm = TRUE))
    expect_identical(.va_cols(at[-c(1L, 7L)]), .va_cols(same[-c(1L, 7L)]))
  }
  # By the polygons, whose IDs the assignment hands back in the class the summary holds.
  at <- .va_at(L, L$held, values = dbl(L$pooled), points = int(L$asg), seeds = NULL)
  want <- .va_at(L, L$held, seeds = NULL)
  expect_type(at$cell_id, "double")
  expect_identical(as.integer(at$cell_id), shift(want$cell_id))
  expect_identical(.va_cols(at[-c(1L, 7L)]), .va_cols(want[-c(1L, 7L)]))
})

test_that("a cell the summary has no row for keeps its ID, in the class of the summary's IDs, and is filled", {
  L <- .va_full()
  # One more seed on a held-out location, for a cell no row was assigned to.
  extra <- sf::st_sf(seed_id = 9999L, method = "by hand", cell_id = 9999L,
                     geometry = sf::st_geometry(L$held)[7L])
  seeds <- rbind(L$seeds[c("seed_id", "method", "cell_id")], extra)
  at <- .va_at(L, L$held, seeds = seeds)
  expect_identical(at$cell_id[7L], 9999L)
  # The locations nearest the new seed: that one, and those around it.
  mine <- at$cell_id == 9999L
  expect_lt(sum(mine), 10L)
  expect_true(all(at$filled[mine]))
  expect_false(anyNA(at$filled_from[mine]))
  expect_identical(at$value[mine],
                   L$pooled$..pooled_resp_z[match(at$filled_from[mine], L$pooled$cell_id)])
  expect_identical(.va_cols(at[!mine, ]), .va_cols(.va_at(L, L$held)[!mine, ]))
  # A summary without the cells' polygons, with IDs as text.
  values <- sf::st_drop_geometry(L$pooled)
  values$cell_id <- as.character(values$cell_id)
  chr <- .va_at(L, L$held, values = values, seeds = seeds)
  expect_identical(chr$cell_id[7L], "9999")
  expect_identical(chr$cell_id, as.character(at$cell_id))
  expect_identical(chr$value, at$value)
  # A cell the summary has no row for has no value to give: a location filled from it
  # gets NA.
  drop <- at$filled_from[7L]
  less <- .va_at(L, L$held, values = L$pooled[L$pooled$cell_id != drop, ], seeds = seeds)
  expect_identical(less$filled_from[7L], drop)
  expect_true(is.na(less$value[7L]))
  # A location in such a cell is not filled when the cell holds scored rows: it has NA.
  busy <- as.integer(names(which.max(table(at$cell_id[!mine]))))
  hit  <- at$cell_id == busy
  expect_gt(sum(hit), 1L)
  less <- .va_at(L, L$held, values = L$pooled[L$pooled$cell_id != busy, ], seeds = seeds)
  expect_identical(less$cell_id, at$cell_id)
  expect_identical(less$filled[hit], rep(FALSE, sum(hit)))
  expect_true(all(is.na(less$value[hit])) && all(is.na(less$filled_from[hit])))
  expect_identical(.va_cols(less[!hit & !mine, ]), .va_cols(at[!hit & !mine, ]))
  # Seeds whose IDs the summary's integer column cannot hold: both ID columns are text,
  # and no ID is lost.
  named <- seeds
  named$cell_id <- as.character(named$cell_id)
  named$cell_id[nrow(named)] <- "new-7"
  chr <- .va_at(L, L$held, seeds = named)
  expect_type(chr$cell_id, "character")
  expect_type(chr$filled_from, "character")
  expect_identical(chr$cell_id, ifelse(mine, "new-7", as.character(at$cell_id)))
  expect_identical(chr$filled_from, as.character(at$filled_from))
  expect_identical(.va_cols(chr[-c(1L, 7L)]), .va_cols(at[-c(1L, 7L)]))
})


# ---- What passes through; other CRSs; empty geometry ------------------------------

test_that("`estimate` and the interval columns follow the summary", {
  L <- .va_full()
  pooled <- .va_at(L, L$held)
  expect_identical(names(pooled), c("cell_id", "value", "estimate", "obs_lo", "obs_hi", "filled",
                                    "filled_from", "dist_nearest", "dist_ratio"))
  expect_identical(unique(pooled$estimate), "pooled")
  expect_identical(attr(pooled, "values_at")$obs_level, 0.90)
  expect_true(all(pooled$obs_lo < pooled$value & pooled$value < pooled$obs_hi))
  plain <- .va_at(L, L$held, values = L$plain)
  expect_identical(names(plain), c("cell_id", "value", "estimate", "filled", "filled_from",
                                   "dist_nearest", "dist_ratio"))
  expect_identical(unique(plain$estimate), "mean")
  expect_null(attr(plain, "values_at")$obs_level)
  expect_true("obs_level" %in% names(attr(plain, "values_at")))
  expect_identical(plain[c("cell_id", "filled", "filled_from", "dist_nearest", "dist_ratio")],
                   pooled[c("cell_id", "filled", "filled_from", "dist_nearest", "dist_ratio")])
  # Pooled values without an interval: no interval columns, and the pooled value is read.
  none <- L$pooled[setdiff(names(L$pooled), c("..obs_lo_resp_z", "..obs_hi_resp_z"))]
  at <- .va_at(L, L$held, values = none)
  expect_identical(names(at), names(plain))
  expect_identical(at$value, pooled$value)
  expect_identical(unique(at$estimate), "pooled")
  expect_null(attr(at, "values_at")$obs_level)
  # A frame that kept the columns and lost its record: the level is not known.
  bare <- L$pooled
  attr(bare, "pooled") <- NULL
  expect_identical(attr(.va_at(L, L$held, values = bare), "values_at")$obs_level, NA_real_)
  expect_identical(.va_at(L, L$held, values = bare)$obs_lo, pooled$obs_lo)
  # The record names the estimate the column holds, with or without the summary's record.
  expect_identical(attr(.va_at(L, L$held, values = bare), "values_at")$estimate, "pooled")
  expect_identical(unique(.va_at(L, L$held, values = bare)$estimate), "pooled")
  expect_identical(attr(plain, "values_at")$estimate, "mean")
})

test_that("a column summarised as response and as predictor is read as the response; a predictor by its own columns", {
  L <- .va_full()
  both <- suppressWarnings(
    summarize_by_cell(L$asg, "z", predictor_vars = c("z", "w"), id_col = "cell_id",
                      cells_sf = L$tess$cells, estimate = "pooled", obs_level = 0.90))
  expect_true(all(c("..pooled_resp_z", "..pooled_pred_z", "..pooled_pred_w", "resp_mean_z",
                    "pred_mean_z", "pred_mean_w", "..obs_lo_pred_w") %in% names(both)))
  want <- .va_at(L, L$held)
  # The predictor's copies moved, so that the two can be told apart.
  moved <- both
  for (nm in c("..pooled_pred_z", "pred_mean_z", "..obs_lo_pred_z", "..obs_hi_pred_z"))
    moved[[nm]] <- moved[[nm]] + 1000
  at <- .va_at(L, L$held, values = moved)
  expect_identical(at$value, want$value)
  expect_identical(at$obs_lo, want$obs_lo)
  # The same for plain means.
  plain <- sf::st_drop_geometry(moved)[c("cell_id", "resp_mean_z", "pred_mean_z", "pred_mean_w")]
  expect_identical(.va_at(L, L$held, values = plain)$value,
                   .va_at(L, L$held, values = L$plain)$value)
  # A predictor: its pooled value and interval, then its mean.
  w <- .va_at(L, L$held, values = both, response = "w")
  src <- match(w$cell_id, both$cell_id)
  expect_identical(w$value, both$..pooled_pred_w[src])
  expect_identical(w$obs_hi, both$..obs_hi_pred_w[src])
  expect_identical(attr(w, "values_at")$column, "w")
  expect_identical(.va_at(L, L$held, values = plain, response = "w")$value,
                   both$pred_mean_w[src])
})

test_that("`newdata` in lon/lat or another CRS is placed where the rows are", {
  L <- .va_gap()
  for (seeds in list(L$seeds, NULL)) {
    want <- .va_at(L, L$held, seeds = seeds)
    for (crs in c(4326, 3857, 32616)) {
      at <- .va_at(L, sf::st_transform(L$held, crs), seeds = seeds)
      expect_identical(at[1:7], want[1:7])
      expect_equal(at$dist_nearest, want$dist_nearest, tolerance = 1e-8)
      expect_identical(attr(at, "values_at"), attr(want, "values_at"))
    }
  }
  # Seeds in lon/lat.
  expect_identical(.va_at(L, L$held, seeds = sf::st_transform(L$seeds, 4326))[1:7],
                   .va_at(L, L$held)[1:7])
})

test_that("a row with an empty geometry is NA in every column, and a far coordinate is no location", {
  L <- .va_full()
  # Row 3 is empty; rows 5, 7, 8 and 9 have a coordinate too large to take a distance
  # from (x, then y alone) or one that is not finite.
  g <- sf::st_geometry(L$held[1:9, ])
  g[[3L]] <- sf::st_point()
  g[[5L]] <- sf::st_point(c(1e200, 5e6))
  g[[7L]] <- sf::st_point(c(5e5, -1e200))
  g[[8L]] <- sf::st_point(c(Inf, 5e6))
  g[[9L]] <- sf::st_point(c(5e5, NaN))
  newdata <- sf::st_sf(k = 1:9, geometry = g)
  for (seeds in list(L$seeds, NULL)) {
    r <- .va_catch(.va_at(L, newdata, seeds = seeds))
    expect_identical(.va_classes(r), character(0))
    at <- r$value
    expect_identical(nrow(at), 9L)
    for (i in c(3L, 5L, 7L, 8L, 9L)) expect_true(all(is.na(unlist(at[i, ]))))
    expect_type(at$estimate, "character")
    expect_type(at$filled, "logical")
    keep <- c(1L, 2L, 4L, 6L)
    want <- .va_at(L, L$held[keep, ], seeds = seeds)
    expect_identical(.va_cols(at[keep, ]), .va_cols(want))
    expect_identical(attr(at, "values_at")$filled, sum(at$filled, na.rm = TRUE))
  }
  # No location at all: no row, and the record still describes the scored rows.
  none <- .va_at(L, L$held[0, ])
  expect_identical(dim(none), c(0L, 9L))
  expect_identical(attr(none, "values_at")$n, 480L)
  expect_identical(attr(none, "values_at")$filled, 0L)
  expect_identical(dim(.va_at(L, L$held[0, ], seeds = NULL)), c(0L, 9L))
})

test_that("a feature that is not a point is read at the point coerce_to_points() gives it", {
  L <- .va_full()
  squares <- sf::st_buffer(L$held[1:30, ], 400, endCapStyle = "SQUARE")
  expect_true(all(sf::st_geometry_type(squares) == "POLYGON"))
  for (seeds in list(L$seeds, NULL))
    expect_identical(.va_at(L, squares, seeds = seeds),
                     .va_at(L, coerce_to_points(squares, "auto"), seeds = seeds))
  # A square is read at its centre, which is the held-out location.
  expect_identical(.va_at(L, squares)$cell_id, .va_at(L, L$held[1:30, ])$cell_id)
  # No column of `newdata` is read: one called like a column of the result, or the ID.
  odd <- L$held[1:30, ]
  odd$cell_id <- "x"; odd$value <- -1; odd$filled <- NA
  expect_identical(.va_at(L, odd, seeds = NULL), .va_at(L, L$held[1:30, ], seeds = NULL))
})

# Two cells made by hand: cell 1 is two squares of side 10, at x from 0 and from 20, and
# cell 2 the square between them; two rows in each square.
.va_parts <- function() .va_keep("parts", {
  square <- function(x0, y0 = 0, side = 10)
    sf::st_polygon(list(cbind(5e5 + x0 + c(0, side, side, 0, 0), 5e6 + y0 + c(0, 0, side, side, 0))))
  cells <- sf::st_sf(cell_id = 1:2,
                     geometry = sf::st_sfc(sf::st_multipolygon(list(square(0), square(20))),
                                           square(10), crs = 32617))
  rows <- .va_points(c(2, 8, 22, 28, 12, 18), rep(5, 6L), cell_id = c(1L, 1L, 1L, 1L, 2L, 2L),
                     z = c(1, 2, 3, 2, 10, 12))
  list(square = square, cells = cells, rows = rows,
       values = summarize_by_cell(rows, "z", id_col = "cell_id", cells_sf = cells,
                                  estimate = "mean"))
})

test_that("a cell of several parts holds the locations in each of them", {
  T <- .va_parts()
  expect_identical(as.character(sf::st_geometry_type(T$values)), c("MULTIPOLYGON", "POLYGON"))
  r <- .va_catch(cell_values_at(.va_points(c(5, 25, 15, 40), c(3, 7, 5, 5)), T$values, T$rows, "z",
                                id_col = "cell_id"))
  expect_identical(.va_classes(r), character(0))
  at <- r$value
  expect_identical(at$cell_id, c(1L, 1L, 2L, NA))
  expect_identical(at$value, c(2, 2, 11, 2))
  expect_identical(at$filled, c(FALSE, FALSE, FALSE, TRUE))
  expect_identical(at$filled_from, c(NA, NA, NA, 1L))
  expect_identical(at$dist_nearest, c(sqrt(13), sqrt(13), 3, 12))
})

test_that("a layer that mixes points with other features is read feature by feature", {
  T <- .va_parts()
  # A point in the first square, a square of side 4 inside the second cell, and a point
  # without coordinates.
  mixed <- sf::st_sf(k = 1:3, geometry = sf::st_sfc(sf::st_point(c(5e5 + 5, 5e6 + 3)),
                                                    T$square(13, 2, 4), sf::st_point(),
                                                    crs = 32617))
  expect_identical(as.character(sf::st_geometry_type(mixed)), c("POINT", "POLYGON", "POINT"))
  seeds <- .va_points(c(5, 15, 25), c(5, 5, 5), cell_id = c(1L, 2L, 3L))
  for (s in list(NULL, seeds)) {
    r <- .va_catch(cell_values_at(mixed, T$values, T$rows, "z", id_col = "cell_id", seeds = s))
    expect_identical(.va_classes(r), character(0))
    at <- r$value
    expect_identical(nrow(at), 3L)
    expect_identical(at$cell_id, c(1L, 2L, NA))
    expect_identical(at$value, c(2, 11, NA))
    expect_identical(at$filled, c(FALSE, FALSE, NA))
    # The square is read at its centre, (15, 4).
    expect_identical(at$dist_nearest, c(sqrt(13), sqrt(10), NA))
  }
})


# ---- By the polygons: the CRS handed over, and what the assignment says ---------------

# The WARN lines the package logs while `expr` runs, in order, and nothing on the console.
.va_logged <- function(env = parent.frame()) {
  log <- new.env(parent = emptyenv())
  log$lines <- character(0)
  testthat::local_mocked_bindings(.sk_log = function(level, msg, raising = FALSE) {
    if (identical(attr(level, "level"), "WARN")) log$lines <- c(log$lines, msg)
    invisible(msg)
  }, .env = env)
  log
}
# The assignment, with what it is handed kept in `seen`.
.va_watch_assignment <- function(env = parent.frame()) {
  seen <- new.env(parent = emptyenv())
  real <- assign_features_to_polygons
  force(real)
  testthat::local_mocked_bindings(assign_features_to_polygons = function(features_sf, polygons_sf, ...) {
    seen$features <- sf::st_crs(features_sf)
    seen$polygons <- sf::st_crs(polygons_sf)
    seen$names    <- names(polygons_sf)
    seen$class    <- class(polygons_sf)
    seen$n        <- nrow(polygons_sf)
    real(features_sf, polygons_sf, ...)
  }, .env = env)
  seen
}

test_that("without a CRS on `newdata` and `points` the locations are placed where the summary measured, not where the polygons' CRS says", {
  L <- .va_full()
  bare <- function(x) sf::st_set_crs(x, NA)
  seen <- .va_watch_assignment()
  want <- .va_at(L, L$held, seeds = NULL)
  expect_identical(seen$features, sf::st_crs(32617))
  # Beside polygons with a CRS, locations with one go to the assignment as they arrived:
  # it transforms them once.
  in_3857 <- .va_at(L, sf::st_transform(L$held, 3857), seeds = NULL)
  expect_identical(seen$features, sf::st_crs(3857))
  expect_identical(seen$polygons, sf::st_crs(32617))
  expect_identical(in_3857[1:7], want[1:7])
  # The working CRS is the record's, and the polygons are in it: no warning of the
  # assignment about a layer the caller did not hand it.
  r <- .va_catch(.va_at(L, bare(L$held), points = bare(L$asg), seeds = NULL))
  expect_identical(.va_classes(r), character(0))
  expect_identical(r$value, want)
  expect_identical(seen$features, sf::st_crs(32617))
  # The polygons in another CRS than the working one: a summary transformed for a map,
  # which keeps its record.
  on_map <- sf::st_transform(L$pooled, 4326)
  expect_identical(attr(on_map, "count_rule")$crs, sf::st_crs(32617))
  want_map <- .va_at(L, L$held, values = on_map, seeds = NULL)
  expect_false(anyNA(want_map$cell_id))
  r <- .va_catch(.va_at(L, bare(L$held), values = on_map, points = bare(L$asg), seeds = NULL))
  expect_identical(.va_classes(r), character(0))
  expect_identical(r$value, want_map)
  expect_identical(seen$features, sf::st_crs(32617))
  expect_identical(seen$polygons, sf::st_crs(4326))
})

test_that("cells without a CRS are read as they were for the rows, whatever CRS `newdata` arrives in", {
  L <- .va_full()
  # 36 square cells over the layer, without a CRS, beside rows with one.
  grid <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(rbind(L$rows["z"], L$held["z"]))), n = c(6, 6))
  units <- sf::st_sf(unit = sprintf("U%02d", seq_along(grid)), geometry = sf::st_set_crs(grid, NA))
  asg <- suppressWarnings(assign_features_to_polygons(L$rows, units, polygon_id_col = "unit",
                                                      quiet = TRUE))
  out <- suppressWarnings(summarize_by_cell(asg, "z", id_col = "unit", cells_sf = units,
                                            estimate = "mean"))
  expect_true(is.na(sf::st_crs(out)))
  expect_identical(sf::st_crs(asg), sf::st_crs(32617))
  seen <- .va_watch_assignment()
  call <- function(newdata) .va_catch(cell_values_at(newdata, out, asg, "z", id_col = "unit"))
  base <- call(L$held)
  expect_false(anyNA(base$value$unit))
  expect_false(any(base$value$filled))
  # The cells of the rows' own coordinates.
  inside <- sf::st_intersects(sf::st_set_crs(sf::st_geometry(L$held), NA), sf::st_geometry(units))
  expect_identical(base$value$unit, units$unit[vapply(inside, function(i) i[1L], integer(1))])
  for (crs in c(4326, 3857, 32616)) {
    r <- call(sf::st_transform(L$held, crs))
    expect_identical(r$value[1:5], base$value[1:5])
    expect_equal(r$value$dist_nearest, base$value$dist_nearest, tolerance = 1e-8)
    # Handed over in the CRS of `points`; the assignment stamps the cells with it and says so.
    expect_identical(seen$features, sf::st_crs(32617))
    expect_true(is.na(seen$polygons))
    expect_length(r$warnings, 1L)
    expect_match(conditionMessage(r$warnings[[1L]]), "harmonize_crs(): `b` has no CRS", fixed = TRUE)
  }
  # Rows in lon/lat beside cells in degrees without a CRS: the CRS of `points`, not the
  # projected one the distances are measured in.
  rows_ll <- sf::st_transform(L$rows, 4326)
  grid_ll <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(sf::st_transform(rbind(L$rows["z"], L$held["z"]),
                                                                        4326))), n = c(6, 6))
  units_ll <- sf::st_sf(unit = sprintf("U%02d", seq_along(grid_ll)),
                        geometry = sf::st_set_crs(grid_ll, NA))
  asg_ll <- suppressWarnings(assign_features_to_polygons(rows_ll, units_ll, polygon_id_col = "unit",
                                                         keep_unassigned = TRUE, quiet = TRUE))
  out_ll <- suppressWarnings(summarize_by_cell(asg_ll, "z", id_col = "unit", cells_sf = units_ll,
                                               estimate = "mean"))
  in_ll  <- suppressWarnings(cell_values_at(sf::st_transform(L$held, 4326), out_ll, asg_ll, "z",
                                            id_col = "unit"))
  expect_identical(seen$features, sf::st_crs(4326))
  in_utm <- suppressWarnings(cell_values_at(L$held, out_ll, asg_ll, "z", id_col = "unit"))
  expect_identical(seen$features, sf::st_crs(4326))
  expect_false(isTRUE(sf::st_is_longlat(attr(in_utm, "values_at")$crs)))
  expect_identical(in_utm[1:5], in_ll[1:5])
  # Nothing with a CRS, the cells neither: the coordinates of all three as they are, and
  # no CRS for the assignment to stamp on anything.  The one warning is
  # ensure_projected()'s, about the rows it took for lon/lat.
  bare <- function(x) sf::st_set_crs(x, NA)
  r <- .va_catch(cell_values_at(bare(sf::st_transform(L$held, 4326)), out_ll, bare(asg_ll), "z",
                                id_col = "unit"))
  expect_true(is.na(seen$features))
  expect_true(is.na(seen$polygons))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]),
               "ensure_projected(): CRS is missing; assuming EPSG:4326", fixed = TRUE)
  expect_identical(r$value$unit, units_ll$unit[vapply(
    sf::st_intersects(bare(sf::st_geometry(sf::st_transform(L$held, 4326))),
                      sf::st_geometry(units_ll)), function(i) i[1L], integer(1))])
  # The same in planar coordinates: no warning at all.
  r <- .va_catch(cell_values_at(bare(L$held), out, bare(asg), "z", id_col = "unit"))
  expect_identical(.va_classes(r), character(0))
  expect_true(is.na(seen$features))
  expect_identical(r$value[1:5], base$value[1:5])
})

test_that("the assignment is handed the cells' polygons and IDs alone: a grouped summary, a row without an ID", {
  L <- .va_full()
  want <- .va_at(L, L$held, seeds = NULL)
  seen <- .va_watch_assignment()
  grouped <- dplyr::group_by(L$pooled, cell_id)
  expect_s3_class(grouped, "grouped_df")
  at <- .va_at(L, L$held, values = grouped, seeds = NULL)
  expect_identical(.va_cols(at[1:7]), .va_cols(want[1:7]))
  expect_identical(seen$names, c("cell_id", "geometry"))
  expect_identical(seen$class, c("sf", "data.frame"))
  expect_identical(seen$n, nrow(L$pooled))
  # A row of the summary without an ID is no cell: its polygon holds no location, and the
  # locations in it are filled from the cell of their nearest scored row, whose value the
  # summary no longer has under that ID.
  gone <- as.integer(names(which.max(table(want$cell_id))))
  values <- L$pooled
  values$cell_id[values$cell_id == gone] <- NA
  at <- .va_at(L, L$held, values = values, seeds = NULL)
  hit <- want$cell_id == gone
  expect_identical(seen$n, nrow(L$pooled) - 1L)
  expect_true(all(is.na(at$cell_id[hit])))
  expect_identical(at$filled, hit)
  expect_identical(.va_cols(at[!hit, ]), .va_cols(want[!hit, ]))
  # No ID at all: no location has a cell, and the assignment is not called.
  seen$n <- NULL
  values$cell_id <- NA_integer_
  at <- .va_at(L, L$held, values = values, seeds = NULL)
  expect_null(seen$n)
  expect_identical(at$filled, rep(TRUE, nrow(at)))
  expect_true(all(is.na(at$cell_id)))
})

test_that("the warning that no location is in a polygon is neither raised nor logged, and every other warning of the assignment is", {
  L <- .va_full()
  bare <- function(x) sf::st_set_crs(x, NA)
  far <- .va_points(c(-2e4, -3e4, 3e5), c(5e4, -3e4, 3e5))
  log <- .va_logged()
  direct <- function()
    .va_catch(assign_features_to_polygons(far, L$pooled, polygon_id_col = "cell_id",
                                          keep_unassigned = TRUE, quiet = TRUE))
  # The assignment raises the warning for these locations, and logs its line.
  expect_identical(.va_classes(direct()), "spatialkit_unassigned")
  expect_length(log$lines, 1L)
  expect_match(log$lines, "none of the 3 feature(s) fall inside any polygon", fixed = TRUE)
  # Inside the call: neither.
  log$lines <- character(0)
  r <- .va_catch(.va_at(L, far, seeds = NULL))
  expect_identical(.va_classes(r), character(0))
  expect_identical(log$lines, character(0))
  expect_identical(r$value$filled, rep(TRUE, 3L))
  # Once the call has returned the warning is the assignment's again.
  expect_null(.answered_state$classes)
  expect_identical(.va_classes(direct()), "spatialkit_unassigned")
  expect_length(log$lines, 1L)
  # Also when the assignment stops: two polygons with one ID.
  twice <- L$pooled
  twice$cell_id[2L] <- twice$cell_id[1L]
  expect_match(conditionMessage(.va_error(.va_at(L, far, values = twice, seeds = NULL))),
               "assign_features_to_polygons(): the polygon ID column 'cell_id'", fixed = TRUE)
  expect_null(.answered_state$classes)
  expect_identical(.va_classes(direct()), "spatialkit_unassigned")
  # Another warning of the assignment reaches the caller, with its line.  Nothing here
  # names a CRS but the polygons, so the assignment stamps theirs on the locations.
  values <- L$pooled
  rule <- attr(values, "count_rule"); rule$crs <- sf::NA_crs_
  attr(values, "count_rule") <- rule
  pooled <- attr(values, "pooled"); pooled$crs <- sf::NA_crs_
  attr(values, "pooled") <- pooled
  log$lines <- character(0)
  r <- .va_catch(.va_at(L, bare(L$held), values = values, points = bare(L$asg), seeds = NULL))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), "harmonize_crs(): `a` has no CRS", fixed = TRUE)
  expect_identical(log$lines, conditionMessage(r$warnings[[1L]]))
  expect_true(is.na(attr(r$value, "values_at")$crs))
  expect_identical(.va_cols(r$value), .va_cols(.va_at(L, L$held, seeds = NULL)))
  # The same with every location outside: that warning passes, the other does not.
  log$lines <- character(0)
  r <- .va_catch(.va_at(L, bare(far), values = values, points = bare(L$asg), seeds = NULL))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), "harmonize_crs(): `a` has no CRS", fixed = TRUE)
  expect_identical(log$lines, conditionMessage(r$warnings[[1L]]))
  expect_identical(r$value$filled, rep(TRUE, 3L))
})

test_that("a span that answers for a class keeps that class back and no other, and closes to what was open before", {
  log <- .va_logged()
  raise <- function(class) .va_classes(.va_catch(.warn_classed(class, "%s is said", class[1L])))
  expect_null(.answered_state$classes)
  expect_identical(raise("spatialkit_one"), "spatialkit_one")
  outer <- .answer_for_classed("spatialkit_one")
  expect_identical(raise("spatialkit_one"), character(0))
  expect_identical(raise("spatialkit_two"), "spatialkit_two")
  # A condition of several classes is kept back when any of them is answered for.
  expect_identical(raise(c("spatialkit_three", "spatialkit_one")), character(0))
  # A span opened inside another adds its class to the outer one's.
  inner <- .answer_for_classed("spatialkit_two")
  expect_identical(raise("spatialkit_two"), character(0))
  expect_identical(raise("spatialkit_one"), character(0))
  inner()
  expect_identical(.answered_state$classes, "spatialkit_one")
  expect_identical(raise("spatialkit_two"), "spatialkit_two")
  expect_identical(raise("spatialkit_one"), character(0))
  outer()
  expect_null(.answered_state$classes)
  expect_identical(raise("spatialkit_one"), "spatialkit_one")
  expect_identical(log$lines, c("spatialkit_one is said", "spatialkit_two is said",
                                "spatialkit_two is said", "spatialkit_one is said"))
  # A warning kept back carries no fields away and returns its text unseen.
  done <- .answer_for_classed("spatialkit_one")
  expect_invisible(.warn_classed("spatialkit_one", "%d", 3L, data = list(n = 3L)))
  expect_identical(.warn_classed("spatialkit_one", "%d", 3L), "3")
  done()
})


# ---- The record of the call ------------------------------------------------------

test_that("attr(, \"values_at\") says how the cells were found and what the distances are measured against", {
  L <- .va_gap()
  at <- .va_at(L, L$held)
  rec <- attr(at, "values_at")
  expect_s3_class(rec, "spatialkit_values_at")
  expect_identical(names(rec), c("assignment", "crs", "column", "estimate", "obs_level", "n",
                                 "locations", "spacing", "filled", "count_rule"))
  expect_identical(rec$assignment, "seeds")
  expect_identical(attr(.va_at(L, L$held, seeds = NULL), "values_at")$assignment, "polygons")
  expect_identical(rec$crs, attr(L$pooled, "count_rule")$crs)
  expect_identical(rec$column, "z")
  expect_identical(rec$estimate, "pooled")
  expect_identical(rec$obs_level, 0.90)
  scored <- is.finite(L$asg$z)
  expect_identical(rec$n, sum(scored))
  expect_identical(rec$locations, sum(scored))
  # The spacing: the median, over the scored locations, of the distance to the nearest other.
  d <- as.matrix(stats::dist(.va_xy(L$asg)[scored, ]))
  diag(d) <- Inf
  expect_equal(rec$spacing, stats::median(apply(d, 1L, min)), tolerance = 1e-12)
  expect_identical(at$dist_ratio, at$dist_nearest / rec$spacing)
  expect_identical(rec$filled, sum(at$filled))
  expect_gt(rec$filled, 0L)
  expect_identical(rec$count_rule, attr(L$pooled, "count_rule"))
  expect_identical(rec$count_rule$check, list(estimate = "match", cells = "match"))
  # Plain means of cells sized for pooled values carry the record with its mismatch.
  expect_identical(attr(.va_at(L, L$held, values = L$plain), "values_at")$count_rule,
                   attr(L$plain, "count_rule"))
  # Without a record: NULL, and the CRS the pooled columns were measured in.
  bare <- L$pooled
  attr(bare, "count_rule") <- NULL
  no_rule <- attr(.va_at(L, L$held, values = bare), "values_at")
  expect_null(no_rule$count_rule)
  expect_true("count_rule" %in% names(no_rule))
  expect_identical(no_rule$crs, attr(L$pooled, "pooled")$crs)
  # The distances are those from the scored rows: a row without a response is not one.
  unscored <- L$asg[!scored, ]
  expect_true(all(.va_at(L, unscored)$dist_nearest > 0))
})

test_that("the working CRS is the record's, then the pooled columns', then the rows' own", {
  L <- .va_full()
  want <- .va_at(L, L$held)
  # The rows and the new locations in another projected CRS than the count was scored in.
  rows <- sf::st_transform(L$asg, 3857)
  held <- sf::st_transform(L$held, 3857)
  at <- cell_values_at(held, L$pooled, rows, "z", id_col = "cell_id", seeds = L$seeds)
  expect_identical(attr(at, "values_at")$crs, sf::st_crs(32617))
  expect_identical(at[1:7], want[1:7])
  expect_equal(at$dist_nearest, want$dist_nearest, tolerance = 1e-8)
  # The record's CRS comes before the pooled columns'.
  values <- L$pooled
  other <- attr(values, "pooled"); other$crs <- sf::st_crs(3857)
  attr(values, "pooled") <- other
  at <- cell_values_at(held, values, rows, "z", id_col = "cell_id", seeds = L$seeds)
  expect_identical(attr(at, "values_at")$crs, sf::st_crs(32617))
  # The record gone: the pooled columns' CRS.
  values <- L$pooled
  attr(values, "count_rule") <- NULL
  at <- cell_values_at(held, values, rows, "z", id_col = "cell_id", seeds = L$seeds)
  expect_identical(attr(at, "values_at")$crs, sf::st_crs(32617))
  expect_equal(at$dist_nearest, want$dist_nearest, tolerance = 1e-8)
  # Both gone: the rows' own CRS, in which a distance is another number (Web Mercator
  # stretches these latitudes by about 1.4).
  attr(values, "pooled") <- NULL
  at <- cell_values_at(held, values, rows, "z", id_col = "cell_id", seeds = L$seeds)
  expect_identical(attr(at, "values_at")$crs, sf::st_crs(3857))
  expect_gt(stats::median(at$dist_nearest / want$dist_nearest), 1.3)
  # A record whose CRS is geographic or missing is not a CRS to measure in.
  for (crs in list(sf::st_crs(4326), sf::NA_crs_)) {
    values <- L$pooled
    rule <- attr(values, "count_rule"); rule$crs <- crs
    attr(values, "count_rule") <- rule
    at <- cell_values_at(held, values, rows, "z", id_col = "cell_id", seeds = L$seeds)
    expect_identical(attr(at, "values_at")$crs, sf::st_crs(32617))
  }
})

test_that("the record prints what was read and what the distances mean", {
  L <- .va_gap()
  at <- .va_at(L, L$held)
  rec <- attr(at, "values_at")
  expect_identical(.va_print(rec), c(
    "Cell values at other locations: z, the pooled value",
    "  cell        : of the nearest seed (the Voronoi cell before its clip to the cover)",
    "  interval    : 90% for ONE OBSERVATION (obs_lo, obs_hi), centred on the value",
    sprintf("  scored rows : %d rows at %d locations, measured in EPSG:32617", rec$n, rec$locations),
    sprintf("  spacing     : %s, the median distance from a scored location to the",
            format(signif(rec$spacing, 4))),
    "                nearest other; dist_ratio is dist_nearest over it",
    sprintf("  filled      : %d locations have no cell with a value and take that of", rec$filled),
    "                the cell of the nearest scored row (filled_from)",
    sprintf("  count       : %d cells, set by size (the summary's record, in `count_rule`)",
            rec$count_rule$cells),
    "Away from the observations no interval here is calibrated."))
  utils::capture.output(returned <- print(rec))
  expect_identical(returned, rec)
  # Plain means by the polygons, nothing filled, no record of the count.
  F <- .va_full()
  values <- F$plain
  attr(values, "count_rule") <- NULL
  rec <- attr(.va_at(F, F$held[1:10, ], values = values, seeds = NULL), "values_at")
  expect_identical(.va_print(rec), c(
    "Cell values at other locations: z, the plain mean",
    "  cell        : the polygon of `values` that holds the location",
    "  scored rows : 480 rows at 480 locations, measured in EPSG:32617",
    sprintf("  spacing     : %s, the median distance from a scored location to the",
            format(signif(rec$spacing, 4))),
    "                nearest other; dist_ratio is dist_nearest over it",
    "  filled      : no location"))
  # One filled location, one scored row, no CRS, no spacing.
  one <- structure(list(assignment = "seeds", crs = sf::NA_crs_, column = "v", estimate = "mean",
                        obs_level = NULL, n = 1L, locations = 1L, spacing = NA_real_,
                        filled = 1L, count_rule = NULL), class = "spatialkit_values_at")
  expect_identical(.va_print(one), c(
    "Cell values at other locations: v, the plain mean",
    "  cell        : of the nearest seed (the Voronoi cell before its clip to the cover)",
    "  scored rows : 1 row at 1 location, coordinates without a CRS",
    "  spacing     : not available: the scored rows sit on one location, so dist_ratio is NA",
    "  filled      : 1 location has no cell with a value and takes that of",
    "                the cell of the nearest scored row (filled_from)"))
  # An interval whose level the frame did not keep.
  one$obs_level <- NA_real_
  expect_true("  interval    : for ONE OBSERVATION (obs_lo, obs_hi), centred on the value" %in%
                .va_print(one))
})

test_that("the print method of the record is registered", {
  registered <- getNamespaceInfo(asNamespace("spatialkit"), "S3methods")
  expect_true(any(registered[, 1L] == "print" & registered[, 2L] == "spatialkit_values_at"))
  expect_true("cell_values_at" %in% getNamespaceExports("spatialkit"))
  expect_identical(names(formals(cell_values_at)),
                   c("newdata", "values", "points", "response_var", "id_col", "seeds"))
  expect_identical(formals(cell_values_at)$id_col, "poly_id")
  expect_null(formals(cell_values_at)$seeds)
})


# ---- Edge cases ------------------------------------------------------------------

test_that("a seed with an NA ID has no cell: the locations nearest it are filled", {
  L <- .va_full()
  want <- .va_at(L, L$held)
  expect_false(any(want$filled))
  # The seed of the cell that holds the most held-out rows.
  gone <- as.integer(names(which.max(table(want$cell_id))))
  seeds <- L$seeds
  seeds$cell_id[seeds$cell_id == gone] <- NA
  at <- .va_at(L, L$held, seeds = seeds)
  hit <- want$cell_id == gone
  expect_gt(sum(hit), 1L)
  expect_true(all(is.na(at$cell_id[hit])))
  expect_identical(at$filled, hit)
  # The cell they are filled from is the one that holds their nearest scored row.
  u <- .va_xy(L$asg)
  nearest <- .va_nearest(.va_xy(L$held)[hit, , drop = FALSE], u)
  expect_identical(at$filled_from[hit], L$asg$cell_id[nearest])
  expect_identical(at$value[hit],
                   L$pooled$..pooled_resp_z[match(at$filled_from[hit], L$pooled$cell_id)])
  expect_identical(.va_cols(at[!hit, ]), .va_cols(want[!hit, ]))
  # A row of the summary without an ID (the rows no cell holds) is no cell either: a
  # location without a cell does not take its value.
  flat <- sf::st_drop_geometry(L$pooled)
  stray <- flat[1L, ]
  stray$cell_id <- NA
  stray$..pooled_resp_z <- 999
  expect_identical(.va_cols(.va_at(L, L$held, values = rbind(flat, stray), seeds = seeds)),
                   .va_cols(at))
  # Several seeds without a cell are no repeated ID.
  seeds$cell_id[1:3] <- NA
  expect_identical(nrow(.va_at(L, L$held, seeds = seeds)), nrow(L$held))
  # A seed with an empty geometry is nearest to nothing.
  moved <- L$seeds
  g <- sf::st_geometry(moved)
  g[[which(moved$cell_id == gone)]] <- sf::st_point()
  moved <- sf::st_set_geometry(moved, g)
  at <- .va_at(L, L$held, seeds = moved)
  expect_false(any(at$cell_id == gone))
  expect_identical(.va_cols(at[!hit, ]), .va_cols(want[!hit, ]))
})

test_that("what `seeds` and `values` must be", {
  L <- .va_full()
  seeds_text <- paste0(
    "cell_values_at(): `seeds` must be an sf POINT layer with a column 'cell_id' that gives ",
    "each seed's cell; for cells from build_tessellation(seeds), set seeds$cell_id <- tess$index.")
  # Two seeds with one ID: both have the cell of another, and both are counted.
  twice <- L$seeds
  twice$cell_id[2L] <- twice$cell_id[1L]
  expect_identical(conditionMessage(.va_error(.va_at(L, L$held, seeds = twice))), paste0(
    "cell_values_at(): 2 seeds have the cell of another in `seeds$cell_id`; each seed stands ",
    "for one cell. Set the column from the tessellation built on these seeds ",
    "(seeds$cell_id <- tess$index)."))
  twice$cell_id[3:4] <- twice$cell_id[1L]
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, seeds = twice))),
               "cell_values_at(): 4 seeds have the cell of another", fixed = TRUE)
  # Two pairs.
  pairs <- L$seeds
  pairs$cell_id[c(2L, 9L)] <- pairs$cell_id[c(1L, 8L)]
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, seeds = pairs))),
               "cell_values_at(): 4 seeds have the cell of another", fixed = TRUE)
  # A plain error, raised without a call.
  e <- .va_error(.va_at(L, L$held, seeds = twice))
  expect_null(e$call)
  expect_identical(class(e), c("simpleError", "error", "condition"))
  # Seeds on one location share its cell and are no error: the same frame as without
  # the repeats, wherever in `seeds` they stand.
  own <- L$seeds[c("seed_id", "method", "cell_id")]
  expect_identical(.va_at(L, L$held, seeds = rbind(own, own[c(5L, 9L), ])), .va_at(L, L$held))
  expect_identical(.va_at(L, L$held, seeds = rbind(own[c(9L, 9L), ], own)), .va_at(L, L$held))
  # A third seed elsewhere with their ID: all three are counted.
  third <- rbind(own, own[5L, ])
  third$cell_id[1L] <- third$cell_id[5L]
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, seeds = third))),
               "cell_values_at(): 3 seeds have the cell of another", fixed = TRUE)
  # Seeds that differ in one coordinate alone are at two locations.
  for (shift in list(c(1, 0), c(0, 1))) {
    beside <- rbind(own, own[5L, ])
    g <- sf::st_geometry(beside)
    g[[nrow(beside)]] <- sf::st_point(sf::st_coordinates(g[5L])[1L, 1:2] + shift)
    expect_match(conditionMessage(.va_error(.va_at(L, L$held, seeds = sf::st_set_geometry(beside, g)))),
                 "cell_values_at(): 2 seeds have the cell of another", fixed = TRUE)
  }
  # A seed without a location and the ID of another is not a seed on its location.
  nowhere <- rbind(own, own[5L, ])
  g <- sf::st_geometry(nowhere)
  g[[nrow(nowhere)]] <- sf::st_point()
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, seeds = sf::st_set_geometry(nowhere, g)))),
               "cell_values_at(): 2 seeds have the cell of another", fixed = TRUE)
  # Seeds without the ID column; not points; not an sf layer.
  bare <- L$seeds[c("seed_id", "method")]
  expect_identical(conditionMessage(.va_error(.va_at(L, L$held, seeds = bare))), seeds_text)
  expect_identical(conditionMessage(.va_error(.va_at(L, L$held, seeds = L$tess$cells))),
                   seeds_text)
  expect_identical(
    conditionMessage(.va_error(.va_at(L, L$held, seeds = sf::st_drop_geometry(L$seeds)))),
    seeds_text)
  expect_identical(conditionMessage(.va_error(.va_at(L, L$held, seeds = 107L))), seeds_text)
  # A summary without the cells' polygons, and no seeds.
  flat <- sf::st_drop_geometry(L$pooled)
  polygons_text <- paste0(
    "cell_values_at(): `values` carries no cell polygons, so a location cannot be placed in a ",
    "cell. Build `values` with `cells_sf`, or pass the seeds the cells were built on as `seeds` ",
    "(seeds$cell_id <- tess$index).")
  expect_identical(conditionMessage(.va_error(.va_at(L, L$held, values = flat, seeds = NULL))),
                   polygons_text)
  # A summary of no rows has no polygon either; with seeds it has no value for any cell.
  expect_identical(
    conditionMessage(.va_error(.va_at(L, L$held, values = L$pooled[0, ], seeds = NULL))),
    polygons_text)
  none <- .va_at(L, L$held[1:5, ], values = L$pooled[0, ])
  expect_true(all(is.na(none$value)))
  expect_identical(none$cell_id, .va_at(L, L$held[1:5, ])$cell_id)
  # A summary whose geometry is not the cells' polygons.
  on_points <- sf::st_set_geometry(flat, sf::st_geometry(L$seeds)[match(flat$cell_id,
                                                                       L$seeds$cell_id)])
  expect_identical(
    conditionMessage(.va_error(.va_at(L, L$held, values = on_points, seeds = NULL))),
    polygons_text)
  # With seeds the summary needs no geometry.
  expect_identical(.va_at(L, L$held, values = flat), .va_at(L, L$held))
  # No value of the column: the median alone, or a column that was not summarised.
  med <- suppressWarnings(summarize_by_cell(L$asg, "z", id_col = "cell_id", agg_funs = median,
                                            estimate = "mean"))
  expect_false("resp_mean_z" %in% names(med))
  expect_identical(conditionMessage(.va_error(.va_at(L, L$held, values = med))), paste0(
    "cell_values_at(): `values` has no cell value of 'z' (neither `..pooled_resp_z` nor ",
    "`resp_mean_z`). Summarise it with summarize_by_cell() first, with the mean in `agg_funs` ",
    "or with estimate = \"pooled\"."))
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, response = "w"))),
               "`values` has no cell value of 'w' (neither `..pooled_resp_w` nor `resp_mean_w`)",
               fixed = TRUE)
  # A column the summary holds as a predictor, without its mean: the predictor's names.
  med <- suppressWarnings(summarize_by_cell(L$asg, "z", predictor_vars = "w", id_col = "cell_id",
                                            agg_funs = median, estimate = "mean"))
  expect_true(all(c("resp_median_z", "pred_median_w") %in% names(med)))
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, values = med, response = "w"))),
               "`values` has no cell value of 'w' (neither `..pooled_pred_w` nor `pred_mean_w`)",
               fixed = TRUE)
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, values = med))),
               "`values` has no cell value of 'z' (neither `..pooled_resp_z` nor `resp_mean_z`)",
               fixed = TRUE)
  # Held as response and as predictor: the response's names.
  both <- suppressWarnings(summarize_by_cell(L$asg, "z", predictor_vars = c("z", "w"),
                                             id_col = "cell_id", agg_funs = median,
                                             estimate = "mean"))
  expect_true(all(c("resp_median_z", "pred_median_z") %in% names(both)))
  expect_match(conditionMessage(.va_error(.va_at(L, L$held, values = both))),
               "`values` has no cell value of 'z' (neither `..pooled_resp_z` nor `resp_mean_z`)",
               fixed = TRUE)
})

test_that("what the other arguments must be", {
  L <- .va_full()
  msg <- function(expr) conditionMessage(.va_error(expr))
  flat <- sf::st_drop_geometry(L$held)
  expect_identical(msg(.va_at(L, flat)), "cell_values_at(): `newdata` must be an sf object.")
  expect_identical(msg(.va_at(L, sf::st_geometry(L$held))),
                   "cell_values_at(): `newdata` must be an sf object.")
  expect_identical(msg(.va_at(L, L$held, points = sf::st_drop_geometry(L$asg))),
                   "cell_values_at(): `points` must be an sf object.")
  expect_match(msg(.va_at(L, L$held, values = list(1))), "`values` must be the summarize_by_cell() result",
               fixed = TRUE)
  for (bad in list(NULL, c("z", "w"), 1, NA_character_, ""))
    expect_match(msg(.va_at(L, L$held, response = bad)),
                 "cell_values_at(): `response_var` must be a single column name", fixed = TRUE)
  expect_match(msg(cell_values_at(L$held, L$pooled, L$asg, "z", id_col = c("a", "b"))),
               "cell_values_at(): `id_col` must be a single column name", fixed = TRUE)
  expect_identical(msg(.va_at(L, L$held, response = "nope")),
                   "cell_values_at(): `points` has no column 'nope'.")
  # The ID column is found as summarize_by_cell() finds it.
  expect_identical(cell_values_at(L$held, L$pooled, L$asg, "z", seeds = L$seeds),
                   .va_at(L, L$held))
  expect_identical(cell_values_at(L$held, L$pooled, L$asg, "z", id_col = "nope", seeds = L$seeds),
                   .va_at(L, L$held))
  none <- L$asg[setdiff(names(L$asg), "cell_id")]
  expect_identical(msg(.va_at(L, L$held, points = none)), paste0(
    "cell_values_at(): could not find an ID column in `points`. Looked for: cell_id, poly_id, ",
    "polygon_id."))
  other <- L$pooled
  names(other)[names(other) == "cell_id"] <- "unit"
  expect_match(msg(.va_at(L, L$held, values = other)),
               "cell_values_at(): `values` has no ID column 'cell_id'", fixed = TRUE)
  # No error carries a call.
  expect_null(.va_error(.va_at(L, flat))$call)
  expect_null(.va_error(.va_at(L, L$held, response = "nope"))$call)
})

test_that("no scored row is an error: no value, no cell, no location, or a column that is not a number", {
  L <- .va_full()
  text <- "cell_values_at(): `points` has no row with a value of 'z', a cell and a location."
  msg <- function(points) conditionMessage(.va_error(.va_at(L, L$held, points = points)))
  rows <- L$asg; rows$z <- NA_real_
  expect_identical(msg(rows), text)
  rows <- L$asg; rows$z <- Inf
  expect_identical(msg(rows), text)
  rows <- L$asg; rows$cell_id <- NA_integer_
  expect_identical(msg(rows), text)
  # A column that is not a number has no value the summary could have read: a string, and
  # the types is.finite() accepts, a logical and a factor.
  rows <- L$asg; rows$z <- as.character(rows$z)
  expect_identical(msg(rows), text)
  rows <- L$asg; rows$z <- rows$z > 0
  expect_identical(msg(rows), text)
  rows <- L$asg; rows$z <- factor(round(rows$z))
  expect_identical(msg(rows), text)
  expect_identical(msg(L$asg[0, ]), text)
  rows <- sf::st_set_geometry(L$asg, sf::st_sfc(rep(list(sf::st_point()), nrow(L$asg)),
                                                 crs = sf::st_crs(L$asg)))
  expect_identical(msg(rows), text)
  # One scored row is enough.
  rows <- L$asg; rows$z[-1L] <- NA
  at <- .va_at(L, L$held[1:3, ], points = rows)
  expect_identical(attr(at, "values_at")$n, 1L)
  expect_identical(attr(at, "values_at")$spacing, NA_real_)
})

test_that("a row of `points` without a usable location is no scored row", {
  L <- .va_full()
  # Rows 1 to 4: an empty point, a y that is not a number, an x that is not a number (sf
  # holds both as empty points), and a y too large to take a distance from.
  rows <- L$asg
  g  <- sf::st_geometry(rows)
  xy <- .va_xy(rows)
  g[[1L]] <- sf::st_point()
  g[[2L]] <- sf::st_point(c(xy[2L, 1L], NaN))
  g[[3L]] <- sf::st_point(c(NaN, xy[3L, 2L]))
  g[[4L]] <- sf::st_point(c(xy[4L, 1L], 1e200))
  rows <- sf::st_set_geometry(rows, g)
  for (seeds in list(L$seeds, NULL)) {
    r <- .va_catch(.va_at(L, L$held, points = rows, seeds = seeds))
    expect_identical(.va_classes(r), character(0))
    expect_identical(attr(r$value, "values_at")$n, nrow(L$asg) - 4L)
    expect_identical(attr(r$value, "values_at")$locations, nrow(L$asg) - 4L)
    expect_identical(r$value, .va_at(L, L$held, points = L$asg[-(1:4), ], seeds = seeds))
  }
})

test_that("`newdata` without a CRS is refused beside rows with one, and runs beside rows without", {
  L <- .va_full()
  bare <- function(x) sf::st_set_crs(x, NA)
  want <- .va_at(L, L$held)
  expect_identical(
    conditionMessage(.va_error(.va_at(L, bare(L$held)))),
    paste0("cell_values_at(): `newdata` has no CRS and `points` has one; set it with ",
           "sf::st_set_crs() so the locations can be placed."))
  expect_identical(
    conditionMessage(.va_error(.va_at(L, bare(L$held), seeds = NULL))),
    paste0("cell_values_at(): `newdata` has no CRS and `points` has one; set it with ",
           "sf::st_set_crs() so the locations can be placed."))
  # Said before the rows are read: also where no row is scored.
  unscored <- L$asg
  unscored$z <- NA_real_
  expect_match(conditionMessage(.va_error(.va_at(L, bare(L$held), points = unscored))),
               "cell_values_at(): `newdata` has no CRS and `points` has one", fixed = TRUE)
  # Neither has one: the coordinates are used as they are.
  r <- .va_catch(.va_at(L, bare(L$held), points = bare(L$asg), seeds = bare(L$seeds)))
  expect_identical(.va_classes(r), character(0))
  expect_identical(r$value, want)
  # Seeds without a CRS beside rows with one cannot be placed either.
  expect_identical(
    conditionMessage(.va_error(.va_at(L, L$held, seeds = bare(L$seeds)))),
    paste0("cell_values_at(): `seeds` has no CRS and `points` has one; set it with ",
           "sf::st_set_crs() so the seeds can be placed."))
  # Nothing with a CRS, and no record of one: the coordinates as they are, and no CRS in
  # the record of the call.
  values <- sf::st_drop_geometry(L$pooled)
  rule <- attr(values, "count_rule"); rule$crs <- sf::NA_crs_
  attr(values, "count_rule") <- rule
  pooled <- attr(values, "pooled"); pooled$crs <- sf::NA_crs_
  attr(values, "pooled") <- pooled
  at <- .va_at(L, bare(L$held), values = values, points = bare(L$asg), seeds = bare(L$seeds))
  expect_true(is.na(attr(at, "values_at")$crs))
  expect_identical(.va_cols(at), .va_cols(want))
  # There a layer that does carry a CRS cannot be placed beside the rows.
  expect_identical(
    conditionMessage(.va_error(.va_at(L, L$held, values = values, points = bare(L$asg),
                                      seeds = bare(L$seeds)))),
    paste0("cell_values_at(): `newdata` has a CRS and `points` has none; set the CRS of ",
           "`points` with sf::st_set_crs() so the locations can be placed."))
  expect_match(
    conditionMessage(.va_error(.va_at(L, bare(L$held), values = values, points = bare(L$asg)))),
    "cell_values_at(): `seeds` has a CRS and `points` has none", fixed = TRUE)
})

test_that("rows without a CRS whose coordinates look like lon/lat are read as the summary read them, and `newdata` with them", {
  L <- .va_keep("lonlat", {
    F <- .va_full()
    rows <- sf::st_transform(F$rows, 4326)
    c(list(rows = rows, held = sf::st_transform(F$held, 4326)), .va_path(rows))
  })
  bare <- function(x) sf::st_set_crs(x, NA)
  want <- .va_at(L, L$held)
  r <- .va_catch(.va_at(L, bare(L$held), points = bare(L$asg), seeds = bare(L$seeds)))
  # ensure_projected() says what it took the rows for; nothing else is raised.
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), "look like lon/lat", fixed = TRUE)
  expect_identical(r$value, want)
  expect_gt(stats::median(r$value$dist_nearest), 100)
  # A layer with a CRS beside them is transformed as usual.
  r <- .va_catch(.va_at(L, L$held, points = bare(L$asg)))
  expect_identical(r$value, want)
  # By the polygons, too.
  r <- .va_catch(.va_at(L, bare(L$held), points = bare(L$asg), seeds = NULL))
  expect_length(r$warnings, 1L)
  expect_match(conditionMessage(r$warnings[[1L]]), "ensure_projected(): input has no CRS",
               fixed = TRUE)
  expect_identical(r$value, .va_at(L, L$held, seeds = NULL))
})

test_that("scored rows on one location: a distance, no spacing and no ratio; the first cell in the order of the IDs", {
  rows <- .va_points(rep(0, 6L), rep(0, 6L), cell_id = c(2L, 1L, 2L, 1L, 2L, 1L),
                     z = c(20, 10, 22, 12, 24, 14))
  values <- .va_means(rows)
  seeds <- .va_points(c(-100, 100, 0), c(0, 0, 1000), cell_id = c(1L, 2L, 3L))
  newdata <- .va_points(c(-60, 70, 0), c(0, 0, 900))
  r <- .va_catch(cell_values_at(newdata, values, rows, "z", id_col = "cell_id", seeds = seeds))
  expect_identical(.va_classes(r), character(0))
  at <- r$value
  expect_identical(at$cell_id, c(1L, 2L, 3L))
  expect_identical(at$value, c(12, 22, 12))
  expect_identical(at$filled, c(FALSE, FALSE, TRUE))
  # Rows of two cells share the location: the cell that comes first in the order of the IDs.
  expect_identical(at$filled_from, c(NA, NA, 1L))
  expect_identical(at$dist_nearest, c(60, 70, 900))
  expect_identical(at$dist_ratio, rep(NA_real_, 3L))
  rec <- attr(at, "values_at")
  expect_identical(rec$n, 6L)
  expect_identical(rec$locations, 1L)
  expect_identical(rec$spacing, NA_real_)
  expect_null(rec$count_rule)
  # Whatever the order of the rows.
  for (o in list(6:1, c(1L, 3L, 5L, 2L, 4L, 6L)))
    expect_identical(cell_values_at(newdata, values, rows[o, ], "z", id_col = "cell_id",
                                    seeds = seeds), at)
  # Strings are ordered as in the C locale, a factor by its levels.
  as_ids <- function(x, map) { x$cell_id <- map[x$cell_id]; x }
  chr <- c("b", "B", "a")
  at_chr <- cell_values_at(newdata, .va_means(as_ids(rows, chr)), as_ids(rows, chr), "z",
                           id_col = "cell_id", seeds = as_ids(seeds, chr))
  expect_identical(at_chr$cell_id, c("b", "B", "a"))
  expect_identical(at_chr$filled_from, c(NA, NA, "B"))
  expect_identical(at_chr$value, c(12, 22, 22))
  fac <- factor(c("u", "v", "w"), levels = c("v", "u", "w"))
  at_fac <- cell_values_at(newdata, .va_means(as_ids(rows, fac)), as_ids(rows, fac), "z",
                           id_col = "cell_id", seeds = as_ids(seeds, fac))
  expect_identical(as.character(at_fac$cell_id), c("u", "v", "w"))
  expect_identical(as.character(at_fac$filled_from), c(NA, NA, "v"))
  expect_true(is.factor(at_fac$filled_from))
})

test_that("among rows of two cells at one location the cell is the first in the C-locale order, whatever the session's collation", {
  # testthat runs a test under the C collation, where sort() gives that order too.  This
  # one runs under a collation that puts "b" before "B".
  suppressWarnings(withr::local_collate("en_US.UTF-8"))
  skip_if_not(identical(sort(c("b", "B", "a")), c("a", "b", "B")),
              "no collation here that orders \"b\" before \"B\"")
  ids <- c("b", "B")
  rows <- .va_points(rep(0, 6L), rep(0, 6L), cell_id = ids[c(2L, 1L, 2L, 1L, 2L, 1L)],
                     z = c(20, 10, 22, 12, 24, 14))
  seeds <- .va_points(c(-100, 100, 0), c(0, 0, 1000), cell_id = c("b", "B", "a"))
  at <- cell_values_at(.va_points(c(-60, 70, 0), c(0, 0, 900)), .va_means(rows), rows, "z",
                       id_col = "cell_id", seeds = seeds)
  expect_identical(at$cell_id, c("b", "B", "a"))
  expect_identical(at$filled_from, c(NA, NA, "B"))
  expect_identical(at$value, c(12, 22, 22))
})

test_that("a cell that holds scored rows is not filled, whatever `values` has for it", {
  L <- .va_full()
  # No value for any cell, as where the summary could not compute the pooled columns:
  # nothing is taken from another cell, and every value is NA.
  values <- L$pooled
  values$..pooled_resp_z <- NA_real_
  values$..obs_lo_resp_z <- NA_real_
  values$..obs_hi_resp_z <- NA_real_
  for (seeds in list(L$seeds, NULL)) {
    want <- .va_at(L, L$held, seeds = seeds)
    at   <- .va_at(L, L$held, values = values, seeds = seeds)
    expect_identical(at$cell_id, want$cell_id)
    expect_identical(at$filled, rep(FALSE, nrow(at)))
    expect_true(all(is.na(at$filled_from)))
    expect_true(all(is.na(at$value)) && all(is.na(at$obs_lo)) && all(is.na(at$obs_hi)))
    expect_identical(at$dist_nearest, want$dist_nearest)
    expect_identical(attr(at, "values_at")$filled, 0L)
    expect_true("  filled      : no location" %in% .va_print(attr(at, "values_at")))
  }
  # A location in no cell is filled all the same, with the NA of the cell it is filled from.
  far <- .va_points(c(-2e4, 3e5), c(5e4, 3e5))
  at <- .va_at(L, far, values = values, seeds = NULL)
  expect_identical(at$filled, c(TRUE, TRUE))
  expect_identical(at$filled_from, .va_at(L, far, seeds = NULL)$filled_from)
  expect_true(all(is.na(at$value)))
})

test_that("the pooled columns of a summary in which no cell holds two values give NA, unfilled", {
  # Twelve rows, each alone in its cell: the summary has plain means and no pooled value.
  rows <- .va_points(seq(0, 1100, by = 100), rep(c(0, 60, 20), 4L), cell_id = 1:12,
                     z = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8))
  r <- .va_catch(summarize_by_cell(rows, "z", id_col = "cell_id", estimate = "pooled"))
  expect_true("spatialkit_pool_unavailable" %in% .va_classes(r))
  values <- r$value
  expect_true(all(is.na(values$..pooled_resp_z)))
  expect_false(anyNA(values$resp_mean_z))
  newdata <- .va_points(c(10, 480, 1090), c(5, 30, 15))
  seeds <- sf::st_sf(cell_id = rows$cell_id, geometry = sf::st_geometry(rows))
  at <- cell_values_at(newdata, values, rows, "z", id_col = "cell_id", seeds = seeds)
  expect_identical(at$cell_id, c(1L, 6L, 12L))
  expect_identical(at$estimate, rep("pooled", 3L))
  expect_identical(at$value, rep(NA_real_, 3L))
  expect_identical(at$filled, rep(FALSE, 3L))
  expect_identical(at$filled_from, rep(NA_integer_, 3L))
  # The plain means of the same frame are there to be read.
  plain <- values[c("cell_id", "resp_mean_z")]
  expect_identical(cell_values_at(newdata, plain, rows, "z", id_col = "cell_id", seeds = seeds)$value,
                   c(3, 9, 8))
})

test_that("a cell whose valued rows have no location keeps its mean, and is filled where it has no value", {
  L <- .va_full()
  want <- .va_at(L, L$held)
  # The rows of the cell that holds the most held-out rows lose their geometry.
  k <- as.integer(names(which.max(table(want$cell_id))))
  rows <- L$asg
  g <- sf::st_geometry(rows)
  for (i in which(rows$cell_id == k)) g[[i]] <- sf::st_point()
  rows <- sf::st_set_geometry(rows, g)
  hit <- want$cell_id == k
  expect_gt(sum(hit), 1L)
  expect_gt(sum(L$asg$cell_id == k), 1L)
  # Plain means: the cell holds no scored row and has a mean, which its locations take.
  plain <- suppressWarnings(summarize_by_cell(rows, "z", id_col = "cell_id",
                                              cells_sf = L$tess$cells, estimate = "mean"))
  mean_k <- plain$resp_mean_z[plain$cell_id == k]
  expect_true(is.finite(mean_k))
  at <- .va_at(L, L$held, values = plain, points = rows)
  expect_identical(at$cell_id, want$cell_id)
  expect_identical(at$filled, rep(FALSE, nrow(at)))
  expect_identical(at$value[hit], rep(mean_k, sum(hit)))
  # No scored row lies in the cell, so its locations are further from the data.
  expect_true(all(at$dist_nearest[hit] > 0))
  expect_identical(attr(at, "values_at")$n, sum(L$asg$cell_id != k))
  # Pooled values: no row of the cell has a location, so it has no pooled value and its
  # locations are filled from another cell.
  pooled <- suppressWarnings(summarize_by_cell(rows, "z", id_col = "cell_id",
                                               cells_sf = L$tess$cells, estimate = "pooled",
                                               obs_level = 0.90))
  expect_true(is.na(pooled$..pooled_resp_z[pooled$cell_id == k]))
  at <- .va_at(L, L$held, values = pooled, points = rows)
  expect_identical(at$filled, hit)
  expect_false(anyNA(at$value))
  expect_true(all(at$filled_from[hit] != k))
})

test_that("two scored locations have a spacing: the distance between them", {
  rows <- .va_points(c(0, 30), c(0, 40), cell_id = 1:2, z = c(1, 2))
  seeds <- sf::st_sf(cell_id = 1:2, geometry = sf::st_geometry(rows))
  for (fnn in c(TRUE, FALSE)) local({
    if (!fnn) .va_no_fnn()
    at <- cell_values_at(.va_points(c(0, 30), c(10, 45)), .va_means(rows), rows, "z",
                         id_col = "cell_id", seeds = seeds)
    expect_identical(attr(at, "values_at")$locations, 2L)
    expect_identical(attr(at, "values_at")$spacing, 50)
    expect_identical(at$dist_nearest, c(10, 5))
    expect_identical(at$dist_ratio, c(0.2, 0.1))
    expect_identical(at$value, c(1, 2))
  })
})

test_that("the record counts rows and locations apart, in the frame and in print", {
  T <- .va_lattice()
  seeds <- .va_points(T$at$x, T$at$y, cell_id = T$cell)
  rec <- attr(cell_values_at(.va_points(13, 21), T$values, T$rows, "z", id_col = "cell_id",
                             seeds = seeds), "values_at")
  expect_identical(rec$n, 60L)
  expect_identical(rec$locations, 30L)
  expect_true("  scored rows : 60 rows at 30 locations, measured in EPSG:32617" %in% .va_print(rec))
})

test_that("an ID column named like a column of the result is refused", {
  T <- .va_lattice()
  for (nm in c("value", "estimate", "obs_lo", "obs_hi", "filled", "filled_from", "dist_nearest",
               "dist_ratio")) {
    rename <- function(x) { names(x)[names(x) == "cell_id"] <- nm; x }
    seeds <- rename(.va_points(T$at$x, T$at$y, cell_id = T$cell))
    e <- .va_error(cell_values_at(.va_points(13, 21), rename(T$values), rename(T$rows), "z",
                                  id_col = nm, seeds = seeds))
    expect_identical(conditionMessage(e), paste0(
      "cell_values_at(): the ID column is named '", nm, "', which is a column of the result ",
      "(value, estimate, obs_lo, obs_hi, filled, filled_from, dist_nearest, dist_ratio). ",
      "Rename it in `points`, `values` and `seeds`."))
    expect_null(e$call)
  }
  # Any other name is the first column of the frame.
  rename <- function(x) { names(x)[names(x) == "cell_id"] <- "values"; x }
  at <- cell_values_at(.va_points(13, 21), rename(T$values), rename(T$rows), "z", id_col = "values",
                       seeds = rename(.va_points(T$at$x, T$at$y, cell_id = T$cell)))
  expect_identical(names(at), c("values", "value", "estimate", "filled", "filled_from",
                                "dist_nearest", "dist_ratio"))
  expect_false(anyDuplicated(names(at)) > 0L)
})


# ---- The bound on the exact searches ------------------------------------------------

test_that("without FNN each exact search is bounded, and the spacing's does not depend on the locations asked for", {
  # 20 scored locations on a 5 x 4 lattice in 4 cells, and the 4 cells' seeds.
  g <- expand.grid(x = seq(0, 400, by = 100), y = seq(0, 300, by = 100))
  cell <- 1L + (g$x > 200) + 2L * (g$y > 100)
  rows <- .va_points(g$x, g$y, cell_id = cell, z = g$x / 100 + g$y / 50)
  values <- .va_means(rows)
  seeds <- .va_points(c(100, 350, 100, 350), c(50, 50, 250, 250), cell_id = 1:4)
  new <- function(k) .va_points(seq(10, by = 37, length.out = k), seq(15, by = 23, length.out = k))
  call <- function(k, seeds_k = seeds)
    cell_values_at(new(k), values, rows, "z", id_col = "cell_id", seeds = seeds_k)
  with_fnn <- if (requireNamespace("FNN", quietly = TRUE)) call(6L) else NULL
  .va_no_fnn()
  free <- call(6L)
  if (!is.null(with_fnn)) expect_identical(free, with_fnn)

  withr::local_options(spatialkit.exact_pairs_max = 100)
  # One new location: 20 pairs and 4 pairs for its two searches, and 400 for the spacing.
  e <- .va_error(call(1L))
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$what, "spacing")
  expect_identical(e$work, 400)
  expect_identical(e$limit, 100)
  expect_null(e$call)
  expect_identical(conditionMessage(e), paste0(
    "cell_values_at(): the spacing of the scored rows (the nearest other of each of 20 scored ",
    "locations) would compare 4e+02 pairs without the FNN package, above the limit of 1e+02 ",
    "(options(spatialkit.exact_pairs_max = )). Install FNN: the result is the same."))
  # By the polygons as well, where no seed is searched for.
  cells <- sf::st_sf(cell_id = 1:4, geometry = sf::st_buffer(sf::st_geometry(seeds), 100,
                                                            endCapStyle = "SQUARE"))
  with_cells <- dplyr::left_join(cells, values, by = "cell_id")
  e <- .va_error(cell_values_at(new(1L), with_cells, rows, "z", id_col = "cell_id"))
  expect_identical(e$what, "spacing")
  # Six new locations: 120 pairs for the nearest scored location, which is named first.
  e <- .va_error(call(6L))
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$what, "nearest scored location")
  expect_identical(e$work, 120)
  expect_match(conditionMessage(e),
               "cell_values_at(): finding the nearest of 20 scored locations for 6 locations would",
               fixed = TRUE)
  # The limit is the most pairs allowed, not the first refused.
  withr::local_options(spatialkit.exact_pairs_max = 400)
  expect_identical(call(6L), free)
  withr::local_options(spatialkit.exact_pairs_max = 399)
  expect_identical(.va_error(call(6L))$what, "spacing")
  # 101 seeds for four locations: 404 pairs, beside 80 and 400.
  withr::local_options(spatialkit.exact_pairs_max = 400)
  many <- rbind(seeds, .va_points(seq(1000, by = 10, length.out = 97), rep(1000, 97),
                                  cell_id = 100L + 1:97))
  e <- .va_error(call(4L, many))
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$what, "nearest seed")
  expect_identical(e$work, 404)
  expect_match(conditionMessage(e),
               "cell_values_at(): finding the nearest of 101 seeds for 4 locations would",
               fixed = TRUE)
  # No pair at all under a limit of 0, but for a call that has none to compare.
  withr::local_options(spatialkit.exact_pairs_max = 0)
  expect_identical(.va_error(call(1L))$what, "nearest scored location")
})

test_that("the bound counts the rows and the seeds that have a location, and no others", {
  # 20 scored locations.  Of 32 rows of `newdata` two have a location, and of 325 seeds 25.
  g <- expand.grid(x = seq(0, 400, by = 100), y = seq(0, 300, by = 100))
  rows <- .va_points(g$x, g$y, cell_id = 1L + (g$x > 200) + 2L * (g$y > 100),
                     z = g$x / 100 + g$y / 50)
  values <- .va_means(rows)
  empty <- function(k) sf::st_sfc(rep(list(sf::st_point()), k), crs = 32617)
  located <- .va_points(c(10, 47), c(15, 38))
  newdata <- rbind(sf::st_sf(geometry = sf::st_geometry(located)), sf::st_sf(geometry = empty(30L)))
  seeds <- rbind(.va_points(c(100, 350, 100, 350), c(50, 50, 250, 250), cell_id = 1:4),
                 .va_points(seq(1000, by = 10, length.out = 21L), rep(1000, 21L),
                            cell_id = 100L + 1:21),
                 sf::st_sf(cell_id = 1000L + 1:300, geometry = empty(300L)))
  call <- function() cell_values_at(newdata, values, rows, "z", id_col = "cell_id", seeds = seeds)
  .va_no_fnn()
  free <- call()
  expect_identical(nrow(free), 32L)
  expect_identical(.va_cols(free[1:2, ]),
                   .va_cols(cell_values_at(located, values, rows, "z", id_col = "cell_id",
                                           seeds = seeds[1:4, ])))
  # 2 x 20 pairs for the nearest scored location, 2 x 25 for the nearest seed, 20 x 20
  # for the spacing: all within 400.
  withr::local_options(spatialkit.exact_pairs_max = 400)
  expect_identical(call(), free)
  # The seeds' search is stopped at 2 x 25 pairs, not at 2 x 325.
  withr::local_options(spatialkit.exact_pairs_max = 45)
  e <- .va_error(call())
  expect_s3_class(e, "spatialkit_work_limit")
  expect_identical(e$what, "nearest seed")
  expect_identical(e$work, 50)
  expect_match(conditionMessage(e),
               "cell_values_at(): finding the nearest of 25 seeds for 2 locations would",
               fixed = TRUE)
  # The scored locations' search at 2 x 20 pairs, not at 32 x 20.
  withr::local_options(spatialkit.exact_pairs_max = 39)
  e <- .va_error(call())
  expect_identical(e$what, "nearest scored location")
  expect_identical(e$work, 40)
  expect_match(conditionMessage(e),
               "cell_values_at(): finding the nearest of 20 scored locations for 2 locations would",
               fixed = TRUE)
})

test_that("with FNN no search is bounded", {
  skip_if_not_installed("FNN")
  g <- expand.grid(x = seq(0, 400, by = 100), y = seq(0, 300, by = 100))
  rows <- .va_points(g$x, g$y, cell_id = 1L + (g$x > 200) + 2L * (g$y > 100),
                     z = g$x / 100 + g$y / 50)
  seeds <- .va_points(c(100, 350, 100, 350), c(50, 50, 250, 250), cell_id = 1:4)
  newdata <- .va_points(c(10, 47, 84), c(15, 38, 61))
  want <- cell_values_at(newdata, .va_means(rows), rows, "z", id_col = "cell_id", seeds = seeds)
  withr::local_options(spatialkit.exact_pairs_max = 0)
  expect_identical(cell_values_at(newdata, .va_means(rows), rows, "z", id_col = "cell_id",
                                  seeds = seeds), want)
})

test_that("the backend is asked for once, and the answer reaches every search", {
  skip_if_not_installed("FNN")
  L <- .va_full()
  asked <- character(0)
  testthat::local_mocked_bindings(.has_backend = function(pkg) {
    asked <<- c(asked, pkg)
    FALSE
  })
  # A kd-tree reached by any search stops the call.
  testthat::local_mocked_bindings(
    get.knn  = function(...) stop("FNN::get.knn() reached"),
    get.knnx = function(...) stop("FNN::get.knnx() reached"), .package = "FNN")
  at <- .va_at(L, L$held)
  expect_identical(asked, "FNN")
  expect_identical(nrow(at), nrow(L$held))
})


# ---- Ties ---------------------------------------------------------------------------

test_that("a location the same distance from two scored locations takes the lower x, then the lower y", {
  T <- .va_lattice()
  # No location has a cell: the one seed stands for none.
  seeds <- .va_points(25, 20, cell_id = NA_integer_)
  # Half-way between two locations along x, along y, between four, and beside the last column.
  q <- data.frame(x = c(5, 0, 5, 25, 45, 50, 13), y = c(0, 5, 5, 15, 40, 35, 21))
  newdata <- .va_points(q$x, q$y)
  u <- cbind(T$at$x, T$at$y)
  first <- .va_nearest(cbind(q$x, q$y), u)
  # The ties are ties: the first six are the same distance from two or four locations.
  tied <- vapply(seq_len(nrow(q)), function(i) {
    d2 <- (q$x[i] - u[, 1L])^2 + (q$y[i] - u[, 2L])^2
    sum(d2 == min(d2))
  }, integer(1))
  expect_identical(tied, c(2L, 2L, 4L, 4L, 2L, 2L, 1L))
  want_cell <- T$cell[first]
  expect_identical(want_cell, c(30L, 30L, 30L, 19L, 6L, 2L, 23L))
  # The other tied location has the lower ID at (45, 40) and comes first in neither the
  # rows nor the IDs elsewhere, so no other rule gives these cells.
  call <- function(rows) cell_values_at(newdata, T$values, rows, "z", id_col = "cell_id",
                                        seeds = seeds)
  at <- call(T$rows)
  expect_identical(at$filled, rep(TRUE, 7L))
  expect_identical(at$cell_id, rep(NA_integer_, 7L))
  expect_identical(at$filled_from, want_cell)
  expect_identical(at$value, 10 * want_cell)
  expect_identical(at$dist_nearest, sqrt((q$x - u[first, 1L])^2 + (q$y - u[first, 2L])^2))
  expect_identical(attr(at, "values_at")$locations, 30L)
  expect_identical(attr(at, "values_at")$spacing, 10)
  # Under permutations of the rows.
  for (seed in 1:5) {
    o <- local({ .va_seed(6210 + seed); sample.int(nrow(T$rows)) })
    expect_identical(call(T$rows[o, ]), at)
  }
  expect_identical(call(T$rows[order(T$rows$cell_id), ]), at)
  expect_identical(call(T$rows[order(-T$rows$cell_id), ]), at)
  # Without FNN.
  local({
    .va_no_fnn()
    expect_identical(call(T$rows), at)
    expect_identical(call(T$rows[order(T$rows$cell_id), ]), at)
  })
  # The lower x before the lower y: without the rows at the corner, a location on the
  # diagonal is the same distance from (0, 10), cell 29, and from (10, 0), cell 25.
  xy <- .va_xy(T$rows)
  no_corner <- T$rows[!(xy[, 1L] == 5e5 & xy[, 2L] == 5e6), ]
  expect_identical(nrow(no_corner), 58L)
  diagonal <- .va_points(c(2, 4.5), c(2, 4.5))
  for (fnn in c(TRUE, FALSE)) local({
    if (!fnn) .va_no_fnn()
    at <- cell_values_at(diagonal, .va_means(no_corner), no_corner, "z", id_col = "cell_id",
                         seeds = seeds)
    expect_identical(at$filled_from, c(29L, 29L))
    expect_identical(at$dist_nearest, sqrt(c(2^2 + 8^2, 4.5^2 + 5.5^2)))
  })
  # Rows of a second cell at a location: the cell that comes first in the order of the IDs.
  more <- rbind(T$rows, .va_points(c(40, 40, 0, 0), c(40, 40, 0, 0),
                                   cell_id = c(3L, 3L, 40L, 40L), z = c(1, 2, 3, 4)))
  values <- .va_means(more)
  for (o in list(seq_len(nrow(more)), rev(seq_len(nrow(more))))) {
    at2 <- cell_values_at(newdata, values, more[o, ], "z", id_col = "cell_id", seeds = seeds)
    expect_identical(at2$filled_from, c(30L, 30L, 30L, 19L, 3L, 2L, 23L))
    expect_identical(at2$value, values$resp_mean_z[match(at2$filled_from, values$cell_id)])
    expect_identical(attr(at2, "values_at")$locations, 30L)
    expect_identical(attr(at2, "values_at")$n, 64L)
  }
})

test_that("a location the same distance from two seeds takes the first of them in `seeds`", {
  T <- .va_lattice()
  # A seed on every location, with the ID of the cell there.
  seeds <- .va_points(T$at$x, T$at$y, cell_id = T$cell)
  q <- data.frame(x = c(5, 0, 5, 25, 45, 50, 13), y = c(0, 5, 5, 15, 40, 35, 21))
  newdata <- .va_points(q$x, q$y)
  u <- cbind(T$at$x, T$at$y)
  first_of <- function(order_of_seeds) vapply(seq_len(nrow(q)), function(i) {
    d2 <- (q$x[i] - u[, 1L])^2 + (q$y[i] - u[, 2L])^2
    T$cell[order_of_seeds[order_of_seeds %in% which(d2 == min(d2))][1L]]
  }, integer(1))
  orders <- list(seq_len(30L), rev(seq_len(30L)), local({ .va_seed(6220); sample.int(30L) }))
  got <- list()
  for (o in orders) {
    at <- cell_values_at(newdata, T$values, T$rows, "z", id_col = "cell_id", seeds = seeds[o, ])
    expect_identical(at$cell_id, first_of(o))
    expect_identical(at$filled, rep(FALSE, 7L))
    expect_identical(at$value, 10 * at$cell_id)
    got[[length(got) + 1L]] <- at$cell_id
    # Under a permutation of the rows, and without FNN.
    expect_identical(cell_values_at(newdata, T$values, T$rows[rev(seq_len(60L)), ], "z",
                                    id_col = "cell_id", seeds = seeds[o, ]), at)
    local({
      .va_no_fnn()
      expect_identical(cell_values_at(newdata, T$values, T$rows, "z", id_col = "cell_id",
                                      seeds = seeds[o, ]), at)
    })
  }
  # The order of the seeds decides, so two orders give two answers at a tie and one elsewhere.
  expect_identical(got[[1L]], c(30L, 30L, 30L, 19L, 6L, 2L, 23L))
  expect_identical(got[[2L]], c(25L, 29L, 24L, 13L, 1L, 1L, 23L))
})


# ---- Coverage at held-out rows of simulated layers ----------------------------------

test_that("the observation interval covers a held-out observation about as often as its level says, and the t interval does not", {
  skip_on_cran()
  # Forty layers, r = 1 to 40.  set.seed(r): 3,000 locations, x then y uniform on a 1e5
  # square.  The response is 3 times a field of 60 cosine waves plus normal noise of
  # standard deviation 1; the field is drawn under set.seed(r + 1000), the noise after
  # it, and then the held-out fifth, sample.int(3000, 600).
  field <- function(xy, seed) {
    set.seed(seed)
    om <- matrix(stats::rnorm(120, sd = 1 / 9000), 60L)
    phase <- stats::runif(60, 0, 2 * pi)
    as.numeric(sqrt(2 / 60) * cos(xy %*% t(om) + rep(phase, each = nrow(xy))) %*%
                 stats::rnorm(60))
  }
  one <- function(r) {
    .va_seed(r)
    xy <- cbind(stats::runif(3000, 0, 1e5), stats::runif(3000, 0, 1e5))
    z  <- 3 * field(xy, r + 1000) + stats::rnorm(3000, sd = 1)
    held <- sample.int(3000L, 600L)
    pts <- sf::st_as_sf(data.frame(x = 5e5 + xy[, 1L], y = 5e6 + xy[, 2L], z = z),
                        coords = c("x", "y"), crs = 32617)
    rows <- pts[-held, ]
    new  <- pts[held, ]
    prof <- resolution_profile(rows, "z", partition = "density", init = "kmeans++", nstart = 3,
                               sac = "auto", seed = 123, estimate = "pooled", obs_per_cell = 4.5)
    sel   <- select_resolution(prof, criterion = "size")
    seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = rows)
    tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    asg   <- assign_features_to_polygons(rows, tess$cells, polygon_id_col = "cell_id",
                                         quiet = TRUE)
    out   <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                               estimate = "pooled", obs_level = 0.90, conf_level = 0.90)
    seeds$cell_id <- tess$index
    at <- cell_values_at(new, out, asg, "z", id_col = "cell_id", seeds = seeds)
    # The t interval of the cell the value came from, read as a range for an observation;
    # a cell without one covers nothing.
    src <- match(ifelse(at$filled, at$filled_from, at$cell_id), out$cell_id)
    lo <- out$..ci_lo_resp_z[src]
    hi <- out$..ci_hi_resp_z[src]
    c(cells = sel$best,
      obs = mean(new$z >= at$obs_lo & new$z <= at$obs_hi),
      ci  = mean(!is.na(lo) & !is.na(hi) & new$z >= lo & new$z <= hi))
  }
  got <- vapply(1:40, one, numeric(3))
  expect_identical(unname(got["cells", ]), rep(533, 40L))
  expect_false(anyNA(got))
  # One run of these seeds gave 0.891 (single layers 0.850 to 0.915) and 0.578 (0.515 to 0.623).
  expect_gt(mean(got["obs", ]), 0.87)
  expect_lt(mean(got["obs", ]), 0.92)
  expect_lt(mean(got["ci", ]), 0.70)
})
