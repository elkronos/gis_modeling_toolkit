# The neighbour searches say where the rounding of the coordinates at their
# offset can have parted a tie.
#
# Two squared distances within .TIE_TOL (1e-7) of each other are a tie in
# every neighbour search.  That holds at any origin while the doubles hold
# the coordinates to about .TIE_TOL / 4 of the distances compared.  An 8 x 8
# lattice at a spacing of 100.1 m moved by 1e11 is stored rounded to
# 1.5e-5 m, nominally tied neighbours lie further apart than .TIE_TOL, and
# which of them a row takes follows the rounding: the k-nearest-neighbour
# weights, every Moran's I, the density partition, the pooling's nearest
# other location and the kriging neighbourhoods moved with the origin, with
# no warning.  Each search now hands the k-th squared distance and the next
# one past its tie to .tie_note(), and the exported call warns once, with
# class "spatialkit_xy_imprecise".  Nothing it returns changes.
#
# Every expectation on a warning compares inputs at least a factor of ten
# inside or outside the band of .tie_lost_by_rounding(), so no last bit that
# moves with the platform decides it.

.tp_lattice_xy <- function(sp, off) {
  g <- expand.grid(i = 0:7, j = 0:7)
  cbind(sp * g$i + off, sp * g$j + off)
}

.tp_sf <- function(xy, crs = 32615) {
  g <- expand.grid(i = 0:7, j = 0:7)
  stopifnot(nrow(xy) == nrow(g))
  sf::st_as_sf(data.frame(z = sin(g$i) + cos(1.3 * g$j), x1 = g$i / 7 + 0.1 * cos(g$j),
                          X = xy[, 1], Y = xy[, 2]),
               coords = c("X", "Y"), crs = crs)
}

.tp_fit <- function(d) {
  new_spatial_fit("tp_lm_fit", stats::lm(z ~ x1, data = sf::st_drop_geometry(d)),
                  z ~ x1, "z", "x1", d)
}
residuals.tp_lm_fit <- function(object, ...) as.numeric(stats::residuals(object$engine))
fitted.tp_lm_fit <- function(object, ...) as.numeric(stats::fitted(object$engine))
registerS3method("residuals", "tp_lm_fit", residuals.tp_lm_fit)
registerS3method("fitted", "tp_lm_fit", fitted.tp_lm_fit)

# The value of `expr` and the spatialkit_xy_imprecise warnings it raised;
# every other warning is muffled.
.tp_catch <- function(expr) {
  w <- list()
  v <- withCallingHandlers(expr, warning = function(cnd) {
    if (inherits(cnd, "spatialkit_xy_imprecise")) w[[length(w) + 1L]] <<- cnd
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

# Closes, when the calling test ends, a watch that a failed expectation
# left open, so that it does not swallow the warnings of the tests after it.
.tp_close_watch <- function(env = parent.frame()) {
  withr::defer({
    .tie_state$open <- FALSE
    .tie_state$said <- FALSE
    .tie_state$spacing <- 0
    .tie_state$dist <- Inf
  }, envir = env)
}

# `expr` run inside a watch of its own, as an exported function runs it:
# its value and whether the watch warned.  A watch that `expr` stops in is
# closed.
.tp_watched <- function(expr) {
  end <- .tie_watch("probe", "xy")
  on.exit(end(end), add = TRUE)
  warned <- FALSE
  v <- withCallingHandlers({
    r <- expr
    end(r)
    r
  }, spatialkit_xy_imprecise = function(cnd) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  list(value = v, warned = warned)
}

# The parted lattices: spacing and offset.
.tp_parted <- list(c(100.1, 1e11), c(100.1, 1e12), c(0.7, 1e10))


test_that(".tie_lost_by_rounding() flags a pair just past the tie line and nothing else", {
  M <- 1e11
  r <- .ulp_of(M)
  a <- 1e4
  line <- a * (1 + .TIE_TOL)
  band <- 8 * r * 2 * sqrt(a)               # 8 r (sqrt(a) + sqrt(b)), b about a
  expect_true(.tie_lost_by_rounding(a, line + band / 10, M))
  expect_false(.tie_lost_by_rounding(a, line + 10 * band, M))
  # A pair within the tie stays tied whatever the rounding did: one-sided.
  expect_false(.tie_lost_by_rounding(a, a * (1 + .TIE_TOL / 10), M))
  expect_false(.tie_lost_by_rounding(a, a, M))
  expect_false(.tie_lost_by_rounding(0, 0, M))
  # A distance that is not finite (none past the tie, or an overflow) is not
  # read.
  expect_false(.tie_lost_by_rounding(a, Inf, M))
  expect_false(.tie_lost_by_rounding(NA_real_, a, M))
  expect_identical(.tie_lost_by_rounding(c(a, a, a), c(line + band / 10, line + 10 * band, a), M),
                   c(TRUE, FALSE, FALSE))
  # The spacing of doubles is that of the largest coordinate: the same pair
  # is clear at an origin of 1e6.
  expect_false(.tie_lost_by_rounding(a, line + band / 10, 1e6))
  # A frame divided by a scale passes the spacing over that scale.
  expect_true(.tie_lost_by_rounding(a / 1e6, (line + band / 10) / 1e6, M, r = r / 1e3))
})

test_that("each neighbour search flags the parted lattice and returns what it returned", {
  for (cs in .tp_parted) {
    xy <- .tp_lattice_xy(cs[1], cs[2])
    u  <- xy[order(xy[, 1], xy[, 2]), ]
    lbl <- sprintf("spacing %g, offset %g", cs[1], cs[2])
    for (fnn in c(TRUE, FALSE)) {
      if (fnn && !requireNamespace("FNN", quietly = TRUE)) next
      got <- .tp_watched(.knn_pairs(xy, 6L, use_fnn = fnn))
      expect_true(got$warned, info = paste(lbl, "pairs, fnn", fnn))
      expect_identical(got$value, .knn_pairs(xy, 6L, use_fnn = fnn), info = lbl)
      rnd <- function() .tie_rounding(u)
      got <- .tp_watched(.knn_exact_rows(u[, 1], u[, 2], seq_len(64), 6L, .TIE_TOL, rnd = rnd()))
      expect_true(got$warned, info = paste(lbl, "exact rows"))
      expect_identical(got$value, .knn_exact_rows(u[, 1], u[, 2], seq_len(64), 6L, .TIE_TOL),
                       info = lbl)
      got <- .tp_watched(.nearest_other(xy, rep(1L, 64), seq_len(64), fnn))
      expect_true(got$warned, info = paste(lbl, "nearest other, fnn", fnn))
      expect_identical(got$value, .nearest_other(xy, rep(1L, 64), seq_len(64), fnn), info = lbl)
      got <- .tp_watched(withr::with_seed(1L, .density_kmeans(xy, 5L, use_fnn = fnn)))
      expect_true(got$warned, info = paste(lbl, "density partition, fnn", fnn))
      expect_identical(got$value, withr::with_seed(1L, .density_kmeans(xy, 5L, use_fnn = fnn)),
                       info = lbl)
    }
    if (requireNamespace("FNN", quietly = TRUE)) {
      got <- .tp_watched(.knn_tied_rows(u, seq_len(64), 6L, .TIE_TOL, rnd = .tie_rounding(u)))
      expect_true(got$warned, info = paste(lbl, "tied rows"))
      expect_identical(got$value, .knn_tied_rows(u, seq_len(64), 6L, .TIE_TOL), info = lbl)
    }
    got <- .tp_watched(.sac_structure_test(xy, sin(xy[, 1] / cs[1]) + cos(xy[, 2] / cs[1])))
    expect_true(got$warned, info = paste(lbl, "structure test"))
    # The centre of the lattice is the same distance from four points at
    # each of several distances: the kriging neighbourhood's tie.
    d2  <- (u[, 1] - (cs[2] + 3.5 * cs[1]))^2 + (u[, 2] - (cs[2] + 3.5 * cs[1]))^2
    got <- .tp_watched(.tie_pick_note(d2, 6L, .tie_rounding(u)))
    expect_true(got$warned, info = paste(lbl, "kriging neighbourhood"))
  }
})

test_that("a search outside a watch reads nothing and warns nothing", {
  xy <- .tp_lattice_xy(100.1, 1e11)
  expect_null(.tie_rounding(xy))
  expect_no_warning(.knn_pairs(xy, 6L, use_fnn = FALSE))
  expect_no_warning(.nearest_other(xy, rep(1L, 64), seq_len(64), FALSE))
})

test_that("residual_morans_i() and estimate_sac_range() warn once on the parted lattice", {
  for (cs in .tp_parted) {
    d   <- .tp_sf(.tp_lattice_xy(cs[1], cs[2]))
    lbl <- sprintf("spacing %g, offset %g", cs[1], cs[2])
    spacing <- format(signif(.ulp_of(cs[2]), 3))
    got <- .tp_catch(residual_morans_i(.tp_fit(d)))
    expect_length(got$warnings, 1L)
    msg <- conditionMessage(got$warnings[[1]])
    expect_match(msg, "residual_morans_i(): the data of `fit` lies so far from the origin of its CRS",
                 fixed = TRUE, info = lbl)
    expect_match(msg, paste0("rounded to ", spacing, ", and some neighbour distances near"),
                 fixed = TRUE, info = lbl)
    expect_match(msg, "Shift the layer nearer the origin", fixed = TRUE, info = lbl)
    expect_identical(got$warnings[[1]]$spacing, .ulp_of(cs[2]), info = lbl)
    # The distance whose tie it reads is a neighbour distance of the
    # lattice: its spacing or a diagonal.
    expect_true(got$warnings[[1]]$distance >= cs[1] * 0.99 &&
                  got$warnings[[1]]$distance <= cs[1] * 3, info = lbl)
    expect_true(is.finite(got$value$observed), info = lbl)
    # The statistic is the one the same call gives with the warning
    # muffled: the warning changes nothing it returns.
    again <- withCallingHandlers(residual_morans_i(.tp_fit(d)),
                                 warning = function(cnd) invokeRestart("muffleWarning"))
    expect_identical(got$value, again, info = lbl)
  }
  skip_if_not_installed("gstat")
  for (cs in .tp_parted) {
    d   <- .tp_sf(.tp_lattice_xy(cs[1], cs[2]))
    got <- .tp_catch(suppressMessages(estimate_sac_range(d, "z")))
    expect_length(got$warnings, 1L)
    expect_match(conditionMessage(got$warnings[[1]]),
                 "estimate_sac_range(): `points_sf` lies so far from the origin of its CRS",
                 fixed = TRUE, info = sprintf("spacing %g, offset %g", cs[1], cs[2]))
  }
})

test_that("the density partition warns once through its exported caller", {
  d <- .tp_sf(.tp_lattice_xy(100.1, 1e11))
  got <- .tp_catch(suppressMessages(voronoi_seeds_kmeans(d, 5, partition = "density")))
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]),
               "voronoi_seeds_kmeans(): `points_sf` lies so far from the origin of its CRS",
               fixed = TRUE)
  ref <- .tp_catch(suppressMessages(voronoi_seeds_kmeans(.tp_sf(.tp_lattice_xy(100.1, 0)), 5,
                                                         partition = "density")))
  expect_length(ref$warnings, 0L)
})

test_that("lattices whose ties hold, and layers with none, raise nothing", {
  # The lattice at 100.1 m up to 1e10 (doubles hold it), the lattice at
  # 250 m on whole metres at 1e12 and 1e15 (exact at every offset), points
  # on 1/8 m and a jittered lattice at 1e12.
  run <- function(xy) {
    d <- .tp_sf(xy)
    .tp_catch(residual_morans_i(.tp_fit(d)))
  }
  ref <- run(.tp_lattice_xy(100.1, 0))
  expect_length(ref$warnings, 0L)
  for (off in c(1e8, 1e10)) {
    got <- run(.tp_lattice_xy(100.1, off))
    expect_length(got$warnings, 0L)
    # The weights are identical; I differs at most by the order in which
    # the same products are summed.
    expect_equal(got$value$observed, ref$value$observed, tolerance = 1e-10,
                 info = sprintf("offset %g", off))
  }
  ref <- run(.tp_lattice_xy(250, 0))
  for (off in c(1e12, 1e15)) {
    got <- run(.tp_lattice_xy(250, off))
    expect_length(got$warnings, 0L)
    expect_equal(got$value$observed, ref$value$observed, tolerance = 1e-10,
                 info = sprintf("250 m lattice, offset %g", off))
  }
  pts <- withr::with_seed(20261106L, cbind(round(8 * stats::runif(64, 0, 4000)) / 8,
                                           round(8 * stats::runif(64, 0, 3000)) / 8))
  jit <- withr::with_seed(20261107L, .tp_lattice_xy(100.1, 0) +
                            stats::runif(128, -5, 5))
  for (xy in list(pts, pts + 1e12, jit + 1e12)) {
    got <- .tp_watched(.knn_pairs(xy, 8L, use_fnn = FALSE))
    expect_false(got$warned)
    expect_false(.tp_watched(.nearest_other(xy, rep(1L, 64), seq_len(64), FALSE))$warned)
  }
  expect_length(run(jit + 1e12)$warnings, 0L)
})

test_that("the nearest-row search of cell_values_at() and predict_surface() holds the lattice", {
  # Its margin (.centres_tol(), at least four units in the last place of the
  # largest coordinate) holds more than the rounding moves a distance, so a
  # point at the centre of a lattice square takes the same point at every
  # offset, and nothing is said.
  g  <- expand.grid(i = 0:7, j = 0:7)
  inner <- g$i < 7 & g$j < 7
  for (sp in c(100.1, 0.7)) {
    ref <- NULL
    for (off in c(0, 1e10, 1e11, 1e12)) {
      xy <- .tp_lattice_xy(sp, off)
      q  <- cbind(sp * (g$i[inner] + 0.5), sp * (g$j[inner] + 0.5)) + off
      for (fnn in c(TRUE, FALSE)) {
        if (fnn && !requireNamespace("FNN", quietly = TRUE)) next
        got <- .tp_watched(.values_at_nearest(q[, 1], q[, 2], xy[, 1], xy[, 2], use_fnn = fnn))
        expect_false(got$warned)
        if (is.null(ref)) ref <- got$value
        expect_identical(got$value, ref, info = sprintf("spacing %g, offset %g", sp, off))
      }
    }
  }
})

test_that("a watch inside another leaves the warning to the outer call", {
  .tp_close_watch()
  xy <- .tp_lattice_xy(100.1, 1e11)
  n <- 0L
  withCallingHandlers({
    outer <- .tie_watch("outer_fn", "a")
    inner <- .tie_watch("inner_fn", "b")
    .knn_pairs(xy, 6L, use_fnn = FALSE)
    inner(NULL)
    outer(NULL)
  }, spatialkit_xy_imprecise = function(cnd) {
    n <<- n + 1L
    expect_match(conditionMessage(cnd), "^outer_fn\\(\\): `a` lies")
    invokeRestart("muffleWarning")
  })
  expect_identical(n, 1L)
  # A call that stops says nothing, and leaves no watch open.
  f <- function() {
    tie_end <- .tie_watch("f", "xy")
    on.exit(tie_end(returnValue(tie_end)), add = TRUE)
    .knn_pairs(xy, 6L, use_fnn = FALSE)
    stop("stopped")
  }
  expect_error(withCallingHandlers(f(), spatialkit_xy_imprecise = function(cnd)
    stop("warned")), "stopped")
  expect_false(isTRUE(.tie_state$open))
})

test_that("the help pages and NEWS say so", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  for (page in c("residual_morans_i.Rd", "estimate_sac_range.Rd"))
    expect_match(read("man", page), "past a tie by less than that rounding", fixed = TRUE,
                 info = page)
  expect_match(read("man", "fit_gwr_model.Rd"), paste0(
    "where a neighbour lies past the edge by less than that rounding the call ",
    "warns, with class \\code{\"spatialkit_xy_imprecise\"}"), fixed = TRUE)
  expect_match(read("NEWS.md"), "past a tie by less than that rounding", fixed = TRUE)
  expect_match(read("NEWS.md"), "the edge of a GWR's boxcar window moved with the origin",
               fixed = TRUE)
})

# The originals of the searches that a test below silences, so that the
# site under test alone decides whether the watch warns.
.tp_quiet_tied_rows <- function() {
  orig <- .knn_tied_rows
  function(u, rows, k, tol = 0, num = seq_len(nrow(u)), rnd = NULL)
    orig(u, rows, k, tol, num)
}

test_that("the kd-tree paths of the searches note a parted tie themselves", {
  skip_if_not_installed("FNN")
  quiet_rows <- .tp_quiet_tied_rows()
  orig_exact <- .nearest_exact
  quiet_exact <- function(x, y, cx, cy, tol = .centres_tol(cx, cy), rnd = NULL)
    orig_exact(x, y, cx, cy, tol)
  g <- expand.grid(i = 0:7, j = 0:7)
  for (cs in .tp_parted) {
    xy  <- .tp_lattice_xy(cs[1], cs[2])
    u   <- xy[order(xy[, 1], xy[, 2]), ]
    lbl <- sprintf("spacing %g, offset %g", cs[1], cs[2])
    # The rows FNN settles on its own, with the tied rows' search silenced.
    local_mocked_bindings(.knn_tied_rows = quiet_rows)
    got <- .tp_watched(.knn_locations(u, 6L, use_fnn = TRUE, rnd = .tie_rounding(u)))
    expect_true(got$warned, info = paste(lbl, ".knn_locations()"))
    expect_identical(got$value, .knn_locations(u, 6L, use_fnn = TRUE), info = lbl)
    # More than 16 centres, in the frame of .kmeans_frame(): the points of
    # the odd columns are as far from the centre on their left as from the
    # one on their right.
    fr  <- .kmeans_frame(xy)
    cen <- fr$z[g$i %% 2 == 0, ]
    local_mocked_bindings(.nearest_exact = quiet_exact)
    got <- .tp_watched(.nearest_centre(fr$z[, 1], fr$z[, 2], cen[, 1], cen[, 2], TRUE,
                                       rnd = .tie_rounding(xy, scale = fr$scale)))
    expect_true(got$warned, info = paste(lbl, ".nearest_centre()"))
    expect_identical(got$value, .nearest_centre(fr$z[, 1], fr$z[, 2], cen[, 1], cen[, 2], TRUE),
                     info = lbl)
    # In the layer's own units the margin of the tie holds more than the
    # rounding moves a distance: nothing to say.
    cen <- xy[g$i %% 2 == 0, ]
    got <- .tp_watched(.nearest_centre(xy[, 1], xy[, 2], cen[, 1], cen[, 2], TRUE,
                                       rnd = .tie_rounding(xy)))
    expect_false(got$warned, info = paste(lbl, ".nearest_centre(), layer units"))
  }
  # The nearest other location of the rows FNN settles on its own.
  for (cs in .tp_parted[c(1, 3)]) {
    xy <- .tp_lattice_xy(cs[1], cs[2])
    u  <- xy[order(xy[, 1], xy[, 2]), ]
    local_mocked_bindings(.knn_tied_rows = quiet_rows)
    got <- .tp_watched(.nearest_location_search(u, seq_len(64), TRUE, rnd = .tie_rounding(u)))
    expect_true(got$warned, info = sprintf("spacing %g, offset %g", cs[1], cs[2]))
    expect_identical(got$value, .nearest_location_search(u, seq_len(64), TRUE))
  }
})

test_that("the structure test notes a parted tie on its kd-tree shortcut", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("Matrix")
  # Scattered points, none tied, and one point whose second and third
  # nearest lie 100.1 m east and 100.1 m north of it.  At 1e11 the two are
  # stored apart by three times the tie and a tenth of the band of the
  # rounding; every other row is clear of a tie, so the k = 2 graph takes
  # the shortcut.  At the origin the two tie, and the graph is
  # .build_knn_weights()'.
  base <- withr::with_seed(20261111L, cbind(stats::runif(40, 0, 3000),
                                            stats::runif(40, 0, 3000)))
  base <- base[sqrt((base[, 1] - 1500)^2 + (base[, 2] - 1500)^2) > 400, ]
  P <- c(1500.3, 1500.7)
  trio <- rbind(P, P + c(30, 40), P + c(100.1, 0), P + c(0, 100.1))
  y <- function(xy) sin(xy[, 1] / 700) + cos(xy[, 2] / 500)
  orig <- .knn_pairs
  local_mocked_bindings(.knn_pairs = function(coords, k, use_fnn, rnd = NULL)
    orig(coords, k, use_fnn, rnd = NULL))
  for (off in c(0, 1e11)) {
    xy <- rbind(base, trio) + off
    p  <- nrow(base) + 1L
    d2 <- sort((xy[-p, 1] - xy[p, 1])^2 + (xy[-p, 2] - xy[p, 2])^2)
    nn <- FNN::get.knn(xy, k = 3L)$nn.dist
    # The shortcut's own test at k = 2 holds for every row at 1e11.
    expect_identical(any(nn[, 3] <= nn[, 2] * sqrt(1 + .TIE_TOL)), off == 0)
    got <- .tp_watched(.sac_structure_test(xy, y(xy)))
    expect_identical(got$warned, off > 0, info = sprintf("offset %g", off))
    expect_identical(got$value, .sac_structure_test(xy, y(xy)))
    if (off > 0)
      expect_true(.tie_lost_by_rounding(d2[2], d2[3], off))
    # Without the point, nothing ties.
    expect_false(.tp_watched(.sac_structure_test(xy[-p, ], y(xy[-p, ])))$warned)
  }
})

test_that("make_folds(method = 'knndm') warns once where its k-means step reads a parted tie", {
  g <- expand.grid(i = 0:19, j = 0:15)
  xy <- cbind(100.1 * g$i + 1e11, 100.1 * g$j + 1e11)
  # Ward's method alone reads no search of its own; the k-means step that
  # reduces a large sample to `ward_max` groups does.
  run <- function(wm) .tp_watched(withr::with_seed(1L, .knndm_folds(
    xy, g_target = seq(20, 200, length.out = 50), k = 5, ward_max = wm)))
  expect_false(run(4000L)$warned)
  got <- run(100L)
  expect_true(got$warned)
  expect_identical(got$value, withr::with_seed(1L, .knndm_folds(
    xy, g_target = seq(20, 200, length.out = 50), k = 5, ward_max = 100L)))
  # End to end: 4200 points, past `ward_max`.
  lat <- function(nx, ny, sp, off) {
    g <- expand.grid(i = 0:(nx - 1), j = 0:(ny - 1))
    sf::st_as_sf(data.frame(x = sp * g$i + off, y = sp * g$j + off, z = sin(g$i / 2)),
                 coords = c("x", "y"), crs = 32615)
  }
  p <- lat(70, 60, 100.1, 1e11); pp <- lat(20, 15, 350.35, 1e11 + 50)
  call <- function() suppressMessages(make_folds(p, k = 5, method = "knndm", seed = 1,
                                                 prediction_points = pp))
  got <- .tp_catch(call())
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]),
               "make_folds(): `points_sf` lies so far from the origin of its CRS",
               fixed = TRUE)
  again <- withCallingHandlers(call(), warning = function(cnd) invokeRestart("muffleWarning"))
  expect_identical(got$value, again)
})

# A 20 x 16 lattice at 100.1 m, cells of 2 x 2 points and four folds of
# blocks, as kriging_adequacy() reads them.
.tp_krig <- function(off, sp = 100.1) {
  g <- expand.grid(i = 0:19, j = 0:15)
  z <- withr::with_seed(3L, sin(g$i * sp / 300 + 0.4) * cos(g$j * sp / 260) +
                          stats::rnorm(nrow(g), 0, 0.4))
  fold <- (g$i %/% 5L) %% 2L * 2L + (g$j %/% 4L) %% 2L + 1L
  pts <- sf::st_as_sf(data.frame(x = off + g$i * sp, y = off + g$j * sp, z = z, fold = fold),
                      coords = c("x", "y"), crs = 32615)
  box <- sf::st_bbox(c(xmin = off - sp / 2, ymin = off - sp / 2,
                       xmax = off + 19.5 * sp, ymax = off + 15.5 * sp),
                     crs = sf::st_crs(32615))
  cells <- sf::st_sf(cell_id = 1:80, geometry = sf::st_make_grid(sf::st_as_sfc(box),
                                                                 n = c(10, 8)))
  list(pts = pts, cells = cells, cell = (g$j %/% 2L) * 10L + g$i %/% 2L + 1L, z = z, g = g)
}

test_that("the kriging neighbourhood of a cell notes a parted tie", {
  skip_if_not_installed("sp")
  for (off in c(0, 1e11)) {
    k <- .tp_krig(off)
    # Each cell's centre is 70.8 m from its four points and 158.3 m from the
    # next eight: the sixth nearest ties with seven others.
    got <- .tp_watched(.cell_neighbourhoods(k$cells, k$pts, k$cell, 6L, 2000L,
                                            skip = rep(FALSE, 80)))
    expect_identical(got$warned, off > 0, info = sprintf("offset %g", off))
    expect_identical(got$value, .cell_neighbourhoods(k$cells, k$pts, k$cell, 6L, 2000L,
                                                     skip = rep(FALSE, 80)))
  }
})

test_that("the count that enlarges a kriging cell notes locations parted from its tie", {
  r <- list(r = .ulp_of(1e11), spacing = .ulp_of(1e11), scale = 1)
  a <- 1e4
  line <- a * (1 + .TIE_TOL)
  band <- 8 * r$r * 2 * sqrt(a)
  near <- line + band / 10
  far  <- line + band * 10
  run <- function(d2, nmax) .tp_watched(.tie_count_note(d2, a, nmax, r))$warned
  # Two own locations at the tie, two parted from it: counted, they would
  # enlarge a cell of nmax 3.
  expect_true(run(c(a, a, near, near, 5e4), 3L))
  # Counted or not, the cell is enlarged (nmax 1), or not (nmax 4).
  expect_false(run(c(a, a, near, near, 5e4), 1L))
  expect_false(run(c(a, a, near, near, 5e4), 4L))
  # Locations well past the tie are not parted from it.
  expect_false(run(c(a, a, far, far, 5e4), 3L))
})

test_that("the enlargement count of a kriging cell is read in .cell_neighbourhoods()", {
  skip_if_not_installed("sp")
  # A cell of the lattice at 1e11 whose four points the rounding stores at
  # two distances from its centre: the nearer are its own, the farther
  # another cell's.  Counted with the tie, nmax own-and-tied locations do
  # not enlarge it; counted with the parted ones too, they would.  The
  # neighbourhood's own tie (.tie_pick_note()) is silenced, so the count
  # alone decides.
  k <- .tp_krig(1e11)
  xy <- sf::st_coordinates(k$pts)
  r <- list(r = .ulp_of(1e11))
  found <- NULL
  for (c in seq_len(80)) {
    cen <- sp::coordinates(sf::as_Spatial(sf::st_set_crs(sf::st_geometry(k$cells)[c],
                                                         sf::NA_crs_)))
    own <- which(k$cell == c)
    d2 <- (xy[own, 1] - cen[1, 1])^2 + (xy[own, 2] - cen[1, 2])^2
    lo <- d2 <= min(d2) * (1 + .TIE_TOL)
    if (any(!lo) && all(.tie_lost_by_rounding(rep(min(d2), sum(!lo)), d2[!lo], r = r$r))) {
      found <- list(c = c, near = own[lo], far = own[!lo])
      break
    }
  }
  expect_false(is.null(found))
  cell <- k$cell
  cell[found$far] <- if (found$c > 1L) found$c - 1L else found$c + 1L
  skip <- seq_len(80) != found$c
  local_mocked_bindings(.tie_pick_note = function(d2, k, rnd) invisible(FALSE))
  got <- .tp_watched(.cell_neighbourhoods(k$cells, k$pts, cell, length(found$near), 2000L,
                                          skip = skip))
  expect_true(got$warned)
  expect_false(got$value$enlarged[found$c])
})

test_that("kriging_adequacy() warns once on the parted lattice", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("sp")
  k0 <- .tp_krig(0)
  sac <- suppressWarnings(suppressMessages(estimate_sac_range(k0$pts, "z")))
  run <- function(off) {
    k <- .tp_krig(off)
    asg <- suppressWarnings(assign_features_to_polygons(k$pts, k$cells,
                                                        polygon_id_col = "cell_id",
                                                        quiet = TRUE))
    function() suppressMessages(kriging_adequacy(asg, "z", k$cells, id_col = "cell_id",
                                                 sac = sac, folds = asg$fold, nmax = 6))
  }
  expect_length(.tp_catch(run(0)())$warnings, 0L)
  call <- run(1e11)
  got <- .tp_catch(call())
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]),
               "kriging_adequacy(): `assigned_points_sf` lies so far from the origin of its CRS",
               fixed = TRUE)
  again <- withCallingHandlers(call(), warning = function(cnd) invokeRestart("muffleWarning"))
  expect_identical(got$value, again)
})

# A 10 x 10 lattice for the GWR.
.tp_gwr <- function(sp, off) {
  g <- expand.grid(i = 0:9, j = 0:9)
  sf::st_as_sf(data.frame(x = sp * g$i + off, y = sp * g$j + off,
                          z = sin(g$i / 2) + cos(1.3 * g$j / 2) + 0.1 * sin(7 * g$i * g$j),
                          x1 = g$i / 10 + 0.3 * cos(g$j), x2 = sin(g$i * g$j)),
               coords = c("x", "y"), crs = 32615)
}

test_that("a GWR's boxcar window warns once where the rounding parts its edge", {
  skip_if_not_installed("sp")
  fit <- function(d, kernel = "boxcar")
    suppressMessages(fit_gwr_model(d, "z", "x1", adaptive = TRUE, bandwidth = 12,
                                   kernel = kernel))
  for (cs in list(c(100.1, 1e11), c(0.7, 1e10))) {
    d <- .tp_gwr(cs[1], cs[2])
    lbl <- sprintf("spacing %g, offset %g", cs[1], cs[2])
    got <- .tp_catch(fit(d))
    expect_length(got$warnings, 1L)
    expect_match(conditionMessage(got$warnings[[1]]),
                 "fit_gwr_model(): `data_sf` lies so far from the origin of its CRS",
                 fixed = TRUE, info = lbl)
    expect_identical(got$warnings[[1]]$spacing, .ulp_of(cs[2]), info = lbl)
    again <- withCallingHandlers(fit(d), warning = function(cnd) invokeRestart("muffleWarning"))
    expect_identical(fitted(got$value), fitted(again), info = lbl)
  }
  # The ties hold: at the origin, and on whole metres at 1e12.  A kernel
  # that weighs a point at its edge 0 has no tie to lose.
  expect_length(.tp_catch(fit(.tp_gwr(100.1, 0)))$warnings, 0L)
  expect_length(.tp_catch(fit(.tp_gwr(250, 1e12)))$warnings, 0L)
  expect_length(.tp_catch(fit(.tp_gwr(100.1, 1e11), "bisquare"))$warnings, 0L)
  # predict(), cv_gwr() and gwr_model_selection() say it once, by their
  # own name.
  d <- .tp_gwr(100.1, 1e11)
  f <- withCallingHandlers(fit(d), warning = function(cnd) invokeRestart("muffleWarning"))
  got <- .tp_catch(predict(f, d, se.fit = TRUE))
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]),
               paste0("predict.gwr_fit(): `newdata` and the data of `object` lie so far ",
                      "from the origin of their CRS that their coordinates are rounded"),
               fixed = TRUE)
  got <- .tp_catch(suppressMessages(cv_gwr(d, "z", "x1", adaptive = TRUE, bandwidth = 12,
                                           kernel = "boxcar", k = 3)))
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]), "^cv_gwr\\(\\): `data_sf` lies")
  got <- .tp_catch(suppressMessages(gwr_model_selection(d, "z", c("x1", "x2"),
                                                        bandwidth = 12, kernel = "boxcar")))
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]), "^gwr_model_selection\\(\\): `data_sf` lies")
})

test_that("a call that says its coordinates are rounded does not say it twice", {
  .tp_close_watch()
  # 120 points over 4000 x 3000 m at 1e15: make_folds() says that the
  # rounding (0.125) is more than a millionth of the extent, and its range
  # estimate reads ties that rounding parts.  One warning of the class.
  d <- withr::with_seed(3L, data.frame(x = stats::runif(120, 0, 4000) + 1e15,
                                       y = stats::runif(120, 0, 3000) + 1e15))
  d$z <- withr::with_seed(4L, sin(d$x / 700) + cos(d$y / 500) + stats::rnorm(120, 0, 0.2))
  p <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32615)
  got <- .tp_catch(suppressMessages(make_folds(p, k = 4, method = "block_kfold",
                                               auto_range = TRUE, response_var = "z", seed = 1)))
  expect_length(got$warnings, 1L)
  expect_match(conditionMessage(got$warnings[[1]]), "more than a millionth of its extent",
               fixed = TRUE)
  # The watch on its own: the extent warning, then a parted tie.
  xy <- .tp_lattice_xy(100.1, 1e11)
  n <- 0L
  withCallingHandlers({
    end <- .tie_watch("f", "xy")
    .warn_xy_imprecise(c(xmin = 1e15, ymin = 1e15, xmax = 1e15 + 10, ymax = 1e15), "f", "xy")
    .knn_pairs(xy, 6L, use_fnn = FALSE)
    expect_false(end(NULL))
  }, spatialkit_xy_imprecise = function(cnd) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  expect_identical(n, 1L)
  # The next call starts afresh.
  expect_true(.tp_watched(.knn_pairs(xy, 6L, use_fnn = FALSE))$warned)
})

test_that("the message names several things in the plural, and a CRS given as such", {
  .tp_close_watch()
  xy <- .tp_lattice_xy(100.1, 1e11)
  msg <- NULL
  withCallingHandlers({
    end <- .tie_watch("select_resolution", "the points of `profile`", plural = TRUE)
    .knn_pairs(xy, 6L, use_fnn = FALSE)
    end(NULL, sf::st_crs(32615))
  }, spatialkit_xy_imprecise = function(cnd) {
    msg <<- conditionMessage(cnd)
    invokeRestart("muffleWarning")
  })
  expect_match(msg, paste0("select_resolution(): the points of `profile` lie so far from the ",
                           "origin of their CRS that their coordinates are rounded to"),
               fixed = TRUE)
})

test_that("ties a forked worker notes reach the watch of the process that forked it", {
  .tp_close_watch()
  skip_on_os("windows")
  xy <- .tp_lattice_xy(100.1, 1e11)
  warned <- FALSE
  withCallingHandlers({
    end <- .tie_watch("f", "xy")
    job <- parallel::mcparallel({
      .knn_pairs(xy, 6L, use_fnn = FALSE)
      .tie_carry()
    })
    carried <- parallel::mccollect(job)[[1]]
    expect_true(carried$spacing > 0)
    .tie_merge(carried)
    end(NULL)
  }, spatialkit_xy_imprecise = function(cnd) {
    warned <<- TRUE
    invokeRestart("muffleWarning")
  })
  expect_true(warned)
  expect_null(.tie_carry())
})
