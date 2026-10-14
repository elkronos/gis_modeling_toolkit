# tests/testthat/test-core-pipeline.R
# ---------------------------------------------------------------------------
# The path from a profile to a summary, through exported functions only:
# resolution_profile(), select_resolution(), get_voronoi_seeds(),
# build_tessellation(), assign_features_to_polygons(), summarize_by_cell().
# The cells built from a selection are the cells the profile scored: every
# row lies in the cell of its nearest scored centre, and the leave-one-out
# error computed on the built and assigned cells is the profile's at the
# chosen level.  On projected and lon/lat input, for both kinds of
# partition, for a count set by a size and one chosen by the leave-one-out
# minimum, and for centres fitted to a subsample.
# Also here: what the check reports on gridded locations, where a row can
# lie exactly on the boundary between two cells, and for features that are
# not points.
#
# Every profile names `partition`, `nstart`, `estimate`, `obs_per_cell` and
# `sac`, every selection its `criterion`, every summary its `estimate`.
# ---------------------------------------------------------------------------

.pp_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.pp_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.pp_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)

.pp_names <- c("version", "rule", "estimate", "cells", "obs_per_cell", "obs_per_cell_asked",
               "response", "n", "n_points", "n_fit", "seed", "features", "scored",
               "partition", "crs", "scored_cells")
.pp_sizes <- c(18, 12, 9, 6, 4.5, 3, 2)

# The layer: seed 7001; `u`, then `v`, 3,000 values each, uniform on 0 to
# 1e5; then the response, sin(u / 2e4) + cos(v / 1.5e4) plus normal noise of
# standard deviation 1.5; the points at (5e5 + u, 5e6 + v) in EPSG:32617.
.pp_store <- new.env(parent = emptyenv())
.pp_layer <- function() {
  if (is.null(.pp_store$pts)) {
    d <- local({
      .pp_seed(7001)
      u <- stats::runif(3000, 0, 1e5)
      v <- stats::runif(3000, 0, 1e5)
      z <- sin(u / 2e4) + cos(v / 1.5e4) + stats::rnorm(3000, 0, 1.5)
      data.frame(x = 5e5 + u, y = 5e6 + v, z = z)
    })
    .pp_store$pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
  }
  .pp_store$pts
}

# One path, by the calls a user writes, with what every path has to show.
.pp_path <- function(pts, prof, criterion, estimate, n_fit = 3000L) {
  sel   <- select_resolution(prof, criterion = criterion)
  seeds <- .pp_catch(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  # The one warning the seeding may raise here: a minimum at an end of the sizes.
  expect_true(all(vapply(seeds$warnings, inherits, logical(1), "spatialkit_bound_chosen")))
  seeds <- seeds$value
  tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                       quiet = TRUE)
  out   <- .pp_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                       estimate = estimate))

  # ---- the record, from the selection to the summary -------------------------
  rec <- sel$count_rule
  expect_identical(names(rec), .pp_names)
  expect_identical(rec$rule, criterion)
  expect_identical(rec$estimate, estimate)
  expect_identical(rec$cells, as.integer(sel$best))
  expect_identical(rec$response, "z")
  expect_identical(rec$n, 3000L)
  expect_identical(rec$n_points, 3000L)
  expect_identical(rec$n_fit, as.integer(n_fit))
  expect_identical(rec$seed, 123)
  expect_identical(rec$features, "points")
  expect_identical(rec$scored$rows, "all")
  expect_identical(rec$scored$n, 3000L)
  expect_identical(rec$scored_cells, NA)
  on_seeds <- rec; on_seeds$scored_cells <- TRUE
  expect_identical(attr(sel$seeds, "count_rule"), on_seeds)
  expect_identical(attr(seeds, "count_rule"), on_seeds)
  expect_identical(tess$params$count_rule, on_seeds)
  expect_identical(attr(tess$cells, "count_rule"), on_seeds)
  expect_identical(attr(asg, "count_rule"), on_seeds)
  got <- attr(out$value, "count_rule")
  expect_identical(got$check, list(estimate = "match", cells = "match"))
  got$check <- NULL
  expect_identical(got, on_seeds)
  expect_length(.pp_of(out, "spatialkit_cells_differ"), 0L)
  expect_length(.pp_of(out, "spatialkit_estimate_mismatch"), 0L)

  # ---- every row in the cell of its nearest scored centre ---------------------
  expect_identical(nrow(seeds), as.integer(sel$best))
  expect_identical(nrow(asg), 3000L)
  cen <- unname(sf::st_coordinates(sel$seeds)[, 1:2])
  expect_identical(sf::st_crs(sel$seeds), rec$crs)
  xy  <- unname(sf::st_coordinates(sf::st_transform(pts, rec$crs))[, 1:2])
  nearest <- .nearest_centre(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L], .has_backend("FNN"))
  expect_false(anyNA(tess$index))
  expect_identical(asg$cell_id, tess$index[nearest])

  # ---- the summary's leave-one-out error is the profile's at that level -------
  at_level <- prof$loo[match(sel$best, prof$levels)]
  expect_identical(rec$scored$loo, at_level)
  fit <- .pool_column(asg$z, asg$cell_id, .pooled_xy(asg, rec$crs)$xy, NULL, want = estimate,
                      use_fnn = .has_backend("FNN"))$fit
  found <- if (identical(estimate, "pooled")) fit$loo_pooled else fit$loo_mean
  expect_true(is.finite(at_level))
  expect_lt(abs(found - at_level) / at_level, 1e-10)
  expect_identical(fit$N, 3000L)
  expect_identical(fit$K, rec$scored$cells)
  if (identical(estimate, "pooled")) {
    own <- attr(out$value, "pooled")$columns$resp_z
    expect_lt(abs(own$loo - at_level) / at_level, 1e-10)
    expect_identical(attr(out$value, "pooled")$crs, rec$crs)
  }
  invisible(list(sel = sel, seeds = seeds, tess = tess, asg = asg, out = out$value))
}


# ---- End to end -------------------------------------------------------------------

test_that("the cells built from a selection reproduce the profile's score, on projected and lon/lat input", {
  pts <- .pp_layer()
  layers <- list(projected = pts, `lon/lat` = sf::st_transform(pts, 4326))
  kinds  <- list(list(partition = "density", nstart = 3), list(partition = "kmeans", nstart = 5))
  paths <- 0L
  for (nm in names(layers)) for (kind in kinds) {
    layer <- layers[[nm]]
    # By a size, for pooled values.
    pooled <- resolution_profile(layer, "z", partition = kind$partition, nstart = kind$nstart,
                                 sac = "auto", seed = 123, estimate = "pooled",
                                 obs_per_cell = 4.5)
    P <- .pp_path(layer, pooled, "size", "pooled")
    expect_identical(P$sel$best, 667L)
    expect_identical(P$sel$count_rule$obs_per_cell_asked, 4.5)
    expect_identical(P$sel$count_rule$partition$kind, kind$partition)
    expect_identical(P$sel$count_rule$partition$nstart, as.integer(kind$nstart))
    # By the lowest leave-one-out error, for plain means.
    plain <- suppressWarnings(
      resolution_profile(layer, "z", partition = kind$partition, nstart = kind$nstart,
                         sac = "auto", seed = 123, estimate = "mean",
                         obs_per_cell = .pp_sizes))
    M <- .pp_path(layer, plain, "loo", "mean")
    expect_true(M$sel$best %in% plain$levels)
    # The seeds come back in the CRS of the points they are for, and the
    # record names the projected CRS the profile worked in.
    expect_identical(sf::st_crs(P$seeds), sf::st_crs(layer))
    expect_false(sf::st_is_longlat(P$sel$count_rule$crs))
    paths <- paths + 2L
  }
  expect_identical(paths, 8L)
})

test_that("centres fitted to a subsample are checked on every row of the layer", {
  skip_if_not_installed("gstat")
  pts <- .pp_layer()
  prof <- suppressWarnings(
    resolution_profile(pts, "z", estimate = "pooled", obs_per_cell = NULL, sample_n = 1500,
                       sac = NULL, partition = "density", nstart = 3, seed = 123))
  expect_identical(attr(prof, "bounds")$n_sample, 1500L)
  S <- .pp_path(pts, prof, "loo", "pooled", n_fit = 1500L)
  rec <- attr(S$out, "count_rule")
  expect_identical(rec$n_fit, 1500L)
  expect_identical(rec$n_points, 3000L)
  expect_match(gsub("\\s+", " ", paste(utils::capture.output(print(rec)), collapse = " ")),
               "fitted to a subsample of 1500 of the 3000 points", fixed = TRUE)
})


# ---- Points on a lattice ----------------------------------------------------------

# 60 x 50 locations at one spacing from one corner, in EPSG:32617; the
# response of seed 9001, drawn after the locations are laid.
.pp_lattice <- function(step, x0, y0) {
  g <- expand.grid(i = seq_len(60L) - 1L, j = seq_len(50L) - 1L)
  z <- local({
    .pp_seed(9001)
    sin(g$i / 12) + cos(g$j / 9) + stats::rnorm(nrow(g), 0, 1.5)
  })
  sf::st_as_sf(data.frame(x = x0 + step * g$i, y = y0 + step * g$j, z = z,
                          row = seq_len(nrow(g))),
               coords = c("x", "y"), crs = 32617)
}

test_that("on gridded locations the check says match exactly when no row was placed in another cell", {
  lattices <- list(`whole metres` = .pp_lattice(100, 5e5, 5e6),
                   `spacing 0.7`  = .pp_lattice(0.7, 5e5 + 0.1, 5e6 + 0.3))
  kinds <- list(list(partition = "kmeans", nstart = 5), list(partition = "density", nstart = 3))
  moved_in <- integer(0)
  for (nm in names(lattices)) for (kind in kinds) {
    pts  <- lattices[[nm]]
    info <- paste(nm, kind$partition)
    prof <- suppressWarnings(
      resolution_profile(pts, "z", partition = kind$partition, nstart = kind$nstart,
                         sac = "auto", seed = 123, estimate = "pooled", obs_per_cell = 4.5))
    sel   <- select_resolution(prof, criterion = "size")
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
    tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    # The profile's placement: each row to its nearest scored centre, a tie
    # to the lower-numbered one.
    cen <- unname(sf::st_coordinates(sel$seeds)[, 1:2])
    xy  <- unname(sf::st_coordinates(pts)[, 1:2])
    scored_in <- tess$index[.nearest_centre(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L], .has_backend("FNN"))]
    # The path's placement: the polygon that holds the row, a row on a shared
    # edge to the smaller polygon.
    asg <- suppressWarnings(assign_features_to_polygons(pts, tess$cells,
                                                        polygon_id_col = "cell_id", quiet = TRUE))
    expect_identical(formals(assign_features_to_polygons)$tie_break,
                     quote(c("smallest_area", "first")))
    placed_in  <- asg$cell_id[match(pts$row, asg$row)]
    unassigned <- sum(is.na(placed_in))
    moved      <- sum(placed_in != scored_in, na.rm = TRUE)
    moved_in[info] <- moved
    # A row placed in another cell than the profile's is a row on a boundary:
    # the centre of the cell it was placed in is as near as its nearest
    # centre, exactly or to the rounding of the coordinates.  Any other row
    # in another cell would be a row the path misplaced.
    elsewhere <- which(!is.na(placed_in) & placed_in != scored_in)
    if (length(elsewhere)) {
      d2 <- function(rows, centre) (xy[rows, 1L] - cen[centre, 1L])^2 +
                                   (xy[rows, 2L] - cen[centre, 2L])^2
      nearest <- match(scored_in[elsewhere], tess$index)
      placed  <- match(placed_in[elsewhere], tess$index)
      gap <- abs(d2(elsewhere, placed) - d2(elsewhere, nearest)) / d2(elsewhere, nearest)
      expect_true(all(gap < 1e-6), info = info)
    }

    r <- .pp_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                     estimate = "pooled"))
    check <- attr(r$value, "count_rule")$check$cells
    w <- .pp_of(r, "spatialkit_cells_differ")
    own <- attr(r$value, "pooled")$columns$resp_z
    at_level <- prof$loo[match(sel$best, prof$levels)]
    if (unassigned == 0L && moved == 0L) {
      expect_identical(check, "match", info = info)
      expect_length(w, 0L)
      expect_lt(abs(own$loo - at_level) / at_level, 1e-10)
      # In any order of the rows: the profile holds them in coordinate
      # order, the summary as they come.
      for (seed in c(1L, 2L)) {
        shuffled <- asg[local({ .pp_seed(seed); sample.int(nrow(asg)) }), ]
        attr(shuffled, "count_rule") <- attr(asg, "count_rule")
        rs <- .pp_catch(summarize_by_cell(shuffled, "z", id_col = "cell_id",
                                          cells_sf = tess$cells, estimate = "pooled"))
        expect_identical(attr(rs$value, "count_rule")$check$cells, "match", info = info)
        loo <- attr(rs$value, "pooled")$columns$resp_z$loo
        expect_lt(abs(loo - at_level) / at_level, 1e-10)
      }
    } else {
      expect_identical(check, "differ", info = info)
      expect_length(w, 1L)
      # A row in no cell is a row fewer; rows in other cells leave the count,
      # and change the occupied cells or, with those unchanged, the error.
      kind_of <- if (unassigned > 0L) "n"
                 else if (!identical(own$cells, sel$count_rule$scored$cells)) "cells"
                 else "loo"
      expect_identical(w[[1L]]$what, kind_of, info = info)
      expect_identical(w[[1L]]$found$n, 3000L - unassigned, info = info)
      if (identical(kind_of, "loo"))
        expect_match(conditionMessage(w[[1L]]),
                     "rows exactly on a boundary between two cells (gridded locations)",
                     fixed = TRUE)
    }
  }
  # Which lattice has rows on a boundary rests on exact ties, which another
  # GEOS or another arithmetic can settle the other way, so the counts are
  # held only on the platform the frozen fixture was recorded on.  There the
  # k-means cells put no row of either lattice on a boundary and the density
  # cells put 5 and 7 of 3,000 there, so both branches above are taken.
  expect_identical(names(moved_in), c("whole metres kmeans", "whole metres density",
                                      "spacing 0.7 kmeans", "spacing 0.7 density"))
  recorded <- readRDS(test_path("fixtures", "core-frozen.rds"))$info
  if (core_frozen_same_platform(core_frozen_platform(), recorded))
    expect_identical(unname(moved_in), c(0L, 5L, 0L, 7L))
})


# ---- Features that are not points -------------------------------------------------

test_that("polygon features are scored at one point and assigned by overlap, and the cells are not compared", {
  # 4,800 squares centred on locations of seed 5105, of side a third of the
  # mean cell width at 4.5 rows a cell.  A square's point is its centre, and
  # the cell that holds the centre holds most of the square unless the
  # square lies across a corner of the cells, which about one square in a
  # thousand does: 10 of these, and 4, 5 and 4 of the layers of seeds 5102
  # to 5104.
  d <- local({
    .pp_seed(5105)
    x <- stats::runif(4800, 0, 1e5); y <- stats::runif(4800, 0, 1e5)
    data.frame(x = 5e5 + x, y = 5e6 + y,
               z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(4800, 0, 1.5))
  })
  side <- sqrt(1e10 / 1067) / 3
  squares <- lapply(seq_len(nrow(d)), function(i) sf::st_polygon(list(cbind(
    d$x[i] + side * c(-0.5, 0.5, 0.5, -0.5, -0.5), d$y[i] + side * c(-0.5, -0.5, 0.5, 0.5, -0.5)))))
  polys <- sf::st_sf(z = d$z, row = seq_len(nrow(d)), geometry = sf::st_sfc(squares, crs = 32617))

  prof <- resolution_profile(polys, "z", partition = "density", nstart = 3, sac = "auto",
                             seed = 123, estimate = "pooled", obs_per_cell = 4.5)
  expect_identical(attr(prof, "estimate")$features, "other")
  sel <- select_resolution(prof, criterion = "size")
  expect_identical(sel$best, 1067L)
  expect_identical(sel$count_rule$features, "other")
  at_points <- coerce_to_points(polys, "auto")
  seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = at_points)
  tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  expect_identical(attr(tess$cells, "count_rule")$scored_cells, TRUE)

  # By largest overlap.  Some squares straddle a boundary, and a few have
  # most of their area in another cell than their point's.
  asg <- assign_features_to_polygons(polys, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  expect_identical(nrow(asg), 4800L)
  # The assigned features carry the record of the cells, points or not.
  expect_true(all(sf::st_geometry_type(asg) == "POLYGON"))
  expect_identical(attr(asg, "count_rule"), attr(tess$cells, "count_rule"))
  expect_gt(sum(lengths(sf::st_intersects(polys, tess$cells)) > 1L), 1000L)
  cen <- unname(sf::st_coordinates(sel$seeds)[, 1:2])
  xy  <- unname(sf::st_coordinates(at_points)[, 1:2])
  scored_in <- tess$index[.nearest_centre(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L], .has_backend("FNN"))]
  elsewhere <- sum(asg$cell_id[match(polys$row, asg$row)] != scored_in)
  expect_gt(elsewhere, 0L)

  r <- .pp_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                   estimate = "pooled"))
  expect_length(.pp_of(r, "spatialkit_cells_differ"), 0L)
  expect_identical(attr(r$value, "count_rule")$check,
                   list(estimate = "match", cells = "not checked: the features are not points"))
  # The summary's own error describes the cells as assigned, which are not
  # the cells as scored wherever a square was assigned elsewhere.
  own <- attr(r$value, "pooled")$columns$resp_z$loo
  expect_true(is.finite(own))
  if (elsewhere > 0L)
    expect_false(isTRUE(all.equal(own, sel$count_rule$scored$loo, tolerance = 1e-8)))

  # The same features assigned by their points are the rows the profile scored.
  by_point <- assign_features_to_polygons(at_points, tess$cells, polygon_id_col = "cell_id",
                                          quiet = TRUE)
  expect_identical(by_point$cell_id[match(polys$row, by_point$row)], scored_in)
  r <- .pp_catch(summarize_by_cell(by_point, "z", id_col = "cell_id", cells_sf = tess$cells,
                                   estimate = "pooled"))
  expect_identical(attr(r$value, "count_rule")$check, list(estimate = "match", cells = "match"))
  expect_length(.pp_of(r, "spatialkit_cells_differ"), 0L)
  at_level <- sel$count_rule$scored$loo
  expect_lt(abs(attr(r$value, "pooled")$columns$resp_z$loo - at_level) / at_level, 1e-10)
})
