# tests/testthat/test-core-gridded-ties.R
# ---------------------------------------------------------------------------
# Points on the edge between two cells, on gridded locations.
# resolution_profile() gives each scored row to its nearest centre in the CRS
# it works in, a row the same distance from two centres to the lower-numbered
# one.  On cells built from the centres it scored (their record of the count
# says scored_cells = TRUE and holds the centres), assign_features_to_polygons()
# gives a point on the edge between two of those cells the cell the profile
# gave it, so the summary's check reads `match` on lattices too, and
# cell_values_at() places such a location in the same cell, with seeds and
# without.  Everything else is placed by the polygons alone: polygons a
# caller brings, cells whose record says scored_cells = FALSE, features that
# are not points, a point farther from the edge than rounding.  A real
# difference between the cells or rows and the ones scored is still reported.
#
# Every profile names `partition`, `nstart`, `estimate`, `obs_per_cell`,
# `sac` and `seed`, every selection its `criterion`, every summary its
# `estimate`.
# ---------------------------------------------------------------------------

# ---- Layers and helpers ------------------------------------------------------

.gt_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

.gt_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.gt_of <- function(r, class) Filter(function(w) inherits(w, class), r$warnings)

.gt_store <- new.env(parent = emptyenv())
.gt_keep <- function(key, expr) {
  if (is.null(.gt_store[[key]])) .gt_store[[key]] <- expr
  .gt_store[[key]]
}

# A lattice: `nx` x `ny` locations at spacing `sp` from (x0, y0), and a
# smooth field plus normal noise of standard deviation 0.4 drawn under seed
# 9001.  `row` numbers the locations.  Sixteen of them, 20 x 15 to 50 x 40 at
# spacings of 1, 30, 250 and 1000, are the layers of the first test.
.gt_lattice <- function(nx, ny, sp, x0 = 4e5, y0 = 4.4e6, crs = 32615) {
  g <- expand.grid(x = x0 + (seq_len(nx) - 1) * sp, y = y0 + (seq_len(ny) - 1) * sp)
  z <- local({
    .gt_seed(9001)
    sin(g$x / (sp * 6)) + cos(g$y / (sp * 5)) + stats::rnorm(nrow(g), 0, 0.4)
  })
  sf::st_as_sf(data.frame(x = g$x, y = g$y, z = z, row = seq_len(nrow(g))),
               coords = c("x", "y"), crs = crs)
}
.gt_dims    <- list(c(20L, 15L), c(30L, 20L), c(40L, 30L), c(50L, 40L))
.gt_spacing <- c(1, 30, 250, 1000)

# The path from a profile to a summary, with its arguments written out: the
# profile, the count by size, the scored centres as seeds, their Voronoi
# cells, the assignment and the summary.
.gt_path <- function(pts, kind, obs_per_cell = 4.5, size = NULL, boundary = NULL, crs = NULL) {
  prof <- suppressWarnings(
    resolution_profile(pts, "z", partition = kind, nstart = 3, sac = "auto", seed = 123,
                       estimate = "pooled", obs_per_cell = obs_per_cell))
  sel   <- select_resolution(prof, criterion = "size", obs_per_cell = size)
  seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts,
                                              boundary = boundary))
  tess  <- build_tessellation(seeds, boundary = boundary, method = "voronoi", crs = crs,
                              quiet = TRUE)
  asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id",
                                       quiet = TRUE)
  out   <- .gt_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                       estimate = "pooled"))
  seeds$cell_id <- tess$index
  list(prof = prof, sel = sel, seeds = seeds, tess = tess, asg = asg, out = out$value,
       warnings = out$warnings)
}
.gt_lattice_path <- function(d, sp, kind)
  .gt_keep(sprintf("%dx%d-%g-%s", d[1L], d[2L], sp, kind),
           .gt_path(.gt_lattice(d[1L], d[2L], sp), kind))

.gt_bare <- function(x) { attr(x, "count_rule") <- NULL; x }
.gt_cell_of <- function(asg, pts, id = "cell_id") asg[[id]][match(pts$row, asg$row)]

# The profile's rule, written out: each row to the centre the least squared
# distance away in the record's CRS, of two as near the lower-numbered one
# (.nearest_exact() replaces the best so far only when a centre is strictly
# nearer), and the cell build_tessellation() built around that centre.
.gt_profile_rule <- function(P, pts) {
  rec <- P$sel$count_rule
  cen <- rec$scored$centres
  xy  <- unname(sf::st_coordinates(sf::st_transform(pts, rec$crs))[, 1:2])
  list(cell = P$tess$index[.nearest_exact(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L])],
       xy = xy, cen = cen)
}

# What the polygons alone give: the same cells without the record.
.gt_polygon_rule <- function(P, pts)
  .gt_cell_of(assign_features_to_polygons(pts, .gt_bare(P$tess$cells),
                                          polygon_id_col = "cell_id", quiet = TRUE), pts)

# Everything a lattice's path has to show.  Returns how many rows the
# polygons alone would have placed in another cell than the profile's.
.gt_expect_scored <- function(P, pts, info) {
  rule   <- .gt_profile_rule(P, pts)
  placed <- .gt_cell_of(P$asg, pts)
  # Every row in the cell the profile scored it in.
  expect_identical(nrow(P$asg), nrow(pts), info = info)
  expect_false(anyNA(placed), info = info)
  expect_identical(placed, rule$cell, info = info)
  # The summary's check reads match, and its error is the profile's.
  expect_identical(attr(P$out, "count_rule")$check, list(estimate = "match", cells = "match"),
                   info = info)
  expect_length(.gt_of(P, "spatialkit_cells_differ"), 0L)
  at_level <- P$sel$count_rule$scored$loo
  own <- attr(P$out, "pooled")$columns$resp_z
  expect_lt(abs(own$loo - at_level) / at_level, 1e-10)
  # cell_values_at() at the rows themselves: their own cells, with seeds and without.
  for (seeds in list(P$seeds, NULL)) {
    at <- .gt_catch(cell_values_at(pts, P$out, P$asg, "z", id_col = "cell_id", seeds = seeds))
    expect_length(at$warnings, 0L)
    expect_identical(at$value$cell_id, placed, info = paste(info, if (is.null(seeds)) "polygons"))
    expect_identical(at$value$filled, rep(FALSE, nrow(pts)), info = info)
  }
  # A row the polygons alone place in another cell is a row on an edge: the
  # centre of that cell is as near as the profile's, exactly or to the
  # rounding of the coordinates.
  old <- .gt_polygon_rule(P, pts)
  expect_false(anyNA(old), info = info)
  elsewhere <- which(old != rule$cell)
  if (length(elsewhere)) {
    d2 <- function(rows, centre) (rule$xy[rows, 1L] - rule$cen[centre, 1L])^2 +
                                 (rule$xy[rows, 2L] - rule$cen[centre, 2L])^2
    nearest <- match(rule$cell[elsewhere], P$tess$index)
    other   <- match(old[elsewhere], P$tess$index)
    gap <- abs(d2(elsewhere, other) - d2(elsewhere, nearest)) / d2(elsewhere, nearest)
    expect_true(all(gap < 1e-6), info = info)
  }
  length(elsewhere)
}

.gt_recorded <- function() readRDS(test_path("fixtures", "core-frozen.rds"))$info
.gt_same_platform <- function() core_frozen_same_platform(core_frozen_platform(), .gt_recorded())


# ---- Sixteen lattices --------------------------------------------------------------

test_that("on lattices at spacings of 1, 30, 250 and 1000 every row is in the cell the profile scored it in, under either kind", {
  moved <- list(density = integer(0), kmeans = integer(0))
  for (kind in c("density", "kmeans")) for (d in .gt_dims) for (sp in .gt_spacing) {
    info <- sprintf("%d x %d at %g, %s", d[1L], d[2L], sp, kind)
    pts <- .gt_lattice(d[1L], d[2L], sp)
    P <- .gt_lattice_path(d, sp, kind)
    expect_identical(P$sel$count_rule$partition$kind, kind)
    expect_identical(P$sel$best, as.integer(round(nrow(pts) / 4.5)), info = info)
    moved[[kind]] <- c(moved[[kind]], .gt_expect_scored(P, pts, info))
  }
  expect_length(moved$density, 16L)
  expect_length(moved$kmeans, 16L)
  # Which rows lie on an edge rests on exact ties, which another GEOS or
  # another arithmetic can settle the other way, so the counts are held only
  # on the platform the frozen fixture was recorded on.  There the polygons
  # alone would place rows of 11 of the 16 lattices in another cell than the
  # profile's under the density kind, and none under k-means, so on those 11
  # the profile's rule decides the rows on an edge.
  if (.gt_same_platform()) {
    expect_identical(moved$density, c(0L, 1L, 0L, 0L, 2L, 3L, 0L, 2L,
                                      3L, 1L, 0L, 1L, 1L, 3L, 3L, 3L))
    expect_identical(moved$kmeans, rep(0L, 16L))
  }
})


# ---- Longitude and latitude, two units, two CRSs ----------------------------------

test_that("on a lattice in longitude and latitude the rows are placed as the profile placed them in its own CRS", {
  layers <- list(
    `30 m in EPSG:4326`  = sf::st_transform(.gt_lattice(30L, 20L, 30), 4326),
    `250 m in EPSG:4326` = sf::st_transform(.gt_lattice(50L, 40L, 250), 4326),
    `degrees` = local({
      g <- expand.grid(i = 0:29, j = 0:19)
      z <- local({ .gt_seed(9001); sin(g$i / 6) + cos(g$j / 5) + stats::rnorm(nrow(g), 0, 0.4) })
      sf::st_as_sf(data.frame(x = -93.6 + 0.003 * g$i, y = 41.1 + 0.003 * g$j, z = z,
                              row = seq_len(nrow(g))), coords = c("x", "y"), crs = 4326)
    }))
  for (nm in names(layers)) for (kind in c("density", "kmeans")) {
    pts <- layers[[nm]]
    info <- paste(nm, kind)
    P <- .gt_path(pts, kind)
    rec <- P$sel$count_rule
    # The profile measured in a projected CRS, and the seeds came back in
    # degrees; the cells are drawn in the profile's CRS.
    expect_false(sf::st_is_longlat(rec$crs), info = info)
    expect_true(sf::st_is_longlat(P$seeds), info = info)
    expect_true(sf::st_crs(P$tess$cells) == rec$crs, info = info)
    .gt_expect_scored(P, pts, info)
    # The seeds, put back in the profile's CRS, are its centres to rounding, so
    # cell_values_at() searches the centres as the record holds them.
    seed_xy <- unname(sf::st_coordinates(sf::st_transform(P$seeds, rec$crs))[, 1:2])
    expect_identical(.values_at_centres(attr(P$out, "count_rule"), seed_xy, rec$crs),
                     list(centres = rec$scored$centres, seed = seq_len(nrow(seed_xy))),
                     info = info)
  }
})

test_that("the same lattices in two units and two CRSs are placed as the profile placed them in each", {
  base <- list(`30 m` = .gt_lattice(30L, 20L, 30, x0 = 4.5e5, y0 = 4.55e6, crs = 26915),
               `250 m` = .gt_lattice(40L, 30L, 250, x0 = 4.5e5, y0 = 4.55e6, crs = 26915))
  crss <- list(`UTM 15N, metres` = sf::st_crs(26915),
               `UTM 15N, US feet` = sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs"),
               `Iowa South, metres` = sf::st_crs(26976),
               `Iowa South, US feet` = sf::st_crs(3418))
  for (b in names(base)) for (cn in names(crss)) for (kind in c("density", "kmeans")) {
    pts <- base[[b]]
    if (!isTRUE(sf::st_crs(pts) == crss[[cn]])) pts <- sf::st_transform(pts, crss[[cn]])
    info <- paste(b, cn, kind)
    P <- .gt_path(pts, kind)
    # The profile worked in the layer's own CRS, in its units.
    expect_true(P$sel$count_rule$crs == crss[[cn]], info = info)
    expect_true(sf::st_crs(P$tess$cells) == crss[[cn]], info = info)
    .gt_expect_scored(P, pts, info)
  }
})

# Midpoints between each centre and its three nearest, in the record's CRS:
# each is the same distance from two centres.
.gt_midpoints <- function(rec) {
  cen <- rec$scored$centres
  m <- do.call(rbind, lapply(seq_len(nrow(cen)), function(i) {
    d <- (cen[, 1L] - cen[i, 1L])^2 + (cen[, 2L] - cen[i, 2L])^2
    j <- order(d)[2:4]
    cbind((cen[i, 1L] + cen[j, 1L]) / 2, (cen[i, 2L] + cen[j, 2L]) / 2)
  }))
  sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L], row = seq_len(nrow(m))),
               coords = c("x", "y"), crs = rec$crs)
}
# At those midpoints: the cell the profile's rule gives, the one the
# assignment gives, and the one cell_values_at() gives with the seeds in
# their own order, reversed, sorted by cell and shuffled.
.gt_expect_midpoints <- function(P, info) {
  rec  <- P$sel$count_rule
  mids <- .gt_midpoints(rec)
  want <- .gt_profile_rule(P, mids)$cell
  expect_false(anyNA(want), info = info)
  placed <- .gt_cell_of(assign_features_to_polygons(mids, P$tess$cells, polygon_id_col = "cell_id",
                                                    keep_unassigned = TRUE, quiet = TRUE), mids)
  expect_identical(placed, want, info = info)
  n <- nrow(P$seeds)
  orders <- list(seq_len(n), rev(seq_len(n)), order(P$seeds$cell_id),
                 local({ .gt_seed(4); sample.int(n) }))
  for (o in orders) {
    at <- cell_values_at(mids, P$out, P$asg, "z", id_col = "cell_id", seeds = P$seeds[o, ])
    expect_identical(at$cell_id, want, info = info)
  }
}

test_that("a lattice in degrees wider than a UTM zone is drawn and placed in the profile's equal-area CRS", {
  # 25 x 11 locations 2.4 degrees apart over the conterminous United States:
  # the profile measures them in Lambert azimuthal equal-area and the seeds
  # come back in degrees. Put back in that CRS they lie off the centres by
  # what the round trip through degrees leaves, which depends on the PROJ
  # build: about five times the margin of an edge on macOS, far less on
  # Windows. Seeds moved by five times the margin are matched as well, so
  # that case is tested on every platform.
  g <- expand.grid(i = 0:24, j = 0:10)
  z <- local({ .gt_seed(9001); sin(g$i / 6) + cos(g$j / 5) + stats::rnorm(nrow(g), 0, 0.4) })
  pts <- sf::st_as_sf(data.frame(x = -124 + 2.4 * g$i, y = 25 + 2.4 * g$j, z = z,
                                 row = seq_len(nrow(g))), coords = c("x", "y"), crs = 4326)
  for (kind in c("density", "kmeans")) {
    P <- suppressWarnings(.gt_path(pts, kind))
    rec <- P$sel$count_rule
    cen <- rec$scored$centres
    expect_match(rec$crs$proj4string, "+proj=laea", fixed = TRUE, info = kind)
    expect_true(sf::st_is_longlat(P$seeds), info = kind)
    # The cells are drawn in the profile's CRS, around the centres as the
    # record holds them: the same cells as seeds placed at the centres there.
    expect_true(sf::st_crs(P$tess$cells) == rec$crs, info = kind)
    at_centres <- sf::st_as_sf(data.frame(x = cen[, 1L], y = cen[, 2L]), coords = c("x", "y"),
                               crs = rec$crs)
    attr(at_centres, "cover") <- attr(P$seeds, "cover")
    at_centres <- .count_rule_put(at_centres, attr(P$seeds, "count_rule"))
    own <- build_tessellation(at_centres, method = "voronoi", quiet = TRUE)
    expect_identical(sf::st_geometry(P$tess$cells), sf::st_geometry(own$cells), info = kind)
    seed_xy <- unname(sf::st_coordinates(sf::st_transform(P$seeds, rec$crs))[, 1:2])
    expect_identical(.values_at_centres(attr(P$out, "count_rule"), seed_xy, rec$crs)$seed,
                     seq_len(nrow(cen)), info = kind)
    tol <- .edge_tol(c(range(cen[, 1L]), range(cen[, 2L])))
    off <- cbind(rep(c(1, -1), length.out = nrow(cen)), rep(c(-1, 1), length.out = nrow(cen)))
    moved <- cen + 5 * tol * off
    expect_gt(max(abs(moved - cen)), tol)
    expect_identical(.values_at_centres(attr(P$out, "count_rule"), moved, rec$crs)$seed,
                     seq_len(nrow(cen)), info = kind)
    .gt_expect_scored(P, pts, kind)
    .gt_expect_midpoints(P, kind)
  }
  # create_voronoi_polygons() draws them there as well.
  P <- suppressWarnings(.gt_path(pts, "density"))
  direct <- create_voronoi_polygons(P$seeds, quiet = TRUE)
  expect_true(sf::st_crs(direct$cells) == P$sel$count_rule$crs)
  expect_identical(sf::st_geometry(direct$cells), sf::st_geometry(P$tess$cells))
})

test_that("seeds placed in the CRS of a boundary are drawn in the profile's CRS; a `crs` the caller names is kept", {
  pts  <- .gt_lattice(40L, 30L, 250)
  bbox <- sf::st_sf(geometry = sf::st_as_sfc(sf::st_bbox(sf::st_buffer(sf::st_as_sfc(sf::st_bbox(pts)),
                                                                       250))))
  for (kind in c("density", "kmeans")) for (bc in c(5070, 4326)) {
    info <- paste(kind, bc)
    P <- suppressWarnings(.gt_path(pts, kind, boundary = sf::st_transform(bbox, bc)))
    expect_true(sf::st_crs(P$seeds) == sf::st_crs(bc), info = info)
    expect_true(sf::st_crs(P$tess$cells) == sf::st_crs(32615), info = info)
    expect_true(sf::st_crs(P$tess$boundary) == sf::st_crs(32615), info = info)
    .gt_expect_scored(P, pts, info)
    .gt_expect_midpoints(P, info)
  }
  # Cells returned in degrees are drawn in the profile's CRS too, and rows in
  # that CRS are joined there, where the profile's rule applies.
  P <- suppressWarnings(.gt_path(pts, "density", crs = 4326))
  expect_true(sf::st_is_longlat(P$tess$cells))
  rule <- .gt_profile_rule(P, pts)
  expect_identical(.gt_cell_of(P$asg, pts), rule$cell)
  expect_identical(attr(P$out, "count_rule")$check, list(estimate = "match", cells = "match"))
  # A projected `crs` is the CRS the cells are drawn in, as asked.
  P <- suppressWarnings(.gt_path(pts, "density", crs = 26915))
  expect_true(sf::st_crs(P$tess$cells) == sf::st_crs(26915))
})

test_that("on a ladder of sizes the record holds the centres of the level chosen, and the rows are placed by them", {
  pts <- .gt_lattice(30L, 20L, 30)
  for (kind in c("density", "kmeans")) for (size in c(18, 9)) {
    info <- paste(kind, size)
    P <- .gt_path(pts, kind, obs_per_cell = c(18, 9, 4.5), size = size)
    levels <- P$prof$levels
    expect_identical(levels, c(33L, 67L, 133L), info = info)
    expect_identical(P$sel$best, as.integer(round(600 / size)), info = info)
    cen <- P$sel$count_rule$scored$centres
    expect_identical(cen, unname(attr(P$prof, "centres")[[as.character(P$sel$best)]]), info = info)
    expect_identical(nrow(cen), P$sel$best, info = info)
    expect_identical(.scored_centres(attr(P$tess$cells, "count_rule")), cen, info = info)
    .gt_expect_scored(P, pts, info)
    .gt_expect_midpoints(P, info)
  }
})


# ---- What the polygons place alone -----------------------------------------------

test_that("a layer on which no row lies on an edge gets the cells it gets without the record", {
  layer <- local({
    .gt_seed(8101)
    x <- stats::runif(600, 0, 1e5); y <- stats::runif(600, 0, 1e5)
    sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y,
                            z = sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(600, 0, 1.5),
                            row = seq_len(600)),
                 coords = c("x", "y"), crs = 32617)
  })
  for (pts in list(layer, sf::st_transform(layer, 4326))) for (kind in c("density", "kmeans")) {
    P <- .gt_path(pts, kind)
    bare <- assign_features_to_polygons(pts, .gt_bare(P$tess$cells), polygon_id_col = "cell_id",
                                        quiet = TRUE)
    expect_identical(.gt_bare(P$asg), bare)
    expect_identical(attr(P$out, "count_rule")$check$cells, "match")
    expect_identical(.gt_cell_of(P$asg, pts), .gt_profile_rule(P, pts)$cell)
  }
  # The lattices on which the polygons and the profile agree: the same cells
  # and the same layer, record apart.
  for (d in .gt_dims[1:2]) for (sp in .gt_spacing) {
    pts <- .gt_lattice(d[1L], d[2L], sp)
    P <- .gt_lattice_path(d, sp, "kmeans")
    if (identical(.gt_polygon_rule(P, pts), .gt_profile_rule(P, pts)$cell))
      expect_identical(.gt_bare(P$asg),
                       assign_features_to_polygons(pts, .gt_bare(P$tess$cells),
                                                   polygon_id_col = "cell_id", quiet = TRUE))
  }
})

test_that("polygons a caller brings keep their rule for points on shared edges", {
  pts <- .gt_lattice(30L, 20L, 30)
  # Squares of 120 m from the lattice's corner, 8 across and 5 up, numbered
  # a row at a time from the lower left: every fourth row and column of the
  # lattice lies on their edges, and their corners on lattice points.
  sq <- list()
  for (j in 0:4) for (i in 0:7)
    sq[[length(sq) + 1L]] <- sf::st_polygon(list(
      cbind(4e5 + 120 * c(i, i + 1, i + 1, i, i), 4.4e6 + 120 * c(j, j, j + 1, j + 1, j))))
  units <- sf::st_sf(GEOID = sprintf("u%03d", seq_along(sq)),
                     geometry = sf::st_sfc(sq, crs = 32615))
  asg <- assign_features_to_polygons(pts, units, polygon_id_col = "GEOID", quiet = TRUE)
  expect_null(attr(asg, "count_rule"))
  expect_identical(nrow(asg), nrow(pts))
  # The rule of `tie_break = "smallest_area"`: of the equal squares a point
  # is on, the lowest, then the leftmost.
  dx <- sf::st_coordinates(pts)[, 1L] - 4e5
  dy <- sf::st_coordinates(pts)[, 2L] - 4.4e6
  col <- floor(dx / 120) - (dx %% 120 == 0 & dx > 0)
  row <- floor(dy / 120) - (dy %% 120 == 0 & dy > 0)
  want <- units$GEOID[row * 8L + col + 1L]
  expect_identical(.gt_cell_of(asg, pts, "GEOID"), want)
  expect_gt(attr(asg, "ties")$n, 100L)
  # The same whatever the order of the squares.
  rev_units <- units[rev(seq_len(nrow(units))), ]
  expect_identical(.gt_cell_of(assign_features_to_polygons(pts, rev_units, polygon_id_col = "GEOID",
                                                           quiet = TRUE), pts, "GEOID"), want)
  # And the summary checks nothing.
  r <- .gt_catch(summarize_by_cell(asg, "z", id_col = "GEOID", cells_sf = units,
                                   estimate = "pooled"))
  expect_null(attr(r$value, "count_rule"))
  expect_length(.gt_of(r, "spatialkit_cells_differ"), 0L)
})

test_that("cells whose record says they are not the scored cells keep the polygons' rule", {
  pts <- .gt_lattice(30L, 20L, 30)
  P <- .gt_lattice_path(c(30L, 20L), 30, "density")
  old <- .gt_polygon_rule(P, pts)
  for (mark in list(FALSE, NA)) {
    cells <- .count_rule_put(P$tess$cells, attr(P$tess$cells, "count_rule"), scored_cells = mark)
    asg <- assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE)
    expect_identical(.gt_cell_of(asg, pts), old)
    expect_identical(attr(asg, "count_rule")$scored_cells, mark)
    expect_identical(.gt_bare(asg), assign_features_to_polygons(pts, .gt_bare(P$tess$cells),
                                                                polygon_id_col = "cell_id",
                                                                quiet = TRUE))
  }
  # Seeds a function placed itself at the count, and a lattice sized by the
  # selection: their records say FALSE, and the polygons place every point.
  own <- suppressWarnings(get_voronoi_seeds(method = "kmeans",
                                            n = local({ s <- P$sel; s$seeds <- NULL; s }),
                                            sample_points = pts, set_seed = 1,
                                            partition = "kmeans", kmeans_nstart = 10))
  expect_identical(attr(own, "count_rule")$scored_cells, FALSE)
  hexes <- suppressWarnings(build_tessellation(pts, boundary = sf::st_sf(geometry = sf::st_as_sfc(
                              sf::st_bbox(pts))), method = "hex", approx_n_cells = P$sel,
                              quiet = TRUE))
  for (cells in list(build_tessellation(own, method = "voronoi", quiet = TRUE)$cells,
                     hexes$cells)) {
    expect_identical(attr(cells, "count_rule")$scored_cells, FALSE)
    expect_identical(.gt_bare(assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id",
                                                          quiet = TRUE)),
                     assign_features_to_polygons(pts, .gt_bare(cells), polygon_id_col = "cell_id",
                                                 quiet = TRUE))
  }
  # The rule changed something on this lattice, so the record is what decided.
  if (.gt_same_platform()) expect_identical(sum(old != .gt_cell_of(P$asg, pts)), 3L)
})


# ---- The rule on cells made by hand ---------------------------------------------------

# Four squares of side 10 from (5e5, 5e6) in EPSG:32617, numbered from the
# lower left a row at a time (1 LL, 2 LR, 3 UL, 4 UR), and a record whose
# centres are `centres` (square middles unless moved).  The record is a real
# one with its centres, count and CRS replaced.
.gt_squares <- function(crs = 32617, x0 = 5e5, y0 = 5e6) {
  sq <- function(i, j) sf::st_polygon(list(rbind(c(i, j), c(i + 10, j), c(i + 10, j + 10),
                                                 c(i, j + 10), c(i, j)) + rep(c(x0, y0), each = 5)))
  sf::st_sf(cell_id = 1:4, geometry = sf::st_sfc(sq(0, 0), sq(10, 0), sq(0, 10), sq(10, 10),
                                                 crs = crs))
}
.gt_middles <- function(x0 = 5e5, y0 = 5e6) cbind(x0 + c(5, 15, 5, 15), y0 + c(5, 5, 15, 15))
.gt_record <- function(centres, crs = sf::st_crs(32617), scored_cells = TRUE) {
  rec <- .gt_lattice_path(c(20L, 15L), 30, "density")$sel$count_rule
  rec$cells <- nrow(centres)
  rec$crs <- crs
  rec$scored$centres <- centres
  rec$scored_cells <- scored_cells
  rec
}
.gt_points <- function(x, y, crs = 32617, x0 = 5e5, y0 = 5e6)
  sf::st_as_sf(data.frame(x = x0 + x, y = y0 + y, row = seq_along(x)), coords = c("x", "y"),
               crs = crs)
# On the edge LL|LR, on LL|UL, on the corner of all four, on LR|UR; inside LL;
# inside UR.
.gt_six <- function(...) .gt_points(c(10, 5, 10, 15, 5, 12), c(5, 10, 10, 10, 5, 17), ...)
.gt_assign <- function(pts, cells, ...)
  assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE, ...)

test_that("on the cells of scored centres a point on an edge goes to the lowest-numbered of the centres it is nearest", {
  sq <- .gt_squares(); pts <- .gt_six()
  # The centres numbered from the upper right: 1 UR, 2 UL, 3 LR, 4 LL.
  cen <- .gt_middles()[c(4, 3, 2, 1), ]
  cells <- .count_rule_put(sq, .gt_record(cen))
  got <- .gt_assign(pts, cells)
  expect_identical(got$cell_id, c(2L, 3L, 4L, 4L, 1L, 4L))
  # The polygons alone: the lowest, then leftmost, of equal squares.
  expect_identical(.gt_assign(pts, sq)$cell_id, c(1L, 1L, 1L, 2L, 1L, 4L))
  # Whatever `tie_break` says, and whatever the order of the squares.
  expect_identical(.gt_assign(pts, cells, tie_break = "first")$cell_id, c(2L, 3L, 4L, 4L, 1L, 4L))
  rev_cells <- .count_rule_put(sq[4:1, ], .gt_record(cen))
  expect_identical(.gt_assign(pts, rev_cells)$cell_id, c(2L, 3L, 4L, 4L, 1L, 4L))
  # The record of the points that matched several polygons is the polygons' own.
  expect_identical(attr(got, "ties")[c("n", "which", "rule")],
                   attr(.gt_assign(pts, sq), "ties")[c("n", "which", "rule")])
  expect_identical(attr(got, "ties")$which, 1:4)
  # Numbered from the lower left the centres give what the polygons give.
  expect_identical(.gt_assign(pts, .count_rule_put(sq, .gt_record(.gt_middles())))$cell_id,
                   c(1L, 1L, 1L, 2L, 1L, 4L))
  # The rows come back as they came, with the cells' record.
  expect_identical(sf::st_geometry(got), sf::st_geometry(pts))
  expect_identical(attr(got, "count_rule"), .gt_record(cen))
})

test_that("a point rounding put a hair inside a cell is on its edge; one a little farther is not", {
  sq <- .gt_squares()
  tol <- .edge_tol(sf::st_bbox(sq))
  expect_equal(tol, 1e-10 * (5e6 + 20))
  # Centres moved right by 2e-5: the edge between their cells lies 2e-5 east
  # of the squares' edge, so a point 1e-5 east of it is nearer the left centre
  # and inside the right square, within the margin of its edge.
  near <- .gt_points(10 + 1e-5, 5)
  moved <- .gt_middles(); moved[, 1L] <- moved[, 1L] + 2e-5
  expect_identical(.gt_assign(near, sq)$cell_id, 2L)
  expect_identical(.gt_assign(near, .count_rule_put(sq, .gt_record(moved)))$cell_id, 1L)
  # Moved 2e-2, with a point 1e-2 east of the edge: nearer the left centre
  # too, and more than ten times the margin inside the right square, so it is
  # not on the edge, and the polygon that holds it decides.
  far <- .gt_points(10 + 1e-2, 5)
  moved[, 1L] <- .gt_middles()[, 1L] + 2e-2
  expect_gt(1e-2, 10 * tol)
  expect_identical(.gt_assign(far, .count_rule_put(sq, .gt_record(moved)))$cell_id, 2L)
})

test_that("a centre that is not alone in one cell gives no cell, and its points keep the polygons' rule", {
  sq <- .gt_squares(); pts <- .gt_six()
  # 1 UR, 2 UL, 3 LR, and 4 moved from LL into UR beside 1: LL holds no centre
  # and UR holds two.
  cen <- .gt_middles()[c(4, 3, 2, 1), ]
  cen[4L, ] <- c(5e5 + 17, 5e6 + 17)
  got <- .gt_assign(pts, .count_rule_put(sq, .gt_record(cen)))
  # LL|LR goes to 3 (LR) and LL|UL to 2 (UL), each alone in its square; the
  # corner and LR|UR are nearest 1, which shares UR, so the polygons decide.
  expect_identical(got$cell_id, c(2L, 3L, 1L, 2L, 1L, 4L))
})

test_that("each margin is the margin: a point a few margins from the other cell, or inside its own, keeps its cell", {
  tol <- .edge_tol(sf::st_bbox(.gt_squares()))
  box <- function(x0, x1) sf::st_polygon(list(cbind(5e5 + c(x0, x1, x1, x0, x0),
                                                    5e6 + c(0, 0, 10, 10, 0))))
  cells <- function(right) sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(box(0, 10), right,
                                                                         crs = 32617))
  # The left centre 5 from the point at (10, 5), the right one 4.9: the
  # profile's rule gives the point the right cell.
  cen <- cbind(5e5 + c(5, 14.9), 5e6 + c(5, 5))
  rec <- .gt_record(cen)
  # On the edge of its own square, five margins short of the right square:
  # not on the edge between the two.
  gap <- cells(box(10 + 5 * tol, 20))
  pt  <- .gt_points(10, 5)
  expect_identical(.gt_assign(pt, gap)$cell_id, 1L)
  expect_identical(.gt_assign(pt, .count_rule_put(gap, rec))$cell_id, 1L)
  # Two margins inside its own square and inside the right one, which
  # overlaps it by three: not on the edge of its own cell.
  overlap <- cells(box(10 - 3 * tol, 20))
  inside  <- .gt_points(10 - 2 * tol, 5)
  expect_identical(.gt_assign(inside, overlap)$cell_id, 1L)
  expect_identical(.gt_assign(inside, .count_rule_put(overlap, rec))$cell_id, 1L)
  # Half a margin inside its own square and a quarter short of the right
  # one: on the edge between them, and the profile's rule decides.
  overlap <- cells(box(10 - tol / 4, 20))
  edge    <- .gt_points(10 - tol / 2, 5)
  expect_identical(.gt_assign(edge, overlap)$cell_id, 1L)
  expect_identical(.gt_assign(edge, .count_rule_put(overlap, rec))$cell_id, 2L)
})

test_that("a centre on the edge of two cells has no cell, and its points keep the polygons' rule", {
  # The squares with LR first: the first square a centre on the edge LL|LR
  # touches is LR.
  sq <- .gt_squares()[c(2L, 1L, 3L, 4L), ]
  # 1 on the edge LL|LR, 2 inside LL, 3 UL, 4 UR.
  cen <- cbind(5e5 + c(10, 2, 5, 15), 5e6 + c(5, 2, 15, 15))
  pt  <- .gt_points(10, 2)
  xy  <- unname(sf::st_coordinates(pt))
  expect_identical(.nearest_exact(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L]), 1L)
  expect_identical(.gt_assign(pt, sq)$cell_id, 1L)
  expect_identical(.gt_assign(pt, .count_rule_put(sq, .gt_record(cen)))$cell_id, 1L)
})

test_that("a centre outside the boundary that clipped its cell gives that cell to its points on an edge", {
  x0 <- 5e5; y0 <- 5e6
  # B (15, 5) numbered 1, A (5, 5) 2, C (5, 15) 3, D (15, 15) 4.  The study
  # area is the left column [0, 10] x [0, 20] and a strip [10, 12] x [-200, 10]:
  # B and D lie outside it, B's cell is the strip, larger than A's square,
  # and D has none.
  cen <- cbind(x0 + c(15, 5, 5, 15), y0 + c(5, 5, 15, 15))
  bnd <- sf::st_sf(geometry = sf::st_union(sf::st_sfc(
    sf::st_polygon(list(cbind(x0 + c(0, 10, 10, 0, 0), y0 + c(0, 0, 20, 20, 0)))),
    sf::st_polygon(list(cbind(x0 + c(10, 12, 12, 10, 10), y0 + c(-200, -200, 10, 10, -200)))),
    crs = 32617)))
  seeds <- .count_rule_put(sf::st_as_sf(data.frame(seed_id = 1:4, x = cen[, 1L], y = cen[, 2L]),
                                        coords = c("x", "y"), crs = 32617), .gt_record(cen))
  tess <- suppressWarnings(build_tessellation(seeds, boundary = bnd, method = "voronoi",
                                              quiet = TRUE))
  expect_identical(tess$index, c(NA, 1L, 2L, NA))
  expect_equal(as.numeric(sf::st_area(tess$cells)), c(100, 100, 420))
  # (10, 5) and (10, 2) are as near B as A: B's cell, the strip.  (11, 5) is
  # in the strip.  The polygons alone give the first two the smaller square.
  pts <- .gt_points(c(10, 10, 11), c(5, 2, 5))
  expect_identical(.gt_assign(pts, tess$cells)$cell_id, c(3L, 3L, 3L))
  expect_identical(.gt_assign(pts, .gt_bare(tess$cells))$cell_id, c(1L, 1L, 3L))
  # Without the strip, B has no cell: a cell that holds a centre is that
  # centre's, so the point keeps the polygons' rule.
  cells_d <- .count_rule_put(.gt_bare(tess$cells)[c(1L, 2L), ], .gt_record(cen))
  expect_identical(.gt_assign(.gt_points(10, 5), cells_d)$cell_id, 1L)
})

test_that("a point on an edge whose nearest centre stands in a third cell keeps the polygons' rule", {
  sq <- .gt_squares()
  # One centre in each square, the one in UR near its lower-left corner, so
  # that the point on the edge LL|LR is nearest it: its cell is 5 away.
  cen <- cbind(5e5 + c(10.5, 2, 18, 2), 5e6 + c(10.5, 18, 2, 2))
  pt <- .gt_points(10, 5)
  xy <- unname(sf::st_coordinates(pt))
  expect_identical(.nearest_exact(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L]), 1L)
  got <- .gt_assign(pt, .count_rule_put(sq, .gt_record(cen)))
  expect_identical(got$cell_id, .gt_assign(pt, sq)$cell_id)
  expect_identical(got$cell_id, 1L)
})

test_that("only points on cells with the scored centres in a CRS the join can measure are placed by the profile's rule", {
  sq <- .gt_squares(); pts <- .gt_six()
  cen <- .gt_middles()[c(4, 3, 2, 1), ]
  polygons <- c(1L, 1L, 1L, 2L, 1L, 4L)
  profile  <- c(2L, 3L, 4L, 4L, 1L, 4L)
  with <- function(rec) .gt_assign(pts, .count_rule_put(sq, rec))$cell_id
  expect_identical(with(.gt_record(cen)), profile)
  # Not the scored cells, or not said to be.
  expect_identical(with(.gt_record(cen, scored_cells = FALSE)), polygons)
  expect_identical(with(.gt_record(cen, scored_cells = NA)), polygons)
  # No centres, too few, or not finite; no scored list.
  rec <- .gt_record(cen); rec$scored$centres <- NULL
  expect_identical(with(rec), polygons)
  rec <- .gt_record(cen); rec$scored$centres <- cen[1:3, ]
  expect_identical(with(rec), polygons)
  rec <- .gt_record(cen); rec$scored$centres[2L, 1L] <- NA
  expect_identical(with(rec), polygons)
  rec <- .gt_record(cen); rec["scored"] <- list(NULL)
  expect_identical(with(rec), polygons)
  # A record in longitude and latitude, or without a CRS beside cells with one.
  expect_identical(with(.gt_record(cen, crs = sf::st_crs(4326))), polygons)
  expect_identical(with(.gt_record(cen, crs = sf::NA_crs_)), polygons)
  # Without a CRS on all three, the coordinates are compared as they are.
  bare <- function(x) sf::st_set_crs(x, NA)
  expect_identical(.gt_assign(bare(pts), .count_rule_put(bare(sq),
                                                         .gt_record(cen, crs = sf::NA_crs_)))$cell_id,
                   profile)
  # Points without a CRS beside a record with one are left to the polygons.
  r <- .gt_catch(.gt_assign(bare(pts), .count_rule_put(sq, .gt_record(cen))))
  expect_identical(r$value$cell_id, polygons)
  # Features that are not points: squares split evenly by an edge, and
  # multipoints, are placed as the polygons place them.
  boxes <- sf::st_sf(row = 1:2, geometry = sf::st_buffer(sf::st_geometry(pts)[c(1L, 4L)], 1,
                                                         endCapStyle = "SQUARE"))
  expect_identical(.gt_assign(boxes, .count_rule_put(sq, .gt_record(cen)))$cell_id,
                   .gt_assign(boxes, sq)$cell_id)
  multi <- sf::st_sf(row = 1:2, geometry = sf::st_sfc(
    sf::st_multipoint(rbind(c(5e5 + 10, 5e6 + 5))), sf::st_multipoint(rbind(c(5e5 + 15, 5e6 + 10))),
    crs = 32617))
  expect_identical(.gt_assign(multi, .count_rule_put(sq, .gt_record(cen)))$cell_id, c(1L, 2L))
  # Cells and points in longitude and latitude: the join is on the sphere, and
  # the polygons decide.
  ll <- function(x) sf::st_transform(x, 4326)
  expect_identical(.gt_assign(ll(pts), .count_rule_put(ll(sq), .gt_record(cen)))$cell_id,
                   .gt_assign(ll(pts), ll(sq))$cell_id)
  # Two squares of 0.01 degree drawn in degrees on either side of the meridian
  # 80 W, five points on it, and a record in EPSG:32617 whose centres stand
  # near the meridian on one side or the other: the sphere decides each
  # point, whichever centre is nearer.
  deg <- function(x0) sf::st_polygon(list(cbind(x0 + c(0, 0.01, 0.01, 0, 0),
                                               45 + c(0, 0, 0.01, 0.01, 0))))
  cells_ll <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(deg(-80.01), deg(-80), crs = 4326))
  on_meridian <- sf::st_as_sf(data.frame(x = -80, y = 45 + 0.01 * (1:5) / 6, row = 1:5),
                              coords = c("x", "y"), crs = 4326)
  sphere <- .gt_assign(on_meridian, cells_ll)$cell_id
  for (at in list(c(-0.001, 0.009), c(-0.009, 0.001))) {
    mid <- sf::st_as_sf(data.frame(x = -80 + at, y = 45.005), coords = c("x", "y"), crs = 4326)
    centres <- unname(sf::st_coordinates(sf::st_transform(mid, 32617))[, 1:2])
    expect_identical(.gt_assign(on_meridian, .count_rule_put(cells_ll, .gt_record(centres)))$cell_id,
                     sphere)
  }
})

test_that("points that arrive in another CRS are measured in the record's, as the profile measured them", {
  sq <- .gt_squares()
  cen <- .gt_middles()[c(4, 3, 2, 1), ]
  cells <- .count_rule_put(sq, .gt_record(cen))
  # The six points and a lattice of 81 more inside the four squares, in degrees.
  g <- expand.grid(x = seq(2, 18, by = 2), y = seq(2, 18, by = 2))
  pts_ll <- sf::st_transform(rbind(.gt_six(), .gt_points(g$x, g$y)), 4326)
  got <- .gt_assign(pts_ll, cells)
  # The profile's rule on their coordinates in EPSG:32617, and the square of
  # the centre it gives.
  xy <- unname(sf::st_coordinates(sf::st_transform(pts_ll, 32617)))
  want <- c(4L, 3L, 2L, 1L)[.nearest_exact(xy[, 1L], xy[, 2L], cen[, 1L], cen[, 2L])]
  expect_identical(got$cell_id, want)
  expect_identical(sf::st_geometry(got), sf::st_geometry(pts_ll))
})


# ---- cell_values_at() ------------------------------------------------------------

test_that("cell_values_at() with the seeds in any order places the rows and the points between centres where the assignment does", {
  for (kind in c("density", "kmeans")) {
    pts <- .gt_lattice(30L, 20L, 30)
    P <- .gt_lattice_path(c(30L, 20L), 30, kind)
    placed <- .gt_cell_of(P$asg, pts)
    n <- nrow(P$seeds)
    for (o in list(rev(seq_len(n)), order(P$seeds$cell_id), local({ .gt_seed(4); sample.int(n) }))) {
      at <- cell_values_at(pts, P$out, P$asg, "z", id_col = "cell_id", seeds = P$seeds[o, ])
      expect_identical(at$cell_id, placed, info = kind)
    }
    .gt_expect_midpoints(P, kind)
  }
})

test_that("seeds are matched to the centres they stand at, whatever their order, and nothing else is", {
  .gt_seed(77)
  cen <- cbind(5e5 + stats::runif(40, 0, 1000), 5e6 + stats::runif(40, 0, 1000))
  o <- sample.int(40)
  for (fnn in unique(c(FALSE, .has_backend("FNN")))) {
    expect_identical(.seeds_as_centres(cen, cen, use_fnn = fnn), 1:40)
    expect_identical(.seeds_as_centres(cen[o, ], cen, use_fnn = fnn)[o], 1:40)
    # Seeds moved by rounding are the centres; seeds a metre away, two seeds
    # at one centre, a subset, or a seed that is not finite, are not.
    expect_identical(.seeds_as_centres(cen + 1e-6, cen, use_fnn = fnn), 1:40)
    moved <- cen; moved[7L, ] <- moved[7L, ] + 1
    expect_null(.seeds_as_centres(moved, cen, use_fnn = fnn))
    twice <- cen; twice[2L, ] <- cen[1L, ]
    expect_null(.seeds_as_centres(twice, cen, use_fnn = fnn))
    expect_null(.seeds_as_centres(cen[-1L, ], cen, use_fnn = fnn))
    nf <- cen; nf[3L, 1L] <- NA
    expect_null(.seeds_as_centres(nf, cen, use_fnn = fnn))
  }
  # With and without FNN, the nearest and the next nearest are the same.
  q <- cbind(5e5 + stats::runif(200, 0, 1000), 5e6 + stats::runif(200, 0, 1000))
  a <- .two_nearest(q[, 1L], q[, 2L], cen[, 1L], cen[, 2L], use_fnn = FALSE)
  expect_identical(a$first, .nearest_exact(q[, 1L], q[, 2L], cen[, 1L], cen[, 2L]))
  expect_true(all(a$d2 >= a$d1))
  if (.has_backend("FNN"))
    expect_identical(.two_nearest(q[, 1L], q[, 2L], cen[, 1L], cen[, 2L], use_fnn = TRUE), a)
})

test_that("cell_values_at() places a location on an edge where the assignment placed such a row, whatever rounding moved the seeds", {
  sq <- .gt_squares()
  cen <- .gt_middles()[c(4, 3, 2, 1), ]
  # Response "z" on the six points and on the middles; a record that speaks
  # for no estimate and no response, so the summary checks nothing.
  rec <- .gt_record(cen); rec$estimate <- NA_character_; rec$response <- NA_character_
  cells <- .count_rule_put(sq, rec)
  rows <- .gt_six()
  rows$z <- c(1, 2, 3, 4, 5, 6)
  asg <- .gt_assign(rows, cells)
  values <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = cells, estimate = "mean")
  expect_identical(attr(values, "count_rule")$check$cells,
                   "not checked: the response is another column")
  # The seeds as the centres, the first of them (UR) moved north by one step
  # of the arithmetic (1e-9 rounds to 9.3e-10 there), so that the corner and
  # the edge LR|UR are a little nearer other seeds than to it.
  at_seeds <- cen
  at_seeds[1L, 2L] <- at_seeds[1L, 2L] + 1e-9
  expect_gt(at_seeds[1L, 2L], cen[1L, 2L])
  seeds <- sf::st_as_sf(data.frame(x = at_seeds[, 1L], y = at_seeds[, 2L],
                                   cell_id = c(4L, 3L, 2L, 1L)),
                        coords = c("x", "y"), crs = 32617)
  newdata <- .gt_six()
  at <- cell_values_at(newdata, values, asg, "z", id_col = "cell_id", seeds = seeds)
  expect_identical(at$cell_id, asg$cell_id)
  expect_identical(at$cell_id, c(2L, 3L, 4L, 4L, 1L, 4L))
  # The seeds are the centres to rounding, in the profile's order.
  seed_xy <- unname(sf::st_coordinates(seeds))
  rule <- attr(values, "count_rule")
  expect_identical(.values_at_centres(rule, seed_xy, sf::st_crs(32617)),
                   list(centres = cen, seed = 1:4))
  # Seeds in any other order are the centres too: the search is over the
  # centres in the profile's order, and each gives the ID of its seed.
  for (o in list(c(2L, 1L, 3L, 4L), 4:1, order(seeds$cell_id), c(3L, 1L, 4L, 2L))) {
    got <- .values_at_centres(rule, seed_xy[o, ], sf::st_crs(32617))
    expect_identical(got$centres, cen)
    expect_identical(o[got$seed], 1:4)
    at_o <- cell_values_at(newdata, values, asg, "z", id_col = "cell_id", seeds = seeds[o, ])
    expect_identical(at_o$cell_id, asg$cell_id, info = paste(o, collapse = " "))
  }
  # A seed farther from its centre than the margin of an edge, by a
  # ten-thousandth of the distance between centres, is still that centre,
  # as seeds that went through another CRS and back are.
  tol <- .edge_tol(sf::st_bbox(sq))
  near <- seed_xy; near[1L, 1L] <- near[1L, 1L] + 1e-3
  expect_gt(1e-3, tol)
  expect_identical(.values_at_centres(rule, near, sf::st_crs(32617))$seed, 1:4)
  # A seed moved by a hundredth of that distance is not, nor are two seeds
  # at one centre, another number of seeds, or seeds in another CRS.
  far <- seed_xy; far[1L, 1L] <- far[1L, 1L] + 0.1
  expect_null(.values_at_centres(rule, far, sf::st_crs(32617)))
  twice <- seed_xy; twice[2L, ] <- seed_xy[1L, ]
  expect_null(.values_at_centres(rule, twice, sf::st_crs(32617)))
  expect_null(.values_at_centres(rule, seed_xy[1:3, ], sf::st_crs(32617)))
  expect_null(.values_at_centres(rule, seed_xy, sf::st_crs(32618)))
  expect_null(.values_at_centres(NULL, seed_xy, sf::st_crs(32617)))
  # Seeds that are not the centres are searched as given, and a tie goes to
  # the seed that comes first.
  far_seeds <- seeds; sf::st_geometry(far_seeds) <- sf::st_geometry(sf::st_as_sf(
    data.frame(x = far[, 1L], y = far[, 2L]), coords = c("x", "y"), crs = 32617))
  o <- c(2L, 1L, 3L, 4L)
  expect_identical(cell_values_at(newdata, values, asg, "z", id_col = "cell_id",
                                  seeds = far_seeds[o, ])$cell_id, c(2L, 3L, 3L, 2L, 1L, 4L))
  # Without seeds the polygons of `values` place the locations with the record.
  at_p <- cell_values_at(newdata, values, asg, "z", id_col = "cell_id")
  expect_identical(at_p$cell_id, asg$cell_id)
  # A summary without the record places them as the polygons alone do.
  bare <- values; attr(bare, "count_rule") <- NULL
  expect_identical(cell_values_at(newdata, bare, asg, "z", id_col = "cell_id")$cell_id,
                   c(1L, 1L, 1L, 2L, 1L, 4L))
  expect_identical(cell_values_at(newdata, bare, asg, "z", id_col = "cell_id",
                                  seeds = seeds)$cell_id, c(2L, 3L, 3L, 2L, 1L, 4L))
})


# ---- The record ----------------------------------------------------------------------

test_that("the record holds the centres the profile scored, as the profile holds them", {
  P <- .gt_lattice_path(c(30L, 20L), 30, "density")
  rec <- P$sel$count_rule
  cen <- rec$scored$centres
  expect_identical(utils::tail(names(rec$scored), 1L), "centres")
  expect_identical(cen, unname(attr(P$prof, "centres")[[as.character(P$sel$best)]]))
  expect_identical(cen, unname(sf::st_coordinates(P$sel$seeds)[, 1:2]))
  expect_identical(dim(cen), c(P$sel$best, 2L))
  # Every layer down the path carries them.
  for (x in list(P$sel$seeds, P$tess$cells, P$asg, P$out))
    expect_identical(attr(x, "count_rule")$scored$centres, cen)
  expect_identical(.scored_centres(attr(P$asg, "count_rule")), cen)
  expect_null(.scored_centres(rec))                       # beside no centres
  # A profile that no longer holds its centres writes none, and the record
  # keeps the field.
  bare <- P$prof; attr(bare, "centres") <- NULL
  r <- .count_rule_new(bare, "size", P$sel$best, scored_cells = TRUE)
  expect_identical(names(r$scored), names(rec$scored))
  expect_null(r$scored$centres)
  expect_null(.scored_centres(r))
  # The k-means kind's centres, as plain numbers.
  K <- .gt_lattice_path(c(30L, 20L), 30, "kmeans")$sel$count_rule$scored$centres
  expect_true(is.double(K) && is.null(dimnames(K)) && ncol(K) == 2L)
})


# ---- A real difference is still reported ------------------------------------------------

test_that("rows dropped, a seed moved, another estimate's cells and a changed value are still reported", {
  pts <- .gt_lattice(30L, 20L, 30)
  P <- .gt_lattice_path(c(30L, 20L), 30, "density")
  sum_of <- function(asg, cells = P$tess$cells, estimate = "pooled")
    .gt_catch(summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = cells,
                                estimate = estimate))
  differ <- function(r) .gt_of(r, "spatialkit_cells_differ")

  # Rows dropped.
  r <- sum_of(P$asg[-(1:30), ])
  expect_length(differ(r), 1L)
  expect_identical(differ(r)[[1L]]$what, "n")
  expect_identical(attr(r$value, "count_rule")$check$cells, "differ")

  # One seed moved one spacing east and one north: its cell is another, and
  # rows nowhere near an edge of the scored cells lie in it.  The record
  # still says these are the scored centres.
  xy <- unname(sf::st_coordinates(P$seeds)[, 1:2])
  k <- 60L
  xy[k, ] <- xy[k, ] + c(30, 30)
  seeds <- sf::st_as_sf(data.frame(seed_id = seq_len(nrow(xy)), x = xy[, 1L], y = xy[, 2L]),
                        coords = c("x", "y"), crs = sf::st_crs(P$seeds))
  attr(seeds, "cover") <- attr(P$seeds, "cover")
  seeds <- .count_rule_put(seeds, attr(P$seeds, "count_rule"))
  expect_true(isTRUE(attr(seeds, "count_rule")$scored_cells))
  tess <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
  asg  <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  rule <- .gt_profile_rule(list(sel = P$sel, tess = tess), pts)
  placed <- .gt_cell_of(asg, pts)
  elsewhere <- which(placed != rule$cell)
  expect_gt(length(elsewhere), 0L)
  d2 <- function(rows, centre) (rule$xy[rows, 1L] - rule$cen[centre, 1L])^2 +
                               (rule$xy[rows, 2L] - rule$cen[centre, 2L])^2
  gap <- abs(d2(elsewhere, match(placed[elsewhere], tess$index)) -
             d2(elsewhere, match(rule$cell[elsewhere], tess$index))) /
    d2(elsewhere, match(rule$cell[elsewhere], tess$index))
  expect_true(any(gap > 1e-3))
  r <- sum_of(asg, tess$cells)
  expect_length(differ(r), 1L)
  expect_true(differ(r)[[1L]]$what %in% c("cells", "loo"))
  expect_identical(attr(r$value, "count_rule")$check$cells, "differ")

  # Rows assigned to the cells of a count chosen for plain means, summarised
  # with the cells of the count chosen for pooled values: the two records
  # differ in the estimate.
  plain <- suppressWarnings(
    resolution_profile(pts, "z", partition = "density", nstart = 3, sac = "auto", seed = 123,
                       estimate = "mean", obs_per_cell = 4.5))
  sel_m <- select_resolution(plain, criterion = "size")
  seeds_m <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel_m, sample_points = pts))
  asg_m <- assign_features_to_polygons(pts, build_tessellation(seeds_m, method = "voronoi",
                                                               quiet = TRUE)$cells,
                                       polygon_id_col = "cell_id", quiet = TRUE)
  r <- sum_of(asg_m, estimate = "mean")
  expect_identical(differ(r)[[1L]]$what, "record")
  expect_identical(differ(r)[[1L]]$field, "estimate")
  # The pooled count's rows summarised as plain means: the estimate's warning.
  r <- sum_of(P$asg, estimate = "mean")
  expect_length(.gt_of(r, "spatialkit_estimate_mismatch"), 1L)
  expect_length(differ(r), 0L)

  # One value changed.
  changed <- P$asg
  changed$z[100L] <- changed$z[100L] + 5
  r <- sum_of(changed)
  expect_length(differ(r), 1L)
  expect_identical(differ(r)[[1L]]$what, "loo")
})
