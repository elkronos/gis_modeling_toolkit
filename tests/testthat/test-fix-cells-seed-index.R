# tests/testthat/test-fix-cells-seed-index.R
# ---------------------------------------------------------------------------
# Seeds outside a boundary that is not convex.
#
# A k-means seed is a mean of observations, so it can lie over a bay or a
# bend of a coast, outside the boundary its cell is clipped to.  `index` gave
# each point the cell it lies in, and such a seed lies in none, so it was NA
# although its cell was built and held observations; the warning told the
# user to pass a boundary that holds it.  With `seeds$cell_id <- tess$index`,
# cell_values_at(seeds = ) put every location nearest that seed in no cell
# and filled it from the cell of its nearest scored row, which is not always
# the cell that holds it.  The Voronoi builders now carry each cell's
# generator through the clip, and `index` gives each seed the cell it
# generates.  Points that are not seeds are indexed by the cell that holds
# them, as before; given to cell_values_at() as seeds with the cell polygons
# in the summary, one with an NA ID outside every cell takes the cell it
# generates there.
# ---------------------------------------------------------------------------

# Every warning an expression raises, muffled, and its value.
.si_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
# The columns of a frame and nothing else: no row names, no record.
.si_cols <- function(d) {
  d <- as.list(d)
  attributes(d) <- list(names = names(d))
  d
}
.si_unassigned <- function(r)
  Filter(function(w) inherits(w, "spatialkit_unassigned"), r$warnings)

# The cell of each location, by the polygons.
.si_holding <- function(newp, cells)
  cells$cell_id[vapply(sf::st_intersects(newp, cells),
                       function(i) if (length(i)) i[1L] else NA_integer_, integer(1))]

# Each cell's own seed: the seed nearest a point inside the cell, which no
# other seed is nearer (a Voronoi cell is the ground nearest its seed).
.si_own_seed <- function(cells, seeds)
  sf::st_nearest_feature(suppressWarnings(sf::st_point_on_surface(sf::st_geometry(cells))),
                         seeds)

.si_notch <- function(crs = 32617)
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(550, 1000), c(550, 300),
    c(450, 300), c(450, 1000), c(0, 1000), c(0, 0)))), crs = crs))
.si_six <- function(crs = 32617)
  sf::st_as_sf(data.frame(x = c(150, 850, 150, 850, 500, 500),
                          y = c(150, 150, 650, 650, 150, 800)),
               coords = c("x", "y"), crs = crs)


test_that("on the North Carolina outline every default seed has its own cell, those over the sounds included", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  nc <- sf::st_transform(sf::st_read(system.file("shape/nc.shp", package = "sf"),
                                     quiet = TRUE), 32617)
  st <- sf::st_sf(geometry = sf::st_union(nc))
  withr::local_seed(21)
  pts <- sf::st_sf(geometry = sf::st_sample(st, 2000))
  xy <- sf::st_coordinates(pts)
  pts$z <- sin(xy[, 1] / 1e5) + cos(xy[, 2] / 8e4) + stats::rnorm(2000, 0, 0.5)
  sel   <- select_resolution(resolution_profile(pts, "z"))
  seeds <- get_voronoi_seeds(boundary = st, method = "kmeans", n = sel, sample_points = pts)
  outside <- lengths(sf::st_intersects(seeds, st)) == 0L
  expect_gt(sum(outside), 0L)                  # the case is there
  b <- .si_catch(build_tessellation(seeds, boundary = st, method = "voronoi", quiet = TRUE))
  expect_length(.si_unassigned(b), 0L)
  tess <- b$value
  expect_identical(nrow(tess$cells), nrow(seeds))
  expect_false(anyNA(tess$index))
  expect_identical(sort(tess$index), seq_len(nrow(tess$cells)))
  expect_identical(tess$params$n_outside, 0L)
  # Each seed's entry is the cell it generates, for the seeds outside as for the rest.
  own <- .si_own_seed(tess$cells, seeds)
  expect_identical(tess$index[own], tess$cells$cell_id)
  expect_true(all(tess$index[outside] %in% tess$cells$cell_id))
  # With `index` alone: every location in the cell that holds it, none filled, and
  # each with its own cell's value.
  asg <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id")
  out <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells, obs_level = 0.9)
  withr::local_seed(5)
  newp <- sf::st_sf(geometry = sf::st_sample(st, 300))
  seeds$cell_id <- tess$index
  at <- cell_values_at(newp, out, asg, "z", id_col = "cell_id", seeds = seeds)
  holds <- .si_holding(newp, tess$cells)
  expect_false(anyNA(holds))
  expect_identical(at$cell_id, holds)
  expect_false(any(at$filled))
  expect_identical(at$value, out$..pooled_resp_z[match(holds, out$cell_id)])
})

test_that("on six notched layers every seed has its own cell and no location is filled, also for the selection's seeds bound with rbind()", {
  U <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(520, 1000), c(520, 200),
    c(480, 200), c(480, 1000), c(0, 1000), c(0, 0)))), crs = 32615))
  n_outside <- 0L
  n_bound_na <- 0L
  for (sd in 1:6) {
    withr::local_seed(sd)
    p  <- sf::st_sf(geometry = sf::st_sample(U, 400))
    p$z <- sf::st_coordinates(p)[, 1] / 500 + stats::rnorm(400, 0, 0.3)
    se <- select_resolution(resolution_profile(p, "z"))
    s  <- get_voronoi_seeds(boundary = U, method = "kmeans", n = se, sample_points = p)
    n_outside <- n_outside + sum(lengths(sf::st_intersects(s, U)) == 0L)
    b  <- .si_catch(build_tessellation(s, boundary = U, method = "voronoi", quiet = TRUE))
    expect_length(.si_unassigned(b), 0L)
    ts <- b$value
    expect_identical(sort(ts$index), seq_len(nrow(ts$cells)), info = paste("layer", sd))
    expect_identical(ts$index[.si_own_seed(ts$cells, s)], ts$cells$cell_id, info = paste("layer", sd))
    as <- assign_features_to_polygons(p, ts$cells, polygon_id_col = "cell_id")
    ou <- suppressWarnings(summarize_by_cell(as, "z", id_col = "cell_id", cells_sf = ts$cells))
    withr::local_seed(100 + sd)
    np <- sf::st_sf(geometry = sf::st_sample(U, 400))
    s$cell_id <- ts$index
    at <- cell_values_at(np, ou, as, "z", id_col = "cell_id", seeds = s)
    expect_identical(at$cell_id, .si_holding(np, ts$cells), info = paste("layer", sd))
    expect_false(any(at$filled), info = paste("layer", sd))
    # The selection's seeds bound with rbind() lose the attributes that make them
    # seeds, so `index` is NA for those outside the square; cell_values_at() gives
    # each of them its cell from the polygons of the summary.
    h  <- floor(nrow(se$seeds) / 2)
    rb <- rbind(se$seeds[seq_len(h), ], se$seeds[(h + 1L):nrow(se$seeds), ])
    tr <- .si_catch(build_tessellation(rb, boundary = U, method = "voronoi", quiet = TRUE))$value
    rb$cell_id <- tr$index
    n_bound_na <- n_bound_na + sum(is.na(rb$cell_id))
    ar <- assign_features_to_polygons(p, tr$cells, polygon_id_col = "cell_id")
    orr <- suppressWarnings(summarize_by_cell(ar, "z", id_col = "cell_id", cells_sf = tr$cells))
    ab <- cell_values_at(np, orr, ar, "z", id_col = "cell_id", seeds = rb)
    expect_identical(ab$cell_id, .si_holding(np, tr$cells), info = paste("layer", sd))
    expect_identical(ab$filled, !(ab$cell_id %in% ar$cell_id), info = paste("layer", sd))
  }
  expect_gt(n_outside, 0L)                     # the case is there
  expect_gt(n_bound_na, 0L)
})

test_that("every kind of seed layer keeps the cell it generates; other points keep the cell that holds them", {
  U <- .si_notch(); six <- .si_six()
  set.seed(4)
  obs <- sf::st_as_sf(data.frame(x = runif(200, 0, 1000), y = runif(200, 0, 300)),
                      coords = c("x", "y"), crs = 32617)
  plain <- suppressWarnings(build_tessellation(six, boundary = U, method = "voronoi", quiet = TRUE))
  expect_identical(which(is.na(plain$index)), 6L)
  notch_cell <- plain$cells$cell_id[lengths(sf::st_intersects(plain$cells, six)) == 0L]
  provided <- get_voronoi_seeds(method = "provided", seeds = six)
  as_kmeans <- provided; as_kmeans$method <- "kmeans"
  as_random <- provided; as_random$method <- "random"
  kinds <- list(
    provided  = list(points = provided),
    kmeans    = list(points = as_kmeans),
    random    = list(points = as_random),
    bound     = list(points = rbind(provided[1:3, ], provided[4:6, ])),
    attribute = list(points = structure(six, cover = obs)),
    kmeans_attr = list(points = structure(six, kmeans = list(k = 6L))),
    n_from_attr = list(points = structure(six, n_from = "selection")),
    cover     = list(points = six, cover = obs))
  for (k in names(kinds)) {
    x <- kinds[[k]]
    r <- .si_catch(build_tessellation(x$points, boundary = U, method = "voronoi",
                                      quiet = TRUE, cover = x$cover))
    expect_length(.si_unassigned(r), 0L)
    expect_identical(r$value$index, replace(plain$index, 6L, notch_cell), info = k)
    expect_identical(r$value$params$n_outside, 0L, info = k)
    # The cells are the cells the same points get when they are not seeds.
    expect_identical(r$value$cells, plain$cells, info = k)
    # create_voronoi_polygons() gives the same index.
    v <- suppressWarnings(create_voronoi_polygons(x$points, boundary = U, quiet = TRUE,
                                                  cover = x$cover))
    expect_identical(v$index, r$value$index, info = k)
  }
  # A MULTIPOINT seed: the smallest cell_id among the cells its vertices generate.
  mp <- sf::st_sf(seed_id = 1:5, method = "provided", geometry = sf::st_sfc(c(
    lapply(1:4, function(i) sf::st_geometry(six)[[i]]),
    list(sf::st_multipoint(rbind(c(500, 150), c(500, 800))))), crs = 32617))
  m <- build_tessellation(mp, boundary = U, method = "voronoi", quiet = TRUE)
  expect_identical(m$index, c(plain$index[1:4], min(plain$index[5L], notch_cell)))
  # So for every pair of the six points joined into one seed, whatever order the
  # cells had before the clip.
  own <- replace(plain$index, 6L, notch_cell)
  g <- sf::st_geometry(six)
  for (i in 1:5) for (j in (i + 1L):6) {
    pair <- sf::st_multipoint(rbind(sf::st_coordinates(g[i]), sf::st_coordinates(g[j])))
    lay <- sf::st_sf(seed_id = 1:5, method = "provided",
                     geometry = sf::st_sfc(c(list(pair), lapply(setdiff(1:6, c(i, j)),
                                                                function(k) g[[k]])),
                                           crs = 32617))
    mi <- build_tessellation(lay, boundary = U, method = "voronoi", quiet = TRUE)
    expect_identical(mi$index[1L], min(own[i], own[j]), info = paste(i, j))
    expect_identical(mi$cells, plain$cells, info = paste(i, j))
  }
  # A vertex far above the square has a cell the clip removes: the seed has the
  # cell of its other vertex, the notch's, though neither lies in a cell.
  far <- sf::st_sf(seed_id = 1:6, method = "provided", geometry = sf::st_sfc(c(
    lapply(1:5, function(k) g[[k]]),
    list(sf::st_multipoint(rbind(c(500, 800), c(500, 2500))))), crs = 32617))
  fr <- .si_catch(build_tessellation(far, boundary = U, method = "voronoi", quiet = TRUE))
  expect_length(.si_unassigned(fr), 0L)
  expect_identical(fr$value$index, own)
  # The same cells but for rounding: the far vertex changes the diagram outside
  # the square only.
  expect_identical(fr$value$cells$cell_id, plain$cells$cell_id)
  apart <- sf::st_sym_difference(sf::st_geometry(fr$value$cells), sf::st_geometry(plain$cells),
                                 by_element = TRUE)
  expect_lt(max(as.numeric(sf::st_area(apart))), 1e-6)
})

test_that("the record of the nearest-cell repair keeps only the seeds whose entry is the cell they were snapped to", {
  U <- .si_notch()
  # The first point 1e-4 below the square, a hair from its own cell; the sixth 1e-4
  # above the notch's floor, a hair from the fifth point's cell, its own cell being
  # the notch's.  Both are within the snapping distance.
  p <- sf::st_as_sf(data.frame(x = c(150, 850, 150, 850, 500, 500),
                               y = c(-1e-4, 150, 650, 650, 299.99995, 300.0001)),
                    coords = c("x", "y"), crs = 32617)
  expect_identical(lengths(sf::st_intersects(p, U)), c(0L, 1L, 1L, 1L, 1L, 0L))
  tp <- .si_catch(build_tessellation(p, boundary = U, method = "voronoi", quiet = TRUE))$value
  expect_identical(tp$params$snapped$which, c(1L, 6L))
  expect_identical(tp$index[6L], tp$index[5L])
  s <- get_voronoi_seeds(method = "provided", seeds = p)
  r <- .si_catch(build_tessellation(s, boundary = U, method = "voronoi", quiet = TRUE))
  expect_length(.si_unassigned(r), 0L)
  ts <- r$value
  expect_identical(ts$cells, tp$cells)
  expect_identical(ts$index[.si_own_seed(ts$cells, s)], ts$cells$cell_id)
  expect_identical(ts$index[1:5], tp$index[1:5])
  # The first seed keeps its entry in the record; the sixth, given its own cell,
  # leaves it.
  expect_identical(ts$params$snapped$n, 1L)
  expect_identical(ts$params$snapped$which, 1L)
  expect_identical(ts$params$snapped$distance, tp$params$snapped$distance[1L])
  expect_identical(ts$params$n_outside, 0L)
})

test_that("a record of the count makes seeds of its own seeds only: observations that carry it keep the cell that holds them", {
  U <- .si_notch(); six <- .si_six()
  withr::local_seed(1)
  pts <- sf::st_sf(geometry = sf::st_sample(U, 600))
  pts$z <- sin(sf::st_coordinates(pts)[, 1L] / 150) + stats::rnorm(600, 0, 0.4)
  sel  <- select_resolution(resolution_profile(pts, "z"))
  tess <- .si_catch(build_tessellation(sel$seeds, boundary = U, method = "voronoi",
                                       quiet = TRUE))$value
  # Thirty observations, some outside the square, assigned to the cells with
  # keep_unassigned = TRUE: the rows carry the cells' record of the count.
  withr::local_seed(9)
  obs <- sf::st_as_sf(data.frame(x = stats::runif(30, 0, 1000), y = stats::runif(30, 0, 1000),
                                 z = 1), coords = c("x", "y"), crs = 32617)
  away <- lengths(sf::st_intersects(obs, U)) == 0L
  expect_gt(sum(away), 0L)
  asg <- .si_catch(assign_features_to_polygons(obs, tess$cells, polygon_id_col = "cell_id",
                                               keep_unassigned = TRUE))$value
  rec <- attr(asg, "count_rule")
  expect_s3_class(rec, "spatialkit_count_rule")
  plain <- .si_catch(build_tessellation(obs, boundary = U, method = "voronoi", quiet = TRUE))
  r <- .si_catch(build_tessellation(asg, boundary = U, method = "voronoi", quiet = TRUE))
  expect_identical(is.na(r$value$index), away)
  expect_identical(r$value$index, plain$value$index)
  expect_identical(r$value$params$n_outside, sum(away))
  w <- .si_unassigned(r)
  expect_length(w, 1L)
  expect_match(conditionMessage(w[[1L]]), sprintf("%d of 30 point(s)", sum(away)), fixed = TRUE)
  v <- .si_catch(create_voronoi_polygons(asg, boundary = U, quiet = TRUE))$value
  expect_identical(v$index, plain$value$index)
  # A layer that carries a record and stands at its centres is seeds without a
  # seed_id column: the six points as the centres of a record of six cells.
  six_rec <- rec
  six_rec$cells <- 6L
  six_rec$crs <- sf::st_crs(32617)
  six_rec$scored$centres <- unname(sf::st_coordinates(six))
  six_rec$scored_cells <- TRUE
  bare <- suppressWarnings(build_tessellation(six, boundary = U, method = "voronoi", quiet = TRUE))
  notch_cell <- bare$cells$cell_id[lengths(sf::st_intersects(bare$cells, six)) == 0L]
  at_c <- .si_catch(build_tessellation(.count_rule_put(six, six_rec), boundary = U,
                                       method = "voronoi", quiet = TRUE))
  expect_length(.si_unassigned(at_c), 0L)
  expect_identical(at_c$value$index, replace(bare$index, 6L, notch_cell))
  # The same record on the six points with one moved 50 m off its centre: they
  # are not its centres, so they are not seeds and the notch point is NA.
  moved <- sf::st_as_sf(data.frame(x = c(200, 850, 150, 850, 500, 500),
                                   y = c(150, 150, 650, 650, 150, 800)),
                        coords = c("x", "y"), crs = 32617)
  off <- .si_catch(build_tessellation(.count_rule_put(moved, six_rec), boundary = U,
                                      method = "voronoi", quiet = TRUE))
  expect_identical(which(is.na(off$value$index)), 6L)
  expect_match(conditionMessage(.si_unassigned(off)[[1L]]), "1 of 6 point(s)", fixed = TRUE)
})

test_that("a seed whose cell the clip removed is NA, and the warning names seeds", {
  set.seed(3)
  pts <- sf::st_as_sf(data.frame(x = runif(300, 0, 1000), y = runif(300, 0, 1000)),
                      coords = c("x", "y"), crs = 32617)
  seeds <- voronoi_seeds_kmeans(pts, k = 12)
  corner <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(400, 0), c(400, 400), c(0, 400), c(0, 0)))), crs = 32617))
  r <- .si_catch(build_tessellation(seeds, boundary = corner, method = "voronoi", quiet = TRUE))
  w <- .si_unassigned(r)
  expect_length(w, 1L)
  tess <- r$value
  # NA exactly for the seeds whose whole Voronoi cell lies outside the boundary.
  vor <- sf::st_collection_extract(sf::st_voronoi(sf::st_union(sf::st_geometry(seeds)),
                                                  envelope = sf::st_geometry(corner)),
                                   "POLYGON", warn = FALSE)
  vor <- sf::st_sfc(vor, crs = 32617)
  own <- vapply(sf::st_intersects(seeds, vor), function(i) i[1L], integer(1))
  reaches <- vapply(own, function(i)
    isTRUE(as.numeric(sf::st_area(sf::st_intersection(vor[i], sf::st_geometry(corner)))) > 0),
    logical(1))
  expect_identical(is.na(tess$index), !reaches)
  n_out <- sum(!reaches)
  expect_gt(n_out, 0L)
  expect_gt(sum(reaches & lengths(sf::st_intersects(seeds, corner)) == 0L), 0L)
  expect_identical(tess$params$n_outside, n_out)
  expect_identical(w[[1L]]$n, n_out)
  expect_identical(w[[1L]]$total, 12L)
  expect_identical(conditionMessage(w[[1L]]), sprintf(paste0(
    "build_tessellation(): %d of 12 seed(s) (%.1f%%) have no cell and are NA in `index`; ",
    "`params$n_outside` holds the count. They lie outside the boundary (grown by `expand`), ",
    "and every part of the study area is nearer another seed, so the clip left nothing of ",
    "their cells. A seed outside the boundary whose cell reaches into it keeps that cell."),
    n_out, 100 * n_out / 12))
  # Every seed with a cell has its own.
  expect_identical(sort(tess$index[!is.na(tess$index)]), seq_len(nrow(tess$cells)))
  # The same points that are not seeds: every one outside the boundary is NA, in
  # the words for points, and the cells are the same.
  bare <- sf::st_sf(geometry = sf::st_geometry(seeds))
  rb <- .si_catch(build_tessellation(bare, boundary = corner, method = "voronoi", quiet = TRUE))
  wb <- .si_unassigned(rb)
  expect_length(wb, 1L)
  expect_match(conditionMessage(wb[[1L]]), "point(s)", fixed = TRUE)
  expect_match(conditionMessage(wb[[1L]]), "Pass a boundary that holds them", fixed = TRUE)
  expect_identical(is.na(rb$value$index), lengths(sf::st_intersects(seeds, corner)) == 0L)
  expect_identical(rb$value$cells, tess$cells)
  # A seed with an empty geometry has no location to index, and is counted with
  # them.  The one emptied had no cell, so the cells in the boundary stay.
  expect_gt(n_out, 1L)
  g <- sf::st_geometry(seeds)
  g[[which(!reaches)[1L]]] <- sf::st_point()
  gone <- sf::st_set_geometry(seeds, g)
  re <- .si_catch(build_tessellation(gone, boundary = corner, method = "voronoi", quiet = TRUE))
  we <- .si_unassigned(re)
  expect_length(we, 1L)
  expect_identical(re$value$index, tess$index)
  expect_identical(re$value$params$n_outside, n_out)
  expect_match(conditionMessage(we[[1L]]), sprintf("%d of 12 seed(s) (", n_out), fixed = TRUE)
  expect_match(conditionMessage(we[[1L]]), paste0(
    "holds the count. 1 of them have empty geometry; the rest lie outside the ",
    "boundary (grown by `expand`), and every part of the study area is nearer ",
    "another seed"), fixed = TRUE)
  # The bare points as the seeds of cell_values_at(): each one outside the
  # boundary whose cell reaches into it takes that cell, since it lies in no cell
  # of the summary and its cell is the one no point names; one whose cell is gone
  # keeps NA.  Every location in the corner has the cell that holds it, and only
  # those nearest a seed with no cell are filled.
  pts$z <- sf::st_coordinates(pts)[, 1L] / 100
  asg <- .si_catch(assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id"))$value
  out <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                           estimate = "mean")
  h <- expand.grid(x = seq(5, 995, 30), y = seq(5, 995, 30))
  newp <- sf::st_as_sf(h, coords = c("x", "y"), crs = 32617)
  bare$cell_id <- rb$value$index
  expect_gt(sum(is.na(bare$cell_id) & reaches), 0L)
  at <- cell_values_at(newp, out, asg, "z", id_col = "cell_id", seeds = bare)
  seeds$cell_id <- tess$index
  expect_identical(.si_cols(at), .si_cols(cell_values_at(newp, out, asg, "z",
                                                         id_col = "cell_id", seeds = seeds)))
  inside <- lengths(sf::st_intersects(newp, corner)) > 0L
  expect_identical(at$cell_id[inside], .si_holding(newp[inside, ], tess$cells))
  gone_near <- is.na(tess$index)[sf::st_nearest_feature(newp, seeds)]
  expect_gt(sum(gone_near), 0L)
  expect_false(any(inside & gone_near))
  # Filled: nearest a seed with no cell, or in a cell that holds no scored row.
  expect_identical(at$filled, gone_near | !(at$cell_id %in% asg$cell_id))
})
