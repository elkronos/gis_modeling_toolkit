# tests/testthat/test-fix-crs-area-refused-probe.R
# ---------------------------------------------------------------------------
# The area distortion of a projected layer is measured over every probe,
# also one s2 refuses in degrees because two of its vertices fall on one
# point there.  .crs_area_error() took the geodesic areas of all its probes
# in one call on their copy in EPSG:4326.  K-means cells built from rows in
# degrees can have two corners a fraction of a nanometre apart; in degrees
# they fall on one point, s2 refuses the ring ("duplicate vertex"), and the
# one call failed, so the figure was NA for the whole layer.  An exact
# repeated vertex, which GEOS accepts, does the same.
# summarize_by_cell(area = TRUE) then refused cells that are valid in their
# own CRS, and ensure_projected(purpose = "area") skipped its check on a
# projected layer that distorts areas.  When the one call fails, each probe
# is now measured on its own, a probe s2 refuses is measured again without
# the vertices that repeat the one before them in degrees, and a probe s2
# refuses even then keeps the figure NA.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
square <- function(x0, y0, w, h = w)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
# The polygons of `cells` s2 refuses once in EPSG:4326.
s2_refused <- function(cells) {
  ll <- sf::st_transform(sf::st_geometry(cells), 4326)
  which(vapply(seq_along(ll), function(i) inherits(tryCatch(
    spatialkit:::.with_s2(sf::st_area(ll[i])), error = function(e) e), "error"), logical(1)))
}
# Polygon `p` in `crs` with a vertex added one to four ulps from its first
# corner: along x on the first edge, or along y on the closing edge.  Either
# way the ring is valid in `crs` (GEOS) and the new vertex can fall on the
# corner in degrees, where s2 refuses the ring.  Which offset does that
# depends on the rounding of the transform, so the first that does is taken;
# when none does, the second vertex is repeated, which GEOS accepts and s2
# refuses too.
with_near_corner <- function(p, crs) {
  m <- unclass(p)[[1]]
  nr <- nrow(m)
  for (k in 1:4) for (axis in 1:2) {
    v <- m[1, ]
    v[axis] <- v[axis] + k * 2^(floor(log2(abs(v[axis]))) - 52)
    q <- if (axis == 1L) sf::st_polygon(list(rbind(m[1, ], v, m[-1, ])))
         else sf::st_polygon(list(rbind(m[-nr, ], v, m[nr, ])))
    s <- sf::st_sfc(q, crs = crs)
    if (isTRUE(sf::st_is_valid(s)) && length(s2_refused(s))) return(q)
  }
  with_repeat(p)
}
# Polygon `p` with its second vertex listed twice.
with_repeat <- function(p) {
  m <- unclass(p)[[1]]
  sf::st_polygon(list(rbind(m[1:2, ], m[2, ], m[-(1:2), ])))
}
# A bowtie: a ring whose edges cross, with lobes of 1.125e7 and 3.125e7 m^2.
bowtie <- function(x0, y0)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1e4, y0 + 1e4), c(x0 + 1e4, y0),
                            c(x0, y0 + 6e3), c(x0, y0))))
# A square of 10 km with a spike of 10 km up from the middle of its top edge.
spike <- function(x0, y0)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 1e4, y0), c(x0 + 1e4, y0 + 1e4),
                            c(x0 + 5e3, y0 + 1e4), c(x0 + 5e3, y0 + 2e4),
                            c(x0 + 5e3, y0 + 1e4), c(x0, y0 + 1e4), c(x0, y0))))
# One polygon in degrees without the vertices whose x and y repeat those of
# the vertex before them, rebuilt from st_coordinates().
drop_repeats <- function(s) {
  xy <- sf::st_coordinates(s)
  multi <- inherits(s[[1]], "MULTIPOLYGON")
  # L1 is the ring; in a MULTIPOLYGON, L2 is the part.
  key <- if (multi) paste(xy[, "L2"], xy[, "L1"]) else as.character(xy[, "L1"])
  rep <- c(FALSE, xy[-1L, "X"] == xy[-nrow(xy), "X"] & xy[-1L, "Y"] == xy[-nrow(xy), "Y"] &
                    key[-1L] == key[-nrow(xy)])
  xy <- xy[!rep, , drop = FALSE]
  key <- key[!rep]
  rings <- lapply(split(seq_len(nrow(xy)), factor(key, unique(key))),
                  function(i) xy[i, c("X", "Y"), drop = FALSE])
  geom <- if (multi) {
    part <- sub(" .*", "", names(rings))
    sf::st_multipolygon(unname(lapply(split(rings, factor(part, unique(part))), unname)))
  } else sf::st_polygon(unname(rings))
  sf::st_sfc(geom, crs = sf::st_crs(s))
}
# The figure of a layer of at most 200 polygons, step by step: the polygons
# densified to a hundredth of the layer's extent, their planar areas, their
# geodesic areas taken one polygon at a time, a polygon s2 refuses measured
# again without its repeated vertices, NA when s2 refuses it even then, and
# the largest relative departure of the ratio from its median.
spread_each <- function(cells) {
  g  <- sf::st_geometry(cells)
  bb <- sf::st_bbox(g)
  p  <- sf::st_segmentize(g, dfMaxLength = max(bb[["xmax"]] - bb[["xmin"]],
                                               bb[["ymax"]] - bb[["ymin"]]) / 100)
  ll <- sf::st_transform(p, 4326)
  area_s2 <- function(s) as.numeric(spatialkit:::.with_s2(sf::st_area(s)))
  geod <- vapply(seq_along(ll), function(i) tryCatch(area_s2(ll[i]), error = function(e)
    tryCatch(area_s2(drop_repeats(ll[i])), error = function(e) NA_real_)), numeric(1))
  if (anyNA(geod)) return(NA_real_)
  r <- as.numeric(sf::st_area(p)) / geod
  max(abs(r / stats::median(r) - 1))
}
# The same with the geodesic areas taken in one call, as before the fix.
spread_one_call <- function(cells) {
  g  <- sf::st_geometry(cells)
  bb <- sf::st_bbox(g)
  p  <- sf::st_segmentize(g, dfMaxLength = max(bb[["xmax"]] - bb[["xmin"]],
                                               bb[["ymax"]] - bb[["ymin"]]) / 100)
  r <- as.numeric(sf::st_area(p)) /
    spatialkit:::.with_s2(as.numeric(sf::st_area(sf::st_transform(p, 4326))))
  max(abs(r / stats::median(r) - 1))
}
# `cells` with each vertex dropped, in the cells' own CRS, whose image in
# degrees is the image of the vertex before it: a layer s2 accepts whole.
drop_in_own_crs <- function(cells) {
  g  <- sf::st_geometry(cells)
  ll <- sf::st_transform(g, 4326)
  for (i in seq_along(g)) {
    m <- unclass(g[[i]])
    d <- unclass(ll[[i]])
    for (j in seq_along(m)) {
      dj <- d[[j]]
      same <- c(FALSE, dj[-1L, 1L] == dj[-nrow(dj), 1L] & dj[-1L, 2L] == dj[-nrow(dj), 2L])
      m[[j]] <- m[[j]][!same, , drop = FALSE]
    }
    g[[i]] <- sf::st_polygon(m)
  }
  sf::st_set_geometry(cells, g)
}
# Rows at the cells' centroids, one per cell, for summarize_by_cell().
rows_in <- function(cells, ids = cells$cell_id)
  sf::st_sf(cell_id = ids, z = seq_along(ids),
            geometry = sf::st_centroid(sf::st_geometry(cells)))

# The layer of BUG-09: 30 x 20 locations 30 m apart in EPSG:32615, a smooth
# field plus noise drawn under seed 9001; the rows in `crs`.  Its k-means
# cells by size, kept in the CRS the profile worked in, with the rows
# assigned to them.
.rp_store <- new.env(parent = emptyenv())
kmeans_cells <- function(crs) {
  key <- as.character(crs)
  if (!is.null(.rp_store[[key]])) return(.rp_store[[key]])
  g <- expand.grid(x = 4e5 + 30 * (0:29), y = 4.4e6 + 30 * (0:19))
  z <- withr::with_seed(9001, sin(g$x / 180) + cos(g$y / 150) + stats::rnorm(600, 0, 0.4),
                        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                        .rng_sample_kind = "Rejection")
  pts <- sf::st_as_sf(data.frame(g, z = z, row = 1:600), coords = c("x", "y"), crs = 32615)
  if (crs != 32615) pts <- sf::st_transform(pts, crs)
  prof  <- quiet(resolution_profile(pts, "z", partition = "kmeans", nstart = 3, sac = "auto",
                                    seed = 123, estimate = "pooled", obs_per_cell = 4.5))
  sel   <- select_resolution(prof, criterion = "size")
  seeds <- quiet(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  tess  <- quiet(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
  asg   <- assign_features_to_polygons(pts, tess$cells, polygon_id_col = "cell_id", quiet = TRUE)
  .rp_store[[key]] <- list(cells = tess$cells, asg = asg)
  .rp_store[[key]]
}

test_that("k-means cells kept in EPSG:32615 for rows in degrees get an area check and a density", {
  k <- kmeans_cells(4326)
  cells <- k$cells
  expect_true(sf::st_crs(cells) == sf::st_crs(32615))
  expect_true(all(sf::st_is_valid(cells)))
  e <- .crs_area_error(cells)
  expect_true(is.finite(e))
  expect_lt(e, .area_error_tol)
  # The spread over every cell, measured step by step.
  expect_identical(e, spread_each(cells))
  # Whatever sf_use_s2() says, and in any row order.
  was <- suppressMessages(sf::sf_use_s2(FALSE))
  withr::defer(suppressMessages(sf::sf_use_s2(was)))
  e_off <- .crs_area_error(cells)
  suppressMessages(sf::sf_use_s2(was))
  expect_identical(e_off, e)
  for (o in withr::with_seed(5, lapply(1:3, function(i) sample(nrow(cells)))))
    expect_identical(.crs_area_error(cells[o, ]), e)
  # summarize_by_cell(area = TRUE) gives every cell its planar area in
  # EPSG:32615 and a density.
  out <- quiet(summarize_by_cell(k$asg, "z", id_col = "cell_id", cells_sf = cells,
                                 estimate = "pooled", area = TRUE))
  expect_identical(attr(out, "area_error"), e)
  expect_identical(nrow(out), nrow(cells))
  expect_true(all(is.finite(out$cell_area) & out$cell_area > 0))
  expect_equal(out$cell_area, as.numeric(sf::st_area(out)))
  expect_equal(out$n_per_area, out$n / out$cell_area)
  expect_equal(sum(out$n), 600)
})

test_that("on those cells the cells s2 refuses in degrees are measured, not left out", {
  # Which cells s2 refuses depends on the rounding of the transforms (here 4
  # of the 133).  Each has two corners that fall on one point in degrees;
  # with one of the two dropped in EPSG:32615, s2 accepts every cell, and
  # the figure is that layer's to within the nanometre the corners are apart.
  cells <- kmeans_cells(4326)$cells
  bad <- s2_refused(cells)
  skip_if(!length(bad), "s2 refuses none of these cells in degrees here")
  fixed <- drop_in_own_crs(cells)
  expect_length(s2_refused(fixed), 0L)
  e <- .crs_area_error(cells)
  expect_true(is.finite(e))
  expect_equal(e, .crs_area_error(fixed), tolerance = 1e-6)
})

test_that("cells whose probes s2 all accepts keep their figure to the last bit", {
  # The same layer with the rows in EPSG:32615: no cell is refused, and the
  # figure is the one-call figure.
  cells <- kmeans_cells(32615)$cells
  expect_length(s2_refused(cells), 0L)
  e <- .crs_area_error(cells)
  expect_identical(e, spread_one_call(cells))
  expect_identical(e, spread_each(cells))
  expect_equal(e, 2.6934256306399007e-06, tolerance = 1e-4)
  # Twenty squares of 1 degree from 30N to 69N in Web Mercator, and 36 squares
  # of 10 km in UTM zone 15.
  merc <- sf::st_sf(cell_id = 1:20, geometry = sf::st_transform(sf::st_sfc(lapply(0:19, function(k)
    square(-100, 30 + 2 * k, 1)), crs = 4326), 3857))
  utm <- sf::st_sf(cell_id = 1:36, geometry = sf::st_sfc(lapply(0:35, function(k)
    square(4e5 + 1e4 * (k %% 6), 4.4e6 + 1e4 * (k %/% 6), 1e4)), crs = 32615))
  for (lay in list(merc, utm)) {
    expect_length(s2_refused(lay), 0L)
    expect_identical(.crs_area_error(lay), spread_one_call(lay))
  }
})

test_that("a valid cell s2 refuses in a CRS that distorts areas is measured, and the distortion reported", {
  # The Web Mercator squares, the eleventh (50N) with a vertex beside its
  # first corner that falls on it in degrees.  The figure is that of all
  # twenty, where the largest ratio is about three times the median; over
  # the other nineteen alone it is not.
  g <- sf::st_transform(sf::st_sfc(lapply(0:19, function(k) square(-100, 30 + 2 * k, 1)),
                                   crs = 4326), 3857)
  g_bad <- g; g_bad[[11]] <- with_near_corner(g[[11]], 3857)
  bad  <- sf::st_sf(cell_id = 1:20, geometry = g_bad)
  full <- sf::st_sf(cell_id = 1:20, geometry = g)
  rest <- sf::st_sf(cell_id = c(1:10, 12:20), geometry = g[-11])
  expect_true(all(sf::st_is_valid(bad)))
  expect_identical(s2_refused(bad), 11L)
  e <- .crs_area_error(bad)
  expect_identical(e, spread_each(bad))
  expect_equal(e, .crs_area_error(full), tolerance = 1e-9)
  expect_gt(abs(e - .crs_area_error(rest)), 0.1)
  expect_gt(e, 2)
  # ensure_projected(purpose = "area") logs the distortion; it logged nothing.
  lines <- capture_spatialkit_log(out <- ensure_projected(bad, purpose = "area"))
  expect_true(sf::st_crs(out) == sf::st_crs(3857))
  expect_true(log_has(lines, sprintf("distorts areas across this extent by up to %.1f%%", 100 * e)))
  # summarize_by_cell(area = TRUE) refuses the cells for their distortion,
  # not because it could not measure it.
  expect_error(quiet(summarize_by_cell(rows_in(full), "z", id_col = "cell_id", cells_sf = bad,
                                       estimate = "mean", area = TRUE)),
               sprintf("distorts areas across the cells by up to %.1f%%", 100 * e), fixed = TRUE)
})

test_that("valid cells s2 refuses are measured however many of a layer there are", {
  # Squares of 10 km in UTM zone 15 with a vertex beside the first corner.
  # Two of two, two of three and one of three refused: each layer gets the
  # figure of its squares without the added vertex.
  sq <- lapply(0:2, function(k) square(4e5 + 1e4 * k, 4.4e6, 1e4))
  cells <- function(...) sf::st_sf(cell_id = seq_along(list(...)),
                                   geometry = sf::st_sfc(list(...), crs = 32615))
  near <- lapply(sq, with_near_corner, crs = 32615)
  one_of_two   <- cells(sq[[1]], near[[2]])
  two_of_three <- cells(near[[1]], near[[2]], sq[[3]])
  one_of_three <- cells(sq[[1]], near[[2]], sq[[3]])
  expect_true(all(sf::st_is_valid(two_of_three)))
  expect_identical(s2_refused(one_of_two), 2L)
  expect_identical(s2_refused(two_of_three), 1:2)
  expect_identical(s2_refused(one_of_three), 2L)
  for (lay in list(one_of_two, two_of_three, one_of_three)) {
    e <- .crs_area_error(lay)
    expect_true(is.finite(e))
    expect_identical(e, spread_each(lay))
    plain <- cells(sq[[1]], sq[[2]], sq[[3]])[seq_len(nrow(lay)), ]
    expect_equal(e, .crs_area_error(plain), tolerance = 1e-6)
  }
  out <- quiet(summarize_by_cell(rows_in(one_of_two), "z", id_col = "cell_id",
                                 cells_sf = one_of_two, estimate = "mean", area = TRUE))
  expect_identical(attr(out, "area_error"), .crs_area_error(one_of_two))
  expect_equal(out$cell_area, rep(1e8, 2), tolerance = 1e-12)
})

test_that("the cells that carry the distortion are measured when they repeat a vertex", {
  # Twenty cells in Web Mercator, all valid: two of half a degree at 30N,
  # and eighteen of 1 degree from 32N to 66N, each with one vertex listed
  # twice.  Left out, the eighteen took the distortion with them (a figure
  # of 7e-15 over the two left, and a density in a CRS whose ratio is six
  # times larger at 66N than at 30N).
  ll <- c(list(square(-100, 30, 0.5), square(-99.5, 30, 0.5)),
          lapply(1:18, function(k) square(-100, 30 + 2 * k, 1)))
  g <- sf::st_transform(sf::st_sfc(ll, crs = 4326), 3857)
  g2 <- g
  for (i in 3:20) g2[[i]] <- with_repeat(g[[i]])
  lay  <- sf::st_sf(cell_id = 1:20, geometry = g2)
  full <- sf::st_sf(cell_id = 1:20, geometry = g)
  expect_true(all(sf::st_is_valid(lay)))
  expect_identical(s2_refused(lay), 3:20)
  e <- .crs_area_error(lay)
  expect_equal(e, .crs_area_error(full), tolerance = 1e-12)
  expect_gt(e, 1.5)
  lines <- capture_spatialkit_log(out <- ensure_projected(lay, purpose = "area"))
  expect_true(log_has(lines, sprintf("distorts areas across this extent by up to %.1f%%", 100 * e)))
  expect_error(quiet(summarize_by_cell(rows_in(full), "z", id_col = "cell_id", cells_sf = lay,
                                       estimate = "mean", area = TRUE)),
               sprintf("distorts areas across the cells by up to %.1f%%", 100 * e), fixed = TRUE)
  # A cell of two parts, one of them with a repeated vertex, and a cell
  # whose repeated vertex differs from the one before it in z alone.
  sq <- lapply(0:3, function(k) square(4e5 + 1e4 * k, 4.4e6, 1e4))
  mp <- sf::st_multipolygon(list(unclass(sq[[1]]), unclass(with_repeat(sq[[3]]))))
  mp0 <- sf::st_multipolygon(list(unclass(sq[[1]]), unclass(sq[[3]])))
  lay_mp  <- sf::st_sf(cell_id = 1:3, geometry = sf::st_sfc(list(mp, sq[[2]], sq[[4]]), crs = 32615))
  lay_mp0 <- sf::st_sf(cell_id = 1:3, geometry = sf::st_sfc(list(mp0, sq[[2]], sq[[4]]), crs = 32615))
  expect_identical(s2_refused(lay_mp), 1L)
  expect_identical(.crs_area_error(lay_mp), .crs_area_error(lay_mp0))
  expect_identical(.crs_area_error(lay_mp), spread_each(lay_mp))
  # A cell with a hole, the hole's ring with a repeated vertex.
  hole  <- unclass(square(4.04e5, 4.404e6, 2e3))[[1]][5:1, ]
  holed <- sf::st_polygon(list(unclass(sq[[1]])[[1]], hole))
  holed_rep <- sf::st_polygon(list(unclass(sq[[1]])[[1]],
                                   unclass(with_repeat(sf::st_polygon(list(hole))))[[1]]))
  lay_h  <- sf::st_sf(cell_id = 1:3, geometry = sf::st_sfc(list(holed_rep, sq[[2]], sq[[3]]),
                                                         crs = 32615))
  lay_h0 <- sf::st_sf(cell_id = 1:3, geometry = sf::st_sfc(list(holed, sq[[2]], sq[[3]]),
                                                         crs = 32615))
  expect_true(all(sf::st_is_valid(lay_h)))
  expect_identical(s2_refused(lay_h), 1L)
  expect_equal(.crs_area_error(lay_h), .crs_area_error(lay_h0), tolerance = 1e-12)
  zsq <- function(k) sf::st_polygon(list(cbind(unclass(sq[[k]])[[1]], 0)))
  zrep <- sf::st_polygon(list(cbind(c(4e5, 4.1e5, 4.1e5, 4.1e5, 4e5, 4e5),
                                    c(4.4e6, 4.4e6, 4.4e6, 4.41e6, 4.41e6, 4.4e6),
                                    c(0, 0, 7, 0, 0, 0))))
  lay_z  <- sf::st_sf(cell_id = 1:3, geometry = sf::st_sfc(list(zrep, zsq(2), zsq(3)),
                                                         crs = 32615))
  expect_identical(quiet(s2_refused(lay_z)), 1L)
  e_z <- quiet(.crs_area_error(lay_z))
  expect_true(is.finite(e_z))
  expect_identical(e_z, quiet(.crs_area_error(sf::st_zm(lay_z))))
})

test_that("a far cell that s2 refuses is measured among more than 200", {
  # 299 squares of 10 km in a Lambert conformal conic around 39N and one
  # 2,200 km north of them, valid, with a vertex beside its first corner.
  # The far cell is always measured, since it reaches furthest north; left
  # out, the figure fell from 14.5 percent to 0.13 percent and the cells got
  # a density.
  lcc <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +datum=WGS84 +units=m"
  cen <- sf::st_transform(sf::st_sfc(sf::st_point(c(-96, 39)), crs = 4326), lcc)[[1]]
  gg <- c(lapply(0:298, function(k) square(cen[1] + 2e4 * (k %% 23), cen[2] + 2e4 * (k %/% 23), 1e4)),
          list(square(cen[1], cen[2] + 2.2e6, 1e4)))
  G <- sf::st_sfc(gg, crs = lcc)
  G2 <- G; G2[[300]] <- with_near_corner(G[[300]], lcc)
  full <- sf::st_sf(cell_id = 1:300, geometry = G)
  bad  <- sf::st_sf(cell_id = 1:300, geometry = G2)
  expect_true(all(sf::st_is_valid(bad)))
  expect_identical(s2_refused(bad), 300L)
  e <- .crs_area_error(bad)
  expect_equal(e, .crs_area_error(full), tolerance = 1e-9)
  expect_gt(e, 0.1)
  lines <- capture_spatialkit_log(out <- ensure_projected(bad, purpose = "area"))
  expect_true(log_has(lines, sprintf("distorts areas across this extent by up to %.1f%%", 100 * e)))
  expect_error(quiet(summarize_by_cell(rows_in(full), "z", id_col = "cell_id", cells_sf = bad,
                                       estimate = "mean", area = TRUE)),
               sprintf("distorts areas across the cells by up to %.1f%%", 100 * e), fixed = TRUE)
})

test_that("a cell s2 refuses for another reason leaves the figure NA and the density refused", {
  # Three squares of 10 km in UTM zone 15 and a fourth cell that is not a
  # simple ring: a bowtie, whose planar area from GEOS (2e7 m^2) is neither
  # lobe nor their sum, or a square with a spike.  Neither is valid in
  # EPSG:32615, and s2 refuses each in degrees with or without repeated
  # vertices.  Left out, the bowtie got a density over 2e7 m^2.
  sq <- lapply(0:3, function(k) square(4e5 + 1e4 * k, 4.4e6, 1e4))
  for (odd in list(bowtie(4.3e5, 4.4e6), spike(4.3e5, 4.4e6))) {
    lay <- sf::st_sf(cell_id = 1:4, geometry = sf::st_sfc(c(sq[1:3], list(odd)), crs = 32615))
    expect_identical(sf::st_is_valid(lay), c(TRUE, TRUE, TRUE, FALSE))
    expect_identical(s2_refused(lay), 4L)
    expect_true(is.na(.crs_area_error(lay)))
    expect_true(is.na(spread_each(lay)))
    expect_error(quiet(summarize_by_cell(rows_in(lay), "z", id_col = "cell_id", cells_sf = lay,
                                         estimate = "mean", area = TRUE)),
                 "could not be measured", fixed = TRUE)
  }
  # A refused cell that the repeats are dropped from does not make up for one
  # that s2 refuses for another reason.
  lay <- sf::st_sf(cell_id = 1:4, geometry = sf::st_sfc(
    list(sq[[1]], with_near_corner(sq[[2]], 32615), sq[[3]], bowtie(4.3e5, 4.4e6)), crs = 32615))
  expect_identical(s2_refused(lay), c(2L, 4L))
  expect_true(is.na(.crs_area_error(lay)))
  # One square and a bowtie: one cell measured, no spread.
  lay <- sf::st_sf(cell_id = 1:2, geometry = sf::st_sfc(list(sq[[1]], bowtie(4.3e5, 4.4e6)),
                                                      crs = 32615))
  expect_true(is.na(.crs_area_error(lay)))
})
