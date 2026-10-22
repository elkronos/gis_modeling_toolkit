# tests/testthat/test-fix-cells-cover.R
# ---------------------------------------------------------------------------
# Voronoi cells built on seeds, with no boundary.
#
# The cells used to be clipped to the convex hull of the seeds themselves
# (grown by 2% of its diagonal).  Seeds lie inside the data they stand for,
# so the cells stopped short of it, and assign_features_to_polygons() then
# removed the uncovered observations with one INFO line in the temp log:
# 72% of 3,000 uniform points at 4 seeds, 25% at 25, 4% at 100, and 15% on
# the profile -> selection -> seeds route.  The seeders now attach the ground
# the seeds stand for as attr(, "cover"), the Voronoi builders derive their
# default boundary from it (or from a `cover` argument), and points in no
# cell are a classed warning with a count on the result.
# ---------------------------------------------------------------------------

.fc_obs <- function(n = 600, seed = 2, side = 1e5, crs = 32617) {
  set.seed(seed)
  p <- sf::st_as_sf(data.frame(x = runif(n, 0, side), y = runif(n, 0, side)),
                    coords = c("x", "y"), crs = crs)
  p$z <- rnorm(n)
  p
}

.fc_square <- function(side = 1e5, crs = 32617, x0 = 0, y0 = 0) {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(x0, y0), c(x0 + side, y0), c(x0 + side, y0 + side), c(x0, y0 + side),
    c(x0, y0)))), crs = crs))
}

# Every warning an expression raises, with its classes, and its value.
.fc_conditions <- function(expr) {
  msgs <- character(0); cls <- list(); conds <- list()
  val <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      cls[[length(cls) + 1L]] <<- class(w)
      conds[[length(conds) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = msgs, classes = cls, conditions = conds)
}
.fc_has_class <- function(res, class)
  vapply(res$classes, function(k) class %in% k, logical(1))
# The profile of the two tests below is scored at three levels and its Cp
# optimum is the first, so a count read off it, or off its selection, comes
# with the warning that a bound chose the count (class
# spatialkit_bound_chosen; test-fix-resolution-bound-warning.R).  It says
# nothing about the cover, so it is counted here and set aside by class.
.fc_bound_chosen <- function(expr) {
  n <- 0L
  val <- withCallingHandlers(expr, spatialkit_bound_chosen = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  expect_identical(n, 1L)
  val
}

# The whole route on a layer of seeds: cells, assignment, summary.  Returns
# what the assertions need; raises nothing itself.
.fc_route <- function(seeds, obs, ...) {
  built <- .fc_conditions(build_tessellation(seeds, method = "voronoi",
                                             quiet = TRUE, ...))
  tess <- built$value
  asg <- .fc_conditions(assign_features_to_polygons(
    obs, tess$cells, polygon_id_col = "cell_id"))
  sm <- suppressWarnings(summarize_by_cell(asg$value, "z", id_col = "cell_id"))
  # The k-means partition itself: each observation's nearest seed, measured
  # in the CRS the cells were built in.
  near <- sf::st_nearest_feature(sf::st_transform(obs, sf::st_crs(tess$cells)),
                                 sf::st_transform(seeds, sf::st_crs(tess$cells)))
  list(tess = tess, built = built, asg = asg, sum_n = sum(sm$n),
       nearest_cell = tess$index[near])
}

.fc_expect_all_assigned <- function(r, obs, info) {
  expect_length(r$built$warnings, 0L)
  expect_length(r$asg$warnings, 0L)
  expect_identical(r$tess$params$n_outside, 0L, info = info)
  expect_identical(nrow(r$asg$value), nrow(obs), info = info)
  expect_identical(attr(r$asg$value, "unassigned")$n, 0L, info = info)
  expect_identical(r$sum_n, nrow(obs), info = info)
  # A Voronoi cell is the ground nearest its seed, so the polygon assignment
  # is the nearest-seed assignment.
  expect_identical(r$asg$value$cell_id, r$nearest_cell, info = info)
}


# --- the seeders attach the ground the seeds stand for -----------------------

test_that("the seeding functions attach a cover", {
  obs <- .fc_obs(200)
  bnd <- .fc_square()

  # k-means on points: the points' convex hull, in their CRS.
  km <- voronoi_seeds_kmeans(obs, k = 6)
  cov <- attr(km, "cover")
  expect_s3_class(cov, "sfc")
  expect_true(sf::st_crs(cov) == sf::st_crs(obs))
  expect_equal(as.numeric(sf::st_area(cov)),
               as.numeric(sf::st_area(sf::st_convex_hull(sf::st_union(obs)))))
  expect_true(all(lengths(sf::st_intersects(obs, cov)) == 1L))

  # Random seeds: the boundary they were drawn in.
  rn <- voronoi_seeds_random(bnd, k = 6, set_seed = 1)
  expect_equal(as.numeric(sf::st_area(attr(rn, "cover"))), 1e10)

  # get_voronoi_seeds(): the boundary when there is one, else the cloud's hull.
  g1 <- get_voronoi_seeds(bnd, method = "random", n = 6, set_seed = 1)
  expect_equal(as.numeric(sf::st_area(attr(g1, "cover"))), 1e10)
  g2 <- get_voronoi_seeds(bnd, method = "kmeans", n = 6, set_seed = 1)
  expect_equal(as.numeric(sf::st_area(attr(g2, "cover"))), 1e10)
  g3 <- get_voronoi_seeds(method = "kmeans", n = 6, sample_points = obs, set_seed = 1)
  expect_equal(as.numeric(sf::st_area(attr(g3, "cover"))),
               as.numeric(sf::st_area(sf::st_convex_hull(sf::st_union(obs)))))
  g4 <- get_voronoi_seeds(bnd, method = "kmeans", n = 6, sample_points = obs,
                          set_seed = 1)
  expect_equal(as.numeric(sf::st_area(attr(g4, "cover"))), 1e10)
  # Provided seeds have a cover only when the call says what they stand for.
  mine <- obs[1:5, ]
  expect_null(attr(get_voronoi_seeds(method = "provided", seeds = mine), "cover"))
  expect_false(is.null(attr(get_voronoi_seeds(bnd, method = "provided", seeds = mine),
                            "cover")))
})

test_that("the cover is kept across the CRS alignment of the seeds", {
  obs <- .fc_obs(200)
  bnd_ll <- sf::st_transform(.fc_square(), 4326)
  # The cloud is projected, the boundary lon/lat: the seeds are moved into
  # the boundary's CRS, and the cover goes with them.
  s <- get_voronoi_seeds(bnd_ll, method = "kmeans", n = 5, sample_points = obs,
                         set_seed = 1)
  expect_true(sf::st_crs(s) == sf::st_crs(4326))
  expect_s3_class(attr(s, "cover"), "sfc")
  # And across what a caller does to them, short of rbind().
  km <- voronoi_seeds_kmeans(obs, k = 5)
  expect_false(is.null(attr(sf::st_transform(km, 4326), "cover")))
  expect_false(is.null(attr(km[1:3, ], "cover")))
  expect_null(attr(rbind(km, km), "cover"))
})


# --- every observation gets a cell -------------------------------------------

test_that("k-means seeds with no boundary: every observation is assigned (projected)", {
  obs <- .fc_obs(600)
  for (k in c(1, 2, 4, 25, 100)) {
    seeds <- suppressWarnings(voronoi_seeds_kmeans(obs, k = k))
    r <- .fc_route(seeds, obs)
    expect_identical(nrow(r$tess$cells), as.integer(k))
    expect_identical(r$tess$params$boundary_source, "seeds' cover")
    .fc_expect_all_assigned(r, obs, info = paste("k =", k))
  }
})

test_that("k-means seeds with no boundary: every observation is assigned (lon/lat)", {
  # The seeds' working CRS is picked from the seeds' own extent, so the cover
  # has to follow them out of the observations' CRS.
  obs <- sf::st_transform(.fc_obs(600), 4326)
  for (k in c(1, 2, 4, 25, 100)) {
    seeds <- suppressWarnings(voronoi_seeds_kmeans(obs, k = k))
    r <- .fc_route(seeds, obs)
    .fc_expect_all_assigned(r, obs, info = paste("lon/lat, k =", k))
  }
})

test_that("k-means seeds on a real lon/lat layer: every feature is assigned", {
  shp <- system.file("shape/nc.shp", package = "sf")
  skip_if(!nzchar(shp), "sf's nc.shp is not available")
  nc <- suppressWarnings(sf::st_centroid(sf::st_geometry(sf::st_read(shp, quiet = TRUE))))
  obs <- sf::st_sf(z = seq_along(nc) / 10, geometry = nc)
  for (k in c(3, 10, 30)) {
    seeds <- voronoi_seeds_kmeans(obs, k = k)
    r <- .fc_route(seeds, obs)
    .fc_expect_all_assigned(r, obs, info = paste("nc, k =", k))
  }
})


# --- a continent-sized lon/lat layer -----------------------------------------

# Points uniform over a box in longitude and latitude, EPSG:4326.
.fc_wide_ll <- function(lon0, lat0, w, h, n, seed) {
  set.seed(seed)
  d <- data.frame(lon = lon0 + runif(n, 0, w), lat = lat0 + runif(n, 0, h))
  d$z <- sin(d$lon * 3) + cos(d$lat * 2) + rnorm(n)
  sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
}

test_that("the cover of a lon/lat cloud is its hull in degrees, whatever sf_use_s2() says", {
  # The cover is moved into the cells' CRS with its edges read as straight
  # lines in degrees (.cover_follow()).  It was the hull s2 builds, whose
  # edges are great circles: over the box of the conterminous United States
  # 24 of 800 points lay outside it read that way, one of them by 130 km
  # in the cells' CRS, more than the 2% margin.
  obs <- .fc_wide_ll(-125, 25, 58, 24, 800, 4)
  cov <- attr(suppressWarnings(voronoi_seeds_kmeans(obs, k = 44, set_seed = 1)), "cover")
  expect_true(sf::st_crs(cov) == sf::st_crs(4326))
  bare <- function(x) sf::st_set_crs(sf::st_geometry(x), NA)
  expect_true(all(lengths(sf::st_intersects(bare(obs), bare(cov))) == 1L))
  # The planar hull of the coordinates, whether or not sf uses s2.
  hull <- sf::st_convex_hull(sf::st_combine(bare(obs)))
  expect_equal(as.numeric(sf::st_area(bare(cov))), as.numeric(sf::st_area(hull)))
  was <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(was)), add = TRUE)
  expect_identical(.cover_of(obs), cov)
})

test_that("k-means seeds on a continent-sized lon/lat layer: every observation is assigned", {
  # Over the box of the conterminous United States one point, (-75.49,
  # 48.98), fell in no cell at 44 and at 178 seeds; over northern Asia 13 of
  # 300 points did at 10 seeds.  Around the North Pole s2's hull, read in
  # degrees, held 19 of 600 points and none of the seeds, so it was set
  # aside (spatialkit_no_cover) and 230 points fell in no cell.
  cases <- list(list(obs = .fc_wide_ll(-125, 25, 58, 24, 800, 4), k = c(44, 178)),
                list(obs = .fc_wide_ll(60, 55, 120, 25, 300, 3), k = 10),
                list(obs = .fc_wide_ll(-180, 70, 360, 19.5, 600, 1), k = 10))
  for (cs in cases) for (k in cs$k) {
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = k,
                                                sample_points = cs$obs, set_seed = 1))
    r <- .fc_route(seeds, cs$obs)
    expect_identical(r$tess$params$boundary_source, "seeds' cover")
    .fc_expect_all_assigned(r, cs$obs, info = paste("k =", k))
  }
})

test_that("a continent-sized lon/lat layer: the summary's check reads match on the default path", {
  obs <- .fc_wide_ll(-125, 25, 58, 24, 800, 4)
  prof  <- suppressWarnings(suppressMessages(resolution_profile(obs, "z", seed = 123)))
  sel   <- suppressWarnings(select_resolution(prof))
  seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = obs))
  tess  <- suppressMessages(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
  asg <- .fc_conditions(assign_features_to_polygons(obs, tess$cells, polygon_id_col = "cell_id"))
  expect_false(any(.fc_has_class(asg, "spatialkit_unassigned")))
  expect_identical(nrow(asg$value), nrow(obs))
  out <- .fc_conditions(summarize_by_cell(asg$value, "z", id_col = "cell_id",
                                          cells_sf = tess$cells))
  expect_false(any(.fc_has_class(out, "spatialkit_cells_differ")))
  chk <- attr(out$value, "count_rule")$check
  expect_identical(chk$cells, "match")
  expect_identical(chk$estimate, "match")
})

test_that("the cover of a cloud in another geographic CRS is its hull in degrees too", {
  # The same 800 points in NAD83 and in ETRS89: s2's hull of either, read in
  # degrees, left the point near the northern edge in no cell at 44 and at
  # 178 seeds, as in WGS 84.
  for (crs in c(4269, 4258)) {
    obs <- sf::st_set_crs(sf::st_set_crs(.fc_wide_ll(-125, 25, 58, 24, 800, 4), NA), crs)
    bare <- function(x) sf::st_set_crs(sf::st_geometry(x), NA)
    hull <- sf::st_convex_hull(sf::st_combine(bare(obs)))
    for (k in c(44, 178)) {
      info <- sprintf("EPSG:%d, k = %d", crs, k)
      seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = k,
                                                  sample_points = obs, set_seed = 1))
      cov <- attr(seeds, "cover")
      expect_true(sf::st_crs(cov) == sf::st_crs(crs), info = info)
      expect_equal(as.numeric(sf::st_area(bare(cov))), as.numeric(sf::st_area(hull)), info = info)
      r <- .fc_route(seeds, obs)
      expect_identical(r$tess$params$boundary_source, "seeds' cover", info = info)
      .fc_expect_all_assigned(r, obs, info = info)
    }
  }
})


# --- a layer across 180 degrees ----------------------------------------------

# Points uniform over a box that starts at lon0 and runs w degrees east, with
# their longitudes put back in [-180, 180), EPSG:4326.
.fc_across_180 <- function(lon0, lat0, w, h, n, seed) {
  set.seed(seed)
  d <- data.frame(lon = ((lon0 + runif(n, 0, w) + 180) %% 360) - 180,
                  lat = lat0 + runif(n, 0, h))
  d$z <- sin(d$lon * 3) + cos(d$lat * 2) + rnorm(n)
  sf::st_as_sf(d, coords = c("lon", "lat"), crs = 4326)
}

test_that("the cover of points on both sides of 180 degrees is their hull with longitudes from 0 to 360", {
  # Taken on the longitudes as given, from -180 to 180, the hull of these
  # points ran the long way round the globe.
  obs <- .fc_across_180(170, -45, 20, 90, 600, 1)
  lon <- sf::st_coordinates(obs)[, 1L]
  expect_true(any(lon < 0) && any(lon > 0))
  cov <- attr(suppressWarnings(voronoi_seeds_kmeans(obs, k = 10, set_seed = 1)), "cover")
  expect_true(sf::st_crs(cov) == sf::st_crs(4326))
  bb <- sf::st_bbox(cov)
  expect_gte(bb[["xmin"]], 170)
  expect_gt(bb[["xmax"]], 180)
  expect_lte(bb[["xmax"]], 190)
  # It holds every point once the longitudes west of 180 are counted past it.
  xy <- sf::st_coordinates(obs)
  xy[xy[, 1L] < 0, 1L] <- xy[xy[, 1L] < 0, 1L] + 360
  past <- sf::st_sfc(lapply(seq_len(nrow(xy)), function(i) sf::st_point(xy[i, 1:2])))
  expect_true(all(lengths(sf::st_intersects(past, sf::st_set_crs(cov, NA))) == 1L))
  # Longitudes that span no more as given stay as given: 0 to 360 over the
  # Pacific, and the conterminous United States.  Others are taken from -180
  # to 180 or from 0 to 360, whichever spans fewer degrees.
  expect_identical(.lon_narrowest(c(160, 170, 200, 210)), c(160, 170, 200, 210))
  expect_identical(.lon_narrowest(c(-125, -100, -67)), c(-125, -100, -67))
  expect_identical(.lon_narrowest(c(170, 175, -175, -170)), c(170, 175, 185, 190))
  expect_equal(.lon_narrowest(c(-100, -90, 30, 40, 170, -175)), c(260, 270, 30, 40, 170, 185))
  expect_equal(.lon_narrowest(c(350, 355, 0, 5, 10)), c(-10, -5, 0, 5, 10))
  expect_identical(.lon_narrowest(c(5, 5, 5)), c(5, 5, 5))
})

test_that("a cover past 180 degrees goes into Web Mercator in two pieces, and into a CRS centred on it in one", {
  cov <- sf::st_sfc(sf::st_polygon(list(rbind(c(170, -45), c(190, -45), c(190, 45),
                                              c(170, 45), c(170, -45)))), crs = 4326)
  pts <- sf::st_sfc(sf::st_point(c(175, 0)), sf::st_point(c(-175, 0)), sf::st_point(c(-171, 40)),
                    sf::st_point(c(171, -40)), crs = 4326)
  # Web Mercator splits the map at 180: the part past it goes to the other
  # side, where the points west of 180 are.
  to <- sf::st_transform(pts, 3857)
  merc <- .cover_follow(cov, pts, to)
  expect_true(sf::st_crs(merc) == sf::st_crs(3857))
  expect_length(sf::st_cast(merc, "POLYGON"), 2L)
  expect_true(all(lengths(sf::st_intersects(to, merc)) == 1L))
  # A projection centred on 180 joins the two sides, and the cover moves as
  # it is.
  to <- sf::st_transform(pts, "+proj=laea +lat_0=0 +lon_0=180 +datum=WGS84 +units=m +no_defs")
  cen <- .cover_follow(cov, pts, to)
  expect_identical(as.character(sf::st_geometry_type(cen)), "POLYGON")
  expect_true(all(lengths(sf::st_intersects(to, cen)) == 1L))
  expect_true(.joins_across_180(sf::st_crs(4326), sf::st_crs(to), c(-45, 45)))
  expect_false(.joins_across_180(sf::st_crs(4326), sf::st_crs(3857), c(-45, 45)))
  # In Mollweide and Equal Earth the 180th meridian is a curve.  The edges
  # the cut makes along it get vertices too, so they follow it: with only
  # their two ends they were chords, the two pieces crossed, and none of
  # these points near 180 was inside.
  near <- sf::st_sfc(sf::st_point(c(179.5, 0)), sf::st_point(c(-179.5, 0)),
                     sf::st_point(c(179.9, 20)), sf::st_point(c(-179.9, -20)), crs = 4326)
  for (crs in c("+proj=moll +datum=WGS84 +units=m +no_defs",
                "+proj=eqearth +datum=WGS84 +units=m +no_defs")) {
    to <- sf::st_transform(near, crs)
    curved <- .cover_follow(cov, near, to)
    expect_false(.joins_across_180(sf::st_crs(4326), sf::st_crs(crs), c(-45, 45)), info = crs)
    expect_length(sf::st_cast(curved, "POLYGON"), 2L)
    expect_true(all(sf::st_is_valid(curved)), info = crs)
    expect_true(all(lengths(sf::st_intersects(to, curved)) == 1L), info = crs)
  }
})

test_that("k-means seeds on layers across 180 degrees: every observation is assigned", {
  # With the cover taken the long way round, 32 of these 600 points fell in
  # no cell at 10 seeds between 170 E and 170 W, and 190 and 50 at 10 and 44
  # seeds between 90 E and 100 W.  Over 190 degrees of longitude, whose seeds
  # span less than 180 and get a projection centred on them, the cover held
  # none of the seeds and was set aside, and 259 fell in no cell.  Over 250
  # degrees the seeds' CRS is Web Mercator, and the cover goes in two pieces.
  # Built in Mollweide (`crs`), whose 180th meridian is a curve, the cover in
  # two pieces left 6 points between 170 E and 170 W in no cell while its
  # edges along 180 had no vertices but their ends.
  moll <- "+proj=moll +datum=WGS84 +units=m +no_defs"
  cases <- list(list(box = c(170, -45, 20, 90), k = 10, ratio = TRUE),
                list(box = c(90, -60, 170, 120), k = c(10, 44)),
                list(box = c(85, -60, 190, 120), k = 10),
                list(box = c(55, -60, 250, 120), k = 10, merc = TRUE),
                list(box = c(170, -45, 20, 90), k = 10, crs = moll))
  for (cs in cases) for (k in cs$k) {
    cs$obs <- .fc_across_180(cs$box[1L], cs$box[2L], cs$box[3L], cs$box[4L], 600, 1)
    info <- sprintf("%g degrees east from %g E, k = %d%s", cs$box[3L], cs$box[1L], k,
                    if (is.null(cs$crs)) "" else ", in Mollweide")
    seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = k,
                                                sample_points = cs$obs, set_seed = 1))
    r <- if (is.null(cs$crs)) .fc_route(seeds, cs$obs) else
      .fc_route(seeds, cs$obs, crs = sf::st_crs(cs$crs))
    expect_identical(r$tess$params$boundary_source, "seeds' cover", info = info)
    .fc_expect_all_assigned(r, cs$obs, info = info)
    if (isTRUE(cs$merc))
      expect_true(sf::st_crs(r$tess$cells) == sf::st_crs(3857), info = info)
    if (!is.null(cs$crs))
      expect_true(sf::st_crs(r$tess$cells) == sf::st_crs(cs$crs), info = info)
    if (isTRUE(cs$ratio)) {
      # The cells hold the points' own hull and little more: the long way
      # round they were more than ten times its area.
      hull <- sf::st_convex_hull(sf::st_union(sf::st_transform(sf::st_geometry(cs$obs),
                                                               sf::st_crs(r$tess$cells))))
      expect_lt(as.numeric(sum(sf::st_area(r$tess$cells)) / sf::st_area(hull)), 1.5)
    }
  }
})

test_that("a layer across 180 degrees: the summary's check reads match on the default path", {
  # 300 points between 90 E and 100 W and 60 S and 60 N: the cover taken the
  # long way round left 14 rows in no cell, and the check read "differ".
  obs <- .fc_across_180(90, -60, 170, 120, 300, 1)
  prof  <- suppressWarnings(suppressMessages(resolution_profile(obs, "z", seed = 123)))
  sel   <- suppressWarnings(select_resolution(prof))
  seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = obs))
  tess  <- suppressMessages(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
  asg <- .fc_conditions(assign_features_to_polygons(obs, tess$cells, polygon_id_col = "cell_id"))
  expect_false(any(.fc_has_class(asg, "spatialkit_unassigned")))
  expect_identical(nrow(asg$value), nrow(obs))
  out <- .fc_conditions(summarize_by_cell(asg$value, "z", id_col = "cell_id",
                                          cells_sf = tess$cells))
  expect_false(any(.fc_has_class(out, "spatialkit_cells_differ")))
  expect_identical(attr(out$value, "count_rule")$check$cells, "match")
})

test_that("observations with no CRS are covered too, in degrees and in metres", {
  # Neither layer has a CRS: seeds and cover are read in one coordinate space,
  # whichever way the builder then interprets the seeds.
  planar <- sf::st_set_crs(.fc_obs(300), NA)
  seeds <- voronoi_seeds_kmeans(planar, k = 9)
  tess <- suppressMessages(build_tessellation(seeds, method = "voronoi", quiet = TRUE))
  expect_identical(tess$params$boundary_source, "seeds' cover")
  expect_true(all(lengths(sf::st_intersects(planar, tess$cells)) >= 1L))

  set.seed(3)
  deg <- sf::st_as_sf(data.frame(x = runif(300, 10, 11), y = runif(300, 50, 51)),
                      coords = c("x", "y"))
  seeds_d <- suppressWarnings(voronoi_seeds_kmeans(deg, k = 9))
  # CRS-less degrees are taken as EPSG:4326 and projected, with a warning;
  # the cover is taken the same way.
  tess_d <- suppressWarnings(suppressMessages(
    build_tessellation(seeds_d, method = "voronoi", quiet = TRUE)))
  deg_p <- sf::st_transform(sf::st_set_crs(deg, 4326), sf::st_crs(tess_d$cells))
  expect_true(all(lengths(sf::st_intersects(deg_p, tess_d$cells)) >= 1L))
})

test_that("seeds drawn in a boundary fill it when the boundary is not passed on", {
  bnd <- .fc_square()
  obs <- .fc_obs(600)
  for (m in c("random", "kmeans")) {
    seeds <- get_voronoi_seeds(bnd, method = m, n = 10, set_seed = 3)
    r <- .fc_route(seeds, obs)
    .fc_expect_all_assigned(r, obs, info = m)
    covered <- sum(as.numeric(sf::st_area(suppressWarnings(
      sf::st_intersection(sf::st_geometry(r$tess$cells), sf::st_geometry(bnd)))))) / 1e10
    # 36% for 10 random seeds and 54% for 10 k-means seeds before.
    expect_gt(covered, 0.9999)
  }
  r <- .fc_route(voronoi_seeds_random(bnd, k = 10, set_seed = 3), obs)
  .fc_expect_all_assigned(r, obs, info = "voronoi_seeds_random")
})

test_that("seeds from a resolution profile or selection hold every observation", {
  skip_if_not_installed("gstat")
  set.seed(11)
  obs <- sf::st_as_sf(data.frame(x = runif(400, 0, 5e4), y = runif(400, 0, 5e4)),
                      coords = c("x", "y"), crs = 32617)
  obs$z <- rnorm(400) + sin(sf::st_coordinates(obs)[, 1] / 8000)
  prof <- suppressWarnings(suppressMessages(
    resolution_profile(obs, "z", levels = c(4, 9, 16), nstart = 5, estimate = NULL,
                       obs_per_cell = NULL, sac = NULL, partition = "kmeans")))
  # The profile keeps the hull of the points it was computed on.
  expect_s3_class(attr(prof, "hull"), "sfc")
  sel <- select_resolution(prof, criterion = "cp")

  # The route the package's warnings prescribe: 15% of rows lost before.
  s1 <- .fc_bound_chosen(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = obs,
                                           partition = "kmeans"))
  .fc_expect_all_assigned(.fc_route(s1, obs), obs, info = "selection + sample_points")
  s2 <- .fc_bound_chosen(voronoi_seeds_kmeans(obs, k = sel, partition = "kmeans"))
  .fc_expect_all_assigned(.fc_route(s2, obs), obs, info = "voronoi_seeds_kmeans(k = selection)")
  # The profile alone: its own hull is the cover.
  s3 <- .fc_bound_chosen(get_voronoi_seeds(method = "kmeans", n = prof, partition = "kmeans"))
  expect_false(is.null(attr(s3, "cover")))
  .fc_expect_all_assigned(.fc_route(s3, obs), obs, info = "profile alone")

  # A selection whose centres carry the profile's hull, given with neither
  # points nor a boundary: get_voronoi_seeds() hands the hull on as the cover.
  # (The hull is put on the centres here; see the next test for what a
  # select_resolution() result carries by itself.)
  sel_h <- sel
  attr(sel_h$seeds, "cover") <- attr(prof, "hull")
  s4 <- .fc_bound_chosen(get_voronoi_seeds(method = "kmeans", n = sel_h, partition = "kmeans"))
  expect_false(is.null(attr(s4, "cover")))
  expect_equal(as.numeric(sf::st_area(attr(s4, "cover"))),
               as.numeric(sf::st_area(attr(prof, "hull"))))
  r4 <- .fc_route(s4, obs)
  expect_identical(r4$tess$params$boundary_source, "seeds' cover")
  .fc_expect_all_assigned(r4, obs, info = "selection carrying the profile's hull")

  # Centres with nothing to say what ground they stand for: the cells stop
  # at the centres' own hull, and that is said, at the build and at the
  # assignment.
  sel_0 <- sel
  attr(sel_0$seeds, "cover") <- NULL
  s5 <- .fc_bound_chosen(get_voronoi_seeds(method = "kmeans", n = sel_0, partition = "kmeans"))
  expect_null(attr(s5, "cover"))
  r5 <- .fc_route(s5, obs)
  expect_true(any(.fc_has_class(r5$built, "spatialkit_no_cover")))
  expect_match(r5$built$warnings[.fc_has_class(r5$built, "spatialkit_no_cover")][1],
               "layer of seeds .* no `boundary` and no cover")
  expect_identical(r5$tess$params$boundary_source, "points")
  expect_true(any(.fc_has_class(r5$asg, "spatialkit_unassigned")))
  # Passing the observations as `cover` repairs it.
  r6 <- .fc_route(s5, obs, cover = obs)
  expect_identical(r6$tess$params$boundary_source, "cover")
  .fc_expect_all_assigned(r6, obs, info = "no cover on the centres, cover = observations")
})

test_that("a select_resolution() result alone holds every observation", {
  skip_if_not_installed("gstat")
  set.seed(11)
  obs <- sf::st_as_sf(data.frame(x = runif(400, 0, 5e4), y = runif(400, 0, 5e4)),
                      coords = c("x", "y"), crs = 32617)
  obs$z <- rnorm(400) + sin(sf::st_coordinates(obs)[, 1] / 8000)
  prof <- suppressWarnings(suppressMessages(
    resolution_profile(obs, "z", levels = c(4, 9, 16), nstart = 5, estimate = NULL,
                       obs_per_cell = NULL, sac = NULL, partition = "kmeans")))
  sel <- select_resolution(prof, criterion = "cp")
  # select_resolution() builds its centres with .profile_seeds() (in
  # R/resolution-profile.R), which puts the profile's hull on them as their
  # cover.  Without it a selection given with neither points nor a boundary
  # had no cover to hand on and its cells stopped at the centres' hull (with
  # two classed warnings; previous test).
  expect_s3_class(attr(sel$seeds, "cover", exact = TRUE), "sfc")
  s <- .fc_bound_chosen(get_voronoi_seeds(method = "kmeans", n = sel, partition = "kmeans"))
  expect_false(is.null(attr(s, "cover")))
  r <- .fc_route(s, obs)
  expect_identical(r$tess$params$boundary_source, "seeds' cover")
  .fc_expect_all_assigned(r, obs, info = "selection alone")
})


# --- the `cover` argument ----------------------------------------------------

test_that("`cover` sets the default boundary: points, a polygon, another CRS, no CRS", {
  obs <- .fc_obs(400)
  seeds <- voronoi_seeds_kmeans(obs, k = 8)
  bare <- seeds
  attr(bare, "cover") <- NULL                       # as rbind() leaves them

  # No cover, seeds: the old cells, and a classed warning.
  res <- .fc_conditions(build_tessellation(bare, method = "voronoi", quiet = TRUE))
  expect_true(any(.fc_has_class(res, "spatialkit_no_cover")))
  expect_gt(sum(lengths(sf::st_intersects(obs, res$value$cells)) == 0L), 50L)

  inside <- function(tess) all(lengths(sf::st_intersects(
    sf::st_transform(obs, sf::st_crs(tess$cells)), tess$cells)) >= 1L)
  # The observations.
  expect_no_warning(t1 <- build_tessellation(bare, method = "voronoi", quiet = TRUE,
                                             cover = obs))
  expect_true(inside(t1))
  expect_identical(t1$params$boundary_source, "cover")
  expect_null(t1$params$cover_outside)
  # The boundary is the observations' hull grown by 2% of its diagonal.
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(obs)))
  bb <- sf::st_bbox(hull)
  grown <- sf::st_buffer(hull, 0.02 * sqrt((bb[["xmax"]] - bb[["xmin"]])^2 +
                                             (bb[["ymax"]] - bb[["ymin"]])^2))
  expect_equal(as.numeric(sf::st_area(t1$boundary)), as.numeric(sf::st_area(grown)),
               tolerance = 1e-8)
  expect_equal(sum(as.numeric(sf::st_area(t1$cells))), as.numeric(sf::st_area(grown)),
               tolerance = 1e-8)
  # A polygon.
  expect_no_warning(t2 <- build_tessellation(bare, method = "voronoi", quiet = TRUE,
                                             cover = .fc_square()))
  expect_true(inside(t2))
  # The observations in another CRS.
  expect_no_warning(t3 <- build_tessellation(bare, method = "voronoi", quiet = TRUE,
                                             cover = sf::st_transform(obs, 4326)))
  expect_true(inside(t3))
  # The observations with no CRS: read in the seeds' coordinates.
  expect_no_warning(t4 <- build_tessellation(bare, method = "voronoi", quiet = TRUE,
                                             cover = sf::st_set_crs(obs, NA)))
  expect_true(inside(t4))
  # create_voronoi_polygons() takes it too.
  expect_no_warning(t5 <- create_voronoi_polygons(bare, cover = obs, quiet = TRUE))
  expect_true(inside(t5))
  expect_identical(t5$params$boundary_source, "cover")
  # The argument wins over the cover the seeds carry.
  half <- .fc_square(side = 5e4)
  t6 <- suppressWarnings(build_tessellation(seeds, method = "voronoi", quiet = TRUE,
                                            cover = sf::st_geometry(half)))
  expect_identical(t6$params$boundary_source, "cover")

  expect_error(build_tessellation(bare, method = "voronoi", cover = "obs"),
               "`cover` must be an sf or sfc point or polygon layer")
})

test_that("a generator outside the cover keeps its ground, and a cover elsewhere is set aside", {
  obs <- .fc_obs(300)
  seeds <- voronoi_seeds_kmeans(obs, k = 8)
  attr(seeds, "cover") <- NULL
  # A cover over the western half: the eastern seeds still get cells.
  west <- .fc_square(side = 5e4)
  t1 <- suppressWarnings(build_tessellation(seeds, method = "voronoi", quiet = TRUE,
                                            cover = west))
  expect_identical(t1$params$n_outside, 0L)
  expect_identical(nrow(t1$cells), 8L)
  # A cover that holds none of them is in another place: refused, by name.
  far <- .fc_square(x0 = 9e6)
  res <- .fc_conditions(build_tessellation(seeds, method = "voronoi", quiet = TRUE,
                                           cover = far))
  hit <- .fc_has_class(res, "spatialkit_no_cover")
  expect_identical(sum(hit), 1L)
  expect_match(res$warnings[hit], "`cover` holds none of the points")
  expect_identical(res$value$params$boundary_source, "points")
})


# --- points in no cell are counted, and warned about -------------------------

test_that("a boundary smaller than the data: the build and the assignment both warn, with counts", {
  obs <- .fc_obs(600)
  seeds <- voronoi_seeds_kmeans(obs, k = 25)
  small <- .fc_square(side = 5e4)
  n_out <- sum(lengths(sf::st_intersects(obs, small)) == 0L)
  expect_gt(n_out, 100L)

  # With the observations as `cover`, build_tessellation() counts the ones
  # the boundary leaves out.
  built <- .fc_conditions(build_tessellation(seeds, boundary = small,
                                             method = "voronoi", quiet = TRUE,
                                             cover = obs))
  hit <- .fc_has_class(built, "spatialkit_unassigned")
  expect_true(any(grepl("`cover` point", built$warnings[hit])))
  expect_identical(built$value$params$cover_outside, n_out)
  expect_identical(built$value$params$boundary_source, "boundary")
  # The seeds whose cells lie wholly outside the boundary are NA in `index`,
  # counted and warned about, as seeds.
  expect_identical(built$value$params$n_outside, sum(is.na(built$value$index)))
  expect_gt(built$value$params$n_outside, 0L)
  expect_true(any(grepl("NA in `index`", built$warnings[hit])))
  expect_true(any(grepl("seed(s)", built$warnings[hit], fixed = TRUE)))

  # assign_features_to_polygons(): an R warning with a class, not a log line.
  asg <- .fc_conditions(assign_features_to_polygons(obs, built$value$cells,
                                                    polygon_id_col = "cell_id"))
  hit_a <- .fc_has_class(asg, "spatialkit_unassigned")
  expect_identical(sum(hit_a), 1L)
  expect_match(asg$warnings[hit_a], sprintf("%d of 600 feature\\(s\\)", n_out))
  expect_identical(asg$conditions[[which(hit_a)]]$n, n_out)
  expect_identical(asg$conditions[[which(hit_a)]]$total, 600L)
  rec <- attr(asg$value, "unassigned")
  expect_identical(rec$n, n_out)
  expect_identical(rec$which, which(lengths(sf::st_intersects(obs, small)) == 0L))
  expect_identical(nrow(asg$value), 600L - n_out)
})

test_that("build_tessellation() counts and warns about points outside every cell", {
  obs <- .fc_obs(300)
  small <- .fc_square(side = 5e4)
  n_out <- sum(lengths(sf::st_intersects(obs, small)) == 0L)
  for (m in c("voronoi", "square", "hex")) {
    res <- .fc_conditions(build_tessellation(obs, boundary = small, method = m,
                                             approx_n_cells = if (m == "voronoi") NULL else 9,
                                             quiet = TRUE))
    hit <- .fc_has_class(res, "spatialkit_unassigned")
    expect_identical(sum(hit), 1L, info = m)
    expect_match(res$warnings[hit], sprintf("%d of 300 point\\(s\\)", n_out), info = m)
    expect_identical(res$value$params$n_outside, n_out, info = m)
    expect_identical(sum(is.na(res$value$index)), n_out, info = m)
    # Caught by class, without matching the message.
    caught <- tryCatch(
      build_tessellation(obs, boundary = small, method = m, quiet = TRUE,
                         approx_n_cells = if (m == "voronoi") NULL else 9),
      spatialkit_unassigned = function(w) w$n)
    expect_identical(caught, n_out, info = m)
  }
  # A lattice checks a point `cover` against its cells as well.
  seeds <- voronoi_seeds_kmeans(obs, k = 5)
  g <- .fc_conditions(build_tessellation(obs[lengths(sf::st_intersects(obs, small)) == 1L, ],
                                         boundary = small, method = "square",
                                         approx_n_cells = 9, quiet = TRUE, cover = obs))
  expect_identical(g$value$params$cover_outside, n_out)
  expect_identical(g$value$params$n_outside, 0L)
  expect_true(any(.fc_has_class(g, "spatialkit_unassigned")))
})


# --- the messages name what the code did -------------------------------------

test_that("the cover check says what set the cells, for the method that ran", {
  obs <- .fc_obs(300)
  small <- .fc_square(side = 5e4)
  inside <- obs[lengths(sf::st_intersects(obs, small)) == 1L, ]
  msg_of <- function(...) {
    res <- .fc_conditions(build_tessellation(..., quiet = TRUE, cover = obs))
    hit <- .fc_has_class(res, "spatialkit_unassigned") & grepl("`cover` point", res$warnings)
    expect_identical(sum(hit), 1L)
    res$warnings[hit]
  }
  # Voronoi cells in a boundary: the boundary set them, and leaving it out
  # would have the cells built to hold the cover.
  v <- msg_of(voronoi_seeds_kmeans(obs, k = 9), boundary = small, method = "voronoi")
  expect_match(v, "`boundary` sets the cells, and it does not hold them", fixed = TRUE)
  expect_match(v, "leave `boundary` out so the cells are built to hold `cover`", fixed = TRUE)
  # A grid is laid over the boundary whatever the cover: that advice was
  # given here too, and a grid cannot be built without a boundary.
  for (m in c("square", "hex")) {
    g <- msg_of(inside, boundary = small, method = m, approx_n_cells = 9)
    expect_match(g, "the grid is laid over `boundary`, and it does not hold them", fixed = TRUE)
    expect_match(g, "`cover` sizes nothing for a grid: pass a boundary that holds the observations",
                 fixed = TRUE)
    expect_false(grepl("leave `boundary` out", g, fixed = TRUE), info = m)
  }
  # Triangles span the hull of the points triangulated.  With no boundary
  # passed at all the message used to blame `boundary` and say to leave it out.
  t0 <- msg_of(inside, method = "triangles")
  expect_match(t0, "the triangles span the convex hull of `points_sf`, and it does not hold them",
               fixed = TRUE)
  expect_match(t0, "`cover` sizes nothing for method = \"triangles\"", fixed = TRUE)
  expect_false(grepl("`boundary` sets the cells", t0, fixed = TRUE))
  expect_false(grepl("leave `boundary` out", t0, fixed = TRUE))
  t1 <- msg_of(inside, boundary = small, method = "triangles")
  expect_match(t1, "convex hull of `points_sf`, clipped to `boundary`, and that does not hold them",
               fixed = TRUE)
  # What the triangle message points to does hold the cover.
  fixed <- build_tessellation(voronoi_seeds_kmeans(inside, k = 9), method = "voronoi",
                              quiet = TRUE, cover = obs)
  expect_true(all(lengths(sf::st_intersects(obs, fixed$cells)) >= 1L))
})

test_that("a cover the seeds carry is called that, not `cover`, through build_tessellation()", {
  obs <- .fc_obs(300)
  seeds <- voronoi_seeds_kmeans(obs, k = 8)
  # The message: no `cover` argument was passed.
  said <- character(0)
  tess <- withCallingHandlers(
    build_tessellation(seeds, method = "voronoi"),
    message = function(m) { said <<- c(said, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("deriving boundary from the cover the seeds carry", said, fixed = TRUE)))
  expect_false(any(grepl("from `cover`", said, fixed = TRUE)))
  expect_identical(tess$params$boundary_source, "seeds' cover")
  expect_null(tess$params$cover_outside)
  # The direct call says the same.
  said2 <- character(0)
  withCallingHandlers(
    create_voronoi_polygons(seeds),
    message = function(m) { said2 <<- c(said2, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("deriving boundary from the cover the seeds carry", said2, fixed = TRUE)))
  # A `cover` that was passed is still called `cover`.
  said3 <- character(0)
  withCallingHandlers(
    build_tessellation(seeds, method = "voronoi", cover = obs),
    message = function(m) { said3 <<- c(said3, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("deriving boundary from `cover`", said3, fixed = TRUE)))

  # The warning for a cover in another place names its source the same way.
  far <- seeds
  attr(far, "cover") <- sf::st_geometry(.fc_square(x0 = 9e6))
  res <- .fc_conditions(build_tessellation(far, method = "voronoi", quiet = TRUE))
  hit <- .fc_has_class(res, "spatialkit_no_cover")
  expect_identical(sum(hit), 1L)
  expect_match(res$warnings[hit], "the cover the seeds carry holds none of the points",
               fixed = TRUE)
  expect_false(grepl("`cover` holds none", res$warnings[hit], fixed = TRUE))
  expect_identical(res$value$params$boundary_source, "points")
  # And with the argument, the argument is named.
  bare <- seeds; attr(bare, "cover") <- NULL
  res2 <- .fc_conditions(build_tessellation(bare, method = "voronoi", quiet = TRUE,
                                            cover = .fc_square(x0 = 9e6)))
  expect_match(res2$warnings[.fc_has_class(res2, "spatialkit_no_cover")],
               "`cover` holds none of the points", fixed = TRUE)
})


# --- one cell per observation is as it was -----------------------------------

test_that("plain one-cell-per-observation Voronoi is unchanged", {
  set.seed(4)
  q <- sf::st_as_sf(data.frame(x = runif(60, 0, 1000), y = runif(60, 0, 1000)),
                    coords = c("x", "y"), crs = 32617)
  expect_no_warning(tess <- build_tessellation(q, method = "voronoi", quiet = TRUE))
  expect_identical(nrow(tess$cells), 60L)
  expect_identical(sort(tess$index), 1:60)
  expect_identical(tess$params$n_outside, 0L)
  expect_identical(tess$params$boundary_source, "points")
  # The boundary is still the points' own hull grown by 2% of its diagonal.
  hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(q)))
  bb <- sf::st_bbox(hull)
  grown <- sf::st_buffer(hull, 0.02 * sqrt((bb[["xmax"]] - bb[["xmin"]])^2 +
                                             (bb[["ymax"]] - bb[["ymin"]])^2))
  expect_equal(as.numeric(sf::st_area(tess$boundary)), as.numeric(sf::st_area(grown)),
               tolerance = 1e-10)
  expect_equal(sum(as.numeric(sf::st_area(tess$cells))), as.numeric(sf::st_area(grown)),
               tolerance = 1e-10)
  # A transect and a single observation raise nothing either.
  line <- sf::st_as_sf(data.frame(x = 5e5 + 0:9 * 10, y = 5e6 + 0:9 * 5),
                       coords = c("x", "y"), crs = 32632)
  expect_no_warning(build_tessellation(line, method = "voronoi", quiet = TRUE))
  expect_no_warning(build_tessellation(q[1, ], method = "voronoi", quiet = TRUE))
})
