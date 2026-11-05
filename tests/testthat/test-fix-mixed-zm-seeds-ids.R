# tests/testthat/test-fix-mixed-zm-seeds-ids.R
# ---------------------------------------------------------------------------
# sf reads a layer in the dimension of its first geometry: st_coordinates()
# of a point layer with an XYZ point first among XY points lays its values
# out of step, and with the XYZ point later stops on its "dimnames";
# st_transform() of a layer with an XYZ geometry first stops with "OGR
# error"; GEOS reads an XYZ ring among XY rings as XY pairs, and s2 drops Z
# only when the first geometry has it.  The seeders read their points and
# their boundary, ensure_stable_poly_id() its polygons, and
# create_grid_polygons_cached() its boundary in x and y: Z and M are dropped
# from the copy read when a layer mixes dimensions.  Each mixed layer is held
# to the same rows in XY, to the bit, warnings included.
#
# The layer: 300 points in EPSG:32617 (seed 7) and the same in EPSG:4326, a
# boundary of four quadrants over them, and 16 grid squares.
# ---------------------------------------------------------------------------

.zi_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)),
       conds = got)
}

.zi_dims <- function(x)
  vapply(sf::st_geometry(x), function(g) class(g)[1L], character(1), USE.NAMES = FALSE)

.zi_store <- new.env(parent = emptyenv())
.zi_fix <- function() {
  if (is.null(.zi_store$fix)) {
    d <- local({
      withr::local_seed(7, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 300
      data.frame(x = 5e5 + stats::runif(n, 0, 1e5), y = 5e6 + stats::runif(n, 0, 1e5),
                 v = stats::rnorm(n))
    })
    pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
    bb <- sf::st_bbox(pts)
    mx <- (bb[[1]] + bb[[3]]) / 2; my <- (bb[[2]] + bb[[4]]) / 2
    q <- function(x0, y0, x1, y1)
      sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
    bnd <- sf::st_sf(id = 1:4, geometry = sf::st_sfc(
      q(bb[[1]], bb[[2]], mx, my), q(mx, bb[[2]], bb[[3]], my),
      q(bb[[1]], my, mx, bb[[4]]), q(mx, my, bb[[3]], bb[[4]]), crs = 32617))
    cells <- sf::st_sf(k = 16:1, geometry = sf::st_make_grid(sf::st_as_sfc(bb), n = c(4, 4)))
    prof <- suppressWarnings(suppressMessages(
      resolution_profile(pts, response_var = "v", levels = c(8, 12))))
    .zi_store$fix <- list(pts = pts, pts_ll = sf::st_transform(pts, 4326), bnd = bnd,
                          bnd_ll = sf::st_transform(bnd, 4326), cells = cells,
                          cells_ll = sf::st_transform(cells, 4326), prof = prof)
  }
  .zi_store$fix
}

# Row `i` of `l` given a third coordinate, Z or M, at its own x and y.
.zi_row <- function(l, i, dim = "XYZ") {
  g  <- sf::st_geometry(l)
  co <- unclass(g[[i]])
  g[[i]] <- switch(class(g[[i]])[2L],
    POINT        = sf::st_point(c(co, 3), dim = dim),
    POLYGON      = sf::st_polygon(lapply(co, function(r) cbind(r, 3)), dim = dim),
    MULTIPOLYGON = sf::st_multipolygon(lapply(co, function(p)
      lapply(p, function(r) cbind(r, 3))), dim = dim))
  sf::st_set_geometry(l, sf::st_sfc(as.list(g), crs = sf::st_crs(l), recompute_bbox = TRUE))
}
.zi_all <- function(l) {
  for (i in seq_len(nrow(l))) l <- .zi_row(l, i)
  l
}

# `f` on the layer with row `i` given Z or M is identical() to `f` on the
# layer, warnings included.
.zi_same <- function(f, layer, rows, dims = c("XYZ", "XYM")) {
  ref <- .zi_catch(f(layer))
  for (d in dims) for (i in rows) {
    mixed <- .zi_row(layer, i, d)
    expect_identical(unique(.zi_dims(mixed)), if (i == 1L) c(d, "XY") else c("XY", d))
    got <- .zi_catch(f(mixed))
    info <- sprintf("row %d %s", i, d)
    expect_identical(got$value, ref$value, info = info)
    expect_identical(got$classes, ref$classes, info = info)
    expect_identical(got$texts, ref$texts, info = info)
  }
  invisible(ref)
}

test_that("the premise: sf reads a mixed point layer in the dimension of its first point", {
  F <- .zi_fix()
  expect_warning(sf::st_coordinates(.zi_row(F$pts, 1)), "not a sub-multiple")
  expect_error(suppressWarnings(sf::st_coordinates(.zi_row(F$pts, 5))), "dimnames")
  expect_identical(sf::st_geometry(sf::st_zm(.zi_row(F$pts, 5))), sf::st_geometry(F$pts))
})

test_that("get_voronoi_seeds(method = 'kmeans') reads mixed sample_points as the same points in XY", {
  F <- .zi_fix()
  for (P in list(F$pts, F$pts_ll))
    .zi_same(function(p) get_voronoi_seeds(method = "kmeans", n = 12, sample_points = p,
                                           set_seed = 1), P, rows = c(1L, 5L))
  .zi_same(function(p) get_voronoi_seeds(F$bnd_ll, method = "kmeans", n = 12,
                                         sample_points = p, set_seed = 1),
           F$pts, rows = c(1L, 5L), dims = "XYZ")
})

test_that("get_voronoi_seeds() with a profile as n reads mixed sample_points as in XY", {
  F <- .zi_fix()
  for (P in list(F$pts, F$pts_ll)) {
    ref <- .zi_same(function(p) get_voronoi_seeds(method = "kmeans", n = F$prof,
                                                  sample_points = p, set_seed = 1),
                    P, rows = c(1L, 5L))
    expect_length(attr(ref$value, "kmeans")$cluster, nrow(P))
  }
  .zi_same(function(p) get_voronoi_seeds(F$bnd, method = "kmeans", n = F$prof,
                                         sample_points = p, set_seed = 1),
           F$pts_ll, rows = c(1L, 5L), dims = "XYZ")
})

test_that("voronoi_seeds_kmeans() reads mixed points as the same points in XY", {
  F <- .zi_fix()
  for (P in list(F$pts, F$pts_ll))
    .zi_same(function(p) voronoi_seeds_kmeans(p, k = 12, set_seed = 1), P, rows = c(1L, 5L))
  .zi_same(function(p) voronoi_seeds_kmeans(p, k = F$prof), F$pts, rows = c(1L, 5L))
})

test_that("the seeders read a mixed boundary as the same boundary in XY", {
  F <- .zi_fix()
  for (B in list(F$bnd, F$bnd_ll)) {
    ref <- .zi_same(function(b) voronoi_seeds_random(b, k = 10, set_seed = 1), B,
                    rows = c(1L, 3L))
    expect_identical(attr(ref$value, "cover"), sf::st_geometry(B))
    ref <- .zi_same(function(b) get_voronoi_seeds(b, method = "random", n = 10, set_seed = 1),
                    B, rows = c(1L, 3L))
    expect_identical(attr(ref$value, "cover"), sf::st_geometry(B))
    .zi_same(function(b) get_voronoi_seeds(b, method = "kmeans", n = 10, set_seed = 1),
             B, rows = c(1L, 3L))
    .zi_same(function(b) get_voronoi_seeds(b, method = "kmeans", n = 10, sample_points = F$pts,
                                           set_seed = 1), B, rows = c(1L, 3L), dims = "XYZ")
  }
  # In degrees an XYZ quadrant after the first gave random seeds outside the
  # boundary; they now fall inside it.
  s <- voronoi_seeds_random(.zi_row(F$bnd_ll, 3), k = 10, set_seed = 1)
  expect_true(all(lengths(sf::st_intersects(s, F$bnd_ll)) > 0L))
})

test_that("ensure_stable_poly_id() keys mixed polygons as the same polygons in XY", {
  F <- .zi_fix()
  for (C in list(F$cells, F$cells_ll)) {
    ref <- .zi_catch(ensure_stable_poly_id(C))
    ref_raw <- .zi_catch(ensure_stable_poly_id(C, make_valid = FALSE))
    for (d in c("XYZ", "XYM")) for (i in c(1L, 5L)) {
      info <- sprintf("row %d %s", i, d)
      mixed <- .zi_row(C, i, d)
      # The repaired copy is returned, in XY, with the warning.
      got <- .zi_catch(ensure_stable_poly_id(mixed))
      expect_identical(got$value, ref$value, info = info)
      expect_identical(got$classes, c(ref$classes, "spatialkit_zm_dropped"), info = info)
      w <- got$conds[[length(got$conds)]]
      expect_identical(w$rows, i, info = info)
      expect_match(conditionMessage(w), "ensure_stable_poly_id(): `polygons_sf` mixes",
                   fixed = TRUE, info = info)
      expect_false(any(grepl("s2 refused", got$texts)), info = info)
      # Without the repair the caller's rows come back with their Z.
      got <- .zi_catch(ensure_stable_poly_id(mixed, make_valid = FALSE))
      expect_identical(got$value$poly_id, ref_raw$value$poly_id, info = info)
      expect_identical(got$value$k, ref_raw$value$k, info = info)
      expect_identical(got$texts, ref_raw$texts, info = info)
      expect_identical(sf::st_geometry(got$value), sf::st_geometry(mixed)[match(got$value$k, mixed$k)],
                       info = info)
    }
  }
})

test_that("ensure_stable_poly_id() keys mixed Voronoi cells in degrees as in XY", {
  F <- .zi_fix()
  seeds <- get_voronoi_seeds(method = "kmeans", n = 12, sample_points = F$pts, set_seed = 1)
  vc <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)$cells
  vc <- sf::st_transform(vc[, setdiff(names(vc), "poly_id")], 4326)
  vc$k <- seq_len(nrow(vc))
  ref <- ensure_stable_poly_id(vc, make_valid = FALSE)
  for (i in c(1L, 5L)) {
    got <- .zi_catch(ensure_stable_poly_id(.zi_row(vc, i), make_valid = FALSE))
    expect_identical(got$value$k, ref$k)
    expect_length(got$classes, 0L)
    got <- .zi_catch(ensure_stable_poly_id(.zi_row(vc, i)))
    expect_identical(got$value$k, ensure_stable_poly_id(vc)$k)
    expect_identical(got$classes, "spatialkit_zm_dropped")
  }
})

test_that("create_grid_polygons_cached() reads a mixed boundary as the same boundary in XY", {
  F <- .zi_fix()
  # The boundary in degrees is projected for the grid; the flattening comes
  # first, so the projection raises no warning about a copy the caller never
  # sees.
  for (B in list(projected = F$bnd, degrees = F$bnd_ll)) {
    env <- new.env()
    ref <- .zi_catch(create_grid_polygons_cached(B, target_cells = 16, cache_env = env))
    expect_length(ls(env), 1L)
    expect_false("spatialkit_zm_dropped" %in% ref$classes)
    for (d in c("XYZ", "XYM")) for (i in c(1L, 3L)) {
      info <- sprintf("%s row %d %s", sf::st_crs(B)$epsg, i, d)
      mixed <- .zi_row(B, i, d)
      got <- .zi_catch(create_grid_polygons_cached(mixed, target_cells = 16,
                                                   cache_env = new.env()))
      expect_identical(got$value, ref$value, info = info)
      expect_identical(got$classes, ref$classes, info = info)
      expect_identical(got$texts, ref$texts, info = info)
      # One cache entry for a mixed boundary and its XY copy.
      expect_identical(suppressWarnings(suppressMessages(
        create_grid_polygons_cached(mixed, target_cells = 16, cache_env = env))),
        ref$value, info = info)
      expect_length(ls(env), 1L)
    }
  }
})

test_that("layers in one dimension take the path they took", {
  F <- .zi_fix()
  # XYZ throughout: the seeders read x and y, as they did.
  ref <- get_voronoi_seeds(method = "kmeans", n = 12, sample_points = F$pts, set_seed = 1)
  expect_identical(get_voronoi_seeds(method = "kmeans", n = 12, sample_points = .zi_all(F$pts),
                                     set_seed = 1), ref)
  expect_identical(voronoi_seeds_kmeans(.zi_all(F$pts), k = 12, set_seed = 1),
                   voronoi_seeds_kmeans(F$pts, k = 12, set_seed = 1))
  bz <- .zi_all(F$bnd)
  expect_identical(.zi_dims(bz), rep("XYZ", 4L))
  s <- voronoi_seeds_random(bz, k = 10, set_seed = 1)
  expect_identical(sf::st_geometry(s), sf::st_geometry(voronoi_seeds_random(F$bnd, k = 10, set_seed = 1)))
  # The cover of a boundary XYZ throughout keeps its Z: nothing flattens it.
  expect_identical(attr(s, "cover"), sf::st_geometry(bz))
  # ensure_stable_poly_id() on cells XYZ throughout: no warning about Z, the
  # XY cells' IDs, and without the repair the caller's own rows.
  cz <- .zi_all(F$cells)
  got <- .zi_catch(ensure_stable_poly_id(cz))
  expect_false("spatialkit_zm_dropped" %in% got$classes)
  expect_identical(got$value$k, ensure_stable_poly_id(F$cells)$k)
  got <- suppressMessages(ensure_stable_poly_id(cz, make_valid = FALSE))
  expect_identical(.zi_dims(got), rep("XYZ", 16L))
  expect_identical(sf::st_geometry(got), sf::st_geometry(cz)[match(got$k, cz$k)])
  # The cache key of a boundary XYZ throughout is the one it had: taken on
  # the boundary with its Z, so it is not the XY boundary's.
  env <- new.env()
  create_grid_polygons_cached(bz, target_cells = 16, cache_env = env)
  expect_identical(ls(env), .cache_key(.project_for_grid(bz, "create_grid_polygons_cached"),
                                       "square", 16))
  create_grid_polygons_cached(F$bnd, target_cells = 16, cache_env = env)
  expect_length(ls(env), 2L)
})

test_that("get_voronoi_seeds(method = 'provided') moves mixed seeds into the boundary's CRS in XY", {
  F <- .zi_fix()
  seeds <- F$pts[1:12, ]
  ref <- .zi_catch(get_voronoi_seeds(F$bnd_ll, method = "provided", seeds = seeds))
  expect_length(ref$classes, 0L)
  for (d in c("XYZ", "XYM")) for (i in c(1L, 5L)) {
    info <- sprintf("row %d %s", i, d)
    mixed <- .zi_row(seeds, i, d)
    got <- .zi_catch(get_voronoi_seeds(F$bnd_ll, method = "provided", seeds = mixed))
    expect_identical(got$value, ref$value, info = info)
    expect_identical(got$classes, "spatialkit_zm_dropped", info = info)
    expect_identical(got$conds[[1L]]$rows, i, info = info)
    expect_match(got$texts, "get_voronoi_seeds(): `seeds` mixes", fixed = TRUE)
    # Seeds that are not moved come back as they were given, Z included.
    got <- .zi_catch(get_voronoi_seeds(F$bnd, method = "provided", seeds = mixed))
    expect_length(got$classes, 0L)
    expect_identical(sf::st_geometry(got$value), sf::st_geometry(mixed), info = info)
  }
})

test_that("the help says the seeders, the IDs and the cached grid read mixed layers in x and y", {
  # The source page from the source tree, the installed page under R CMD
  # check, which runs the tests where there is no man/.
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  txt <- page("get_voronoi_seeds.Rd")
  expect_match(txt, paste0("A cloud that mixes XY and XYZ (or XYM) points is read without Z and ",
                           "M before it is projected or clustered"), fixed = TRUE)
  expect_match(txt, paste0("seeds that mix XY and XYZ (or XYM) points are then returned in XY, ",
                           "with a warning of class \"spatialkit_zm_dropped\""), fixed = TRUE)
  expect_match(txt, "defining the sampling area. One that mixes XY and XYZ (or XYM) geometries is read without Z and M",
               fixed = TRUE)
  expect_match(page("voronoi_seeds_kmeans.Rd"),
               "Points that mix XY and XYZ (or XYM) are read without Z and M", fixed = TRUE)
  expect_match(page("voronoi_seeds_random.Rd"),
               "One that mixes XY and XYZ (or XYM) geometries is read without Z and M", fixed = TRUE)
  txt <- page("ensure_stable_poly_id.Rd")
  expect_match(txt, "are keyed without Z and M", fixed = TRUE)
  expect_match(txt, paste0("are repaired in XY and returned in XY, with a warning of class ",
                           "\"spatialkit_zm_dropped\". With FALSE the rows come back with the ",
                           "geometry they arrived with, Z and M included."), fixed = TRUE)
  expect_match(page("create_grid_polygons_cached.Rd"),
               "shares its cache entry with the same boundary in XY", fixed = TRUE)
})
