# tests/testthat/test-fix-mixed-zm-tessellation.R
# ---------------------------------------------------------------------------
# sf reads a layer in the dimension of its first geometry: st_coordinates()
# of a point layer with an XYZ point after XY ones stops on its "dimnames",
# st_transform() and st_dimension() of one with an XYZ point first stop with
# "OGR error", GEOS stops on the WKB of such a layer and on any geometry
# that carries M.  The tessellation builders read their points, boundary
# and cover in x and y when a layer mixes dimensions or carries M; so do
# the seeders and the cached grid read their boundary.  ensure_stable_poly_id()
# repairs polygons that carry Z or M in XY and says so.  Each layer is held
# to the same rows in XY, to the bit, warnings included.
#
# The layer: 300 points in EPSG:32617 (seed 7) and the same in EPSG:4326,
# 20 k-means seeds, four quadrants over the points, their bounding box, and
# 16 grid squares.
# ---------------------------------------------------------------------------

.zt_catch <- function(expr) {
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

.zt_dims <- function(x)
  vapply(sf::st_geometry(x), function(g) class(g)[1L], character(1), USE.NAMES = FALSE)

.zt_store <- new.env(parent = emptyenv())
.zt_fix <- function() {
  if (is.null(.zt_store$fix)) {
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
    box <- sf::st_sf(id = 1L, geometry = sf::st_sfc(q(bb[[1]], bb[[2]], bb[[3]], bb[[4]]),
                                                   crs = 32617))
    cells <- sf::st_sf(k = 16:1, geometry = sf::st_make_grid(sf::st_as_sfc(bb), n = c(4, 4)))
    seeds <- suppressWarnings(suppressMessages(
      get_voronoi_seeds(method = "kmeans", n = 20, sample_points = pts, set_seed = 1)))
    .zt_store$fix <- list(pts = pts, pts_ll = sf::st_transform(pts, 4326), bnd = bnd,
                          bnd_ll = sf::st_transform(bnd, 4326), box = box,
                          box_ll = sf::st_transform(box, 4326), cells = cells,
                          cells_ll = sf::st_transform(cells, 4326), seeds = seeds,
                          seeds_ll = sf::st_transform(seeds, 4326))
  }
  .zt_store$fix
}

# Row `i` of `l` given a third coordinate (or a third and a fourth for
# XYZM) at its own x and y, keeping the layer's attributes.
.zt_lift <- function(g, dim) {
  extra <- if (dim == "XYZM") c(3, 4) else 3
  co <- unclass(g)
  switch(class(g)[2L],
    POINT   = sf::st_point(c(co, extra), dim = dim),
    POLYGON = sf::st_polygon(lapply(co, function(r)
      cbind(r, matrix(extra, nrow(r), length(extra), byrow = TRUE))), dim = dim))
}
.zt_row <- function(l, i, dim = "XYZ") {
  g <- sf::st_geometry(l)
  g[[i]] <- .zt_lift(g[[i]], dim)
  sf::st_geometry(l) <- sf::st_sfc(as.list(g), crs = sf::st_crs(l), recompute_bbox = TRUE)
  l
}
.zt_all <- function(l, dim = "XYZ") {
  g <- lapply(sf::st_geometry(l), .zt_lift, dim = dim)
  sf::st_geometry(l) <- sf::st_sfc(g, crs = sf::st_crs(l))
  l
}

# `f` on the layer with row `i` given Z or M is identical() to `f` on the
# layer, warnings included.
.zt_same <- function(f, layer, rows, dims = c("XYZ", "XYM"), what = "") {
  ref <- .zt_catch(f(layer))
  for (d in dims) for (i in rows) {
    mixed <- .zt_row(layer, i, d)
    expect_identical(unique(.zt_dims(mixed)), if (i == 1L) c(d, "XY") else c("XY", d))
    got <- .zt_catch(f(mixed))
    info <- sprintf("%s row %d %s", what, i, d)
    expect_identical(got$value, ref$value, info = info)
    expect_identical(got$classes, ref$classes, info = info)
    expect_identical(got$texts, ref$texts, info = info)
  }
  invisible(ref)
}

test_that("the premise: the builders' readers stop on a mixed point layer or on M", {
  F <- .zt_fix()
  expect_error(suppressWarnings(sf::st_coordinates(.zt_row(F$seeds, 5))), "dimnames")
  expect_error(suppressWarnings(sf::st_transform(.zt_row(F$pts, 1), 4326)), "OGR error")
  expect_error(sf::st_union(.zt_all(F$box, "XYM")), "XYM")
  expect_identical(sf::st_zm(.zt_row(F$seeds, 5)), F$seeds)
  expect_identical(attr(sf::st_zm(.zt_row(F$seeds, 1)), "cover"), attr(F$seeds, "cover"))
})

test_that("create_voronoi_polygons() and build_tessellation() read mixed seeds as in XY", {
  F <- .zt_fix()
  for (S in list(F$seeds, F$seeds_ll)) {
    tag <- format(sf::st_crs(S)$epsg)
    ref <- .zt_same(function(s) create_voronoi_polygons(s, quiet = TRUE), S,
                    rows = c(1L, 5L), what = paste("cvp", tag))
    expect_identical(ref$value$params$boundary_source, "seeds' cover")
    .zt_same(function(s) build_tessellation(s, method = "voronoi", quiet = TRUE), S,
             rows = c(1L, 5L), what = paste("voronoi", tag))
  }
})

test_that("build_tessellation() grids, triangles and clip_target_for() read mixed points as in XY", {
  F <- .zt_fix()
  for (k in c("pts", "pts_ll")) {
    P <- F[[k]]
    B <- if (k == "pts") F$box else F$box_ll
    for (m in c("hex", "square", "triangles"))
      .zt_same(function(p) build_tessellation(p, boundary = B, method = m,
                                              approx_n_cells = if (m != "triangles") 20,
                                              quiet = TRUE),
               P, rows = c(1L, 5L), what = paste(m, k))
    # In degrees the points are projected: flattened first, so
    # ensure_projected() raises no warning about an XY copy that is not
    # returned.
    ref <- .zt_same(function(p) clip_target_for(p, quiet = TRUE), P, rows = c(1L, 5L),
                    what = paste("clip", k))
    expect_false("spatialkit_zm_dropped" %in% ref$classes)
    .zt_same(function(p) clip_target_for(p, boundary = B, expand = 0.05, quiet = TRUE),
             P, rows = c(1L, 5L), what = paste("clip boundary", k))
  }
})

test_that("a mixed cover gives the XY cover's cells and no false warning", {
  F <- .zt_fix()
  for (CV in list(F$bnd, F$bnd_ll)) {
    tag <- format(sf::st_crs(CV)$epsg)
    ref <- .zt_same(function(cv) create_voronoi_polygons(F$seeds, cover = cv, quiet = TRUE),
                    CV, rows = c(1L, 3L), what = paste("cvp cover", tag))
    expect_identical(ref$value$params$boundary_source, "cover")
    expect_false("spatialkit_no_cover" %in% ref$classes)
    ref <- .zt_same(function(cv) build_tessellation(F$seeds, method = "voronoi", cover = cv,
                                                    quiet = TRUE),
                    CV, rows = c(1L, 3L), what = paste("voronoi cover", tag))
    expect_identical(ref$value$params$boundary_source, "cover")
    # The cover the seeds carry.
    with_cover <- function(cv) { s <- F$seeds; attr(s, "cover") <- sf::st_geometry(cv); s }
    ref <- .zt_same(function(cv) create_voronoi_polygons(with_cover(cv), quiet = TRUE),
                    CV, rows = c(1L, 3L), what = paste("cvp seeds' cover", tag))
    expect_identical(ref$value$params$boundary_source, "seeds' cover")
    .zt_same(function(cv) build_tessellation(with_cover(cv), method = "voronoi", quiet = TRUE),
             CV, rows = c(1L, 3L), what = paste("voronoi seeds' cover", tag))
  }
  # A point cover checked against cells built in a boundary.
  ref <- .zt_same(function(cv) create_voronoi_polygons(F$seeds, boundary = F$box, cover = cv,
                                                       quiet = TRUE),
                  F$pts, rows = c(1L, 5L), what = "cvp point cover")
  expect_identical(ref$value$params$cover_outside, 0L)
  ref <- .zt_same(function(cv) build_tessellation(F$pts, boundary = F$box, method = "hex",
                                                  approx_n_cells = 20, cover = cv,
                                                  quiet = TRUE),
                  F$pts, rows = c(1L, 5L), what = "hex point cover")
  expect_identical(ref$value$params$cover_outside, 0L)
})

test_that("points and a cover that carry M throughout are read as in XY", {
  F <- .zt_fix()
  for (d in c("XYM", "XYZM")) {
    ref <- .zt_catch(create_voronoi_polygons(F$seeds, quiet = TRUE))
    got <- .zt_catch(create_voronoi_polygons(.zt_all(F$seeds, d), quiet = TRUE))
    expect_identical(got[1:3], ref[1:3], info = d)
    ref <- .zt_catch(build_tessellation(F$pts, boundary = F$box, method = "hex",
                                        approx_n_cells = 16, quiet = TRUE))
    got <- .zt_catch(build_tessellation(.zt_all(F$pts, d), boundary = F$box, method = "hex",
                                        approx_n_cells = 16, quiet = TRUE))
    expect_identical(got[1:3], ref[1:3], info = d)
    ref <- .zt_catch(create_voronoi_polygons(F$seeds, cover = F$bnd, quiet = TRUE))
    got <- .zt_catch(create_voronoi_polygons(F$seeds, cover = .zt_all(F$bnd, d), quiet = TRUE))
    expect_identical(got[1:3], ref[1:3], info = d)
    expect_identical(got$value$params$boundary_source, "cover")
  }
})

test_that("a boundary that carries M throughout is read as the same boundary in XY", {
  F <- .zt_fix()
  calls <- list(
    grid    = function(b) create_grid_polygons(b, target_cells = 16),
    cached  = function(b) create_grid_polygons_cached(b, target_cells = 16, cache_env = new.env()),
    random  = function(b) get_voronoi_seeds(b, method = "random", n = 10, set_seed = 1),
    kmeans  = function(b) get_voronoi_seeds(b, method = "kmeans", n = 10, set_seed = 1),
    vsr     = function(b) voronoi_seeds_random(b, k = 10, set_seed = 1),
    voronoi = function(b) create_voronoi_polygons(F$seeds[1:15, ], boundary = b, quiet = TRUE),
    tess_v  = function(b) build_tessellation(F$seeds[1:15, ], method = "voronoi", boundary = b,
                                             quiet = TRUE),
    tess_h  = function(b) build_tessellation(F$pts, boundary = b, method = "hex",
                                             approx_n_cells = 16, quiet = TRUE),
    clip    = function(b) clip_target_for(F$pts, boundary = b, quiet = TRUE))
  for (k in c("box", "box_ll")) {
    B <- F[[k]]
    for (nm in names(calls)) {
      ref <- .zt_catch(calls[[nm]](B))
      expect_false("spatialkit_zm_dropped" %in% ref$classes)
      for (d in c("XYM", "XYZM")) {
        b <- .zt_all(B, d)
        expect_identical(.zt_dims(b), d)
        got <- .zt_catch(calls[[nm]](b))
        info <- paste(nm, k, d)
        expect_identical(got$value, ref$value, info = info)
        expect_identical(got$classes, ref$classes, info = info)
        expect_identical(got$texts, ref$texts, info = info)
      }
    }
    # One cache entry for the boundary in XY and with M.
    env <- new.env()
    ref <- suppressWarnings(suppressMessages(
      create_grid_polygons_cached(B, target_cells = 16, cache_env = env)))
    for (d in c("XYM", "XYZM")) {
      expect_identical(suppressWarnings(suppressMessages(
        create_grid_polygons_cached(.zt_all(B, d), target_cells = 16, cache_env = env))), ref)
      expect_length(ls(env), 1L)
    }
  }
  # The seeders' cover is the boundary read in XY.
  s <- voronoi_seeds_random(.zt_all(F$box, "XYM"), k = 10, set_seed = 1)
  expect_identical(attr(s, "cover"), sf::st_geometry(F$box))
})

test_that("a boundary XYZ throughout takes the path it took", {
  F <- .zt_fix()
  bz <- .zt_all(F$box)
  # The cells and the clip target are the XY boundary's, with no warning,
  # for the boundary projected and in degrees.
  calls <- list(function(b) create_grid_polygons(b, target_cells = 16),
                function(b) create_grid_polygons_cached(b, target_cells = 16,
                                                        cache_env = new.env()),
                function(b) create_voronoi_polygons(F$seeds[1:15, ], boundary = b, quiet = TRUE),
                function(b) build_tessellation(F$seeds[1:15, ], method = "voronoi",
                                               boundary = b, quiet = TRUE),
                function(b) build_tessellation(F$pts, boundary = b, method = "hex",
                                               approx_n_cells = 16, quiet = TRUE),
                function(b) clip_target_for(F$pts, boundary = b, quiet = TRUE))
  for (k in c("box", "box_ll")) for (i in seq_along(calls)) {
    info <- paste(k, i)
    ref <- .zt_catch(calls[[i]](F[[k]]))
    got <- .zt_catch(calls[[i]](.zt_all(F[[k]])))
    expect_identical(got$value, ref$value, info = info)
    expect_identical(got$classes, ref$classes, info = info)
    expect_false("spatialkit_zm_dropped" %in% got$classes, info = info)
  }
  # The seeders' cover keeps its Z: a boundary XYZ throughout is not flattened.
  s <- voronoi_seeds_random(bz, k = 10, set_seed = 1)
  expect_identical(attr(s, "cover"), sf::st_geometry(bz))
  # Points and a cover XYZ throughout give the XY cells.
  expect_identical(create_voronoi_polygons(.zt_all(F$seeds), quiet = TRUE),
                   create_voronoi_polygons(F$seeds, quiet = TRUE))
  expect_identical(create_voronoi_polygons(F$seeds, cover = .zt_all(F$bnd), quiet = TRUE),
                   create_voronoi_polygons(F$seeds, cover = F$bnd, quiet = TRUE))
  # The helper hands a layer XY or XYZ throughout back as it came.
  expect_identical(.xy_if_mixed_or_m(bz), bz)
  expect_identical(.xy_if_mixed_or_m(F$box), F$box)
  expect_identical(.xy_if_mixed_or_m(.zt_all(F$box, "XYM")), F$box)
})

test_that("ensure_stable_poly_id() repairs polygons that carry Z or M throughout in XY, with a warning", {
  F <- .zt_fix()
  cz <- .zt_all(F$cells)
  ce <- cz
  g <- sf::st_geometry(ce)
  g[[5]] <- sf::st_polygon()
  sf::st_geometry(ce) <- sf::st_sfc(as.list(g), crs = sf::st_crs(ce))
  expect_identical(unique(.zt_dims(ce)), c("XYZ", "XY"))
  layers <- list(xyz = cz, xyz_ll = .zt_all(F$cells_ll), xym = .zt_all(F$cells, "XYM"),
                 xyzm = .zt_all(F$cells, "XYZM"), xyz_empty = ce)
  for (nm in names(layers)) {
    L <- layers[[nm]]
    rows <- which(.zt_dims(L) != "XY")
    ref <- .zt_catch(ensure_stable_poly_id(sf::st_zm(L)))
    expect_length(ref$classes, 0L)
    got <- .zt_catch(ensure_stable_poly_id(L))
    expect_identical(got$value, ref$value, info = nm)
    expect_identical(unique(.zt_dims(got$value)), "XY", info = nm)
    expect_identical(got$classes, "spatialkit_zm_dropped", info = nm)
    expect_identical(got$conds[[1L]]$rows, rows, info = nm)
    expect_identical(got$texts, sprintf(paste0(
      "ensure_stable_poly_id(): %d of the %d geometries of `polygons_sf` carry Z or M ",
      "coordinates, and the repair with sf::st_make_valid() is made in x and y, so the ",
      "repaired polygons are returned in XY. Pass make_valid = FALSE to have the rows ",
      "returned with their Z and M, or drop them with sf::st_zm() first to silence this."),
      length(rows), nrow(L)), info = nm)
    # Without the repair the caller's rows come back with their Z or M.
    raw <- .zt_catch(ensure_stable_poly_id(L, make_valid = FALSE))
    expect_false("spatialkit_zm_dropped" %in% raw$classes)
    expect_identical(raw$value$k, ref$value$k, info = nm)
    expect_identical(sf::st_geometry(raw$value), sf::st_geometry(L)[match(raw$value$k, L$k)],
                     info = nm)
  }
  # A mixed layer keeps the mixed layer's warning, once.
  got <- .zt_catch(ensure_stable_poly_id(.zt_row(F$cells, 5)))
  expect_identical(got$classes, "spatialkit_zm_dropped")
  expect_match(got$texts, "ensure_stable_poly_id(): `polygons_sf` mixes", fixed = TRUE)
  expect_identical(got$value, ensure_stable_poly_id(F$cells))
  # A layer XY throughout raises no such warning.
  expect_false("spatialkit_zm_dropped" %in% .zt_catch(ensure_stable_poly_id(F$cells))$classes)
})

test_that("ensure_stable_poly_id() raises no warning for empty polygons in Z or M", {
  # An empty polygon after the first has no coordinates, so nothing is
  # dropped from it: polygons XY with an empty POLYGON Z and an empty
  # POLYGON M among them give the value of their copy in XY, with no
  # warning.  An XYZ layer with such an empty polygon counts only the rows
  # that hold Z.
  F <- .zt_fix()
  put_empty <- function(l, i, dim) {
    g <- sf::st_geometry(l)
    for (j in i) g[[j]] <- sf::st_polygon(dim = dim)
    sf::st_geometry(l) <- sf::st_sfc(as.list(g), crs = sf::st_crs(l))
    l
  }
  for (cells in list(F$cells, F$cells_ll)) {
    e <- put_empty(put_empty(cells, 3, "XYZ"), 7, "XYM")
    expect_identical(unique(.zt_dims(e)), c("XY", "XYZ", "XYM"))
    ref <- .zt_catch(ensure_stable_poly_id(sf::st_zm(e)))
    got <- .zt_catch(ensure_stable_poly_id(e))
    expect_identical(got$value, ref$value)
    expect_identical(got$classes, ref$classes)
    expect_false("spatialkit_zm_dropped" %in% got$classes)
    expect_identical(unique(.zt_dims(got$value)), "XY")
  }
  ez <- put_empty(.zt_all(F$cells), 9, "XYM")
  got <- .zt_catch(ensure_stable_poly_id(ez))
  expect_identical(got$classes, "spatialkit_zm_dropped")
  expect_identical(got$conds[[1L]]$rows, setdiff(1:16, 9L))
  expect_match(got$texts, "15 of the 16 geometries", fixed = TRUE)
  expect_identical(got$value, ensure_stable_poly_id(sf::st_zm(ez)))
})

test_that("the help says the builders read their layers in x and y", {
  # The source page from the source tree, the installed page under R CMD
  # check, which runs the tests where there is no man/.
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  for (nm in c("clip_target_for.Rd", "create_voronoi_polygons.Rd", "build_tessellation.Rd"))
    expect_match(page(nm), "Points that mix XY and XYZ (or XYM), or that carry M, are read in x and y",
                 fixed = TRUE, info = nm)
  for (nm in c("create_voronoi_polygons.Rd", "build_tessellation.Rd"))
    expect_match(page(nm), paste0("A cover, given or attached, that mixes XY and XYZ (or XYM) ",
                                  "geometries, or that carries M, is read in x and y"),
                 fixed = TRUE, info = nm)
  # A boundary whose geometries carry M is read as ?spatialkit-geometry
  # says, in x and y.
  for (nm in c("clip_target_for.Rd", "create_voronoi_polygons.Rd", "create_grid_polygons.Rd",
               "build_tessellation.Rd", "get_voronoi_seeds.Rd", "voronoi_seeds_random.Rd"))
    expect_match(page(nm), "mixed dimensions are handled as \"Boundaries\" in spatialkit-geometry",
                 fixed = TRUE, info = nm)
  expect_match(page("spatialkit-geometry.Rd"),
               "A boundary that mixes dimensions or carries M is read in x and y.", fixed = TRUE)
  expect_match(page("spatialkit-geometry.Rd"), "GEOS refuses M", fixed = TRUE)
  expect_match(page("create_grid_polygons_cached.Rd"),
               "or whose geometries carry M (XYM or XYZM), is read without Z and M", fixed = TRUE)
  expect_match(page("ensure_stable_poly_id.Rd"),
               "and polygons that all carry Z or M, are repaired in XY", fixed = TRUE)
  expect_match(page("ensure_stable_poly_id.Rd"),
               "An empty polygon in Z or M has no coordinates to drop and raises no warning",
               fixed = TRUE)
})
