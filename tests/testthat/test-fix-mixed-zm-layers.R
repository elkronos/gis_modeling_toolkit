# tests/testthat/test-fix-mixed-zm-layers.R
# ---------------------------------------------------------------------------
# sf reads a layer in the dimension of its first geometry: s2 drops Z and M
# only when the first geometry carries them, and GEOS and st_transform() read
# an XYZ ring among XY rows as XY pairs.  rbind() of an XY and an XYZ layer
# gives such a layer without a word.  The join, the cells' areas, the CRS
# chooser and the tessellation builders' `boundary` now hand sf a copy of a
# mixed layer without Z and M; the rows and cells returned keep their own
# geometry; ensure_projected() and harmonize_crs() return a transformed copy
# of a mixed layer in XY with a warning of class "spatialkit_zm_dropped".
#
# The layer: six 0.1-degree squares near (-94.3, 39.6) in EPSG:4326, one
# point in each.  sf::st_zm() of the mixed squares is identical() to the XY
# squares, so every result on a mixed layer is held to the XY layer's.
# ---------------------------------------------------------------------------

.mz_catch <- function(expr) {
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

.mz_dims <- function(x)
  vapply(sf::st_geometry(x), function(g) class(g)[1L], character(1), USE.NAMES = FALSE)

.mz_sq <- function(x0, y0, d = 0.1)
  rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d), c(x0, y0))

.mz_fix <- function() {
  base <- lapply(0:5, function(i) sf::st_polygon(list(.mz_sq(-94.3 + 0.1 * i, 39.6))))
  squares <- sf::st_sf(cell_id = 1:6, geometry = sf::st_sfc(base, crs = 4326))
  xy <- cbind(-94.3 + 0.1 * (0:5) + 0.04, 39.645)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], row = 1:6, z = 1:6),
                      coords = c("x", "y"), crs = 4326)
  # Row `i` of the squares given a third coordinate, Z or M, at 7.
  with3 <- function(i, dim = "XYZ") {
    g <- base
    g[[i]] <- sf::st_polygon(list(cbind(.mz_sq(-94.3 + 0.1 * (i - 1), 39.6), 7)), dim = dim)
    sf::st_sf(cell_id = 1:6, geometry = sf::st_sfc(g, crs = 4326))
  }
  mixed <- with3(3)
  allz <- sf::st_sf(cell_id = 1:6, geometry = sf::st_sfc(lapply(0:5, function(i)
    sf::st_polygon(list(cbind(.mz_sq(-94.3 + 0.1 * i, 39.6), 7)))), crs = 4326))
  # Patches around each point; patch 2 XYZ.
  pg <- lapply(1:6, function(i) sf::st_polygon(list(.mz_sq(xy[i, 1] - 0.01, xy[i, 2] - 0.01, 0.02))))
  patches <- sf::st_sf(row = 1:6, geometry = sf::st_sfc(pg, crs = 4326))
  pg[[2]] <- sf::st_polygon(list(cbind(.mz_sq(xy[2, 1] - 0.01, xy[2, 2] - 0.01, 0.02), 3)))
  patches_z <- sf::st_sf(row = 1:6, geometry = sf::st_sfc(pg, crs = 4326))
  # Points; point `i` XYZ.
  pts_z <- function(i) {
    g <- lapply(1:6, function(k) if (k == i) sf::st_point(c(xy[k, ], 3)) else sf::st_point(xy[k, ]))
    sf::st_sf(row = 1:6, z = 1:6, geometry = sf::st_sfc(g, crs = 4326))
  }
  list(squares = squares, pts = pts, with3 = with3, mixed = mixed, allz = allz,
       patches = patches, patches_z = patches_z, pts_z = pts_z, xy = xy)
}

.mz_join <- function(f, p, ...)
  .mz_catch(assign_features_to_polygons(f, p, "cell_id", quiet = TRUE, ...))

test_that("the premise: sf reads a mixed layer in the dimension of its first geometry", {
  F <- .mz_fix()
  expect_identical(.mz_dims(F$mixed), c("XY", "XY", "XYZ", "XY", "XY", "XY"))
  expect_identical(lengths(suppressMessages(sf::st_intersects(F$pts, F$mixed))),
                   c(1L, 1L, 0L, 2L, 2L, 2L))
  expect_identical(sf::st_geometry(sf::st_zm(F$mixed)), sf::st_geometry(F$squares))
})

test_that("points are joined to cells that mix XY and XYZ as to the same cells in XY", {
  F <- .mz_fix()
  ref <- .mz_join(F$pts, F$squares)
  expect_identical(ref$value$cell_id[order(ref$value$row)], 1:6)
  bound <- rbind(sf::st_sf(cell_id = c(1L, 2L, 4L, 5L, 6L),
                           geometry = sf::st_geometry(F$squares)[-3]),
                 sf::st_sf(cell_id = 3L, geometry = sf::st_geometry(F$mixed)[3]))
  expect_identical(.mz_dims(bound), c("XY", "XY", "XY", "XY", "XY", "XYZ"))
  cases <- list(direct = F$mixed, rbind = bound, reversed = F$mixed[6:1, ],
                xym = F$with3(3, "XYM"), first = F$with3(1), last = F$with3(6))
  for (nm in names(cases)) {
    r <- .mz_join(F$pts, cases[[nm]])
    expect_identical(r$classes, character(0), info = nm)
    expect_identical(r$value$cell_id[order(r$value$row)], 1:6, info = nm)
    expect_identical(attr(r$value, "ties")$n, 0L, info = nm)
    expect_identical(sf::st_geometry(r$value), sf::st_geometry(F$pts), info = nm)
  }
  # The same cells in the same order give the XY cells' result to the bit.
  expect_identical(.mz_join(F$pts, F$mixed)$value, ref$value)
  expect_identical(.mz_join(F$pts, F$with3(3, "XYM"))$value, ref$value)
})

test_that("an XYZ cell whose ring repeats a vertex in x and y among XY cells gets BUG-25's error", {
  F <- .mz_fix()
  m <- cbind(.mz_sq(-94.1, 39.6), 0)
  m <- rbind(m[1:2, ], c(m[2, 1:2], 5), m[3:5, ])
  g <- sf::st_geometry(F$squares); g[[3]] <- sf::st_polygon(list(m))
  rep3 <- sf::st_set_geometry(F$squares, g)
  # sf calls every row of the layer valid, and row 3 alone not.
  expect_true(all(suppressMessages(sf::st_is_valid(rep3))))
  expect_false(suppressMessages(sf::st_is_valid(rep3[3, ])))
  e <- tryCatch(assign_features_to_polygons(F$pts, rep3, "cell_id", quiet = TRUE),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_cells_degenerate")
  expect_identical(as.integer(e$rows), 3L)
})

test_that("features that mix XY and XYZ are joined as in XY and come back with their own geometry", {
  F <- .mz_fix()
  for (largest in c(TRUE, FALSE)) {
    ref <- .mz_join(F$patches, F$squares, largest = largest)
    r <- .mz_join(F$patches_z, F$squares, largest = largest)
    expect_identical(r$classes, ref$classes, info = largest)
    expect_identical(r$value$cell_id[order(r$value$row)], 1:6, info = largest)
    expect_identical(attr(r$value, "ties")$n, 0L, info = largest)
    expect_identical(sf::st_geometry(r$value), sf::st_geometry(F$patches_z), info = largest)
    expect_identical(.mz_dims(r$value)[2], "XYZ")
    expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(ref$value))
  }
  # Points that mix XY and XYZ, the XYZ point first or second.
  for (i in 1:2) {
    r <- .mz_join(F$pts_z(i), F$squares)
    expect_identical(r$classes, character(0))
    expect_identical(r$value$cell_id[order(r$value$row)], 1:6)
    expect_identical(sf::st_geometry(r$value), sf::st_geometry(F$pts_z(i)))
  }
})

test_that("projected joins with a layer that mixes XY and XYZ run as with the layer in XY", {
  F <- .mz_fix()
  sq32 <- sf::st_transform(F$squares, 32615)
  gu <- sf::st_geometry(sq32)
  gu[[3]] <- sf::st_polygon(list(cbind(gu[[3]][[1]], 7)))
  cells32 <- sf::st_set_geometry(sq32, gu)
  pu <- sf::st_transform(F$pts, 32615)
  cases <- list(
    both_32615   = list(pu, cells32, pu, sq32),
    lonlat_cells = list(pu, F$mixed, pu, F$squares),
    patches      = list(F$patches_z, sq32, F$patches, sq32),
    xyz_first    = list(F$pts_z(1), sq32, sf::st_zm(F$pts_z(1)), sq32))
  for (nm in names(cases)) {
    k <- cases[[nm]]
    r <- .mz_join(k[[1]], k[[2]])
    expect_identical(r$classes, character(0), info = nm)
    expect_identical(r$value$cell_id[order(r$value$row)], 1:6, info = nm)
    expect_identical(sf::st_geometry(r$value), sf::st_geometry(k[[1]]), info = nm)
    ref <- .mz_join(k[[3]], k[[4]])
    expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(ref$value),
                     info = nm)
  }
})

test_that("features without a CRS that mix dimensions keep them unless they are moved", {
  F <- .mz_fix()
  sq32 <- sf::st_transform(F$squares, 32615)
  # Coordinates in metres: the polygons' CRS is stamped on, and Z kept.
  p32 <- sf::st_coordinates(sf::st_transform(F$pts, 32615))
  g <- lapply(1:6, function(k) if (k == 2) sf::st_point(c(p32[k, ], 3)) else sf::st_point(p32[k, ]))
  pz32 <- sf::st_sf(row = 1:6, geometry = sf::st_sfc(g))
  r <- .mz_join(pz32, sq32)
  expect_false("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(r$value$cell_id, 1:6)
  expect_identical(sf::st_geometry(r$value), sf::st_set_crs(sf::st_geometry(pz32), 32615))
  # Neither layer has a CRS.
  r <- .mz_join(pz32, sf::st_set_crs(sq32, NA))
  expect_identical(r$classes, character(0))
  expect_identical(sf::st_geometry(r$value), sf::st_geometry(pz32))
  # Degrees taken as lon/lat and moved into the polygons' CRS: XY, with the warning.
  pz <- sf::st_set_crs(F$pts_z(2), NA)
  r <- .mz_join(pz, sq32)
  expect_true("spatialkit_zm_dropped" %in% r$classes)
  expect_match(r$texts[r$classes == "spatialkit_zm_dropped"],
               "assign_features_to_polygons(): `features_sf` mixes geometries with and without Z or M",
               fixed = TRUE)
  expect_identical(r$value$cell_id, 1:6)
  expect_identical(unique(.mz_dims(r$value)), "XY")
  ref <- suppressWarnings(assign_features_to_polygons(
    sf::st_zm(pz), sq32, "cell_id", quiet = TRUE))
  expect_identical(sf::st_geometry(r$value), sf::st_geometry(ref))
})

test_that("summarize_by_cell(area = TRUE) measures cells that mix XY and XYZ in x and y", {
  F <- .mz_fix()
  sumy <- function(pts, cells) {
    a <- assign_features_to_polygons(pts, cells, "cell_id", quiet = TRUE)
    .mz_catch(summarize_by_cell(a, "z", id_col = "cell_id", cells_sf = cells,
                                estimate = "mean", area = TRUE))
  }
  ref <- sumy(F$pts, F$squares)
  r <- sumy(F$pts, F$mixed)
  expect_identical(r$classes, ref$classes)
  expect_identical(r$value$cell_area, ref$value$cell_area)
  expect_equal(r$value$cell_area / 1e6, rep(95.2, 6), tolerance = 1e-3)
  expect_identical(r$value$n, rep(1L, 6))
  expect_identical(attr(r$value, "area_error"), attr(ref$value, "area_error"))
  expect_identical(sf::st_geometry(r$value), sf::st_geometry(F$mixed))
  # Projected cells: the XYZ cell among XY ones was refused as not measurable.
  sq32 <- sf::st_transform(F$squares, 32615)
  gu <- sf::st_geometry(sq32)
  gu[[3]] <- sf::st_polygon(list(cbind(gu[[3]][[1]], 7)))
  cells32 <- sf::st_set_geometry(sq32, gu)
  pu <- sf::st_transform(F$pts, 32615)
  ref <- sumy(pu, sq32)
  r <- sumy(pu, cells32)
  expect_identical(r$classes, ref$classes)
  expect_identical(r$value$cell_area, ref$value$cell_area)
  expect_identical(attr(r$value, "area_error"), attr(ref$value, "area_error"))
  expect_true(is.finite(attr(r$value, "area_error")))
  expect_identical(.mz_dims(r$value)[3], "XYZ")
})

test_that("ensure_projected() chooses from a mixed layer in XY and returns it in XY with a warning", {
  F <- .mz_fix()
  ref <- .mz_catch(ensure_projected(F$squares))
  r <- .mz_catch(ensure_projected(F$mixed))
  expect_identical(r$classes, c("spatialkit_zm_dropped", ref$classes))
  expect_match(r$texts[1], paste0(
    "ensure_projected(): `x` mixes geometries with and without Z or M coordinates ",
    "(1 of its 6 geometries carry them)"), fixed = TRUE)
  expect_identical(r$conds[[1]]$rows, 3L)
  expect_identical(sf::st_crs(r$value)$epsg, 32615L)
  expect_identical(r$value, ref$value)
  # The XYM row too, and an sfc.
  expect_identical(suppressWarnings(ensure_projected(F$with3(3, "XYM"))), ref$value)
  expect_identical(suppressWarnings(ensure_projected(sf::st_geometry(F$mixed))),
                   suppressWarnings(ensure_projected(sf::st_geometry(F$squares))))
  # A target CRS.
  ref <- .mz_catch(ensure_projected(F$squares, target_crs = 32615))
  r <- .mz_catch(ensure_projected(F$mixed, target_crs = 32615))
  expect_identical(r$classes, c("spatialkit_zm_dropped", ref$classes))
  expect_identical(r$value, ref$value)
  # Points with the XYZ point first stopped sf with "OGR error".
  r <- .mz_catch(ensure_projected(F$pts_z(1)))
  expect_true("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(r$value, suppressWarnings(ensure_projected(sf::st_zm(F$pts_z(1)))))
  # Without a CRS, taken as lon/lat.
  r <- .mz_catch(ensure_projected(sf::st_set_crs(F$mixed, NA)))
  expect_true("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(r$value, suppressWarnings(ensure_projected(sf::st_set_crs(F$squares, NA))))
  r <- .mz_catch(ensure_projected(sf::st_set_crs(F$mixed, NA), target_crs = 32615))
  expect_true("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(r$value, suppressWarnings(
    ensure_projected(sf::st_set_crs(F$squares, NA), target_crs = 32615)))
  # Projected input is returned as it came, and its area distortion measured.
  m32 <- sf::st_transform(sf::st_zm(F$mixed), 32615)
  gu <- sf::st_geometry(m32); gu[[3]] <- sf::st_polygon(list(cbind(gu[[3]][[1]], 7)))
  m32 <- sf::st_set_geometry(m32, gu)
  r <- .mz_catch(ensure_projected(m32, purpose = "area"))
  expect_identical(r$classes, character(0))
  expect_identical(r$value, m32)
  expect_true(is.finite(.crs_area_error(.xy_if_mixed(m32))))
})

test_that("harmonize_crs() returns a transformed mixed layer in XY with a warning", {
  F <- .mz_fix()
  pu <- sf::st_transform(F$pts, 32615)
  ref <- harmonize_crs(F$squares, pu, prefer = "b")
  r <- .mz_catch(harmonize_crs(F$mixed, pu, prefer = "b"))
  expect_identical(r$classes, "spatialkit_zm_dropped")
  expect_match(r$texts, "harmonize_crs(): `a` mixes geometries", fixed = TRUE)
  expect_identical(r$value, ref)
  r <- .mz_catch(harmonize_crs(pu, F$mixed, prefer = "a"))
  expect_match(r$texts, "harmonize_crs(): `b` mixes geometries", fixed = TRUE)
  expect_identical(r$value$b, harmonize_crs(pu, F$squares)$b)
  r <- .mz_catch(harmonize_crs(F$mixed, pu, target_crs = 32615))
  expect_identical(r$classes, "spatialkit_zm_dropped")
  expect_identical(r$value$a, ref$a)
  # A layer that is not transformed comes back as it came.
  r <- .mz_catch(harmonize_crs(F$mixed, F$pts, prefer = "b"))
  expect_identical(r$classes, character(0))
  expect_identical(r$value$a, F$mixed)
  r <- .mz_catch(harmonize_crs(F$mixed, pu, prefer = "a"))
  expect_identical(r$classes, character(0))
  expect_identical(r$value$a, F$mixed)
  # Without a CRS: the lon/lat test reads the layer in XY, so a mixed layer
  # in degrees is taken as lon/lat and moved, as the same layer in XY is.
  expect_true(.looks_like_lonlat(sf::st_set_crs(F$mixed, NA))$lonlat)
  ref <- .mz_catch(harmonize_crs(sf::st_set_crs(F$squares, NA), pu))
  r <- .mz_catch(harmonize_crs(sf::st_set_crs(F$mixed, NA), pu))
  expect_identical(r$classes, c(ref$classes, "spatialkit_zm_dropped"))
  expect_identical(r$value, ref$value)
  # Points with the XYZ point first.
  r <- .mz_catch(harmonize_crs(F$pts_z(1), pu, prefer = "b"))
  expect_identical(r$classes, "spatialkit_zm_dropped")
  expect_identical(r$value$a, harmonize_crs(sf::st_zm(F$pts_z(1)), pu, prefer = "b")$a)
})

test_that("the tessellation builders read a mixed boundary as the same boundary in XY", {
  F <- .mz_fix()
  seeds <- local({
    withr::local_seed(1, .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    sf::st_as_sf(data.frame(x = stats::runif(12, -94.29, -93.71),
                            y = stats::runif(12, 39.61, 39.69)),
                 coords = c("x", "y"), crs = 4326)
  })
  calls <- list(
    grid    = function(b) create_grid_polygons(b, target_cells = 12, type = "square"),
    gridsfc = function(b) create_grid_polygons(sf::st_geometry(b), target_cells = 12,
                                               type = "square"),
    voronoi = function(b) create_voronoi_polygons(seeds, boundary = b),
    tess_v  = function(b) build_tessellation(seeds, method = "voronoi", boundary = b),
    tess_h  = function(b) build_tessellation(seeds, method = "hex", boundary = b,
                                             approx_n_cells = 12),
    clip    = function(b) clip_target_for(seeds, boundary = b, quiet = TRUE))
  for (nm in names(calls)) {
    ref <- .mz_catch(calls[[nm]](F$squares))
    for (b in list(F$mixed, F$with3(1), F$with3(3, "XYM"))) {
      r <- .mz_catch(calls[[nm]](b))
      expect_identical(r$classes, ref$classes, info = nm)
      expect_identical(r$texts, ref$texts, info = nm)
      expect_identical(r$value, ref$value, info = nm)
    }
  }
  expect_identical(nrow(create_grid_polygons(F$mixed, target_cells = 12, type = "square",
                                             quiet = TRUE)), 16L)
})

test_that("layers in one dimension are read and transformed as before", {
  F <- .mz_fix()
  # Every cell XYZ: the join places every point, and ensure_projected()
  # keeps Z without the warning.
  r <- .mz_join(F$pts, F$allz)
  expect_identical(r$value$cell_id[order(r$value$row)], 1:6)
  r <- .mz_catch(ensure_projected(F$allz))
  expect_false("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(unique(.mz_dims(r$value)), "XYZ")
  expect_identical(sf::st_crs(r$value)$epsg, 32615L)
  # A CRS-less XYZ point layer keeps its Z through the join.
  p32 <- sf::st_coordinates(sf::st_transform(F$pts, 32615))
  pz <- sf::st_sf(row = 1:6, geometry = sf::st_sfc(lapply(1:6, function(k)
    sf::st_point(c(p32[k, ], k)))))
  r <- .mz_join(pz, sf::st_transform(F$squares, 32615))
  expect_identical(unique(.mz_dims(r$value)), "XYZ")
  expect_identical(r$value$cell_id, 1:6)
  # The helpers hand a layer in one dimension back as it came.
  expect_identical(.xy_if_mixed_warn(F$allz, "f", "x"), F$allz)
  expect_identical(.xy_if_mixed_warn(F$squares, "f", "x"), F$squares)
  expect_no_warning(.xy_if_mixed_warn(F$allz, "f", "x"))
  r <- .mz_catch(.xy_if_mixed_warn(F$mixed, "f", "x"))
  expect_identical(r$classes, "spatialkit_zm_dropped")
  expect_identical(r$value, sf::st_zm(F$mixed))
})

test_that("the help says mixed layers are read in x and y", {
  # The source page from the source tree, the installed page under R CMD
  # check, which runs the tests where there is no man/.
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  txt <- page("assign_features_to_polygons.Rd")
  expect_match(txt, paste0("When either layer mixes dimensions, the join runs on a copy of it ",
                           "without Z and M, and the rows come back with the geometry they ",
                           "arrived with, Z and M included."), fixed = TRUE)
  txt <- page("summarize_by_cell.Rd")
  expect_match(txt, paste0("Cells that mix XY and XYZ (or XYM) geometries are measured, and ",
                           "their spread taken, without Z and M"), fixed = TRUE)
  for (nm in c("ensure_projected.Rd", "harmonize_crs.Rd"))
    expect_match(page(nm), "spatialkit_zm_dropped", fixed = TRUE)
  # A boundary that mixes dimensions is read as ?spatialkit-geometry says,
  # in x and y.
  for (nm in c("clip_target_for.Rd", "create_voronoi_polygons.Rd",
               "create_grid_polygons.Rd", "build_tessellation.Rd"))
    expect_match(page(nm), "mixed dimensions are handled as \"Boundaries\" in spatialkit-geometry",
                 fixed = TRUE, info = nm)
  expect_match(page("spatialkit-geometry.Rd"),
               "A boundary that mixes dimensions or carries M is read in x and y.", fixed = TRUE)
})

test_that("the scored-cells path joins mixed cells and mixed points as in XY", {
  # Cells built in lon/lat from a profile's centres carry the record that
  # moves the join into the profile's CRS: the cells are moved from their
  # copy in XY, and the points on an edge are read from theirs.
  g <- expand.grid(x = 4e5 + (0:19) * 30, y = 4.4e6 + (0:14) * 30)
  z <- local({
    withr::local_seed(9001, .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    sin(g$x / 180) + cos(g$y / 150) + stats::rnorm(nrow(g), 0, 0.4)
  })
  pts <- sf::st_as_sf(data.frame(x = g$x, y = g$y, z = z, row = seq_len(nrow(g))),
                      coords = c("x", "y"), crs = 32615)
  prof <- suppressWarnings(resolution_profile(pts, "z", partition = "kmeans", nstart = 3,
                                              sac = "auto", seed = 123, estimate = "pooled",
                                              obs_per_cell = 4.5))
  sel <- select_resolution(prof, criterion = "size", obs_per_cell = 4.5)
  seeds <- suppressWarnings(get_voronoi_seeds(method = "kmeans", n = sel, sample_points = pts))
  cells <- build_tessellation(seeds, method = "voronoi", crs = 4326, quiet = TRUE)$cells
  expect_true(sf::st_is_longlat(cells))
  expect_true(.scored_crs(.count_rule_get(cells)) == sf::st_crs(32615))
  ref <- .mz_join(pts, cells)
  expect_false(anyNA(ref$value$cell_id))
  # Cell 2 given Z, its record kept.
  gc <- sf::st_geometry(cells)
  gc[[2]] <- sf::st_polygon(lapply(gc[[2]], function(r) cbind(r, 5)))
  mc <- sf::st_set_geometry(cells, gc)
  attr(mc, "count_rule") <- attr(cells, "count_rule")
  expect_identical(.mz_dims(mc)[1:3], c("XY", "XYZ", "XY"))
  expect_true(.scored_crs(.count_rule_get(mc)) == sf::st_crs(32615))
  r <- .mz_join(pts, mc)
  expect_identical(r$classes, ref$classes)
  expect_identical(r$value, ref$value)
  # Points that mix XY and XYZ, the XYZ point first or in row 5.
  for (i in c("1", "5")) {
    gp <- sf::st_geometry(pts)
    k <- as.integer(i)
    gp[[k]] <- sf::st_point(c(as.numeric(gp[[k]]), 2))
    mp <- sf::st_set_geometry(pts, sf::st_sfc(lapply(gp, identity), crs = 32615))
    r <- .mz_join(mp, cells)
    expect_identical(r$classes, ref$classes, info = i)
    expect_identical(sf::st_drop_geometry(r$value), sf::st_drop_geometry(ref$value), info = i)
    expect_identical(attr(r$value, "ties"), attr(ref$value, "ties"), info = i)
    expect_identical(sf::st_geometry(r$value), sf::st_geometry(mp), info = i)
  }
})

test_that("ensure_projected(purpose = \"area\") measures a mixed projected layer's distortion in XY", {
  sq5 <- function(y0) sf::st_polygon(list(.mz_sq(-100, y0, 5)))
  m <- sf::st_transform(sf::st_sf(id = 1:6, geometry = sf::st_sfc(
    lapply(30 + 5 * (0:5), sq5), crs = 4326)), 3857)
  gm <- sf::st_geometry(m)
  gm[[2]] <- sf::st_polygon(lapply(gm[[2]], function(r) cbind(r, 1)))
  mixed <- sf::st_set_geometry(m, gm)
  e <- .crs_area_error(m)
  expect_gt(e, .area_error_tol)
  line <- sprintf("distorts areas across this extent by up to %.1f%%", 100 * e)
  lines_xy <- capture_spatialkit_log(out_xy <- ensure_projected(m, purpose = "area"))
  expect_true(log_has(lines_xy, line))
  r <- .mz_catch(lines <- capture_spatialkit_log(out <- ensure_projected(mixed, purpose = "area")))
  expect_identical(r$classes, character(0))
  expect_true(log_has(lines, line))
  expect_identical(out, mixed)
  expect_identical(out_xy, m)
})

test_that("the warning names the rows that carry Z or M, XYM rows included", {
  F <- .mz_fix()
  for (dim in c("XYZ", "XYM")) {
    r <- .mz_catch(ensure_projected(F$with3(3, dim)))
    w <- r$conds[r$classes == "spatialkit_zm_dropped"]
    expect_length(w, 1L)
    expect_identical(w[[1]]$rows, 3L, info = dim)
    expect_match(conditionMessage(w[[1]]), "(1 of its 6 geometries carry them)",
                 fixed = TRUE, info = dim)
  }
  # A layer that mixes XYZ and XYM with no XY row: every row is counted.
  g <- lapply(0:5, function(i) sf::st_polygon(list(cbind(.mz_sq(-94.3 + 0.1 * i, 39.6), 7)),
                                              dim = if (i == 2) "XYM" else "XYZ"))
  zm <- sf::st_sf(cell_id = 1:6, geometry = sf::st_sfc(g, crs = 4326))
  r <- .mz_catch(ensure_projected(zm))
  w <- r$conds[r$classes == "spatialkit_zm_dropped"]
  expect_identical(w[[1]]$rows, 1:6)
  expect_identical(r$value, suppressWarnings(ensure_projected(F$squares)))
})

test_that("an empty geometry of another dimension after the first leaves Z in place", {
  F <- .mz_fix()
  pu <- sf::st_transform(F$pts, 32615)
  # sf transforms an XYZ layer with an empty XY row after the first with Z.
  ez <- rbind(F$allz, sf::st_sf(cell_id = 7L, geometry = sf::st_sfc(sf::st_polygon(), crs = 4326)))
  expect_identical(.mz_dims(ez), c(rep("XYZ", 6), "XY"))
  ref <- suppressWarnings(ensure_projected(F$allz))
  r <- .mz_catch(ensure_projected(ez))
  expect_false("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(unique(.mz_dims(r$value)[1:6]), "XYZ")
  expect_identical(sf::st_crs(r$value)$epsg, 32615L)
  expect_identical(sf::st_geometry(r$value)[1:6], sf::st_geometry(ref))
  r <- .mz_catch(harmonize_crs(ez, pu, prefer = "b"))
  expect_identical(r$classes, character(0))
  expect_identical(sf::st_geometry(r$value$a)[1:6],
                   sf::st_geometry(harmonize_crs(F$allz, pu, prefer = "b")$a))
  # The empty XY row first: sf reads the rest as XY pairs, so the layer is
  # handed over in XY with the warning.
  ez1 <- rbind(sf::st_sf(cell_id = 7L, geometry = sf::st_sfc(sf::st_polygon(), crs = 4326)),
               F$allz)
  r <- .mz_catch(ensure_projected(ez1))
  w <- r$conds[r$classes == "spatialkit_zm_dropped"]
  expect_identical(w[[1]]$rows, 2:7)
  expect_identical(sf::st_geometry(r$value)[2:7], sf::st_geometry(sf::st_zm(ref)))
  # An empty point holds NA coordinates, which sf misreads: XYZ points with
  # an empty XY point after them stopped sf::st_transform() with "OGR error",
  # and are handed over in XY with the warning.
  ep <- rbind(F$pts_z(1)[, "row"], sf::st_sf(row = 7L, geometry = sf::st_sfc(sf::st_point(), crs = 4326)))
  gp <- sf::st_geometry(ep)
  for (k in 2:6) gp[[k]] <- sf::st_point(c(F$xy[k, ], 3))
  ep <- sf::st_set_geometry(ep, gp)
  expect_identical(.mz_dims(ep), c(rep("XYZ", 6), "XY"))
  expect_error(suppressWarnings(sf::st_transform(ep, 32615)))
  r <- .mz_catch(ensure_projected(ep))
  w <- r$conds[r$classes == "spatialkit_zm_dropped"]
  expect_identical(w[[1]]$rows, 1:6)
  expect_identical(unique(.mz_dims(r$value)), "XY")
  # Features without a CRS, in degrees, moved into the polygons' CRS: the
  # rows keep their Z, and the empty row is left out as it has no location.
  pz <- sf::st_sf(row = 1:7, geometry = sf::st_sfc(c(
    lapply(1:6, function(k) sf::st_polygon(list(cbind(
      .mz_sq(F$xy[k, 1] - 0.01, F$xy[k, 2] - 0.01, 0.02), k)))),
    list(sf::st_polygon()))))
  expect_identical(.mz_dims(pz), c(rep("XYZ", 6), "XY"))
  sq32 <- sf::st_transform(F$squares, 32615)
  r <- .mz_join(pz, sq32)
  expect_false("spatialkit_zm_dropped" %in% r$classes)
  expect_identical(r$value$row, 1:6)
  expect_identical(r$value$cell_id, 1:6)
  expect_identical(unique(.mz_dims(r$value)), "XYZ")
  ref <- .mz_join(pz[1:6, ], sq32)
  expect_identical(sf::st_geometry(r$value), sf::st_geometry(ref$value))
})

test_that("the lon/lat test of a layer without a CRS leaves empty geometries out", {
  # sf::st_coordinates() stops on a polygon layer with an empty polygon, so
  # a layer in degrees a fraction of a degree across was taken as planar.
  # An empty point, whose coordinates are NA, did no harm, and does none.
  F <- .mz_fix()
  sq <- sf::st_set_crs(F$squares, NA)
  expect_true(.looks_like_lonlat(sq)$lonlat)
  e_poly <- rbind(sq, sf::st_sf(cell_id = 7L, geometry = sf::st_sfc(sf::st_polygon())))
  expect_true(.looks_like_lonlat(e_poly)$lonlat)
  expect_true(.looks_like_lonlat(sf::st_geometry(e_poly))$lonlat)
  pts <- sf::st_set_crs(F$pts[, "row"], NA)
  e_pt <- rbind(sf::st_sf(row = 0L, geometry = sf::st_sfc(sf::st_point())), pts)
  expect_true(.looks_like_lonlat(e_pt)$lonlat)
  # ensure_projected() takes the layer as lon/lat, as the layer without the
  # empty row.
  ref <- suppressWarnings(ensure_projected(sq))
  r <- suppressWarnings(ensure_projected(e_poly))
  expect_identical(sf::st_crs(r)$epsg, 32615L)
  expect_identical(sf::st_geometry(r)[1:6], sf::st_geometry(ref))
  # A planar layer in metres is still taken as planar.
  m <- sf::st_set_crs(sf::st_transform(F$squares, 32615), NA)
  m <- rbind(m, sf::st_sf(cell_id = 7L, geometry = sf::st_sfc(sf::st_polygon())))
  expect_false(.looks_like_lonlat(m)$lonlat)
})
