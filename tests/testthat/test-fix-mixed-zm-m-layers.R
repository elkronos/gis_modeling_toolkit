# tests/testthat/test-fix-mixed-zm-m-layers.R
# ---------------------------------------------------------------------------
# GEOS does not read M: it stops with "GEOS does not support XYM or XYZM
# geometries" on any geometry that carries one.  assign_features_to_polygons()
# joined a layer whose geometries all carry M as it came, and so stopped
# wherever the join runs with GEOS; cell_values_at() stopped in that join and
# in its search of the cells no seed names, and kriging_adequacy() in the
# repair of such cells.  The join and the cells are defined in x and y and M
# is read nowhere, so these layers are now read in x and y, with no warning,
# and every result is held to the same rows in XY, to the bit.  The rows
# assign_features_to_polygons() returns keep their own geometry, M included.
#
# The layer: 300 points in EPSG:32617 (seed 6102), the Voronoi cells of 20
# k-means seeds of 240 of them, their pooled summary, and squares of 6 km
# around 60 of the points, which straddle the cells' edges.
# ---------------------------------------------------------------------------

.zmm_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)))
}

.zmm_store <- new.env(parent = emptyenv())
.zmm_layer <- function() {
  if (is.null(.zmm_store$layer)) {
    d <- local({
      withr::local_seed(6102, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 300
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      z <- sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 0.6)
      te <- sample.int(n, 60)
      list(d = data.frame(x = 5e5 + x, y = 5e6 + y, z = z, ..r = seq_len(n)), te = te)
    })
    pts  <- sf::st_as_sf(d$d, coords = c("x", "y"), crs = 32617)
    rows <- pts[-d$te, ]
    seeds <- get_voronoi_seeds(method = "kmeans", n = 20, sample_points = rows,
                               set_seed = 1)
    tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    asg   <- assign_features_to_polygons(rows, tess$cells, polygon_id_col = "cell_id",
                                         quiet = TRUE)
    pooled <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                estimate = "pooled")
    seeds$cell_id <- tess$index
    held <- pts[d$te, ]
    held$..r <- seq_len(nrow(held))
    xy <- sf::st_coordinates(held)
    squares <- sf::st_sf(..r = held$..r, z = held$z, geometry = sf::st_sfc(
      lapply(seq_len(nrow(xy)), function(i) {
        x0 <- xy[i, 1] - 3000; y0 <- xy[i, 2] - 3000
        sf::st_polygon(list(rbind(c(x0, y0), c(x0 + 6000, y0), c(x0 + 6000, y0 + 6000),
                                  c(x0, y0 + 6000), c(x0, y0))))
      }), crs = 32617))
    .zmm_store$layer <- list(held = held, squares = squares, cells = tess$cells,
                             seeds = seeds, asg = asg, pooled = pooled)
  }
  .zmm_store$layer
}

# Every geometry of `l` given M at 3 (dim = "XYM"), or Z at 3 and M at 4
# (dim = "XYZM"), or Z at 3 (dim = "XYZ"), at its own x and y.  Empty
# geometries become empty geometries of that dimension.
.zmm_all <- function(l, dim) {
  extra <- switch(dim, XYM = 3, XYZ = 3, XYZM = c(3, 4))
  add <- function(r) do.call(cbind, c(list(r), as.list(extra)))
  na  <- rep(NA_real_, 2L + length(extra))
  g <- lapply(sf::st_geometry(l), function(s) {
    co <- unclass(s)
    switch(class(s)[2L],
      POINT = if (sf::st_is_empty(s)) sf::st_point(na, dim = dim)
              else sf::st_point(c(co, extra), dim = dim),
      POLYGON = sf::st_polygon(lapply(co, add), dim = dim),
      MULTIPOLYGON = sf::st_multipolygon(lapply(co, function(p) lapply(p, add)), dim = dim))
  })
  sf::st_set_geometry(l, sf::st_sfc(g, crs = sf::st_crs(l), recompute_bbox = TRUE))
}

.zmm_dims <- function(x)
  unique(vapply(sf::st_geometry(x), function(g) class(g)[1L], character(1)))

.zmm_join <- function(f, p, ...)
  .zmm_catch(assign_features_to_polygons(f, p, "cell_id", quiet = TRUE,
                                         keep_unassigned = TRUE, ...))

# The result of a join on rows that carry M against the same rows in XY:
# identical once the geometry is put aside, warnings included, and the rows
# returned carry `geom`, the caller's geometry, row for row.
.zmm_same_join <- function(got, ref, geom, info) {
  expect_identical(sf::st_drop_geometry(got$value), sf::st_drop_geometry(ref$value),
                   info = info)
  for (a in c("unassigned", "ties", "count_rule"))
    expect_identical(attr(got$value, a), attr(ref$value, a), info = paste(info, a))
  expect_identical(got$classes, ref$classes, info = info)
  expect_identical(got$texts, ref$texts, info = info)
  expect_false("spatialkit_zm_dropped" %in% got$classes, info = info)
  expect_identical(sf::st_geometry(got$value), geom[got$value$..r], info = info)
  expect_identical(sf::st_crs(got$value), sf::st_crs(ref$value), info = info)
}

test_that("the premise: GEOS stops on geometries that carry M, and sf transforms them with M", {
  L <- .zmm_layer()
  m <- .zmm_all(L$held, "XYM")
  expect_identical(.zmm_dims(m), "XYM")
  expect_error(sf::st_intersects(m, L$cells), "GEOS does not support XYM or XYZM")
  moved <- sf::st_transform(m, 4326)
  expect_identical(.zmm_dims(moved), "XYM")
  expect_identical(sf::st_geometry(sf::st_zm(moved)),
                   sf::st_geometry(sf::st_transform(L$held, 4326)))
})

test_that("features_sf whose geometries carry M are joined as the same rows in XY", {
  L <- .zmm_layer()
  for (crs in c(32617L, 4326L)) for (kind in c("points", "squares")) {
    base <- if (kind == "points") L$held else L$squares
    if (crs != 32617L) base <- sf::st_transform(base, crs)
    for (lg in c(TRUE, FALSE)) {
      ref <- .zmm_join(base, L$cells, largest = lg)
      for (dim in c("XYM", "XYZM")) {
        f <- .zmm_all(base, dim)
        .zmm_same_join(.zmm_join(f, L$cells, largest = lg), ref, sf::st_geometry(f),
                       sprintf("%s %s in EPSG:%d, largest = %s", kind, dim, crs, lg))
      }
    }
  }
})

test_that("polygons_sf whose geometries carry M are read as the same polygons in XY", {
  L <- .zmm_layer()
  cells_ll <- sf::st_transform(L$cells, 4326)
  for (cells in list(L$cells, cells_ll)) for (kind in c("points", "squares")) {
    base <- if (kind == "points") L$held else L$squares
    for (fc in c(32617L, 4326L)) {
      f <- if (fc == 32617L) base else sf::st_transform(base, fc)
      for (lg in c(TRUE, FALSE)) {
        ref <- .zmm_join(f, cells, largest = lg)
        for (dim in c("XYM", "XYZM"))
          .zmm_same_join(.zmm_join(f, .zmm_all(cells, dim), largest = lg), ref,
                         sf::st_geometry(f),
                         sprintf("cells %s in EPSG:%s, %s in EPSG:%d, largest = %s", dim,
                                 sf::st_crs(cells)$epsg, kind, fc, lg))
      }
    }
  }
})

test_that("features without a CRS that carry M come back with their own geometry", {
  L <- .zmm_layer()
  # Given the cells' CRS: the coordinates are metres and are only stamped.
  bare <- sf::st_set_crs(L$squares, NA)
  ref <- .zmm_join(bare, L$cells)
  expect_identical(sf::st_geometry(ref$value),
                   sf::st_set_crs(sf::st_geometry(bare), 32617)[ref$value$..r])
  for (dim in c("XYM", "XYZM")) {
    f <- .zmm_all(bare, dim)
    .zmm_same_join(.zmm_join(f, L$cells), ref,
                   sf::st_set_crs(sf::st_geometry(f), 32617), paste("stamped", dim))
  }
  # In degrees: taken as EPSG:4326 and moved into the cells' CRS, M kept.
  for (base in list(L$held, L$squares)) {
    bare_ll <- sf::st_set_crs(sf::st_transform(base, 4326), NA)
    ref <- .zmm_join(bare_ll, L$cells)
    for (dim in c("XYM", "XYZM")) {
      f <- .zmm_all(bare_ll, dim)
      got <- .zmm_join(f, L$cells)
      moved <- sf::st_transform(sf::st_set_crs(sf::st_geometry(f), 4326), 32617)
      info <- paste("moved from lon/lat", dim, class(sf::st_geometry(base))[1L])
      .zmm_same_join(got, ref, moved, info)
      expect_identical(.zmm_dims(got$value), dim, info = info)
      expect_identical(sf::st_coordinates(sf::st_zm(got$value)),
                       sf::st_coordinates(ref$value), info = info)
    }
  }
})

test_that("layers that carry M join with GEOS when neither has a CRS", {
  L <- .zmm_layer()
  bare_c <- sf::st_set_crs(L$cells, NA)
  for (base in list(L$held, L$squares)) {
    bare <- sf::st_set_crs(base, NA)
    for (lg in c(TRUE, FALSE)) {
      ref <- .zmm_join(bare, bare_c, largest = lg)
      for (dim in c("XYM", "XYZM")) {
        f <- .zmm_all(bare, dim)
        info <- sprintf("no CRS, %s %s, largest = %s", class(sf::st_geometry(base))[1L],
                        dim, lg)
        .zmm_same_join(.zmm_join(f, bare_c, largest = lg), ref, sf::st_geometry(f), info)
        .zmm_same_join(.zmm_join(bare, .zmm_all(bare_c, dim), largest = lg), ref,
                       sf::st_geometry(bare), paste("cells", info))
      }
    }
  }
})

test_that("layers that carry M join with GEOS in degrees with s2 off", {
  # The join then runs with GEOS on lon/lat taken as planar.  Areas in
  # degrees with s2 off need lwgeom, which is not a dependency: the ties are
  # broken by "first", which reads no area, and the largest overlap of
  # squares runs only where sf can measure an area in degrees with s2 off.
  L <- .zmm_layer()
  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  cells_ll <- sf::st_transform(L$cells, 4326)
  has_area <- tryCatch({
    sf::st_area(cells_ll[1L, ])
    TRUE
  }, error = function(e) FALSE)
  for (base in list(L$held, L$squares)) {
    f0 <- sf::st_transform(base, 4326)
    for (lg in c(TRUE, FALSE)) {
      if (lg && identical(base, L$squares) && !has_area) next
      ref <- .zmm_join(f0, cells_ll, largest = lg, tie_break = "first")
      for (dim in c("XYM", "XYZM")) {
        f <- .zmm_all(f0, dim)
        info <- sprintf("s2 off, %s %s, largest = %s", class(sf::st_geometry(base))[1L],
                        dim, lg)
        .zmm_same_join(.zmm_join(f, cells_ll, largest = lg, tie_break = "first"), ref,
                       sf::st_geometry(f), info)
        .zmm_same_join(.zmm_join(f0, .zmm_all(cells_ll, dim), largest = lg,
                                 tie_break = "first"), ref,
                       sf::st_geometry(f0), paste("cells", info))
      }
    }
  }
})

test_that("an empty point among points that carry M is joined as the same empty point in XY", {
  L <- .zmm_layer()
  base <- L$held
  g <- sf::st_geometry(base)
  g[[4]] <- sf::st_point()
  base <- sf::st_set_geometry(base, sf::st_sfc(as.list(g), crs = 32617,
                                               recompute_bbox = TRUE))
  ref <- .zmm_join(base, L$cells)
  expect_true(4L %in% attr(ref$value, "unassigned")$which)
  for (dim in c("XYM", "XYZM")) {
    f <- .zmm_all(base, dim)
    expect_identical(.zmm_dims(f), dim)
    .zmm_same_join(.zmm_join(f, L$cells), ref, sf::st_geometry(f), paste("projected", dim))
    # Without a CRS, in degrees: moved with the empty point.
    bare_ll <- sf::st_set_crs(sf::st_transform(base, 4326), NA)
    fl <- .zmm_all(bare_ll, dim)
    .zmm_same_join(.zmm_join(fl, L$cells), .zmm_join(bare_ll, L$cells),
                   sf::st_transform(sf::st_set_crs(sf::st_geometry(fl), 4326), 32617),
                   paste("moved from lon/lat", dim))
  }
})

test_that("features that are all empty polygons in M are joined as the same empty polygons in XY", {
  # GEOS refuses an empty polygon in M as well.
  L <- .zmm_layer()
  expect_error(sf::st_is_valid(sf::st_sfc(sf::st_polygon(dim = "XYM"), crs = 32617)),
               "GEOS does not support XYM or XYZM")
  base <- sf::st_sf(..r = 1:3, geometry = sf::st_sfc(rep(list(sf::st_polygon()), 3),
                                                     crs = 32617))
  for (lg in c(TRUE, FALSE)) {
    ref <- .zmm_join(base, L$cells, largest = lg)
    expect_identical(attr(ref$value, "unassigned")$which, 1:3)
    for (dim in c("XYM", "XYZM")) {
      f <- sf::st_set_geometry(base, sf::st_sfc(rep(list(sf::st_polygon(dim = dim)), 3),
                                                crs = 32617))
      expect_identical(.zmm_dims(f), dim)
      .zmm_same_join(.zmm_join(f, L$cells, largest = lg), ref, sf::st_geometry(f),
                     paste("empty", dim, "largest =", lg))
    }
  }
})

test_that("layers XYZ throughout are joined as they came and keep their Z", {
  L <- .zmm_layer()
  for (base in list(L$held, L$squares)) {
    ref <- .zmm_join(base, L$cells)
    f <- .zmm_all(base, "XYZ")
    got <- .zmm_join(f, L$cells)
    .zmm_same_join(got, ref, sf::st_geometry(f), "XYZ")
    expect_identical(.zmm_dims(got$value), "XYZ")
    got_p <- .zmm_join(base, .zmm_all(L$cells, "XYZ"))
    .zmm_same_join(got_p, ref, sf::st_geometry(base), "cells XYZ")
  }
})

.zmm_at <- function(L, newdata, values = L$pooled, seeds = NULL)
  .zmm_catch(cell_values_at(newdata, values, L$asg, "z", id_col = "cell_id",
                            seeds = seeds))

test_that("cell_values_at() reads newdata and values that carry M as the same rows in XY", {
  L <- .zmm_layer()
  nd <- L$held[1:30, ]
  # Point newdata without seeds, placed by the polygons of `values`.
  ref <- .zmm_at(L, nd)
  for (dim in c("XYM", "XYZM")) {
    got <- .zmm_at(L, .zmm_all(nd, dim))
    expect_identical(got, ref, info = paste("newdata", dim))
  }
  # `values` whose polygons carry M, without seeds.
  for (dim in c("XYM", "XYZM")) {
    got <- .zmm_at(L, nd, values = .zmm_all(L$pooled, dim))
    expect_identical(got, ref, info = paste("values", dim))
    expect_false("spatialkit_zm_dropped" %in% got$classes)
  }
  # A seed without an ID inside its cell, and one moved beyond every cell,
  # whose cell then has no seed and is offered to it (.values_at_free_cells()).
  inside <- L$seeds
  inside$cell_id[3] <- NA
  beyond <- inside
  g <- sf::st_geometry(beyond)
  g[3] <- g[3] + c(1e6, 0)
  beyond <- sf::st_set_geometry(beyond, sf::st_set_crs(g, 32617))
  for (sd in list(inside, beyond)) {
    ref <- .zmm_at(L, nd, seeds = sd)
    for (dim in c("XYM", "XYZM")) {
      got <- .zmm_at(L, nd, values = .zmm_all(L$pooled, dim), seeds = sd)
      expect_identical(got, ref, info = paste("values with seeds", dim))
    }
  }
})

test_that("kriging_adequacy() reads cells_sf that carry M as the same cells in XY", {
  L <- .zmm_layer()
  sac <- suppressWarnings(estimate_sac_range(L$asg, "z", cutoff = 1))
  for (cells in list(L$cells, sf::st_transform(L$cells, 4326))) {
    ref <- .zmm_catch(kriging_adequacy(L$asg, "z", cells, id_col = "cell_id", sac = sac,
                                       seed = 1))
    for (dim in c("XYM", "XYZM")) {
      got <- .zmm_catch(kriging_adequacy(L$asg, "z", .zmm_all(cells, dim), id_col = "cell_id",
                                         sac = sac, seed = 1))
      expect_identical(got, ref, info = sprintf("cells %s in EPSG:%s", dim,
                                                sf::st_crs(cells)$epsg))
    }
  }
})

test_that("the help says layers that carry M are read in x and y", {
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
  expect_match(txt, paste0("The same holds, with no warning, for a layer whose geometries carry ",
                           "M (XYM or XYZM), which GEOS does not read"), fixed = TRUE)
  expect_match(txt, paste0("A layer whose geometries mix dimensions or carry M (XYM or XYZM) is ",
                           "joined in x and y, and the rows returned keep their Z and M"),
               fixed = TRUE)
  expect_match(txt, paste0("Polygons whose geometries mix dimensions or carry M (XYM or XYZM) ",
                           "are read in x and y."), fixed = TRUE)
  txt <- page("cell_values_at.Rd")
  expect_match(txt, paste0("A layer whose geometries mix dimensions or carry M (XYM or XYZM) is ",
                           "read in x and y"), fixed = TRUE)
  expect_match(txt, paste0("Its polygons are read in x and y when they mix dimensions or carry ",
                           "M (XYM or XYZM)."), fixed = TRUE)
  expect_match(txt, paste0("Polygons of values whose geometries carry M are read in x and y ",
                           "as well."), fixed = TRUE)
  expect_match(page("kriging_adequacy.Rd"),
               "So are cells whose geometries carry M (XYM or XYZM), which GEOS does not read.",
               fixed = TRUE)
})
