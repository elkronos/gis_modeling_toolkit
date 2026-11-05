# tests/testthat/test-fix-mixed-zm-points.R
# ---------------------------------------------------------------------------
# sf reads a point layer in the dimension of its first point:
# st_coordinates() of a layer that mixes XY and XYZ points lays the values out
# of step (XYZ first) or stops on its "dimnames" (XY first), and
# st_transform() of one with an XYZ point first stops with "OGR error".
# rbind() of an XY and an XYZ layer builds such a layer without a word.
# cell_values_at() reads `newdata`, `seeds` and `points`, and
# summarize_by_cell(estimate = "pooled") reads its rows, in x and y: Z and M
# are dropped from the copy read when a layer mixes dimensions.  Each mixed
# layer is held to the same rows in XY, to the bit, warnings included.
#
# The layer: 300 locations (seed 6102), 60 held out, 20 k-means Voronoi
# cells, pooled values; the working CRS is EPSG:32617.
# ---------------------------------------------------------------------------

.zm_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)))
}

.zm_store <- new.env(parent = emptyenv())
.zm_layer <- function() {
  if (is.null(.zm_store$layer)) {
    d <- local({
      withr::local_seed(6102, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 300
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      z <- sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 0.6)
      te <- sample.int(n, 60)
      list(d = data.frame(x = 5e5 + x, y = 5e6 + y, z = z), te = te)
    })
    pts  <- sf::st_as_sf(d$d, coords = c("x", "y"), crs = 32617)
    rows <- pts[-d$te, ]
    seeds <- get_voronoi_seeds(method = "kmeans", n = 20, sample_points = rows, seed = 1)
    tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    asg   <- assign_features_to_polygons(rows, tess$cells, polygon_id_col = "cell_id",
                                         quiet = TRUE)
    pooled <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                estimate = "pooled")
    seeds$cell_id <- tess$index
    .zm_store$layer <- list(held = pts[d$te, ][1:8, ], tess = tess, seeds = seeds,
                            asg = asg, pooled = pooled)
  }
  .zm_store$layer
}

# Row `i` of `l` given a third coordinate at its own x and y: Z, or M.
.zm_row <- function(l, i, dim = "XYZ") {
  g  <- sf::st_geometry(l)
  xy <- sf::st_coordinates(l)[i, 1:2]
  g[i] <- sf::st_sfc(sf::st_point(c(xy, 3), dim = dim), crs = sf::st_crs(l))
  sf::st_set_geometry(l, g)
}

.zm_at <- function(L, newdata, points = L$asg, seeds = L$seeds)
  .zm_catch(cell_values_at(newdata, L$pooled, points, "z", id_col = "cell_id",
                           seeds = seeds))

test_that("the premise: sf reads a mixed point layer in the dimension of its first point", {
  xy_first  <- sf::st_sfc(sf::st_point(c(1, 2)), sf::st_point(c(3, 4, 5)),
                          sf::st_point(c(6, 7)), crs = 32617)
  xyz_first <- sf::st_sfc(sf::st_point(c(3, 4, 5)), sf::st_point(c(1, 2)),
                          sf::st_point(c(6, 7)), crs = 32617)
  expect_error(suppressWarnings(sf::st_coordinates(xy_first)))
  got <- suppressWarnings(sf::st_coordinates(xyz_first))
  expect_false(identical(unname(got[, 1:2]), rbind(c(3, 4), c(1, 2), c(6, 7))))
  # GDAL prints "OGR: Not enough data" on the way.
  utils::capture.output(expect_error(sf::st_transform(xyz_first, 32616)))
})

test_that("newdata that mixes XY and XYZ points gives the frame of the same rows in XY", {
  L <- .zm_layer()
  for (by in c("seeds", "polygons")) {
    sd <- if (by == "seeds") L$seeds else NULL
    for (epsg in c(32617, 32616)) {
      nd  <- if (epsg == 32617) L$held else sf::st_transform(L$held, epsg)
      ref <- .zm_at(L, nd, seeds = sd)
      for (i in c(1L, 3L)) for (dim in c("XYZ", "XYM")) {
        lab <- sprintf("%s, EPSG:%d, row %d %s", by, epsg, i, dim)
        got <- .zm_at(L, .zm_row(nd, i, dim), seeds = sd)
        expect_identical(got$value, ref$value, label = lab)
        expect_identical(got$classes, ref$classes, label = lab)
        expect_identical(got$texts, ref$texts, label = lab)
      }
    }
  }
})

test_that("a mixed row of newdata beyond 1e150 still has NA in every column", {
  L <- .zm_layer()
  nd <- sf::st_transform(L$held, 32616)
  g <- sf::st_geometry(nd)
  g[2] <- sf::st_sfc(sf::st_point(c(5e5, 1e200)), crs = 32616)
  far <- sf::st_set_geometry(nd, g)
  for (sd in list(L$seeds, NULL)) {
    ref <- .zm_at(L, far, seeds = sd)
    for (i in c(1L, 3L)) {
      got <- .zm_at(L, .zm_row(far, i), seeds = sd)
      expect_identical(got$value, ref$value)
      expect_identical(got$classes, ref$classes)
      expect_identical(got$texts, ref$texts)
    }
    mixed_far <- far
    g2 <- sf::st_geometry(mixed_far)
    g2[2] <- sf::st_sfc(sf::st_point(c(5e5, 1e200, 3)), crs = 32616)
    mixed_far <- sf::st_set_geometry(mixed_far, g2)
    got <- .zm_at(L, mixed_far, seeds = sd)
    expect_identical(got$value, ref$value)
    expect_true(all(is.na(unlist(got$value[2, setdiff(names(got$value), "cell_id")]))))
  }
})

test_that("points that mix XY and XYZ points give the frame and the summary of the same rows in XY", {
  L <- .zm_layer()
  for (sd in list(L$seeds, NULL)) {
    ref <- .zm_at(L, L$held, seeds = sd)
    for (i in c(1L, 5L)) for (dim in c("XYZ", "XYM")) {
      got <- .zm_at(L, L$held, points = .zm_row(L$asg, i, dim), seeds = sd)
      expect_identical(got$value, ref$value)
      expect_identical(got$classes, ref$classes)
      expect_identical(got$texts, ref$texts)
    }
  }
  sm <- function(a) .zm_catch(summarize_by_cell(a, "z", id_col = "cell_id",
                                                cells_sf = L$tess$cells, estimate = "pooled"))
  ref <- sm(L$asg)
  for (i in c(1L, 5L)) for (dim in c("XYZ", "XYM")) {
    got <- sm(.zm_row(L$asg, i, dim))
    expect_identical(got$value, ref$value)
    expect_identical(attr(got$value, "pooled"), attr(ref$value, "pooled"))
    expect_identical(got$classes, ref$classes)
    expect_identical(got$texts, ref$texts)
  }
})

test_that("seeds that mix XY and XYZ points give the frame of the same seeds in XY", {
  L <- .zm_layer()
  ref <- .zm_at(L, L$held)
  for (i in c(1L, 5L)) for (dim in c("XYZ", "XYM")) {
    got <- .zm_at(L, L$held, seeds = .zm_row(L$seeds, i, dim))
    expect_identical(got$value, ref$value)
    expect_identical(got$classes, ref$classes)
    expect_identical(got$texts, ref$texts)
  }
  # Seeds in another CRS, as the frame of the same seeds in XY there.
  s16 <- sf::st_transform(L$seeds, 32616)
  ref <- .zm_at(L, L$held, seeds = s16)
  for (i in c(1L, 5L)) {
    got <- .zm_at(L, L$held, seeds = .zm_row(s16, i))
    expect_identical(got$value, ref$value)
    expect_identical(got$classes, ref$classes)
  }
})

test_that("layers in one dimension take the path they took", {
  L <- .zm_layer()
  allz <- function(l) {
    xy <- sf::st_coordinates(l)[, 1:2]
    sf::st_set_geometry(l, sf::st_sfc(lapply(seq_len(nrow(xy)), function(i)
      sf::st_point(c(xy[i, ], 3))), crs = sf::st_crs(l)))
  }
  nd  <- sf::st_transform(L$held, 32616)
  ref <- .zm_at(L, nd)
  got <- .zm_at(L, allz(nd))
  expect_identical(got$value, ref$value)
  expect_identical(got$classes, ref$classes)
  got <- .zm_at(L, L$held, points = allz(L$asg))
  expect_identical(got$value, .zm_at(L, L$held)$value)
  # The helper hands a layer in one dimension back as it came.
  z <- allz(L$held)
  expect_identical(.xy_if_mixed(z), z)
  expect_identical(.xy_if_mixed(L$held), L$held)
  expect_identical(.geom_dims(.zm_row(L$held, 2)), c("XY", "XYZ"))
  expect_identical(.geom_dims(.xy_if_mixed(.zm_row(L$held, 2, "XYM"))), "XY")
})

# Row `i` of `l` replaced by `geom`.
.zm_set <- function(l, i, geom) {
  g <- sf::st_geometry(l)
  g[i] <- sf::st_sfc(geom, crs = sf::st_crs(l))
  sf::st_set_geometry(l, g)
}

# Every point of `l` given a Z of 3 at its own x and y.
.zm_all <- function(l) {
  xy <- sf::st_coordinates(l)[, 1:2]
  sf::st_set_geometry(l, sf::st_sfc(lapply(seq_len(nrow(xy)), function(i)
    sf::st_point(c(xy[i, ], 3))), crs = sf::st_crs(l)))
}

# A square of side 2h around (x0, y0), with a Z of 3 when `z`.
.zm_square <- function(x0, y0, h, z = FALSE) {
  m <- rbind(c(x0 - h, y0 - h), c(x0 + h, y0 - h), c(x0 + h, y0 + h),
             c(x0 - h, y0 + h), c(x0 - h, y0 - h))
  if (z) m <- cbind(m, 3)
  sf::st_polygon(list(m))
}

.zm_sum <- function(L, a)
  .zm_catch(summarize_by_cell(a, "z", id_col = "cell_id", cells_sf = L$tess$cells,
                              estimate = "pooled"))

.zm_same <- function(got, ref, lab, pooled = FALSE) {
  expect_identical(got$value, ref$value, label = lab)
  if (pooled)
    expect_identical(attr(got$value, "pooled"), attr(ref$value, "pooled"), label = lab)
  expect_identical(got$classes, ref$classes, label = lab)
  expect_identical(got$texts, ref$texts, label = lab)
}

test_that("an XYZ layer with a polygon among its points is read in x and y once reduced to points", {
  # coerce_to_points(, "auto") keeps the Z of a POINT Z and gives the point
  # on the surface of a POLYGON Z in XY: the point layer read mixes dimensions
  # although the layer given does not.
  L <- .zm_layer()
  xy <- sf::st_coordinates(L$held)
  for (i in c(1L, 2L)) {
    nd_z  <- .zm_set(.zm_all(L$held), i, .zm_square(xy[i, 1], xy[i, 2], 50, z = TRUE))
    nd_xy <- .zm_set(L$held, i, .zm_square(xy[i, 1], xy[i, 2], 50))
    expect_identical(.geom_dims(nd_z), "XYZ")
    for (sd in list(L$seeds, NULL)) {
      lab <- sprintf("newdata, polygon at row %d, %s", i,
                     if (is.null(sd)) "polygons" else "seeds")
      .zm_same(.zm_at(L, nd_z, seeds = sd), .zm_at(L, nd_xy, seeds = sd), lab)
    }
  }
  axy <- sf::st_coordinates(L$asg)
  for (i in c(1L, 5L)) {
    p_z  <- .zm_set(.zm_all(L$asg), i, .zm_square(axy[i, 1], axy[i, 2], 1, z = TRUE))
    p_xy <- .zm_set(L$asg, i, .zm_square(axy[i, 1], axy[i, 2], 1))
    for (sd in list(L$seeds, NULL)) {
      lab <- sprintf("points, polygon at row %d, %s", i,
                     if (is.null(sd)) "polygons" else "seeds")
      .zm_same(.zm_at(L, L$held, points = p_z, seeds = sd),
               .zm_at(L, L$held, points = p_xy, seeds = sd), lab)
    }
    .zm_same(.zm_sum(L, p_z), .zm_sum(L, p_xy),
             sprintf("pooled summary, polygon at row %d", i), pooled = TRUE)
  }
})

test_that("an XYZ line layer with one empty line is read in x and y once reduced to points", {
  # coerce_to_points(, "auto") gives an empty LINESTRING Z an empty XY point,
  # among the POINT Z points of the other lines; an empty XY line makes the
  # layer given mix dimensions already.
  L <- .zm_layer()
  axy <- sf::st_coordinates(L$asg)
  lines <- function(z) lapply(seq_len(nrow(axy)), function(i) {
    m <- rbind(axy[i, ] - c(1, 0), axy[i, ] + c(1, 0))
    sf::st_linestring(if (z) cbind(m, 3) else m)
  })
  empty <- list(XYZ = sf::st_linestring(matrix(numeric(0), 0, 3), dim = "XYZ"),
                XY = sf::st_linestring())
  for (i in c(1L, 7L)) for (e in names(empty)) {
    l_z  <- sf::st_set_geometry(L$asg, sf::st_sfc(lines(TRUE), crs = sf::st_crs(L$asg)))
    l_xy <- sf::st_set_geometry(L$asg, sf::st_sfc(lines(FALSE), crs = sf::st_crs(L$asg)))
    l_z  <- .zm_set(l_z, i, empty[[e]])
    l_xy <- .zm_set(l_xy, i, sf::st_linestring())
    if (e == "XYZ") expect_identical(.geom_dims(l_z), "XYZ")
    lab <- sprintf("empty %s line at row %d", e, i)
    .zm_same(.zm_sum(L, l_z), .zm_sum(L, l_xy), paste("pooled summary,", lab),
             pooled = TRUE)
    .zm_same(.zm_at(L, L$held, points = l_z), .zm_at(L, L$held, points = l_xy),
             paste("cell_values_at(),", lab))
  }
})

test_that("an XYZ layer with an empty XY point is read as the same rows in XY", {
  # sf::st_point() is an empty XY point, and so is a row without geometry
  # after rbind(): a layer of XYZ points with one of them mixes dimensions.
  L <- .zm_layer()
  for (i in c(1L, 4L)) {
    nd_z  <- .zm_set(.zm_all(L$held), i, sf::st_point())
    nd_xy <- .zm_set(L$held, i, sf::st_point())
    expect_identical(.geom_dims(nd_z), if (i == 1L) c("XY", "XYZ") else c("XYZ", "XY"))
    for (sd in list(L$seeds, NULL)) {
      lab <- sprintf("newdata, empty point at row %d, %s", i,
                     if (is.null(sd)) "polygons" else "seeds")
      .zm_same(.zm_at(L, nd_z, seeds = sd), .zm_at(L, nd_xy, seeds = sd), lab)
    }
  }
  for (i in c(1L, 7L)) {
    p_z  <- .zm_set(.zm_all(L$asg), i, sf::st_point())
    p_xy <- .zm_set(L$asg, i, sf::st_point())
    .zm_same(.zm_sum(L, p_z), .zm_sum(L, p_xy),
             sprintf("pooled summary, empty point at row %d", i), pooled = TRUE)
    .zm_same(.zm_at(L, L$held, points = p_z), .zm_at(L, L$held, points = p_xy),
             sprintf("cell_values_at(), empty point at row %d", i))
  }
})

test_that("seeds with one seed repeated at its own location and a mixed row give the frame of the same seeds in XY", {
  # One ID on one location is allowed; the check of shared IDs reads the
  # seeds in x and y.
  L <- .zm_layer()
  s2 <- rbind(L$seeds, L$seeds[3, ])
  ref <- .zm_at(L, L$held, seeds = s2)
  expect_s3_class(ref$value, "data.frame")
  for (i in c(1L, 5L)) for (dim in c("XYZ", "XYM")) {
    .zm_same(.zm_at(L, L$held, seeds = .zm_row(s2, i, dim)), ref,
             sprintf("repeated seed, row %d %s", i, dim))
  }
})

test_that("?cell_values_at and ?summarize_by_cell say the layers are read in x and y", {
  # The source page from the source tree, the installed page under R CMD
  # check, which runs the tests where there is no man/.
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  txt <- page("cell_values_at.Rd")
  expect_match(txt, paste0("newdata, seeds and points are read in x and y: a layer whose ",
                           "points, or the points taken from its other features, do not all ",
                           "have the same dimensions (some with a Z or M coordinate) is read ",
                           "without its Z and M"), fixed = TRUE)
  txt <- page("summarize_by_cell.Rd")
  expect_match(txt, paste0("The distances are taken in x and y: a layer whose features, or the ",
                           "points taken from them, do not all have the same dimensions (some ",
                           "with a Z or M coordinate) is read without its Z and M."), fixed = TRUE)
})
