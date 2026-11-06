# tests/testthat/test-fix-mixed-zm-summaries.R
# ---------------------------------------------------------------------------
# sf reads a layer in the dimension of its first geometry: st_transform() of
# a layer with an XYZ point first among XY points stops with "OGR error",
# st_coordinates() lays its values out of step or stops on its "dimnames",
# and GEOS reads an XYZ ring among XY rings as XY pairs.  The variogram
# design effect of summarize_by_cell() and kriging_adequacy() read their
# points, and kriging_adequacy() its cells, in x and y: Z and M are dropped
# from the copy read when a layer mixes dimensions, with no warning, since
# neither function returns that copy.  Each mixed layer is held to the same
# rows in XY, to the bit, warnings included.
#
# The layer: 400 points in EPSG:32617 (seed 7), 16 grid cells over them, and
# the variogram estimate_sac_range() fits to them.
# ---------------------------------------------------------------------------

.zs_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value,
       classes = vapply(got, function(w) class(w)[1L], character(1)),
       texts = vapply(got, conditionMessage, character(1)))
}

.zs_store <- new.env(parent = emptyenv())
.zs_layer <- function() {
  if (is.null(.zs_store$layer)) {
    d <- local({
      withr::local_seed(7, .rng_kind = "Mersenne-Twister",
                        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
      n <- 400
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      cx <- stats::runif(30, 0, 1e5); cy <- stats::runif(30, 0, 1e5)
      z <- rowSums(sapply(1:30, function(k)
        exp(-((x - cx[k])^2 + (y - cy[k])^2) / (2 * 8e3^2)))) + stats::rnorm(n, 0, 0.1)
      data.frame(x = 5e5 + x, y = 5e6 + y, z = z)
    })
    pts <- sf::st_as_sf(d, coords = c("x", "y"), crs = 32617)
    cells <- sf::st_sf(cell_id = 1:16,
                       geometry = sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(pts)), n = c(4, 4)))
    asg <- assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id", quiet = TRUE)
    sac <- suppressWarnings(estimate_sac_range(pts, "z", cutoff = 1))
    sac_bare <- sac
    attr(sac_bare, "crs") <- NULL
    .zs_store$layer <- list(asg = asg, cells = cells, sac = sac, sac_bare = sac_bare)
  }
  .zs_store$layer
}

# Row `i` of `l` given a third coordinate, Z or M, at its own x and y.
.zs_row <- function(l, i, dim = "XYZ") {
  g  <- sf::st_geometry(l)
  co <- unclass(g[[i]])
  g[[i]] <- switch(class(g[[i]])[2L],
    POINT      = sf::st_point(c(co, 3), dim = dim),
    MULTIPOINT = sf::st_multipoint(cbind(co, 3), dim = dim),
    POLYGON    = sf::st_polygon(lapply(co, function(r) cbind(r, 3)), dim = dim))
  sf::st_set_geometry(l, sf::st_sfc(as.list(g), crs = sf::st_crs(l), recompute_bbox = TRUE))
}
.zs_all <- function(l) {
  for (i in seq_len(nrow(l))) l <- .zs_row(l, i)
  l
}

.zs_deff <- function(a, sac)
  .zs_catch(summarize_by_cell(a, "z", id_col = "cell_id", estimate = "mean",
                              deff = "variogram", sac = sac))
.zs_krige <- function(a, cells, sac = NULL)
  .zs_catch(kriging_adequacy(a, "z", cells, id_col = "cell_id", sac = sac, seed = 1))

test_that("the premise: a mixed point layer stops sf::st_transform() or garbles its coordinates", {
  L <- .zs_layer()
  expect_error(sf::st_transform(.zs_row(L$asg, 1L), 4326), "OGR error")
  expect_error(suppressWarnings(sf::st_coordinates(.zs_row(L$asg, 5L))), "dimnames")
  expect_true(is.finite(as.numeric(L$sac)))
})

test_that("summarize_by_cell(deff = \"variogram\") reads mixed points as the same points in XY", {
  L <- .zs_layer()
  for (sac in list(L$sac, L$sac_bare, NULL)) {
    ref <- .zs_deff(L$asg, sac)
    for (dim in c("XYZ", "XYM")) for (i in c(1L, 5L)) {
      got <- .zs_deff(.zs_row(L$asg, i, dim), sac)
      expect_identical(got, ref, info = sprintf("row %d %s, sac %s", i, dim,
                                                if (is.null(sac)) "NULL" else "given"))
    }
  }
})

test_that("summarize_by_cell(deff = \"variogram\") reads mixed polygon features with no warning about Z", {
  L <- .zs_layer()
  polys <- sf::st_set_geometry(L$asg, sf::st_buffer(sf::st_geometry(L$asg), 10))
  ref <- .zs_deff(polys, L$sac)
  expect_false("spatialkit_zm_dropped" %in% ref$classes)
  for (i in c(1L, 5L)) {
    got <- .zs_deff(.zs_row(polys, i), L$sac)
    expect_false("spatialkit_zm_dropped" %in% got$classes, info = sprintf("row %d", i))
    expect_identical(got, ref, info = sprintf("row %d", i))
  }
})

test_that("kriging_adequacy() reads mixed points as the same points in XY", {
  L <- .zs_layer()
  for (sac in list(L$sac, NULL)) {
    ref <- .zs_krige(L$asg, L$cells, sac)
    for (dim in c("XYZ", "XYM")) for (i in c(1L, 5L)) {
      got <- .zs_krige(.zs_row(L$asg, i, dim), L$cells, sac)
      expect_identical(got, ref, info = sprintf("row %d %s, sac %s", i, dim,
                                                if (is.null(sac)) "NULL" else "given"))
    }
  }
})

test_that("kriging_adequacy() reads mixed multipoints with no warning about Z", {
  L <- .zs_layer()
  mp <- sf::st_cast(L$asg, "MULTIPOINT")
  ref <- .zs_krige(mp, L$cells, L$sac)
  for (i in c(1L, 5L)) {
    got <- .zs_krige(.zs_row(mp, i), L$cells, L$sac)
    expect_false("spatialkit_zm_dropped" %in% got$classes, info = sprintf("row %d", i))
    expect_identical(got, ref, info = sprintf("row %d", i))
  }
})

test_that("kriging_adequacy() reads cells_sf that mixes XY and XYZ as the same cells in XY", {
  L <- .zs_layer()
  for (cells in list(L$cells, sf::st_transform(L$cells, 4326))) {
    ref <- .zs_krige(L$asg, cells)
    for (dim in c("XYZ", "XYM")) for (i in c(1L, 5L)) {
      got <- .zs_krige(L$asg, .zs_row(cells, i, dim))
      expect_identical(got, ref, info = sprintf("cells in EPSG:%s, row %d %s",
                                                sf::st_crs(cells)$epsg, i, dim))
    }
  }
})

test_that("layers in one dimension take the path they took", {
  L <- .zs_layer()
  # The flattening hands a layer in one dimension back as it came.
  expect_identical(.xy_if_mixed(.zs_all(L$asg)), .zs_all(L$asg))
  expect_identical(.xy_if_mixed(.zs_all(L$cells)), .zs_all(L$cells))
  # Points XYZ throughout for the design effect, and cells XYZ throughout,
  # give what the same layers in XY give.  (Points XYZ throughout are handed
  # to gstat with their Z, as before.)
  expect_identical(.zs_deff(.zs_all(L$asg), L$sac), .zs_deff(L$asg, L$sac))
  expect_identical(.zs_krige(L$asg, .zs_all(L$cells), L$sac),
                   .zs_krige(L$asg, L$cells, L$sac))
})

# Rows `rows` of the XYZ layer `l` made MULTIPOINT Z (the point and one 1 m
# away) or POLYGON Z (a 5 m triangle on the point), all at the same height.
.zs_retype <- function(l, rows, type) {
  lst <- as.list(sf::st_geometry(l))
  for (j in rows) {
    co <- unclass(lst[[j]])
    lst[[j]] <- switch(type,
      MULTIPOINT = sf::st_multipoint(rbind(co, co + c(1, 1, 0)), dim = "XYZ"),
      POLYGON    = sf::st_polygon(list(rbind(co, co + c(5, 0, 0), co + c(5, 5, 0), co)),
                                  dim = "XYZ"))
  }
  sf::st_set_geometry(l, sf::st_sfc(lst, crs = sf::st_crs(l)))
}

test_that("a layer XYZ throughout whose multipoints or polygons give XY points is read in x and y", {
  L <- .zs_layer()
  allz <- .zs_all(L$asg)
  mpz <- .zs_retype(allz, c(2L, 7L, 11L), "MULTIPOINT")
  pgz <- .zs_retype(allz, c(2L, 7L, 11L), "POLYGON")
  # The layers are XYZ throughout, and the points coerce_to_points() takes
  # from them are not: those of the multipoints and polygons are XY.
  for (l in list(mpz, pgz)) {
    expect_identical(unique(.geom_dims(l)), "XYZ")
    expect_setequal(.geom_dims(suppressWarnings(coerce_to_points(l, "auto"))),
                    c("XYZ", "XY"))
  }
  for (sac in list(L$sac, L$sac_bare)) {
    for (l in list(MULTIPOINT = mpz, POLYGON = pgz)) {
      ref <- .zs_deff(sf::st_zm(l), sac)
      got <- .zs_deff(l, sac)
      expect_false("spatialkit_zm_dropped" %in% got$classes)
      expect_identical(got, ref)
    }
  }
  ref <- .zs_krige(sf::st_zm(mpz), L$cells, L$sac)
  got <- .zs_krige(mpz, L$cells, L$sac)
  expect_false("spatialkit_zm_dropped" %in% got$classes)
  expect_identical(got, ref)
})

test_that("kriging_adequacy() kriges points XYZ throughout in x and y, as their XY copy", {
  # gstat::krige() kriged these points in three dimensions with a variogram
  # fitted in x and y, so the heights moved the result.  The points are now
  # read in x and y on every layer, Z and M dropped, as the variogram is
  # fitted, so the result is the XY one, warnings included.
  L <- .zs_layer()
  h <- local({
    withr::local_seed(11, .rng_kind = "Mersenne-Twister",
                      .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    stats::runif(nrow(L$asg), 0, 1000)
  })
  lst <- as.list(sf::st_geometry(L$asg))
  for (j in seq_along(lst)) lst[[j]] <- sf::st_point(c(unclass(lst[[j]]), h[j]), dim = "XYZ")
  ptsz <- sf::st_set_geometry(L$asg, sf::st_sfc(lst, crs = sf::st_crs(L$asg)))
  expect_identical(unique(.geom_dims(ptsz)), "XYZ")
  got <- .zs_krige(ptsz, L$cells, L$sac)
  ref <- .zs_krige(L$asg, L$cells, L$sac)
  expect_true(any(is.finite(ref$value$kr_var)))
  expect_false("spatialkit_zm_dropped" %in% got$classes)
  expect_identical(got, ref)
})

test_that("?kriging_adequacy and ?summarize_by_cell say mixed layers are read in x and y", {
  # The source page from the source tree, the installed page under R CMD
  # check, which runs the tests where there is no man/.
  page <- function(name) {
    file <- test_path("..", "..", "man", name)
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[name]]
    expect_s3_class(rd, "Rd")
    gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  }
  txt <- page("kriging_adequacy.Rd")
  # The points of every layer, mixed or not, are read in x and y.
  expect_match(txt, paste0("Any Z or M dimension is dropped, on every layer, because ",
                           "gstat::krige() uses every coordinate dimension of its locations and ",
                           "the variogram is fitted in x and y"), fixed = TRUE)
  expect_match(txt, paste0("Cells that do not all have the same dimensions are read in x and ",
                           "y, without their Z and M."), fixed = TRUE)
  txt <- page("summarize_by_cell.Rd")
  expect_match(txt, paste0("The distances are taken in x and y: a layer whose features, or the ",
                           "points taken from them, do not all have the same dimensions (some ",
                           "with a Z or M coordinate) is read without its Z and M. Exponential"),
               fixed = TRUE)
})
