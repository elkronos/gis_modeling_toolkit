# BUG-191: create_grid_polygons_cached() projected a boundary in degrees with
# .project_for_grid() and handed that copy to create_grid_polygons(), which
# then transformed it to the `crs` in `...`.  The uncached builder transforms
# the boundary once, so with `crs` given every corner of the cached grid
# differed from the uncached grid's in the last bits and no cell matched on
# geometry.  With `crs` given the builder is now handed the caller's
# boundary; the projected copy is kept for the cache key, which is unchanged.

.cg_s2 <- function(on, env = parent.frame()) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
.cg_q <- function(x) suppressWarnings(suppressMessages(x))

.cg_nc8 <- function() {
  sf::st_transform(sf::st_read(system.file("shape/nc.shp", package = "sf"),
                               quiet = TRUE)[1:8, ], 4326)
}

# The key as .cache_key() wrote it before the change, rebuilt from its parts
# on the boundary .project_for_grid() gives; `s2` is the setting hashed for a
# grid returned in degrees (BUG-185), NULL for any other grid.
.cg_old_key <- function(b, type, target_cells, ..., s2 = NULL) {
  bp <- .cg_q(spatialkit:::.project_for_grid(b, "create_grid_polygons_cached"))
  geom_hash <- digest::digest(sf::st_as_binary(spatialkit:::.boundary_union(bp)))
  dots <- list(...)
  if (length(dots) && !is.null(names(dots))) dots <- dots[order(names(dots))]
  payload <- list(geom_hash = geom_hash, crs = sf::st_crs(bp)$wkt,
                  target_cells = target_cells, args = dots,
                  version = spatialkit:::.spatialkit_version())
  if (!is.null(s2)) payload$s2 <- s2
  paste0("spatialkit_grid::", type, "::", digest::digest(payload))
}

.cg_calls <- list(
  cellsize_4326  = list(cellsize = 0.1, crs = 4326, quiet = TRUE),
  target_4326    = list(target_cells = 40, crs = 4326, quiet = TRUE),
  cellsize_32617 = list(cellsize = 10000, crs = 32617, quiet = TRUE),
  target_32617   = list(target_cells = 40, crs = 32617, quiet = TRUE))

test_that("with `crs` given, the cached grid of a lon/lat boundary is the uncached grid", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  b <- .cg_nc8()
  for (s in c(TRUE, FALSE)) for (nm in names(.cg_calls)) {
    .cg_s2(s)
    a <- .cg_calls[[nm]]
    lbl <- paste(nm, "s2", s)
    env <- new.env()
    g <- .cg_q(do.call(create_grid_polygons_cached, c(list(b), a, list(cache_env = env))))
    u <- .cg_q(ensure_stable_poly_id(do.call(create_grid_polygons, c(list(b), a))))
    expect_identical(g, u, label = lbl)
    # Every cell matches on geometry, as the help's example matches them.
    expect_identical(sf::st_as_text(sf::st_geometry(g)),
                     sf::st_as_text(sf::st_geometry(u)), label = lbl)
    # A second call is a hit on the one entry.
    expect_identical(.cg_q(do.call(create_grid_polygons_cached,
                                   c(list(b), a, list(cache_env = env)))), g, label = lbl)
    expect_length(ls(env), 1L)
  }
})

test_that("`crs` by position, by partial name and as a layer is read as create_grid_polygons() reads it", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  b <- .cg_nc8()
  utm <- sf::st_transform(b[1, ], 32617)
  for (s in c(TRUE, FALSE)) {
    .cg_s2(s)
    u <- .cg_q(ensure_stable_poly_id(create_grid_polygons(b, cellsize = 10000,
                                                          crs = 32617, quiet = TRUE)))
    by_pos <- .cg_q(create_grid_polygons_cached(b, NULL, "square", 10000, NULL, TRUE,
                                                32617, TRUE, cache_env = new.env()))
    by_part <- .cg_q(create_grid_polygons_cached(b, cellsize = 10000, cr = 32617,
                                                 quiet = TRUE, cache_env = new.env()))
    by_layer <- .cg_q(create_grid_polygons_cached(b, cellsize = 10000, crs = utm,
                                                  quiet = TRUE, cache_env = new.env()))
    expect_identical(by_pos, u, label = paste("position, s2", s))
    expect_identical(by_part, u, label = paste("partial name, s2", s))
    expect_identical(by_layer, u, label = paste("layer, s2", s))
  }
})

test_that("without `crs` the cached grid of a lon/lat boundary is still the uncached grid", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  b <- .cg_nc8()
  for (s in c(TRUE, FALSE)) {
    .cg_s2(s)
    for (a in list(list(target_cells = 40, quiet = TRUE),
                   list(target_cells = 30, type = "hex", quiet = TRUE))) {
      g <- .cg_q(do.call(create_grid_polygons_cached, c(list(b), a, list(cache_env = new.env()))))
      u <- .cg_q(ensure_stable_poly_id(do.call(create_grid_polygons, c(list(b), a))))
      expect_identical(g, u)
    }
  }
})

test_that("the cache keys are the keys they were before the change", {
  skip_if_not(file.exists(system.file("shape/nc.shp", package = "sf")))
  b <- .cg_nc8()
  calls <- c(.cg_calls, list(
    no_crs_square = list(target_cells = 40, quiet = TRUE),
    no_crs_hex    = list(target_cells = 30, type = "hex", quiet = TRUE),
    no_crs_cell   = list(cellsize = 10000, quiet = TRUE)))
  for (s in c(TRUE, FALSE)) for (nm in names(calls)) {
    .cg_s2(s)
    a <- calls[[nm]]
    tc <- a$target_cells
    type <- if (is.null(a$type)) "square" else a$type
    rest <- a[setdiff(names(a), c("target_cells", "type"))]
    in_degrees <- identical(a$crs, 4326)
    want <- do.call(.cg_old_key, c(list(b, type, tc), rest,
                                   list(s2 = if (in_degrees) s else NULL)))
    env <- new.env()
    .cg_q(do.call(create_grid_polygons_cached, c(list(b), a, list(cache_env = env))))
    expect_identical(ls(env), want, label = paste(nm, "s2", s))
  }
})
