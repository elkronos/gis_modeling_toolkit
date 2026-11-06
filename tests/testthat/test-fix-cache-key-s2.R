# BUG-185: create_grid_polygons_cached() keyed a grid returned in degrees
# without the session's sf::sf_use_s2() setting.  Such a grid is clipped by
# s2 or by GEOS (when laid in degrees) and its cells are repaired by
# ensure_stable_poly_id() with s2's rebuild or GEOS's, so create_grid_polygons()
# gives one grid per setting, and the cache handed back the grid of whichever
# setting filled it first.  The setting now goes into the key of a grid
# returned in degrees, and into no other key.

local_s2 <- function(on, env = parent.frame()) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
quiet <- function(x) suppressWarnings(suppressMessages(x))

nc8 <- function() {
  sf::st_transform(sf::st_read(system.file("shape/nc.shp", package = "sf"),
                               quiet = TRUE)[1:8, ], 4326)
}

# The grid create_grid_polygons() gives in the session, numbered as the cached
# builder numbers it.
uncached <- function(b, args) {
  quiet(ensure_stable_poly_id(do.call(create_grid_polygons, c(list(b), args))))
}
cached <- function(b, args, env) {
  quiet(do.call(create_grid_polygons_cached, c(list(b), args, list(cache_env = env))))
}

# The key as .cache_key() wrote it before the change, rebuilt here from its
# parts, so a key that should be unchanged is checked against the old formula
# and not against the new code.
old_key <- function(bnd, type, target_cells, ...) {
  geom_hash <- digest::digest(sf::st_as_binary(spatialkit:::.boundary_union(bnd)))
  dots <- list(...)
  if (length(dots) && !is.null(names(dots))) dots <- dots[order(names(dots))]
  paste0("spatialkit_grid::", type, "::",
         digest::digest(list(geom_hash = geom_hash, crs = sf::st_crs(bnd)$wkt,
                             target_cells = target_cells, args = dots,
                             version = spatialkit:::.spatialkit_version())))
}

test_that("a grid returned in degrees is cached once per sf_use_s2() setting, in either order", {
  # A projected boundary: create_grid_polygons_cached() hands it to
  # create_grid_polygons() as it is, so the two builders transform it alike.
  b <- sf::st_transform(nc8(), 32617)
  calls <- list(
    square_deg = list(cellsize = 0.1, crs = 4326, quiet = TRUE),
    hex_deg    = list(type = "hex", cellsize = 0.1, crs = 4326, quiet = TRUE),
    target_deg = list(target_cells = 40, crs = 4326, quiet = TRUE))
  for (nm in names(calls)) {
    a <- calls[[nm]]
    local_s2(TRUE);  u_on  <- uncached(b, a)
    local_s2(FALSE); u_off <- uncached(b, a)
    # The two settings give two grids, so one cached grid cannot serve both.
    expect_false(identical(u_on, u_off), label = nm)
    for (first in c(TRUE, FALSE)) {
      env <- new.env()
      local_s2(first)
      g1 <- cached(b, a, env)
      local_s2(!first)
      g2 <- cached(b, a, env)
      expect_identical(g1, if (first) u_on else u_off, label = paste(nm, first))
      expect_identical(g2, if (first) u_off else u_on, label = paste(nm, !first))
      expect_length(ls(env), 2L)
      # A second call under each setting is a hit on its own entry.
      expect_identical(cached(b, a, env), g2)
      local_s2(first)
      expect_identical(cached(b, a, env), g1)
      expect_length(ls(env), 2L)
    }
  }
})

test_that("a lon/lat boundary laid in degrees gets the grid of the session's setting", {
  b <- nc8()
  a <- list(cellsize = 0.1, crs = 4326, quiet = TRUE)
  # create_grid_polygons_cached() hands create_grid_polygons() the boundary
  # after .project_for_grid(), so that is the boundary its uncached grid is
  # built on.
  bp <- quiet(spatialkit:::.project_for_grid(b, "create_grid_polygons_cached"))
  local_s2(TRUE);  u_on  <- uncached(bp, a)
  local_s2(FALSE); u_off <- uncached(bp, a)
  expect_false(identical(u_on, u_off))
  for (first in c(TRUE, FALSE)) {
    env <- new.env()
    local_s2(first);  g1 <- cached(b, a, env)
    local_s2(!first); g2 <- cached(b, a, env)
    expect_identical(g1, if (first) u_on else u_off)
    expect_identical(g2, if (first) u_off else u_on)
    expect_length(ls(env), 2L)
  }
})

test_that("a grid laid and returned in a projected CRS keeps its key under either setting", {
  ll <- nc8()
  pr <- sf::st_transform(ll, 32617)
  cases <- list(
    list(b = pr, args = list(cellsize = 10000, quiet = TRUE)),
    list(b = pr, args = list(target_cells = 30, type = "hex", quiet = TRUE)),
    list(b = ll, args = list(target_cells = 30, quiet = TRUE)),
    list(b = ll, args = list(cellsize = 10000, crs = 32617, quiet = TRUE)),
    list(b = ll, args = list(n = 5, crs = sf::st_crs(32617), quiet = TRUE)))
  for (i in seq_along(cases)) {
    b <- cases[[i]]$b
    a <- cases[[i]]$args
    bp <- quiet(spatialkit:::.project_for_grid(b, "create_grid_polygons_cached"))
    tc <- a$target_cells
    type <- if (is.null(a$type)) "square" else a$type
    rest <- a[setdiff(names(a), c("target_cells", "type"))]
    want <- do.call(old_key, c(list(bp, type, tc), rest))
    env <- new.env()
    for (on in c(TRUE, FALSE)) {
      local_s2(on)
      g <- cached(b, a, env)
      expect_identical(ls(env), want, label = paste("case", i, "s2", on))
    }
  }
})

test_that("an argument passed by position or as a layer is read as create_grid_polygons() reads it", {
  b <- sf::st_transform(nc8(), 32617)
  layer4326 <- sf::st_transform(b[1, ], 4326)
  env <- new.env()
  local_s2(TRUE)
  # `cellsize` by position, `crs` as a layer.
  g_on <- quiet(create_grid_polygons_cached(b, NULL, "square", 0.1, crs = layer4326,
                                            quiet = TRUE, cache_env = env))
  local_s2(FALSE)
  g_off <- quiet(create_grid_polygons_cached(b, NULL, "square", 0.1, crs = layer4326,
                                             quiet = TRUE, cache_env = env))
  expect_length(ls(env), 2L)
  expect_identical(g_off, uncached(b, list(cellsize = 0.1, crs = 4326, quiet = TRUE)))
  local_s2(TRUE)
  expect_identical(g_on, uncached(b, list(cellsize = 0.1, crs = 4326, quiet = TRUE)))
})

test_that("`crs` passed by position or by partial name is cached once per setting, in either order", {
  # create_grid_polygons() matches `crs` by position (after `cellsize`, `n`
  # and `clip`) and by partial name (`cr =`), so the key reads it the same
  # way: a key that looked only for a named `crs` held one entry for both
  # settings and handed the second setting the first one's grid.
  b <- sf::st_transform(nc8(), 32617)
  want_args <- list(cellsize = 0.1, crs = 4326, quiet = TRUE)
  local_s2(TRUE);  u_on  <- uncached(b, want_args)
  local_s2(FALSE); u_off <- uncached(b, want_args)
  expect_false(identical(u_on, u_off))
  calls <- list(
    by_position = function(env)
      create_grid_polygons_cached(b, NULL, "square", 0.1, NULL, TRUE, 4326, TRUE,
                                  cache_env = env),
    by_partial_name = function(env)
      create_grid_polygons_cached(b, cellsize = 0.1, cr = 4326, quiet = TRUE,
                                  cache_env = env))
  for (nm in names(calls)) {
    for (first in c(TRUE, FALSE)) {
      env <- new.env()
      local_s2(first);  g1 <- quiet(calls[[nm]](env))
      local_s2(!first); g2 <- quiet(calls[[nm]](env))
      expect_length(ls(env), 2L)
      expect_identical(g1, if (first) u_on else u_off, label = paste(nm, first))
      expect_identical(g2, if (first) u_off else u_on, label = paste(nm, !first))
    }
  }
})

test_that("?create_grid_polygons_cached says a grid in degrees is cached once per setting", {
  file <- test_path("..", "..", "man", "create_grid_polygons_cached.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["create_grid_polygons_cached.Rd"]]
  expect_s3_class(rd, "Rd")
  page <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(page, "is cached once per sf::sf_use_s2() setting", fixed = TRUE)
})
