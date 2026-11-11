# tests/testthat/test-fix-stable-ids-warnings.R
# ---------------------------------------------------------------------------
# ensure_stable_poly_id(method = "surface_point") takes the points on the
# surface of the copy in EPSG:4326 its sort key is measured on.  sf warned
# "st_point_on_surface may not give correct results for longitude/latitude
# data" on every call, in its own words, about a key the help says is a
# sort order in degrees.  The text is muffled at that one call
# (.key_surface_points()), as .pooled_xy() muffles it, and the points are
# the same to the bit, so no ID moves.
# ---------------------------------------------------------------------------

.sw_text <- "st_point_on_surface may not give correct results"

# Six 1 km squares, 3 columns by 2 rows, in UTM 17N, with a column `name`.
.sw_layer <- function() {
  sq <- function(x0, y0, s = 1000)
    sf::st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s),
                              c(x0, y0 + s), c(x0, y0))))
  cells <- list()
  for (j in 0:1) for (i in 0:2)
    cells[[length(cells) + 1L]] <- sq(6e5 + 1000 * i, 3.99e6 + 1000 * j)
  sf::st_sf(name = letters[1:6], geometry = sf::st_sfc(cells, crs = 32617))
}

# The same call with sf's lon/lat text muffled by hand, and nothing else.
.sw_muffled <- function(expr) {
  withCallingHandlers(expr, warning = function(w)
    if (grepl(.sw_text, conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning"))
}

.sw_cases <- function() {
  lay <- .sw_layer()
  list(projected = lay,
       lonlat    = sf::st_transform(lay, 4326),
       no_crs    = sf::st_set_crs(lay, NA),
       one_row   = lay[4, ],
       bare_sfc  = sf::st_geometry(lay))
}

test_that("surface_point raises no warning on projected, lon/lat, CRS-less and one-row layers", {
  for (nm in names(.sw_cases())) {
    x <- .sw_cases()[[nm]]
    expect_no_warning(out <- ensure_stable_poly_id(x, method = "surface_point"))
    ref <- .sw_muffled(ensure_stable_poly_id(x, method = "surface_point"))
    expect_identical(out, ref, info = nm)
  }
})

test_that("surface_point gives the IDs it gave: the order of the squares is unchanged", {
  out <- ensure_stable_poly_id(.sw_layer(), method = "surface_point")
  # Columns west to east, south before north within a column, as the
  # centroid key gives on the same squares.
  expect_identical(out$name, c("a", "d", "b", "e", "c", "f"))
  expect_identical(out$poly_id, 1:6)
  cen <- ensure_stable_poly_id(.sw_layer(), method = "centroid")
  expect_identical(out$name, cen$name)
})

test_that("the points on the surface are sf's own, bit for bit", {
  g <- sf::st_geometry(sf::st_transform(.sw_layer(), 4326))
  expect_identical(
    .key_surface_points(g),
    suppressWarnings(sf::st_point_on_surface(.drop_empty_parts(g))))
})

test_that("other warnings of the call still reach the caller", {
  # A row that is not a polygon is dropped with the package's own warning.
  lay <- .sw_layer()
  mixed <- sf::st_sf(name = c(lay$name, "pt"),
                     geometry = c(sf::st_geometry(lay),
                                  sf::st_sfc(sf::st_point(c(6e5, 3.99e6)), crs = 32617)))
  expect_warning(ensure_stable_poly_id(mixed, method = "surface_point"),
                 "dropping 1 non-polygon row")
})
