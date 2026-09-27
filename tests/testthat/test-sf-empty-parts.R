# sf 1.0-18 changed st_is_empty(): a multi-part geometry now counts as empty
# when its FIRST part is empty, so MULTILINESTRING (EMPTY, (0 0, 10 0)) is
# "empty" there and not in earlier sf.  st_cast() consults it, and
# coerce_to_points() returned an EMPTY POINT for such a feature on current
# sf, losing its real part.  The newer rule is emulated here so the test
# bites on any installed sf.

.sf118_is_empty <- function(x) {
  g <- sf::st_geometry(x)
  vapply(unclass(g), function(item) {
    if (inherits(item, "POINT")) return(all(is.na(item)))
    if (length(item) == 0L) return(TRUE)
    if (is.list(item)) {
      f <- item[[1]]
      return(length(f) == 0L || (is.list(f) && length(f[[1]]) == 0L))
    }
    FALSE
  }, logical(1))
}

test_that("a line whose first part is empty keeps its real part under sf >= 1.0-18", {
  local_mocked_bindings(st_is_empty = .sf118_is_empty, .package = "sf")
  seg <- function(x0, y0, x1, y1) rbind(c(x0, y0), c(x1, y1))
  e2  <- matrix(numeric(0), 0, 2)

  x <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    sf::st_multilinestring(list(e2, seg(0, 0, 10, 0))),
    sf::st_multilinestring(list(e2, seg(3, 4, 3, 4))), crs = 32617))
  expect_true(all(sf::st_is_empty(x)))          # the emulation is in force
  got <- coerce_to_points(x, "auto")
  expect_equal(unname(sf::st_coordinates(got)), rbind(c(5, 0), c(3, 4)))

  ln <- function(x0, y0, d) rbind(c(x0, y0), c(x0 + d, y0), c(x0 + 2 * d, y0))
  ll <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(
    sf::st_multilinestring(list(ln(-80, 35, 0.1))),
    sf::st_multilinestring(list(e2, ln(-79, 35.5, 0.1))),
    sf::st_multilinestring(list(ln(-78, 36, 0.1))), crs = 4326))
  got <- suppressWarnings(coerce_to_points(ll, "auto"))
  expect_equal(unname(sf::st_coordinates(got)[2, ]), c(-78.9, 35.5),
               tolerance = 1e-4)
})
