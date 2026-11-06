# BUG-195: ensure_stable_poly_id(make_valid = TRUE) on a layer in degrees
# with s2 on repaired the rows s2 or GEOS calls invalid with s2's rebuild,
# which does not split crossing edges.  A bowtie, a figure eight and a
# MULTIPOLYGON whose parts overlap came back invalid under s2 and GEOS, where
# with s2 off GEOS returns them valid; st_area() stopped on them and their
# projected copy stopped GEOS's overlay.  Such a row now gets GEOS's repair
# of its coordinates, as with s2 off.

.ce_s2 <- function(on, env = parent.frame()) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
.ce_q <- function(x) suppressWarnings(suppressMessages(x))
.ce_sq <- function(x0, y0, w) rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w),
                                    c(x0, y0 + w), c(x0, y0))
.ce_ok <- sf::st_polygon(list(.ce_sq(-80, 30, 0.3)))
.ce_cases <- list(
  bowtie = sf::st_polygon(list(rbind(c(-92, 39), c(-91.5, 39.5), c(-91.5, 39),
                                     c(-92, 39.5), c(-92, 39)))),
  overlapping_parts = sf::st_multipolygon(list(list(.ce_sq(-90, 39, 0.5)),
                                               list(.ce_sq(-89.7, 39.2, 0.5)))),
  figure_eight = sf::st_polygon(list(rbind(c(-88, 39), c(-87.5, 39), c(-88, 39.5),
                                           c(-87.5, 39.5), c(-88, 39)))))
.ce_layer <- function(bad) sf::st_sf(v = 1:2, geometry = sf::st_sfc(bad, .ce_ok, crs = 4326))

# The call under one setting: the rows in input order and their IDs.
.ce_run <- function(x, s2) local({
  .ce_s2(s2)
  o <- .ce_q(ensure_stable_poly_id(x))
  list(g = sf::st_geometry(o)[order(o$v)], ids = o$poly_id[order(o$v)])
})
.ce_s2_valid <- function(g) local({ .ce_s2(TRUE); sf::st_is_valid(g) })
.ce_geos_valid <- function(g) sf::st_is_valid(sf::st_set_crs(g, NA))
.ce_s2_area <- function(g) local({ .ce_s2(TRUE); as.numeric(sf::st_area(g)) })

test_that("a row with crossing edges comes back valid under s2 and GEOS with s2 on", {
  for (nm in names(.ce_cases)) {
    x <- .ce_layer(.ce_cases[[nm]])
    # The input row is invalid under both engines.
    expect_identical(.ce_s2_valid(sf::st_geometry(x)), c(FALSE, TRUE), info = nm)
    expect_identical(.ce_geos_valid(sf::st_geometry(x)), c(FALSE, TRUE), info = nm)
    on  <- .ce_run(x, TRUE)
    off <- .ce_run(x, FALSE)
    for (r in list(on, off)) {
      expect_identical(.ce_s2_valid(r$g[1]), TRUE, info = nm)
      expect_identical(.ce_geos_valid(r$g[1]), TRUE, info = nm)
    }
    # The row s2 on returns is the one s2 off returns: the same area.
    expect_equal(.ce_s2_area(on$g[1]), .ce_s2_area(off$g[1]), tolerance = 1e-6, info = nm)
    expect_identical(on$g[[1]], off$g[[1]], info = nm)
    # Its projected copy is cut by its bounding box without error.
    p <- sf::st_transform(on$g[1], 32615)
    cut <- sf::st_intersection(p, sf::st_as_sfc(sf::st_bbox(p)))
    expect_equal(as.numeric(sf::st_area(cut)), as.numeric(sf::st_area(p)),
                 tolerance = 1e-9, info = nm)
    # The valid square comes back as it arrived.
    expect_identical(on$g[[2]], .ce_ok, info = nm)
    # The IDs are the same under both settings and those given before the change.
    expect_identical(on$ids, off$ids, info = nm)
    expect_identical(on$ids, c(1L, 2L), info = nm)
  }
})

test_that("the repair helper gives GEOS's repair of the coordinates to rows s2's rebuild leaves invalid", {
  .ce_s2(TRUE)
  g <- sf::st_sfc(c(.ce_cases, list(.ce_ok)), crs = 4326)
  h <- spatialkit:::.make_valid_rows(g)
  geos <- sf::st_set_crs(suppressWarnings(sf::st_make_valid(sf::st_set_crs(g, NA))), 4326)
  for (i in 1:3) expect_identical(h[[i]], geos[[i]], info = names(.ce_cases)[i])
  expect_identical(h[[4]], .ce_ok)
  expect_identical(sf::st_crs(h), sf::st_crs(g))
  x <- sf::st_sf(v = 1:4, geometry = g)
  expect_identical(sf::st_geometry(spatialkit:::.make_valid_rows(x)), h)
})

# With s2 on, such a row used to be keyed from the sort copy's own repair,
# whose second pass drops the overlap of overlapping parts, so its ID could
# differ from the one s2 off gives.  It is now keyed on GEOS's repair under
# both settings.  A 1 degree square and a 0.6 degree square that overlap,
# as one MULTIPOLYGON, beside a row of 111 small squares whose longitudes
# fall on both sides of its centroid: before, the row took ID 52 with s2 on
# and 55 with s2 off.
test_that("a row with overlapping parts takes the same ID with s2 on and off", {
  ov <- sf::st_multipolygon(list(list(.ce_sq(-90, 39, 1)),
                                 list(rbind(c(-89.2, 39.3), c(-88.6, 39.3), c(-88.6, 39.9),
                                            c(-89.2, 39.9), c(-89.2, 39.3)))))
  thin <- lapply(seq(-89.9, -88.8, by = 0.01),
                 function(x0) sf::st_polygon(list(.ce_sq(x0, 45, 0.001))))
  x <- sf::st_sf(v = seq_len(length(thin) + 1L),
                 geometry = sf::st_sfc(c(list(ov), thin), crs = 4326))
  on  <- .ce_run(x, TRUE)
  off <- .ce_run(x, FALSE)
  expect_identical(on$ids, off$ids)
  expect_identical(on$ids[1], 55L)
  expect_identical(on$g[[1]], off$g[[1]])
  expect_equal(.ce_s2_area(on$g[1]), .ce_s2_area(off$g[1]), tolerance = 1e-12)
})

# A polygon valid on the plane whose edges cross as great circles: s2 calls
# it invalid, its rebuild leaves it so, and GEOS's repair leaves the crossing
# in place.  The repair says so, and the row is returned as GEOS repaired it,
# which is the row s2 off returns.
test_that("a row neither repair makes valid under both checks is returned as GEOS repaired it, with a warning", {
  arc <- sf::st_polygon(list(rbind(c(-80, 70), c(80, 70), c(80, 71), c(5, 72),
                                   c(-5, 72), c(-80, 71), c(-80, 70))))
  sq <- sf::st_polygon(list(.ce_sq(100, 10, 1)))
  g <- sf::st_sfc(arc, sq, crs = 4326)
  expect_identical(.ce_s2_valid(g[1]), FALSE)
  expect_identical(.ce_geos_valid(g[1]), TRUE)
  geos <- suppressWarnings(sf::st_make_valid(sf::st_set_crs(g[1], NA)))[[1]]
  .ce_s2(TRUE)
  msg <- "could not make 1 row\\(s\\) valid under both s2 and GEOS \\(rows 1\\)"
  w <- character()
  h <- withCallingHandlers(spatialkit:::.make_valid_rows(g),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  expect_length(w, 1L)
  expect_match(w, msg)
  expect_identical(h[[1]], geos)
  expect_identical(h[[2]], sq)
  # Through ensure_stable_poly_id(): the repair warns once (the sort copy's
  # refusal warns as well), and the rows and IDs are those s2 off gives.
  x <- sf::st_sf(v = 1:2, geometry = g)
  w <- character()
  o <- withCallingHandlers(ensure_stable_poly_id(x),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  expect_identical(sum(grepl(msg, w)), 1L)
  o <- o[order(o$v), ]
  expect_identical(sf::st_geometry(o)[[1]], geos)
  expect_identical(sf::st_geometry(o)[[2]], sq)
  off <- .ce_run(x, FALSE)
  expect_identical(o$poly_id, off$ids)
  expect_identical(sf::st_geometry(o)[[1]], off$g[[1]])
})
