# tests/testthat/test-fix-crs-score-geometry.R
# ---------------------------------------------------------------------------
# .crs_distance_error() scores a lon/lat layer of fewer than 40 features on
# their representative points plus the vertices of their outline, with
# points added along each edge.  What it read of the geometry depended on
# things that say nothing about the projection:
#   * a MULTIPOINT was filled between its member points as if they were
#     edges, so the score moved with the order of the members;
#   * sf::st_coordinates() refuses a layer that mixes geometry types or
#     holds GEOMETRYCOLLECTIONs, so such a layer's outline was not read
#     (a single study-area box as a GEOMETRYCOLLECTION scored NA for every
#     candidate), and across 180 degrees the wrap test of the chooser
#     stopped ensure_projected() with sf's bare message;
#   * the ring key read Z and M, so a POLYGON Z whose vertices differ in Z
#     was not filled, and Z or M on a line or polygon made every candidate
#     NA.
# The layer is now scored in x and y, read feature by feature when it mixes
# types or holds collections, and nothing is filled between points.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
box <- function(x0, y0, w, h)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
lay <- function(...) sf::st_sf(id = seq_along(list(...)), geometry = sf::st_sfc(list(...), crs = 4326))
lay_of <- function(gs) sf::st_sf(id = seq_along(gs), geometry = sf::st_sfc(gs, crs = 4326))
# A MULTIPOLYGON or MULTILINESTRING whose first part is empty, which
# sf::st_is_empty() calls empty.
empty_first <- function(s, type = "MULTIPOLYGON")
  structure(c(list(if (type == "MULTIPOLYGON") list() else matrix(numeric(0), 0, 2)), list(unclass(s))),
            class = c("XY", type, "sfg"))
points_of <- function(m)
  lay_of(lapply(seq_len(nrow(m)), function(i) sf::st_point(m[i, ])))
score <- function(x, code) .crs_distance_error(x, sf::st_crs(code))
# The log lines of ensure_projected(), without their time stamps, the
# result and its crs_choice.
projected <- function(x) {
  r <- NULL
  l <- capture_spatialkit_log(r <- quiet(ensure_projected(x)))
  list(x = r, log = sub("^[A-Z]+ \\[[^]]*\\] ", "", l), choice = attr(r, "crs_choice"))
}

mp_layers <- function() {
  set.seed(53)
  lapply(1:3, function(k) {
    m1 <- cbind(stats::runif(5, -110, -95), stats::runif(5, 30, 45))
    m2 <- cbind(stats::runif(5, -94, -85), stats::runif(5, 30, 45))
    list(a = lay(sf::st_multipoint(m1), sf::st_multipoint(m2)),
         b = lay(sf::st_multipoint(m1[5:1, ]), sf::st_multipoint(m2[c(2, 4, 1, 5, 3), ])),
         p = points_of(rbind(m1, m2)), m1 = m1, m2 = m2)
  })
}

# ---- BUG-49: MULTIPOINT ----------------------------------------------------

test_that("a MULTIPOINT layer scores as its member points, in any member order", {
  # Points were added on the segments between consecutive members of each
  # MULTIPOINT: two MULTIPOINTs of five points scored 0.75% in UTM zone 14
  # in one member order and 0.78% in another.
  for (L in mp_layers()) for (code in c(32614, 5070)) {
    ea <- score(L$a, code); eb <- score(L$b, code); ep <- score(L$p, code)
    expect_true(is.finite(ep))
    expect_lt(abs(ea - eb), 1e-12)
    expect_lt(abs(ea - ep), 1e-12)
    expect_lt(abs(eb - ep), 1e-12)
  }
  # One MULTIPOINT, as st_union() of a point layer gives.
  m1 <- mp_layers()[[1]]$m1
  for (code in c(32614, 5070)) {
    ep <- score(points_of(m1), code)
    expect_lt(abs(score(lay(sf::st_multipoint(m1)), code) - ep), 1e-12)
    expect_lt(abs(score(lay(sf::st_multipoint(m1[c(3, 1, 5, 2, 4), ])), code) - ep), 1e-12)
  }
})

test_that("ensure_projected() chooses and logs the same for a MULTIPOINT layer in any member order", {
  L  <- mp_layers()[[1]]
  pa <- projected(L$a); pb <- projected(L$b); pp <- projected(L$p)
  expect_identical(pa$log, pb$log)
  expect_identical(pa$log, pp$log)
  expect_identical(pa$choice$name, pb$choice$name)
  expect_identical(pa$choice$chosen, pb$choice$chosen)
  expect_lt(max(abs(pa$choice$distance_error - pb$choice$distance_error)), 1e-12)
  expect_lt(max(abs(pa$choice$distance_error - pp$choice$distance_error)), 1e-12)
  expect_true(sf::st_crs(pa$x) == sf::st_crs(pb$x))
  # A layer whose CRS changed with the member order: Albers in one order,
  # a Lambert azimuthal CRS in the other.
  ms <- list(cbind(c(-90.9626110740393, -76.6553423834325, -90.4804468274048),
                   c(42.7435264007139, 37.1770066637421, 41.3754164936988)),
             cbind(c(-80.1999502693829, -90.5653981710625, -90.9396888390483),
                   c(39.6546434664189, 33.2920624776002, 35.4919825264752)))
  fa <- projected(lay_of(lapply(ms, sf::st_multipoint)))
  fb <- projected(lay_of(lapply(ms, function(m) sf::st_multipoint(m[c(3, 1, 2), ]))))
  expect_identical(fa$choice$name[fa$choice$chosen], fb$choice$name[fb$choice$chosen])
  expect_true(sf::st_crs(fa$x) == sf::st_crs(fb$x))
  expect_identical(fa$log, fb$log)
})

test_that("points inside a layer of mixed types are not filled between either", {
  L <- mp_layers()[[1]]
  q <- c(-100, 37)
  mixed_a <- lay(sf::st_multipoint(L$m1), sf::st_point(q), sf::st_multipoint(L$m2))
  mixed_b <- lay(sf::st_multipoint(L$m1[5:1, ]), sf::st_point(q), sf::st_multipoint(L$m2[c(3, 5, 1, 2, 4), ]))
  expect_s3_class(sf::st_geometry(mixed_a), "sfc_GEOMETRY")
  p <- points_of(rbind(L$m1, q, L$m2))
  for (code in c(32614, 5070)) {
    ep <- score(p, code)
    expect_lt(abs(score(mixed_a, code) - ep), 1e-12)
    expect_lt(abs(score(mixed_b, code) - ep), 1e-12)
  }
  # And the same members in a GEOMETRYCOLLECTION.
  gc <- lay(sf::st_geometrycollection(list(sf::st_multipoint(L$m1[c(2, 4, 1, 3, 5), ]),
                                           sf::st_point(q))),
            sf::st_multipoint(L$m2))
  for (code in c(32614, 5070))
    expect_lt(abs(score(gc, code) - score(p, code)), 1e-12)
})

test_that("MULTIPOINT layers of 40 distinct locations or more score as before", {
  # Measured on one point per feature; casting to POINT would change this.
  set.seed(7)
  big <- lay_of(lapply(1:45, function(j)
    sf::st_multipoint(cbind(stats::runif(3, -110, -85), stats::runif(3, 30, 45)))))
  expect_identical(score(big, 5070), as.numeric("0x1.3c83cd0b0504p-7"))
  expect_identical(score(big, 32614), as.numeric("0x1.b09e66399c6p-7"))
})

# ---- BUG-123: mixed types and collections ----------------------------------

test_that("a layer of mixed types or of collections scores and chooses as its cast", {
  b1 <- box(-120, 30, 20, 17); b2 <- box(-95, 30, 20, 17); c1 <- box(-124, 25, 57, 24)
  hole <- sf::st_polygon(list(rbind(c(-120, 30), c(-80, 30), c(-80, 48), c(-120, 48), c(-120, 30)),
                              rbind(c(-110, 35), c(-100, 35), c(-100, 40), c(-110, 40), c(-110, 35))))
  ring <- lapply(0:17, function(k) box(-180 + 20 * k, -75, 4, 3))
  ring_mixed <- ring; ring_mixed[[1]] <- sf::st_multipolygon(list(ring[[1]]))
  cases <- list(
    # Two boxes, one listed as a MULTIPOLYGON: scored on one pair between
    # their interior points.
    two_boxes = list(lay(b1, sf::st_multipolygon(list(b2))),
                     lay(sf::st_multipolygon(list(b1)), sf::st_multipolygon(list(b2)))),
    # One study-area box as a GEOMETRYCOLLECTION: every candidate NA and
    # UTM zone 15 kept.
    collection = list(lay(sf::st_geometrycollection(list(c1))), lay(c1)),
    # Two boxes in one GEOMETRYCOLLECTION against one MULTIPOLYGON.
    collection2 = list(lay(sf::st_geometrycollection(list(b1, b2))),
                       lay(sf::st_multipolygon(list(b1, b2)))),
    # Across 180: ensure_projected() stopped with sf's bare message.
    fiji = list(lay(box(176, -20, 3, 4), sf::st_multipolygon(list(box(-179.5, -20, 3, 4)))),
                lay(sf::st_multipolygon(list(box(176, -20, 3, 4))),
                    sf::st_multipolygon(list(box(-179.5, -20, 3, 4))))),
    # Round the South Pole.
    polar = list(lay_of(ring_mixed), lay_of(lapply(ring, function(b) sf::st_multipolygon(list(b))))),
    # A polygon with a hole beside a MULTIPOLYGON, and a POLYGON beside a
    # MULTIPOLYGON of two parts: no edge joins two rings or two parts.
    hole = list(lay(hole, sf::st_multipolygon(list(box(-75, 30, 5, 5)))),
                lay(sf::st_multipolygon(list(hole)), sf::st_multipolygon(list(box(-75, 30, 5, 5))))),
    two_parts = list(lay(box(-120, 30, 10, 10), sf::st_multipolygon(list(box(-100, 35, 8, 8), box(-85, 30, 9, 9)))),
                     lay(sf::st_multipolygon(list(box(-120, 30, 10, 10))),
                         sf::st_multipolygon(list(box(-100, 35, 8, 8), box(-85, 30, 9, 9))))),
    # A GEOMETRYCOLLECTION whose first member is empty, which
    # sf::st_is_empty() calls empty: every candidate scored NA.
    empty_member = list(lay(sf::st_geometrycollection(list(sf::st_linestring(), c1))), lay(c1)))
  for (nm in names(cases)) {
    a <- cases[[nm]][[1]]; b <- cases[[nm]][[2]]
    pa <- NULL
    expect_no_error(pa <- projected(a))
    pb <- projected(b)
    expect_identical(pa$choice, pb$choice, label = nm)
    expect_identical(pa$log, pb$log, label = nm)
    expect_true(sf::st_crs(pa$x) == sf::st_crs(pb$x), label = nm)
    expect_true(all(is.finite(pa$choice$distance_error)), label = nm)
    for (code in c(32614, 5070, 3857))
      expect_identical(score(a, code), score(b, code), label = paste(nm, code))
  }
  # A builder that projects a boundary across 180 through ensure_projected().
  fiji <- cases$fiji[[1]]
  g <- NULL
  expect_no_error(g <- quiet(create_grid_polygons(fiji, target_cells = 12)))
  g_cast <- quiet(create_grid_polygons(cases$fiji[[2]], target_cells = 12))
  expect_identical(nrow(g), nrow(g_cast))
  expect_true(sf::st_crs(g) == sf::st_crs(g_cast))
})

test_that("a polygon beside a line is scored on the vertices of both", {
  # The box alone scores 0.067 in UTM zone 14; read with the line, the
  # layer scored on two interior points and gave 0.008.
  bx <- box(-120, 30, 45, 17)
  ln <- sf::st_linestring(rbind(c(-118, 44), c(-80, 33)))
  mixed <- lay(bx, ln)
  # The points it is scored on, built here: the two representative points,
  # the six vertices, and six points along each of the five edges (40 / 6
  # vertices), each edge filled from its end with the smaller longitude.
  edge_pts <- function(a, b, k = 6) {
    if (a[1] > b[1] || (a[1] == b[1] && a[2] > b[2])) { t <- a; a <- b; b <- t }
    f <- seq_len(k) / (k + 1)
    cbind(a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2]))
  }
  v  <- unclass(bx)[[1]]
  lv <- unclass(ln)
  rep_pts <- sf::st_coordinates(suppressWarnings(sf::st_point_on_surface(sf::st_sfc(bx, ln, crs = 4326))))
  pts <- rbind(rep_pts[, 1:2], v[1:4, ], lv,
               do.call(rbind, lapply(1:4, function(i) edge_pts(v[i, ], v[i + 1, ]))),
               edge_pts(lv[1, ], lv[2, ]))
  expect_identical(nrow(pts), 38L)
  for (code in c(32614, 5070, 3857)) {
    ref <- score(points_of(pts), code)
    expect_equal(score(mixed, code), ref, tolerance = 1e-12, label = paste(code))
    # The same features as GEOMETRYCOLLECTIONs, and in the other order.
    expect_equal(score(lay(sf::st_geometrycollection(list(bx)), sf::st_geometrycollection(list(ln))), code),
                 ref, tolerance = 1e-12)
    expect_equal(score(lay(ln, bx), code), ref, tolerance = 1e-12)
  }
})

test_that("the vertex reader gives sf::st_coordinates() on a layer of one type, and the cast's x and y otherwise", {
  g1 <- sf::st_sfc(box(-110, 30, 12, 15), box(-98, 30, 12, 15), crs = 4326)
  expect_identical(.vertex_coords(g1), sf::st_coordinates(g1))
  gm <- sf::st_sfc(box(-110, 30, 12, 15), sf::st_multipolygon(list(box(-98, 30, 12, 15))), crs = 4326)
  vm <- .vertex_coords(gm)
  expect_identical(unname(vm[, 1:2]), unname(sf::st_coordinates(sf::st_cast(gm, "MULTIPOLYGON"))[, 1:2]))
  expect_identical(colnames(vm), c("X", "Y", "L1", "L2", "L3", "L4"))
  # Empty members are left out; Z is dropped.
  gz <- sf::st_sfc(sf::st_geometrycollection(list(sf::st_point(c(1, 2, 3)), sf::st_linestring(matrix(numeric(0), 0, 3), dim = "XYZ"))),
                   sf::st_point(c(4, 5, 6)), crs = 4326)
  vz <- .vertex_coords(gz)
  expect_identical(unname(vz[, 1:2]), rbind(c(1, 2), c(4, 5)))
  # A layer that mixes dimensions as well as types is read in x and y.
  pz <- sf::st_polygon(list(cbind(rbind(c(176, -20), c(179, -20), c(179, -16), c(176, -16), c(176, -20)),
                                  c(1, 2, 3, 4, 1))), dim = "XYZ")
  gd <- sf::st_sfc(pz, box(-179.5, -20, 3, 4), sf::st_multipolygon(list(box(-170, -20, 2, 2))), crs = 4326)
  expect_identical(unname(.vertex_coords(gd)[, 1:2]),
                   unname(sf::st_coordinates(sf::st_cast(sf::st_zm(gd), "MULTIPOLYGON"))[, 1:2]))
  # A leaf whose first part is empty is read.
  ge <- sf::st_sfc(box(-110, 30, 12, 15), empty_first(box(-98, 30, 12, 15)), crs = 4326)
  expect_identical(unname(.vertex_coords(ge)[, 1:2]), unname(sf::st_coordinates(g1)[, 1:2]))
})

test_that("a feature whose first part is empty is scored on its other parts", {
  # sf::st_is_empty() calls such a feature empty, and it was left out: a
  # layer of one such box scored NA for every candidate.
  c1 <- box(-124, 25, 57, 24)
  ln <- sf::st_linestring(rbind(c(-120, 30), c(-100, 45), c(-80, 35)))
  cases <- list(mpoly = list(lay(empty_first(c1)), lay(sf::st_multipolygon(list(c1)))),
                mline = list(lay(empty_first(ln, "MULTILINESTRING")), lay(sf::st_multilinestring(list(unclass(ln))))))
  for (nm in names(cases)) {
    a <- cases[[nm]][[1]]; b <- cases[[nm]][[2]]
    expect_true(sf::st_is_empty(sf::st_geometry(a))[1L])
    for (code in c(32614, 5070, 3857)) {
      expect_true(is.finite(score(a, code)), label = paste(nm, code))
      expect_identical(score(a, code), score(b, code), label = paste(nm, code))
    }
    pa <- projected(a); pb <- projected(b)
    expect_identical(pa$choice, pb$choice, label = nm)
    expect_identical(pa$log, pb$log, label = nm)
  }
})

test_that("a layer in degrees with no CRS is read as lon/lat whatever its geometry types", {
  # sf::st_coordinates() refuses a layer of mixed types, and the refusal
  # made a small layer in degrees look planar: it was returned with no CRS
  # where the same boxes as MULTIPOLYGONs were projected.
  b1 <- box(10.123, 50.456, 0.3, 0.2); b2 <- box(10.5, 50.456, 0.3, 0.2)
  mp <- function(...) sf::st_sfc(lapply(list(...), function(b) sf::st_multipolygon(list(b))))
  cases <- list(
    mixed = list(sf::st_sfc(b1, sf::st_multipolygon(list(b2))), mp(b1, b2)),
    collection = list(sf::st_sfc(sf::st_geometrycollection(list(b1))), mp(b1)),
    # A feature whose first part is empty, which sf::st_is_empty() calls
    # empty: the layer had no coordinates left to read.
    empty_first = list(sf::st_sfc(empty_first(b1)), mp(b1)))
  for (cs in cases) {
    s <- cs[[1]]; cast <- cs[[2]]
    expect_true(.looks_like_lonlat(s)$lonlat)
    expect_true(.looks_like_lonlat(sf::st_sf(id = seq_along(s), geometry = s))$lonlat)
    a <- quiet(ensure_projected(sf::st_sf(id = seq_along(s), geometry = s)))
    b <- quiet(ensure_projected(sf::st_sf(id = seq_along(cast), geometry = cast)))
    expect_true(sf::st_crs(a) == sf::st_crs(b))
    expect_false(is.na(sf::st_crs(a)))
    expect_identical(attr(a, "crs_assumed"), attr(b, "crs_assumed"))
  }
})

test_that("the vertex-mean centre reads a layer of mixed types", {
  # s2 refuses a bowtie even after the repair, and the centre falls back to
  # the mean of the vertices, which sf::st_coordinates() refused on a layer
  # of mixed types: the centre was NA and the layer went to EPSG:3857.
  bow <- sf::st_polygon(list(rbind(c(-110, 30), c(-90, 45), c(-90, 30), c(-110, 45), c(-110, 30))))
  tri <- sf::st_polygon(list(rbind(c(-100, 32), c(-95, 32), c(-95, 40), c(-100, 32))))
  m  <- sf::st_sfc(bow, sf::st_multipolygon(list(tri)), crs = 4326)
  mc <- sf::st_cast(m, "MULTIPOLYGON")
  expect_s3_class(m, "sfc_GEOMETRY")
  ca <- .lonlat_centre(m); cb <- .lonlat_centre(mc)
  expect_true(all(is.finite(ca)))
  expect_equal(ca, cb, tolerance = 1e-12)
  a <- projected(sf::st_sf(id = 1:2, geometry = m)); b <- projected(sf::st_sf(id = 1:2, geometry = mc))
  expect_true(sf::st_crs(a$x) == sf::st_crs(b$x))
  expect_false(sf::st_crs(a$x) == sf::st_crs(3857))
  expect_identical(a$log, b$log)
})

# ---- BUG-124: Z and M ------------------------------------------------------

test_that("a layer with Z or M scores, chooses and logs as its XY copy", {
  zbox <- sf::st_polygon(list(cbind(rbind(c(-120, 30), c(-75, 30), c(-75, 47), c(-120, 47), c(-120, 30)),
                                    c(0, 10, 20, 30, 0))), dim = "XYZ")
  zline <- sf::st_linestring(rbind(c(-124, 30, 10), c(-90, 32, 20), c(-70, 44, 30)), dim = "XYZ")
  mbox <- sf::st_polygon(list(cbind(rbind(c(-124, 25), c(-67, 25), c(-67, 49), c(-124, 49), c(-124, 25)),
                                    c(1, 2, 3, 4, 1))), dim = "XYM")
  zmbox <- sf::st_polygon(list(cbind(rbind(c(-124, 25), c(-67, 25), c(-67, 49), c(-124, 49), c(-124, 25)),
                                     c(5, 6, 7, 8, 5), c(1, 2, 3, 4, 1))), dim = "XYZM")
  zmulti <- sf::st_multipoint(cbind(c(-110, -100, -90, -80), c(32, 44, 35, 41), c(1, 2, 3, 4)), dim = "XYZ")
  cases <- list(zbox = lay(zbox), zline = lay(zline), mbox = lay(mbox), zmbox = lay(zmbox),
                zmixed = lay(zbox, zline), zmulti = lay(zmulti))
  for (nm in names(cases)) {
    a <- cases[[nm]]; b <- sf::st_zm(a, drop = TRUE, what = "ZM")
    pa <- projected(a); pb <- projected(b)
    expect_identical(pa$choice, pb$choice, label = nm)
    expect_identical(pa$log, pb$log, label = nm)
    expect_true(all(is.finite(pa$choice$distance_error)), label = nm)
    for (code in c(32614, 5070))
      expect_identical(score(a, code), score(b, code), label = paste(nm, code))
  }
})

# ---- Layers the change does not reach --------------------------------------

test_that("layers of one geometry type in x and y score as before", {
  # Recorded before the change, as hexadecimal doubles, which read back
  # exactly.
  set.seed(3)
  pts <- lay_of(lapply(1:60, function(i) sf::st_point(c(stats::runif(1, -120, -75), stats::runif(1, 30, 47)))))
  hole <- sf::st_polygon(list(rbind(c(-120, 30), c(-80, 30), c(-80, 48), c(-120, 48), c(-120, 30)),
                              rbind(c(-110, 35), c(-100, 35), c(-100, 40), c(-110, 40), c(-110, 35))))
  L <- list(
    polys = lay(box(-110, 30, 12, 15), box(-98, 30, 12, 15)),
    one = lay(box(-124, 25, 57, 24)),
    hole = lay(hole),
    lines = lay(sf::st_linestring(rbind(c(-120, 30), c(-100, 45), c(-80, 35))),
                sf::st_linestring(rbind(c(-118, 44), c(-80, 33)))),
    mline = lay(sf::st_multilinestring(list(rbind(c(-120, 30), c(-100, 45)), rbind(c(-90, 32), c(-75, 44))))),
    mpoly = lay(sf::st_multipolygon(list(box(-120, 30, 10, 10), box(-95, 35, 15, 8))),
                sf::st_multipolygon(list(box(-105, 42, 6, 5)))),
    fiji = lay(box(176, -20, 3, 4), box(-179.5, -20, 3, 4)),
    pts = pts,
    pz = lay_of(lapply(1:30, function(i) sf::st_point(c(-120 + i, 30 + i / 2, i)))))
  pins <- list(
    polys = c("0x1.24c5058a6a2p-6", "0x1.3dd7a441239p-7", "0x1.aa54ab19b5e58p-2"),
    one   = c("0x1.1890edd49d908p-3", "0x1.6ce269f2b0cp-6", "0x1.1e5a7bdf06264p-1"),
    hole  = c("0x1.96f26ed8c8ep-5", "0x1.97be825f248p-7", "0x1.049511137fe3cp-1"),
    lines = c("0x1.6ef5d3124ab2p-5", "0x1.3a77c1f98b8cp-7", "0x1.9c50de140d7ecp-2"),
    mline = c("0x1.79cc996e428ap-5", "0x1.37f81fe82334p-7", "0x1.9f0b226970cccp-2"),
    mpoly = c("0x1.9c522d24c1aep-5", "0x1.3ef8ceefcc2p-7", "0x1.dadef504fc2fp-2"),
    fiji  = c("0x1.2e4c450f272f2p+1", "0x1.a4b6e5535fcb4p-2", "0x1.fa7abe28ff40ep+7"),
    pts   = c("0x1.574dbb751c0ep-5", "0x1.3a800be54d2p-7", "0x1.d26a219662278p-2"),
    pz    = c("0x1.65f401a77e14p-5", "0x1.0fa289ec73b8p-8", "0x1.a11f2b8c12da4p-2"))
  pins <- lapply(pins, as.numeric)
  for (nm in names(L))
    expect_identical(vapply(c(32614, 5070, 3857), function(code) score(L[[nm]], code), numeric(1)),
                     pins[[nm]], label = nm)
})
