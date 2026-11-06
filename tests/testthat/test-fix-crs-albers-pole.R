# tests/testthat/test-fix-crs-albers-pole.R
# ---------------------------------------------------------------------------
# BUG-204.  On the comparison path and the purpose = "area" path the CRS
# chooser builds an Albers candidate with standard parallels ymin + span / 6
# and ymax - span / 6, and on a box with no height it falls back to lat - 5
# and lat + 5 around the centre on the sphere.  A line along one parallel
# that spans 180 degrees of longitude runs over the pole on the sphere, so
# its centre is the pole and the fallback parallels are 85 and 95.  The cone
# test passed, sf::st_crs() stopped with "invalid crs", and ensure_projected()
# stopped with it, with either purpose.  The same held for any centre within
# 5 degrees of a pole on a box with no height.
#
# The Albers candidate is now built only when both standard parallels lie
# between -90 and 90; otherwise it is left out of the comparison, as a
# degenerate cone is.  Layers whose parallels lie in range get the choice
# they got.
# ---------------------------------------------------------------------------

.cap_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.cap_line <- function(m)
  sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_linestring(m), crs = 4326))
# The lines of the item: centre on the pole on the sphere.
.cap_pole_lines <- function() list(
  line20  = rbind(c(0, 20), c(180, 20)),
  line45  = rbind(c(0, 45), c(180, 45)),
  line60  = rbind(c(0, 60), c(180, 60)),
  edge70  = rbind(c(-10, 70), c(170, 70)),
  south30 = rbind(c(0, -30), c(180, -30)))
# Centre within 5 degrees of a pole but not on it.
.cap_near_lines <- function() list(
  north80 = rbind(c(0, 80), c(160, 80)),    # centre (80, 88.2)
  south84 = rbind(c(0, -84), c(120, -84)))  # centre (60, -87.0)
# The rule before the fix: the cone test alone.
.cap_cone_only <- function(lat1, lat2) {
  cone_n <- (sin(lat1 * pi / 180) + sin(lat2 * pi / 180)) / 2
  is.finite(cone_n) && abs(cone_n) >= 0.05
}

test_that(".albers_parallels_ok() keeps the cone test and refuses parallels beyond a pole", {
  expect_true(.albers_parallels_ok(35, 45))
  expect_true(.albers_parallels_ok(80, 90))
  expect_true(.albers_parallels_ok(-90, -80))
  expect_true(.albers_parallels_ok(79.6, 89.6))
  # Parallels close together, as on a box of small but nonzero height.
  expect_true(.albers_parallels_ok(40, 40.1))
  expect_true(.albers_parallels_ok(-60.05, -60))
  expect_false(.albers_parallels_ok(85, 95))
  expect_false(.albers_parallels_ok(-95, -85))
  expect_false(.albers_parallels_ok(85, 90.000001))
  expect_false(.albers_parallels_ok(-90.000001, -85))
  # The cone test as before: parallels across the equator degenerate.
  expect_false(.albers_parallels_ok(-10, 10))
  expect_false(.albers_parallels_ok(-2, 4))
  expect_false(.albers_parallels_ok(NA_real_, 40))
  expect_false(.albers_parallels_ok(NaN, NaN))
})

test_that("a lon/lat line whose centre is at or near a pole is projected, with either purpose", {
  lines <- c(.cap_pole_lines(), .cap_near_lines())
  for (nm in names(lines)) for (p in c("distance", "area")) {
    x <- .cap_line(lines[[nm]])
    lab <- sprintf("%s, purpose = \"%s\"", nm, p)
    pick <- .cap_quiet(.pick_local_projected_crs(x, p))
    expect_false(is.null(pick$candidates), label = lab)
    expect_false(any(grepl("^Albers", pick$candidates$name)), label = lab)
    expect_identical(sum(pick$candidates$chosen), 1L, label = lab)
    expect_identical(pick$crs$input, pick$candidates$crs[pick$candidates$chosen],
                     label = lab)
    # The Lambert azimuthal candidate centred on the data places the line
    # and is the one chosen.
    expect_match(pick$candidates$name[pick$candidates$chosen],
                 "^Lambert azimuthal equal-area centred on", label = lab)
    out <- .cap_quiet(ensure_projected(x, purpose = p))
    expect_identical(sf::st_crs(out), pick$crs, label = lab)
    expect_identical(attr(out, "crs_choice"), pick$candidates, label = lab)
    xy <- sf::st_coordinates(out)
    expect_identical(nrow(xy), 2L, label = lab)
    expect_true(all(is.finite(xy[, c("X", "Y")])), label = lab)
  }
  # On the pole itself: the Lambert azimuthal CRS centred on that pole.
  for (nm in names(.cap_pole_lines())) {
    pick <- .cap_quiet(.pick_local_projected_crs(.cap_line(.cap_pole_lines()[[nm]])))
    expect_match(pick$crs$input,
                 if (nm == "south30") "+lat_0=-90.000000 " else "+lat_0=90.000000 ",
                 fixed = TRUE, label = nm)
  }
})

test_that("layers whose Albers parallels lie between -90 and 90 get the choice they got", {
  set.seed(1)
  st39 <- sf::st_as_sf(data.frame(x = stats::runif(39, -120, -75),
                                  y = stats::runif(39, 30, 47)),
                       coords = c("x", "y"), crs = 4326)
  conus <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-122, 26), c(-70, 26), c(-70, 48), c(-122, 48), c(-122, 26)))), crs = 4326))
  nc <- sf::st_transform(
    sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE), 4326)
  anti <- sf::st_as_sf(data.frame(x = c(-179, 179), y = c(-10, 10)),
                       coords = c("x", "y"), crs = 4326)
  ant <- sf::st_as_sf(data.frame(x = seq(-170, 170, by = 20), y = -75),
                      coords = c("x", "y"), crs = 4326)
  layers <- list(
    st39 = st39, conus = conus, nc = nc, anti = anti, ant = ant,
    flat40 = .cap_line(rbind(c(-100, 40), c(-80, 40))),     # fallback 35.4, 45.4
    near70 = .cap_line(rbind(c(0, 70), c(150, 70))),        # fallback 79.6, 89.6
    ctrl170 = .cap_line(rbind(c(0, 20), c(170, 20))))
  albers_seen <- character(0)
  for (nm in names(layers)) for (p in c("distance", "area")) {
    x <- layers[[nm]]
    lab <- sprintf("%s, purpose = \"%s\"", nm, p)
    got_pick <- .cap_quiet(.pick_local_projected_crs(x, p))
    got_out  <- .cap_quiet(ensure_projected(x, purpose = p))
    local({
      local_mocked_bindings(.albers_parallels_ok = .cap_cone_only)
      want_pick <- .cap_quiet(.pick_local_projected_crs(x, p))
      want_out  <- .cap_quiet(ensure_projected(x, purpose = p))
      expect_identical(got_pick, want_pick, label = lab)
      expect_identical(sf::st_crs(got_out), sf::st_crs(want_out), label = lab)
      expect_identical(attr(got_out, "crs_choice"), attr(want_out, "crs_choice"),
                       label = lab)
      expect_identical(sf::st_geometry(got_out), sf::st_geometry(want_out), label = lab)
    })
    if (any(grepl("^Albers", got_pick$candidates$name)))
      albers_seen <- c(albers_seen, lab)
  }
  # The comparison above covers layers with an Albers candidate on both paths,
  # among them the fallback parallels of a box with no height.
  for (nm in c("st39", "conus", "nc", "flat40", "near70"))
    for (p in c("distance", "area"))
      expect_true(sprintf("%s, purpose = \"%s\"", nm, p) %in% albers_seen,
                  label = sprintf("%s %s has an Albers candidate", nm, p))
})

test_that("a point layer near a pole is projected with purpose = \"area\"", {
  # A point layer, or one that stays near its zone's central meridian,
  # reaches only the purpose = "area" path.  With a centre within 5 degrees
  # of a pole that path stopped with "invalid crs" where purpose =
  # "distance" kept the UTM zone.
  pt <- function(x, y) sf::st_as_sf(data.frame(x = x, y = y),
                                    coords = c("x", "y"), crs = 4326)
  layers <- list(
    p89   = pt(10, 89),               # fallback parallels 84, 94
    p855  = pt(10, 85.5),             # 80.5, 90.5
    pm87  = pt(10, -87),              # -92, -82
    two87 = pt(c(0, 2), c(87, 87)))   # 82, 92
  for (nm in names(layers)) {
    x <- layers[[nm]]
    pick <- .cap_quiet(.pick_local_projected_crs(x, "area"))
    expect_identical(pick$candidates$name,
                     sprintf("Lambert azimuthal equal-area centred on (%.1f, %.1f)",
                             mean(sf::st_coordinates(x)[, 1]),
                             sf::st_coordinates(x)[1, 2]),
                     label = nm)
    expect_true(pick$candidates$chosen, label = nm)
    out <- .cap_quiet(ensure_projected(x, purpose = "area"))
    expect_identical(sf::st_crs(out), pick$crs, label = nm)
    expect_match(sf::st_crs(out)$input, "^\\+proj=laea ", label = nm)
    xy <- sf::st_coordinates(out)
    expect_identical(nrow(xy), nrow(x), label = nm)
    expect_true(all(is.finite(xy)), label = nm)
    # purpose = "distance" keeps the UTM zone, as it did.
    expect_match(sf::st_crs(.cap_quiet(ensure_projected(x)))$input,
                 "^EPSG:32[67]", label = nm)
  }
  # A parallel of exactly 90 is in range: the Albers row stays.
  pick85 <- .cap_quiet(.pick_local_projected_crs(pt(10, 85), "area"))
  expect_identical(pick85$candidates$name,
                   c("Lambert azimuthal equal-area centred on (10.0, 85.0)",
                     "Albers equal-area (lat_1=80.0, lat_2=90.0, lon_0=10.0)"))
})
