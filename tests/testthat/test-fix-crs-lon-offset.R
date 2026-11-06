# tests/testthat/test-fix-crs-lon-offset.R
# ---------------------------------------------------------------------------
# BUG-210.  The CRS chooser takes the UTM zone from the centre on the
# sphere and measured how far the box reaches from that zone's central
# meridian as max(abs(lon_min - cm), abs(lon_max - cm)), not wrapped.  For
# a layer across a pole the centre's longitude, and so the zone, can lie
# far outside the box, and the comparison path's warning said two lines
# from -10 to 170 at 70 and 75 degrees north "reach 257.0 deg from the
# central meridian of UTM zone 16".  No longitude lies more than 180
# degrees from a meridian.
#
# The offset is now the largest distance on the circle from the meridian to
# any longitude of the box (.lon_offset_from()), and the 5-degree gate and
# both warnings read it.  The choice is unchanged except for a box within 5
# degrees of the meridian across the antimeridian, which the rule assigns
# to the zone.
# ---------------------------------------------------------------------------

.clo_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.clo_line <- function(...)
  sf::st_sf(id = seq_along(list(...)),
            geometry = sf::st_sfc(lapply(list(...), sf::st_linestring), crs = 4326))
.clo_pts <- function(x, y)
  sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = 4326)
# The figure before the fix.
.clo_old_offset <- function(lon_min, lon_max, cm)
  max(abs(lon_min - cm), abs(lon_max - cm))
# The chooser's output and its log lines without their time stamps.
.clo_run <- function(x, purpose = "distance") {
  pick <- NULL
  lines <- capture_spatialkit_log(
    pick <- .clo_quiet(.pick_local_projected_crs(x, purpose)))
  list(pick = pick, lines = sub("^[A-Z]+ \\[[^]]*\\] ", "", lines))
}

test_that(".lon_offset_from() gives the largest offset on the circle, never above 180", {
  # The meridian opposite -87 (93) lies in the box.
  expect_identical(.lon_offset_from(-10, 170, -87), 180)
  expect_identical(.lon_offset_from(100, 170, -87), 173)
  expect_identical(.lon_offset_from(-100, -80, -87), 13)
  # Across the antimeridian: 178 lies 5 degrees from -177, not above the gate.
  expect_identical(.lon_offset_from(178, 180, -177), 5)
  expect_identical(.lon_offset_from(181, 183, -177), 2)
  # A box that starts or ends on the opposite meridian reaches 180.
  expect_identical(.lon_offset_from(93, 100, -87), 180)
  expect_identical(.lon_offset_from(80, 93, -87), 180)
  # Where the old figure was at most 180 and the opposite meridian lies
  # outside the box, the figure is the old one.
  expect_identical(.lon_offset_from(-124, -67, -93), 31)
  expect_identical(.lon_offset_from(0, 180, 93), 93)
  # Against a brute-force maximum over the box's longitudes.
  set.seed(210)
  for (i in 1:300) {
    a <- stats::runif(1, -180, 180)
    b <- a + stats::runif(1, 0, 180)
    cm <- 6 * sample(1:60, 1) - 183
    lon <- c(seq(a, b, length.out = 2001), cm + 180 + 360 * (-2:2))
    lon <- lon[lon >= a & lon <= b]
    d <- abs(lon - cm) %% 360
    off <- .lon_offset_from(a, b, cm)
    expect_true(off >= 0 && off <= 180)
    expect_equal(off, max(pmin(d, 360 - d)))
  }
  # A box or meridian that is not finite gives a figure that is not finite.
  expect_false(is.finite(.lon_offset_from(NA_real_, NA_real_, -87)))
  expect_false(is.finite(.lon_offset_from(NaN, 10, -87)))
  expect_false(is.finite(.lon_offset_from(0, 10, NA_real_)))
})

test_that("a layer across the pole is said to reach 180 degrees from its zone's meridian", {
  layers <- list(
    two_lines = .clo_line(rbind(c(-10, 70), c(170, 70)), rbind(c(-10, 75), c(170, 75))),
    one_line  = .clo_line(rbind(c(-10, 70), c(170, 70))))
  for (nm in names(layers)) {
    got <- .clo_run(layers[[nm]])
    msg <- grep("deg from the central meridian", got$lines, value = TRUE)
    expect_length(msg, 1L)
    expect_match(msg, "reaches 180.0 deg from the central meridian of UTM zone 16 (180.0 deg of longitude in total)",
                 fixed = TRUE, label = nm)
    figs <- as.numeric(regmatches(msg, gregexpr("[0-9]+\\.[0-9](?= deg)", msg, perl = TRUE))[[1L]])
    expect_true(length(figs) >= 1L && all(figs <= 180), label = nm)
    # The CRS, the candidates and their scores are the ones before the fix.
    local({
      local_mocked_bindings(.lon_offset_from = .clo_old_offset)
      want <- .clo_run(layers[[nm]])
      expect_identical(got$pick, want$pick, label = nm)
      expect_match(grep("deg from the central meridian", want$lines, value = TRUE),
                   "reaches 257.0 deg", fixed = TRUE, label = nm)
    })
    expect_match(got$pick$candidates$name[got$pick$candidates$chosen],
                 "^Lambert azimuthal equal-area centred on \\(-90\\.0, 90\\.0\\)", label = nm)
  }
})

test_that("layers whose offset was at most 180 get the choice and the log lines they got", {
  set.seed(1)
  conus <- sf::st_as_sf(data.frame(lon = stats::runif(50, -124, -67),
                                   lat = stats::runif(50, 25, 49)),
                        coords = c("lon", "lat"), crs = 4326)
  set.seed(1)
  st39 <- .clo_pts(stats::runif(39, -120, -75), stats::runif(39, 30, 47))
  nc <- sf::st_transform(
    sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE), 4326)
  layers <- list(
    conus = conus, st39 = st39, nc = nc,
    # The layers of test-fix-crs-albers-pole.R.
    line20  = .clo_line(rbind(c(0, 20), c(180, 20))),
    line45  = .clo_line(rbind(c(0, 45), c(180, 45))),
    line60  = .clo_line(rbind(c(0, 60), c(180, 60))),
    south30 = .clo_line(rbind(c(0, -30), c(180, -30))),
    north80 = .clo_line(rbind(c(0, 80), c(160, 80))),
    south84 = .clo_line(rbind(c(0, -84), c(120, -84))),
    near70  = .clo_line(rbind(c(0, 70), c(150, 70))),
    ctrl170 = .clo_line(rbind(c(0, 20), c(170, 20))),
    flat40  = .clo_line(rbind(c(-100, 40), c(-80, 40))),
    anti = .clo_pts(c(-179, 179), c(-10, 10)),
    ant  = .clo_pts(seq(-170, 170, by = 20), -75),
    p180 = .clo_pts(rep(180, 5), seq(10, 20, length.out = 5)),
    pm180 = .clo_pts(rep(-180, 5), seq(10, 20, length.out = 5)),
    wake = .clo_pts(c(-78.7, -78.5), c(35.7, 35.9)),
    wide = .clo_pts(c(-100, -92), c(38, 42)))
  for (nm in names(layers)) for (p in c("distance", "area")) {
    lab <- sprintf("%s, purpose = \"%s\"", nm, p)
    got <- .clo_run(layers[[nm]], p)
    local({
      local_mocked_bindings(.lon_offset_from = .clo_old_offset)
      want <- .clo_run(layers[[nm]], p)
      expect_identical(got$pick, want$pick, label = lab)
      expect_identical(got$lines, want$lines, label = lab)
    })
  }
  # The CONUS layer still logs the switch to Albers.
  got <- .clo_run(conus)
  expect_match(got$pick$crs$input, "+proj=aea", fixed = TRUE)
  expect_true(log_has(got$lines, "reaches 30\\.2 deg from the central meridian of UTM zone 15"))
})

test_that("a box within 5 degrees of the zone's meridian across the antimeridian keeps the zone", {
  # Longitudes 181 and 183 are -179 and -177: 2 degrees from the central
  # meridian of UTM zone 1.  The old figure, 360, sent the layer to the
  # comparison, which chose an Albers CRS.
  beyond <- .clo_pts(c(181, 183), c(10, 11))
  inside <- .clo_pts(c(-179, -177), c(10, 11))
  got <- .clo_run(beyond)
  expect_identical(got$pick$crs, sf::st_crs(32601))
  expect_identical(got$pick$candidates$name, "UTM zone 1")
  expect_false(log_has(got$lines, "central meridian of UTM zone"))
  expect_identical(.clo_run(inside)$pick$crs, sf::st_crs(32601))
  out <- .clo_quiet(ensure_projected(beyond))
  expect_identical(sf::st_crs(out), sf::st_crs(32601))
  expect_true(all(is.finite(sf::st_coordinates(out))))
  local({
    local_mocked_bindings(.lon_offset_from = .clo_old_offset)
    old <- .clo_run(beyond)
    expect_true(log_has(old$lines, "reaches 360\\.0 deg from the central meridian of UTM zone 1"))
    expect_match(old$pick$crs$input, "+proj=aea", fixed = TRUE)
  })
  # Longitudes -181 and -183 are 179 and 177: 2 degrees from the central
  # meridian of UTM zone 60.
  beyond_w <- .clo_pts(c(-181, -183), c(10, 11))
  inside_w <- .clo_pts(c(179, 177), c(10, 11))
  got <- .clo_run(beyond_w)
  expect_identical(got$pick$crs, sf::st_crs(32660))
  expect_identical(got$pick$candidates$name, "UTM zone 60")
  expect_false(log_has(got$lines, "central meridian of UTM zone"))
  expect_identical(.clo_run(inside_w)$pick$crs, sf::st_crs(32660))
  out <- .clo_quiet(ensure_projected(beyond_w))
  expect_identical(sf::st_crs(out), sf::st_crs(32660))
  expect_true(all(is.finite(sf::st_coordinates(out))))
  local({
    local_mocked_bindings(.lon_offset_from = .clo_old_offset)
    old <- .clo_run(beyond_w)
    expect_true(log_has(old$lines, "reaches 360\\.0 deg from the central meridian of UTM zone 60"))
    expect_match(old$pick$crs$input, "+proj=aea", fixed = TRUE)
  })
})

test_that("a layer with longitudes beyond 180 wider than 5 degrees is said to reach at most 180 degrees", {
  # Longitudes 170, 185 and 200 are 170, -175 and -160: the farthest lies
  # 17 degrees from the central meridian of UTM zone 1 (-177).  The old
  # figure was 377.  The layer goes to the comparison either way, and the
  # choice is the same.
  wide <- .clo_pts(c(170, 185, 200), c(10, 12, 14))
  got <- .clo_run(wide)
  msg <- grep("deg from the central meridian", got$lines, value = TRUE)
  expect_length(msg, 1L)
  expect_match(msg, "reaches 17.0 deg from the central meridian of UTM zone 1 ", fixed = TRUE)
  local({
    local_mocked_bindings(.lon_offset_from = .clo_old_offset)
    want <- .clo_run(wide)
    expect_identical(got$pick, want$pick)
    expect_match(grep("deg from the central meridian", want$lines, value = TRUE),
                 "reaches 377.0 deg from the central meridian of UTM zone 1 ", fixed = TRUE)
  })
})
