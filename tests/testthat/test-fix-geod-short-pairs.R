# tests/testthat/test-fix-geod-short-pairs.R
# ---------------------------------------------------------------------------
# .geod_distance() stopped Vincenty's iteration once a pass moved lambda by
# less than 1e-12 rad.  On an east-west pair whose longitude difference is
# under about 1e-12 / (f cos^2 U) rad the first pass already moved it by
# less, the loop ended, and sigma, taken from the starting lambda, came out
# short by about f cos^2 U: 0.34% at the equator, 0.20% at latitude 40, on
# pairs of a few millimetres; pairs out to about a metre stopped after the
# second pass, up to 1.1e-5 short.  .crs_distance_error() takes the worst
# ratio over its pairs, so one such pair raised the score of every
# candidate.  The change is now measured against the pair's longitude
# difference, taken the short way round.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
points_of <- function(m)
  sf::st_sf(id = seq_len(nrow(m)), geometry = sf::st_sfc(lapply(seq_len(nrow(m)), function(i)
    sf::st_point(m[i, ])), crs = 4326))
a_wgs <- 6378137; f_wgs <- 1 / 298.257223563; e2_wgs <- f_wgs * (2 - f_wgs)

test_that("east-west pairs of a few millimetres agree with the ellipsoid's local metric", {
  # -3.35e-3 at the equator, -1.97e-3 at latitude 40, -8.4e-4 at 60 and
  # -1.0e-4 at 80; the pairs that stopped on the second pass, -1.1e-5.
  for (lat in c(0, 40, 60, 80)) {
    N <- a_wgs / sqrt(1 - e2_wgs * sin(lat * pi / 180)^2)
    for (gap in c(1e-8, 3e-8, 1e-7, 1e-6, 1e-4)) {
      g <- (-104 + gap) - (-104)
      ref <- N * cos(lat * pi / 180) * g * pi / 180
      d <- .geod_distance(-104, lat, -104 + gap, lat)
      expect_lt(abs(d / ref - 1), 1e-9, label = sprintf("%g degrees east at latitude %g", gap, lat))
    }
  }
})

test_that("north-south pairs agree with the meridian's local metric and stop on the first pass", {
  # A pair along a meridian has a longitude difference of 0, so its test is
  # 0 <= 0; with `<` it would run 200 passes and fall back to the sphere.
  for (lat in c(0, 40, 60, 80)) {
    for (gap in c(1e-8, 1e-6, 1e-4)) {
      phi <- (lat + gap / 2) * pi / 180
      M <- a_wgs * (1 - e2_wgs) / (1 - e2_wgs * sin(phi)^2)^1.5
      ref <- M * ((lat + gap) - lat) * pi / 180
      d <- .geod_distance(-104, lat, -104, lat + gap)
      expect_lt(abs(d / ref - 1), 1e-5, label = sprintf("%g degrees north at latitude %g", gap, lat))
    }
  }
})

test_that("longer pairs do not move", {
  # Recorded before the change.
  lon1 <- c(-113.4424, 72.8547, 26.3975, -119.5013, 159.7822, 159.651, -133.5028, 120.0416,
            -11.5133, 17.9941, 18.9627, -93.9979)
  lat1 <- c(44.2873, -54.2606, -16.102, 60.1032, 80.9877, -46.6097, -9.3824, -72.2535, 27.5228,
            -19.1166, 57.2712, -59.4148)
  lon2 <- c(-131.7697, 71.5075, -15.6929, -136.6538, 215.2995, 115.5356, -192.2531, 79.7987,
            25.7098, 62.2574, 20.6765, -78.7343)
  lat2 <- c(71.8416, -71.471, -2.7239, 32.1408, 89, -62.8488, -40.1757, -89, 63.0461, 4.4725, 89,
            -71.4877)
  before <- c(3221871.5613050521, 1919289.2578919767, 4840153.932565351, 3355676.4135925518,
              947745.37881610333, 3271687.6232934371, 6701416.9154408248, 1897619.8502407966,
              4777806.6472374462, 5506095.236437263, 3540207.8857926973, 1510659.9232936671)
  expect_equal(.geod_distance(lon1, lat1, lon2, lat2), before, tolerance = 1e-10)
})

test_that("nearly antipodal pairs fall back to the sphere as before", {
  # The first three converge after about 190 passes under a test relative
  # to a longitude difference near pi, looser than 1e-12 rad; the test stays
  # in radians from 1 rad of longitude up.
  pairs <- rbind(c(28.3263, 32.5437, 208.3050, -33.0132), c(-49.9813, 27.1479, 130.2741, -27.6363),
                 c(49.7471, -7.3395, 230.2526, 6.8660), c(0, 0, 179.7, 0.2), c(10, 20, -170.2, -20.1))
  rad <- pi / 180
  for (i in seq_len(nrow(pairs))) {
    p <- pairs[i, ]
    la1 <- p[2] * rad; la2 <- p[4] * rad; dl <- (p[3] - p[1]) * rad
    sph <- 6371008.8 * acos(pmin(1, pmax(-1, sin(la1) * sin(la2) + cos(la1) * cos(la2) * cos(dl))))
    expect_identical(.geod_distance(p[1], p[2], p[3], p[4]), sph, label = paste("pair", i))
  }
})

test_that("a point a few millimetres east of another leaves the score of the layer as it was", {
  # 30 points between latitudes 33 and 41 score 0.000398 in UTM zone 13;
  # with a 31st point 2e-8 degrees (1.8 mm) east of one, 0.001965, which
  # ensure_projected() reported in `crs_choice`.  Between latitudes 58 and
  # 62, 0.000400 and 0.000548 with a point 5e-8 degrees (3 mm) east of one.
  set.seed(7)
  pts <- cbind(stats::runif(30, -107.5, -102.5), stats::runif(30, 33, 41))
  plus <- rbind(pts, pts[5, ] + c(2e-8, 0))
  e0 <- .crs_distance_error(points_of(pts), sf::st_crs(32613))
  expect_lt(abs(.crs_distance_error(points_of(plus), sf::st_crs(32613)) - e0), 1e-6)
  expect_lt(abs(.crs_distance_error(points_of(plus), sf::st_crs(5070)) -
                .crs_distance_error(points_of(pts), sf::st_crs(5070))), 1e-6)
  ch0 <- attr(quiet(ensure_projected(points_of(pts))), "crs_choice")
  ch  <- attr(quiet(ensure_projected(points_of(plus))), "crs_choice")
  expect_identical(ch$name[ch$chosen], "UTM zone 13")
  expect_lt(abs(ch$distance_error[ch$chosen] - ch0$distance_error[ch0$chosen]), 1e-6)
  set.seed(9)
  pts60 <- cbind(stats::runif(30, -107.5, -102.5), stats::runif(30, 58, 62))
  plus60 <- rbind(pts60, pts60[5, ] + c(5e-8, 0))
  for (code in c(32613, 3035))
    expect_lt(abs(.crs_distance_error(points_of(plus60), sf::st_crs(code)) -
                  .crs_distance_error(points_of(pts60), sf::st_crs(code))), 1e-6,
              label = paste("latitudes 58 to 62 in EPSG", code))
})

test_that("east-west pairs across longitude 180 agree with the same pairs away from it", {
  # The test was scaled on the longitude difference as written, close to
  # 2 pi for a short pair across 180, so the first pass stopped such a pair
  # as before and it came out 0.34% short at the equator and 0.20% at
  # latitude 40.  A pair whose longitudes differ by more than 180 degrees is
  # now measured the short way round.
  for (lat in c(0, 40)) {
    N <- a_wgs / sqrt(1 - e2_wgs * sin(lat * pi / 180)^2)
    for (gap in c(1e-8, 2e-8, 5e-8, 1e-7, 1e-6)) {
      lon1 <- 180 - gap / 2; lon2 <- -180 + gap / 2
      g <- (180 - lon1) + (lon2 + 180)
      ref <- N * cos(lat * pi / 180) * g * pi / 180
      lab <- sprintf("%g degrees across 180 at latitude %g", gap, lat)
      expect_lt(abs(.geod_distance(lon1, lat, lon2, lat) / ref - 1), 1e-9, label = lab)
      expect_lt(abs(.geod_distance(lon2, lat, lon1, lat) / ref - 1), 1e-9, label = paste(lab, "listed west first"))
    }
  }
  # A gap of an odd number of steps of 2^-45 degrees, the spacing of
  # doubles near 180: written as a difference of about 360 degrees it
  # rounds to an even number of steps, 3e-6 of a 1e-8 degree gap.
  lon1 <- 180 - 176000 * 2^-45; lon2 <- -180 + 176001 * 2^-45
  ref <- a_wgs * 352001 * 2^-45 * pi / 180
  expect_lt(abs(.geod_distance(lon1, 0, lon2, 0) / ref - 1), 1e-9)
  expect_lt(abs(.geod_distance(lon2, 0, lon1, 0) / ref - 1), 1e-9)
  # Longer pairs across 180 measure as the same pairs shifted away from it.
  expect_equal(.geod_distance(c(170, -175), c(-17, 40), c(-170, 175), c(-16, 41)),
               .geod_distance(c(-10, 5), c(-17, 40), c(10, -5), c(-16, 41)), tolerance = 1e-12)
  expect_identical(.geod_distance(180, 40, -180, 40), 0)
})

test_that("a pair of stations a millimetre apart across 180 leaves the score of the layer as it was", {
  # 30 stations around Fiji score 0.000992 for the Lambert azimuthal CRS
  # chosen; with two points 1.2e-8 degrees (1.3 mm) apart either side of
  # 180, 0.003083, and with the same pair at longitude 178, 0.000989.
  set.seed(4)
  lon <- stats::runif(30, 175, 185); lon <- ifelse(lon > 180, lon - 360, lon)
  pts <- cbind(lon, stats::runif(30, -20, -14))
  across <- rbind(pts, c(180 - 6e-9, -17), c(-180 + 6e-9, -17))
  at178  <- rbind(pts, c(178 - 6e-9, -17), c(178 + 6e-9, -17))
  chosen <- function(m) {
    ch <- attr(quiet(ensure_projected(points_of(m))), "crs_choice")
    ch$distance_error[ch$chosen]
  }
  e0 <- chosen(pts)
  expect_lt(abs(chosen(across) - e0), 1e-5)
  expect_lt(abs(chosen(at178) - e0), 1e-5)
})
