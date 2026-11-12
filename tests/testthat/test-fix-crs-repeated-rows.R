# tests/testthat/test-fix-crs-repeated-rows.R
# ---------------------------------------------------------------------------
# .crs_distance_error() scores a candidate CRS on the pairs of a sample of
# points.  Its counts read rows, not places.  A polygon layer gave one point
# per row, added its outline only below 40 rows and thinned the outline,
# then not densified, above 160 vertex rows: one polygon repeated over many
# rows (a long-format layer) scored by the number of repeats, and from 40
# rows every candidate scored NA and the UTM zone was kept.  A point layer
# took its sample of 40 over rows, so repeat visits to one station pushed
# other stations out of it and moved the scores and the CRS chosen.  Copies
# of one location now count once, and the outline is thinned on its distinct
# vertices.
# ---------------------------------------------------------------------------

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
box <- function(x0, y0, w, h)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h), c(x0, y0))))
oct <- function(lon0, lat0, r) {
  a <- seq(0, 2 * pi, length.out = 9)[-9]
  m <- cbind(lon0 + r * cos(a) / cos(lat0 * pi / 180), lat0 + r * sin(a))
  sf::st_polygon(list(rbind(m, m[1, ])))
}
lay <- function(gs) sf::st_sf(id = seq_along(gs), geometry = sf::st_sfc(gs, crs = 4326))
points_of <- function(m)
  sf::st_sf(id = seq_len(nrow(m)), geometry = sf::st_sfc(lapply(seq_len(nrow(m)), function(i)
    sf::st_point(m[i, ])), crs = 4326))
stations <- function() {
  set.seed(1)
  cbind(stats::runif(39, -120, -75), stats::runif(39, 30, 47))
}

test_that("one polygon repeated over many rows scores, and gets the CRS, of one row", {
  # The box over the conterminous United States got a Lambert azimuthal CRS
  # with "1.76% against the zone's 10.87%" on one row, an Albers CRS on 32
  # rows, "1.28% against the zone's 8.33%" on 33 to 39 rows (the outline
  # thinned to its corners), and from 40 rows NA for every candidate, with
  # UTM zone 15 kept.
  shapes <- list(conus = box(-122, 26, 52, 22), box = box(-110, 32, 14, 12), octagon = oct(-98, 39, 9))
  for (nm in names(shapes)) {
    one <- lay(list(shapes[[nm]]))
    x1 <- quiet(ensure_projected(one))
    e1 <- lapply(c(32613, 5070), function(code) .crs_distance_error(one, sf::st_crs(code)))
    for (k in c(2, 20, 33, 40, 100)) {
      xk <- lay(rep(list(shapes[[nm]]), k))
      lab <- paste(nm, "repeated", k, "times")
      pk <- quiet(ensure_projected(xk))
      expect_true(sf::st_crs(pk) == sf::st_crs(x1), label = lab)
      expect_identical(attr(pk, "crs_choice")$distance_error, attr(x1, "crs_choice")$distance_error,
                       label = lab)
      expect_identical(lapply(c(32613, 5070), function(code) .crs_distance_error(xk, sf::st_crs(code))),
                       e1, label = lab)
    }
  }
})

test_that("copies of one polygon listed from other vertices or the other way score as one row", {
  # A 12-gon listed from each of its vertices in both directions: no two
  # rows are exact copies, so only the thinning on distinct vertices reads
  # them as one outline.  From 13 rows the 12 distinct vertices were
  # thinned as if there were more than 160, and from 40 rows every
  # candidate scored NA.
  a12 <- seq(0, 2 * pi, length.out = 13)[-13]
  m <- cbind(-103 + 7 * cos(a12) / cos(38 * pi / 180), 38 + 7 * sin(a12))
  ring <- function(s, rev = FALSE) {
    k <- ((0:11 + s) %% 12) + 1
    mm <- m[k, ]
    if (rev) mm <- mm[12:1, ]
    sf::st_polygon(list(rbind(mm, mm[1, ])))
  }
  variants <- unlist(lapply(0:11, function(s) list(ring(s), ring(s, TRUE))), recursive = FALSE)
  one <- lay(list(ring(0)))
  for (code in c(32613, 5070)) {
    e1 <- .crs_distance_error(one, sf::st_crs(code))
    expect_true(is.finite(e1))
    for (k in c(12, 13, 24, 40, 100))
      expect_identical(.crs_distance_error(lay(rep(variants, length.out = k)), sf::st_crs(code)), e1,
                       label = paste(k, "rows in EPSG", code))
  }
})

test_that("copies interleaved with other features score as the features once", {
  three <- list(box(-122, 30, 6, 5), box(-100, 40, 7, 6), box(-80, 28, 8, 9))
  once <- lay(three)
  for (code in c(32613, 5070)) {
    e1 <- .crs_distance_error(once, sf::st_crs(code))
    expect_identical(.crs_distance_error(lay(rep(three, each = 20)), sf::st_crs(code)), e1)
    expect_identical(.crs_distance_error(lay(rep(three, 20)), sf::st_crs(code)), e1)
  }
  x <- NULL
  l <- capture_spatialkit_log(x <- quiet(ensure_projected(lay(rep(list(box(-122, 26, 52, 22)), 40)))))
  expect_true(log_has(l, "1.76% against the zone's 10.87%"))
  expect_false(log_has(l, "could not be scored"))
  expect_match(attr(x, "crs_choice")$name[attr(x, "crs_choice")$chosen], "^Lambert azimuthal")
})

test_that("repeat visits to one station leave the scores, the CRS and the log as for each station once", {
  # 39 stations scored UTM zone 14 at 5.07% and Albers at 0.61%; with
  # station 1 listed 100 times, 4.59% and 0.45%, and listed 1,001 times,
  # 2.32% and 0.35%, the 40 rows sampled holding 4 of the 39 stations.
  st <- stations()
  x1 <- NULL
  l1 <- capture_spatialkit_log(x1 <- quiet(ensure_projected(points_of(st))))
  expect_true(log_has(l1, "0.61% against the zone's 5.07%"))
  e1 <- lapply(c(32614, 5070), function(code) .crs_distance_error(points_of(st), sf::st_crs(code)))
  layers <- list(x2 = rbind(st, st[1, ]), x100 = rbind(st, st[rep(1, 99), ]),
                 x1001 = rbind(st, st[rep(1, 1000), ]), each3 = st[rep(1:39, each = 3), ])
  for (nm in names(layers)) {
    xk <- NULL
    lk <- capture_spatialkit_log(xk <- quiet(ensure_projected(points_of(layers[[nm]]))))
    expect_true(sf::st_crs(xk) == sf::st_crs(x1), label = nm)
    expect_identical(attr(xk, "crs_choice"), attr(x1, "crs_choice"), label = nm)
    expect_true(log_has(lk, "0.61% against the zone's 5.07%"), label = nm)
    expect_identical(lapply(c(32614, 5070), function(code)
      .crs_distance_error(points_of(layers[[nm]]), sf::st_crs(code))), e1, label = nm)
  }
})

test_that("repeat visits to stations on a grid, which share longitudes, score as each station once", {
  # Copies of one location are found through a sort on x, then y.  On a
  # grid many stations share a longitude, so a sort on x alone would leave
  # copies apart and count them again.
  grid <- as.matrix(expand.grid(lon = seq(-118, -76, by = 6), lat = seq(30, 45, by = 3)))
  dimnames(grid) <- NULL
  expect_identical(nrow(grid), 48L)
  once <- points_of(grid)
  twice <- points_of(rbind(grid, grid[48:1, ]))
  inter <- points_of(rbind(grid[1:20, ], grid[rep(30, 50), ], grid[21:48, ], grid[rep(30, 50), ]))
  for (code in c(32614, 5070)) {
    e1 <- .crs_distance_error(once, sf::st_crs(code))
    expect_true(is.finite(e1))
    expect_identical(.crs_distance_error(twice, sf::st_crs(code)), e1, label = paste("twice in EPSG", code))
    expect_identical(.crs_distance_error(inter, sf::st_crs(code)), e1,
                     label = paste("one station 101 times in EPSG", code))
  }
  ch1 <- attr(quiet(ensure_projected(once)), "crs_choice")
  expect_identical(attr(quiet(ensure_projected(twice)), "crs_choice"), ch1)
  expect_identical(attr(quiet(ensure_projected(inter)), "crs_choice"), ch1)
})

test_that("layers with no repeated location score as before", {
  # The figures recorded before the change, to the last bits the change of
  # the geodesic's stopping test moves (test-fix-geod-short-pairs.R).
  lon <- NULL; lat <- NULL; r <- NULL
  blob <- function(lon0, lat0, r, nv, wig) {
    aa <- seq(0, 2 * pi, length.out = nv + 1L)[-(nv + 1L)]
    rr <- r * (1 + wig * sin(5 * aa + lon0) + 0.5 * wig * cos(11 * aa + lat0))
    m  <- cbind(lon0 + rr * cos(aa) / cos(lat0 * pi / 180), lat0 + rr * sin(aa))
    sf::st_polygon(list(rbind(m, m[1, ])))
  }
  set.seed(49); lon <- stats::runif(5, -118, -72); lat <- stats::runif(5, 27, 48); r <- stats::runif(5, 1, 3)
  outlines5 <- lay(lapply(1:5, function(i) blob(lon[i], lat[i], r[i], 200, 0.15)))
  grid4 <- lay(list(box(-110, 35, 3.7, 2.9), box(-106.3, 35, 3.7, 2.9),
                    box(-110, 37.9, 3.7, 2.9), box(-106.3, 37.9, 3.7, 2.9)))
  three <- lay(list(box(-122, 30, 6, 5), box(-100, 40, 7, 6), box(-80, 28, 8, 9)))
  set.seed(7)
  squares300 <- lay(lapply(1:300, function(i) box(stats::runif(1, -122, -70), stats::runif(1, 26, 48), 0.25, 0.25)))
  st39 <- points_of(stations())
  pins <- list(
    list(lay(list(box(-122, 26, 52, 22))), 0.16343244114212907, 0.015629578061724381),
    list(grid4, 0.0021318696389771574, 0.009771734310373148),
    list(three, 0.13713877145650955, 0.00956240728386204),
    list(outlines5, 0.031026601372790896, 0.0094795167675296721),
    list(squares300, 0.13828976658735592, 0.0096177977396322278),
    list(st39, 0.082827299527691745, 0.00968393880054641))
  for (k in seq_along(pins)) {
    expect_equal(.crs_distance_error(pins[[k]][[1L]], sf::st_crs(32613)), pins[[k]][[2L]],
                 tolerance = 1e-8, label = paste("layer", k, "in EPSG:32613"))
    expect_equal(.crs_distance_error(pins[[k]][[1L]], sf::st_crs(5070)), pins[[k]][[3L]],
                 tolerance = 1e-8, label = paste("layer", k, "in EPSG:5070"))
  }
  ch <- attr(quiet(ensure_projected(three)), "crs_choice")
  expect_match(ch$name[ch$chosen], "^Albers")
  expect_equal(ch$distance_error, c(0.10037393466258582, 0.025961369931918332, 0.0079429376127582163),
               tolerance = 1e-8)
})

test_that("the help says copies of one location count once and pairs under 1 mm are left out", {
  # The source page from the source tree, the installed page under R CMD check,
  # which runs the tests where there is no man/.
  file <- test_path("..", "..", "man", "ensure_projected.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[["ensure_projected.Rd"]]
  expect_s3_class(rd, "Rd")
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("comparing planar with geodesic distances over the pairs at least ",
                           "1 mm apart."), fixed = TRUE)
  expect_match(txt, paste0("plus the vertices of its outline when fewer than 40 distinct points ",
                           "remain, so a single study-area polygon is scored too, however many ",
                           "rows repeat it; copies of one location count once, and rows without ",
                           "a finite location are left out."), fixed = TRUE)
  expect_false(grepl("when it has fewer than 40 features", txt, fixed = TRUE))
})
