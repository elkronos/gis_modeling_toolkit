# tests/testthat/test-fix-far-rows-moved.R
# ---------------------------------------------------------------------------
# A coordinate beyond 1e150 marks a row without a location.  The bound is read
# on the coordinates as given, before the move to the working CRS: a move
# between CRSs can bring such a coordinate back finite (a northing of 1e200 in
# EPSG:32616 lands at (87304, 5766268) in EPSG:32617, one in EPSG:3857 at
# (5e5, 9997965)), and a transform of it is not a place the caller gave.
# Rows of `newdata` in cell_values_at(), and rows of `points` there and in
# summarize_by_cell(estimate = "pooled"), are each held to what the same row
# emptied gives.
#
# The layer is that of test-core-values-at.R: 600 locations (seed 6102), 120
# held out, the density kind from 3 k-means++ starts under seed 123, the count
# from a size of 4.5, and pooled values.  The working CRS is EPSG:32617.
# ---------------------------------------------------------------------------

.fm_seed <- function(seed, env = parent.frame())
  withr::local_seed(seed, .local_envir = env, .rng_kind = "Mersenne-Twister",
                    .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")

# A call with its warnings kept and muffled.
.fm_catch <- function(expr) {
  got <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    got[[length(got) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = got)
}
.fm_classes <- function(r) vapply(r$warnings, function(w) class(w)[1L], character(1))
.fm_texts <- function(r) vapply(r$warnings, conditionMessage, character(1))
# The columns of a frame and nothing else: no row names, no record.
.fm_cols <- function(d) {
  d <- as.list(d)
  attributes(d) <- list(names = names(d))
  d
}

.fm_store <- new.env(parent = emptyenv())
.fm_layer <- function() {
  if (is.null(.fm_store$layer)) {
    d <- local({
      .fm_seed(6102)
      n <- 600
      x <- stats::runif(n, 0, 1e5); y <- stats::runif(n, 0, 1e5)
      z <- sin(x / 2e4) + cos(y / 1.5e4) + stats::rnorm(n, 0, 0.6)
      te <- sample.int(n, 120)
      list(d = data.frame(x = 5e5 + x, y = 5e6 + y, z = z), te = te)
    })
    pts  <- sf::st_as_sf(d$d, coords = c("x", "y"), crs = 32617)
    rows <- pts[-d$te, ]
    prof <- suppressWarnings(
      resolution_profile(rows, "z", partition = "density", init = "kmeans++", nstart = 3,
                         sac = "auto", seed = 123, estimate = "pooled", obs_per_cell = 4.5))
    sel   <- select_resolution(prof, criterion = "size")
    seeds <- get_voronoi_seeds(method = "kmeans", n = sel, sample_points = rows)
    tess  <- build_tessellation(seeds, method = "voronoi", quiet = TRUE)
    asg   <- assign_features_to_polygons(rows, tess$cells, polygon_id_col = "cell_id",
                                         quiet = TRUE)
    pooled <- summarize_by_cell(asg, "z", id_col = "cell_id", cells_sf = tess$cells,
                                estimate = "pooled")
    seeds$cell_id <- tess$index
    .fm_store$layer <- list(held = pts[d$te, ], tess = tess, seeds = seeds, asg = asg,
                            pooled = pooled)
  }
  .fm_store$layer
}

# `layer` in EPSG `epsg`, its row `i` replaced by the point `at` (an empty one
# when `at` is NULL).
.fm_put <- function(layer, epsg, at, i = 1L) {
  l <- sf::st_transform(layer, epsg)
  g <- sf::st_geometry(l)
  g[i] <- sf::st_sfc(if (is.null(at)) sf::st_point() else sf::st_point(at), crs = epsg)
  sf::st_set_geometry(l, g)
}

.fm_at <- function(L, newdata, points = L$asg, seeds = L$seeds)
  cell_values_at(newdata, L$pooled, points, "z", id_col = "cell_id", seeds = seeds)


# The far rows of the tests below: a northing of 1e200 in EPSG:32616 and in
# EPSG:3857, and an easting of 1e200 in EPSG:3413 (polar stereographic),
# which lands at (500000, -9997965) in EPSG:32617.
.fm_far <- list(list(32616, c(5e5, 1e200)), list(3857, c(5e5, 1e200)),
                list(3857, c(5e5, -1e200)), list(3413, c(1e200, 0)))

test_that("the moves of the item bring a coordinate beyond 1e150 back finite", {
  for (cs in .fm_far) {
    far <- sf::st_sfc(sf::st_point(cs[[2L]]), crs = cs[[1L]])
    xy  <- sf::st_coordinates(sf::st_transform(far, 32617))[1L, 1:2]
    expect_true(all(is.finite(xy)), info = cs[[1L]])
    expect_true(all(abs(xy) < 1e150), info = cs[[1L]])
  }
})

test_that("a row of `newdata` beyond 1e150 in another CRS than the working one has NA in every column", {
  L <- .fm_layer()
  nd0 <- L$held[1:6, ]
  cases <- c(.fm_far, list(list(32616, c(5e5, -1e200))))
  for (cs in cases) {
    epsg <- cs[[1L]]; at <- cs[[2L]]
    tag <- sprintf("EPSG:%d, (%g, %g)", epsg, at[1L], at[2L])
    far <- .fm_put(nd0, epsg, at)
    emp <- .fm_put(nd0, epsg, NULL)
    for (seeds in list(L$seeds, NULL)) {
      how <- paste(tag, if (is.null(seeds)) "polygons" else "seeds")
      r <- .fm_catch(.fm_at(L, far, seeds = seeds))
      e <- .fm_catch(.fm_at(L, emp, seeds = seeds))
      expect_true(all(is.na(unlist(as.data.frame(r$value)[1L, ]))), info = how)
      expect_identical(r$value, e$value, info = how)
      expect_identical(.fm_classes(r), .fm_classes(e), info = how)
      # The other rows are as without the far row.
      rest <- .fm_at(L, far[-1L, ], seeds = seeds)
      expect_identical(lapply(.fm_cols(r$value), `[`, -1L), .fm_cols(rest), info = how)
    }
  }
})

test_that("a row of `newdata` without a location in the working CRS itself, or only once moved there, has none still", {
  L <- .fm_layer()
  nd0 <- L$held[1:6, ]
  for (seeds in list(L$seeds, NULL)) {
    # In the working CRS itself, as before.
    far <- .fm_put(nd0, 32617, c(5e5, 1e200))
    emp <- .fm_put(nd0, 32617, NULL)
    r <- .fm_at(L, far, seeds = seeds)
    expect_true(all(is.na(unlist(as.data.frame(r)[1L, ]))))
    expect_identical(r, .fm_at(L, emp, seeds = seeds))
    # Metres stamped EPSG:4326: finite as given, off the globe once moved.
    ll <- sf::st_transform(nd0, 4326)
    g <- sf::st_geometry(ll)
    g[1L] <- sf::st_sfc(sf::st_point(c(5e5, 5e6)), crs = 4326)
    stamped <- sf::st_set_geometry(ll, g)
    r <- .fm_at(L, stamped, seeds = seeds)
    expect_true(all(is.na(unlist(as.data.frame(r)[1L, ]))))
    expect_identical(r, .fm_at(L, .fm_put(nd0, 4326, NULL), seeds = seeds))
  }
})

test_that("a row of `points` beyond 1e150 in another CRS than the working one is a row without a location in cell_values_at()", {
  L <- .fm_layer()
  for (cs in .fm_far) {
    epsg <- cs[[1L]]
    far <- .fm_put(L$asg, epsg, cs[[2L]])
    emp <- .fm_put(L$asg, epsg, NULL)
    for (seeds in list(L$seeds, NULL)) {
      how <- sprintf("EPSG:%d, (%g, %g), %s", epsg, cs[[2L]][1L], cs[[2L]][2L],
                     if (is.null(seeds)) "polygons" else "seeds")
      r <- .fm_catch(.fm_at(L, L$held, points = far, seeds = seeds))
      e <- .fm_catch(.fm_at(L, L$held, points = emp, seeds = seeds))
      expect_identical(r$value, e$value, info = how)
      expect_identical(.fm_texts(r), .fm_texts(e), info = how)
      expect_identical(attr(r$value, "values_at")$n, 479L, info = how)
    }
  }
  # The reader itself.
  far <- .fm_put(L$asg, 32616, c(5e5, 1e200))
  xy <- .pooled_xy(far, sf::st_crs(32617))$xy
  expect_true(all(is.na(xy[1L, ])))
  expect_true(all(is.finite(xy[-1L, ])))
  expect_identical(xy[-1L, ], .pooled_xy(sf::st_transform(L$asg[-1L, ], 32616),
                                         sf::st_crs(32617))$xy)
  # In the working CRS itself, as before.
  far17 <- .fm_put(L$asg, 32617, c(5e5, 1e200))
  expect_identical(.fm_at(L, L$held, points = far17),
                   .fm_at(L, L$held, points = .fm_put(L$asg, 32617, NULL)))
})

test_that("summarize_by_cell(estimate = \"pooled\") leaves a row of `points` beyond 1e150 in another CRS out of the pooled columns, with the warning", {
  L <- .fm_layer()
  expect_false(is.null(.count_rule_get(sf::st_transform(L$asg, 32616))))
  run <- function(l) .fm_catch(summarize_by_cell(l, "z", id_col = "cell_id",
                                                 cells_sf = L$tess$cells, estimate = "pooled"))
  skipped <- paste0("summarize_by_cell(): 1 row(s) in 1 cell(s) have no location (an empty ",
                    "point, or a coordinate that is not finite or lies beyond 1e150). They ",
                    "count in `n` and the plain columns and take no part in the pooled ",
                    "columns.")
  # A layer in a geographic CRS: the text also names a latitude beyond 90
  # (BUG-223), which the row in metres stamped EPSG:4326 has.
  skipped_ll <- paste0("summarize_by_cell(): 1 row(s) in 1 cell(s) have no location (an empty ",
                       "point, a coordinate that is not finite or lies beyond 1e150, or a ",
                       "latitude beyond 90 degrees). They count in `n` and the plain columns ",
                       "and take no part in the pooled columns.")
  cases <- list(list(32616, c(5e5, 1e200)), list(3857, c(5e5, -1e200)),
                list(3413, c(1e200, 0)), list(32617, c(5e5, 1e200)),
                # Metres stamped EPSG:4326: finite as given, off the globe once moved.
                list(4326, c(5e5, 5e6)))
  for (cs in cases) {
    tag <- sprintf("EPSG:%d", cs[[1L]])
    r <- run(.fm_put(L$asg, cs[[1L]], cs[[2L]]))
    e <- run(.fm_put(L$asg, cs[[1L]], NULL))
    expect_identical(r$value, e$value, info = tag)
    expect_identical(.fm_classes(r), .fm_classes(e), info = tag)
    expect_identical(.fm_texts(r), .fm_texts(e), info = tag)
    expect_true("spatialkit_pool_rows_skipped" %in% .fm_classes(r), info = tag)
    expect_true((if (cs[[1L]] == 4326) skipped_ll else skipped) %in% .fm_texts(r), info = tag)
    expect_identical(attr(r$value, "pooled")$columns$resp_z$skipped, 1L, info = tag)
    expect_identical(attr(r$value, "pooled")$columns$resp_z$n, 479L, info = tag)
  }
})

# 300 locations in degrees near (-84.2, 40.2), in 25 cells named by hand, and
# no record of a count: summarize_by_cell() chooses the CRS from the layer.
.fm_lonlat <- function() {
  .fm_seed(11)
  n <- 300
  lon <- stats::runif(n, -84.5, -84.0); lat <- stats::runif(n, 40, 40.5)
  cid <- paste0(floor((lon + 84.5) * 6), "_", floor((lat - 40) * 6))
  z <- sin(lon * 10) + cos(lat * 8) + stats::rnorm(n, 0, 0.5)
  sf::st_as_sf(data.frame(lon, lat, z, cell_id = cid), coords = c("lon", "lat"), crs = 4326)
}

# `layer` with its row `i` replaced by the point `at`, or by an empty point,
# in the CRS it has.
.fm_set <- function(layer, at, i = 1L) {
  g <- sf::st_geometry(layer)
  g[i] <- sf::st_sfc(if (is.null(at)) sf::st_point() else sf::st_point(at),
                     crs = sf::st_crs(layer))
  sf::st_set_geometry(layer, g)
}

test_that("a row of `points` beyond 1e150 as given takes no part in choosing the CRS of the pooled summary", {
  ll <- .fm_lonlat()
  expect_null(.count_rule_get(ll))
  run <- function(l) .fm_catch(summarize_by_cell(l, "z", id_col = "cell_id",
                                                 estimate = "pooled"))
  nc <- sf::st_set_crs(ll, NA)
  # In EPSG:4326 a longitude of 1e200 put the centre of the projection off
  # the globe, and sf refused it; a latitude of 1e200 raised sf's warnings on
  # the bounding box.  Without a CRS, such a row took the layer out of the
  # lon/lat envelope, and the other rows were read as planar degrees.
  cases <- list(list("EPSG:4326", ll, c(1e200, 40.2)), list("EPSG:4326", ll, c(-84.2, 1e200)),
                list("no CRS", nc, c(-84.2, 1e200)), list("no CRS", nc, c(-1e200, 40.2)))
  for (cs in cases) {
    tag <- sprintf("%s, (%g, %g)", cs[[1L]], cs[[3L]][1L], cs[[3L]][2L])
    r <- run(.fm_set(cs[[2L]], cs[[3L]]))
    e <- run(.fm_set(cs[[2L]], NULL))
    expect_identical(r$value, e$value, info = tag)
    expect_identical(.fm_classes(r), .fm_classes(e), info = tag)
    expect_identical(.fm_texts(r), .fm_texts(e), info = tag)
    expect_false(is.na(attr(r$value, "pooled")$crs), info = tag)
    expect_identical(attr(r$value, "pooled")$columns$resp_z$skipped, 1L, info = tag)
  }
  # The CRS-less degrees are taken for lon/lat, with the warning that says so.
  r <- run(.fm_set(nc, c(-84.2, 1e200)))
  expect_true(any(grepl("assuming EPSG:4326", .fm_texts(r), fixed = TRUE)))
  expect_false(isTRUE(sf::st_is_longlat(attr(r$value, "pooled")$crs)))
})

test_that("a row of `points` beyond 1e150 as given in a layer without a CRS leaves the others read as lon/lat in cell_values_at()", {
  L <- .fm_layer()
  bare <- function(l) sf::st_set_crs(sf::st_transform(l, 4326), NA)
  pts <- bare(L$asg)
  nd  <- bare(L$held)
  for (seeds in list(bare(L$seeds), NULL)) {
    how <- if (is.null(seeds)) "polygons" else "seeds"
    for (at in list(c(-80, 1e200), c(1e200, 45))) {
      tag <- sprintf("%s, (%g, %g)", how, at[1L], at[2L])
      r <- .fm_catch(.fm_at(L, nd, points = .fm_set(pts, at), seeds = seeds))
      e <- .fm_catch(.fm_at(L, nd, points = .fm_set(pts, NULL), seeds = seeds))
      expect_identical(r$value, e$value, info = tag)
      expect_identical(.fm_texts(r), .fm_texts(e), info = tag)
      expect_true(any(grepl("look like lon/lat", .fm_texts(r), fixed = TRUE)), info = tag)
      expect_identical(attr(r$value, "values_at")$n, 479L, info = tag)
    }
  }
})

test_that("a far row of points with three coordinates is emptied in its own dimensions", {
  ll <- .fm_lonlat()
  xyz <- sf::st_as_sf(data.frame(sf::st_coordinates(ll), Z = 0, z = ll$z),
                      coords = c("X", "Y", "Z"), crs = 4326)
  g <- sf::st_geometry(xyz)
  far <- g; far[1L] <- sf::st_sfc(sf::st_point(c(1e200, 40.2, 0)), crs = 4326)
  emp <- g; emp[1L] <- sf::st_sfc(sf::st_point(rep(NA_real_, 3L)), crs = 4326)
  r <- .pooled_xy(sf::st_set_geometry(xyz, far))
  e <- .pooled_xy(sf::st_set_geometry(xyz, emp))
  # The emptied row comes out NaN and the far one NA: both are rows without a
  # location to every caller, which reads is.finite().
  expect_identical(r$crs, e$crs)
  expect_identical(r$xy[-1L, ], e$xy[-1L, ])
  expect_true(all(is.na(r$xy[1L, ])))
  expect_true(all(is.na(e$xy[1L, ])))
  expect_true(all(is.finite(r$xy[-1L, ])))
})

test_that("the warning's entry in the help names the bound", {
  # The source page from the source tree, the installed page under R CMD check,
  # which runs the tests where there is no man/.
  read_rd <- function(page) {
    file <- test_path("..", "..", "man", page)
    if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[[page]]
  }
  rd <- read_rd("summarize_by_cell.Rd")
  expect_s3_class(rd, "Rd")
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("rows with a value and a cell have no location: an empty point, ",
                           "a coordinate that is not finite or lies beyond 1e150, as given or ",
                           "in the CRS the distances are measured in, or in a geographic CRS a ",
                           "vertex at a latitude beyond 90 degrees."), fixed = TRUE)
  rd <- read_rd("cell_values_at.Rd")
  expect_s3_class(rd, "Rd")
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("not finite or lies beyond 1e150 as given or in the working CRS, ",
                           "has NA in every column"), fixed = TRUE)
})
