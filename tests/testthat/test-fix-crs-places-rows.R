# tests/testthat/test-fix-crs-places-rows.R
# ---------------------------------------------------------------------------
# BUG-193.  PROJ gives no coordinates for a point a UTM zone cannot place
# (about 90 degrees or more from its central meridian on the equator).  The
# CRS chooser scored a candidate on the points it could place and left the
# others out, so on LINESTRING (0 0, 180 0) it kept UTM zone 46, which holds
# neither end of the line.  ensure_projected() returned the line empty, and
# five points on the equator with two of them empty, without a word, and
# coerce_to_points(x, "auto") handed the empty line to sf::st_line_sample(),
# which segfaulted and ended the R session.
#
# A candidate is now used only if it places every located vertex of the
# layer; a sampled point it cannot place makes its score NA; a move out of
# lon/lat that leaves a located row empty, or without some of its vertices,
# warns with class spatialkit_rows_unplaced; and coerce_to_points() never
# samples a line its temporary projection did not place, but refuses it
# with class spatialkit_features_degenerate, site "projection".
#
# The cases that crashed run in a child R process, so that a relapse fails
# a test and does not end the test run.
# ---------------------------------------------------------------------------

.cpr_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
.cpr_line <- function()
  sf::st_sf(v = 1, geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(180, 0))),
                                         crs = 4326))
.cpr_equator <- function()
  sf::st_as_sf(data.frame(id = 1:5, x = seq(0, 180, by = 45), y = 0),
               coords = c("x", "y"), crs = 4326)

# Run `code` in a fresh R process that loads this same copy of spatialkit,
# and return its exit status, its output and what it saved as `res`.
.cpr_child <- function(code) {
  ns_path   <- getNamespaceInfo(asNamespace("spatialkit"), "path")
  installed <- file.exists(file.path(ns_path, "Meta", "package.rds"))
  load_line <- if (installed) {
    sprintf("suppressPackageStartupMessages(library(spatialkit, lib.loc = %s))",
            deparse(dirname(ns_path)))
  } else {
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))", deparse(ns_path))
  }
  script <- tempfile(fileext = ".R")
  out    <- tempfile(fileext = ".rds")
  on.exit(unlink(c(script, out)), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               load_line,
               "suppressMessages(sf::sf_use_s2(TRUE))",
               "q <- function(e) suppressWarnings(suppressMessages(e))",
               "try_ <- function(expr) tryCatch(expr, error = function(e) e)",
               "line <- sf::st_sf(v = 1, geometry = sf::st_sfc(",
               "  sf::st_linestring(rbind(c(0, 0), c(180, 0))), crs = 4326))",
               "res <- list()",
               code,
               sprintf("saveRDS(res, %s)", deparse(out))),
             script)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                     c("--vanilla", shQuote(script)),
                                     stdout = TRUE, stderr = TRUE, env = "R_TESTS="))
  status <- attr(output, "status")
  list(status = if (is.null(status)) 0L else as.integer(status),
       output = output,
       res    = if (file.exists(out)) readRDS(out) else NULL)
}
.cpr_expect_ran <- function(child) {
  expect_identical(child$status, 0L,
                   label = paste0("child R exit status; output:\n",
                                  paste(utils::tail(child$output, 15L), collapse = "\n")))
  expect_false(is.null(child$res))
}
.cpr_xy <- function(p) unname(sf::st_coordinates(p)[1L, 1:2])


test_that("ensure_projected() places LINESTRING (0 0, 180 0) in a CRS that holds both ends", {
  x <- .cpr_line()
  lines <- capture_spatialkit_log(p <- .cpr_quiet(ensure_projected(x)))
  expect_false(sf::st_is_empty(p)[1L])
  expect_true(all(is.finite(sf::st_coordinates(p)[, 1:2])))
  expect_identical(nrow(sf::st_coordinates(p)), 2L)
  expect_false(sf::st_crs(p) == sf::st_crs(32646))
  ch <- attr(p, "crs_choice")
  expect_identical(ch$name[1L], "UTM zone 46")
  expect_true(is.na(ch$distance_error[1L]))
  expect_false(ch$chosen[1L])
  expect_identical(sum(ch$chosen), 1L)
  expect_true(log_has(lines, "UTM zone 46 cannot place every location of the layer"))
  # The chooser alone gives the same.
  tr <- .cpr_quiet(.pick_local_projected_crs(x))
  expect_true(tr$crs == sf::st_crs(p))
  expect_identical(tr$candidates, ch)
})

test_that("five points on the equator from 0 to 180 come back with every row placed", {
  x <- .cpr_equator()
  expect_no_warning(p <- suppressMessages(ensure_projected(x)))
  expect_identical(sum(sf::st_is_empty(p)), 0L)
  expect_true(all(is.finite(sf::st_coordinates(p))))
  ch <- attr(p, "crs_choice")
  expect_true(is.na(ch$distance_error[ch$name == "UTM zone 46"]))
  expect_false(ch$chosen[ch$name == "UTM zone 46"])
  expect_false(sf::st_crs(p) == sf::st_crs(32646))
})

test_that("a point the CRS cannot place makes the distance score NA", {
  x <- .cpr_equator()
  expect_true(is.na(.crs_distance_error(x, sf::st_crs(32646))))
  # Without the two points the zone cannot place, it is scored.
  expect_true(is.finite(.crs_distance_error(x[2:4, ], sf::st_crs(32646))))
  expect_true(is.finite(.crs_distance_error(x, sf::st_crs(
    "+proj=laea +lat_0=0 +lon_0=90 +datum=WGS84 +units=m +no_defs"))))
})

test_that("a candidate that cannot place the layer is not used on the area path", {
  # The test of placement is replaced by one under which the Lambert
  # azimuthal candidate places nothing.
  sq <-sf::st_as_sf(data.frame(x = c(-110, -70, -110, -70), y = c(30, 30, 48, 48)),
                     coords = c("x", "y"), crs = 4326)
  base <- .cpr_quiet(.pick_local_projected_crs(sq, purpose = "area"))
  expect_identical(nrow(base$candidates), 2L)
  local_mocked_bindings(.places_vertices = function(xy, from, crs)
    !grepl("+proj=laea", crs$input, fixed = TRUE))
  lines <- capture_spatialkit_log(tr <- .cpr_quiet(.pick_local_projected_crs(sq, purpose = "area")))
  expect_true(grepl("+proj=aea", tr$crs$input, fixed = TRUE))
  expect_true(is.na(tr$candidates$distance_error[1L]))
  expect_identical(tr$candidates$chosen, c(FALSE, TRUE))
  expect_identical(tr$candidates$distance_error[2L], base$candidates$distance_error[2L])
  expect_true(log_has(lines, "Lambert azimuthal .* cannot place every location of the layer"))
})

test_that("with no candidate that places the layer, the global fallback is used", {
  sq <- sf::st_as_sf(data.frame(x = c(-110, -70, -110, -70), y = c(30, 30, 48, 48)),
                     coords = c("x", "y"), crs = 4326)
  local_mocked_bindings(.places_vertices = function(xy, from, crs) FALSE)
  for (p in c("distance", "area")) {
    lines <- capture_spatialkit_log(tr <- .cpr_quiet(.pick_local_projected_crs(sq, purpose = p)))
    expect_null(tr$candidates)
    if (p == "distance") expect_true(tr$crs == sf::st_crs(3857))
    else expect_true(grepl("eqearth", tr$crs$input, fixed = TRUE))
    expect_true(log_has(lines, "no candidate projection can place every location of the layer"))
  }
})

test_that("layers every candidate places choose and score as before", {
  # The result with the placement test switched off is the result of the
  # chooser before it: for layers every candidate places they are identical.
  set.seed(1)
  st39 <- sf::st_as_sf(data.frame(x = stats::runif(39, -120, -75), y = stats::runif(39, 30, 47)),
                       coords = c("x", "y"), crs = 4326)
  box <- function(x0, y0, w, h)
    sf::st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + h), c(x0, y0 + h),
                              c(x0, y0))))
  conus <- sf::st_sf(id = 1, geometry = sf::st_sfc(box(-122, 26, 52, 22), crs = 4326))
  nc <- sf::st_transform(sf::st_read(system.file("shape/nc.shp", package = "sf"),
                                     quiet = TRUE), 4326)
  fiji <- sf::st_sf(id = 1, geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(177, -19), c(-178, -19), c(-178, -15), c(177, -15), c(177, -19)))), crs = 4326))
  set.seed(2)
  ant <- sf::st_as_sf(data.frame(x = seq(-175, 175, by = 10), y = stats::runif(36, -80, -66)),
                      coords = c("x", "y"), crs = 4326)
  lat10 <- sf::st_as_sf(data.frame(x = seq(0, 180, by = 20), y = 10),
                        coords = c("x", "y"), crs = 4326)
  layers <- list(st39 = st39, conus = conus, nc = nc, fiji = fiji, antarctic = ant,
                 lat10 = lat10)
  run <- function() lapply(layers, function(x) lapply(c("distance", "area"), function(p) {
    tr <- .cpr_quiet(.pick_local_projected_crs(x, purpose = p))
    ep <- .cpr_quiet(ensure_projected(x, purpose = p))
    list(tr$crs, tr$candidates, sf::st_geometry(ep), attr(ep, "crs_choice"))
  }))
  now <- run()
  local_mocked_bindings(.places_vertices = function(xy, from, crs) TRUE)
  expect_identical(now, run())
  # Each of them has a candidate scored (the comparison or area path, or the
  # one projection on the other paths).
  expect_true(all(vapply(now, function(r) is.data.frame(r[[1L]][[2L]]), logical(1))))
})

test_that("the located vertices are the distinct ones, and one that shares a coordinate with another counts", {
  pts <- function(m) sf::st_as_sf(data.frame(x = m[, 1L], y = m[, 2L]),
                                  coords = c("x", "y"), crs = 4326)
  v <- .located_vertices(pts(rbind(c(0, 0), c(90, 0), c(180, 0), c(90, 10), c(0, 0))))
  expect_identical(v, rbind(c(0, 0), c(90, 0), c(90, 10), c(180, 0)))
  # Left out: a coordinate that is not finite, beyond 1e150, or a latitude
  # beyond 90, which no CRS places.
  v <- .located_vertices(pts(rbind(c(1, 2), c(Inf, 0), c(5, 1e200), c(-96, 95), c(40, -96))))
  expect_identical(v, rbind(c(1, 2)))
  # UTM zone 46 cannot place (180, 0), which shares its latitude with
  # (90, 0), the vertex sorted before it: the zone is not placed, scores NA
  # and is named.
  x <- pts(rbind(c(30, 5), c(90, 0), c(180, 0)))
  expect_false(.places_vertices(.located_vertices(x), sf::st_crs(4326), sf::st_crs(32646)))
  lines <- capture_spatialkit_log(tr <- .cpr_quiet(.pick_local_projected_crs(x)))
  expect_identical(tr$candidates$name[1L], "UTM zone 46")
  expect_true(is.na(tr$candidates$distance_error[1L]))
  expect_false(tr$candidates$chosen[1L])
  expect_true(log_has(lines, "UTM zone 46 cannot place every location of the layer"))
})

test_that("a row at a latitude beyond 90 does not send the other rows to the global fallback", {
  # No CRS places such a row, the global fallback included, so it takes no
  # part in the test of placement: the stations keep a Lambert azimuthal CRS
  # (0.88%) and do not get Web Mercator or Equal Earth.  The row comes back
  # empty, and the warning names it.
  set.seed(1)
  st39 <- sf::st_as_sf(data.frame(x = stats::runif(39, -120, -75), y = stats::runif(39, 30, 47)),
                       coords = c("x", "y"), crs = 4326)
  zone <- .cpr_quiet(.pick_local_projected_crs(st39))$candidates
  for (bad in list(c(-96, 95), c(40, -96))) {
    x <- rbind(st39, sf::st_as_sf(data.frame(x = bad[1L], y = bad[2L]),
                                  coords = c("x", "y"), crs = 4326))
    for (p in c("distance", "area")) {
      tag <- sprintf("(%g, %g), %s", bad[1L], bad[2L], p)
      tr <- .cpr_quiet(.pick_local_projected_crs(x, purpose = p))
      expect_false(is.null(tr$candidates), info = tag)
      expect_true(grepl("+proj=laea", tr$crs$input, fixed = TRUE), info = tag)
      expect_true(all(is.finite(tr$candidates$distance_error)), info = tag)
      if (p == "distance")
        expect_identical(tr$candidates$distance_error[1L], zone$distance_error[1L], info = tag)
      w <- NULL
      e <- withCallingHandlers(suppressMessages(ensure_projected(x, purpose = p)),
        warning = function(cnd) {
          if (inherits(cnd, "spatialkit_rows_unplaced")) w <<- cnd
          invokeRestart("muffleWarning")
        })
      expect_true(sf::st_crs(e) == tr$crs, info = tag)
      expect_identical(w$rows, 40L, info = tag)
      expect_identical(which(sf::st_is_empty(e)), 40L, info = tag)
    }
  }
})

test_that("the pooled readers report a row the move cannot place in their own words, without the warning", {
  x <- .cpr_equator()
  expect_no_warning(r <- suppressMessages(.pooled_xy(x, sf::st_crs(32646))))
  expect_identical(which(is.na(r$xy[, 1L])), c(1L, 5L))
  expect_true(all(is.finite(r$xy[2:4, ])))
  # The warning is given again after the call.
  expect_warning(suppressMessages(ensure_projected(x, target_crs = 32646)),
                 class = "spatialkit_rows_unplaced")
})

test_that("ensure_projected() warns when its move out of lon/lat leaves a located row unplaced", {
  x <- .cpr_equator()
  w <- NULL
  p <- withCallingHandlers(
    suppressMessages(ensure_projected(x, target_crs = 32646)),
    spatialkit_rows_unplaced = function(cnd) { w <<- cnd; invokeRestart("muffleWarning") })
  expect_s3_class(w, "spatialkit_rows_unplaced")
  expect_identical(w$rows, c(1L, 5L))
  expect_identical(w$n, 5L)
  expect_match(conditionMessage(w), "rows 1 and 5", fixed = TRUE)
  expect_match(conditionMessage(w), "EPSG:32646", fixed = TRUE)
  expect_identical(which(sf::st_is_empty(p)), c(1L, 5L))
  # The same from a CRS-less layer moved to the target, and from the CRS the
  # chooser returns, with and without a CRS (the chooser is replaced by one
  # that returns UTM zone 46, as it did before).
  x0 <- sf::st_set_crs(x, NA)
  expect_warning(expect_warning(suppressMessages(ensure_projected(x0, target_crs = 32646)),
                                class = "spatialkit_rows_unplaced",
                                regexp = "rows 1 and 5 of its 5"),
                 regexp = "look like lon/lat")
  local({
    local_mocked_bindings(.pick_local_projected_crs = function(x, purpose = "distance")
      list(crs = sf::st_crs(32646), candidates = NULL))
    expect_warning(suppressMessages(ensure_projected(x)),
                   class = "spatialkit_rows_unplaced", regexp = "rows 1 and 5 of its 5")
    expect_warning(expect_warning(suppressMessages(ensure_projected(x0)),
                                  class = "spatialkit_rows_unplaced",
                                  regexp = "rows 1 and 5 of its 5"),
                   regexp = "assuming EPSG:4326")
  })
  # A line that loses a vertex is named too: GDAL drops the vertex it cannot
  # move.
  l <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(80, 0), c(100, 0))),
    sf::st_linestring(rbind(c(10, 0), c(90, 0), c(170, 0))), crs = 4326))
  expect_warning(suppressMessages(ensure_projected(l, target_crs = 32646)),
                 class = "spatialkit_rows_unplaced", regexp = "row 2 of its 2")
  # A row without a location as given is emptied, as before, without it;
  # so is a move between two lon/lat CRSs, and a layer every row of which
  # is placed.
  z <- sf::st_as_sf(data.frame(x = c(9.1, Inf, 9.15), y = c(48.7, 48.75, 48.8)),
                    coords = c("x", "y"), crs = 4326)
  expect_no_warning(suppressMessages(ensure_projected(z)))
  expect_no_warning(suppressMessages(ensure_projected(x, target_crs = 4269)))
  expect_no_warning(suppressMessages(ensure_projected(x)))
  expect_no_warning(suppressMessages(ensure_projected(.cpr_line())))
})

test_that("coerce_to_points() takes the midpoint of LINESTRING (0 0, 180 0) and does not crash", {
  skip_on_cran()
  child <- .cpr_child(c(
    "res$auto <- q(coerce_to_points(line, 'auto'))",
    "res$mid  <- q(coerce_to_points(line, 'line_midpoint'))",
    "mls <- sf::st_sf(v = 1, geometry = sf::st_sfc(sf::st_multilinestring(",
    "  list(rbind(c(0, 0), c(180, 0)))), crs = 4326))",
    "res$mls  <- q(coerce_to_points(mls, 'auto'))",
    # A line with no location as given gets an empty point, as an empty line.
    "inf <- sf::st_sf(v = 1:2, geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(Inf, 0))),",
    "  sf::st_linestring(rbind(c(9, 48), c(9.2, 48.1))), crs = 4326))",
    "res$inf  <- q(coerce_to_points(inf, 'auto'))",
    "res$inf_mid <- q(coerce_to_points(inf, 'line_midpoint'))"))
  .cpr_expect_ran(child)
  r <- child$res
  for (nm in c("auto", "mid", "mls")) {
    expect_false(sf::st_is_empty(r[[nm]])[1L], label = nm)
    expect_lt(max(abs(.cpr_xy(r[[nm]]) - c(90, 0))), 1e-6)
  }
  for (nm in c("inf", "inf_mid")) {
    expect_identical(sf::st_is_empty(r[[nm]]), c(TRUE, FALSE), label = nm)
  }
})

test_that("prep_model_data() and make_folds() run on a layer of that line and two points", {
  skip_on_cran()
  child <- .cpr_child(c(
    "x <- sf::st_sf(z = c(1, 2, 3), a = c(0.5, 1.5, 1),",
    "  geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(180, 0))),",
    "    sf::st_point(c(30, 5)), sf::st_point(c(150, -5)), crs = 4326))",
    "res$prep <- q(prep_model_data(x, 'z', 'a'))",
    "res$folds <- q(make_folds(x, k = 3, method = 'random_kfold', seed = 1))"))
  .cpr_expect_ran(child)
  expect_identical(nrow(child$res$prep), 3L)
  expect_false(any(sf::st_is_empty(child$res$prep)))
  expect_true(all(is.finite(sf::st_coordinates(child$res$prep))))
  expect_false(is.null(child$res$folds))
})

test_that("a line the temporary projection does not place is refused, not sampled", {
  skip_on_cran()
  # ensure_projected() is replaced by a move into UTM zone 46, which empties
  # the line, as the chooser did before.
  child <- .cpr_child(c(
    "mock <- function(x, ...) sf::st_transform(x, 32646)",
    "testthat::local_mocked_bindings(ensure_projected = mock, .package = 'spatialkit')",
    "two <- sf::st_sf(z = c(1, 2, 3), a = c(0.5, 1.5, 1),",
    "  geometry = sf::st_sfc(sf::st_point(c(95, 5)), sf::st_linestring(rbind(c(0, 0), c(180, 0))),",
    "    sf::st_point(c(91, -5)), crs = 4326))",
    "mls <- sf::st_sf(v = 1, geometry = sf::st_sfc(sf::st_multilinestring(",
    "  list(rbind(c(0, 0), c(180, 0)))), crs = 4326))",
    "res$auto <- try_(q(coerce_to_points(line, 'auto')))",
    "res$mid  <- try_(q(coerce_to_points(line, 'line_midpoint')))",
    "res$mls  <- try_(q(coerce_to_points(mls, 'auto')))",
    "res$two  <- try_(q(coerce_to_points(two, 'auto')))",
    "res$prep <- try_(q(prep_model_data(two, 'z', 'a')))",
    "res$folds <- try_(q(make_folds(two, k = 3, method = 'random_kfold', seed = 1)))",
    # The remedies take the point without the projection.
    "res$cen <- q(coerce_to_points(line, 'centroid'))",
    "res$nop <- q(coerce_to_points(line, 'auto', tmp_project = FALSE))"))
  .cpr_expect_ran(child)
  r <- child$res
  for (nm in c("auto", "mid", "mls")) {
    e <- r[[nm]]
    expect_s3_class(e, "spatialkit_features_degenerate")
    expect_identical(e$site, "projection")
    expect_identical(e$rows, 1L)
    expect_identical(e$n, 1L)
    expect_identical(e$crs, "EPSG:32646")
    expect_match(conditionMessage(e), "coerce_to_points(): the temporary projection to EPSG:32646",
                 fixed = TRUE)
    expect_match(conditionMessage(e), "row 1 of the 1 features of `x`", fixed = TRUE)
    expect_match(conditionMessage(e), "mode = \"centroid\", or tmp_project = FALSE", fixed = TRUE)
  }
  expect_identical(r$two$rows, 2L)
  expect_match(conditionMessage(r$two), "row 2 of the 3 features of `x`", fixed = TRUE)
  expect_s3_class(r$prep, "spatialkit_features_degenerate")
  expect_identical(r$prep$site, "projection")
  expect_identical(r$prep$rows, 2L)
  expect_match(conditionMessage(r$prep), "prep_model_data(): the temporary projection",
               fixed = TRUE)
  expect_match(conditionMessage(r$prep), "row 2 of the 3 features of `data_sf`", fixed = TRUE)
  expect_match(conditionMessage(r$prep), "Pass pointize = \"centroid\"", fixed = TRUE)
  expect_s3_class(r$folds, "spatialkit_features_degenerate")
  expect_identical(r$folds$site, "projection")
  expect_match(conditionMessage(r$folds), "make_folds(): the temporary projection", fixed = TRUE)
  expect_match(conditionMessage(r$folds), "`points_sf`", fixed = TRUE)
  expect_match(conditionMessage(r$folds), "coerce_to_points(points_sf, \"centroid\")",
               fixed = TRUE)
  for (nm in c("cen", "nop")) {
    expect_lt(max(abs(.cpr_xy(r[[nm]]) - c(90, 0))), 1e-6)
  }
})

test_that(".pointize_refused() raises a projection refusal again in the caller's words", {
  e <- tryCatch(.raise_unplaced_lines(2L, 3L, "coerce_to_points", "`x`", "EPSG:32646",
                                      "Pass mode = \"centroid\"."),
                error = function(e) e)
  r <- tryCatch(.pointize_refused(e, "estimate_sac_range", arg = "`points_sf`", at = c(4L, 7L, 9L),
                                  n = 10L),
                error = function(e) e)
  expect_s3_class(r, "spatialkit_features_degenerate")
  expect_identical(r$site, "projection")
  expect_identical(r$rows, 7L)
  expect_identical(r$n, 10L)
  expect_identical(r$crs, "EPSG:32646")
  expect_match(conditionMessage(r), "estimate_sac_range(): the temporary projection to EPSG:32646",
               fixed = TRUE)
  expect_match(conditionMessage(r), "row 7 of the 10 features of `points_sf`", fixed = TRUE)
  # Raised again by a caller of that caller (select_features_forward()).
  r2 <- tryCatch(.pointize_refused(r, "select_features_forward", arg = "`train_sf`",
                                   at = 11:20, n = 20L),
                 error = function(e) e)
  expect_identical(r2$rows, 17L)
  expect_identical(r2$site, "projection")
  expect_match(conditionMessage(r2), "row 17 of the 20 features of `train_sf`", fixed = TRUE)
})

test_that("?ensure_projected and ?coerce_to_points state the placement rule", {
  rd <- function(topic) {
    path <- testthat::test_path("..", "..", "man", paste0(topic, ".Rd"))
    skip_if_not(file.exists(path), "the Rd sources are not here")
    paste(readLines(path), collapse = " ")
  }
  ep <- gsub("\\s+", " ", rd("ensure_projected"))
  expect_match(ep, "must place every location of the layer", fixed = TRUE)
  expect_match(ep, "spatialkit_rows_unplaced", fixed = TRUE)
  cp <- gsub("\\s+", " ", rd("coerce_to_points"))
  expect_match(cp, "site = \"projection\"", fixed = TRUE)
  expect_match(cp, "tmp_project = FALSE", fixed = TRUE)
})
