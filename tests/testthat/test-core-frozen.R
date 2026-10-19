# tests/testthat/test-core-frozen.R
# ---------------------------------------------------------------------------
# The cell pipeline returns what it returned on the baseline.
#
# fixtures/core-frozen.rds holds, for a fixed list of calls
# (helper-core-frozen.R), what each returned, printed, raised and logged on
# the baseline code and where it left the random-number stream, with the
# data the calls read and the platform they ran on.
# dev/baseline-core-frozen.R wrote it.  The calls are run again here, twice:
# as they are, and with FNN and Matrix switched off.
#
# On the platform of the fixture (the same R, sf, GEOS, PROJ, optional
# packages and linear algebra library) every record has to be identical().
# Anywhere else k-means differs in its last bits and Voronoi polygons and the
# text of a CRS change with GEOS, PROJ and sf, so names, classes, attribute
# names, the classes of the conditions and every text with its numbers taken
# out have to be identical, numbers equal to 1e-12 and coordinates equal to
# 1e-8 of the extent.  That second comparison is made on every platform:
# what is identical passes it.
#
# With the environment variable SPATIALKIT_FROZEN_STRICT set to "true" the
# identical() comparison fails where it would be skipped.  Set it wherever
# the fixture's platform is the one at hand, so that a changed package
# cannot switch the comparison off unseen.
# ---------------------------------------------------------------------------

.cf_baseline <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- readRDS(test_path("fixtures", "core-frozen.rds"))
    cache
  }
})

# The calls run once for the file.
.cf_replay <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- core_frozen_run(.cf_baseline()$data)
    cache
  }
})

# The calls once more, with FNN and Matrix switched off inside the package.
# FNN itself is made to stop when called, so a search that still reaches it
# shows as an error in the call's record.
.cf_replay_bare <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      if (requireNamespace("FNN", quietly = TRUE)) {
        reached <- function(...) stop("FNN was called with FNN switched off.", call. = FALSE)
        local_mocked_bindings(get.knn = reached, get.knnx = reached, .package = "FNN")
      }
      cache <<- core_frozen_run(.cf_baseline()$data,
                                without = .cf_baseline()$without$packages)
    }
    cache
  }
})

# What the replay without FNN and Matrix is held to: the records, but for the
# calls the baseline answers differently without the two.
.cf_baseline_bare <- function() {
  b <- .cf_baseline()
  rec <- b$records
  rec[names(b$without$records)] <- b$without$records
  rec
}

.cf_strict <- function()
  isTRUE(as.logical(Sys.getenv("SPATIALKIT_FROZEN_STRICT", unset = "false")))

# A comparison that cannot be made here: skipped, or failed when strict.
.cf_unmet <- function(why) {
  if (.cf_strict())
    fail(paste(why, "SPATIALKIT_FROZEN_STRICT is set, so the comparison has to run."))
  else skip(why)
}

# gstat decides what a call returns (without it a profile has no variogram),
# so without it there is nothing to compare.  FNN and Matrix do not: the
# records are the same without them, but for those the fixture keeps apart.
.cf_packages_missing <- function() {
  had  <- .cf_baseline()$info$packages[c("gstat", "dplyr", "s2")]
  have <- core_frozen_platform()$packages[names(had)]
  names(had)[!is.na(had) & is.na(have)]
}

# The parts of two runs that have to be identical on the baseline's platform.
.cf_expect_identical <- function(now, baseline) {
  for (nm in names(baseline)) {
    for (part in c("value", "conditions", "output", "printed", "log", "rng"))
      expect_identical(now[[nm]][[part]], baseline[[nm]][[part]],
                       label = sprintf("the %s of `%s`", part, nm),
                       expected.label = "the baseline's")
  }
}

.cf_expect_no_differences <- function(now, baseline) {
  expect_identical(names(now), names(baseline))
  same <- !nzchar(.cf_platform_differs())
  for (nm in names(baseline)) {
    d <- core_frozen_check(now[[nm]], baseline[[nm]], same)
    expect(length(d) == 0L,
           paste(c(sprintf("`%s` is not what it was on the baseline.", nm), d),
                 collapse = "\n"))
  }
}

# The platform at hand against the fixture's, in words; "" when they agree.
.cf_platform_differs <- function() {
  info <- .cf_baseline()$info
  here <- core_frozen_platform()
  if (core_frozen_same_platform(here, info)) return("")
  keys <- core_frozen_platform_keys()
  off  <- keys[!mapply(identical, here[keys], info[keys])]
  paste0("Not the platform the baseline was recorded on (", paste(off, collapse = ", "),
         " differ).")
}


test_that("the fixture names its baseline, its platform and the calls of the helper", {
  skip_on_cran()
  base <- .cf_baseline()
  info <- base$info
  expect_match(info$tree, "^[0-9a-f]{40}$")
  expect_identical(names(info$code), c("R", "DESCRIPTION", "NAMESPACE"))
  for (field in c("platform", "r", "sf", "geos", "proj"))
    expect_true(is.character(info[[field]]) && nzchar(info[[field]]), info = field)
  expect_named(info$packages, c("gstat", "FNN", "Matrix", "dplyr", "s2"))
  expect_named(base$data, c("pts", "clus", "big", "ll", "twice", "lattice"))
  # The calls are the ones the records were taken with, word for word.  A
  # call that names an argument it did not name then no longer pins that
  # argument's default.  The text is replaced only where the defaults
  # change, by dev/baseline-core-frozen.R --calls-only.
  expect_identical(names(base$records), names(core_frozen_calls()))
  expect_identical(core_frozen_call_text(), base$calls,
                   label = "the text of the helper's calls",
                   expected.label = "the text the fixture was recorded with")
  for (nm in names(base$records))
    expect_named(base$records[[nm]],
                 c("value", "conditions", "output", "printed", "log", "rng"), info = nm)
  # The records the baseline gives without FNN and Matrix, where they differ.
  expect_identical(base$without$packages, c("FNN", "Matrix"))
  expect_true(all(names(base$without$records) %in% names(base$records)))
  for (nm in names(base$without$records))
    expect_false(identical(base$without$records[[nm]], base$records[[nm]]), info = nm)
})

# On the fixture's platform every call is compared in full; elsewhere by what
# holds on any platform (core_frozen_portable_differences()).
test_that("every call returns the baseline's names, classes, numbers, geometries and texts", {
  skip_on_cran()
  gone <- .cf_packages_missing()
  if (length(gone)) {
    .cf_unmet(paste("The baseline was recorded with", paste(gone, collapse = ", "),
                    "installed."))
  } else {
    .cf_expect_no_differences(.cf_replay(), .cf_baseline()$records)
  }
})

test_that("on the baseline's platform every value, condition, print, log line and stream is identical", {
  skip_on_cran()
  gone <- .cf_packages_missing()
  off  <- .cf_platform_differs()
  if (length(gone)) {
    .cf_unmet(paste("The baseline was recorded with", paste(gone, collapse = ", "),
                    "installed."))
  } else if (nzchar(off)) {
    .cf_unmet(off)
  } else {
    .cf_expect_identical(.cf_replay(), .cf_baseline()$records)
  }
})

test_that("without FNN and Matrix every call returns what the baseline returns without them", {
  skip_on_cran()
  gone <- .cf_packages_missing()
  if (length(gone)) {
    .cf_unmet(paste("The baseline was recorded with", paste(gone, collapse = ", "),
                    "installed."))
  } else {
    now <- .cf_replay_bare()
    baseline <- .cf_baseline_bare()
    .cf_expect_no_differences(now, baseline)
    # No call reached FNN: its error would be in the record.
    reached <- vapply(now, function(r) any(vapply(r$conditions, function(c)
      grepl("FNN was called", c$message, fixed = TRUE), logical(1))), logical(1))
    expect_identical(names(now)[reached], character(0))
    # The identical() comparison where the platform allows it; its absence
    # elsewhere is reported by the test above.
    if (!nzchar(.cf_platform_differs())) .cf_expect_identical(now, baseline)
  }
})

test_that("the fixed list covers the calls the baseline is for", {
  skip_on_cran()
  rec <- .cf_baseline()$records
  cls <- function(nm) class(rec[[nm]]$value)[1L]
  raised <- function(nm) unlist(lapply(rec[[nm]]$conditions, function(c) c$class[1L]))
  # Eleven kinds of profile and three more, each read by every criterion.
  profiles <- grep("^profile_", names(rec), value = TRUE)
  expect_true(all(c("profile_default", "profile_geometry", "profile_levels", "profile_split",
                    "profile_missing", "profile_min_cell_n", "profile_n_levels",
                    "profile_density", "profile_sampled", "profile_logical",
                    "profile_no_value") %in% profiles))
  for (p in profiles) {
    expect_identical(cls(p), "resolution_profile", info = p)
    expect_false(is.null(rec[[p]]$printed), info = p)
    short <- sub("^profile_", "", p)
    expect_true(all(sprintf("select_%s_%s", short,
                            c("cp", "reliability", "elbow", "moran_z")) %in% names(rec)),
                info = p)
    expect_true(sprintf("summary_%s", short) %in% names(rec), info = p)
  }
  # What makes each profile its own case.
  b <- function(p) attr(rec[[p]]$value, "bounds")
  expect_true(is.na(attr(rec$profile_geometry$value, "variable")))
  expect_true(b("profile_levels")$levels_given)
  expect_false(is.null(attr(rec$profile_split$value, "split")))
  expect_lt(b("profile_missing")$n_response, b("profile_missing")$n)
  expect_identical(b("profile_min_cell_n")$min_cell_n, 5L)
  expect_identical(attr(rec$profile_density$value, "partition")$kind, "density")
  expect_null(attr(rec$profile_default$value, "partition"))
  expect_lt(b("profile_sampled")$n_sample, b("profile_sampled")$n)
  expect_identical(b("profile_no_value")$n_response, 0L)
  expect_true(b("profile_predictors")$predictors)
  # Every criterion returns a selection on some profile and stops on another.
  for (cr in c("cp", "reliability", "elbow", "moran_z_diagnostic")) {
    got <- vapply(grep(sprintf("^select_.*_%s$", cr), names(rec), value = TRUE), cls,
                  character(1))
    expect_true("resolution_selection" %in% got, info = cr)
    expect_true("NULL" %in% got, info = cr)
  }
  expect_identical(cls("select_column_subset"), "resolution_selection")
  expect_identical(rec$select_column_subset_without$conditions[[1L]]$kind, "error")
  expect_identical(cls("select_column_subset_cp"), "resolution_selection")
  # Seeds from a number, a selection and a profile, of both kinds; the other
  # two methods; cells; assignments.
  for (nm in c("seeds_number", "seeds_selection", "seeds_profile", "seeds_density",
               "seeds_random", "seeds_provided", "kmeans_number", "kmeans_selection",
               "kmeans_profile", "kmeans_density"))
    expect_identical(cls(nm), "sf", info = nm)
  expect_identical(attr(rec$seeds_number$value, "kmeans")$nstart, 10L)
  expect_identical(attr(rec$seeds_density$value, "kmeans")$nstart, 10L)
  expect_identical(nrow(rec$seeds_selection$value), rec$select_default_cp$value$best)
  expect_identical(nrow(rec$kmeans_profile$value), rec$select_default_cp$value$best)
  for (nm in grep("^cells_", names(rec), value = TRUE))
    expect_s3_class(rec[[nm]]$value$cells, "sf")
  expect_true(anyNA(rec$assigned_kept$value$cell_id))
  expect_lt(nrow(rec$assigned_dropped$value), nrow(rec$assigned_number$value))
  # The summaries.
  v <- function(nm) rec[[nm]]$value
  expect_false(inherits(v("values_default"), "sf"))
  expect_identical(attr(v("values_kish"), "deff_applied")$method, "kish")
  expect_true(all(v("values_sac")$deff_method == "variogram"))
  expect_equal(v("values_deff_2")[["..se_resp_z"]],
               sqrt(2) * v("values_default")[["..se_resp_z"]])
  expect_true(all(c("..ci_lo_resp_z", "..ci_hi_resp_z") %in% names(v("values_conf_level"))))
  expect_s3_class(v("values_cells_sf"), "sf")
  expect_true(all(c("cell_area", "n_per_area") %in% names(v("values_area"))))
  expect_true("resp_median_z" %in% names(v("values_median")))
  expect_true(any(grepl("x1", names(v("values_predictors")))) &&
                any(grepl("x2", names(v("values_predictors")))))
  expect_true("spatialkit_nonfinite" %in% raised("values_nonfinite"))
  expect_true(1L %in% v("values_sparse")$n)          # a one-row region
  expect_true(anyNA(v("values_sparse")$n))           # an empty cell
  expect_true(anyNA(v("values_unassigned")$cell_id))
  expect_identical(v("values_data_frame")$resp_mean_z, v("values_default")$resp_mean_z)
  # Conditions of every kind were recorded.
  kinds <- unlist(lapply(rec, function(r) vapply(r$conditions, `[[`, character(1), "kind")))
  expect_true(all(c("warning", "message", "error") %in% kinds))
  expect_gt(sum(lengths(lapply(rec, `[[`, "log"))), 0L)

  # ---- The calls added after the first list ----------------------------------
  # A count read off a profile: by its elbow when it is geometry-only, the
  # classed error when `cp` is missing beside another criterion, the plain
  # error when nothing is finite, and the two classed warnings.
  for (f in c("seeds", "kmeans", "hex")) {
    expect_false(is.null(v(sprintf("count_%s_geometry", f))), info = f)
    expect_identical(raised(sprintf("count_%s_no_cp", f)), "spatialkit_no_default_criterion",
                     info = f)
    expect_identical(raised(sprintf("count_%s_no_value", f)), "simpleError", info = f)
    expect_true("spatialkit_no_structure" %in% raised(sprintf("count_%s_noise", f)), info = f)
  }
  expect_equal(nrow(v("count_seeds_geometry")), v("select_geometry_elbow")$best)
  expect_true("spatialkit_bound_chosen" %in% raised("count_seeds_levels"))
  expect_identical(nrow(v("count_seeds_vector")), 12L)        # the first of c(12, 20, 8)
  expect_equal(nrow(v("count_seeds_candidates")), v("levels_default")[[1L]])
  expect_identical(nrow(v("count_hex_selection")$cells), nrow(v("cells_hex_profile")$cells))
  expect_true("spatialkit_partition_unused" %in% raised("count_seeds_selection_density_named"))
  # The fields a classed condition carries are part of its record.
  fields <- function(nm, class) Filter(function(c) c$class[1L] == class,
                                       rec[[nm]]$conditions)[[1L]]$fields
  expect_named(fields("count_seeds_levels", "spatialkit_bound_chosen"),
               c("edge", "criterion", "n"))
  expect_named(fields("count_seeds_noise", "spatialkit_no_structure"),
               c("criterion", "n", "one_cell"))
  expect_true(is.finite(fields("count_seeds_noise", "spatialkit_no_structure")$one_cell))
  # Where a call leaves the stream: only the calls that draw from the
  # caller's stream by design move it.
  moved <- names(rec)[!vapply(rec, function(r) r$rng$untouched, logical(1))]
  expect_setequal(moved, c("seeds_number", "seeds_density", "seeds_boundary", "seeds_random",
                           "kmeans_unseeded", "random_number"))
  for (nm in moved) expect_type(rec[[nm]]$rng$state, "integer")
  # Longitude and latitude, repeated locations, a lattice, a constant column.
  crs_of <- function(nm) sf::st_crs(v(nm))$input
  for (nm in c("seeds_lonlat", "seeds_lonlat_density", "kmeans_lonlat", "kmeans_lonlat_density",
               "assigned_lonlat"))
    expect_identical(crs_of(nm), "EPSG:4326", info = nm)
  # The cells of such seeds are drawn in a projected CRS, and the summary
  # joined to them is in it.
  expect_identical(sf::st_crs(v("cells_lonlat")$cells)$input, "EPSG:32632")
  expect_identical(crs_of("values_lonlat"), "EPSG:32632")
  expect_identical(nrow(v("seeds_more_than_locations")), 150L)
  expect_identical(nrow(v("kmeans_more_than_locations")), 150L)
  expect_match(rec$seeds_more_than_locations$conditions[[1L]]$message, "clamping")
  expect_true("spatialkit_constant_response" %in% raised("profile_constant"))
  expect_identical(cls("profile_lattice"), "resolution_profile")
  # Regions the caller brings, with an ID of each class.
  expect_type(v("values_units_character")$GEOID, "character")
  expect_type(v("values_units_integer")$id_int, "integer")
  expect_type(v("values_units_double")$id_dbl, "double")
  expect_s3_class(v("values_units_factor")$id_fac, "factor")
  expect_identical(sum(!is.na(v("values_id_number_text")$n)), 30L)
  # The summary's other present paths.
  expect_identical(v("values_id_found")$resp_mean_z, v("values_default")$resp_mean_z)
  expect_identical(vapply(rec$values_not_quiet$conditions, `[[`, character(1), "kind"),
                   "message")
  expect_match(rec$values_other_cells$conditions[[1L]]$message, "match no `cells_sf")
  a <- v("values_area_missing")
  expect_false(identical(a$n, a$cell_weight))        # rows against rows with a value
  expect_equal(sum(a$n_per_area * a$cell_area), 400)
  expect_identical(nrow(v("values_one_region")), 1L)
  expect_true(all(v("values_one_row_regions")$n == 1L))
  expect_equal(sum(v("values_one_row")$n, na.rm = TRUE), 1)
  expect_true(all(is.na(v("values_no_rows")$n)))
})

test_that("the comparison for another platform passes the last bits and nothing more", {
  skip_on_cran()
  rec <- .cf_baseline()$records
  same <- function(a, b) length(core_frozen_differences(a, b)) == 0L
  bare <- function(value) list(value = value, conditions = list(), output = character(0),
                               printed = NULL, log = character(0),
                               rng = list(untouched = TRUE, state = NULL))

  # A number: a relative 1e-14 passes, 1e-9 does not.
  prof <- rec$profile_default
  near <- prof; near$value$wss <- near$value$wss * (1 + 1e-14)
  far  <- prof; far$value$wss[3L] <- far$value$wss[3L] * (1 + 1e-9)
  expect_true(same(prof, prof))
  expect_true(same(near, prof))
  expect_false(same(far, prof))
  # The tolerance is 1e-12: every element moved by a relative 2e-13 passes,
  # by 5e-12 it does not.
  near <- prof; near$value$wss <- near$value$wss * (1 + 2e-13)
  far  <- prof; far$value$wss  <- far$value$wss * (1 + 5e-12)
  expect_true(same(near, prof))
  expect_false(same(far, prof))
  # A number inside an attribute, and an integer.
  moved <- prof; attr(moved$value, "bounds")$area <- attr(moved$value, "bounds")$area * 1.001
  expect_false(same(moved, prof))
  moved <- prof; moved$value$levels[2L] <- moved$value$levels[2L] + 1L
  expect_false(same(moved, prof))
  # A missing value where there was a number.
  moved <- prof; moved$value$cp[4L] <- NA
  expect_false(same(moved, prof))
  # Names, classes, attribute names and values that are not numbers.
  moved <- prof; names(moved$value)[2L] <- "WSS"
  expect_false(same(moved, prof))
  moved <- prof; class(moved$value) <- "data.frame"
  expect_false(same(moved, prof))
  moved <- prof; attr(moved$value, "extra") <- 1
  expect_false(same(moved, prof))
  moved <- prof; attr(moved$value, "variable") <- "residuals"
  expect_false(same(moved, prof))
  moved <- prof; moved$value$levels <- as.numeric(moved$value$levels)
  expect_false(same(moved, prof))
  # The levels of a factor.
  fac <- rec$values_units_factor
  expect_s3_class(fac$value$id_fac, "factor")
  moved <- fac; levels(moved$value$id_fac)[3L] <- "other"
  expect_false(same(moved, fac))

  # A geometry: the tolerance is 1e-8 of the extent.  A vertex moved by
  # 1e-10 or 2e-9 of it passes, by 5e-8 or 1e-6 it does not, and neither
  # does a cell more or fewer.
  cells <- rec$cells_number
  shift <- function(r, by, cell = 1L, vertex = 2L) {
    g <- sf::st_geometry(r$value$cells)
    ring <- g[[cell]][[1L]]
    ring[vertex, 1L] <- ring[vertex, 1L] + by
    if (vertex == 1L) ring[nrow(ring), 1L] <- ring[1L, 1L]
    g[[cell]] <- sf::st_polygon(list(ring))
    sf::st_geometry(r$value$cells) <- g
    r
  }
  extent_of <- function(r) {
    bb <- sf::st_bbox(r$value$cells)
    max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  }
  extent <- extent_of(cells)
  expect_true(same(shift(cells, 1e-10 * extent), cells))
  expect_true(same(shift(cells, 2e-9 * extent), cells))
  expect_false(same(shift(cells, 5e-8 * extent), cells))
  expect_false(same(shift(cells, 1e-6 * extent), cells))
  fewer <- cells; fewer$value$cells <- fewer$value$cells[-1L, ]
  expect_false(same(fewer, cells))
  # A ring that starts at another vertex, or runs the other way, is the same
  # polygon.
  turned <- cells
  g <- sf::st_geometry(turned$value$cells)
  ring <- g[[1L]][[1L]]
  open <- ring[-nrow(ring), , drop = FALSE]
  open <- open[c(3:nrow(open), 1:2), , drop = FALSE]
  g[[1L]] <- sf::st_polygon(list(rbind(open, open[1L, ])))
  sf::st_geometry(turned$value$cells) <- g
  expect_true(same(turned, cells))
  g[[1L]] <- sf::st_polygon(list(ring[nrow(ring):1, , drop = FALSE]))
  sf::st_geometry(turned$value$cells) <- g
  expect_true(same(turned, cells))
  # Both at once, with every vertex moved inside the tolerance.
  nudged <- ring[nrow(ring):1, , drop = FALSE]
  nudged <- nudged[-nrow(nudged), , drop = FALSE]
  nudged <- nudged[c(2:nrow(nudged), 1L), , drop = FALSE] + 1e-10 * extent
  g[[1L]] <- sf::st_polygon(list(rbind(nudged, nudged[1L, ])))
  sf::st_geometry(turned$value$cells) <- g
  expect_true(same(turned, cells))
  # Two vertices that share the smallest x, as every cell of a hexagon grid
  # has and every cell cut by a rectangular boundary: a last bit on one of
  # them must not put the two rings out of step.
  hex <- rec$cells_hex_profile
  ring <- sf::st_geometry(hex$value$cells)[[1L]][[1L]]
  open <- ring[-nrow(ring), , drop = FALSE]
  tied <- which(open[, 1L] == min(open[, 1L]))
  expect_gte(length(tied), 2L)
  lowest <- tied[which.min(open[tied, 2L])]
  hex_extent <- extent_of(hex)
  expect_true(same(shift(hex, 1e-10 * hex_extent, vertex = lowest), hex))
  expect_true(same(shift(hex, -1e-10 * hex_extent, vertex = lowest), hex))
  expect_false(same(shift(hex, 1e-6 * hex_extent, vertex = lowest), hex))
  # The text of a CRS is not compared, its name is.
  other <- cells
  crs <- sf::st_crs(other$value$cells)
  crs$wkt <- paste(crs$wkt, " ")
  attr(sf::st_geometry(other$value$cells), "crs") <- crs
  expect_true(same(other, cells))
  crs$input <- "EPSG:32633"
  attr(sf::st_geometry(other$value$cells), "crs") <- crs
  expect_false(same(other, cells))
  # A CRS that is not on a geometry column: the one a profile keeps for its
  # centres.
  expect_s3_class(attr(prof$value, "centres_crs"), "crs")
  moved <- prof
  crs <- attr(moved$value, "centres_crs")
  crs$wkt <- paste(crs$wkt, " ")
  attr(moved$value, "centres_crs") <- crs
  expect_true(same(moved, prof))
  crs$input <- "EPSG:32633"
  attr(moved$value, "centres_crs") <- crs
  expect_false(same(moved, prof))

  # Empty geometries are compared as empty.
  holes <- bare(sf::st_sfc(sf::st_point(), sf::st_polygon(), sf::st_point(c(1, 2))))
  expect_true(same(holes, holes))
  filled <- holes
  filled$value[[1L]] <- sf::st_point(c(0, 0))
  expect_false(same(filled, holes))

  # Conditions: when each was raised, its kind and its class count, and so
  # does its text, but for the numbers in it.
  warned <- rec$profile_missing
  expect_gt(length(warned$conditions), 0L)
  expect_match(warned$conditions[[1L]]$message, "[0-9]")
  renumbered <- warned
  renumbered$conditions[[1L]]$message <- gsub("[0-9]+", "7", warned$conditions[[1L]]$message)
  expect_true(same(renumbered, warned))
  reworded <- warned; reworded$conditions[[1L]]$message <- "other words"
  expect_false(same(reworded, warned))
  reclassed <- warned; reclassed$conditions[[1L]]$class <- c("simpleWarning", "warning", "condition")
  expect_false(same(reclassed, warned))
  silent <- warned; silent$conditions <- list()
  expect_false(same(silent, warned))
  later <- warned; later$conditions[[1L]]$when <- "print"
  expect_false(same(later, warned))
  # The fields of a classed condition: their names, and their numbers to the
  # tolerance.
  bound <- rec$count_seeds_levels
  i <- which(vapply(bound$conditions, function(c) c$class[1L], character(1)) ==
               "spatialkit_bound_chosen")
  expect_length(i, 1L)
  moved <- bound; moved$conditions[[i]]$fields$n <- moved$conditions[[i]]$fields$n + 1L
  expect_false(same(moved, bound))
  moved <- bound; moved$conditions[[i]]$fields$edge <- "another edge"
  expect_false(same(moved, bound))
  moved <- bound; moved$conditions[[i]]$fields$n <- NULL
  expect_false(same(moved, bound))
  flat <- rec$count_seeds_noise
  i <- which(vapply(flat$conditions, function(c) c$class[1L], character(1)) ==
               "spatialkit_no_structure")
  near <- flat; near$conditions[[i]]$fields$one_cell <- near$conditions[[i]]$fields$one_cell * (1 + 1e-14)
  far  <- flat; far$conditions[[i]]$fields$one_cell  <- far$conditions[[i]]$fields$one_cell * (1 + 1e-9)
  expect_true(same(near, flat))
  expect_false(same(far, flat))
  # A call that stopped against one that returned.
  stopped <- rec$select_default_elbow
  expect_null(stopped$value)
  expect_false(same(stopped, rec$select_default_cp))

  # The lines logged, written and printed: their number and their text but
  # for the numbers.
  logged <- rec$count_seeds_noise
  expect_gt(length(logged$log), 0L)
  moved <- logged; moved$log <- c(moved$log, "WARN one more caution")
  expect_false(same(moved, logged))
  moved <- logged; moved$log <- character(0)
  expect_false(same(moved, logged))
  moved <- logged; moved$log[1L] <- sub("^WARN", "ERROR", moved$log[1L])
  expect_false(same(moved, logged))
  moved <- logged; moved$log[1L] <- gsub("[0-9]+", "7", moved$log[1L])
  expect_true(same(moved, logged))
  moved <- prof; moved$output <- "written by the call"
  expect_false(same(moved, prof))
  moved <- prof; moved$printed <- c(moved$printed, "  estimate    : pooled")
  expect_false(same(moved, prof))
  moved <- prof; moved$printed[1L] <- "Another heading"
  expect_false(same(moved, prof))
  moved <- prof; moved$printed <- gsub("[0-9]+", "7", moved$printed)
  expect_true(same(moved, prof))
  moved <- prof; moved$printed <- NULL
  expect_false(same(moved, prof))
  expect_identical(core_frozen_mask(c("a  1,234.5e-3 b", "x = NaN, -2", " NA ")),
                   c("a # b", "x = NA, #", "NA"))
  expect_null(core_frozen_mask(NULL))

  # The stream: left as it was, or moved.
  drew <- rec$seeds_number
  expect_false(drew$rng$untouched)
  moved <- drew; moved$rng <- list(untouched = TRUE, state = NULL)
  expect_false(same(moved, drew))
  moved <- prof; moved$rng <- drew$rng
  expect_false(same(moved, prof))
  # How far a drawing call moves it is held by the identical() comparison only.
  moved <- drew; moved$rng$state[5L] <- moved$rng$state[5L] + 1L
  expect_true(same(moved, drew))
})

test_that("two platform records are the same platform only when every key agrees", {
  here <- core_frozen_platform()
  expect_true(core_frozen_same_platform(here, here))
  expect_identical(core_frozen_platform_keys(),
                   c("platform", "r", "sf", "geos", "proj", "packages", "blas", "lapack"))
  for (key in core_frozen_platform_keys()) {
    other <- here
    other[[key]] <- if (key == "packages") replace(here$packages, "s2", "0.0.0") else "another"
    expect_false(core_frozen_same_platform(other, here), info = key)
  }
  for (pkg in names(here$packages)) {
    other <- here
    other$packages[[pkg]] <- NA_character_
    expect_false(core_frozen_same_platform(other, here), info = pkg)
  }
  # GDAL reads and writes files; no call of the list does.
  other <- here; other$gdal <- "another"
  expect_true(core_frozen_same_platform(other, here))

  # The switch that turns a skipped comparison into a failed one.
  withr::with_envvar(c(SPATIALKIT_FROZEN_STRICT = "true"), {
    expect_true(.cf_strict())
    expect_failure(.cf_unmet("Not the platform."), "SPATIALKIT_FROZEN_STRICT")
  })
  withr::with_envvar(c(SPATIALKIT_FROZEN_STRICT = NA), {
    expect_false(.cf_strict())
    expect_condition(.cf_unmet("Not the platform."), class = "skip")
  })
})

test_that("the recorder keeps what a call returns, prints, raises and logs, and puts the session back", {
  skip_on_cran()
  calls <- list(
    warns   = quote({ warning("w1"); 1 }),
    says    = quote({ message("m1"); 2 }),
    stops   = quote(stop("e1")),
    writes  = quote({ cat("out\n"); 3 }),
    logs    = quote({ .log_warn("caution %d", 5L); .log_info("a note"); 4 }),
    once    = quote({ .log_warn_once("core-frozen-test", "said once"); 5 }),
    again   = quote({ .log_warn_once("core-frozen-test", "said once"); warns + 5 }),
    draws   = quote(stats::runif(2)),
    redraws = quote(stats::runif(2)),
    seeded  = quote({ cleanup <- .with_seed(4L); x <- stats::runif(1); cleanup(); x }),
    classed = quote(.warn_classed("spatialkit_core_frozen_test", "classed %s", "text",
                                  data = list(n = 3L, where = "here"))),
    inside  = quote((function() warning("from a function"))()),
    quiet   = quote(warning("no call", call. = FALSE)),
    shown   = quote(structure(list(best = 7L), class = "resolution_selection")))
  data <- .cf_baseline()$data
  once <- get(".warned_once", envir = asNamespace("spatialkit"))
  assign("core-frozen-sentinel", TRUE, envir = once)
  withr::defer(rm(list = intersect(c("core-frozen-sentinel", "core-frozen-test"),
                                   ls(once, all.names = TRUE)), envir = once))
  before <- list(options = options()[c("width", "digits", "scipen", "OutDec",
                                       "useFancyQuotes", "warn")],
                 collate = Sys.getlocale("LC_COLLATE"),
                 language = Sys.getenv("LANGUAGE", unset = NA),
                 once = sort(ls(once, all.names = TRUE)),
                 threshold = logger::log_threshold(namespace = "spatialkit", index = 2),
                 appender = logger::log_appender(namespace = "spatialkit", index = 2),
                 kind = RNGkind())
  withr::local_seed(99)
  expected_next <- stats::runif(1)
  set.seed(99)
  rec <- core_frozen_run(data, calls)
  expect_identical(stats::runif(1), expected_next)      # the stream was not touched

  after <- list(options = options()[c("width", "digits", "scipen", "OutDec",
                                      "useFancyQuotes", "warn")],
                collate = Sys.getlocale("LC_COLLATE"),
                language = Sys.getenv("LANGUAGE", unset = NA),
                once = sort(ls(once, all.names = TRUE)),
                threshold = logger::log_threshold(namespace = "spatialkit", index = 2),
                appender = logger::log_appender(namespace = "spatialkit", index = 2),
                kind = RNGkind())
  expect_identical(after, before)

  cond <- function(nm) rec[[nm]]$conditions
  for (nm in names(rec))
    expect_named(rec[[nm]], c("value", "conditions", "output", "printed", "log", "rng"),
                 info = nm)
  expect_identical(rec$warns$value, 1)
  expect_identical(cond("warns")[[1L]][c("when", "kind", "class", "message", "fields")],
                   list(when = "call", kind = "warning",
                        class = c("simpleWarning", "warning", "condition"),
                        message = "w1", fields = structure(list(), names = character(0))))
  expect_identical(cond("says")[[1L]][c("kind", "message")],
                   list(kind = "message", message = "m1\n"))
  expect_null(rec$stops$value)
  expect_identical(cond("stops")[[1L]][c("kind", "message")],
                   list(kind = "error", message = "e1"))
  expect_identical(rec$writes$output, "out")
  expect_identical(rec$writes$value, 3)
  # WARN lines are kept without their time stamp; INFO lines are not kept.
  expect_identical(rec$logs$log, "WARN caution 5")
  # What is logged once per session is logged by every call that reaches it.
  expect_length(rec$once$log, 1L)
  expect_identical(rec$again$log, rec$once$log)
  # A call reads the results before it; every call starts from the same seed.
  expect_identical(rec$again$value, 6)
  expect_identical(rec$redraws$value, rec$draws$value)
  # Where each call left the stream: a call that draws moves it, and to the
  # same state each time; one that draws under its own seed and puts the
  # stream back does not.
  expect_true(rec$warns$rng$untouched)
  expect_null(rec$warns$rng$state)
  expect_false(rec$draws$rng$untouched)
  expect_type(rec$draws$rng$state, "integer")
  expect_identical(rec$redraws$rng, rec$draws$rng)
  expect_true(rec$seeded$rng$untouched)
  # A classed warning keeps its class and its fields, and is logged as well
  # as raised.
  expect_identical(cond("classed")[[1L]]$class[1L], "spatialkit_core_frozen_test")
  expect_identical(cond("classed")[[1L]]$fields, list(n = 3L, where = "here"))
  expect_identical(rec$classed$log, "WARN classed text")
  # The function a condition's call names, NA when it carries none.
  expect_identical(cond("classed")[[1L]]$call, NA_character_)
  expect_identical(cond("quiet")[[1L]]$call, NA_character_)
  expect_match(cond("inside")[[1L]]$call, "^\\(?function")
  expect_identical(cond("says")[[1L]]$call, "message")
  # A selection is printed; a number is not.
  expect_false(is.null(rec$shown$printed))
  expect_null(rec$draws$printed)
})

test_that("switching FNN and Matrix off reaches every function that asks for them", {
  skip_on_cran()
  off <- core_frozen_without(c("FNN", "Matrix"))
  expect_true(all(c(".partition_centres", ".nearest_centre", ".knn_locations", ".sac_knn",
                    ".build_knn_weights") %in% names(off)))
  asks <- function(f) any(grepl("requireNamespace\\(\"(FNN|Matrix)\"", deparse(f)))
  ns <- asNamespace("spatialkit")
  for (nm in names(off)) {
    if (nm == ".has_backend") next
    expect_true(asks(get(nm, envir = ns)), info = nm)
    expect_false(asks(off[[nm]]), info = nm)
    expect_identical(environment(off[[nm]]), environment(get(nm, envir = ns)), info = nm)
  }
  # No function that asks is left out.
  every <- Filter(function(nm) { f <- get(nm, envir = ns); is.function(f) && asks(f) },
                  ls(ns, all.names = TRUE))
  expect_setequal(setdiff(names(off), ".has_backend"), every)
  expect_identical(formals(off$.nearest_centre)$use_fnn, FALSE)
  expect_identical(formals(off$.build_knn_weights)$use_matrix, FALSE)
  # Only the packages named: with FNN alone, Matrix is still asked for.
  fnn_only <- core_frozen_without("FNN")
  expect_identical(formals(fnn_only$.build_knn_weights)$use_fnn, FALSE)
  expect_true(is.call(formals(fnn_only$.build_knn_weights)$use_matrix))
  # A run without them leaves the package's own functions in place.
  before <- get(".sac_knn", envir = ns)
  rec <- core_frozen_run(.cf_baseline()$data,
                         list(asked = quote(deparse(.sac_knn)),
                              three = quote(.sac_knn(cbind(c(0, 1, 5), c(0, 0, 0)), 1L)$index[, 1L])),
                         without = c("FNN", "Matrix"))
  expect_false(any(grepl("requireNamespace", rec$asked$value)))
  expect_identical(rec$three$value, c(2L, 1L, 2L))
  expect_identical(get(".sac_knn", envir = ns), before)
})
