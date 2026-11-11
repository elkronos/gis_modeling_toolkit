# tests/testthat/test-fix-s2-pointize-callers.R
# ---------------------------------------------------------------------------
# The refusal of coerce_to_points(x, "auto") raised in the caller's words.
#
# In longitude and latitude with s2 on, coerce_to_points(x, "auto") takes the
# centroid of a GEOMETRYCOLLECTION with s2, and stops with the error of class
# `spatialkit_features_degenerate` when s2 refuses a polygon in it even
# without its repeated vertices (a bowtie, a figure eight).  The functions
# that reduce the caller's layer that way let the error through as it came,
# naming coerce_to_points(), an argument `x` and a `mode` the caller's
# function does not take (BUG-180).  Each now raises it in its own name,
# naming its argument and that argument's rows, with the class kept.  The
# same holds for the refusal of a centroid GEOS puts outside its feature in
# a projected CRS (review finding R-fb180-1), which area_of_applicability()
# and predict_surface() raise in their words too.
# ---------------------------------------------------------------------------

.pc_s2 <- function(on = TRUE, env = parent.frame()) {
  old <- suppressMessages(sf::sf_use_s2(on))
  withr::defer(suppressMessages(sf::sf_use_s2(old)), envir = env)
}
.pc_quiet <- function(expr) {
  r <- NULL
  capture_spatialkit_log(r <- suppressWarnings(suppressMessages(expr)))
  r
}
.pc_sq <- function(x0, y0, d = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0 + d), c(x0, y0 + d),
                            c(x0, y0))))
.pc_rep <- function(x0, y0, d = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0, y0 + d), c(x0, y0))))
.pc_bow <- function(x0, y0, d = 0.1)
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + d, y0 + d), c(x0 + d, y0), c(x0, y0 + d),
                            c(x0, y0))))
.pc_fig8 <- function(x0, y0, d = 0.1) {
  m <- d / 2
  sf::st_polygon(list(rbind(c(x0, y0), c(x0 + m, y0 + m), c(x0 + d, y0), c(x0 + d, y0 + d),
                            c(x0 + m, y0 + m), c(x0, y0 + d), c(x0, y0))))
}
# 64 squares of 0.1 degree at (-94, 39) in EPSG:4326 with a response, two
# predictors and the cell ID of each square; row 3 a collection holding the
# polygon `maker` draws, or the plain layer with maker = NULL.
.pc_layer <- function(maker = NULL) {
  withr::local_seed(1)
  g <- list(); k <- 0L
  for (j in 1:8) for (i in 1:8) {
    k <- k + 1L
    g[[k]] <- .pc_sq(-94 + (i - 1) * 0.1, 39 + (j - 1) * 0.1)
  }
  x <- sf::st_sf(v = stats::rnorm(k), x1 = stats::rnorm(k), x2 = stats::rnorm(k),
                 poly_id = seq_len(k), geometry = sf::st_sfc(g, crs = 4326))
  x$v <- x$v + x$x1
  if (!is.null(maker))
    sf::st_geometry(x)[3] <- sf::st_geometrycollection(list(maker(-93.8, 39)))
  x
}
# The points, cells, summary and folds the calls are given beside the layer.
.pc_kit <- function() {
  x <- .pc_layer()
  pts <- .pc_quiet(prep_model_data(x, "v", c("x1", "x2")))
  pts$poly_id <- x$poly_id
  list(x = x, pts = pts,
       vals = .pc_quiet(summarize_by_cell(pts, "v", cells_sf = x[, "poly_id"])),
       folds = .pc_quiet(make_folds(pts, k = 4, method = "random_kfold", seed = 1)))
}
# Each call on `bad` in the argument its error names.
.pc_calls <- function(bad, K) {
  rec <- local({
    pp <- .pc_quiet(ensure_projected(coerce_to_points(K$x, "auto")))
    prof <- .pc_quiet(resolution_profile(pp, "v", obs_per_cell = 4.5, nstart = 1))
    .pc_quiet(select_resolution(prof, criterion = "size"))$count_rule
  })
  with_rec <- bad
  attr(with_rec, "count_rule") <- rec
  incomplete <- bad
  incomplete$x2[1] <- NA
  sff <- function(layer, ...)
    select_features_forward(layer, "v", c("x1", "x2"), fit_fn = fit_rf_model, k = 3,
                            method = "random_kfold", ...)
  list(
    list("fold_separation", "data_sf", function() fold_separation(K$folds, bad)),
    list("fold_separation", "prediction_points",
         function() fold_separation(K$folds, K$pts, prediction_points = bad)),
    list("make_folds", "points_sf",
         function() make_folds(bad, k = 4, method = "random_kfold", seed = 1)),
    list("make_folds", "prediction_points",
         function() make_folds(K$pts, k = 4, method = "nndm", seed = 1,
                               prediction_points = bad)),
    list("make_folds", "prediction_points",
         function() make_folds(K$pts, k = 4, method = "knndm", seed = 1,
                               prediction_points = bad)),
    list("estimate_sac_range", "points_sf", function() estimate_sac_range(bad, "v")),
    list("resolution_profile", "data_sf",
         function() resolution_profile(bad, "v", "x1", levels = c(2, 4))),
    list("determine_optimal_levels", "data_sf",
         function() determine_optimal_levels(bad, max_levels = 4)),
    list("cell_values_at", "newdata",
         function() cell_values_at(bad, values = K$vals, points = K$pts, response_var = "v")),
    list("cell_values_at", "points",
         function() cell_values_at(K$pts, values = K$vals, points = bad, response_var = "v")),
    list("summarize_by_cell", "assigned_points_sf",
         function() summarize_by_cell(bad, "v", cells_sf = K$x[, "poly_id"],
                                      deff = "variogram")),
    list("summarize_by_cell", "assigned_points_sf",
         function() summarize_by_cell(bad, "v", cells_sf = K$x[, "poly_id"],
                                      estimate = "pooled")),
    list("summarize_by_cell", "assigned_points_sf",
         function() summarize_by_cell(with_rec, "v", cells_sf = K$x[, "poly_id"])),
    list("select_features_forward", "train_sf", function() sff(bad)),
    list("select_features_forward", "train_sf", function() sff(bad, select_on = "split")),
    list("select_features_forward", "train_sf", function() sff(incomplete)),
    list("select_features_forward", "train_sf",
         function() sff(incomplete, select_on = "split")))
}

test_that("each caller raises the refusal of a collection in its own name and argument", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("gstat")
  .pc_s2(TRUE)
  K <- .pc_kit()
  for (maker in list(.pc_bow, .pc_fig8)) {
    bad <- .pc_layer(maker)
    for (cl in .pc_calls(bad, K)) {
      fn <- cl[[1L]]; arg <- cl[[2L]]
      e <- tryCatch(.pc_quiet(cl[[3L]]()), error = function(e) e)
      label <- sprintf("%s() on `%s`", fn, arg)
      expect_s3_class(e, "spatialkit_features_degenerate")
      expect_identical(e$rows, 3L, label = label)
      expect_identical(e$n, 64L, label = label)
      msg <- conditionMessage(e)
      expect_true(startsWith(msg, paste0(fn, "(): s2, ")), label = label)
      expect_true(grepl(sprintf("refuses 1 of the 64 features of `%s` (row 3; ", arg), msg,
                        fixed = TRUE), label = label)
      expect_true(grepl(sprintf("Transform `%s` to a projected CRS", arg), msg, fixed = TRUE),
                  label = label)
      expect_true(grepl(sprintf("coerce_to_points(%s, \"point_on_surface\")", arg), msg,
                        fixed = TRUE), label = label)
      for (no in c("coerce_to_points() again", "`x`", "mode ="))
        expect_false(grepl(no, msg, fixed = TRUE), label = paste(label, no))
      # s2's own error is kept as the parent.
      expect_match(conditionMessage(e$parent), "Loop 0 is not valid")
    }
  }
})

test_that("the remedies the refusal names run, or refuse the bowtie by name", {
  skip_if_not_installed("gstat")
  .pc_s2(TRUE)
  K <- .pc_kit()
  runs <- function(layer) {
    f <- .pc_quiet(make_folds(layer, k = 4, method = "random_kfold", seed = 1))
    expect_length(f$folds, 4L)
    expect_no_error(.pc_quiet(fold_separation(f, layer)))
    expect_no_error(.pc_quiet(determine_optimal_levels(layer, max_levels = 4)))
    expect_no_error(.pc_quiet(cell_values_at(layer, values = K$vals, points = K$pts,
                                             response_var = "v")))
  }
  for (maker in list(.pc_bow, .pc_fig8)) {
    bad <- .pc_layer(maker)
    runs(.pc_quiet(coerce_to_points(bad, "point_on_surface")))
  }
  # Moved to a projected CRS, the figure eight gets GEOS's centroid, which
  # lies in it.  The bowtie's lobes then differ in area, and the centroid
  # GEOS takes lies far outside it: the call is refused, naming its row,
  # where it ran on that point without a word.  Repaired, it runs.
  runs(.pc_quiet(ensure_projected(.pc_layer(.pc_fig8))))
  moved <- .pc_quiet(ensure_projected(.pc_layer(.pc_bow)))
  e <- tryCatch(.pc_quiet(make_folds(moved, k = 4, method = "random_kfold", seed = 1)),
                error = function(e) e)
  expect_s3_class(e, "spatialkit_features_degenerate")
  expect_identical(e$site, "geos_centroid")
  expect_identical(e$rows, 3L)
  runs(.pc_quiet(sf::st_make_valid(moved)))
})

test_that("a failure of the inner folds other than the refusal keeps its plain wrap", {
  skip_if_not_installed("ranger")
  K <- .pc_kit()
  e <- tryCatch(.pc_quiet(select_features_forward(K$pts, "v", c("x1", "x2"),
                                                  fit_fn = fit_rf_model, k = 3,
                                                  block_size = -1)),
                error = function(e) e)
  expect_s3_class(e, "simpleError")
  expect_false(inherits(e, "spatialkit_features_degenerate"))
  expect_match(conditionMessage(e),
               "^select_features_forward\\(\\): could not build the inner CV folds, so no ")
  expect_match(conditionMessage(e), "make_folds(): `block_size` must be", fixed = TRUE)
})

test_that("coerce_to_points() and prep_model_data() keep their own words", {
  .pc_s2(TRUE)
  bad <- .pc_layer(.pc_bow)
  e <- tryCatch(.pc_quiet(coerce_to_points(bad, "auto")), error = function(e) e)
  expect_match(conditionMessage(e), "^coerce_to_points\\(\\): .* of `x` \\(row 3; ")
  expect_match(conditionMessage(e), "and call coerce_to_points() again", fixed = TRUE)
  expect_match(conditionMessage(e), "Or pass mode = \"point_on_surface\"", fixed = TRUE)
  for (pz in c("auto", "centroid")) {
    e <- tryCatch(.pc_quiet(prep_model_data(bad, "v", "x1", pointize = pz)),
                  error = function(e) e)
    expect_match(conditionMessage(e), "^prep_model_data\\(\\): .* of `data_sf` \\(row 3; ")
    # Row 3 is a collection, which pointize = "auto" takes the centroid of
    # and refuses again, so "auto" is not offered under "centroid" either
    # (BUG-190).
    expect_match(conditionMessage(e), "Or pass pointize = \"point_on_surface\", which takes",
                 fixed = TRUE)
    expect_false(grepl("pointize = \"auto\"", conditionMessage(e), fixed = TRUE))
    expect_false(grepl("reduce `data_sf` to points first", conditionMessage(e), fixed = TRUE))
  }
})

test_that("layers s2 accepts, and projected layers, go through the callers as before", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("gstat")
  .pc_s2(TRUE)
  K <- .pc_kit()
  layers <- list(plain = .pc_layer(),
                 repeated = .pc_layer(.pc_rep),
                 projected = .pc_quiet(ensure_projected(.pc_layer(.pc_fig8))))
  run_all <- function(bad) lapply(.pc_calls(bad, K), function(cl) {
    withr::local_seed(42)
    tryCatch(.pc_quiet(cl[[3L]]()), error = function(e) conditionMessage(e))
  })
  now <- lapply(layers, run_all)
  # The same calls with the layer reduced by coerce_to_points() as it was
  # before: the wrapper changes nothing where s2 accepts the layer.
  testthat::local_mocked_bindings(.pointize_auto = function(x, ...)
    coerce_to_points(x, "auto"))
  was <- lapply(layers, run_all)
  for (nm in names(layers)) {
    expect_identical(now[[nm]], was[[nm]], label = nm)
    # Every call ran.  (In the projected layer GEOS reads the figure eight,
    # which s2 refuses in degrees, and its centroid lies in it.)
    expect_false(any(vapply(now[[nm]], is.character, logical(1))), label = nm)
  }
})

# The calls of .pc_calls() and those of area_of_applicability() and
# predict_surface(), on `bad`, a layer in a projected CRS.
.pc_calls_projected <- function(bad, K) {
  pts <- .pc_quiet(ensure_projected(K$pts))
  fit <- .pc_quiet(fit_rf_model(pts, "v", c("x1", "x2"), num_trees = 20))
  fit_xy <- .pc_quiet(fit_rf_model(pts, "v", c("x1", "x2"), num_trees = 20,
                                   include_coords = TRUE))
  c(.pc_calls(bad, K), list(
    list("area_of_applicability", "newdata",
         function() area_of_applicability(bad, train_sf = pts,
                                           predictor_vars = c("x1", "x2"),
                                           coords = "yes")),
    list("area_of_applicability", "train_sf",
         function() area_of_applicability(pts, train_sf = bad,
                                           predictor_vars = c("x1", "x2"),
                                           coords = "yes")),
    list("area_of_applicability", "newdata",
         function() area_of_applicability(bad, model = fit_xy)),
    list("predict_surface", "grid", function() predict_surface(fit, grid = bad))))
}

test_that("a centroid GEOS puts outside its feature is refused in each caller's words", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("gstat")
  .pc_s2(TRUE)
  K <- .pc_kit()
  # The collection of row 3 holds the bowtie moved to the CRS
  # ensure_projected() picks: its lobes differ in area there, and GEOS's
  # centroid of it lies far outside it.
  bad <- .pc_quiet(ensure_projected(.pc_layer(.pc_bow)))
  for (cl in .pc_calls_projected(bad, K)) {
    fn <- cl[[1L]]; arg <- cl[[2L]]
    e <- tryCatch(.pc_quiet(cl[[3L]]()), error = function(e) e)
    label <- sprintf("%s() on `%s`", fn, arg)
    expect_s3_class(e, "spatialkit_features_degenerate")
    expect_identical(e$site, "geos_centroid", label = label)
    expect_identical(e$rows, 3L, label = label)
    expect_identical(e$n, 64L, label = label)
    msg <- conditionMessage(e)
    expect_true(startsWith(msg, sprintf(
      "%s(): the centroid GEOS takes of 1 of the 64 features of `%s` lies outside the feature's bounding box (row 3: \"Self-intersection[",
      fn, arg)), label = label)
    expect_true(grepl(sprintf("Reduce `%s` to points first, for example with coerce_to_points(%s, \"point_on_surface\")",
                              arg, arg), msg, fixed = TRUE), label = label)
    expect_true(grepl("repair the features with sf::st_make_valid(), and run the call again.",
                      msg, fixed = TRUE), label = label)
    for (no in c("coerce_to_points() again", "`x`", "mode =", "Transform", "s2"))
      expect_false(grepl(no, msg, fixed = TRUE), label = paste(label, no))
  }
})
