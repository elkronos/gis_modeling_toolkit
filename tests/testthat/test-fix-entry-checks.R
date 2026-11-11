# tests/testthat/test-fix-entry-checks.R
# ---------------------------------------------------------------------------
# Inputs whose parts contradict each other are refused, or made consistent,
# where they are read, instead of giving a wrong result further on.
#
# - cv_spatial() scored a fold's predictions against `response_var` whatever
#   response the fold's fit modelled: a fit_fn that fits price, called with
#   response_var = "slope", was scored as a model of slope on every fold.
# - kriging_adequacy() kriges each location once, so folds whose test sets
#   shared a location kept the last fold's prediction there, and the
#   statistic moved with the order of the folds.
# - kriging_adequacy() read a vector of fold labels against the rows kept,
#   so one label per row of `assigned_points_sf` was refused as soon as a
#   row was dropped for a missing response or no location, with a message
#   that said the labels must have one entry per row.
# - create_voronoi_polygons() and build_tessellation() handed on a seeds
#   layer's record of its count unchanged when the seeds were fewer or more
#   than its count, so the cells claimed to be the cells the profile scored and
#   summarize_by_cell() warned with causes that did not happen.
# ---------------------------------------------------------------------------

.ec_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

# --- cv_spatial(): a fold whose fit models another response ----------------

.ec_site <- function(n = 80) {
  set.seed(1)
  site <- sf::st_as_sf(data.frame(x = 5e5 + stats::runif(n, 0, 1000),
                                  y = 5e6 + stats::runif(n, 0, 1000),
                                  elev = stats::rnorm(n)),
                       coords = c("x", "y"), crs = 32632)
  site$price <- 10 + 0.01 * (sf::st_coordinates(site)[, 1] - 5e5) +
    2 * site$elev + stats::rnorm(n)
  site$slope <- stats::rnorm(n)
  site
}

test_that("cv_spatial() fails a fold whose fit models another response, naming both", {
  site <- .ec_site()
  registerS3method("predict", "ec_lm_fit", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))
  })
  # The ?new_spatial_fit example's fit_fn: it names its response itself.
  lm_fit <- function(train_sf)
    new_spatial_fit("ec_lm_fit", stats::lm(price ~ elev, sf::st_drop_geometry(train_sf)),
                    price ~ elev, "price", "elev", train_sf)
  ok <- .ec_quiet(cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1))
  expect_identical(as.integer(ok$n_folds_succeeded), 3L)

  # Called for slope, its predictions of price were scored against slope and
  # every fold was `ok`.
  msg <- paste0("fit_fn() returned a fit of 'price' (its `response_var`), ",
                "where the call's `response_var` is 'slope'")
  w <- NULL
  cv <- logger::with_log_threshold(
    withCallingHandlers(
      suppressMessages(cv_spatial(site, "slope", "elev", fit_fn = lm_fit, k = 3, seed = 1)),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }),
    threshold = logger::FATAL, namespace = "spatialkit", index = 2)
  expect_identical(as.integer(cv$n_folds_succeeded), 0L)
  expect_identical(cv$fold_status$status, rep("error", 3L))
  expect_identical(cv$fold_status$message, rep(msg, 3L))
  expect_true(all(is.na(cv$overall$RMSE)))
  # The all-folds-failed warning gives the cause.
  expect_true(any(grepl("all folds failed", w, fixed = TRUE) &
                  grepl(paste0("First error: ", msg), w, fixed = TRUE)))
})

test_that("cv_spatial() fails a fold whose fit has no response_var, and says so", {
  site <- .ec_site(40)
  registerS3method("predict", "ec_odd_fit", function(object, newdata = NULL, ...)
    rep(0, nrow(newdata)))
  odd <- function(train_sf) {
    f <- new_spatial_fit("ec_odd_fit", NULL, price ~ elev, "price", "elev", train_sf)
    f$response_var <- NULL
    f
  }
  cv <- .ec_quiet(cv_spatial(site, "price", "elev", fit_fn = odd, k = 3, seed = 1))
  expect_identical(unique(cv$fold_status$message),
                   paste0("fit_fn() returned a fit whose `response_var` is not one ",
                          "column name, where the call's `response_var` is 'price'"))
})

test_that("cv_spatial() reads a named response_var by the column it names", {
  site <- .ec_site()
  registerS3method("predict", "ec_lm_fit", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))
  })
  lm_fit <- function(train_sf)
    new_spatial_fit("ec_lm_fit", stats::lm(price ~ elev, sf::st_drop_geometry(train_sf)),
                    price ~ elev, "price", "elev", train_sf)
  # c(r = "price") ran every fold before the response check; the check read
  # its names as a difference and failed every fold, saying the fit of
  # 'price' was not a fit of 'price'.
  plain <- .ec_quiet(cv_spatial(site, "price", "elev", fit_fn = lm_fit, k = 3, seed = 1))
  named <- .ec_quiet(cv_spatial(site, c(r = "price"), "elev", fit_fn = lm_fit, k = 3,
                                seed = 1))
  expect_identical(as.integer(named$n_folds_succeeded), 3L)
  expect_identical(named$overall$RMSE, plain$overall$RMSE)
  expect_true(is.finite(named$overall$RMSE))
  # A named response the fit does not model still fails every fold.
  other <- .ec_quiet(cv_spatial(site, c(r = "slope"), "elev", fit_fn = lm_fit, k = 3,
                                seed = 1))
  expect_identical(as.integer(other$n_folds_succeeded), 0L)
  expect_identical(unique(other$fold_status$message),
                   paste0("fit_fn() returned a fit of 'price' (its `response_var`), ",
                          "where the call's `response_var` is 'slope'"))
})

test_that("select_features_forward() does not score a hold-out fit of another response", {
  site <- .ec_site(120)
  registerS3method("predict", "ec_lm_fit", function(object, newdata = NULL, ...) {
    if (is.null(newdata)) newdata <- object$data_sf
    as.numeric(stats::predict(object$engine, sf::st_drop_geometry(newdata)))
  })
  # A fit_fn whose response depends on the rows it is given: slope on the
  # sweep's folds, price on the whole selection half, where the hold-out
  # score is taken.  Its fit of price was scored against slope.
  n_sel <- NULL
  fit_fn <- function(tr, vars) {
    resp <- if (!is.null(n_sel) && nrow(tr) == n_sel) "price" else "slope"
    f <- stats::reformulate(if (length(vars)) vars else "1", resp)
    new_spatial_fit("ec_lm_fit", stats::lm(f, sf::st_drop_geometry(tr)), f, resp, vars, tr)
  }
  split <- .spatial_half_split(site, seed = 1, caller = "test")
  n_sel <- length(split$selection)
  w <- character(0)
  fs <- logger::with_log_threshold(
    withCallingHandlers(
      suppressMessages(select_features_forward(site, "slope", "elev", fit_fn = fit_fn, k = 3,
                                               seed = 1, select_on = "split",
                                               method = "random_kfold", tol = -Inf,
                                               rule = "tol")),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }),
    threshold = logger::FATAL, namespace = "spatialkit", index = 2)
  expect_identical(fs$selected, "elev")
  expect_identical(fs$score_holdout, NA_real_)
  expect_true(any(grepl(paste0("the hold-out score could not be computed (fit_fn() returned ",
                               "a fit of 'price' (its `response_var`), where the call's ",
                               "`response_var` is 'slope')"), w, fixed = TRUE)))
})

# --- kriging_adequacy(): the folds supplied ----------------------------------

skip_if_not_installed("gstat")

# 200 points over 4 cells, and a variogram fitted once, so each call kriges
# and cross-validates only.
.ec_ka <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(1)
    n <- 200
    p <- sf::st_as_sf(data.frame(x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000)),
                      coords = c("x", "y"), crs = 32617)
    p$z <- sin(sf::st_coordinates(p)[, 1] / 200) + stats::rnorm(n, 0, 0.3)
    bb <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000), crs = 32617)
    cells <- sf::st_sf(cell_id = 1:4, geometry = sf::st_make_grid(sf::st_as_sfc(bb), n = c(2, 2)))
    a <- .ec_quiet(assign_features_to_polygons(p, cells, polygon_id_col = "cell_id",
                                               quiet = TRUE))
    sac <- .ec_quiet(estimate_sac_range(a, "z"))
    cache <<- list(a = a, cells = cells, sac = sac, lab = rep(1:4, length.out = n))
    cache
  }
})
.ec_run <- function(a, folds) {
  L <- .ec_ka()
  .ec_quiet(kriging_adequacy(a, "z", L$cells, id_col = "cell_id", sac = L$sac,
                             folds = folds))
}
.ec_splits <- function(lab) list(folds = lapply(sort(unique(lab)), function(f)
  list(train = which(lab != f), test = which(lab == f))))

test_that("kriging_adequacy() refuses folds whose test sets share a location, before fitting", {
  L <- .ec_ka()
  n <- nrow(L$a)
  good <- .ec_splits(L$lab)$folds
  shared <- lapply(1:4, function(i) {
    te <- union(which(L$lab == i), 1:50)
    list(train = setdiff(seq_len(n), te), test = te)
  })
  # Kriged four times, rows 1 to 50 kept the last fold's prediction, and the
  # statistic changed with the order of the folds.  The refusal comes before
  # the variogram is estimated.
  local_mocked_bindings(estimate_sac_range = function(...) stop("estimated"))
  for (f in list(shared, rev(shared)))
    expect_error(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id",
                                  folds = list(folds = f)),
                 paste0("^kriging_adequacy\\(\\): the location of row [0-9]+ is held ",
                        "out by both fold 1 and fold 2 of `folds`; each location is ",
                        "kriged once, so a location held out twice would keep the ",
                        "last fold's prediction\\. Pass folds whose test sets are ",
                        "disjoint, such as a make_folds\\(\\) result\\.$"))
  # The row named is one both folds hold out.
  e <- tryCatch(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id",
                                 folds = list(folds = shared)), error = identity)
  r <- as.integer(sub(".*the location of row ([0-9]+) .*", "\\1", conditionMessage(e)))
  expect_true(r %in% shared[[1]]$test && r %in% shared[[2]]$test)
  # One row in the test sets of two folds is enough.
  twice <- good
  twice[[2]]$test <- c(twice[[2]]$test, twice[[1]]$test[1])
  twice[[2]]$train <- setdiff(twice[[2]]$train, twice[[1]]$test[1])
  expect_error(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id",
                                folds = list(folds = twice)),
               sprintf("the location of row %d is held out by both fold 1 and fold 2",
                       twice[[1]]$test[1]), fixed = TRUE)
  # The folds are named by their positions in `folds`, not next to each
  # other: a row of fold 1 also in fold 3's test set, and one of fold 2 in
  # fold 4's.
  for (ab in list(c(1L, 3L), c(2L, 4L))) {
    f <- good
    r1 <- f[[ab[1]]]$test[1]
    f[[ab[2]]]$test <- c(f[[ab[2]]]$test, r1)
    f[[ab[2]]]$train <- setdiff(f[[ab[2]]]$train, r1)
    expect_error(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id",
                                  folds = list(folds = f)),
                 sprintf("the location of row %d is held out by both fold %d and fold %d of `folds`",
                         r1, ab[1], ab[2]), fixed = TRUE)
  }
})

test_that("kriging_adequacy() gives disjoint folds the same statistic in either order", {
  L <- .ec_ka()
  good <- .ec_splits(L$lab)
  a <- .ec_run(L$a, good)
  b <- .ec_run(L$a, list(folds = rev(good$folds)))
  expect_identical(attr(a, "cv"), attr(b, "cv"))
  expect_identical(attr(a, "cv")$n_pred, nrow(L$a))
})

test_that("kriging_adequacy() reads one fold label per row of assigned_points_sf when a row is dropped", {
  L <- .ec_ka()
  b <- L$a
  b$z[3] <- NA
  # One label per row of the layer was refused once row 3 was dropped.
  full <- .ec_run(b, L$lab)
  kept <- .ec_run(b, L$lab[-3])
  expect_identical(attr(full, "cv"), attr(kept, "cv"))
  expect_identical(attr(full, "cv")$n_pred, 199L)
  expect_identical(sf::st_drop_geometry(full), sf::st_drop_geometry(kept))
  # The same as the layer the caller cut row 3 from.
  expect_identical(attr(full, "cv"), attr(.ec_run(L$a[-3, ], L$lab[-3]), "cv"))
  # By position in the layer as given: the rows shuffled with their labels.
  set.seed(7)
  o <- sample(nrow(b))
  expect_identical(attr(.ec_run(b[o, ], L$lab[o]), "cv"), attr(full, "cv"))
  # The labels of dropped rows are not used.
  lab2 <- L$lab
  lab2[3] <- 99L
  expect_identical(attr(.ec_run(b, lab2), "cv"), attr(full, "cv"))
  # Rows dropped for no location are read the same way.
  e <- L$a
  sf::st_geometry(e)[5] <- sf::st_point()
  expect_identical(attr(.ec_run(e, L$lab), "cv"),
                   attr(.ec_run(L$a[-5, ], L$lab[-5]), "cv"))
})

test_that("kriging_adequacy() with no row dropped reads a label vector as the same splits", {
  L <- .ec_ka()
  lab_cv <- attr(.ec_run(L$a, L$lab), "cv")
  expect_identical(lab_cv, attr(.ec_run(L$a, .ec_splits(L$lab)), "cv"))
  expect_identical(lab_cv$k, 4L)
  expect_identical(lab_cv$n_pred, 200L)
})

test_that("kriging_adequacy() refuses a label vector of another length, naming its length and the rows, before fitting", {
  L <- .ec_ka()
  b <- L$a
  b$z[3] <- NA
  local_mocked_bindings(estimate_sac_range = function(...) stop("estimated"))
  expect_error(kriging_adequacy(b, "z", L$cells, id_col = "cell_id", folds = L$lab[1:150]),
               paste0("kriging_adequacy(): `folds` has 150 labels; give one per row ",
                      "of `assigned_points_sf` (200), or one per row kept (199, after ",
                      "dropping 1 with a missing response or no location), or a ",
                      "make_folds() result."), fixed = TRUE)
  expect_error(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id", folds = L$lab[1:199]),
               paste0("kriging_adequacy(): `folds` has 199 labels; give one per row ",
                      "of `assigned_points_sf` (200), or a make_folds() result."),
               fixed = TRUE)
  expect_error(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id", folds = 1L),
               "`folds` has 1 label; give one per row", fixed = TRUE)
  # A vector longer than the layer is refused too, with a row dropped or not;
  # it is not read by position.
  expect_error(kriging_adequacy(L$a, "z", L$cells, id_col = "cell_id",
                                folds = c(L$lab, 1L)),
               paste0("kriging_adequacy(): `folds` has 201 labels; give one per row ",
                      "of `assigned_points_sf` (200), or a make_folds() result."),
               fixed = TRUE)
  expect_error(kriging_adequacy(b, "z", L$cells, id_col = "cell_id", folds = c(L$lab, 1L)),
               paste0("kriging_adequacy(): `folds` has 201 labels; give one per row ",
                      "of `assigned_points_sf` (200), or one per row kept (199, after ",
                      "dropping 1 with a missing response or no location), or a ",
                      "make_folds() result."), fixed = TRUE)
  # The helper that reads the labels takes the two lengths only.
  pts <- L$a[1:10, ]
  expect_identical(.fold_labels_for(rep(1:2, 5), pts), rep(1:2, 5))
  expect_identical(.fold_labels_for(rep(1:2, 6), pts, given = c(1:4, 6:11), n_given = 12L),
                   rep(1:2, 6)[c(1:4, 6:11)])
  for (m in c(9L, 11L, 13L))
    expect_error(.fold_labels_for(rep(1:2, length.out = m), pts, given = c(1:4, 6:11),
                                  n_given = 12L),
                 "a vector of fold labels with one per row of `assigned_points_sf`",
                 fixed = TRUE)
})

test_that("?kriging_adequacy states the length of a label vector and refuses shared test sets", {
  file <- test_path("..", "..", "man", "kriging_adequacy.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["kriging_adequacy.Rd"]]
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  for (s in c(
    "or a vector of fold labels, one per row of assigned_points_sf (the labels of rows dropped for a missing response or location are not used; one label per row kept is read as well)",
    "So are folds whose test sets share a location: each location is kriged once, under one fold."))
    expect_match(txt, s, fixed = TRUE)
})

# --- create_voronoi_polygons(): seeds fewer than their record's count --------

.ec_sel <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(4)
    n <- 600
    xy <- data.frame(x = 5e5 + stats::runif(n, 0, 5000), y = 5e6 + stats::runif(n, 0, 5000))
    xy$z <- sin(xy$x / 800) + cos(xy$y / 900) + stats::rnorm(n, sd = 0.3)
    pts <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
    prof <- .ec_quiet(resolution_profile(pts, "z", levels = c(4, 8, 16, 32), seed = 1,
                                         estimate = "mean"))
    sel <- .ec_quiet(select_resolution(prof))
    bnd <- sf::st_as_sf(sf::st_buffer(sf::st_convex_hull(sf::st_union(pts)), 1))
    cache <<- list(pts = pts, sel = sel, bnd = bnd)
    cache
  }
})

.ec_summary <- function(pts, cells) {
  asg <- .ec_quiet(assign_features_to_polygons(pts, cells, polygon_id_col = "cell_id"))
  classes <- character(0)
  v <- logger::with_log_threshold(
    withCallingHandlers(summarize_by_cell(asg, "z", cells_sf = cells),
                        warning = function(w) {
                          classes <<- c(classes, class(w)[1L])
                          invokeRestart("muffleWarning")
                        }),
    threshold = logger::FATAL, namespace = "spatialkit", index = 2)
  list(asg = asg, classes = classes, check = attr(v, "count_rule")$check)
}

test_that("seeds fewer than their record's count give cells whose record says they are not the scored cells", {
  S <- .ec_sel()
  n_sel <- nrow(S$sel$seeds)
  expect_identical(as.integer(.count_rule_get(S$sel$seeds)$cells), n_sel)
  half <- S$sel$seeds[seq_len(n_sel / 2), ]
  # The record survives row subsetting, and says scored_cells TRUE.
  expect_true(isTRUE(.count_rule_get(half)$scored_cells))
  for (route in c("build_tessellation", "create_voronoi_polygons")) {
    tess <- .ec_quiet(if (route == "build_tessellation")
      build_tessellation(half, boundary = S$bnd, method = "voronoi", quiet = TRUE)
      else create_voronoi_polygons(half, boundary = S$bnd, quiet = TRUE))
    rc <- .count_rule_get(tess$cells)
    expect_identical(nrow(tess$cells), as.integer(n_sel / 2))
    expect_identical(as.integer(rc$cells), n_sel)
    expect_false(rc$scored_cells)
    expect_false(tess$params$count_rule$scored_cells)
    # The record as get_voronoi_seeds(method = "provided") gives these seeds.
    prov <- get_voronoi_seeds(method = "provided", seeds = half)
    expect_identical(rc, .count_rule_get(prov))
    # Every row is assigned, and the summary says why it does not check the
    # cells, without naming causes that did not happen.
    s <- .ec_summary(S$pts, tess$cells)
    expect_identical(sum(is.na(s$asg$cell_id)), 0L)
    expect_false("spatialkit_cells_differ" %in% s$classes)
    expect_identical(s$check$cells, "not checked: other cells were built from this count")
  }
  # The cells are those of the same seeds through get_voronoi_seeds().
  a <- .ec_quiet(build_tessellation(half, boundary = S$bnd, method = "voronoi", quiet = TRUE))
  b <- .ec_quiet(build_tessellation(get_voronoi_seeds(method = "provided", seeds = half),
                                    boundary = S$bnd, method = "voronoi", quiet = TRUE))
  expect_identical(sf::st_geometry(a$cells), sf::st_geometry(b$cells))
})

test_that("seeds more than their record's count give cells whose record says they are not the scored cells", {
  S <- .ec_sel()
  n_sel <- nrow(S$sel$seeds)
  # The selection's seeds and one more, at the centre of the boundary, with
  # the record put back as the caller who binds them would keep it.
  extra <- S$sel$seeds[1, ]
  sf::st_geometry(extra) <- sf::st_centroid(sf::st_geometry(S$bnd))
  more <- rbind(S$sel$seeds, extra)
  attr(more, "count_rule") <- attr(S$sel$seeds, "count_rule")
  expect_true(isTRUE(.count_rule_get(more)$scored_cells))
  for (route in c("build_tessellation", "create_voronoi_polygons")) {
    tess <- .ec_quiet(if (route == "build_tessellation")
      build_tessellation(more, boundary = S$bnd, method = "voronoi", quiet = TRUE)
      else create_voronoi_polygons(more, boundary = S$bnd, quiet = TRUE))
    rc <- .count_rule_get(tess$cells)
    expect_identical(nrow(tess$cells), n_sel + 1L)
    expect_identical(as.integer(rc$cells), n_sel)
    expect_false(rc$scored_cells)
    expect_false(tess$params$count_rule$scored_cells)
    expect_identical(rc, .count_rule_get(get_voronoi_seeds(method = "provided", seeds = more)))
  }
})

test_that("seeds cut down and moved to another CRS are tessellated in theirs, as provided seeds are", {
  S <- .ec_sel()
  # The profile scored its seeds in EPSG:32632; seeds that are not those
  # centres are not moved back into it.
  half <- sf::st_transform(S$sel$seeds[seq_len(nrow(S$sel$seeds) / 2), ], 3857)
  expect_false(is.null(.count_rule_get(half)))
  a <- .ec_quiet(build_tessellation(half, boundary = S$bnd, method = "voronoi", quiet = TRUE))
  b <- .ec_quiet(build_tessellation(get_voronoi_seeds(method = "provided", seeds = half),
                                    boundary = S$bnd, method = "voronoi", quiet = TRUE))
  v <- .ec_quiet(create_voronoi_polygons(half, boundary = S$bnd, quiet = TRUE))
  expect_identical(sf::st_crs(a$cells), sf::st_crs(3857))
  expect_identical(sf::st_crs(a$cells), sf::st_crs(b$cells))
  expect_identical(sf::st_crs(v$cells), sf::st_crs(b$cells))
  expect_identical(sf::st_geometry(a$cells), sf::st_geometry(b$cells))
})

test_that("a selection's seeds as they come still give the scored cells", {
  S <- .ec_sel()
  tess <- .ec_quiet(build_tessellation(S$sel$seeds, boundary = S$bnd, method = "voronoi",
                                       quiet = TRUE))
  rc <- .count_rule_get(tess$cells)
  expect_true(rc$scored_cells)
  expect_identical(rc, .count_rule_get(S$sel$seeds))
  s <- .ec_summary(S$pts, tess$cells)
  expect_identical(s$check$cells, "match")
})

test_that("?create_voronoi_polygons says seeds fewer or more than the count are not the scored centres", {
  file <- test_path("..", "..", "man", "create_voronoi_polygons.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file)
        else tools::Rd_db("spatialkit")[["create_voronoi_polygons.Rd"]]
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("So are seeds fewer or more than the record's count, such as a ",
                           "selection's seeds cut down to some of its rows"), fixed = TRUE)
})
