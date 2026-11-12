# tests/testthat/test-fix-aoa-coord-forest.R
# ---------------------------------------------------------------------------
# area_of_applicability() and a forest fitted with include_coords = TRUE.
#
# Such a forest had location measured only through its "..x" and "..y"
# predictor columns.  Under weights = "auto" those columns take the forest's
# own importance, so a forest that barely splits on its coordinates gave
# location almost no weight, and points well beyond the training extent read
# as inside.  The forest now gets the location index a GWR and a Bayesian
# fit get, and a point is inside only when it is inside on both indices.
# The predictor index is the one it was: the coordinates stay in it, weighted
# as before.
# ---------------------------------------------------------------------------

cf_pts <- function(m, x0, x1) {
  d <- data.frame(x = stats::runif(m, x0, x1), y = stats::runif(m, 0, 1000),
                  a1 = stats::rnorm(m), a2 = stats::rnorm(m))
  d$z <- 2 * d$a1 + d$a2 + 0.5 * sin(d$x / 150) + stats::rnorm(m, 0, 0.3)
  sf::st_as_sf(d, coords = c("x", "y"), crs = 32632)
}

# The class stand-in test-fix-aoa-rf-location.R uses for a GWR.
cf_fake_gwr <- function(train, vars)
  structure(list(data_sf = train, predictor_vars = vars, response_var = "z",
                 info = list(), n = nrow(train)),
            class = c("gwr_fit", "spatial_fit"))

cf_forest <- function(seed = 1) {
  set.seed(seed)
  trn <- cf_pts(300, 0, 1000)
  suppressWarnings(suppressMessages(
    fit_rf_model(trn, "z", c("a1", "a2"), num_trees = 200,
                 include_coords = TRUE, seed = seed)))
}


test_that("a forest fitted with coordinates puts points beyond the extent outside", {
  skip_if_not_installed("ranger")
  rf <- cf_forest(1)
  imp <- pmax(rf$info$importance, 0)
  # The case of the bug: the forest gives its coordinates little importance.
  expect_lt(sum(imp[c("..x", "..y")]) / sum(imp), 0.1)

  set.seed(2)
  far <- cf_pts(200, 1500, 2500)
  a <- area_of_applicability(far, model = rf)
  # On the predictor index alone most of these points read as inside; that
  # was the whole answer before.
  expect_gt(mean(a$aoa$DI <= a$threshold), 0.3)
  # With location measured on its own, none is inside.
  expect_identical(a$n_inside, 0L)
  expect_identical(a$location$n_outside, 200L)
  expect_identical(a$params$coords, "location index and predictors")
  expect_true(all(c("..x", "..y", "DI", "DI_location", "AOA") %in%
                    names(a$aoa)))
  expect_identical(a$aoa$AOA,
                   a$aoa$DI <= a$threshold &
                     a$aoa$DI_location <= a$location$threshold)
  expect_identical(a$n_outside_predictors,
                   sum(a$aoa$DI > a$threshold, na.rm = TRUE))
  expect_identical(a$n_new, a$n_inside + a$n_outside + a$n_na)

  # Inside the extent nearly every point is still inside.
  set.seed(3)
  within <- area_of_applicability(cf_pts(200, 0, 1000), model = rf)
  expect_gt(mean(within$aoa$AOA), 0.75)

  # It is the location index a GWR gets on the same training rows.
  g <- area_of_applicability(far, model = cf_fake_gwr(rf$data_sf, c("a1", "a2")))
  expect_identical(a$aoa$DI_location, g$aoa$DI_location)
  expect_identical(a$location$threshold, g$location$threshold)
  expect_identical(a$location$normalizer, g$location$normalizer)
  expect_identical(a$location$train_DI, g$location$train_DI)
})

test_that("the forest's predictor index is the one it was", {
  skip_if_not_installed("ranger")
  rf <- cf_forest(4)
  set.seed(5)
  nd <- cf_pts(100, 0, 2500)
  a <- area_of_applicability(nd, model = rf)
  expect_identical(a$predictor_vars, c("a1", "a2", "..x", "..y"))
  expect_identical(names(a$weights), c("a1", "a2", "..x", "..y"))
  # The same index built by hand: the coordinates as two ordinary columns,
  # weighted by the forest's importance, and no location index.
  xy_cols <- function(s) {
    crd <- sf::st_coordinates(s)
    s$..x <- crd[, 1]; s$..y <- crd[, 2]
    s
  }
  hand <- area_of_applicability(xy_cols(nd), train_sf = xy_cols(rf$data_sf),
                                predictor_vars = c("a1", "a2", "..x", "..y"),
                                weights = pmax(rf$info$importance, 0))
  expect_null(hand$location)
  expect_identical(a$aoa$DI, hand$aoa$DI)
  expect_identical(a$threshold, hand$threshold)
  expect_identical(a$train_DI, hand$train_DI)
  expect_identical(a$weights, hand$weights)
  expect_identical(a$n_outside_predictors, sum(hand$aoa$DI > hand$threshold))
})

test_that("location counts whatever weight the coordinates get", {
  skip_if_not_installed("ranger")
  rf <- cf_forest(6)
  set.seed(7)
  nd <- rbind(cf_pts(25, 0, 1000), cf_pts(25, 1500, 2500))
  beyond <- seq_len(50) > 25
  for (w in list("auto", NULL, c(a1 = 1, a2 = 1, ..x = 0, ..y = 0))) {
    a <- area_of_applicability(nd, model = rf, weights = w)
    expect_identical(sum(a$aoa$AOA[beyond]), 0L)
  }
  # A threshold the caller passes sets the predictor index alone.
  big <- area_of_applicability(nd, model = rf, threshold = 1e6)
  expect_true(big$params$threshold_supplied)
  expect_identical(big$n_outside_predictors, 0L)
  expect_identical(sum(big$aoa$AOA[beyond]), 0L)
  # Under the forest's own blocked folds too.
  fo <- suppressWarnings(suppressMessages(
    make_folds(rf$data_sf, k = 5, method = "block_kfold", seed = 1)))
  af <- suppressWarnings(area_of_applicability(nd, model = rf, folds = fo))
  expect_identical(sum(af$aoa$AOA[beyond]), 0L)
  expect_gt(sum(af$aoa$AOA[!beyond]), 20L)
})

test_that("coords = \"yes\" is the default for a forest, and \"no\" leaves location out", {
  skip_if_not_installed("ranger")
  rf <- cf_forest(8)
  set.seed(9)
  nd <- cf_pts(60, 0, 2500)
  a   <- area_of_applicability(nd, model = rf)
  yes <- area_of_applicability(nd, model = rf, coords = "yes")
  expect_identical(yes$aoa$DI, a$aoa$DI)
  expect_identical(yes$aoa$DI_location, a$aoa$DI_location)
  lines <- capture_spatialkit_log(
    off <- area_of_applicability(nd, model = rf, coords = "no"))
  expect_true(log_has(lines, "coords = \"no\", so location is not measured"))
  expect_true(log_has(lines, "a forest fitted with include_coords = TRUE"))
  expect_null(off$location)
  expect_identical(off$predictor_vars, c("a1", "a2"))
  expect_identical(off$params$coords, "none")
})

test_that("the location index of a forest reads the layer the predictor index reads", {
  skip_if_not_installed("ranger")
  rf <- cf_forest(10)
  set.seed(11)
  nd <- cf_pts(30, 0, 1000)
  # A CRS on one side only is reconciled once, with one warning.
  nocrs <- sf::st_set_crs(nd, NA)
  n_warn <- 0L
  a <- withCallingHandlers(area_of_applicability(nocrs, model = rf),
                           warning = function(w) {
                             n_warn <<- n_warn + 1L
                             invokeRestart("muffleWarning")
                           })
  expect_identical(n_warn, 1L)
  expect_equal(a$aoa$DI_location,
               area_of_applicability(nd, model = rf)$aoa$DI_location)
  # Polygons are measured at the points the predictor index reads.
  cells <- sf::st_sf(a1 = stats::rnorm(9), a2 = stats::rnorm(9),
                     geometry = sf::st_make_grid(rf$data_sf, n = c(3, 3)))
  ap <- suppressWarnings(area_of_applicability(cells, model = rf))
  crd <- sf::st_coordinates(ap$aoa)
  expect_identical(unname(crd[, 1]), ap$aoa$..x)
  expect_identical(
    ap$aoa$DI_location,
    area_of_applicability(ap$aoa, model = cf_fake_gwr(rf$data_sf, c("a1", "a2")))$aoa$DI_location)
  # A data.frame has no location to measure.
  expect_error(area_of_applicability(sf::st_drop_geometry(nd), model = rf),
               "coordinates as predictors.*coords = \"no\"")
  # Training data in longitude/latitude has no map distance.
  expect_error(suppressWarnings(
    area_of_applicability(nd, model = rf,
                          train_sf = sf::st_transform(rf$data_sf, 4326))),
    "geographic CRS")
})

test_that("print() names both ways location is measured for a forest", {
  skip_if_not_installed("ranger")
  rf <- cf_forest(12)
  set.seed(13)
  a <- area_of_applicability(cf_pts(40, 1500, 2500), model = rf)
  out <- utils::capture.output(print(a))
  expect_true(any(grepl("location    : measured on its own", out, fixed = TRUE)))
  expect_true(any(grepl("and as the predictors ..x and ..y", out, fixed = TRUE)))
  expect_false(any(grepl("location    : measured as the predictors", out,
                         fixed = TRUE)))
  expect_true(any(grepl("on location", out, fixed = TRUE)))
  # plot() counts the two indices apart, as for a GWR.
  skip_if_not_installed("ggplot2")
  p <- plot(a)
  expect_match(p$labels$subtitle,
               sprintf("%d outside on location; %d outside the AOA (either)",
                       a$location$n_outside, a$n_outside), fixed = TRUE)
  expect_match(p$labels$caption, "location is a second index", fixed = TRUE)
  expect_no_error(ggplot2::ggplot_build(plot(a, type = "histogram")))
})

test_that("only a forest's print() names the coordinates as predictors", {
  set.seed(23)
  trn <- cf_pts(150, 0, 1000)
  out <- utils::capture.output(print(
    area_of_applicability(cf_pts(20, 0, 2000), model = cf_fake_gwr(trn, c("a1", "a2")))))
  expect_true(any(grepl("location    : measured on its own", out, fixed = TRUE)))
  expect_false(any(grepl("..x", out, fixed = TRUE)))
  # An object made before the fix (params$coords "predictors", no location)
  # still prints what it measured.
  old <- area_of_applicability(cf_pts(20, 0, 1000), train_sf = trn,
                               predictor_vars = c("a1", "a2"))
  old$params$coords <- "predictors"
  out_old <- utils::capture.output(print(old))
  expect_true(any(grepl("location    : measured as the predictors ..x and ..y",
                        out_old, fixed = TRUE)))
  expect_false(any(grepl("measured on its own", out_old, fixed = TRUE)))
})

test_that("the forest's location index takes the folds and repeat rule a GWR's takes", {
  skip_if_not_installed("ranger")
  set.seed(21)
  s <- cf_pts(75, 0, 1000)
  trn <- s[rep(seq_len(75), each = 2), ]          # two visits per site
  trn$a1 <- stats::rnorm(150)
  trn$z  <- 2 * trn$a1 + trn$a2 + stats::rnorm(150, 0, 0.3)
  rf <- suppressWarnings(suppressMessages(
    fit_rf_model(trn, "z", c("a1", "a2"), num_trees = 100,
                 include_coords = TRUE, seed = 1)))
  g  <- cf_fake_gwr(rf$data_sf, c("a1", "a2"))
  set.seed(22)
  nd <- rbind(cf_pts(30, 0, 1000), cf_pts(30, 1000, 1400))
  nd$a1[3] <- NA                                  # missing predictor, in the extent
  nd$a1[60] <- NA                                 # missing predictor, far beyond it
  sf::st_geometry(nd)[60] <- sf::st_point(c(3000, 500))
  sf::st_geometry(nd)[5] <- sf::st_point()        # no location
  # Repeated locations, no folds.
  a  <- suppressWarnings(area_of_applicability(nd, model = rf))
  ag <- suppressWarnings(area_of_applicability(nd, model = g))
  expect_gt(a$location$n_duplicates, 0L)
  expect_identical(a$location, ag$location)
  expect_identical(a$aoa$DI_location, ag$aoa$DI_location)
  # A row with a missing predictor is outside when it is outside on location,
  # and unclassified when it is inside there; a row with no location is
  # unclassified.
  expect_true(is.na(a$aoa$DI[3]) && is.na(a$aoa$DI[60]))
  expect_lte(a$aoa$DI_location[3], a$location$threshold)
  expect_gt(a$aoa$DI_location[60], a$location$threshold)
  expect_identical(a$aoa$AOA[c(3L, 5L, 60L)], c(NA, NA, FALSE))
  expect_true(is.na(a$aoa$DI_location[5]))
  expect_identical(a$location$n_outside,
                   sum(a$aoa$DI_location > a$location$threshold, na.rm = TRUE))
  expect_identical(a$n_na, sum(is.na(a$aoa$AOA)))
  expect_identical(a$n_new, a$n_inside + a$n_outside + a$n_na)
  # Block folds.
  fo <- suppressWarnings(suppressMessages(
    make_folds(rf$data_sf, k = 5, method = "block_kfold", seed = 1)))
  af  <- suppressWarnings(area_of_applicability(nd, model = rf, folds = fo))
  agf <- suppressWarnings(area_of_applicability(nd, model = g, folds = fo))
  expect_false(identical(af$location$threshold, a$location$threshold))
  expect_identical(af$location, agf$location)
  expect_identical(af$aoa$DI_location, agf$aoa$DI_location)
})

test_that("the help says what changes for a forest", {
  root <- testthat::test_path("..", "..")
  if (!(file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man"))))
    skip("the help's source is not here to read")
  read <- function(...) gsub("[[:space:]]+", " ",
                             paste(readLines(file.path(root, ...), warn = FALSE), collapse = " "))
  rd <- read("man", "area_of_applicability.Rd")
  # Under coords = "no" the forest's `aoa` is newdata as passed.
  expect_match(rd, paste0("For a forest fitted with \\code{include_coords = TRUE}, unless ",
                          "\\code{coords = \"no\"}, it is the object the computation ran on"),
               fixed = TRUE)
})
