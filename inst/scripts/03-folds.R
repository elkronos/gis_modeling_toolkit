# =============================================================================
# 03  Fold schemes: what each one's score estimates
# =============================================================================
#   source(system.file("scripts", "03-folds.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()
bnd <- tour_boundary(pts)

rmse_of <- function(folds, label) {
  cv <- cv_spatial(pts, "z", "elev", fit_fn = tour_fit, folds = folds)
  r  <- cv$overall$RMSE[1L]
  cat(sprintf("  %-28s k=%d  RMSE %.3f\n", label, folds$k, r))
  invisible(list(cv = cv, rmse = r))
}

step("03.1", "The same model, two fold schemes")
# The learner here is a cubic trend surface, so it can memorise where things are.
# Random folds leave a test point's neighbours in the training set, so the
# score is the error of predicting among the samples, where the memorised
# surface holds. Block folds hold out whole blocks, so the score is the error
# of predicting into gaps a block wide, where it does not.
f_rand  <- make_folds(pts, k = 5, method = "random_kfold", seed = 1)
f_block <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                      boundary = bnd, block_size = 250)
r1 <- rmse_of(f_rand,  "random_kfold")
r2 <- rmse_of(f_block, "block_kfold (250 m)")
cat(sprintf("  ratio %.2f: how much worse the model predicts into a gap a block\n",
            r2$rmse / r1$rmse),
    "  wide than among the samples\n", sep = "")
cat("  Neither number is the error of every map. Report the one whose folds\n",
    "  match where the map will predict: random for filling in among the\n",
    "  samples, blocked for predicting into gaps a block wide.\n", sep = "")

# fold_separation() says which one matches a map of the study area: it sets
# the distances from the map's locations to the data beside the held-out
# distances of each scheme.
set.seed(3)
map_pts <- sf::st_sample(bnd, 900, type = "regular")
reading <- function(f)
  attr(fold_separation(f, pts, prediction_points = map_pts), "target")$reading
cat(sprintf("  against a map of the study area: random folds read \"%s\", blocks \"%s\"\n",
            reading(f_rand), reading(f_block)))
cat("  \"matched\" is the error of that map; \"farther\" tends to read high for it\n",
    "  and \"closer\" low.\n", sep = "")

step("03.2", "A block size from the autocorrelation range")
# `auto_range = TRUE` fits a variogram and uses the estimated range as the
# minimum block size, so blocks are at least as wide as the autocorrelation
# and the score is the error of a map that has to reach that far from the data.
f_auto <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                     boundary = bnd, auto_range = TRUE, response_var = "z")
cat(sprintf("  estimated range %.0f -> block_size %.0f, %d blocks over %d folds\n",
            f_auto$params$sac_range, f_auto$params$block_size,
            f_auto$params$blocks_used, f_auto$k))
r3 <- rmse_of(f_auto, "block_kfold (auto_range)")
if (r3$rmse > 3 * r2$rmse) {
  cat(sprintf("  Far worse again: with %d blocks each fold hides about %.0f percent\n",
              f_auto$params$blocks_used, 100 / f_auto$k))
  cat("  of the map, and a cubic trend surface extrapolated that far returns\n",
      "  values well outside the data. Keep the blocks, and stop predicting\n",
      "  where the model has no support; script 07 measures where that is.\n",
      sep = "")
}

step("03.3", "Draw the folds")
if (!skip_without("ggplot2", "the fold maps")) {
  look_for("random folds: every colour is mixed in everywhere. Nothing is held ",
           "out except rows, so a test point's neighbours are in the training set.")
  show_plot(plot_folds(f_rand, pts, boundary = bnd),  "03-folds-random.png")
  look_for("block folds: every point inside a grid cell shares one colour, so ",
           "a fold holds out whole blocks. A fold can still be scattered across ",
           "the map, because make_folds() gives each block to the fold with the ",
           "fewest points so far, wherever the block lies.")
  show_plot(plot_folds(f_block, pts, boundary = bnd), "03-folds-block.png")
}

step("03.4", "Two schemes for specific questions")
# buffered_loo: hold out one point and everything within `buffer` of it. The
# score is the error of predicting a location at least `buffer` from every
# training sample, and it costs one model fit per row.
f_loo <- make_folds(pts, k = 5, method = "buffered_loo", seed = 1, buffer = 240)
cat(sprintf("  buffered_loo: %d folds; fold 1 trains on %d of %d rows\n",
            f_loo$k, length(f_loo$folds[[1]]$train), nrow(pts)))

# nndm: matches the train-to-test distance distribution of the CV to the
# distance distribution of the prediction task. Use it when you know where you
# will predict: here, a regular grid over the area rather than the sample
# locations, the usual case for a map.
set.seed(7)
grid_pts <- sf::st_sf(geometry = sf::st_sample(bnd, 200, type = "regular"))
f_nndm <- make_folds(pts, k = 5, method = "nndm", seed = 1,
                     prediction_points = grid_pts)
cat(sprintf("  nndm: prediction points sit a median %.0f m from the nearest sample;\n",
            f_nndm$params$target_median))
cat(sprintf("        the folds reproduce that at %.0f m\n",
            f_nndm$params$realised_median))
cat("  Hand it the points you will predict on; other points set a different task.\n")

step("03.5", "Read k off the folds, do not assume it")
# When fewer than k blocks hold data, make_folds() returns fewer folds, and a
# CV you write up as 5-fold is not one. It happened in 03.2: make_folds()
# warned (class "spatialkit_k_reduced", also logged as the WARN line there) and
# recorded the number of folds it built in `$k`. Read `$k` off the folds; it
# stays right where a warning is suppressed or missed.
cat(sprintf("  03.1 asked for k = 5 and got k = %d\n", f_block$k))
cat(sprintf("  03.2 asked for k = 5 and got k = %d\n", f_auto$k))

# At a larger size make_folds() stops rather than return a single fold.
bad <- try(make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                      boundary = bnd, block_size = 600), silent = TRUE)
if (inherits(bad, "try-error"))
  cat("  at block_size = 600 it stops:\n    ",
      conditionMessage(attr(bad, "condition")), "\n", sep = "")
