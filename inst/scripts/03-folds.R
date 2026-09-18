# =============================================================================
# 03  Fold schemes: why random k-fold flatters a spatial model
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
# memorised surface is still valid there and the score looks good. Block folds
# hold out whole regions, so it is not.
f_rand  <- make_folds(pts, k = 5, method = "random_kfold", seed = 1)
f_block <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                      boundary = bnd, block_size = 250)
r1 <- rmse_of(f_rand,  "random_kfold")
r2 <- rmse_of(f_block, "block_kfold (250 units)")
cat(sprintf("  ratio %.2fx -- the optimism random folds buy you\n", r2$rmse / r1$rmse))
cat("  Neither number is the truth. Report the one that matches how the model\n",
    "  will be used: random for gap-filling inside the sampled area, blocked\n",
    "  for predicting somewhere new.\n", sep = "")

step("03.2", "Let the data choose the block size")
# `auto_range = TRUE` fits a variogram and uses the estimated range as the
# minimum block size, so blocks are at least as wide as the autocorrelation
# they are supposed to break.
f_auto <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                     boundary = bnd, auto_range = TRUE, response_var = "z")
cat(sprintf("  estimated range %.0f -> block_size %.0f, %d blocks over %d folds\n",
            f_auto$params$sac_range, f_auto$params$block_size,
            f_auto$params$blocks_used, f_auto$k))
r3 <- rmse_of(f_auto, "block_kfold (auto_range)")
if (r3$rmse > 3 * r2$rmse) {
  cat(sprintf("  Far worse again, and not a bug: %d blocks means each fold hides\n",
              f_auto$params$blocks_used))
  cat(sprintf("  about %.0f%% of the map, and a cubic trend surface asked to\n",
              100 / f_auto$k))
  cat("  extrapolate that far produces nonsense. The answer is not a gentler\n",
      "  fold scheme; it is to stop predicting where the model has no support.\n",
      "  Script 07 measures how far that is.\n", sep = "")
}

step("03.3", "Draw the folds")
if (!skip_without("ggplot2", "the fold maps")) {
  look_for("random folds: every colour is mixed in everywhere. Nothing is held ",
           "out except rows, so a test point's neighbours are in the training set.")
  show_plot(plot_folds(f_rand, pts, boundary = bnd),  "03-folds-random.png")
  look_for("block folds: every point inside a grid cell shares one colour, so ",
           "a fold holds out whole blocks. The fold itself is still scattered ",
           "across the map -- blocks are assigned to folds, not halves of the ",
           "map to folds.")
  show_plot(plot_folds(f_block, pts, boundary = bnd), "03-folds-block.png")
}

step("03.4", "Two schemes for specific questions")
# buffered_loo: hold out one point and everything within `buffer` of it. The
# strictest answer to 'how well do I predict an unsampled location', and it
# costs one model fit per row.
f_loo <- make_folds(pts, k = 5, method = "buffered_loo", seed = 1, buffer = 240)
cat(sprintf("  buffered_loo: %d folds; fold 1 trains on %d of %d rows\n",
            f_loo$k, length(f_loo$folds[[1]]$train), nrow(pts)))

# nndm: matches the train-to-test distance distribution of the CV to the
# distance distribution of the actual prediction task. Use it when you know
# where you will predict -- here, on a regular grid rather than at the sample
# locations, which is the usual mapping case.
set.seed(7)
grid_pts <- sf::st_sf(geometry = sf::st_sample(bnd, 200, type = "regular"))
f_nndm <- make_folds(pts, k = 5, method = "nndm", seed = 1,
                     prediction_points = grid_pts)
cat(sprintf("  nndm: prediction points sit %.0f units from the nearest sample;\n",
            f_nndm$params$target_median))
cat(sprintf("        the folds reproduce that at %.0f units\n",
            f_nndm$params$realised_median))
cat("  Hand it the SAME points you will predict on, or it matches the wrong task.\n")

step("03.5", "Read k off the folds, do not assume it")
# When fewer than k blocks contain data, k is reduced to match, and the CV you
# write up as 5-fold is not one. This already happened in 03.2. It arrives as a
# log line and nothing else: make_folds() raises no R condition for it, so
# tryCatch(warning =) will not catch it and suppressWarnings() will not hide
# it. Reading `$k` off the folds is the only reliable check.
cat(sprintf("  03.1 asked for k = 5 and got k = %d\n", f_block$k))
cat(sprintf("  03.2 asked for k = 5 and got k = %d\n", f_auto$k))

# Push it further and the package refuses instead of returning one useless fold.
bad <- try(make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                      boundary = bnd, block_size = 600), silent = TRUE)
if (inherits(bad, "try-error"))
  cat("  at block_size = 600 it stops:\n    ",
      conditionMessage(attr(bad, "condition")), "\n", sep = "")
