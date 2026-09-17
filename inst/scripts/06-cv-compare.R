# =============================================================================
# 06  Comparing models without fooling yourself
# =============================================================================
#   source(system.file("scripts", "06-cv-compare.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()
bnd <- tour_boundary(pts)

# Three candidates, no optional dependencies. `tour_fit` (a cubic trend surface
# plus elevation) comes from _common.R; the other two are defined here.
mean_fit <- function(train_sf, ...) {
  new_spatial_fit(subclass = "mean_fit",
                  engine = stats::lm(z ~ 1, data = sf::st_drop_geometry(train_sf)),
                  formula = z ~ 1, response_var = "z",
                  predictor_vars = character(0), data_sf = train_sf)
}
elev_fit <- function(train_sf, ...) {
  new_spatial_fit(subclass = "elev_fit",
                  engine = stats::lm(z ~ elev, data = sf::st_drop_geometry(train_sf)),
                  formula = z ~ elev, response_var = "z",
                  predictor_vars = "elev", data_sf = train_sf)
}
predict.mean_fit <- predict.elev_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  as.numeric(stats::predict(object$engine, newdata = sf::st_drop_geometry(newdata)))
}
registerS3method("predict", "mean_fit", predict.mean_fit)
registerS3method("predict", "elev_fit", predict.elev_fit)

step("06.1", "One set of folds for every candidate")
# This is the part people skip. Build the folds ONCE and pass them to every
# run, or the differences you read off the table are partly fold luck.
folds <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                    boundary = bnd, block_size = 250)
cat(sprintf("  %d blocked folds, block_size 250, seed 1\n", folds$k))

runs <- list(
  "intercept only" = cv_spatial(pts, "z", character(0), fit_fn = mean_fit, folds = folds),
  "elevation"      = cv_spatial(pts, "z", "elev",      fit_fn = elev_fit, folds = folds),
  "trend surface"  = cv_spatial(pts, "z", "elev",      fit_fn = tour_fit, folds = folds)
)
tab <- do.call(rbind, lapply(names(runs), function(nm)
  data.frame(model = nm, RMSE = runs[[nm]]$overall$RMSE[1L],
             MAE = runs[[nm]]$overall$MAE[1L], R2 = runs[[nm]]$overall$R2[1L])))
print(tab, row.names = FALSE, digits = 4)

step("06.2", "Compare fold by fold, not mean to mean")
# The same folds means the per-fold scores are PAIRED, so you can look at the
# differences directly instead of comparing two averages with overlapping spread.
per_fold <- sapply(runs, function(r) r$fold_metrics$RMSE)
rownames(per_fold) <- paste("fold", runs[[1]]$fold_metrics$fold)
print(round(per_fold, 3))
d <- per_fold[, "elevation"] - per_fold[, "intercept only"]
cat(sprintf("\n  elevation minus intercept, per fold: %s\n",
            paste(sprintf("%+.3f", d), collapse = "  ")))
cat(sprintf("  better on %d of %d folds; mean gain %.3f, spread %.3f\n",
            sum(d < 0), length(d), -mean(d), stats::sd(d)))
if (all(d < 0)) {
  cat("  A clean sweep: every fold agrees, so the ranking does not depend on\n",
      "  which part of the map was held out.\n", sep = "")
} else {
  cat("  It wins on average but not everywhere, so the ranking does depend on\n",
      "  which part of the map was held out. Say how many folds, not just the\n",
      "  mean.\n", sep = "")
}

# The trend surface is the counter-example: worst of the three on average, yet
# it wins outright on several individual folds.
d2 <- per_fold[, "trend surface"] - per_fold[, "intercept only"]
cat(sprintf("  trend surface beats the intercept on %d of %d folds and loses the\n",
            sum(d2 < 0), length(d2)))
cat(sprintf("  other %d badly enough to finish last overall (worst fold %.2f).\n",
            sum(d2 >= 0), max(per_fold[, "trend surface"])))

step("06.3", "Draw it")
if (!skip_without("ggplot2", "the metric plot")) {
  look_for("one point per fold, with the pooled value marked. A model whose ",
           "folds are scattered wide has not been measured precisely, whatever ",
           "its mean says.")
  show_plot(plot_cv_metrics(runs[["trend surface"]], metric = "RMSE"),
            "06-cv-metrics.png")
}

step("06.4", "The built-in backends")
# compare_models_cv() drives the package's own model families (RF, GWR,
# Bayesian) through one interface and one set of folds.
avail <- c(RF = have("ranger"), GWR = have("GWmodel"), Bayesian = have("brms"))
cat("  backends available here:",
    paste(names(avail)[avail], collapse = ", "), "\n")
if (any(avail[c("RF", "GWR")])) {
  cmp <- compare_models_cv(pts, "z", c("elev", "noise"),
                           models = names(avail)[avail & names(avail) != "Bayesian"],
                           folds = folds, quiet = TRUE, seed = 1)
  print(cmp$overall, row.names = FALSE, digits = 4)
  if (!skip_without("ggplot2", "the backend comparison plot")) {
    look_for("several models side by side on the same folds. Overlapping fold ",
             "spreads mean the ranking is not stable; a clear separation means ",
             "it is.")
    show_plot(plot_cv_metrics(cmp, metric = "RMSE"), "06-backends.png")
  }
  cat("  Bayesian is left out here because it is slow; script 10 runs it and\n",
      "  shows what its extra output buys.\n", sep = "")
} else {
  cat("  SKIPPED: install.packages(c(\"ranger\", \"GWmodel\"))\n")
}

step("06.5", "The score that picked the winner cannot also grade it")
# Pick the best of three models by cross-validated RMSE and you have spent that
# CV on choosing. Reporting the same number as the chosen model's performance
# is the in-sample mistake again, one level up.
best <- tab$model[which.min(tab$RMSE)]
cat(sprintf("  best by CV: %s at RMSE %.3f\n", best, min(tab$RMSE)))
cat(sprintf("  spread of the three: %.3f to %.3f\n", min(tab$RMSE), max(tab$RMSE)))
cat("  The honest write-up either (a) reports the selection, then measures the\n",
    "  winner on data held out from the whole comparison, or (b) reports all\n",
    "  three numbers and lets the reader see how close the race was.\n",
    "  `resolution_profile(select_on = 'split')` in script 02 is the same idea\n",
    "  applied to a resolution choice.\n", sep = "")
