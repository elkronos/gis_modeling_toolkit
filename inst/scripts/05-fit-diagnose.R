# =============================================================================
# 05  Fit a model, then argue with it
# =============================================================================
#   source(system.file("scripts", "05-fit-diagnose.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()

step("05.1", "Wrap a model so the rest of the package can use it")
# `new_spatial_fit()` is the adapter. Give it your fitted engine, the formula,
# the response and predictor names, and the training layer. A `predict` method
# is enough for the cv_* functions; add `fitted` and `residuals` and the
# in-sample diagnostics and plots work too. Leave `fitted` out and
# `model_metrics()` names the method it wants -- which is why script 09.5
# registers one. `_common.R` builds `tour_fit` this way, methods included, in
# about twenty-five lines.
fit <- tour_fit(pts)
print(summary(fit))

step("05.2", "In-sample numbers are not a score")
ins <- model_metrics(fit)
cv  <- cv_spatial(pts, "z", "elev", fit_fn = tour_fit,
                  folds = make_folds(pts, k = 5, method = "block_kfold",
                                     seed = 1, boundary = tour_boundary(pts),
                                     block_size = 250))
cat(sprintf("  in-sample RMSE   %.3f\n", ins$RMSE))
cat(sprintf("  blocked CV RMSE  %.3f  (%.1fx worse)\n",
            cv$overall$RMSE[1L], cv$overall$RMSE[1L] / ins$RMSE))
cat("  The first number says how well the model reproduces what it was shown.\n",
    "  Only the second is an estimate of anything you care about.\n", sep = "")

step("05.3", "What is left in the residuals")
# If the residuals are still clustered in space, something structured is
# missing: a predictor, a trend term, or a correlated error structure.
mi <- residual_morans_i(fit)
print(mi)
cat(sprintf("\n  I = %.3f, p = %.2g: %s\n", mi$observed, mi$p_value,
            if (mi$p_value < 0.05)
              "the residuals are still spatially organised" else
              "no detectable structure left in the residuals"))
cat("  Structure left over means a predictor, a trend term or a correlated\n",
    "  error structure is missing. Run the same statistic again after adding\n",
    "  one: a drop toward zero is your evidence that it did real work.\n", sep = "")

# A second look, with a different weights matrix. Moran's I is a statement
# about the neighbourhood you defined, so changing k changes the answer.
for (kk in c(4L, 8L, 16L)) {
  m <- residual_morans_i(fit, k = kk)
  cat(sprintf("  k = %2d neighbours -> I = %.3f, z = %5.2f\n", kk, m$observed, m$z))
}
cat("  Report the k you used. An unstated k is an unreproducible statistic.\n")

step("05.4", "Three plots, three questions")
if (!skip_without("ggplot2", "the diagnostic plots")) {
  look_for("residuals in space: patches of one colour mean unmodelled ",
           "structure. Salt and pepper is what you want.")
  show_plot(plot(fit, type = "residuals"), "05-residuals.png")

  look_for("observed against predicted: points should sit on the 1:1 line. ",
           "A flatter-than-1:1 cloud means the model is shrinking toward the ",
           "mean, which is what over-smoothing looks like.")
  show_plot(plot(fit, type = "observed_predicted"), "05-obs-pred.png")

  if (have("gstat")) {
    look_for("the residual variogram: a curve that still rises with distance ",
             "is the same message as 05.3, drawn. A flat line means what is ",
             "left is noise.")
    show_plot(plot(fit, type = "variogram"), "05-resid-variogram.png")
  } else {
    cat("  SKIPPED (residual variogram): install.packages(\"gstat\")\n")
  }
}

step("05.5", "Percentage error on a response that crosses zero")
# `z` here is centred near zero, so dividing by it is meaningless: a residual of
# 0.1 next to an observation of 0.001 is a 10,000% error. The metric is reported
# because it is conventional, not because it is always applicable.
cat(sprintf("  response range: %.2f to %.2f, %d observations within 0.1 of zero\n",
            min(pts$z), max(pts$z), sum(abs(pts$z) < 0.1)))
cat(sprintf("  MAPE  %8.1f%%   <- driven by the near-zero observations\n", ins$MAPE))
cat(sprintf("  SMAPE %8.1f%%   <- bounded, but still not interpretable here\n", ins$SMAPE))
cat(sprintf("  RMSE  %8.3f     <- in the units of z, which is what to report\n", ins$RMSE))
cat("  Percentage errors need a response with a meaningful zero and no values\n",
    "  near it. Rainfall and population qualify; anomalies and log-ratios do not.\n",
    sep = "")

step("05.6", "Is aggregating to cells even worth it?")
# Cell means are a summary. `kriging_adequacy()` asks whether a geostatistical
# model of the same data would beat that summary, which is the question behind
# 'should I report cell means or a surface'.
if (!skip_without("gstat", "the adequacy check")) {
  bnd   <- tour_boundary(pts)
  seeds <- get_voronoi_seeds(boundary = bnd, method = "kmeans", n = 20,
                             sample_points = pts, set_seed = 1)
  tess  <- build_tessellation(seeds, boundary = bnd, method = "voronoi", quiet = TRUE)
  asg   <- assign_features_to_polygons(pts, tess$cells)
  ka    <- kriging_adequacy(asg, "z", cells_sf = tess$cells, seed = 1)
  print(ka)

  # Two independent questions, answered by two different lines of that print.
  zv    <- attr(ka, "cv")$zscore_var
  moved <- sum(abs(ka$kr_shift) > 1, na.rm = TRUE)
  cat(sprintf("\n  are the kriging uncertainties honest?  z-score variance %.2f\n", zv))
  cat(sprintf("     %s\n",
    if (zv > 0.7 && zv < 1.4) "yes, near 1: the standard errors mean what they say"
    else if (zv >= 1.4) "no, above 1: the model understates its own error"
    else "no, below 1: the model overstates its own error"))
  cat(sprintf("  would kriging move the cell means?     %d of %d cells by > 1 SE\n",
              moved, nrow(ka)))
  cat(sprintf("     %s\n", if (moved == 0)
    "no: report the plain cell means, kriging buys nothing here" else
    "yes: the plain means are biased by where the samples fell inside each cell"))
  cat("  The interesting case is honest uncertainties AND large shifts: that is\n",
      "  when the extra machinery earns its place in the methods section.\n",
      sep = "")
}
