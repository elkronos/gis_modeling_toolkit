# =============================================================================
# 09  Geographically weighted regression: coefficients that vary over space
# =============================================================================
#   source(system.file("scripts", "09-gwr.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

if (skip_without(c("GWmodel", "sp"), "the whole of script 09")) {
  cat("  GWR needs GWmodel, which needs sp. Nothing else in the tour does.\n")
} else {

pts <- tour_points()

step("09.1", "Fit one")
# GWR fits a separate regression at every observation, weighting the others by
# distance. `adaptive = TRUE` (the default) sets the bandwidth in NUMBER OF
# NEIGHBOURS rather than distance, which is what you want when sampling density
# varies across the study area.
gw <- fit_gwr_model(pts, "z", c("elev", "slope"))
print(summary(gw))
cat(sprintf("\n  bandwidth: %s %s, kernel %s\n", gw$info$bandwidth,
            if (isTRUE(gw$info$adaptive)) "nearest neighbours" else "CRS units",
            gw$info$kernel))
if (isTRUE(gw$info$bandwidth_is_fallback))
  cat("  (that bandwidth is a fallback: the search did not converge)\n")

step("09.2", "The output is a coefficient per point")
cf <- coef(gw)
cat(sprintf("  %d rows x %d coefficients\n", nrow(cf), ncol(cf)))
print(round(apply(cf, 2, function(v) stats::quantile(v, c(0, 0.5, 1), na.rm = TRUE)), 3))
cat("  A coefficient whose sign changes across the map is the finding. It is\n",
    "  also the thing to be most careful about: local collinearity produces\n",
    "  sign flips that mean nothing.\n", sep = "")

step("09.3", "Check for local collinearity before believing any of it")
# Two predictors that are globally independent can be nearly collinear inside a
# small neighbourhood. Where that happens the local coefficients are unstable
# and their signs are arbitrary.
cat(sprintf("  local condition number above threshold: %d of %d fits\n",
            gw$info$n_local_collinear, nrow(cf)))
cat(sprintf("  locally singular: %d, non-finite coefficients: %d\n",
            gw$info$n_local_singular, sum(gw$info$nonfinite_coef)))
if (gw$info$n_local_collinear > 0)
  cat("  Widen the bandwidth or drop a predictor before reading the maps.\n")

step("09.4", "Map the coefficients")
if (!skip_without("ggplot2", "the coefficient maps")) {
  look_for("whether a coefficient surface has structure or looks like static. ",
           "Smooth patches are a spatially varying relationship; salt and ",
           "pepper is the bandwidth being too narrow.")
  show_plot(plot(gw, type = "coefficients", term = "elev"), "09-coef-elev.png")
  look_for("the same for slope. Compare the two: a predictor whose coefficient ",
           "barely moves does not need GWR at all.")
  show_plot(plot(gw, type = "coefficients", term = "slope"), "09-coef-slope.png")
}

step("09.5", "Is it actually better than a global fit?")
# GWR fits one regression per point, so its in-sample error is nearly always the
# lower of the two. The comparison that means something is cross-validated, on
# blocked folds.
bnd   <- tour_boundary(pts)
folds <- make_folds(pts, k = 5, method = "block_kfold", seed = 1,
                    boundary = bnd, block_size = 300)
gl_fit <- function(train_sf, ...) {
  new_spatial_fit(subclass = "gl_fit",
                  engine = stats::lm(z ~ elev + slope,
                                     data = sf::st_drop_geometry(train_sf)),
                  formula = z ~ elev + slope, response_var = "z",
                  predictor_vars = c("elev", "slope"), data_sf = train_sf)
}
predict.gl_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  as.numeric(stats::predict(object$engine, newdata = sf::st_drop_geometry(newdata)))
}
# cv_spatial() only ever needs `predict`. model_metrics() on the fit ITSELF
# needs `fitted` as well, and says so rather than guessing.
fitted.gl_fit <- function(object, ...) stats::fitted(object$engine)
registerS3method("predict", "gl_fit", predict.gl_fit)
registerS3method("fitted",  "gl_fit", fitted.gl_fit)

cv_g <- cv_gwr(pts, "z", c("elev", "slope"), folds = folds)
cv_l <- cv_spatial(pts, "z", c("elev", "slope"), fit_fn = gl_fit, folds = folds)
cat(sprintf("  in-sample:  GWR %.3f   global lm %.3f\n",
            model_metrics(gw)$RMSE, model_metrics(gl_fit(pts))$RMSE))
cat(sprintf("  blocked CV: GWR %.3f   global lm %.3f\n",
            cv_g$overall$RMSE[1L], cv_l$overall$RMSE[1L]))
gap <- cv_g$overall$RMSE[1L] / cv_l$overall$RMSE[1L] - 1
if (abs(gap) < 0.02) {
  cat("  Out of sample the two are within a couple of percent of each other, so\n",
      "  the whole in-sample gain was local fitting that did not generalise. On\n",
      "  this map the relationship does not vary enough to pay for one\n",
      "  regression per point.\n", sep = "")
} else if (gap < 0) {
  cat(sprintf("  GWR is %.0f%% better out of sample, so it earns its place.\n", -100 * gap))
} else {
  cat(sprintf("  GWR is %.0f%% worse out of sample: the local fits were memorising.\n",
              100 * gap))
}

step("09.6", "Picking predictors for a GWR")
# gwr_model_selection() ranks subsets by AICc. It is fast because it reuses one
# bandwidth for every candidate model, and that is also its limitation.
ms <- gwr_model_selection(pts, "z", c("elev", "slope", "noise"), quiet = TRUE)
print(ms)
if (!skip_without("ggplot2", "the model selection plot")) {
  look_for("each dot is one model, the line joins the best at each size, and ",
           "the coloured dot is the winner. A line that rises means the extra ",
           "predictors only cost you; a flat stretch means those sizes are ",
           "indistinguishable, so take the smallest.")
  show_plot(plot(ms), "09-model-selection.png")
}
cat("  The ranking is in-sample and all models share one bandwidth, so treat it\n",
    "  as a shortlist. Script 08 scores a shortlist properly.\n", sep = "")

step("09.7", "Before you read a coefficient map as a causal map")
# A GWR coefficient surface describes how the local regression had to bend to
# fit the data there. Omitted variables, edge effects and local collinearity all
# show up as spatial variation in the coefficients, and none of them is the
# effect of the predictor changing across the region.
cat("  Before writing 'the effect of elevation is stronger in the north':\n")
cat("   - check 09.3 for local collinearity\n")
cat("   - check whether the pattern survives a wider bandwidth\n")
cat("   - check whether an omitted predictor has the same spatial pattern\n")
cat("   - check the edges, where the kernel is one-sided by construction\n")

}
