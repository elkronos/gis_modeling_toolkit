# =============================================================================
# 07  From a fit to a map, and where the map stops meaning anything
# =============================================================================
#   source(system.file("scripts", "07-surface-aoa.R", package = "spatialkit"))
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

pts <- tour_points()
bnd <- tour_boundary(pts)

# Train on the western half only. `slope` runs west to east in this fixture, so
# the eastern half is genuinely new ground in predictor space, not just on the
# map -- which is the situation the area of applicability exists to detect.
west <- pts[sf::st_coordinates(pts)[, 1] < 5e5 + 500, ]   # x runs 5e5 to 5e5 + 1000
ws_fit <- function(train_sf, ...) {
  new_spatial_fit(subclass = "ws_fit",
                  engine = stats::lm(z ~ elev + slope,
                                     data = sf::st_drop_geometry(train_sf)),
                  formula = z ~ elev + slope, response_var = "z",
                  predictor_vars = c("elev", "slope"), data_sf = train_sf)
}
predict.ws_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  as.numeric(stats::predict(object$engine, newdata = sf::st_drop_geometry(newdata)))
}
registerS3method("predict", "ws_fit", predict.ws_fit)

fit <- ws_fit(west)

step("07.1", "The default grid covers the training data, not the boundary")
# This catches people out. The grid is built from the bounding box of the data
# the model was fitted on. A `boundary` only clips that grid; it cannot extend
# it past where the training points reach.
cat(sprintf("  trained on %d of %d points, slope %.2f to %.2f\n",
            nrow(west), nrow(pts), min(west$slope), max(west$slope)))
srf_default <- predict_surface(fit, boundary = bnd, n_cells = 2000,
                               covariates = pts)
cat(sprintf("  boundary spans x = %.0f to %.0f\n",
            sf::st_bbox(bnd)["xmin"], sf::st_bbox(bnd)["xmax"]))
cat(sprintf("  grid spans     x = %.0f to %.0f  <- the training extent\n",
            min(sf::st_coordinates(srf_default)[, 1]),
            max(sf::st_coordinates(srf_default)[, 1])))
cat("  Hand it a `grid` of your own to predict anywhere else.\n")

step("07.2", "Predicting onto ground the model has not seen")
# A plain point grid over the whole boundary. `covariates` supplies the
# predictor values by nearest neighbour, so the layer you pass has to cover the
# area you are predicting on, not just the area you trained on.
grid <- sf::st_sf(geometry = sf::st_make_grid(bnd, n = c(60, 60), what = "centers"))
grid <- grid[lengths(sf::st_intersects(grid, bnd)) > 0L, ]
srf  <- predict_surface(fit, grid = grid, covariates = pts)
cat(sprintf("  %d prediction points, slope %.2f to %.2f (training stopped at %.2f)\n",
            nrow(srf), min(srf$slope), max(srf$slope), max(west$slope)))

step("07.3", "Which of those predictions are extrapolations")
aoa <- area_of_applicability(srf, model = fit)
print(aoa)
cat(sprintf("\n  %.0f%% of the map is inside the AOA. The rest is the model\n",
            100 * aoa$n_inside / (aoa$n_inside + aoa$n_outside)))
cat("  answering a question it was never asked. It will still return numbers.\n")

step("07.4", "Draw the surface and the mask")
if (!skip_without("ggplot2", "the surface maps")) {
  srf$inside <- aoa$aoa$AOA
  look_for("the predicted surface. Nothing about it looks wrong in the east, ",
           "which is the problem: an extrapolation looks exactly like a ",
           "prediction.")
  show_plot(
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = srf, ggplot2::aes(colour = .data$.pred), size = 0.8) +
      ggplot2::geom_sf(data = west, shape = 1, size = 0.7, colour = "grey20") +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(title = "Predicted z", subtitle = "circles: the training points") +
      ggplot2::theme_minimal(),
    "07-surface.png")

  look_for("the same map with the AOA mask. The boundary between the two ",
           "follows predictor space, not the edge of the training block.")
  show_plot(
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = srf, ggplot2::aes(colour = .data$inside), size = 0.8) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
                                   name = "inside AOA") +
      ggplot2::labs(title = "Area of applicability") +
      ggplot2::theme_minimal(),
    "07-aoa.png")

  look_for("the DI distribution: training points against prediction points, ",
           "with the threshold marked. A prediction set whose curve sits well ",
           "right of the training curve is extrapolating even where it passes.")
  show_plot(plot(aoa), "07-aoa-ecdf.png")
}

step("07.5", "Getting it out of R as a raster")
# `predict_surface()` returns POINTS, one per grid node, because that is what a
# model predicts at. Most GIS software wants a raster.
cat("  the surface is a",
    as.character(unique(sf::st_geometry_type(srf))), "layer\n")
out <- if (nzchar(tour_out)) tour_out else tempdir()
if (have("stars")) {
  # The grid was 60 x 60 over the boundary, so that is the cell size to ask for.
  bb <- sf::st_bbox(bnd)
  dx <- as.numeric(bb["xmax"] - bb["xmin"]) / 60
  dy <- as.numeric(bb["ymax"] - bb["ymin"]) / 60
  r  <- stars::st_rasterize(srf[".pred"], dx = dx, dy = dy)
  f <- file.path(out, "07-surface.tif")
  stars::write_stars(r, f)
  cat("  wrote", f, "\n")
} else {
  cat("  SKIPPED (GeoTIFF): install.packages(\"stars\"), or rasterise with\n",
      "  terra::rasterize() if you already have terra.\n", sep = "")
}
# The vector route, for anyone who would rather keep the cells:
f2 <- file.path(out, "07-surface.gpkg")
sf::st_write(srf[".pred"], f2, delete_dsn = TRUE, quiet = TRUE)
cat("  wrote", f2, " (GeoPackage, opens in QGIS)\n")

step("07.6", "Mask first, then summarise")
# Summary statistics over the whole surface include the extrapolated part.
inside <- aoa$aoa$AOA
cat(sprintf("  mean of the whole surface      %.3f\n", mean(srf$.pred)))
cat(sprintf("  mean inside the AOA only       %.3f\n", mean(srf$.pred[inside])))
cat(sprintf("  they differ by %.1f%% of the inside-AOA standard deviation\n",
            100 * abs(mean(srf$.pred) - mean(srf$.pred[inside])) /
              stats::sd(srf$.pred[inside])))
cat("  Do it in that order. Any number you quote from a surface -- a mean, a\n",
    "  total, a share above some threshold -- inherits whatever the model made\n",
    "  up outside its applicability.\n", sep = "")
