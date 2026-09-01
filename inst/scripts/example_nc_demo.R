#!/usr/bin/env Rscript
# =============================================================================
# spatialkit demo - North Carolina with synthetic data
# =============================================================================
#
# Author:  Justin Chase <jchase.msu@gmail.com>
# Package: spatialkit (https://github.com/elkronos/gis_modeling_toolkit)
# Licence: MIT (see the LICENSE file shipped with the package)
# Version: this script ships with the package; check
#          packageVersion("spatialkit") for the version you are running.
#
# A self-contained example showing how spatialkit handles different
# tessellations (Voronoi, hex grid, square grid, Delaunay triangles), cell
# level aggregation, spatially blocked cross-validation, and model fitting.
#
# The boundary comes from sf's built-in nc.shp - no external files needed.
#
# SETUP
#   1. Install the package:
#        remotes::install_github("elkronos/gis_modeling_toolkit")
#      or, from a local clone:
#        devtools::install("path/to/spatialkit")
#
#   2. Install what this script uses beyond the hard dependencies
#      (sf, dplyr, logger and digest install automatically):
#        install.packages(c("ggplot2", "patchwork"))  # plots and the panel
#        install.packages("ranger")                   # random forest backend
#        install.packages("gstat")                    # autocorrelation range
#        install.packages("geometry")                 # Delaunay triangulation
#        install.packages(c("sp", "GWmodel"))         # GWR backend (optional)
#      Anything missing is skipped with a message rather than failing.
#
#   3. Run it:
#        source(system.file("scripts", "example_nc_demo.R", package = "spatialkit"))
#
# OUTPUT
#   PNGs are written to a per-session temporary directory by default, so
#   sourcing this script never scatters files into your working directory.
#   The path is printed at the start and at the end.
#
#   To write somewhere permanent instead, set the output directory BEFORE
#   sourcing:
#        Sys.setenv(SPATIALKIT_DEMO_OUTPUT = "~/spatialkit-demo")
#        source(system.file("scripts", "example_nc_demo.R", package = "spatialkit"))
# =============================================================================

library(spatialkit)
library(sf)
library(dplyr)

has_ggplot  <- requireNamespace("ggplot2",  quietly = TRUE)
has_ranger  <- requireNamespace("ranger",   quietly = TRUE)
has_gstat   <- requireNamespace("gstat",    quietly = TRUE)
has_geom    <- requireNamespace("geometry", quietly = TRUE)
has_gwmodel <- requireNamespace("GWmodel",  quietly = TRUE) &&
               requireNamespace("sp",       quietly = TRUE)

if (!has_ggplot)
  stop("This demo draws maps and needs ggplot2. install.packages(\"ggplot2\")",
       call. = FALSE)
library(ggplot2)

set.seed(42)

# --- Output directory: tempdir() unless the user asked for somewhere else ----
out_dir <- Sys.getenv("SPATIALKIT_DEMO_OUTPUT", unset = "")
if (!nzchar(out_dir)) out_dir <- file.path(tempdir(), "spatialkit-demo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save_plot <- function(file, plot, width = 10, height = 6) {
  ggsave(file.path(out_dir, file), plot, width = width, height = height, dpi = 200)
  cat(sprintf("  wrote %s\n", file.path(out_dir, file)))
}
cat("Output directory:", out_dir, "\n")
cat("  (set SPATIALKIT_DEMO_OUTPUT before sourcing to change it)\n\n")

# =============================================================================
# 1. NORTH CAROLINA BOUNDARY
# =============================================================================
# sf ships nc.shp (100 county polygons). Dissolve to a state outline and
# project to NAD83 / NC State Plane (ftUS), EPSG:2264, so distances are planar
# rather than angular.  Projected does not mean metric: 2264's unit is the US
# SURVEY FOOT, so every distance, bandwidth and block size below - and the
# autocorrelation range printed in section 6 - is in feet, not metres.

nc_counties <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_boundary <- nc_counties |>
  st_union() |>
  st_transform(2264) |>
  st_as_sf()

cat("Loaded North Carolina boundary (EPSG:2264)\n")

# =============================================================================
# 2. SYNTHETIC OBSERVATIONS
# =============================================================================

n_points <- 300

pts_raw    <- st_sample(nc_boundary, size = n_points, type = "random")
pts_coords <- st_coordinates(pts_raw)
x_coords   <- pts_coords[, 1]
y_coords   <- pts_coords[, 2]

# Predictor: elevation (gradient west to east + noise)
elevation <- scale(x_coords)[, 1] * 500 + rnorm(n_points, 3000, 400)

# Predictor: population density (higher near two fake cities)
city1 <- c(1530000, 550000)   # Charlotte-ish in EPSG:2264
city2 <- c(2150000, 750000)   # Raleigh-ish
dist_to_city <- pmin(
  sqrt((x_coords - city1[1])^2 + (y_coords - city1[2])^2),
  sqrt((x_coords - city2[1])^2 + (y_coords - city2[2])^2)
)
pop_density <- pmax(exp(-dist_to_city / 400000) * 5000 + rnorm(n_points, 200, 100), 10)

# Response: predictor effects PLUS a spatial field neither predictor explains.
# That last term is what makes random cross-validation optimistic in section 6.
spatial_field <- 20 * sin(x_coords / 250000) * cos(y_coords / 250000)

y_response <- 50 +
  0.01  * elevation +
  0.005 * pop_density +
  spatial_field +
  rnorm(n_points, 0, 5)

points_sf <- st_sf(
  y           = y_response,
  elevation   = elevation,
  pop_density = pop_density,
  geometry    = pts_raw
)

cat(sprintf("Generated %d synthetic observation points\n", n_points))

response_var   <- "y"
predictor_vars <- c("elevation", "pop_density")

# =============================================================================
# 3. FOUR TESSELLATION TYPES
# =============================================================================

cat("\n--- Building tessellations ---\n")

# 3a. Voronoi from ~40 k-means seeds.  Cells are built around the SEEDS, not
# the observations: one cell per observation is a nearest-neighbour
# interpolation rather than an aggregation, and leaves no within-cell
# variation to compute a standard error from.
seeds <- get_voronoi_seeds(
  boundary      = nc_boundary,
  sample_points = points_sf,
  method        = "kmeans",
  n             = 40
)
tess_voronoi <- build_tessellation(
  seeds, boundary = nc_boundary,
  method = "voronoi", clip = TRUE, quiet = TRUE
)
cat(sprintf("  Voronoi: %d cells\n", nrow(tess_voronoi$cells)))

# 3b / 3c. Hex and square grids (~50 cells).  Both keep `cell_id`.
tess_hex <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "hex", approx_n_cells = 50, clip = TRUE, quiet = TRUE
)
cat(sprintf("  Hex grid: %d cells\n", nrow(tess_hex$cells)))

tess_square <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "square", approx_n_cells = 50, clip = TRUE, quiet = TRUE
)
cat(sprintf("  Square grid: %d cells\n", nrow(tess_square$cells)))

# 3d. Delaunay triangles.  Without 'geometry' this does NOT error: it falls
# back to sf::st_triangulate() on the point set and logs a warning.  Guard on
# the package so the two paths are not silently conflated.
tess_tri <- NULL
if (has_geom) {
  tess_tri <- build_tessellation(
    points_sf, boundary = nc_boundary,
    method = "triangles", clip = TRUE, quiet = TRUE
  )
  cat(sprintf("  Delaunay triangles: %d cells\n", nrow(tess_tri$cells)))
} else {
  cat("  Delaunay triangles: skipped (install.packages(\"geometry\"))\n")
}

# =============================================================================
# 4. CELL-LEVEL AGGREGATION AND CHOROPLETHS
# =============================================================================
# summarize_by_cell() does the aggregation: mean, sd and se for the response
# and every predictor, plus n and cell_weight.  There is no need to write the
# group_by()/summarise() by hand, and doing so loses the design-effect
# machinery shown below.

make_choropleth <- function(tess, boundary, points, title = NULL,
                            fill_col = "resp_mean_y",
                            legend_title = "Mean y", palette = "viridis") {
  cells    <- tess$cells
  assigned <- assign_features_to_polygons(points, cells, polygon_id_col = "cell_id")
  stats_df <- summarize_by_cell(assigned, response_var = "y", id_col = "cell_id")
  cells    <- left_join(cells,
                        as.data.frame(stats_df)[, c("cell_id", "resp_mean_y")],
                        by = "cell_id")

  plot_tessellation_map(
    tessellation_sf = cells,
    boundary        = boundary,
    fill_col        = fill_col,
    palette         = palette,
    tile_alpha      = 0.9,
    outline_col     = "white",
    outline_size    = 0.3,
    boundary_col    = "grey20",
    boundary_size   = 0.8,
    legend_title    = legend_title,
    title           = title,
    subtitle        = sprintf("%d cells  |  %d observations",
                              nrow(cells), nrow(points)),
    caption         = "spatialkit demo - synthetic North Carolina data"
  )
}

cat("\n--- Generating choropleth maps ---\n")

p_vor <- make_choropleth(tess_voronoi, nc_boundary, points_sf,
                         "Voronoi tessellation - mean response")
p_hex <- make_choropleth(tess_hex, nc_boundary, points_sf,
                         "Hexagonal grid - mean response")
p_sq  <- make_choropleth(tess_square, nc_boundary, points_sf,
                         "Square grid - mean response")

save_plot("01_voronoi_choropleth.png", p_vor)
save_plot("02_hex_choropleth.png",     p_hex)
save_plot("03_square_choropleth.png",  p_sq)

if (!is.null(tess_tri)) {
  save_plot("04_delaunay_choropleth.png",
            make_choropleth(tess_tri, nc_boundary, points_sf,
                            "Delaunay triangulation - mean response"))
}

# Design-effect-corrected standard errors.  The ..se_* columns are IID at the
# default deff = 1, which is anticonservative when points inside a cell are
# spatially correlated.
assigned_vor <- assign_features_to_polygons(points_sf, tess_voronoi$cells,
                                            polygon_id_col = "cell_id")
cells_iid  <- summarize_by_cell(assigned_vor, response_var = response_var,
                                predictor_vars = predictor_vars,
                                id_col = "cell_id")
cells_kish <- summarize_by_cell(assigned_vor, response_var = response_var,
                                predictor_vars = predictor_vars,
                                id_col = "cell_id", deff = "kish")
deff <- attr(cells_kish, "deff_applied")
cat(sprintf("\n  Design effect: method = %s | ICC(response) = %.3f | median deff = %.2f\n",
            deff$method, deff$icc_resp, stats::median(deff$deff, na.rm = TRUE)))
cat(sprintf("  Response SE inflation: median %.2fx\n",
            stats::median(cells_kish[["..se_resp_y"]] / cells_iid[["..se_resp_y"]],
                          na.rm = TRUE)))

# =============================================================================
# 5. SIDE-BY-SIDE COMPARISON PANEL
# =============================================================================

cat("\n--- Building comparison panel ---\n")
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_panel <- (p_vor | p_hex | p_sq) +
    plot_annotation(
      title    = "spatialkit - tessellation comparison (North Carolina)",
      subtitle = sprintf("%d observations  |  3 tessellation methods", n_points),
      theme    = theme(
        plot.title    = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11, colour = "grey40")
      )
    )
  save_plot("05_comparison_panel.png", p_panel, width = 18, height = 6)
} else {
  cat("  Install 'patchwork' for the side-by-side panel: install.packages('patchwork')\n")
}

# =============================================================================
# 6. SPATIAL CROSS-VALIDATION
# =============================================================================
# The point of the package: the fold scheme, not the model, decides whether
# the reported score describes the task you actually have.

cat("\n--- Spatial cross-validation ---\n")

if (has_gstat) {
  sac <- estimate_sac_range(points_sf, response_var = response_var,
                            predictor_vars = predictor_vars)
  if (is.na(sac)) {
    cat("  Autocorrelation range: not identified (",
        attr(sac, "rejected_reason"), ")\n", sep = "")
  } else {
    cat(sprintf("  Autocorrelation range: %.0f ft\n", as.numeric(sac)))
  }
} else {
  cat("  Skipping estimate_sac_range(): install.packages(\"gstat\")\n")
}

folds_random  <- make_folds(points_sf, k = 5, method = "random_kfold", seed = 42)
folds_blocked <- make_folds(points_sf, k = 5, method = "block_kfold",  seed = 42)
cat(sprintf("  Folds built: random k = %d, blocked k = %d\n",
            folds_random$k, folds_blocked$k))

if (requireNamespace("patchwork", quietly = TRUE)) {
  save_plot("06_folds_random_vs_blocked.png",
            plot_folds(folds_random,  points_sf, boundary = nc_boundary) +
            plot_folds(folds_blocked, points_sf, boundary = nc_boundary),
            width = 14, height = 6)
}

if (has_ranger) {
  # include_coords = TRUE hands the forest the coordinates, so it can memorise
  # the training surface.  That is the failure mode the default (FALSE) exists
  # to prevent, and the one random folds cannot see.
  rf_args <- list(include_coords = TRUE, num_trees = 300, seed = 1)
  cv_random  <- do.call(cv_rf, c(list(points_sf, response_var, predictor_vars,
                                      folds = folds_random),  rf_args))
  cv_blocked <- do.call(cv_rf, c(list(points_sf, response_var, predictor_vars,
                                      folds = folds_blocked), rf_args))
  cat(sprintf("  random_kfold : R2 = %.3f  RMSE = %.3f\n",
              cv_random$overall$R2, cv_random$overall$RMSE))
  cat(sprintf("  block_kfold  : R2 = %.3f  RMSE = %.3f   <- the one to report\n",
              cv_blocked$overall$R2, cv_blocked$overall$RMSE))
} else {
  cat("  Skipping the CV contrast: install.packages(\"ranger\")\n")
}

# =============================================================================
# 7. FIT A MODEL AND MAP ITS RESIDUALS
# =============================================================================

cat("\n--- Fitting models ---\n")

rf_fit <- NULL
if (has_ranger) {
  rf_fit <- fit_rf_model(points_sf, response_var = response_var,
                         predictor_vars = predictor_vars)
  cat("  Random forest fitted.\n")
  # NOTE: fitted() on an rf_fit returns OUT-OF-BAG predictions, so summary()
  # reports out-of-bag metrics, not in-sample ones, and is not comparable with
  # a gwr_fit or bayesian_fit.  coef() on an rf_fit errors by design - a forest
  # has no coefficients; use $info$importance.
  print(rf_fit$info$importance)
  mi <- residual_morans_i(rf_fit)
  if (!is.null(mi))
    cat(sprintf("  Residual Moran's I = %.4f (z = %.2f, p = %.3g)\n",
                mi$observed, mi$z, mi$p_value))
  save_plot("07_rf_residuals.png", plot(rf_fit, type = "residuals"))
} else {
  cat("  Skipping the random forest: install.packages(\"ranger\")\n")
}

gwr_fit <- NULL
if (has_gwmodel) {
  gwr_fit <- tryCatch(
    fit_gwr_model(data_sf = points_sf, response_var = response_var,
                  predictor_vars = predictor_vars,
                  adaptive = TRUE, kernel = "bisquare"),
    error = function(e) {
      cat("  GWR skipped:", conditionMessage(e), "\n"); NULL
    }
  )
  if (!is.null(gwr_fit)) {
    met <- model_metrics(gwr_fit)     # a spatial_fit stores $info, not metrics
    cat(sprintf("  GWR: bandwidth %.1f | in-sample R2 %.3f | RMSE %.3f\n",
                gwr_fit$info$bandwidth, met$R2, met$RMSE))
  }
} else {
  cat("  Skipping GWR: install.packages(c(\"sp\", \"GWmodel\"))\n")
}

# Aggregate absolute model error onto each tessellation.
fit_for_error <- if (!is.null(rf_fit)) rf_fit else gwr_fit
if (!is.null(fit_for_error)) {
  cat("\n--- Mapping model error per tessellation ---\n")
  pts_err <- points_sf
  pts_err$abs_error <- abs(as.numeric(residuals(fit_for_error)))

  for (tess_info in list(
    list(tess = tess_voronoi, name = "voronoi", label = "Voronoi"),
    list(tess = tess_hex,     name = "hex",     label = "Hex grid"),
    list(tess = tess_square,  name = "square",  label = "Square grid")
  )) {
    cells <- tess_info$tess$cells
    asgn  <- assign_features_to_polygons(pts_err, cells, polygon_id_col = "cell_id")
    errs  <- summarize_by_cell(asgn, response_var = "abs_error", id_col = "cell_id")
    cells <- left_join(cells,
                       as.data.frame(errs)[, c("cell_id", "resp_mean_abs_error")],
                       by = "cell_id")

    p <- plot_tessellation_map(
      tessellation_sf = cells,
      boundary        = nc_boundary,
      fill_col        = "resp_mean_abs_error",
      palette         = "magma",
      tile_alpha      = 0.9,
      outline_col     = "white",
      outline_size    = 0.3,
      boundary_col    = "grey20",
      boundary_size   = 0.8,
      legend_title    = "Mean |residual|",
      title           = sprintf("Model residuals - %s", tess_info$label),
      subtitle        = sprintf("Absolute residual aggregated to %d cells",
                                nrow(cells)),
      caption         = "spatialkit demo - synthetic North Carolina data"
    )
    save_plot(sprintf("08_error_%s.png", tess_info$name), p)
  }
}

# =============================================================================
# 8. PREDICTION SURFACE AND AREA OF APPLICABILITY
# =============================================================================

if (!is.null(fit_for_error)) {
  cat("\n--- Prediction surface ---\n")
  surf <- predict_surface(fit_for_error, n_cells = 3000,
                          covariates = points_sf, boundary = nc_boundary)
  aoa  <- area_of_applicability(surf, model = fit_for_error, folds = folds_blocked)
  cat(sprintf("  %d grid cells | inside AOA: %d | outside: %d | undetermined: %d\n",
              nrow(surf), aoa$n_inside, aoa$n_outside, aoa$n_na))

  # AOA is NA where a predictor was missing or non-finite, and `!NA` is NA,
  # which R silently skips in a subscripted assignment.  Test for TRUE.
  inside <- aoa$aoa$AOA %in% TRUE
  surf$.pred_masked <- ifelse(inside, surf$.pred, NA_real_)

  p_surf <- ggplot() +
    geom_sf(data = surf, aes(colour = .pred_masked), size = 0.6) +
    geom_sf(data = nc_boundary, fill = NA, colour = "grey20") +
    scale_colour_viridis_c(name = "Predicted y", na.value = "grey85") +
    theme_void() +
    ggtitle("Predicted surface (extrapolations blanked out)")
  save_plot("09_prediction_surface.png", p_surf)
}

# =============================================================================
# DONE
# =============================================================================

cat("\n==============================================================\n")
cat("  All outputs written to:\n    ", out_dir, "\n", sep = "")
print(sort(basename(list.files(out_dir, pattern = "\\.png$", full.names = TRUE))))
cat("==============================================================\n")
