#!/usr/bin/env Rscript
# =============================================================================
# Worked example: North Carolina
# =============================================================================
#
# Author:  Justin Chase <jchase.msu@gmail.com>
# Package: spatialkit (https://github.com/elkronos/gis_modeling_toolkit)
# Licence: MIT (see the LICENSE file shipped with the package)
# Version: this script ships with the package; packageVersion("spatialkit")
#          gives the version you are running.
#
# The steps of vignette("spatialkit_nc_demo") as one script, in the same
# order: a simulated field over North Carolina, four tessellations, cell means
# and their standard errors, a random forest and a GWR fit, random folds
# against blocks, and a predicted surface with its area of applicability. The
# article explains each step; this script prints the numbers and writes the
# figures. The state outline comes from the nc.shp file bundled with sf, so
# the script needs no other files.
#
# SETUP
#   1. Install the package:
#        install.packages("spatialkit")
#      or the development version from GitHub:
#        remotes::install_github("elkronos/gis_modeling_toolkit")
#
#   2. Install the optional packages the script uses:
#        install.packages(c("ggplot2", "patchwork"))  # figures and panels
#        install.packages("geometry")                 # Delaunay triangles
#        install.packages("ranger")                   # random forest
#        install.packages("sp")                       # GWR
#        install.packages("gstat")                    # the variogram
#      The script stops without ggplot2. Without one of the others it skips
#      the steps that need it and says so; without patchwork it writes each
#      map of a panel to a file of its own.
#
#   3. Run it:
#        source(system.file("scripts", "example_nc_demo.R", package = "spatialkit"))
#
# OUTPUT
#   PNGs go to a temporary folder by default, so sourcing this script writes
#   nothing to your working directory. The script prints the folder at the
#   start and at the end. To write somewhere permanent, set the folder before
#   sourcing:
#        Sys.setenv(SPATIALKIT_DEMO_OUTPUT = "~/spatialkit-demo")
#        source(system.file("scripts", "example_nc_demo.R", package = "spatialkit"))
# =============================================================================

has_ggplot    <- requireNamespace("ggplot2",   quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
has_geom      <- requireNamespace("geometry",  quietly = TRUE)
has_ranger    <- requireNamespace("ranger",    quietly = TRUE)
has_gstat     <- requireNamespace("gstat",     quietly = TRUE)
has_gwr       <- requireNamespace("sp",        quietly = TRUE)

if (!has_ggplot)
  stop("This script draws maps and needs ggplot2: install.packages(\"ggplot2\")",
       call. = FALSE)

library(spatialkit)
library(sf)
library(ggplot2)

set.seed(42)

# Output folder: a temporary one unless SPATIALKIT_DEMO_OUTPUT names another.
out_dir <- Sys.getenv("SPATIALKIT_DEMO_OUTPUT", unset = "")
if (!nzchar(out_dir)) out_dir <- file.path(tempdir(), "spatialkit-demo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save_plot <- function(file, plot, width = 7, height = 3) {
  ggsave(file.path(out_dir, file), plot, width = width, height = height, dpi = 200)
  cat(sprintf("  wrote %s\n", file.path(out_dir, file)))
}
skipping <- function(what, pkg)
  cat(sprintf("  Skipping %s: install.packages(\"%s\")\n", what, pkg))
cat("Output folder:", out_dir, "\n")
cat("  (set SPATIALKIT_DEMO_OUTPUT before sourcing to change it)\n")

# =============================================================================
# 1. The state outline, in feet
# =============================================================================
# The county boundaries shipped with sf, dissolved into one state outline and
# projected to NAD83 / North Carolina (ftUS), EPSG:2264, so that distances are
# planar. Projected does not mean metric: every distance and block size below
# is in US survey feet, the unit the CRS names (the GWR bandwidth is a count
# of neighbours).

cat("\n--- The state outline ---\n")
nc_counties <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_boundary <- nc_counties |>
  st_union() |>
  st_transform(2264) |>
  st_as_sf()
cat("  North Carolina, EPSG:2264 (US survey feet)\n")

# =============================================================================
# 2. A simulated field with known structure
# =============================================================================
# Three hundred points inside the outline carry two predictors and a response.
# `elevation` rises from west to east, plus noise; `pop_density` decays with
# the distance from two invented cities; `y` is the sum of the two predictors'
# effects, a spatial field that neither predictor carries, and noise. The
# field is structure that nearby observations share, so a model's error
# depends on how far from the data it predicts.

cat("\n--- A simulated field ---\n")
n_points <- 300

pts_raw    <- st_sample(nc_boundary, size = n_points, type = "random")
pts_coords <- st_coordinates(pts_raw)
x_coords   <- pts_coords[, 1]
y_coords   <- pts_coords[, 2]

elevation <- scale(x_coords)[, 1] * 500 + rnorm(n_points, 3000, 400)

city1 <- c(1530000, 550000)   # near Charlotte in EPSG:2264
city2 <- c(2150000, 750000)   # near Raleigh
dist_to_city <- pmin(
  sqrt((x_coords - city1[1])^2 + (y_coords - city1[2])^2),
  sqrt((x_coords - city2[1])^2 + (y_coords - city2[2])^2)
)
pop_density <- pmax(exp(-dist_to_city / 400000) * 5000 + rnorm(n_points, 200, 100), 10)

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
cat(sprintf("  %d simulated observations\n", n_points))

save_plot("01_observations.png",
          ggplot() +
            geom_sf(data = nc_boundary, fill = "grey95", colour = "black") +
            geom_sf(data = points_sf, aes(colour = y), size = 1.1) +
            scale_colour_viridis_c(name = "Response (y)") +
            theme_void() +
            ggtitle("Simulated observations"))

# =============================================================================
# 3. Four tessellations
# =============================================================================
# vignette("cells") says which kind of cell suits which purpose.

cat("\n--- Four tessellations ---\n")

# Voronoi cells grow around seed points, not around the observations: one
# cell per observation would be a nearest-neighbour interpolation, with no
# spread inside a cell to compute a standard error from. get_voronoi_seeds()
# clusters the observations first, each weighted by its local density, so
# cells are small where the observations are dense.
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
# A k-means centre is a mean of observations, so it can lie outside an outline
# that is not convex. Its cell, the part of the state nearer to it than to any
# other seed, is still built inside the outline, and tess_voronoi$index holds
# the cell_id each seed generates.
outside <- lengths(st_intersects(seeds, nc_boundary)) == 0
cat(sprintf("  seeds outside the state outline: %d of %d; their cells: %s\n",
            sum(outside), nrow(seeds),
            if (any(outside)) paste(tess_voronoi$index[outside], collapse = ", ")
            else "none"))

# Hex and square grids of about 50 cells. All four methods return cells with
# a cell_id column; the grid methods keep poly_id beside it.
tess_hex <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "hex", approx_n_cells = 50, clip = TRUE, quiet = TRUE
)
tess_square <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "square", approx_n_cells = 50, clip = TRUE, quiet = TRUE
)

# Delaunay triangles. Without geometry, method = "triangles" falls back to
# sf::st_triangulate() with a logged warning; the guard below keeps the script
# on the geometry::delaunayn() path the article shows. The corners of the
# triangles are the observations, so most triangles hold none.
tess_tri <- NULL
if (has_geom) {
  tess_tri <- build_tessellation(
    points_sf, boundary = nc_boundary,
    method = "triangles", clip = TRUE, quiet = TRUE
  )
  tri_assigned <- assign_features_to_polygons(points_sf, tess_tri$cells,
                                              polygon_id_col = "cell_id")
  cat(sprintf("  %d of %d triangles hold a point\n",
              length(unique(na.omit(tri_assigned$cell_id))), nrow(tess_tri$cells)))
} else {
  skipping("the Delaunay triangles", "geometry")
}

cat("  cells per method:\n")
print(sapply(Filter(Negate(is.null),
                    list(Voronoi = tess_voronoi, Hex = tess_hex,
                         Square = tess_square, Triangles = tess_tri)),
             function(tess) nrow(tess$cells)))

# =============================================================================
# 4. Cell means and their standard errors
# =============================================================================
# summarize_by_cell() reduces the points to one row per cell: the count n; the
# mean, standard deviation and standard error of the response and of every
# predictor; a cell_weight; and beside each mean a pooled value, the cell's
# mean pulled toward the mean of the observations in the 8 nearest cells
# (nearest by the centroid of their observations).

cat("\n--- Cell means and their standard errors ---\n")
assigned <- assign_features_to_polygons(points_sf, tess_voronoi$cells,
                                        polygon_id_col = "cell_id")

cell_stats <- summarize_by_cell(
  assigned,
  response_var   = "y",
  predictor_vars = c("elevation", "pop_density"),
  id_col         = "cell_id"
)
print(head(as.data.frame(cell_stats)[, c("cell_id", "n", "resp_mean_y",
                                         "..se_resp_y", "..pooled_resp_y")]))

# At the default deff = 1 the ..se_* columns treat the points in a cell as
# independent. That is the standard error of each cell's own mean, the value a
# map of the cells shows, when the points are spread through the cell. As a
# standard error for the population mean it rests on deff = 1, which is
# anticonservative when points inside a cell are spatially correlated;
# deff = "kish" widens it by a design effect estimated from that correlation.
cell_kish <- summarize_by_cell(
  assigned,
  response_var   = "y",
  predictor_vars = c("elevation", "pop_density"),
  id_col         = "cell_id",
  deff           = "kish"
)
deff <- attr(cell_kish, "deff_applied")
cat(sprintf("  method = %s | ICC(response) = %.3f | median design effect = %.2f\n",
            deff$method, deff$icc_resp, stats::median(deff$deff, na.rm = TRUE)))
# How much the Kish standard error of the response exceeds the default one.
# An estimated design effect also corrects the within-cell standard deviation,
# so the ratio can exceed the square root of the design effect.
inflation <- cell_kish$`..se_resp_y` / cell_stats$`..se_resp_y`
cat(sprintf("  Kish standard error over the default: median %.2f (square root of the median design effect %.2f)\n",
            stats::median(inflation, na.rm = TRUE),
            sqrt(stats::median(deff$deff, na.rm = TRUE))))

# Maps of the cell means. summarize_by_cell(cells_sf = ) returns the summary
# on the cells themselves, ready for plot_tessellation_map().
on_cells <- function(tess) {
  assigned <- assign_features_to_polygons(points_sf, tess$cells,
                                          polygon_id_col = "cell_id")
  summarize_by_cell(assigned, response_var = "y", id_col = "cell_id",
                    cells_sf = tess$cells)
}
# plot_tessellation_map() sets a fill scale of its own; the one added here
# replaces it, and suppressMessages() keeps ggplot2's note on that quiet.
map_means <- function(cells, label, limits = NULL) suppressMessages(
  plot_tessellation_map(
    tessellation_sf = cells,
    boundary        = nc_boundary,
    fill_col        = "resp_mean_y",
    outline_col     = "white",
    outline_size    = 0.3,
    boundary_col    = "grey20",
    boundary_size   = 0.8,
    title           = sprintf("%s, %d cells", label, nrow(cells))
  ) +
    scale_fill_viridis_c(name = "Mean y", limits = limits, na.value = "grey60")
)

# The Voronoi, hex and square maps share one colour scale, so a colour means
# the same value on every map. Grey cells hold no observation.
means <- lapply(list(Voronoi = tess_voronoi, Hex = tess_hex, Square = tess_square),
                on_cells)
lims  <- range(unlist(lapply(means, `[[`, "resp_mean_y")), na.rm = TRUE)
if (has_patchwork) {
  library(patchwork)
  save_plot("02_cell_means.png",
            (map_means(means$Voronoi, "Voronoi", lims) /
             map_means(means$Hex,     "Hex",     lims) /
             map_means(means$Square,  "Square",  lims)) +
              plot_layout(guides = "collect"),
            height = 7)
} else {
  for (nm in names(means))
    save_plot(sprintf("02_cell_means_%s.png", tolower(nm)),
              map_means(means[[nm]], nm, lims))
}
# The triangle map, on its own colour scale.
if (!is.null(tess_tri))
  save_plot("03_cell_means_triangles.png",
            map_means(on_cells(tess_tri), "Delaunay triangles"))

# =============================================================================
# 5. Models fitted on the points
# =============================================================================
# The models are fitted on the 300 points, not on the cell means. The cells can
# still serve as cross-validation blocks:
#   make_folds(points_sf, k = 5, method = "block_kfold", blocks = tess_voronoi$cells)

cat("\n--- Models fitted on the points ---\n")
rf_fit <- NULL
if (has_ranger) {
  rf_fit <- fit_rf_model(points_sf, response_var = "y",
                         predictor_vars = c("elevation", "pop_density"))
  print(rf_fit)
  # summary() on an rf_fit reports out-of-bag metrics, because fitted() on an
  # rf_fit returns out-of-bag predictions; a gwr_fit or bayesian_fit reports
  # in-sample ones, so the two are not comparable. A forest has no
  # coefficients, so coef() on an rf_fit errors; the importance is in $info.
  print(summary(rf_fit))
  print(rf_fit$info$importance)

  # Structure in the residual map is spatial autocorrelation the model has not
  # captured, and residual_morans_i() puts a number on it.
  save_plot("04_rf_residuals.png", plot(rf_fit, type = "residuals"), height = 3.4)
  mi <- residual_morans_i(rf_fit)
  cat(sprintf("  Moran's I of the forest's residuals = %.4f (z = %.2f, p = %.3g)\n",
              mi$observed, mi$z, mi$p_value))
} else {
  skipping("the random forest and every step that uses it", "ranger")
}

gwr_fit <- NULL
if (has_gwr) {
  gwr_fit <- fit_gwr_model(
    data_sf        = points_sf,
    response_var   = "y",
    predictor_vars = c("elevation", "pop_density"),
    adaptive       = TRUE,
    kernel         = "bisquare"
  )
  gwr_met <- model_metrics(gwr_fit)
  # adaptive = TRUE: the bandwidth is a number of nearest neighbours, not feet.
  cat(sprintf("  GWR bandwidth: %.0f neighbours | in-sample R2: %.3f | RMSE: %.3f\n",
              gwr_fit$info$bandwidth, gwr_met$R2, gwr_met$RMSE))
  # GWR fits one weighted regression at each location, so each coefficient
  # has a value at every point.
  save_plot("05_gwr_elevation.png",
            plot(gwr_fit, type = "coefficients", term = "elevation"), height = 3.4)
  # cn, the condition index that includes the intercept, is above 30 where the
  # local intercept is poorly determined (here because elevation sits far from
  # zero). The slopes are judged by cn_slopes, which centres the predictors in
  # each window, and count as collinear above 30.
  lc <- gwr_fit$info$local_collinearity
  cat(sprintf("  windows with a poorly determined intercept (cn above 30): %d of %d; largest cn_slopes %.1f\n",
              sum(lc$cn > 30, na.rm = TRUE), nrow(lc),
              max(lc$cn_slopes, na.rm = TRUE)))
} else {
  skipping("GWR", "sp")
}

# =============================================================================
# 6. Scoring the model on spatial folds
# =============================================================================
# A random fold leaves each held-out point's neighbours in the training data; a
# blocked fold holds out whole blocks of the state.

cat("\n--- Scoring the model on spatial folds ---\n")

# estimate_sac_range() fits a variogram to the residuals of y on the
# predictors and returns the autocorrelation range in CRS units, here US
# survey feet. It returns NA when the variogram never levels off within the
# distances it was fitted over.
if (has_gstat) {
  sac <- estimate_sac_range(points_sf, response_var = "y",
                            predictor_vars = c("elevation", "pop_density"))
  if (is.na(sac)) {
    cat("  autocorrelation range: not identified (",
        attr(sac, "rejected_reason"), ")\n", sep = "")
  } else {
    cat(sprintf("  autocorrelation range: %.0f ft\n", as.numeric(sac)))
  }
  save_plot("06_variogram.png", plot(sac), height = 3.6)
} else {
  skipping("the autocorrelation range", "gstat")
}

# With no range identified, the block size comes from the simulation. The
# field changes sign about every 785,000 ft (pi times 250,000), and 400,000 ft is
# about half that, the distance from a zero of the wave to its peak. On your
# own data, where the field is unknown, cv_block_size_sweep() gives the
# cross-validated error at each block size
# (vignette("spatial-cross-validation")).
folds_random  <- make_folds(points_sf, k = 5, method = "random_kfold", seed = 42)
folds_blocked <- make_folds(points_sf, k = 5, method = "block_kfold",
                            block_size = 400000, seed = 42)
cat(sprintf("  folds: random k = %d, blocked k = %d over %d blocks\n",
            folds_random$k, folds_blocked$k, folds_blocked$params$blocks_used))

if (has_patchwork) {
  save_plot("07_folds.png",
            (plot_folds(folds_random,  points_sf, boundary = nc_boundary) /
             plot_folds(folds_blocked, points_sf, boundary = nc_boundary)) +
              plot_layout(guides = "collect"),
            height = 6)
} else {
  save_plot("07_folds_random.png",
            plot_folds(folds_random,  points_sf, boundary = nc_boundary))
  save_plot("07_folds_blocked.png",
            plot_folds(folds_blocked, points_sf, boundary = nc_boundary))
}

if (has_ranger) {
  # include_coords = TRUE gives the forest the coordinates, which lets it fit
  # the training surface by location, so its error grows quickly with the
  # distance from the data. The same predictors and settings on both fold
  # schemes show how much.
  cv_random  <- cv_rf(points_sf, "y", c("elevation", "pop_density"),
                      folds = folds_random, include_coords = TRUE,
                      num_trees = 300, seed = 1)
  cv_blocked <- cv_rf(points_sf, "y", c("elevation", "pop_density"),
                      folds = folds_blocked, include_coords = TRUE,
                      num_trees = 300, seed = 1)
  cat(sprintf("  random_kfold : R2 = %.3f  RMSE = %.3f\n",
              cv_random$overall$R2, cv_random$overall$RMSE))
  cat(sprintf("  block_kfold  : R2 = %.3f  RMSE = %.3f\n",
              cv_blocked$overall$R2, cv_blocked$overall$RMSE))
}

# Which of the two to report depends on where the map will predict. The
# surface in section 7 covers the state these points were sampled across at
# random. fold_separation() compares the distances from such a map's locations
# to the data with the held-out distances of each scheme; map_pts, a regular
# sample of the state, stands in for that map's grid.
map_pts <- st_sample(nc_boundary, size = 2000, type = "regular")
reading <- function(f)
  attr(fold_separation(f, points_sf, prediction_points = map_pts), "target")$reading
cat(sprintf("  against a map of the state: random folds read \"%s\", blocks \"%s\"\n",
            reading(folds_random), reading(folds_blocked)))
cat("  (\"matched\": the held-out points sit about as far from their training\n",
    "  data as the map's locations sit from the data; \"farther\": they sit\n",
    "  farther, so that score tends to read high for this map)\n", sep = "")

# compare_models_cv() scores each backend on the folds passed in, so the
# scores are comparable; a backend whose package is missing is dropped with a
# message. The comparison uses the random folds, the scheme fold_separation()
# reads as "matched" for this map.
if (has_ranger) {
  cmp <- compare_models_cv(
    points_sf, "y", c("elevation", "pop_density"),
    models   = c("RF", "GWR"),
    folds    = folds_random,
    rf_args  = list(num_trees = 300)
  )
  print(cmp$overall[, c("model", "RMSE", "MAE", "R2")], row.names = FALSE)
}

# =============================================================================
# 7. A surface, and where it applies
# =============================================================================
# predict_surface() builds a regular grid over the training extent, takes each
# grid point's covariates from the nearest observation, clips the grid to the
# outline and predicts. The surface uses rf_fit, the forest on the two
# predictors alone.

if (has_ranger) {
  cat("\n--- A surface, and where it applies ---\n")
  surf <- predict_surface(rf_fit, n_cells = 3000,
                          covariates = points_sf, boundary = nc_boundary)

  # area_of_applicability() marks where the cross-validated score applies.
  # Pass it the folds the score came from, here the random ones.
  aoa <- area_of_applicability(surf, model = rf_fit, folds = folds_random)
  cat(sprintf("  inside: %d | outside: %d | undetermined: %d | DI threshold %.3f\n",
              aoa$n_inside, aoa$n_outside, aoa$n_na, aoa$threshold))

  # To blank the map outside the area. AOA is NA where a predictor is missing;
  # count those as outside.
  surf$.pred_masked <- ifelse(aoa$aoa$AOA %in% TRUE, surf$.pred, NA_real_)
  save_plot("08_surface.png",
            ggplot() +
              geom_sf(data = surf, aes(colour = .pred_masked), size = 0.6) +
              geom_sf(data = nc_boundary, fill = NA, colour = "grey20") +
              scale_colour_viridis_c(name = "Predicted y", na.value = "grey85") +
              theme_void() +
              ggtitle("Predicted surface, blank outside the area of applicability"))

  # predict_surface(covariates = points_sf) copies predictor values from the
  # nearest observation, so every grid point's dissimilarity index is 0 and
  # none is outside. Covariates from elsewhere (a raster, another survey, a
  # future scenario) can fall outside the training range, and the index marks
  # them:
  extreme <- st_sf(
    elevation   = c(mean(points_sf$elevation), max(points_sf$elevation) * 4),
    pop_density = c(mean(points_sf$pop_density), max(points_sf$pop_density) * 4),
    geometry    = st_geometry(points_sf)[1:2]
  )
  print(st_drop_geometry(area_of_applicability(extreme, model = rf_fit,
                                               folds = folds_random)$aoa)[, c("DI", "AOA")])
}

# =============================================================================
# Done
# =============================================================================

cat("\n==============================================================\n")
cat("  Figures written to:\n    ", out_dir, "\n", sep = "")
print(sort(basename(list.files(out_dir, pattern = "\\.png$", full.names = TRUE))))
cat("==============================================================\n")
