#!/usr/bin/env Rscript
# =============================================================================
# spatialkit Demo — North Carolina with Synthetic Data
# =============================================================================
#
# A self-contained example showing how spatialkit handles different
# tessellations (Voronoi, hex grid, square grid, Delaunay triangles)
# and model types (GWR) on synthetic spatial data.
#
# The boundary comes from sf's built-in nc.shp — no external files needed.
#
# SETUP:
#   1. Install the spatialkit package:
#        devtools::install("path/to/spatialkit")
#   2. Install required packages:
#        install.packages(c("sf", "dplyr", "ggplot2", "logger", "digest"))
#        install.packages(c("sp", "GWmodel"))     # for GWR models
#        install.packages("geometry")              # for Delaunay triangulation
#   3. Source or run this script:
#        source("example_nc_demo.R")
#
# OUTPUT:
#   Saves PNG images to ./output/ showing each tessellation + model combo.
# =============================================================================

library(spatialkit)
library(sf)
library(dplyr)
library(ggplot2)

set.seed(42)
dir.create("output", showWarnings = FALSE)

# =============================================================================
# 1. LOAD NORTH CAROLINA BOUNDARY
# =============================================================================
# The sf package ships with nc.shp (100 county polygons). We dissolve them
# into a single state outline and project to NAD83 / NC State Plane (ftUS).

nc_counties <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_boundary <- nc_counties |>
  st_union() |>
  st_transform(2264) |>
  st_as_sf()

cat("✓ Loaded North Carolina boundary (EPSG:2264)\n")

# =============================================================================
# 2. GENERATE FAKE SPATIAL DATA POINTS
# =============================================================================

n_points <- 300

# Sample random points inside the state boundary
pts_raw <- st_sample(nc_boundary, size = n_points, type = "random")
pts_coords <- st_coordinates(pts_raw)
x_coords <- pts_coords[, 1]
y_coords <- pts_coords[, 2]

# Predictor: elevation (gradient west→east + noise)
elevation <- scale(x_coords)[, 1] * 500 + rnorm(n_points, 3000, 400)

# Predictor: population density (higher near Charlotte & Raleigh)
city1 <- c(1530000, 550000)   # Charlotte-ish in EPSG:2264
city2 <- c(2150000, 750000)   # Raleigh-ish in EPSG:2264
dist_to_city <- pmin(
  sqrt((x_coords - city1[1])^2 + (y_coords - city1[2])^2),
  sqrt((x_coords - city2[1])^2 + (y_coords - city2[2])^2)
)
pop_density <- exp(-dist_to_city / 400000) * 5000 + rnorm(n_points, 200, 100)
pop_density <- pmax(pop_density, 10)

# Response: spatially varying function of predictors
y_response <- 50 +
  0.01  * elevation +
  0.005 * pop_density +
  2.0   * sin(x_coords / 300000) * cos(y_coords / 300000) +
  rnorm(n_points, 0, 5)

points_sf <- st_sf(
  y           = y_response,
  elevation   = elevation,
  pop_density = pop_density,
  geometry    = pts_raw
)

cat(sprintf("✓ Generated %d fake observation points\n", n_points))

# =============================================================================
# 3. BUILD FOUR TESSELLATION TYPES
# =============================================================================

cat("\n--- Building tessellations ---\n")

# 3a. Voronoi (k-means seeds)
seeds <- get_voronoi_seeds(
  boundary      = nc_boundary,
  sample_points = points_sf,
  method        = "kmeans",
  n             = 40
)
tess_voronoi <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "voronoi", clip = TRUE, quiet = TRUE
)
cat(sprintf("  Voronoi: %d cells\n", nrow(tess_voronoi$cells)))

# 3b. Hexagonal grid (~50 cells)
tess_hex <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "hex", approx_n_cells = 50, clip = TRUE, quiet = TRUE
)
cat(sprintf("  Hex grid: %d cells\n", nrow(tess_hex$cells)))

# 3c. Square grid (~50 cells)
tess_square <- build_tessellation(
  points_sf, boundary = nc_boundary,
  method = "square", approx_n_cells = 50, clip = TRUE, quiet = TRUE
)
cat(sprintf("  Square grid: %d cells\n", nrow(tess_square$cells)))

# 3d. Delaunay triangles
tess_tri <- tryCatch(
  build_tessellation(
    points_sf, boundary = nc_boundary,
    method = "triangles", clip = TRUE, quiet = TRUE
  ),
  error = function(e) {
    cat("  Triangles: skipped —", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(tess_tri))
  cat(sprintf("  Delaunay triangles: %d cells\n", nrow(tess_tri$cells)))

# =============================================================================
# 4. CHOROPLETH HELPER
# =============================================================================

make_choropleth <- function(tess, boundary, points, fill_var = "y",
                            palette = "viridis", title = NULL,
                            legend_title = "Mean Response (y)",
                            subtitle_extra = NULL) {
  cells <- tess$cells
  id_col <- if ("cell_id" %in% names(cells)) "cell_id" else "poly_id"
  if (!id_col %in% names(cells)) {
    cells$cell_id <- seq_len(nrow(cells))
    id_col <- "cell_id"
  }

  assigned <- assign_features_to_polygons(points, cells, polygon_id_col = id_col)

  cell_summary <- assigned |>
    st_drop_geometry() |>
    group_by(.data[[id_col]]) |>
    summarise(
      fill_value = mean(.data[[fill_var]], na.rm = TRUE),
      n_obs      = n(),
      .groups    = "drop"
    )

  cells <- left_join(cells, cell_summary, by = id_col)

  sub <- subtitle_extra %||%
    sprintf("%d cells  |  %d observations", nrow(cells), nrow(points))

  plot_tessellation_map(
    tessellation_sf = cells,
    boundary        = boundary,
    fill_col        = "fill_value",
    palette         = palette,
    tile_alpha      = 0.9,
    outline_col     = "white",
    outline_size    = 0.3,
    boundary_col    = "grey20",
    boundary_size   = 0.8,
    legend_title    = legend_title,
    title           = title,
    subtitle        = sub,
    caption         = "spatialkit demo — synthetic North Carolina data"
  )
}

# =============================================================================
# 5. PLOT CHOROPLETHS
# =============================================================================

cat("\n--- Generating choropleth maps ---\n")

p_vor <- make_choropleth(tess_voronoi, nc_boundary, points_sf,
                          title = "Voronoi Tessellation — Mean Response")
p_hex <- make_choropleth(tess_hex, nc_boundary, points_sf,
                          title = "Hexagonal Grid — Mean Response")
p_sq  <- make_choropleth(tess_square, nc_boundary, points_sf,
                          title = "Square Grid — Mean Response")

ggsave("output/01_voronoi_choropleth.png",  p_vor, width = 10, height = 6, dpi = 200)
ggsave("output/02_hex_choropleth.png",      p_hex, width = 10, height = 6, dpi = 200)
ggsave("output/03_square_choropleth.png",   p_sq,  width = 10, height = 6, dpi = 200)

if (!is.null(tess_tri)) {
  p_tri <- make_choropleth(tess_tri, nc_boundary, points_sf,
                            title = "Delaunay Triangulation — Mean Response")
  ggsave("output/04_delaunay_choropleth.png", p_tri, width = 10, height = 6, dpi = 200)
}

cat("  Saved choropleth plots to output/\n")

# =============================================================================
# 6. FIT GWR MODEL
# =============================================================================

response_var   <- "y"
predictor_vars <- c("elevation", "pop_density")

cat("\n--- Fitting GWR model ---\n")
gwr_fit <- tryCatch({
  fit <- fit_gwr_model(
    data_sf        = points_sf,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    adaptive       = TRUE,
    kernel         = "bisquare"
  )
  cat("  ✓ GWR model fitted successfully\n")
  cat(sprintf("    Bandwidth: %.1f | R²: %.3f\n",
              fit$info$bandwidth, fit$metrics$r_squared))
  fit
}, error = function(e) {
  cat("  ✗ GWR skipped:", conditionMessage(e), "\n")
  NULL
})

# =============================================================================
# 7. MAP GWR RESIDUALS
# =============================================================================

if (!is.null(gwr_fit)) {

  pts_with_gwr <- points_sf
  pts_with_gwr$gwr_fitted   <- as.numeric(fitted(gwr_fit))
  pts_with_gwr$gwr_residual <- as.numeric(residuals(gwr_fit))
  pts_with_gwr$abs_error    <- abs(pts_with_gwr$gwr_residual)

  for (tess_info in list(
    list(tess = tess_voronoi, name = "voronoi", label = "Voronoi"),
    list(tess = tess_hex,     name = "hex",     label = "Hex Grid"),
    list(tess = tess_square,  name = "square",  label = "Square Grid")
  )) {
    cells <- tess_info$tess$cells
    id_col <- if ("cell_id" %in% names(cells)) "cell_id" else "poly_id"
    if (!id_col %in% names(cells)) {
      cells$cell_id <- seq_len(nrow(cells))
      id_col <- "cell_id"
    }

    asgn <- assign_features_to_polygons(pts_with_gwr, cells, polygon_id_col = id_col)

    cell_err <- asgn |>
      st_drop_geometry() |>
      group_by(.data[[id_col]]) |>
      summarise(mean_abs_error = mean(abs_error, na.rm = TRUE), .groups = "drop")

    cells <- left_join(cells, cell_err, by = id_col)

    p <- plot_tessellation_map(
      tessellation_sf = cells,
      boundary        = nc_boundary,
      fill_col        = "mean_abs_error",
      palette         = "magma",
      tile_alpha      = 0.9,
      outline_col     = "white",
      outline_size    = 0.3,
      boundary_col    = "grey20",
      boundary_size   = 0.8,
      legend_title    = "Mean |Residual|",
      title           = sprintf("GWR Residuals — %s", tess_info$label),
      subtitle        = sprintf("Abs. residual aggregated to %d cells", nrow(cells)),
      caption         = "spatialkit demo — GWR on synthetic NC data"
    )
    fname <- sprintf("output/05_gwr_error_%s.png", tess_info$name)
    ggsave(fname, p, width = 10, height = 6, dpi = 200)
    cat(sprintf("  Saved %s\n", fname))
  }
}

# =============================================================================
# 8. SIDE-BY-SIDE COMPARISON PANEL
# =============================================================================

cat("\n--- Building comparison panel ---\n")
tryCatch({
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)

    p_panel <- (p_vor | p_hex | p_sq) +
      plot_annotation(
        title    = "spatialkit — Tessellation Comparison (North Carolina)",
        subtitle = sprintf("%d observations  |  3 tessellation methods", n_points),
        theme    = theme(
          plot.title    = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 11, color = "grey40")
        )
      )
    ggsave("output/06_comparison_panel.png", p_panel, width = 18, height = 6, dpi = 200)
    cat("  Saved comparison panel to output/06_comparison_panel.png\n")
  } else {
    cat("  Install 'patchwork' for the side-by-side panel: install.packages('patchwork')\n")
  }
}, error = function(e) {
  cat("  Panel generation skipped:", conditionMessage(e), "\n")
})

# =============================================================================
# 9. CROSS-VALIDATION
# =============================================================================

cat("\n--- Cross-validation (GWR, 5-fold) ---\n")
tryCatch({
  cv_results <- cv_gwr(
    data_sf        = points_sf,
    response_var   = response_var,
    predictor_vars = predictor_vars,
    k              = 5,
    adaptive       = TRUE
  )
  cat(sprintf("  CV RMSE: %.3f  |  CV R²: %.3f\n",
              cv_results$summary$rmse, cv_results$summary$r_squared))
}, error = function(e) {
  cat("  CV skipped:", conditionMessage(e), "\n")
})

# =============================================================================
# DONE
# =============================================================================

cat("\n══════════════════════════════════════════════════════════\n")
cat("  All outputs saved to ./output/\n")
cat("  Files produced:\n")
cat("    01_voronoi_choropleth.png  — Voronoi filled by mean(y)\n")
cat("    02_hex_choropleth.png      — Hex grid filled by mean(y)\n")
cat("    03_square_choropleth.png   — Square grid filled by mean(y)\n")
cat("    04_delaunay_choropleth.png — Delaunay filled by mean(y)\n")
cat("    05_gwr_error_*.png         — GWR residuals per tessellation\n")
cat("    06_comparison_panel.png    — Side-by-side 3-panel view\n")
cat("══════════════════════════════════════════════════════════\n")
