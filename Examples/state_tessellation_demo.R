suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -- Load toolkit --------------------------------------------------------------
script_candidates <- c("R/gis_modeling_toolkit.R", "../R/gis_modeling_toolkit.R", "gis_modeling_toolkit.R")
script_path <- script_candidates[file.exists(script_candidates)][1]
if (is.na(script_path)) stop("Cannot find gis_modeling_toolkit.R in R/, ../R/, or current dir")
suppressWarnings(suppressMessages(source(script_path)))

# -- Data: North Carolina as a demo basemap -----------------------------------
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_boundary <- sf::st_as_sf(sf::st_union(nc))

# ✅ Project first (for proper distances/areas/label points)
nc_boundary <- ensure_projected(nc_boundary)
nc          <- sf::st_transform(nc, sf::st_crs(nc_boundary))

# Sample points inside the projected boundary
set.seed(42)
pts <- sf::st_sample(nc_boundary, size = 300, type = "random", exact = TRUE) |>
  sf::st_sf()

# -- Synthesize a response and fit a simple model ------------------------------
xy <- sf::st_coordinates(pts)
pts$xs <- as.numeric(scale(xy[,1]))
pts$ys <- as.numeric(scale(xy[,2]))
set.seed(42)
pts$resp <- 2 + 0.7*pts$xs - 0.4*pts$ys + 0.2*(pts$xs*pts$ys) + rnorm(nrow(pts), sd = 0.5)

glm_fit      <- lm(resp ~ xs + ys + I(xs*ys), data = sf::st_drop_geometry(pts))
pts$pred_glm <- as.numeric(predict(glm_fit, newdata = sf::st_drop_geometry(pts)))

# -- Build tessellations (UPDATED API) ----------------------------------------
levels <- 8L

# Voronoi seeds via k-means, then build Voronoi cells
seeds_v <- get_voronoi_seeds(boundary = nc_boundary, method = "kmeans", n = levels, set_seed = 42)
bt_v <- build_tessellation(
  points_sf = seeds_v,
  method    = "voronoi",
  boundary  = nc_boundary,
  clip      = TRUE
)

# Hex and Square grids (approximate number of cells)
bt_h <- build_tessellation(points_sf = pts, method = "hex",
                           boundary = nc_boundary, approx_n_cells = levels)
bt_s <- build_tessellation(points_sf = pts, method = "square",
                           boundary = nc_boundary, approx_n_cells = levels)

cells_v <- bt_v$cells
cells_h <- bt_h$cells
cells_s <- bt_s$cells

# -- Assign points to polygons (NEW helper uses polygon_id_col) ----------------
to_poly <- function(polys, pts) {
  # Normalize Voronoi 'cell_id' -> 'poly_id' for joins/labels
  if (!"poly_id" %in% names(polys) && "cell_id" %in% names(polys)) {
    polys$poly_id <- polys$cell_id
  }
  assign_features_to_polygons(
    features_sf    = pts,
    polygons_sf    = polys,
    polygon_id_col = "poly_id",
    predicate      = sf::st_within
  )
}

ap_v <- to_poly(cells_v, pts)
ap_h <- to_poly(cells_h, pts)
ap_s <- to_poly(cells_s, pts)

# -- Mapping helper (UPDATED joins/labels) ------------------------------------
make_outcome_map <- function(polys, assigned_pts, title, boundary, basemap) {
  # Ensure polygons have poly_id (Voronoi may only have cell_id)
  if (!"poly_id" %in% names(polys) && "cell_id" %in% names(polys)) {
    polys$poly_id <- polys$cell_id
  }
  # Summarize per cell (mean of predicted outcome)
  sum_tbl <- summarize_by_cell(
    assigned_points_sf = assigned_pts,
    response_var       = "pred_glm"
  )
  polys_f <- dplyr::left_join(polys, sum_tbl, by = "poly_id")
  
  labs <- suppressWarnings(sf::st_point_on_surface(polys_f))
  labs <- labs[!is.na(polys_f$mean_response), , drop = FALSE]
  
  ggplot() +
    geom_sf(data = polys_f, aes(fill = mean_response), color = "white", linewidth = 0.2) +
    geom_sf(data = boundary, fill = NA, color = "black",  linewidth = 0.3) +
    geom_sf(data = basemap,  fill = NA, color = "grey80", linewidth = 0.2) +
    geom_sf_text(data = labs, aes(label = poly_id), size = 3, check_overlap = TRUE) +
    # keep a scale here; we'll override with one shared scale below
    scale_fill_viridis_c(name = "Mean predicted\noutcome", option = "C", direction = -1) +
    coord_sf(crs = sf::st_crs(boundary)) +
    labs(title = title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Build maps for each tessellation
p_v <- make_outcome_map(
  polys        = cells_v,
  assigned_pts = ap_v,
  title        = sprintf("Voronoi (k = %d): mean predicted outcome", levels),
  boundary     = nc_boundary,
  basemap      = nc
)

p_h <- make_outcome_map(
  polys        = cells_h,
  assigned_pts = ap_h,
  title        = sprintf("Hex Grid (≈%d cells): mean predicted outcome", levels),
  boundary     = nc_boundary,
  basemap      = nc
)

p_s <- make_outcome_map(
  polys        = cells_s,
  assigned_pts = ap_s,
  title        = sprintf("Square Grid (≈%d cells): mean predicted outcome", levels),
  boundary     = nc_boundary,
  basemap      = nc
)

# ---- One legend only: apply a SINGLE shared fill scale to all panels ---------
# Use a global range so Patchwork can collect to exactly one legend.
global_limits <- range(pts$pred_glm, na.rm = TRUE)

common_theme <- theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
)

combined <- (p_v | p_h | p_s) +
  plot_layout(guides = "collect") &
  common_theme &
  scale_fill_viridis_c(
    name = "Mean predicted\noutcome",
    option = "C",
    direction = -1,
    limits = global_limits
  ) &
  theme(legend.position = "bottom")

ggplot2::ggsave("nc_tessellations_triptych.png", combined,
                width = 18, height = 6, dpi = 240, bg = "white")

cat("Saved: nc_tessellations_triptych.png\n")
