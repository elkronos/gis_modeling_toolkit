suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2)
})

# -- Load toolkit --------------------------------------------------------------
script_candidates <- c("R/gis_modeling_toolkit.R", "../R/gis_modeling_toolkit.R", "gis_modeling_toolkit.R")
script_path <- script_candidates[file.exists(script_candidates)][1]
if (is.na(script_path)) stop("Cannot find gis_modeling_toolkit.R in R/, ../R/, or current dir")
source(script_path)

# -- Data: North Carolina as a demo basemap -----------------------------------
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_boundary <- sf::st_union(nc)                               # single polygon boundary

set.seed(42)
pts <- sf::st_sample(nc_boundary, size = 300, type = "random", exact = TRUE) |> sf::st_sf()

# Align CRSs to boundary for sensible distances/areas
pts         <- ensure_projected(pts, nc_boundary)
nc_boundary <- ensure_projected(nc_boundary)
nc          <- sf::st_transform(nc, sf::st_crs(nc_boundary))

# -- Build tessellations -------------------------------------------------------
levels <- 8  # change if you want a different number of cells

bt_v <- build_tessellation(pts, levels = levels, method = "voronoi",
                           boundary = nc_boundary, seeds = "kmeans")
bt_h <- build_tessellation(pts, levels = levels, method = "hex",
                           boundary = nc_boundary, seeds = "kmeans")
bt_s <- build_tessellation(pts, levels = levels, method = "square",
                           boundary = nc_boundary, seeds = "kmeans")

# -- Plots (labels show polygon IDs) ------------------------------------------
p_v <- plot_tessellation_map(
  bt_v[[as.character(levels)]]$polygons,
  bt_v[[as.character(levels)]]$data,
  title   = sprintf("Voronoi (k = %d)", levels),
  boundary = nc_boundary,
  basemap  = nc,
  label    = "poly_id"
)

p_h <- plot_tessellation_map(
  bt_h[[as.character(levels)]]$polygons,
  bt_h[[as.character(levels)]]$data,
  title   = sprintf("Hex Grid (≈%d cells)", levels),
  boundary = nc_boundary,
  basemap  = nc,
  label    = "poly_id"
)

p_s <- plot_tessellation_map(
  bt_s[[as.character(levels)]]$polygons,
  bt_s[[as.character(levels)]]$data,
  title   = sprintf("Square Grid (≈%d cells)", levels),
  boundary = nc_boundary,
  basemap  = nc,
  label    = "poly_id"
)

# -- Combine to a single white-background image --------------------------------
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

common_theme <- ggplot2::theme(
  plot.background  = ggplot2::element_rect(fill = "white", color = NA),
  panel.background = ggplot2::element_rect(fill = "white", color = NA),
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank()
)

combined <- (p_v | p_h | p_s) + patchwork::plot_layout(guides = "collect") & common_theme

ggplot2::ggsave("nc_tessellations_triptych.png", combined,
                width = 18, height = 6, dpi = 240, bg = "white")
cat("Saved: nc_tessellations_triptych.png\n")
