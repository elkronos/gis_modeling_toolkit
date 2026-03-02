# =============================================================================
# Run this in RStudio to render the spatialkit vignette
# =============================================================================
# Place spatialkit_nc_demo.Rmd into:
#   C:/repos/gis_modeling_toolkit/vignettes/
#
# Then run this script (or just the lines you need).
# =============================================================================

# -- 1. Install spatialkit from your local repo (one-time) --------------------
devtools::install("C:/repos/gis_modeling_toolkit")

# -- 2. Install any missing dependencies (one-time) --------------------------
pkgs <- c("sf", "dplyr", "ggplot2", "logger", "digest",
           "sp", "GWmodel",     # GWR models
           "geometry",           # Delaunay triangulation
           "knitr", "rmarkdown") # vignette rendering

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

# Optional but nice for the 3-panel comparison at the end:
# install.packages("patchwork")

# -- 3. Render the vignette to HTML ------------------------------------------
rmarkdown::render(
  input       = "C:/repos/gis_modeling_toolkit/vignettes/spatialkit_nc_demo.Rmd",
  output_dir  = "C:/repos/gis_modeling_toolkit/vignettes",
  output_file = "spatialkit_nc_demo.html"
)

# This opens the result in your browser / RStudio viewer automatically.
# All plots are embedded inline in the HTML — no separate PNG files needed.

# -- Alternative: build ALL vignettes via devtools ----------------------------
# devtools::build_vignettes("C:/repos/gis_modeling_toolkit")
