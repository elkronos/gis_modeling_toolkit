# =============================================================================
# Development helper: render the spatialkit vignette to HTML
#
# Run from the package root directory (the one containing DESCRIPTION).
# This is a dev-only convenience script; it is .Rbuildignore'd and not
# shipped with the package.  During R CMD build/check the vignette is
# built automatically via VignetteBuilder: knitr.
# =============================================================================

if (!file.exists("DESCRIPTION"))
  stop("Run this script from the package root directory.")

# -- 1. Install spatialkit from the local source tree (one-time) --------------
# devtools::install(".")

# -- 2. Install any missing dependencies (one-time) ---------------------------
pkgs <- c("sf", "dplyr", "ggplot2", "logger", "digest",
          "sp", "GWmodel",      # GWR models
          "geometry",           # Delaunay triangulation
          "knitr", "rmarkdown") # vignette rendering

to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

# Optional but nice for the 3-panel comparison at the end:
# install.packages("patchwork")

# -- 3. Render the vignette to HTML (output goes to docs/) --------------------
rmarkdown::render(
  input       = file.path("vignettes", "spatialkit_nc_demo.Rmd"),
  output_dir  = "docs",
  output_file = "spatialkit_nc_demo.html"
)
