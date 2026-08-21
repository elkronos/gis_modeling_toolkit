# =============================================================================
# Development helper: render the spatialkit vignette to HTML
#
#   Rscript dev/render_vignette.R
#
# Run from the package root directory (the one containing DESCRIPTION).
# This is a dev-only convenience script; it is .Rbuildignore'd and not
# shipped with the package.  During R CMD build/check the vignette is
# built automatically via VignetteBuilder: knitr.
# =============================================================================

if (!file.exists("DESCRIPTION"))
  stop("Run this script from the package root directory.")

# -- 1. Load spatialkit FROM THE WORKING TREE --------------------------------
# The vignette's first chunk calls library(spatialkit).  library() attaches an
# already-loaded namespace rather than re-reading the library path, so calling
# load_all() here means the vignette renders against the source tree.  Without
# it the vignette silently documents the INSTALLED (released) package -- the
# same trap that once cost a five-hour baseline run.
if (!requireNamespace("pkgload", quietly = TRUE))
  stop("pkgload is required so the vignette renders against the working tree. ",
       "install.packages('pkgload')", call. = FALSE)
pkgload::load_all(".", quiet = TRUE)
cat("spatialkit loaded from: ",
    tryCatch(getNamespaceInfo(asNamespace("spatialkit"), "path"),
             error = function(e) NA_character_), "\n", sep = "")

# -- 2. Dependencies ---------------------------------------------------------
# 'sp' and 'GWmodel' are Suggests, and the vignette degrades gracefully when
# they are absent -- but it degrades by SKIPPING the GWR fit and the
# cross-validation, which are exactly the sections whose metric output this
# render exists to verify.  A render without them would show no metric lines
# and look identical to a broken accessor.  So they are required here even
# though they are optional for the package.
required <- c("sf", "dplyr", "ggplot2", "logger", "digest",
              "knitr", "rmarkdown",
              "sp", "GWmodel")
optional <- c("geometry",   # true Delaunay; otherwise falls back to the hull
              "patchwork")  # the 3-panel side-by-side comparison

.have <- function(p) vapply(p, requireNamespace, logical(1), quietly = TRUE)
miss_req <- required[!.have(required)]
miss_opt <- optional[!.have(optional)]

if (length(miss_req) || length(miss_opt)) {
  # Rscript starts with repos = c(CRAN = "@CRAN@"), so install.packages()
  # aborts with "trying to use CRAN without setting a mirror" before it
  # installs anything.  Set one only if the user has not already.
  r <- getOption("repos")
  if (is.null(r) || is.na(r["CRAN"]) || !nzchar(r["CRAN"]) ||
      identical(unname(r["CRAN"]), "@CRAN@"))
    options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (length(miss_req)) {
  message("Installing required packages: ", paste(miss_req, collapse = ", "))
  install.packages(miss_req)
  still <- miss_req[!.have(miss_req)]
  if (length(still))
    stop("Could not install: ", paste(still, collapse = ", "),
         ". Install these by hand, then re-run.", call. = FALSE)
}

if (length(miss_opt)) {
  message("Installing optional packages: ", paste(miss_opt, collapse = ", "))
  try(install.packages(miss_opt), silent = TRUE)
  still <- miss_opt[!.have(miss_opt)]
  if (length(still))
    message("NOTE: unavailable, those sections will degrade or skip: ",
            paste(still, collapse = ", "))
}

# -- 3. Render the vignette to HTML (output goes to docs/) --------------------
rmarkdown::render(
  input       = file.path("vignettes", "spatialkit_nc_demo.Rmd"),
  output_dir  = "docs",
  output_file = "spatialkit_nc_demo.html"
)

# -- 4. Verify the metric lines actually produced output ---------------------
# sprintf() returns character(0) if any argument has length zero, and
# cat(character(0)) prints nothing -- which is how the wrong accessors went
# unnoticed through a CRAN release.  Grep the render rather than trusting it.
html <- readLines(file.path("docs", "spatialkit_nc_demo.html"),
                  warn = FALSE, encoding = "UTF-8")
html <- paste(html, collapse = "\n")
ok <- TRUE
for (pat in c("Bandwidth: [0-9]", "CV RMSE: [0-9]")) {
  if (!grepl(pat, html)) {
    ok <- FALSE
    cat(sprintf("FAIL: no rendered output matching /%s/\n", pat))
  }
}
if (ok) {
  cat("\nPASS: both metric lines rendered with numbers.\n")
} else {
  cat("\nThe metric lines are still empty. Either an accessor is wrong or the",
      "\nGWR/CV chunks were skipped. Check the render before committing.\n")
}
