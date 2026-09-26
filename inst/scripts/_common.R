# =============================================================================
# Shared setup for the numbered spatialkit tour scripts (01-, 02-, ...)
# =============================================================================
#
# Sourced by each of the numbered scripts. Nothing here is exported or supported
# API; it exists so the ten of them do not each repeat the simulated fixture and
# the plotting plumbing.
#
# WHERE FIGURES GO
#   By default every figure is drawn to the active graphics device, so that
#   sourcing a script in RStudio or R.app puts the plots in front of you.
#   To write them to files instead, set an output directory before sourcing:
#     Sys.setenv(SPATIALKIT_TOUR_OUTPUT = "~/spatialkit-tour")
#
# PACING
#   Between figures the scripts pause when the session is interactive. Set
#     Sys.setenv(SPATIALKIT_TOUR_PAUSE = "no")
#   to run straight through.
# =============================================================================

suppressPackageStartupMessages({
  library(spatialkit)
  library(sf)
})

have <- function(pkg) requireNamespace(pkg, quietly = TRUE)

tour_out <- local({
  d <- Sys.getenv("SPATIALKIT_TOUR_OUTPUT", unset = "")
  if (nzchar(d)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
})

tour_pause <- function() {
  if (identical(Sys.getenv("SPATIALKIT_TOUR_PAUSE"), "no")) return(invisible())
  if (!interactive()) return(invisible())
  invisible(readline("  [enter] "))
}

# Announce what the next figure is for. A plot with no stated expectation is
# not a check; it is a picture you will nod at.
look_for <- function(...) cat("\n  LOOK FOR:", paste0(..., collapse = ""), "\n")

# Draw to the device, or write a PNG when an output directory was set.
show_plot <- function(p, file, width = 9, height = 6) {
  if (nzchar(tour_out) && have("ggplot2")) {
    path <- file.path(tour_out, file)
    ggplot2::ggsave(path, p, width = width, height = height, dpi = 200)
    cat("  wrote", path, "\n")
  } else {
    print(p)
  }
  tour_pause()
  invisible(p)
}

step <- function(n, title) cat(sprintf("\n=== %s  %s\n", n, title))

skip_without <- function(pkgs, what) {
  missing <- pkgs[!vapply(pkgs, have, logical(1))]
  if (length(missing)) {
    cat(sprintf("  SKIPPED (%s): install.packages(c(%s))\n", what,
                paste(sprintf('"%s"', missing), collapse = ", ")))
    return(TRUE)
  }
  FALSE
}

# --- The fixture the whole tour runs on --------------------------------------
# A simulated exponential field on a 1000-unit square. The structure is known,
# so every diagnostic has something true to find: the covariance is exp(-d/80),
# which puts the true effective range at 240 units.
#
# Do not expect a fitted variogram to recover that number. Scripts 02 and 04
# report about 357 on `z`; the field on its own, with no nugget at all, still
# fits at about 315. Most of the gap is what fitting a single realisation over
# a limited span of lags does, and the nugget that `elev` and the measurement
# noise contribute widens it further. A fitted range is a working number here,
# not a measurement of the 240 -- which is worth knowing before anyone "fixes"
# the fixture to make the two agree.
#
# The columns: `z` is the response, built from `elev` plus the field. `noise` is
# a decoy. `slope` is a west-to-east gradient that does NOT enter `z` at all --
# but this particular field realisation happens to drift west to east too, so
# slope predicts z anyway. Script 08 measures that accident and shows what a
# variable search does with it.
tour_points <- function(n = 400, a = 80, seed = 42) {
  set.seed(seed)
  xy <- data.frame(x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$elev  <- stats::rnorm(n)                    # a real predictor
  xy$noise <- stats::rnorm(n)                    # a decoy
  xy$field <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% stats::rnorm(n))
  xy$z <- 0.8 * xy$elev + xy$field + stats::rnorm(n, sd = 0.3)
  # The west-to-east gradient. Holding out the eastern half of the map also
  # holds out the top of this variable's range, which is what makes the area of
  # applicability in script 07 have something to find. Scripts 01 to 06 ignore
  # the column.
  xy$slope <- 2 * xy$x / 1000 + stats::rnorm(n, sd = 0.2)
  # Move the square into UTM zone 32N, the CRS it is stamped with: x from
  # 500000 and y from 5000000. At 0 to 1000 it would sit on the equator west
  # of the zone. Every result is computed in planar units, so the shift only
  # changes the coordinates the scripts print.
  xy$x <- xy$x + 5e5
  xy$y <- xy$y + 5e6
  sf::st_as_sf(xy[, c("x", "y", "elev", "noise", "slope", "z")],
               coords = c("x", "y"), crs = 32632)
}

tour_boundary <- function(pts) clip_target_for(pts, expand = 0.02, quiet = TRUE)

# A learner with no optional dependencies: a cubic trend surface in the
# coordinates. It is deliberately capable of memorising location, which is what
# makes the fold scheme visible in script 03.
tour_fit <- function(train_sf, ...) {
  d <- sf::st_drop_geometry(train_sf)
  d$x <- sf::st_coordinates(train_sf)[, 1]
  d$y <- sf::st_coordinates(train_sf)[, 2]
  new_spatial_fit(
    subclass       = "tour_fit",
    engine         = stats::lm(z ~ elev + poly(x, 3) * poly(y, 3), data = d),
    formula        = z ~ elev + poly(x, 3) * poly(y, 3),
    response_var   = "z",
    predictor_vars = "elev",
    data_sf        = train_sf
  )
}
predict.tour_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(stats::fitted(object$engine))
  d <- sf::st_drop_geometry(newdata)
  d$x <- sf::st_coordinates(newdata)[, 1]
  d$y <- sf::st_coordinates(newdata)[, 2]
  as.numeric(stats::predict(object$engine, newdata = d))
}
residuals.tour_fit <- function(object, ...) stats::residuals(object$engine)
fitted.tour_fit    <- function(object, ...) stats::fitted(object$engine)
for (.g in c("predict", "residuals", "fitted"))
  registerS3method(.g, "tour_fit", get(paste0(.g, ".tour_fit")))
rm(.g)

if (!nzchar(tour_out))
  cat("Figures go to the active device.",
      "Set SPATIALKIT_TOUR_OUTPUT to write PNGs instead.\n")
