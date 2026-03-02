#' Plot a tessellation map with optional boundary, seeds, and features
#'
#' Builds a layered ggplot2 map of polygon tessellations and optional
#' overlays for a study boundary, seed points, and additional features.
#'
#' @param tessellation_sf An sf POLYGON/MULTIPOLYGON layer. Required.
#' @param boundary Optional sf/sfc polygon outline layer.
#' @param seeds_sf Optional sf/sfc point layer of seed locations.
#' @param features_sf Optional sf/sfc layer of additional features.
#' @param fill_col Column name in tessellation_sf to map to fill. NULL = no fill.
#' @param palette Viridis palette name. Default "viridis".
#' @param na_fill Fill for NA values. Default "grey90".
#' @param tile_alpha Alpha for filled polygons. Default 0.9.
#' @param outline_col,outline_size Tessellation outline aesthetics.
#' @param features_col,features_size Feature overlay aesthetics.
#' @param seeds_col,seeds_size Seed point aesthetics.
#' @param boundary_col,boundary_size Boundary outline aesthetics.
#' @param labels Logical; draw per-cell labels. Default FALSE.
#' @param label_col Column for label text. Default "grid_id".
#' @param label_size Label text size. Default 2.7.
#' @param legend Logical; show fill legend. Default TRUE.
#' @param legend_title Optional legend title.
#' @param theme A ggplot2 theme. Default theme_void().
#' @param target_crs Optional CRS for plotting.
#' @param title,subtitle,caption Plot annotations.
#' @param xlim,ylim Optional numeric vectors of length 2 for coordinate limits
#'   (in the plot CRS). Default NULL (auto).
#' @param expand Logical; expand plot area slightly beyond data limits.
#'   Default TRUE.
#' @return A ggplot2 object.
#' @export
plot_tessellation_map <- function(tessellation_sf,
                                  boundary = NULL,
                                  seeds_sf = NULL,
                                  features_sf = NULL,
                                  fill_col = NULL,
                                  palette = "viridis",
                                  na_fill = "grey90",
                                  tile_alpha = 0.9,
                                  outline_col = "white",
                                  outline_size = 0.2,
                                  features_col = "#333333",
                                  features_size = 0.5,
                                  seeds_col = "#1f77b4",
                                  seeds_size = 1.5,
                                  boundary_col = "#111111",
                                  boundary_size = 0.6,
                                  labels = FALSE,
                                  label_col = "grid_id",
                                  label_size = 2.7,
                                  legend = TRUE,
                                  legend_title = NULL,
                                  theme = ggplot2::theme_void(),
                                  target_crs = NULL,
                                  title = NULL,
                                  subtitle = NULL,
                                  caption = NULL,
                                  xlim = NULL,
                                  ylim = NULL,
                                  expand = TRUE) {

  # --- validation ---
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_tessellation_map() requires package 'ggplot2'. Install it with install.packages('ggplot2').", call. = FALSE)
  if (missing(tessellation_sf) || is.null(tessellation_sf))
    stop("plot_tessellation_map(): 'tessellation_sf' is required.")
  .assert_sf(tessellation_sf, c("POLYGON", "MULTIPOLYGON"), "tessellation_sf")

  # --- pick plot CRS ---
  plot_crs <- if (!is.null(target_crs)) sf::st_crs(target_crs) else sf::st_crs(tessellation_sf)
  if (is.na(plot_crs))
    .log_warn("plot_tessellation_map(): CRS is undefined; drawing layers as-is.")

  # --- coerce & transform helper ---
  maybe_transform <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "sfc")) x <- sf::st_as_sf(x)
    if (!inherits(x, c("sf", "sfc")))
      stop("plot_tessellation_map(): all layers must be sf/sfc when provided.")
    if (is.na(plot_crs)) return(x)
    x_crs <- sf::st_crs(x)
    if (is.na(x_crs)) return(x)
    if (x_crs == plot_crs) return(x)
    sf::st_transform(x, plot_crs)
  }

  tess <- maybe_transform(tessellation_sf)
  bnd  <- maybe_transform(boundary)
  sds  <- maybe_transform(seeds_sf)
  fea  <- maybe_transform(features_sf)

  # --- build ggplot ---
  p <- ggplot2::ggplot()

  # Fill mapping
  has_fill <- !is.null(fill_col) && fill_col %in% names(tess)
  if (has_fill) {
    tess$`..__fill__` <- tess[[fill_col]]
    p <- p + ggplot2::geom_sf(
      data = tess,
      mapping = ggplot2::aes(fill = `..__fill__`),
      color = outline_col,
      linewidth = outline_size,
      alpha = tile_alpha
    )
  } else {
    p <- p + ggplot2::geom_sf(
      data = tess, color = outline_col, linewidth = outline_size, fill = NA
    )
  }

  # Features below boundary
  if (!is.null(fea) && nrow(fea) > 0L && !all(sf::st_is_empty(fea)))
    p <- p + ggplot2::geom_sf(data = fea, color = features_col,
                               linewidth = features_size, size = features_size)

  # Seeds
  if (!is.null(sds) && nrow(sds) > 0L && !all(sf::st_is_empty(sds)))
    p <- p + ggplot2::geom_sf(data = sds, color = seeds_col, size = seeds_size)

  # Boundary outline
  if (!is.null(bnd) && nrow(bnd) > 0L && !all(sf::st_is_empty(bnd)))
    p <- p + ggplot2::geom_sf(data = bnd, fill = NA, color = boundary_col,
                               linewidth = boundary_size)

  # Labels
  if (isTRUE(labels)) {
    if (!label_col %in% names(tess)) {
      .log_warn("plot_tessellation_map(): label_col '%s' not found; skipping labels.",
                label_col)
    } else if (nrow(tess) > 0L && !all(sf::st_is_empty(tess))) {
      centers <- suppressWarnings(sf::st_point_on_surface(tess))
      centers$`..__lab__` <- tess[[label_col]]
      p <- p + ggplot2::geom_sf_text(
        data = centers, ggplot2::aes(label = `..__lab__`), size = label_size
      )
    }
  }

  # Fill scale
  if (has_fill) {
    lab <- legend_title %||% fill_col
    is_cont <- is.numeric(tess$`..__fill__`)
    if (is_cont) {
      p <- p + ggplot2::scale_fill_viridis_c(
        option = palette, na.value = na_fill, name = lab,
        guide = if (legend) "colourbar" else "none"
      )
    } else {
      p <- p + ggplot2::scale_fill_viridis_d(
        option = palette, na.value = na_fill, name = lab,
        guide = if (legend) "legend" else "none", drop = FALSE
      )
    }
  }

  # Coord, theme, labs
  coord_args <- list(expand = expand)
  if (!is.na(plot_crs)) coord_args$crs <- plot_crs
  if (!is.null(xlim))   coord_args$xlim <- xlim
  if (!is.null(ylim))   coord_args$ylim <- ylim

  p <- p +
    do.call(ggplot2::coord_sf, coord_args) +
    theme +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption)

  p
}
