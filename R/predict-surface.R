# =============================================================================
# Prediction to a surface
# =============================================================================

#' Build a regular point grid covering a bounding box
#'
#' Cell centres, not corners, so every point lies inside the cell it represents.
#'
#' @param bb An \code{sf} bbox.
#' @param crs Target CRS for the result.
#' @param cell_size Edge length in CRS units, or NULL to derive from n_cells.
#' @param n_cells Approximate number of cells when \code{cell_size} is NULL.
#' @return An sf POINT layer with columns \code{..grid_x}, \code{..grid_y}.
#' @keywords internal
#' @noRd
.make_prediction_grid <- function(bb, crs, cell_size = NULL, n_cells = 10000L) {
  w <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0)
    stop("predict_surface(): the fitted data has no finite extent.", call. = FALSE)

  if (is.null(cell_size)) {
    n_cells <- max(1L, as.integer(n_cells))
    cell_size <- sqrt((w * h) / n_cells)
  }
  if (!is.finite(cell_size) || cell_size <= 0)
    stop("predict_surface(): `cell_size` must be a positive finite number.",
         call. = FALSE)

  # seq(from, to, by) ERRORS with "wrong sign in 'by' argument" when from > to,
  # rather than returning length 0 -- so a cell_size wider than the extent used
  # to abort with an internal-looking message, and the length-0 guards below
  # were unreachable.  Build the axis explicitly instead.
  .axis <- function(lo, hi) {
    lo <- as.numeric(lo); hi <- as.numeric(hi)
    n <- floor((hi - lo) / cell_size)
    if (!is.finite(n) || n < 1L) return(lo + (hi - lo) / 2)   # one centred cell
    lo + cell_size / 2 + seq.int(0L, n - 1L) * cell_size
  }
  xs <- .axis(bb[["xmin"]], bb[["xmax"]])
  ys <- .axis(bb[["ymin"]], bb[["ymax"]])

  g <- expand.grid(..grid_x = xs, ..grid_y = ys, KEEP.OUT.ATTRS = FALSE)
  out <- sf::st_as_sf(g, coords = c("..grid_x", "..grid_y"), crs = crs,
                      remove = FALSE)
  attr(out, "cell_size") <- cell_size
  out
}


#' Predict a fitted spatial model onto a regular grid
#'
#' Builds a prediction surface over the extent of the training data (or over a
#' grid you supply), predicts in chunks, and returns an \code{sf} layer.
#'
#' \code{predict()} on a \code{spatial_fit} requires \code{newdata} to be
#' constructed by hand, which makes the most common downstream task -- produce
#' a map -- more work than it should be.  This wraps the grid construction,
#' covariate join, chunking and CRS handling.
#'
#' Prediction over a grid is embarrassingly parallel in the sense that rows do
#' not interact, so it is chunked: for \code{bayesian_fit} the posterior draw
#' matrix is \code{n_draws x n_newdata}, which will exhaust memory on a fine
#' grid long before the fit itself would.
#'
#' @param object A \code{spatial_fit} (e.g. from \code{fit_gwr_model()} or
#'   \code{fit_bayesian_spatial_model()}).
#' @param grid Optional \code{sf} POINT layer to predict onto.  When
#'   \code{NULL}, a regular grid is built over the training extent.  Must have
#'   at least one row.
#' @param cell_size Grid resolution in CRS units.  Ignored when \code{grid} is
#'   supplied; when \code{NULL}, derived from \code{n_cells}.
#' @param n_cells Approximate cell count used to derive \code{cell_size}.
#'   Default 10000.  Also ignored when \code{grid} is supplied -- the grid you
#'   pass is used verbatim.
#' @param boundary Optional polygonal \code{sf}/\code{sfc}; grid points outside
#'   it are dropped.
#' @param covariates Optional \code{sf} layer carrying the model's predictors.
#'   Required when the model has predictors and \code{grid} does not already
#'   contain them.  Values are taken from the nearest feature.
#' @param chunk_size Rows per prediction call. Default 5000.
#' @param se Logical; also return a standard-error/posterior-SD column where the
#'   backend supports it.  Default FALSE.
#' @param ... Passed to \code{predict()}.
#' @return An \code{sf} POINT layer with a \code{.pred} column (and
#'   \code{.pred_se} when \code{se = TRUE} and available).  For an
#'   auto-generated grid the resolution is attached as attribute
#'   \code{"cell_size"}.  For a user-supplied \code{grid} it is only whatever
#'   \code{"cell_size"} attribute that object already carried -- usually
#'   \code{NULL}, and \code{NULL} for certain if the grid had to be
#'   re-projected, since \code{st_transform()} does not preserve custom
#'   attributes.  The resolution of a grid you built is not this function's to
#'   infer.
#' @family prediction
#' @examples
#' \donttest{
#' # Any spatial_fit works here; a forest keeps the example free of the
#' # optional GWR/Stan backends.
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   pts <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$price <- 10 + 0.01 * st_coordinates(pts)[, 1] + 2 * pts$elev + rnorm(n)
#'   fit  <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
#'   surf <- predict_surface(fit, n_cells = 500, covariates = pts)
#'   surf[".pred"]
#'   # Check where that surface is extrapolating before mapping it.
#'   area_of_applicability(surf, model = fit)
#' }
#' }
#' @export
predict_surface <- function(object, grid = NULL, cell_size = NULL,
                            n_cells = 10000L, boundary = NULL,
                            covariates = NULL, chunk_size = 5000L,
                            se = FALSE, ...) {
  if (!inherits(object, "spatial_fit"))
    stop("predict_surface(): `object` must be a spatial_fit.", call. = FALSE)

  train <- object$data_sf
  if (!inherits(train, "sf"))
    stop("predict_surface(): the fit carries no training geometry.", call. = FALSE)
  target_crs <- sf::st_crs(train)

  # ---- grid ----------------------------------------------------------------
  if (is.null(grid)) {
    grid <- .make_prediction_grid(sf::st_bbox(train), target_crs,
                                  cell_size = cell_size, n_cells = n_cells)
  } else {
    if (!inherits(grid, "sf"))
      stop("predict_surface(): `grid` must be an sf object.", call. = FALSE)
    # The boundary-clip path below guards its own empty result; without this
    # an empty `grid` reaches seq(1L, 0L, by = chunk_size) and aborts with
    # "wrong sign in 'by' argument", which names nothing the caller passed.
    if (nrow(grid) == 0L)
      stop("predict_surface(): `grid` has no rows; there is nothing to ",
           "predict onto.", call. = FALSE)
    grid <- ensure_projected(grid, target_crs = .crs_or_null(target_crs))
  }
  res <- attr(grid, "cell_size")

  # ---- clip ----------------------------------------------------------------
  if (!is.null(boundary)) {
    bnd <- ensure_projected(sf::st_geometry(boundary), target_crs = .crs_or_null(target_crs))
    keep <- lengths(sf::st_intersects(grid, bnd)) > 0L
    grid <- grid[keep, , drop = FALSE]
    if (nrow(grid) == 0L)
      stop("predict_surface(): no grid points fall inside `boundary`.",
           call. = FALSE)
  }

  # ---- covariates ----------------------------------------------------------
  preds <- object$predictor_vars
  missing_preds <- setdiff(preds, names(grid))
  if (length(missing_preds) > 0L) {
    if (is.null(covariates))
      stop("predict_surface(): the model uses predictor(s) ",
           paste(sQuote(missing_preds), collapse = ", "),
           " which are absent from the grid. Supply `covariates`, or a `grid` ",
           "that already carries them.", call. = FALSE)
    if (!inherits(covariates, "sf"))
      stop("predict_surface(): `covariates` must be an sf object.", call. = FALSE)
    # st_nearest_feature() against an empty layer fails inside sf with a
    # message that names neither argument.
    if (nrow(covariates) == 0L)
      stop("predict_surface(): `covariates` has no rows, so there is no ",
           "nearest feature to take predictor values from.", call. = FALSE)

    cov_missing <- setdiff(missing_preds, names(covariates))
    if (length(cov_missing) > 0L)
      stop("predict_surface(): `covariates` lacks column(s) ",
           paste(sQuote(cov_missing), collapse = ", "), ".", call. = FALSE)

    covariates <- .replay_crs_assumption(covariates, train, "predict_surface")
    covariates <- ensure_projected(covariates, target_crs = .crs_or_null(target_crs))
    nn  <- sf::st_nearest_feature(grid, covariates)
    cdf <- sf::st_drop_geometry(covariates)[nn, missing_preds, drop = FALSE]
    for (cn in missing_preds) grid[[cn]] <- cdf[[cn]]
  }

  # ---- chunked prediction --------------------------------------------------
  n <- nrow(grid)
  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq(1L, n, by = chunk_size)

  preds_vec <- rep(NA_real_, n)
  se_vec    <- if (isTRUE(se)) rep(NA_real_, n) else NULL
  se_ok     <- isTRUE(se)

  # With se = TRUE the point predictions can be read off the same draw matrix
  # instead of drawing the posterior a second time per chunk -- but only when
  # the caller left the summary at its default, since predict() returns
  # colMeans(draws) for summary = "mean" and the column medians otherwise.
  # Deriving unconditionally would silently change .pred for anyone passing
  # summary = "median".
  #
  # Test a PREFIX match, not `.dots$summary`.  `$` on a list partial-matches
  # only in the other direction, so list(summ = "median")$summary is NULL --
  # but R's argument matching does partial-match a supplied name onto a
  # formal, so predict_surface(fit, se = TRUE, summ = "median") is legal and
  # arrives at predict() as summary = "median".  The old test missed it and
  # .pred silently became colMeans(d): exactly what the paragraph above says
  # must not happen.
  .dots <- list(...)
  dot_names <- names(.dots)
  supplied_summary <- !is.null(dot_names) &&
    any(nzchar(dot_names) & startsWith("summary", dot_names))
  derive_pred <- se_ok && !supplied_summary

  for (s in starts) {
    e   <- min(s + chunk_size - 1L, n)
    idx <- s:e
    blk <- grid[idx, , drop = FALSE]

    got_pred <- FALSE
    if (se_ok) {
      d <- try(stats::predict(object, newdata = blk, draws = TRUE, ...),
               silent = TRUE)
      if (inherits(d, "try-error") || !is.matrix(d)) {
        se_ok <- FALSE          # backend has no draw interface; report once
      } else {
        se_vec[idx] <- apply(d, 2L, stats::sd)
        if (derive_pred) {
          preds_vec[idx] <- colMeans(d)
          got_pred <- TRUE
        }
      }
    }

    if (!got_pred) {
      p <- try(stats::predict(object, newdata = blk, ...), silent = TRUE)
      if (inherits(p, "try-error"))
        stop("predict_surface(): prediction failed on rows ", s, "-", e, ": ",
             as.character(p), call. = FALSE)
      preds_vec[idx] <- as.numeric(p)
    }
  }

  if (isTRUE(se) && !se_ok)
    .log_warn(paste0("predict_surface(): `se = TRUE` but this backend does not ",
                     "expose posterior draws; returning predictions only."))

  grid$.pred <- preds_vec
  if (isTRUE(se) && se_ok) grid$.pred_se <- se_vec

  attr(grid, "cell_size") <- res
  grid
}
