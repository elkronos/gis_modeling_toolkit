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
# Ceiling on the number of prediction cells.  A surface of more than a few
# million points is not something anyone plots or writes out; a request for one
# is a units mistake, and refusing it in milliseconds beats failing to allocate
# after minutes.
.surface_max_cells <- 5e6

.make_prediction_grid <- function(bb, crs, cell_size = NULL, n_cells = 10000L) {
  w <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  h <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0)
    stop("predict_surface(): the fitted data has no finite extent.", call. = FALSE)

  if (is.null(cell_size)) {
    # as.integer() silently returns NA above 2^31, and the error below then
    # blamed `cell_size` -- an argument the caller did not pass.  Validate the
    # argument the caller actually gave.
    if (!is.numeric(n_cells) || length(n_cells) != 1L || !is.finite(n_cells) ||
        n_cells < 1)
      stop("predict_surface(): `n_cells` must be a single positive finite number.",
           call. = FALSE)
    if (n_cells > .surface_max_cells)
      stop(sprintf(paste0("predict_surface(): `n_cells` = %s is above the %s ",
                          "cells this function will build."),
                   format(n_cells, big.mark = ",", scientific = FALSE),
                   format(.surface_max_cells, big.mark = ",",
                          scientific = FALSE)), call. = FALSE)
    cell_size <- sqrt((w * h) / n_cells)
  }
  if (!is.numeric(cell_size) || length(cell_size) != 1L ||
      !is.finite(cell_size) || cell_size <= 0)
    stop("predict_surface(): `cell_size` must be a positive finite number.",
         call. = FALSE)

  # Refuse an absurd grid before expand.grid() builds it.  The classic units
  # mistake -- a cell_size in kilometres on a metre CRS -- asked for 3.7e10
  # cells and R aborted trying to allocate 139.6 GB (or swapped the machine on
  # a bigger one).  create_grid_polygons() guards the identical mistake and
  # names the CRS units in its message; this builder is the sibling that did
  # not.
  n_est <- ceiling(w / cell_size) * ceiling(h / cell_size)
  if (is.finite(n_est) && n_est > .surface_max_cells) {
    unit_lbl <- tryCatch({
      u <- sf::st_crs(crs)$units_gdal
      if (is.null(u) || is.na(u) || !nzchar(u)) "CRS units" else u
    }, error = function(e) "CRS units")
    stop(sprintf(paste0("predict_surface(): a cell size of %s on an extent of ",
                        "%s x %s would produce about %s cells, above the %s ",
                        "this function will build. Check that `cell_size` is in ",
                        "the CRS units of the fitted data (%s), or raise ",
                        "`n_cells` instead."),
                 format(signif(cell_size, 6)), format(signif(w, 4)),
                 format(signif(h, 4)),
                 format(n_est, big.mark = ",", scientific = FALSE),
                 format(.surface_max_cells, big.mark = ",", scientific = FALSE),
                 unit_lbl),
         call. = FALSE)
  }

  # seq(from, to, by) ERRORS with "wrong sign in 'by' argument" when from > to,
  # rather than returning length 0 -- so a cell_size wider than the extent used
  # to abort with an internal-looking message, and the length-0 guards below
  # were unreachable.  Build the axis explicitly instead.
  .axis <- function(lo, hi) {
    lo <- as.numeric(lo); hi <- as.numeric(hi)
    # A relative tolerance on the floor.  When the extent is an exact multiple
    # of the cell size the ratio can land a hair below the integer -- 0.3 / 0.1
    # is 2.9999999999999996 -- and a whole column went missing, leaving a
    # cell-wide strip uncovered: 1197 of 10000 random squares at the default
    # n_cells came out 99 columns wide instead of 100.
    r <- (hi - lo) / cell_size
    n <- floor(r + sqrt(.Machine$double.eps) * max(1, r))
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
#' constructed by hand, which makes the most common downstream task (produce
#' a map) more work than it should be.  This wraps the grid construction,
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
#'   at least one row.  It is brought into the fit's CRS first: a CRS-less grid
#'   is given the interpretation the training data got (the assumption recorded
#'   on the fit), with a warning, and then reprojected.  Otherwise a CRS-less
#'   grid can land thousands of kilometres from the covariates and every cell
#'   takes the same nearest feature.  A grid of polygons
#'   (\code{\link{create_grid_polygons}()} output, say) is reduced to one
#'   representative point per cell, as \code{\link{coerce_to_points}()} does,
#'   so covariates are taken at the location predicted for; \code{boundary}
#'   then keeps the cells whose point falls inside it.
#' @param cell_size Grid resolution in CRS units.  Ignored when \code{grid} is
#'   supplied; when \code{NULL}, derived from \code{n_cells}.  A value that
#'   would produce more than 5,000,000 cells is refused, naming the implied
#'   count and the CRS units.  The usual cause is a value in the wrong unit.  A
#'   \code{cell_size} wider than the extent yields a single centred cell.
#' @param n_cells Approximate cell count used to derive \code{cell_size}.
#'   Default 10000.  Must be a single positive finite number and at most
#'   5,000,000; anything else is an error.  Also ignored when \code{grid} is
#'   supplied.  The grid you pass is used verbatim.
#' @param boundary Optional polygonal \code{sf}/\code{sfc}; grid points outside
#'   it are dropped.  Put through the same CRS replay and reprojection as
#'   \code{grid}.
#' @param covariates Optional \code{sf} layer carrying the model's predictors.
#'   Required when the model has predictors and \code{grid} does not already
#'   contain them.  Values are taken from the nearest feature.
#' @param chunk_size Rows per prediction call. Default 5000.  A pure
#'   performance knob for the GWR and random-forest backends, whose rows do not
#'   interact.  For a \code{bayesian_fit} it is also that, \emph{provided} the
#'   grid stays inside the training extent.  Beyond it the GP boundary has to
#'   grow and predictions depend on which rows share the call; see
#'   \code{\link{predict.bayesian_fit}}.
#' @param se Logical; also return a standard-error/posterior-SD column where the
#'   backend supports it.  Default FALSE.  For a \code{bayesian_fit} this is
#'   the SD of the posterior draws \code{predict()} returns, and those are of
#'   the expected value by default (\code{type = "epred"}): the uncertainty
#'   of the mean surface, not of a new observation, which also carries the
#'   observation noise.  For the predictive SD, the one that goes with
#'   prediction intervals and \code{cv_bayes()}'s calibration, pass
#'   \code{type = "predict"} as well.
#' @param ... Passed to \code{predict()}, e.g. \code{type = "predict"} for a
#'   \code{bayesian_fit}.  Not \code{draws}, which this function sets itself
#'   and refuses here.
#' @return An \code{sf} POINT layer with a \code{.pred} column (and
#'   \code{.pred_se} when \code{se = TRUE} and available; one a supplied
#'   \code{grid} already carried, from an earlier surface, is removed
#'   otherwise).  For an
#'   auto-generated grid the resolution is attached as attribute
#'   \code{"cell_size"}.  For a user-supplied \code{grid} it is only whatever
#'   \code{"cell_size"} attribute that object already carried.  That is usually
#'   \code{NULL}, and \code{NULL} for certain if the grid had to be
#'   re-projected, since \code{st_transform()} does not preserve custom
#'   attributes.  The resolution of a grid you built is not this function's to
#'   infer.
#' @family prediction
#' @examples
#' # Any spatial_fit works here; a forest keeps the example free of the
#' # optional GWR/Stan backends.
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   pts <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   pts$price <- 10 + 0.01 * (st_coordinates(pts)[, 1] - 5e5) +
#'     2 * pts$elev + rnorm(n)
#'   fit  <- fit_rf_model(pts, "price", "elev", num_trees = 100, seed = 1)
#'   surf <- predict_surface(fit, n_cells = 500, covariates = pts)
#'   print(surf[".pred"])        # one prediction per grid cell, as an sf layer
#'   # Check where that surface is extrapolating before mapping it.  The grid
#'   # took its covariates from the nearest observation, so here nothing is
#'   # outside; a grid with its own covariate raster is where this bites.
#'   area_of_applicability(surf, model = fit)
#' }
#' @export
predict_surface <- function(object, grid = NULL, cell_size = NULL,
                            n_cells = 10000L, boundary = NULL,
                            covariates = NULL, chunk_size = 5000L,
                            se = FALSE, ...) {
  if (!inherits(object, "spatial_fit"))
    stop("predict_surface(): `object` must be a spatial_fit.", call. = FALSE)

  # `draws` is this function's to set: it asks for the draw matrix itself when
  # se = TRUE.  Passed through `...` it reached predict() too, and a backend
  # that honours it returned an n_draws x n matrix that as.numeric() flattened
  # into .pred column by column -- cell 1's draws, then cell 2's -- with only
  # a cryptic length warning; with se = TRUE the duplicated argument failed
  # inside try() and was reported as a backend without draws.  A prefix
  # counts, since R's argument matching would complete it.
  dot_nms <- names(list(...))
  if (!is.null(dot_nms) && any(nzchar(dot_nms) & startsWith("draws", dot_nms)))
    stop("predict_surface(): `draws` cannot be passed through `...`; ",
         "predict_surface() requests the posterior draws itself when se = TRUE ",
         "and returns their SD as .pred_se. For the draw matrix, call ",
         "predict(object, newdata = grid, draws = TRUE).", call. = FALSE)

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
    # The same replay `covariates` gets, and for exactly the same reason: a
    # CRS-less grid whose OWN bounding box does not trip the lon/lat heuristic
    # had the fit's projected CRS stamped onto raw degrees, landing 5,300 km
    # from the correctly-reprojected covariates.  st_nearest_feature() then
    # handed every cell the same covariate row and the whole surface collapsed
    # to one constant -- silently, with no error anywhere.
    grid <- .replay_crs_assumption(grid, train, "predict_surface", "grid")
    grid <- ensure_projected(grid, target_crs = .crs_or_null(target_crs))
    # A polygon grid -- create_grid_polygons() output, say -- was used as it
    # was.  st_nearest_feature() then gave each cell whichever covariate point
    # inside it the spatial index returned first, not the one at its centre,
    # while predict() pointized the cell by itself, so the covariates and the
    # location predicted at no longer matched: predictions off by up to 2.6 on
    # a 0-30 response, and a row shuffle of `covariates` moved them by up to
    # 4.7.  Reduce it to representative points first, as
    # area_of_applicability() does, so the surface is the POINT layer the
    # manual promises.
    if (!all(sf::st_geometry_type(grid, by_geometry = TRUE) == "POINT"))
      grid <- coerce_to_points(grid, "auto")
  }
  res <- attr(grid, "cell_size")

  # ---- clip ----------------------------------------------------------------
  if (!is.null(boundary)) {
    bnd <- .replay_crs_assumption(sf::st_geometry(boundary), train,
                                  "predict_surface", "boundary")
    bnd <- ensure_projected(bnd, target_crs = .crs_or_null(target_crs))
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

    covariates <- .replay_crs_assumption(covariates, train, "predict_surface",
                                        "covariates")
    covariates <- ensure_projected(covariates, target_crs = .crs_or_null(target_crs))
    nn  <- sf::st_nearest_feature(grid, covariates)
    cdf <- sf::st_drop_geometry(covariates)[nn, missing_preds, drop = FALSE]
    for (cn in missing_preds) grid[[cn]] <- cdf[[cn]]
  }

  # ---- chunked prediction --------------------------------------------------
  n <- nrow(grid)
  # Inf is the natural way to ask for "one chunk, do not split", and it is
  # exactly what as.integer() turns into NA -- after which seq(by = NA) fails
  # with "invalid '(to - from)/by'", naming nothing the caller passed.
  .check_scalar(chunk_size, "chunk_size", "predict_surface", min = 1,
                max = .Machine$integer.max, what = "a single positive number")
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
      # One value per row, or the assignment below recycles or truncates
      # whatever came back into .pred without a word that means anything.
      if (length(p) != length(idx))
        stop(sprintf(paste0("predict_surface(): predict() returned %d value(s) ",
                            "for the %d rows %d-%d; expected one per row. Check ",
                            "what the backend's predict() returns for the ",
                            "arguments passed through `...`."),
                     length(p), length(idx), s, e), call. = FALSE)
      preds_vec[idx] <- as.numeric(p)
    }
  }

  if (isTRUE(se) && !se_ok)
    .log_warn(paste0("predict_surface(): `se = TRUE` but this backend does not ",
                     "expose posterior draws; returning predictions only."))

  grid$.pred <- preds_vec
  # A grid that is an earlier surface carries that model's .pred_se.  It was
  # kept whenever this call did not replace it -- beside the new .pred, and
  # even after the log said "returning predictions only".
  grid$.pred_se <- if (isTRUE(se) && se_ok) se_vec else NULL

  attr(grid, "cell_size") <- res
  grid
}
