#' How far the held-out points sit from the training data
#'
#' Blocked cross-validation exists to put distance between a test point and
#' the training points that could tell you its value.  Whether it succeeded is
#' a measurement, and \code{\link{make_folds}()} only offers a proxy for it:
#' the block size compared against the estimated autocorrelation range,
#' reported as a warning.  That comparison is about the \emph{design}.  This
#' is about the \emph{result}: for every held-out point, the distance to its
#' nearest training point, summarised per fold.
#'
#' The distances do not say by themselves which way a score errs.  A
#' cross-validated error reads low when the held-out points sit closer to
#' their training data than the locations the model will predict at, and high
#' when they sit farther, so the same held-out distances can go with either,
#' depending on the target.  Random folds on evenly sited data, scored for a
#' map of the sampled area, are close to unbiased although nearly every
#' held-out point has a training neighbour inside the range; range-sized
#' blocks on the same data read high; random folds on clustered data scored
#' for a map of the whole region read low.  So without a prediction target
#' the function reports the distances and draws no conclusion.  Give it
#' \code{prediction_points} and it sets the held-out distances beside the
#' distances from those points to the same training data, and words its
#' reading by the signed gap between the two.
#'
#' Nothing is estimated: the distances come from the geometry, and the range,
#' when one is shown, is the one the folds already carry (\code{make_folds(
#' auto_range = TRUE)} records it) or the one you pass as \code{sac}.
#'
#' @section Reading the gap:
#' With a target, the held-out distances and the prediction distances are
#' compared as two distributions.  Each fold's prediction distances are
#' measured to that fold's own training rows, so that a five-fold design is
#' not read as too far only because each fold trains on four fifths of the
#' sites; the distances to all sites are reported beside them.  (With more
#' than 50 folds, the leave-one-out schemes, every training set is nearly the
#' whole sample and the prediction points are measured against all sites
#' once.)  Folds count by their share of the held-out points.
#'
#' Two numbers summarise the comparison: the signed difference of the means,
#' held-out minus prediction (the signed area between the two cumulative
#' distributions), and the Wasserstein distance \eqn{W_1} (the unsigned
#' area).  When a range is available both are also given on the correlation
#' scale, after each distance \eqn{d} is replaced by \eqn{1 - \rho(d)}, with
#' \eqn{\rho} the correlation the range implies (see \code{within_range}
#' under Value).  That scale counts a difference between two distances well
#' inside the range and ignores one between two distances beyond it.
#'
#' The reading is taken on the correlation scale when there is a range, with
#' a tolerance of 0.1, and otherwise on the mean distance, with a tolerance
#' of a quarter of the prediction mean.  (When the prediction points are
#' measured to all sites and every one lies on a site, as when the sites
#' themselves are passed as \code{prediction_points} to a leave-one-out
#' design, the prediction mean and with it that tolerance
#' are 0: held-out distances above 0 then read \code{"farther"}, and
#' held-out distances that are all 0 as well read \code{"matched"}.  With
#' 50 folds or fewer the points are measured to each fold's own training
#' rows, which leave out that fold's held-out sites, so the mean is above 0
#' and the quarter rule gives the reading.)
#' \describe{
#'   \item{\code{"closer"}}{The gap is below minus the tolerance and at least
#'     half of \eqn{W_1} in size: the held-out points sit closer to their
#'     training data than the prediction points do, and the cross-validated
#'     error tends to read low for that target.}
#'   \item{\code{"farther"}}{The gap is above the tolerance and at least
#'     half of \eqn{W_1}: the held-out points sit farther, and the error
#'     tends to read high.}
#'   \item{\code{"matched"}}{The gap is within the tolerance and
#'     \eqn{W_1} within twice the tolerance.}
#'   \item{\code{"mixed"}}{The two distributions cross, closer at some
#'     distances and farther at others (the gap is under half of
#'     \eqn{W_1}, or it is within the tolerance while \eqn{W_1} is not), so
#'     the distances give no single direction.}
#' }
#' The reading is a statement about nearest-neighbour distances, not a
#' measured bias.  A learner that draws on several neighbours can still score
#' low on clustered samples whose nearest-neighbour distances match.
#'
#' @param folds A \code{\link{make_folds}()} result, or its \code{$folds}
#'   element (a list of \code{train}/\code{test} splits).
#' @param data_sf The layer the folds were built on.  Row identifiers are
#'   matched through \code{..row_id} when the layer carries one, and by row
#'   position otherwise, which is what \code{make_folds()} and every
#'   \code{cv_*()} do.  As in \code{cv_*()}, a \code{make_folds()} result
#'   whose recorded rows sit at other locations here (folds built on another
#'   layer, such as the points before
#'   \code{\link{assign_features_to_polygons}()} dropped some) is refused;
#'   the location check is skipped when one of the two layers is POINT and
#'   the other is not.
#' @param sac Optional: an \code{\link{estimate_sac_range}()} result or a
#'   single number.  Defaults to the range the folds carry, if any.  Supplying
#'   one adds the \code{within_range} column and, with a prediction target,
#'   the gap on the correlation scale.  A
#'   range that records its CRS (the folds' own, or an
#'   \code{estimate_sac_range()} result) is compared with distances measured
#'   in that CRS, whatever CRS \code{data_sf} is in.  A bare number is taken
#'   to be in the units the distances are otherwise measured in: those of
#'   \code{data_sf} if it is projected, and for geographic (lon/lat) input
#'   metres, in the CRS \code{\link{ensure_projected}()} chooses (as for
#'   \code{make_folds()}'s \code{block_size}), not degrees.  A \code{units}
#'   object is refused.
#' @param prediction_points Optional \code{sf} or \code{sfc} layer of the
#'   locations the model will predict at: the grid from
#'   \code{\link{predict_surface}()}, or a grid over the region of interest
#'   (\code{sf::st_make_grid(region, what = "centers")}).  A non-POINT layer
#'   is reduced to one representative point per feature, as in
#'   \code{make_folds(method = "nndm")}, so pass grid cells or points, not
#'   one polygon of the region.  It is brought into the CRS the distances are
#'   measured in.  Default \code{NULL}: no target and no reading, except for
#'   NNDM and kNNDM folds (\code{method = "nndm"} and \code{"knndm"}), which
#'   carry the distances from the prediction points they were built for and
#'   are compared against those, measured to all sites.  The carried
#'   distances are lengths in the CRS the folds were built in
#'   (\code{folds$params$crs}).  When the distances here are measured in
#'   another CRS (\code{data_sf} is in another projected CRS, or \code{sac}
#'   records another), the two are not compared, \code{print()} says so, and
#'   passing \code{prediction_points} has the target measured here.
#' @return A data.frame of class \code{fold_separation}, one row per fold:
#'   \code{fold} (the fold's number: for the \code{$folds} of a
#'   \code{cv_*()} result, the \code{fold_id} its \code{fold_metrics} use,
#'   which differs from the list position once a fold has been dropped),
#'   \code{n_train}, \code{n_test}, \code{n_blocks} (\code{NA}
#'   for a scheme with no blocks), \code{min_dist} and \code{median_dist}
#'   (distance from a held-out point to its nearest training point, in the
#'   units of the CRS the \code{crs} attribute names), \code{within_range},
#'   \code{target_median_dist} (the median distance from the prediction
#'   points to the fold's training data) and \code{dist_gap} (the fold's mean
#'   held-out distance minus its mean prediction distance); the last two are
#'   \code{NA} without a target, and \code{print()} shows them under the
#'   headings \code{target} and \code{gap}.
#'
#'   \code{within_range} is a correlation-weighted share, \code{NA} without a
#'   range: the mean over the fold's held-out points of \eqn{\rho(d)}, the
#'   correlation the range leaves between a point and its nearest training
#'   point at distance \eqn{d}.  It is 1 for a point on top of training data
#'   and falls towards 0 beyond the range, so a point at a twentieth of the
#'   range counts for more than one at nine tenths of it.  \eqn{\rho(d)} is
#'   \eqn{\exp(-3 d / range)}, the exponential model whose effective range is
#'   the range given (0.05 at the range), or the spherical model when the
#'   range comes from an \code{estimate_sac_range()} fit that was spherical;
#'   the \code{range_model} attribute says which.  Before this version the
#'   column was the share of held-out points closer than the range, which
#'   counted every such point the same.
#'
#'   Attributes: \code{method}, \code{sac_range}, \code{range_model},
#'   \code{crs} (the CRS the distances were measured in),
#'   \code{n_unknown_ids} (the number of distinct rows the folds name that
#'   \code{data_sf} does not have, each counted once however many folds
#'   mention it, as in the \code{cv_*()} field of the same name),
#'   \code{pooled} (a list with the number, minimum, median and mean of the
#'   held-out distances over all folds, and the pooled \code{within_range})
#'   and \code{target}: \code{NULL} without a prediction target, otherwise a
#'   list with \code{n_prediction_points}, \code{source}, \code{basis}
#'   (\code{"each fold's training data"} or \code{"all sites"}),
#'   \code{held_out_mean}, \code{held_out_median}, \code{target_mean},
#'   \code{target_median}, \code{target_all_mean}, \code{target_all_median}
#'   (the prediction distances to all sites), \code{mean_diff}, \code{w1},
#'   \code{cor_gap}, \code{cor_w1} (\code{NA} without a range), \code{scale}
#'   (\code{"correlation"} or \code{"distance"}), \code{tolerance},
#'   \code{reading} (see \strong{Reading the gap}) and
#'   \code{nndm_out_of_reach} (for NNDM and kNNDM folds, the
#'   \code{params$target_out_of_reach} that \code{make_folds()} recorded;
#'   \code{NA} for other schemes).  \code{not_compared} is \code{NULL}, or,
#'   for NNDM or kNNDM folds whose carried target was not compared because the
#'   distances here are measured in another CRS, one sentence naming the two
#'   CRSs; \code{target} is then \code{NULL} and \code{print()} gives the
#'   sentence where it would otherwise say that no target is given.
#' @family cross-validation
#' @seealso \code{\link{make_folds}()} for the fold schemes and the block
#'   sizing this measures the outcome of; \code{\link{cv_block_size_sweep}()}
#'   for choosing a block size by cross-validated error instead.
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 200
#' pts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000)),
#'   coords = c("x", "y"), crs = 32632
#' )
#'
#' # Random folds put a training point close to every held-out one.
#' random  <- make_folds(pts, k = 4, method = "random_kfold", seed = 1)
#' print(fold_separation(random, pts, sac = 200))
#'
#' # Blocked folds hold out whole neighbourhoods, so the distances grow.
#' blocked <- make_folds(pts, k = 4, method = "block_kfold",
#'                       block_size = 250, seed = 1)
#' print(fold_separation(blocked, pts, sac = 200))
#'
#' # Which of the two suits a map of the sampled square?  Compare each with
#' # the distances from the cells of that map to the training data.
#' cells <- st_make_grid(pts, n = c(20, 20), what = "centers")
#' print(fold_separation(random, pts, sac = 200, prediction_points = cells))
#' print(fold_separation(blocked, pts, sac = 200, prediction_points = cells))
#' @export
fold_separation <- function(folds, data_sf, sac = NULL, prediction_points = NULL) {
  if (!inherits(data_sf, "sf"))
    stop("fold_separation(): `data_sf` must be an sf object.", call. = FALSE)
  splits <- if (is.list(folds) && !is.null(folds$folds)) folds$folds else folds
  if (!is.list(splits) || !length(splits) ||
      !all(vapply(splits, function(s) is.list(s) && !is.null(s$test), logical(1))))
    stop("fold_separation(): `folds` must be a make_folds() result or its ",
         "`$folds` element (a list of train/test splits).", call. = FALSE)
  # A bare split list has no `params`; `$` on it would partial-match nothing
  # and return NULL, which is what the lookups below rely on.
  params <- if (is.list(folds) && !is.null(folds$folds) && is.list(folds$params))
    folds$params else NULL

  # Same row-identity convention as make_folds() and every cv_*(): the stamped
  # column when there is one, row position otherwise.  Getting this wrong
  # would silently measure the distance between the wrong pairs of points.
  ids <- if ("..row_id" %in% names(data_sf)) data_sf[["..row_id"]] else
    seq_len(nrow(data_sf))
  # Which is what happened when the folds came from another layer: built on
  # prep_model_data()'s 200 points and measured on the 195 that
  # assign_features_to_polygons() kept, every row after the first dropped one
  # was paired with its neighbour's folds, and within_range read 1.0 against
  # 0.23-0.48 on the right layer, with nothing said.  cv_*() refuse such folds
  # through the probe make_folds() records; so does this.
  .check_fold_provenance(folds, data_sf, "fold_separation")

  # as.numeric() strips a `units` object to its number in whatever unit it
  # was written in, so 20 km became a range of 20 compared against metres.
  if (inherits(sac, "units"))
    stop("fold_separation(): `sac` must be a plain number in the units of the ",
         "CRS the distances are measured in (metres for lon/lat input); got ",
         format(sac), ".", call. = FALSE)
  sac_obj <- if (!is.null(sac)) sac else params$sac_range
  sac_val <- if (!is.null(sac_obj)) suppressWarnings(as.numeric(sac_obj)[1L]) else NA_real_
  if (!length(sac_val) || !is.finite(sac_val) || sac_val <= 0) sac_val <- NA_real_

  # A range is a length in the CRS it was measured in.  The one the folds
  # carry is in the folds' CRS (params$crs), and an estimate_sac_range()
  # result records its own; measuring the distances in whatever CRS data_sf
  # happened to arrive in compared a metre range with foot distances, off by
  # a factor of 3.3.  So when the range says where it is from, the distances
  # are measured there.  A bare number says nothing, and stays in data_sf's
  # (projected) units.
  rng_crs <- NULL
  if (is.finite(sac_val)) {
    cr <- if (!is.null(sac)) attr(sac, "crs")
          else attr(params$sac_range, "crs") %||% params$crs
    cr <- if (is.null(cr)) NULL
          else tryCatch(sf::st_crs(cr), error = function(e) NULL)
    if (!is.null(cr) && !is.na(cr) && !isTRUE(sf::st_is_longlat(cr)))
      rng_crs <- cr
  }

  pts <- data_sf
  if (!all(sf::st_geometry_type(pts, by_geometry = TRUE) == "POINT"))
    pts <- coerce_to_points(pts, "auto")
  pts <- if (is.null(rng_crs)) ensure_projected(pts) else
    .transform_or_stamp(pts, rng_crs, what = "data_sf", caller = "fold_separation")
  pts <- sf::st_zm(pts, drop = TRUE, what = "ZM")
  crs_lbl <- .fold_crs_label(pts)

  # A distance is only meaningful between finite coordinates.
  xy <- sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  ok <- stats::complete.cases(xy) & is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  if (!any(ok))
    stop("fold_separation(): `data_sf` has no usable coordinates.", call. = FALSE)

  # Each fold's rows in this layer.  A row the folds name and the layer does
  # not have appears once in one fold's test set and once in every other
  # fold's training set, so the unknown IDs are collected and counted once
  # each: five missing rows read as 15, 25 or 50 at k = 3, 5 or 10 (and as
  # 996 over 200 buffered leave-one-out folds) when the mentions were summed,
  # where the cv_*() field of the same name says 5.
  unknown_ids <- vector("list", length(splits))
  idx <- lapply(seq_along(splits), function(j) {
    s  <- splits[[j]]
    te <- match(s$test,  ids)
    tr <- match(s$train, ids)
    unknown_ids[[j]] <<- c(s$test[is.na(te)], s$train[is.na(tr)])
    list(te = te[!is.na(te) & ok[te]], tr = tr[!is.na(tr) & ok[tr]])
  })
  unknown <- length(unique(unlist(unknown_ids, use.names = FALSE)))

  # ---- The prediction target, when there is one ----------------------------
  # Supplied points are measured here.  NNDM and kNNDM folds carry the
  # distances from the prediction points they were built for, to all sites;
  # those serve when nothing is supplied and they were measured in the CRS
  # used here.
  pred <- NULL; d_all <- NULL; target_source <- NA_character_
  target_note <- NULL
  # The two distance-matched methods, and the name each goes by in a sentence.
  matched <- is.character(folds$method) && length(folds$method) == 1L &&
    folds$method %in% c("nndm", "knndm")
  matched_name <- if (matched && identical(folds$method, "knndm")) "kNNDM" else "NNDM"
  if (!is.null(prediction_points)) {
    pred <- .fold_sep_prediction_points(prediction_points, pts)
    sites <- sort(unique(unlist(lapply(idx, function(i) c(i$te, i$tr)),
                                use.names = FALSE)))
    if (length(sites)) {
      d_all <- .nn_dist_to(pred, pts[sites, ])
      d_all <- d_all[is.finite(d_all)]
    }
    target_source <- "prediction_points"
  } else if (matched && is.numeric(params$target_distances)) {
    if (identical(params$crs, crs_lbl)) {
      d_all <- params$target_distances[is.finite(params$target_distances)]
      target_source <- sprintf("the prediction points the %s folds were built for",
                               matched_name)
    } else {
      # The carried distances are lengths in the CRS the folds were built in.
      # Measured here in another one (`data_sf` arrived in it, or `sac`
      # records it), the held-out distances cannot be set beside them, and
      # the target was dropped without a word: the print then said "No
      # prediction target is given" of folds that were built for one.  The
      # reason is kept and printed in that sentence's place.
      crs_chr <- function(v) if (is.character(v) && length(v) == 1L && !is.na(v)) v
                             else "a layer with no recorded CRS"
      target_note <- sprintf(paste0(
        "These %s folds carry the distances from their %d prediction ",
        "points to the data in %s, and the distances here are measured in ",
        "%s, so the two are not compared."),
        matched_name,
        as.integer(params$n_prediction_points %||% length(params$target_distances)),
        crs_chr(params$crs), crs_chr(crs_lbl))
      .log_info(paste0("fold_separation(): %s Pass `prediction_points` to ",
                       "have the target measured here."), target_note)
    }
  }
  if (!length(d_all)) d_all <- NULL
  per_fold_target <- !is.null(pred) && !is.null(d_all) && length(splits) <= 50L

  cor_fn <- .fold_sep_cor_fn(sac_val, sac_obj)
  rho <- cor_fn$fn
  fold_blocks <- params$fold_blocks
  d_cv <- vector("list", length(splits)); d_tg <- vector("list", length(splits))

  rows <- lapply(seq_along(splits), function(j) {
    s  <- splits[[j]]
    te <- idx[[j]]$te; tr <- idx[[j]]$tr
    d <- if (length(te) && length(tr)) .nn_dist_to(pts[te, ], pts[tr, ]) else numeric(0)
    d <- d[is.finite(d)]
    d_cv[[j]] <<- d
    # What the fold's own model would face at the prediction points.
    tg <- if (is.null(d_all) || !length(tr)) numeric(0)
          else if (per_fold_target) {
            v <- .nn_dist_to(pred, pts[tr, ]); v[is.finite(v)]
          } else d_all
    if (per_fold_target) d_tg[[j]] <<- tg
    data.frame(
      # A cv_*() result's $folds carries each split's index in the original
      # fold list as fold_id, and its fold_metrics are labelled by it.  After
      # a fold is dropped the list position no longer matches, and labelling
      # by position paired each fold's error with another fold's distances.
      fold         = as.integer(s$fold_id %||% j),
      n_train      = length(tr),
      n_test       = length(te),
      n_blocks     = if (is.list(fold_blocks) && length(fold_blocks) >= j)
                       length(fold_blocks[[j]]) else NA_integer_,
      min_dist     = if (length(d)) min(d) else NA_real_,
      median_dist  = if (length(d)) stats::median(d) else NA_real_,
      within_range = if (length(d) && !is.null(rho)) mean(rho(d)) else NA_real_,
      target_median_dist = if (length(tg)) stats::median(tg) else NA_real_,
      dist_gap     = if (length(d) && length(tg)) mean(d) - mean(tg) else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (unknown > 0L)
    .log_info(paste0("fold_separation(): the folds name %d row(s) `data_sf` ",
                     "does not have, which were skipped. Are these the folds ",
                     "this layer was split into?"), unknown)

  cv_all <- unlist(d_cv, use.names = FALSE)
  pooled <- list(
    n            = length(cv_all),
    min_dist     = if (length(cv_all)) min(cv_all) else NA_real_,
    median_dist  = if (length(cv_all)) stats::median(cv_all) else NA_real_,
    mean_dist    = if (length(cv_all)) mean(cv_all) else NA_real_,
    within_range = if (length(cv_all) && !is.null(rho)) mean(rho(cv_all)) else NA_real_)

  target <- NULL
  if (!is.null(d_all) && length(cv_all)) {
    if (per_fold_target) {
      # Folds count by their share of the held-out points, as they do in a
      # pooled cross-validation score.
      n_cv <- lengths(d_cv); n_tg <- lengths(d_tg)
      use  <- which(n_cv > 0L & n_tg > 0L)
      tg_all <- unlist(d_tg[use], use.names = FALSE)
      tg_w   <- rep(n_cv[use] / n_tg[use], times = n_tg[use])
    } else {
      tg_all <- d_all; tg_w <- rep(1, length(d_all))
    }
    if (length(tg_all)) {
      cmp <- .fold_sep_compare(cv_all, tg_all, tg_w)
      cmp_cor <- if (is.null(rho)) list(diff = NA_real_, w1 = NA_real_)
                 else .fold_sep_compare(1 - rho(cv_all), 1 - rho(tg_all), tg_w)
      on_cor <- !is.null(rho)
      tol  <- if (on_cor) .fold_sep_tol_cor else .fold_sep_tol_dist
      # On the distance scale the tolerance is a share of the prediction
      # mean, so the gap is read in distance units against that share of the
      # mean.  It used to be divided by the mean first.  Prediction points on
      # top of the sites (a model scored where it was sampled) have a
      # prediction mean of 0: the gap was then Inf, or 0 / 0, the reading NA,
      # and the print said no target had been given.  With every prediction
      # distance 0 the tolerance is 0 and the held-out distances lie wholly
      # on one side: a held-out mean above 0 reads "farther", and one of 0
      # reads "matched".
      reading <- if (on_cor) .fold_sep_reading(cmp_cor$diff, cmp_cor$w1, tol)
                 else .fold_sep_reading(cmp$diff, cmp$w1, tol * cmp$mean_b)
      target <- list(
        n_prediction_points = if (!is.null(pred)) nrow(pred) else length(d_all),
        source = target_source,
        basis  = if (per_fold_target) "each fold's training data" else "all sites",
        held_out_mean = cmp$mean_a, held_out_median = stats::median(cv_all),
        target_mean = cmp$mean_b,
        target_median = if (per_fold_target) .fold_sep_wquantile(tg_all, tg_w, 0.5)
                        else stats::median(tg_all),
        target_all_mean = mean(d_all), target_all_median = stats::median(d_all),
        mean_diff = cmp$diff, w1 = cmp$w1,
        cor_gap = cmp_cor$diff, cor_w1 = cmp_cor$w1,
        scale = if (on_cor) "correlation" else "distance",
        tolerance = tol,
        reading = reading,
        # What make_folds(method = "nndm") or "knndm" found about its own
        # target: TRUE when the prediction points are nearer than
        # leave-one-out reaches, which no fold design can match.
        nndm_out_of_reach = if (matched)
          isTRUE(params$target_out_of_reach) else NA)
    }
  }

  structure(out, class = c("fold_separation", "data.frame"),
            method = folds$method %||% "supplied splits",
            sac_range = sac_val, range_model = cor_fn$model,
            n_unknown_ids = unknown, crs = crs_lbl,
            pooled = pooled, target = target,
            # Why a target the folds carry is not in `target` (NULL when
            # there is nothing to explain).  Not named `target_...`: attr()
            # matches partially, so attr(x, "target") would return this
            # sentence whenever `target` itself is NULL.
            not_compared = if (is.null(target)) target_note else NULL)
}


# Tolerances of the reading (see "Reading the gap" in ?fold_separation).
# Over 60 draws of random five-fold splits on evenly sited points (150 and
# 300 of them, ranges of 30 and 100 km on boxes 300 to 800 km wide), compared
# with a grid over the same box, the gap ran from -0.064 to +0.025 on the
# correlation scale and from -0.16 to +0.08 of the prediction mean on the
# distance scale, and every draw read "matched".  Random folds on clustered
# points against the same grids gave -0.15 or below and -0.44 or below (32
# of 32 "closer").
.fold_sep_tol_cor  <- 0.10
.fold_sep_tol_dist <- 0.25


#' Bring the prediction points into the CRS the distances are measured in
#'
#' @param prediction_points An sf/sfc layer, as passed.
#' @param pts The data points, already in the measuring CRS.
#' @return An sf POINT layer with finite coordinates.
#' @keywords internal
#' @noRd
.fold_sep_prediction_points <- function(prediction_points, pts) {
  if (inherits(prediction_points, "sfc"))
    prediction_points <- sf::st_sf(geometry = prediction_points)
  if (!inherits(prediction_points, "sf"))
    stop("fold_separation(): `prediction_points` must be an sf or sfc layer of ",
         "the locations to be predicted; got an object of class ",
         class(prediction_points)[1L], ".", call. = FALSE)
  if (!nrow(prediction_points))
    stop("fold_separation(): `prediction_points` has no rows.", call. = FALSE)
  if (!all(sf::st_geometry_type(prediction_points, by_geometry = TRUE) == "POINT")) {
    .log_info(paste0("fold_separation(): `prediction_points` are not POINTs; ",
                     "using one representative point per feature so the ",
                     "target distances are point-to-point."))
    prediction_points <- coerce_to_points(prediction_points, "auto")
  }
  if (is.null(.crs_or_null(pts))) {
    if (!is.na(sf::st_crs(prediction_points)))
      stop("fold_separation(): `data_sf` has no CRS and `prediction_points` has ",
           "one, so the distance between them cannot be measured. Set the CRS ",
           "of `data_sf` with sf::st_set_crs().", call. = FALSE)
    pred <- prediction_points
  } else {
    pred <- .transform_or_stamp(prediction_points, sf::st_crs(pts),
                                what = "prediction_points",
                                caller = "fold_separation")
  }
  pred <- sf::st_zm(pred, drop = TRUE, what = "ZM")
  empty <- sf::st_is_empty(pred)
  if (any(empty)) pred <- pred[!empty, , drop = FALSE]
  xy <- if (nrow(pred)) sf::st_coordinates(pred)[, 1:2, drop = FALSE] else matrix(0, 0, 2)
  keep <- is.finite(xy[, 1L]) & is.finite(xy[, 2L])
  if (!any(keep))
    stop("fold_separation(): `prediction_points` has no usable coordinates.",
         call. = FALSE)
  pred[keep, , drop = FALSE]
}


#' The correlation a range leaves at a distance
#'
#' The weight behind \code{within_range} and the correlation-scale gap: 1 at
#' distance 0, falling towards 0 beyond the range.  The range is an effective
#' range, so the exponential form is \code{exp(-3 d / range)}.  A range that
#' carries a single-structure spherical \code{variogram_model} (an
#' \code{estimate_sac_range()} fit) uses the spherical form instead; any
#' other model, and a bare number, is read as exponential.
#'
#' @param range A positive number, or \code{NA}.
#' @param sac_obj The object the range came from, for its model attribute.
#' @return A list with \code{fn} (\code{NULL} without a range) and
#'   \code{model} (\code{"Exp"}, \code{"Sph"} or \code{NA}).
#' @keywords internal
#' @noRd
.fold_sep_cor_fn <- function(range, sac_obj = NULL) {
  if (!is.finite(range) || range <= 0) return(list(fn = NULL, model = NA_character_))
  vm <- tryCatch(attr(sac_obj, "variogram_model"), error = function(e) NULL)
  fam <- "Exp"
  if (is.data.frame(vm) && "model" %in% names(vm)) {
    st <- setdiff(as.character(vm$model), "Nug")
    if (length(st) == 1L && identical(st, "Sph")) fam <- "Sph"
  }
  fn <- if (identical(fam, "Sph"))
    function(d) ifelse(d >= range, 0, 1 - 1.5 * (d / range) + 0.5 * (d / range)^3)
  else function(d) exp(-3 * d / range)
  list(fn = fn, model = fam)
}


#' Quantiles of a weighted sample (inverse of the weighted ECDF)
#'
#' @param x Values.
#' @param w Non-negative weights, one per value.
#' @param p Probabilities.
#' @return The smallest value whose cumulative weight share reaches each
#'   \code{p}.
#' @keywords internal
#' @noRd
.fold_sep_wquantile <- function(x, w, p) {
  o  <- order(x); x <- x[o]
  cw <- cumsum(w[o]); cw <- cw / cw[length(cw)]
  x[pmin(length(x), findInterval(p, cw, left.open = TRUE) + 1L)]
}


#' Signed mean difference and Wasserstein distance of two distance samples
#'
#' \code{diff} is \code{mean(a)} minus the weighted mean of \code{b}: the
#' signed area between the two cumulative distributions.  \code{w1} is the
#' unsigned area, the Wasserstein-1 distance, from the two quantile functions
#' on 2000 probability midpoints.  \code{w1} is never below \code{|diff|};
#' the two are equal when one distribution lies wholly to one side of the
#' other, and \code{w1} is the larger when they cross.
#'
#' @param a Unweighted values (the held-out distances).
#' @param b Values with weights \code{wb} (the prediction distances).
#' @return A list: \code{mean_a}, \code{mean_b}, \code{diff}, \code{w1}.
#' @keywords internal
#' @noRd
.fold_sep_compare <- function(a, b, wb = rep(1, length(b))) {
  p  <- (seq_len(2000L) - 0.5) / 2000
  qa <- .fold_sep_wquantile(a, rep(1, length(a)), p)
  qb <- .fold_sep_wquantile(b, wb, p)
  mean_a <- mean(a); mean_b <- stats::weighted.mean(b, wb)
  d <- mean_a - mean_b
  list(mean_a = mean_a, mean_b = mean_b, diff = d,
       w1 = max(mean(abs(qa - qb)), abs(d)))
}


#' Word a gap between held-out and prediction distances
#'
#' @param gap Signed gap, held-out minus prediction, on the reading's scale:
#'   the correlation scale, or distance units.
#' @param w1 The Wasserstein distance on the same scale.
#' @param tol The tolerance on that scale.  In distance units it is a share
#'   of the prediction mean and can be 0, when every prediction point lies
#'   on a site.
#' @return \code{"closer"}, \code{"farther"}, \code{"matched"},
#'   \code{"mixed"}, or \code{NA} when an argument is not finite.
#' @keywords internal
#' @noRd
.fold_sep_reading <- function(gap, w1, tol) {
  if (!is.finite(gap) || !is.finite(w1) || !is.finite(tol)) return(NA_character_)
  # Two identical distributions are matched whatever the tolerance.  With a
  # tolerance of 0 (prediction points on the sites, on the distance scale)
  # the test below would otherwise name a side for a gap of exactly 0.
  if (gap == 0 && w1 == 0) return("matched")
  # A direction needs a net shift past the tolerance that is also most of
  # the difference between the two distributions.  When they cross, the
  # signed area is the remainder of two opposite ones and names no side.
  if (abs(gap) >= tol && abs(gap) >= 0.5 * w1)
    return(if (gap < 0) "closer" else "farther")
  if (abs(gap) < tol && w1 < 2 * tol) return("matched")
  "mixed"
}


#' @export
print.fold_separation <- function(x, ...) {
  sac <- attr(x, "sac_range"); crs <- attr(x, "crs"); tg <- attr(x, "target")
  if (is.null(sac)) sac <- NA_real_
  if (is.null(crs)) crs <- NA_character_
  cat(sprintf("Fold separation: %s, %d fold(s), %d held-out point(s)%s\n",
              attr(x, "method") %||% "supplied splits", nrow(x),
              sum(x$n_test, na.rm = TRUE),
              if (is.na(crs)) "" else sprintf(" (%s)", crs)))
  # The unit of the CRS, not its identifier: "in EPSG:32617 units" named no
  # unit.  The CRS itself is on the line above.
  if (is.finite(sac)) {
    unit <- if (is.na(crs)) NA_character_ else .crs_unit_label(crs)
    cat(sprintf("  autocorrelation range: %s%s\n", format(sac, digits = 4),
                if (is.na(crs)) ""
                else if (!is.na(unit)) sprintf(" (in %s, like the distances)", unit)
                else sprintf(" (in the units of %s, like the distances)", crs)))
  }
  if (is.list(tg))
    cat(sprintf("  prediction target: %d point(s)%s, measured to %s\n",
                as.integer(tg$n_prediction_points),
                if (identical(tg$source, "prediction_points")) ""
                else sprintf(" (%s)", tg$source), tg$basis))
  num <- function(v) format(signif(v, 4))
  df <- as.data.frame(x)
  df$min_dist    <- signif(df$min_dist, 4)
  df$median_dist <- signif(df$median_dist, 4)
  # The two target columns are shown under shorter headings (`target` for
  # target_median_dist, `gap` for dist_gap), so that the table of a blocked
  # design with a range and a target fits an 80-column line.
  tgt <- df[["target_median_dist"]]; gap <- df[["dist_gap"]]
  df$target_median_dist <- NULL; df$dist_gap <- NULL
  if (!is.null(tgt) && !all(is.na(tgt))) df$target <- signif(tgt, 4)
  if (!is.null(gap) && !all(is.na(gap))) df$gap <- signif(gap, 4)
  # A mean correlation, so two decimals, not a percentage of points.
  wr <- df$within_range; df$within_range <- NULL
  if (!is.null(wr) && !all(is.na(wr)))
    df$within_range <- ifelse(is.na(wr), "NA", sprintf("%.2f", wr))
  if (all(is.na(df$n_blocks))) df$n_blocks <- NULL
  # A leave-one-out design has one row per point.  Printing hundreds of
  # one-point folds pushed the summary below off the screen; the object
  # itself keeps every row.
  shown <- if (nrow(df) > 20L) 10L else nrow(df)
  cat("\n"); print(df[seq_len(shown), , drop = FALSE], row.names = FALSE)
  if (shown < nrow(df))
    cat(sprintf("  ... and %d more fold(s); as.data.frame() shows every row\n",
                nrow(df) - shown))
  cat("\n")

  say <- function(...) cat(strwrap(paste0(...), width = 74, prefix = "  "), sep = "\n")
  # An object from before `pooled` was recorded still prints.
  pl <- attr(x, "pooled")
  if (!is.list(pl))
    pl <- list(min_dist = suppressWarnings(min(x$min_dist, na.rm = TRUE)),
               median_dist = NA_real_,
               within_range = stats::weighted.mean(x$within_range, x$n_test, na.rm = TRUE))
  # Every fold can be unmeasurable at once: no test rows anywhere, or none
  # with a usable coordinate.  Say there is nothing to report.
  if (!is.finite(pl$min_dist)) {
    say("No held-out point could be measured against its training data, so ",
        "there is no separation to report. Check n_test above: a fold with ",
        "no scorable rows is a finding about the fold scheme, not about the ",
        "model.")
    return(invisible(x))
  }
  say(if (is.finite(pl$median_dist))
        sprintf("Held-out points sit a median of %s from their nearest training point; the closest is %s away.",
                num(pl$median_dist), num(pl$min_dist))
      else sprintf("The closest held-out point is %s from its nearest training point.",
                   num(pl$min_dist)),
      if (is.finite(sac) && is.finite(pl$within_range))
        sprintf(paste0(" Under the range of %s that leaves a mean correlation ",
                       "of %.2f between a held-out point and its nearest ",
                       "training point (%s)."),
                format(sac, digits = 4), pl$within_range,
                if (identical(attr(x, "range_model"), "Sph"))
                  "the spherical correlation, averaged over the hold-out: 1 at distance 0, 0 at the range"
                else "exp(-3 d / range) averaged over the hold-out: 1 at distance 0, 0.05 at the range")
      else "")
  say(.fold_separation_advice(tg, attr(x, "method"), attr(x, "not_compared")))
  invisible(x)
}


#' The closing lines of print.fold_separation(): the target and the reading
#'
#' Without a prediction target the distances carry no direction, so nothing
#' is concluded and no remedy is named.  With one, the two sets of distances
#' are stated with their signed gap, and the reading is worded as a statement
#' about distances and the way a score tends to read for that target.
#'
#' "No prediction target is given" is said only when none exists.  A target
#' the folds carry that could not be compared here gets its reason
#' (\code{note}) in that sentence's place, and a target that was measured but
#' has no reading gets its numbers and a statement that none is drawn.  The
#' no-target sentence used to stand for all three, so it could appear under
#' a header naming the target, or for NNDM folds built for one.
#' @param target \code{attr(x, "target")}: a list, or \code{NULL}.
#' @param method The fold method (\code{attr(x, "method")}).
#' @param note \code{attr(x, "not_compared")}: why a target that exists is
#'   not in \code{target}, or \code{NULL}.
#' @return Character(1).
#' @keywords internal
#' @noRd
.fold_separation_advice <- function(target, method = NULL, note = NULL) {
  if (!is.list(target)) {
    if (is.character(note) && length(note) == 1L && !is.na(note) && nzchar(note))
      return(paste(note, "Pass `prediction_points` to have the target measured here."))
    return(paste("No prediction target is given, so these are distances and",
                 "no conclusion about the score: a cross-validated error reads",
                 "low when the hold-out sits closer to its training data than",
                 "the locations to be predicted do, and high when it sits",
                 "farther. Pass `prediction_points` to compare the two."))
  }
  rdg <- target$reading
  if (!(is.character(rdg) && length(rdg) == 1L && !is.na(rdg))) rdg <- NA_character_
  num <- function(v) format(signif(v, 4))
  # The sign is taken from the number as it prints, so a gap that rounds to
  # zero is not shown as "-0.00".
  signed <- function(v, fmt = NULL) {
    if (!is.finite(v)) return(format(v))
    shown <- if (is.null(fmt)) num(abs(v)) else sprintf(fmt, abs(v))
    zero  <- isTRUE(suppressWarnings(as.numeric(shown)) == 0)
    paste0(if (zero) "" else if (v > 0) "+" else "-", shown)
  }
  facts <- sprintf(paste0(
    "The prediction points sit a median of %s from %s. Held-out minus ",
    "prediction: %s in mean distance (W1 %s)%s."),
    num(target$target_median),
    if (identical(target$basis, "all sites")) "the nearest of all sites"
    else sprintf("the same training data (%s from the nearest of all sites)",
                 num(target$target_all_median)),
    signed(target$mean_diff), num(target$w1),
    if (is.finite(target$cor_gap %||% NA_real_))
      sprintf("; %s on the correlation scale (W1 %.2f)",
              signed(target$cor_gap, "%.2f"), target$cor_w1)
    else "")
  # Measured, and not readable: the numbers stand, and nothing is concluded
  # from them.  (Not reachable from fold_separation() as it stands, where
  # every gap is finite; kept so that an NA reading can never print as "no
  # target".)
  if (is.na(rdg))
    return(paste(facts, "The gap is not a finite number on the reading's",
                 "scale, so no reading is drawn from these distances."))
  reading <- switch(rdg,
    closer = paste(
      "The held-out points sit closer to their training data than the",
      "prediction points do, so the cross-validated error tends to read low",
      "for these prediction points."),
    farther = paste0(
      "The held-out points sit farther from their training data than the ",
      "prediction points do, so the cross-validated error tends to read high ",
      "for these prediction points.",
      if (identical(method, "nndm") && isTRUE(target$nndm_out_of_reach))
        paste0(" NNDM can only lengthen a held-out distance, so it cannot ",
               "match prediction points this near to the data.")
      else if (identical(method, "knndm") && isTRUE(target$nndm_out_of_reach))
        paste0(" A held-out point is never nearer to its training data than ",
               "to its nearest other sample, so kNNDM cannot match prediction ",
               "points this near to the data.")
      else ""),
    matched = paste(
      "The held-out points sit about as far from their training data as the",
      "prediction points do. Matched nearest-neighbour distances do not make",
      "a score unbiased by themselves: a learner that draws on several",
      "neighbours can still read low on clustered samples."),
    mixed = paste(
      "The two distributions cross: the held-out points sit closer than the",
      "prediction points at some distances and farther at others, so the",
      "distances give no single direction for the score."),
    "")
  paste(facts, reading)
}


#' The linear unit of a CRS, spelt for a sentence
#'
#' From a CRS label as \code{.fold_crs_label()} writes it (an
#' \code{AUTHORITY:CODE}, a proj string or a WKT), through
#' \code{sf::st_crs()$units_gdal}.
#' @param crs A CRS label.
#' @return Character(1), \code{NA} when the unit cannot be read.
#' @keywords internal
#' @noRd
.crs_unit_label <- function(crs) {
  u <- tryCatch(suppressWarnings(sf::st_crs(crs)$units_gdal),
                error = function(e) NULL)
  if (!(is.character(u) && length(u) == 1L && !is.na(u) && nzchar(u)))
    return(NA_character_)
  switch(u, metre = "metres", kilometre = "kilometres", foot = "feet",
         "US survey foot" = "US survey feet", degree = "degrees", u)
}


#' Refuse folds built on another layer, as the cv_*() functions do
#'
#' \code{.check_fold_probe()} for the fold consumers that are not
#' \code{cv_*()} (\code{fold_separation()}, \code{kriging_adequacy()}):
#' \code{data_sf} gets \code{..row_id} by position when it has none, as they
#' match fold entries.  When the probe was taken on POINT geometry and
#' \code{data_sf} is not, or the reverse, the locations cannot be compared (a
#' polygon's probe point is its centroid, a pointized copy's is whatever
#' point \code{coerce_to_points()} chose), so the check is skipped with an
#' INFO line rather than refusing a layer that works today.  A bare split list,
#' a label vector, or folds without a probe pass unchecked.
#' @param folds Whatever the caller was given as \code{folds}.
#' @param data_sf The layer the caller applies them to, as passed.
#' @param caller Name for the messages.
#' @return \code{invisible(NULL)}; called for the error.
#' @keywords internal
#' @noRd
.check_fold_provenance <- function(folds, data_sf, caller) {
  probe <- tryCatch(folds$params$row_probe, error = function(e) NULL)
  if (is.null(probe) || !inherits(data_sf, "sf") || nrow(data_sf) == 0L)
    return(invisible(NULL))
  if (!("..row_id" %in% names(data_sf)))
    data_sf$..row_id <- seq_len(nrow(data_sf))
  now_points <- all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT")
  if (is.logical(probe$points) && length(probe$points) == 1L &&
      !is.na(probe$points) && !identical(probe$points, now_points)) {
    .log_info(paste0("%s(): the supplied `folds` were built on %s geometry and ",
                     "this layer has %s geometry, so their row locations cannot ",
                     "be compared; skipping the provenance check."),
              caller, if (probe$points) "POINT" else "non-POINT",
              if (now_points) "POINT" else "non-POINT")
    return(invisible(NULL))
  }
  .check_fold_probe(folds, data_sf, caller)
}
