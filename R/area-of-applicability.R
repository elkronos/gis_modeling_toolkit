# =============================================================================
# Area of Applicability (Meyer & Pebesma 2021)
# =============================================================================
#
# Design note
# -----------
# Everything here is a pure function of numeric matrices.  There is no model
# fitting and no optional backend beyond FNN, which is injectable so the dense
# fallback is reachable on a machine that has FNN installed.  That is
# deliberate: the whole method is arithmetic, so all of it should be testable
# against hand-computed values rather than only against itself.
# =============================================================================


#' Extract a numeric predictor matrix
#'
#' Categorical predictors are refused rather than silently dummy-coded: a
#' one-hot column has no meaningful standard deviation, so scaling it puts an
#' arbitrary number into a Euclidean distance.
#'
#' @keywords internal
#' @noRd
.aoa_matrix <- function(data, vars, what = "data") {
  df <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else
    as.data.frame(data)
  missing_v <- setdiff(vars, names(df))
  if (length(missing_v) > 0L)
    stop(sprintf("area_of_applicability(): predictor(s) absent from %s: %s",
                 what, paste(sQuote(missing_v), collapse = ", ")),
         call. = FALSE)
  sub <- df[, vars, drop = FALSE]
  # Logicals are admitted: storage.mode() below gives them the natural 0/1
  # coding, whose standard deviation IS meaningful.  fit_rf_model() fits a
  # logical predictor, cv_rf() cross-validates it and predict() predicts with
  # it, so refusing the same model here would be self-inconsistent.  Only
  # factor/character predictors -- which would need an arbitrary one-hot
  # scaling -- are refused.
  bad <- vars[!vapply(sub, function(v) is.numeric(v) || is.logical(v),
                      logical(1))]
  if (length(bad) > 0L)
    stop(sprintf(paste0("area_of_applicability(): predictor(s) in %s are not ",
                        "numeric or logical: %s. The dissimilarity index is a ",
                        "Euclidean distance in scaled predictor space and has ",
                        "no definition for categorical predictors; recode them ",
                        "yourself if you have a scaling you can defend."),
                 what, paste(sQuote(bad), collapse = ", ")), call. = FALSE)
  m <- as.matrix(sub)
  storage.mode(m) <- "double"
  colnames(m) <- vars
  m
}


#' Attach coordinates as ordinary predictor columns
#'
#' Mirrors the \code{..x}/\code{..y} naming used by \code{.rf_frame()}, so a
#' forest fitted with \code{include_coords = TRUE} is measured in the same
#' space it was fitted in.
#'
#' @keywords internal
#' @noRd
.aoa_add_coords <- function(data_sf) {
  crd <- sf::st_coordinates(data_sf)
  if (ncol(crd) < 2L)
    stop("area_of_applicability(): could not extract X/Y coordinates.",
         call. = FALSE)
  data_sf[["..x"]] <- as.numeric(crd[, 1L])
  data_sf[["..y"]] <- as.numeric(crd[, 2L])
  data_sf
}


#' Centring and scaling derived from the TRAINING data only
#'
#' @return A list with \code{center}, \code{scale} and \code{keep} (predictors
#'   with usable variance).
#' @keywords internal
#' @noRd
.aoa_scaling <- function(X_train, tol = .Machine$double.eps^0.5) {
  ctr <- colMeans(X_train)
  scl <- apply(X_train, 2L, stats::sd)
  keep <- is.finite(ctr) & is.finite(scl) & scl > tol
  names(keep) <- colnames(X_train)
  list(center = ctr, scale = scl, keep = keep)
}


#' @keywords internal
#' @noRd
.aoa_apply_scaling <- function(X, sc) {
  k <- names(sc$keep)[sc$keep]
  if (length(k) == 0L)
    stop("area_of_applicability(): no predictor has usable variance in the ",
         "training data.", call. = FALSE)
  Z <- X[, k, drop = FALSE]
  # sweep() passes STATS through array(), which keeps names on atomic vectors,
  # so a named centring vector can leave a stray `names` attribute sitting on
  # the matrix alongside its dimnames.  Guard both ends.
  Z <- sweep(Z, 2L, unname(sc$center[k]), "-")
  Z <- sweep(Z, 2L, unname(sc$scale[k]),  "/")
  names(Z) <- NULL
  Z
}


#' Validate and align predictor weights
#'
#' The dissimilarity index is invariant to the overall scale of the weights --
#' both the nearest-neighbour distance and the normalising mean pairwise
#' distance carry the same factor -- so users need not normalise importance
#' values before passing them.  They are rescaled to mean 1 here purely so the
#' recorded values are readable.
#'
#' @param fill_vars Character vector of predictors that this package appended
#'   to \code{vars} itself (the \code{"..x"}/\code{"..y"} coordinate columns of
#'   a coordinate-using model).  The caller has never seen them, so a missing
#'   weight for one is filled in rather than raising an error, and unnamed
#'   weights may be one-per-\emph{user-visible}-predictor.
#' @keywords internal
#' @noRd
.aoa_weight_vector <- function(weights, vars, fill_vars = character(0)) {
  if (is.null(weights)) return(stats::setNames(rep(1, length(vars)), vars))
  if (!is.numeric(weights))
    stop("area_of_applicability(): `weights` must be numeric.", call. = FALSE)
  fill_vars    <- intersect(fill_vars, vars)
  user_vars    <- setdiff(vars, fill_vars)
  if (is.null(names(weights))) {
    # An unnamed vector may be one value per predictor the CALLER knows about;
    # the coordinate columns appended for an include_coords model are filled in
    # below rather than demanded from someone who never heard of them.
    if (length(weights) == length(vars)) {
      names(weights) <- vars
    } else if (length(fill_vars) > 0L && length(weights) == length(user_vars)) {
      names(weights) <- user_vars
    } else {
      stop(sprintf(paste0("area_of_applicability(): unnamed `weights` must ",
                          "have one value per predictor (%d%s), got %d. Name ",
                          "them to be unambiguous."),
                   length(vars),
                   if (length(fill_vars) > 0L)
                     sprintf(" including the coordinate columns, or %d without them",
                             length(user_vars)) else "",
                   length(weights)), call. = FALSE)
    }
  }
  missing_w <- setdiff(vars, names(weights))
  # A model fitted with include_coords = TRUE has "..x"/"..y" appended to
  # predictor_vars BEFORE this point, so a user who supplied importances for
  # the model's real predictors used to get a hard error naming two columns
  # they never heard of.  Default them instead to the mean of the weights that
  # WERE supplied (1 when none were), i.e. location counts as much as a typical
  # predictor -- which is the honest default for a model that splits on it.
  fill_missing <- intersect(missing_w, fill_vars)
  if (length(fill_missing) > 0L) {
    known <- weights[intersect(user_vars, names(weights))]
    known <- known[is.finite(known) & known >= 0]
    weights[fill_missing] <- if (length(known)) mean(known) else 1
    missing_w <- setdiff(missing_w, fill_missing)
  }
  if (length(missing_w) > 0L)
    stop("area_of_applicability(): `weights` has no entry for: ",
         paste(sQuote(missing_w), collapse = ", "), call. = FALSE)
  w <- weights[vars]
  if (anyNA(w) || any(!is.finite(w)) || any(w < 0))
    stop("area_of_applicability(): `weights` must be finite and non-negative. ",
         "Permutation importance is slightly negative for predictors that do ",
         "not help, so pass pmax(importance, 0).", call. = FALSE)
  if (all(w == 0))
    stop("area_of_applicability(): all `weights` are zero.", call. = FALSE)
  w * (length(w) / sum(w))
}


#' Minimum distance from each query row to any data row
#'
#' @param exclude Optional integer vector, one per query row, giving a row of
#'   \code{data} to exclude (used for the dense leave-one-out training
#'   reference). Supplying it forces the dense path; the FNN path has its own
#'   self-excluding entry point (\code{get.knn}) and never needs it.
#' @keywords internal
#' @noRd
.aoa_min_dist <- function(query, data,
                          use_fnn = requireNamespace("FNN", quietly = TRUE),
                          chunk_size = NULL, exclude = NULL) {
  nq <- nrow(query); nd <- nrow(data)
  if (nd < 1L)
    stop(".aoa_min_dist(): no reference rows.", call. = FALSE)
  if (nq < 1L) return(numeric(0))

  if (isTRUE(use_fnn) && is.null(exclude)) {
    d <- FNN::get.knnx(data = data, query = query, k = 1L)$nn.dist
    return(as.numeric(d[, 1L]))
  }

  # Dense fallback, chunked so the query-by-reference matrix stays bounded
  # regardless of grid size.  ||a-b||^2 = |a|^2 + |b|^2 - 2a.b loses precision
  # for near-coincident points, but the index is a ratio of that distance to
  # the mean pairwise training distance, so the absolute error is negligible
  # against the quantity it is divided by.
  if (is.null(chunk_size))
    chunk_size <- max(1L, min(10000L, as.integer(floor(4e6 / max(1L, nd)))))
  out  <- numeric(nq)
  dsq  <- rowSums(data^2)
  for (s in seq.int(1L, nq, by = chunk_size)) {
    e  <- min(s + chunk_size - 1L, nq)
    Q  <- query[s:e, , drop = FALSE]
    d2 <- outer(rowSums(Q^2), dsq, "+") - 2 * tcrossprod(Q, data)
    if (!is.null(exclude)) {
      ex <- exclude[s:e]
      ok <- !is.na(ex)
      if (any(ok)) d2[cbind(which(ok), as.integer(ex[ok]))] <- Inf
    }
    out[s:e] <- sqrt(pmax(apply(d2, 1L, min), 0))
  }
  out
}


#' Mean pairwise distance among training points (the DI normaliser)
#'
#' @keywords internal
#' @noRd
.aoa_normalizer <- function(X, max_n = 5000L, seed = 123L) {
  n <- nrow(X)
  if (n < 2L)
    stop("area_of_applicability(): at least two training rows are needed.",
         call. = FALSE)
  subsampled <- FALSE
  if (is.finite(max_n) && n > max_n) {
    cleanup <- .with_seed(seed)
    on.exit(cleanup(), add = TRUE)
    X <- X[sort(sample.int(n, as.integer(max_n))), , drop = FALSE]
    subsampled <- TRUE
    .log_info("area_of_applicability(): normaliser estimated from a %d-point subsample of %d training points.",
              as.integer(max_n), n)
  }
  m <- mean(as.numeric(stats::dist(X)))
  if (!is.finite(m) || m <= 0)
    stop("area_of_applicability(): the mean pairwise training distance is ",
         "zero or undefined; the training predictors carry no spread.",
         call. = FALSE)
  list(value = m, subsampled = subsampled, n_used = nrow(X))
}


#' Normalise a fold specification into train/test position lists
#'
#' Accepts a \code{make_folds()} result, a bare list of \code{train}/\code{test}
#' splits, or a plain vector of fold labels.  The \code{train} slot is used as
#' given rather than reconstructed as "everything else", so buffered and NNDM
#' folds contribute the training set that was actually available -- which is
#' the point of having built them.
#'
#' @param row_ids Optional vector of \code{..row_id} values, one per training
#'   row, in training-row order.  When supplied, the \code{train}/\code{test}
#'   entries of a list-shaped \code{folds} are treated as \code{..row_id}
#'   VALUES and resolved to row positions with \code{match()}; when
#'   \code{NULL} they are treated as positions, as before.  Fold labels are
#'   always positional -- they are one label per row by definition.
#' @keywords internal
#' @noRd
.aoa_fold_splits <- function(folds, n, row_ids = NULL) {
  if (is.null(folds)) return(NULL)

  if (is.atomic(folds) && !is.list(folds)) {
    if (length(folds) != n)
      stop(sprintf(paste0("area_of_applicability(): `folds` has %d labels but ",
                          "the training data has %d rows."),
                   length(folds), n), call. = FALSE)
    # droplevels() matters: a factor subset from a larger data set keeps its
    # unused levels, and each one would otherwise become an empty fold that
    # inflates the reported fold count and reaches .aoa_min_dist() with no
    # test rows.
    f <- droplevels(as.factor(folds))
    if (anyNA(f))
      stop("area_of_applicability(): `folds` contains missing labels.",
           call. = FALSE)
    if (nlevels(f) < 2L)
      stop("area_of_applicability(): `folds` must define at least two ",
           "non-empty folds.", call. = FALSE)
    # Fall through to the shared validation below rather than returning here,
    # so label-built splits get the same checks as hand-built ones.
    sp <- lapply(levels(f), function(lv) {
      te <- which(f == lv)
      list(test = te, train = setdiff(seq_len(n), te))
    })
    # which() already returns positions, so there is nothing to resolve.
    row_ids <- NULL
  } else {
    sp <- if (is.list(folds) && !is.null(folds$folds)) folds$folds else folds
  }

  if (!is.list(sp) || length(sp) == 0L)
    stop("area_of_applicability(): `folds` is not a recognised fold object.",
         call. = FALSE)
  shape_ok <- vapply(sp, function(z)
    is.list(z) && !is.null(z$test) && !is.null(z$train), logical(1))
  if (!all(shape_ok))
    stop("area_of_applicability(): `folds` must be a make_folds() result, a ",
         "list of train/test splits, or a vector of fold labels.",
         call. = FALSE)

  # Every make_folds() branch emits ..row_id VALUES in its train/test slots, as
  # its @return documents, and those coincide with row POSITIONS only when the
  # training data had no pre-existing ..row_id column.  cv_gwr(), cv_bayes()
  # and cv_spatial() all stamp one before prep_model_data(), so the documented
  # workflow -- "pass the same make_folds() result you passed to cv_spatial()"
  # -- hands us IDs.  Indexing X with them computed the training reference
  # distances for the WRONG rows whenever the IDs happened to land inside 1:n,
  # silently shifting the threshold.  Resolve, do not index.
  .resolve <- function(v, j, slot) {
    if (is.null(row_ids)) return(as.integer(v))
    pos <- match(v, row_ids)
    if (anyNA(pos))
      stop(sprintf(paste0("area_of_applicability(): fold %d's %s set names %d ",
                          "..row_id value(s) that are not in the training data ",
                          "(first: %s). `folds` was built on a different data ",
                          "set, or on rows that have since been removed."),
                   j, slot, sum(is.na(pos)),
                   format(v[is.na(pos)][1L])), call. = FALSE)
    as.integer(pos)
  }

  seen <- integer(n)
  out  <- vector("list", length(sp))
  for (j in seq_along(sp)) {
    te <- .resolve(sp[[j]]$test,  j, "test")
    tr <- .resolve(sp[[j]]$train, j, "train")
    if (anyNA(te) || anyNA(tr) || any(te < 1L | te > n) || any(tr < 1L | tr > n))
      stop(sprintf(paste0("area_of_applicability(): fold %d refers to rows ",
                          "outside 1:%d. `folds` was probably built on a ",
                          "different data set, or on data carrying its own ",
                          "..row_id."), j, n), call. = FALSE)
    if (length(te) == 0L)
      stop(sprintf("area_of_applicability(): fold %d has an empty test set.",
                   j), call. = FALSE)
    if (length(tr) == 0L)
      stop(sprintf("area_of_applicability(): fold %d has an empty training set.",
                   j), call. = FALSE)
    # A row in both slots would match itself at distance zero, dragging the
    # whole threshold to 0 and reporting almost everything as outside the AOA.
    # Silent, and the wrong direction, so it has to be fatal.
    both <- intersect(te, tr)
    if (length(both) > 0L)
      stop(sprintf(paste0("area_of_applicability(): fold %d has %d row(s) in ",
                          "both its train and test slots. Each would be its ",
                          "own nearest neighbour at distance zero, collapsing ",
                          "the threshold."), j, length(both)), call. = FALSE)
    clash <- te[seen[te] > 0L]
    if (length(clash) > 0L)
      stop(sprintf(paste0("area_of_applicability(): row %d appears in the test ",
                          "set of both fold %d and fold %d. The training ",
                          "reference distance would be ambiguous."),
                   clash[1L], seen[clash[1L]], j), call. = FALSE)
    seen[te] <- j
    out[[j]] <- list(test = te, train = tr)
  }
  if (any(seen == 0L))
    stop(sprintf(paste0("area_of_applicability(): %d training row(s) are in no ",
                        "test fold, so they have no reference distance."),
                 sum(seen == 0L)), call. = FALSE)
  out
}


#' Reference nearest-neighbour distances for the training data
#'
#' @keywords internal
#' @noRd
.aoa_train_dist <- function(X, splits = NULL,
                            use_fnn = requireNamespace("FNN", quietly = TRUE),
                            chunk_size = NULL) {
  n <- nrow(X)
  if (n < 2L)
    stop("area_of_applicability(): at least two training rows are needed.",
         call. = FALSE)

  if (is.null(splits)) {
    if (isTRUE(use_fnn))
      return(as.numeric(FNN::get.knn(X, k = 1L)$nn.dist[, 1L]))
    return(.aoa_min_dist(X, X, use_fnn = FALSE, chunk_size = chunk_size,
                         exclude = seq_len(n)))
  }

  out <- rep(NA_real_, n)
  for (sp in splits)
    out[sp$test] <- .aoa_min_dist(X[sp$test, , drop = FALSE],
                                  X[sp$train, , drop = FALSE],
                                  use_fnn = use_fnn, chunk_size = chunk_size)
  if (anyNA(out))
    stop("area_of_applicability(): some training rows received no reference ",
         "distance.", call. = FALSE)
  out
}


#' Outlier-removed maximum of the training dissimilarity index
#'
#' The threshold is the largest training DI that is not an upper outlier by the
#' usual rule, i.e. the largest value at or below \code{Q3 + 1.5 * IQR} -- the
#' "(outlier-removed) maximum" of Meyer & Pebesma (2021).  Note that CAST
#' obtains it via \code{grDevices::boxplot.stats()}, which uses Tukey's hinges
#' rather than the type-7 quantiles \code{stats::quantile()} returns; the two
#' agree closely but not exactly.
#'
#' @keywords internal
#' @noRd
.aoa_threshold <- function(di) {
  d <- di[is.finite(di)]
  if (length(d) == 0L)
    stop("area_of_applicability(): no finite training DI values.", call. = FALSE)
  fence <- stats::quantile(d, 0.75, names = FALSE) + 1.5 * stats::IQR(d)
  inside <- d[d <= fence]
  if (length(inside) == 0L) return(max(d))
  max(inside)
}


# -----------------------------------------------------------------------------
# Exported
# -----------------------------------------------------------------------------

#' Area of applicability of a spatial prediction model
#'
#' Computes the dissimilarity index (DI) of Meyer & Pebesma (2021) for a set of
#' prediction locations and flags those that fall inside the model's area of
#' applicability (AOA) -- the region of predictor space where the model's
#' cross-validated performance estimate can be expected to hold.
#'
#' @section Why a map alone is not enough:
#' A fitted model will return a number for any location you hand it, including
#' locations whose predictor values look nothing like anything it was trained
#' on. Those predictions are extrapolations dressed as interpolations, and a
#' cross-validation score says nothing about them, because the held-out folds
#' were drawn from the same predictor distribution as the training data. The
#' AOA marks where the score applies.
#'
#' @section How it is computed:
#' Predictors are centred and scaled using the training data's own means and
#' standard deviations, then optionally weighted by variable importance. For a
#' prediction point \eqn{p}, the DI is the distance to its nearest training
#' point in that space, divided by the mean pairwise distance among training
#' points. The same quantity is computed for the training data itself, using
#' each point's nearest neighbour \emph{outside its own cross-validation fold},
#' and the threshold is the largest training DI that is not an upper outlier.
#' Prediction points at or below that threshold are inside the AOA.
#'
#' The DI is invariant to the overall scale of \code{weights}: the numerator
#' and the normaliser carry the same factor. Importance values can be passed
#' as-is.
#'
#' @section The fold scheme changes the answer, and should:
#' With \code{folds = NULL} the training reference is each point's nearest
#' neighbour anywhere in the training data, which for clustered data is very
#' close, giving a small threshold and a conservative AOA. Passing the folds
#' you actually validated with makes the reference distances larger and the AOA
#' correspondingly wider. That is not a loophole -- the AOA is defined relative
#' to a performance estimate, and a spatially blocked estimate is a claim about
#' predicting further away. Pass the same \code{make_folds()} result you passed
#' to \code{\link{cv_spatial}}. Buffered and NNDM folds use the training set
#' they actually left available, not merely "everything outside the fold".
#'
#' @section Limitations:
#' Predictors must be numeric; categorical variables are refused rather than
#' silently dummy-coded. Predictors with zero variance in the training data are
#' dropped, and a prediction point taking a different value there is a form of
#' extrapolation this index cannot express. Without \code{weights} every
#' predictor counts equally, which overstates dissimilarity along directions
#' the model barely uses.
#'
#' @param newdata Prediction locations: an \code{sf} object (typically from
#'   \code{\link{predict_surface}}) or a data.frame, carrying the predictor
#'   columns.
#' @param model A fitted \code{spatial_fit}, supplying the training data and
#'   predictor names. Optional if \code{train_sf} and \code{predictor_vars} are
#'   given directly, which lets this be used with any model.
#' @param train_sf Training data, if not taken from \code{model}.
#' @param predictor_vars Predictor names, if not taken from \code{model}.
#' @param weights Optional named numeric vector of predictor importances. Any
#'   positive scale works. When \code{model} was fitted with
#'   \code{include_coords = TRUE}, the coordinates \code{"..x"} and
#'   \code{"..y"} are measured as predictors too (see below), and you are not
#'   expected to supply importances for them: any you leave out default to the
#'   mean of the weights you did supply, so location counts about as much as a
#'   typical predictor. Naming them explicitly overrides that. An unnamed
#'   vector may have one value per predictor either with or without the two
#'   coordinate columns.
#' @param folds Cross-validation folds: a \code{\link{make_folds}} result, a
#'   list of \code{train}/\code{test} splits, or a vector of fold labels with
#'   one entry per training row. Default \code{NULL} (plain nearest neighbour).
#' @param threshold Optional numeric override for the DI threshold.
#' @param normalizer_max_n Subsample the training data to this many points when
#'   computing the mean pairwise distance, which is quadratic. Default 5000.
#' @param seed Seed for that subsample. Default 123.
#' @param chunk_size Query rows per distance block on the dense path. Default
#'   \code{NULL} (chosen from the training size).
#' @param use_fnn Use \pkg{FNN} for nearest-neighbour search when available.
#'   Exposed so the dense fallback can be tested.
#'
#' @return An object of class \code{aoa}: a list with \code{aoa} (the input
#'   \code{newdata} with numeric \code{DI} and logical \code{AOA} columns
#'   added), \code{threshold}, \code{train_DI}, \code{normalizer},
#'   \code{weights}, \code{predictor_vars}, \code{dropped_vars}, counts, and
#'   \code{params}.
#'
#' @references
#' Meyer, H. and Pebesma, E. (2021). Predicting into unknown space? Estimating
#' the area of applicability of spatial prediction models. \emph{Methods in
#' Ecology and Evolution} 12(9), 1620--1633. \doi{10.1111/2041-210X.13650}
#'
#' @seealso \code{\link{predict_surface}} to build the grid,
#'   \code{\link{make_folds}} for the fold scheme.
#' @family cross-validation
#' @examples
#' \donttest{
#' library(sf)
#' set.seed(1)
#' n <- 120
#' train <- st_as_sf(
#'   data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
#'              a = rnorm(n), b = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' train$z <- 2 * train$a - train$b + rnorm(n, 0, 0.3)
#'
#' # Prediction points, some of them well outside the training predictor range
#' newpts <- st_as_sf(
#'   data.frame(x = runif(50, 0, 1000), y = runif(50, 0, 1000),
#'              a = c(rnorm(40), rnorm(10, 8)), b = rnorm(50)),
#'   coords = c("x", "y"), crs = 32632
#' )
#'
#' res <- area_of_applicability(newpts, train_sf = train,
#'                              predictor_vars = c("a", "b"))
#' res
#' table(res$aoa$AOA)
#' }
#' @export
area_of_applicability <- function(newdata, model = NULL, train_sf = NULL,
                                  predictor_vars = NULL, weights = NULL,
                                  folds = NULL, threshold = NULL,
                                  normalizer_max_n = 5000L, seed = 123L,
                                  chunk_size = NULL,
                                  use_fnn = requireNamespace("FNN",
                                                             quietly = TRUE)) {
  if (!is.null(model)) {
    if (!inherits(model, "spatial_fit"))
      stop("area_of_applicability(): `model` must be a spatial_fit.",
           call. = FALSE)
    if (is.null(train_sf))       train_sf       <- model$data_sf
    if (is.null(predictor_vars)) predictor_vars <- model$predictor_vars
  }
  if (is.null(train_sf) || is.null(predictor_vars))
    stop("area_of_applicability(): supply `model`, or both `train_sf` and ",
         "`predictor_vars`.", call. = FALSE)
  if (is.null(newdata) || is.null(nrow(newdata)) || nrow(newdata) == 0L)
    stop("area_of_applicability(): `newdata` has no rows.", call. = FALSE)
  if (is.null(nrow(train_sf)) || nrow(train_sf) < 2L)
    stop("area_of_applicability(): at least two training rows are needed to ",
         "define a reference distance.", call. = FALSE)

  predictor_vars <- unique(as.character(predictor_vars))
  if (length(predictor_vars) == 0L)
    stop("area_of_applicability(): no predictors.", call. = FALSE)

  # An explicit use_fnn = TRUE on a machine without FNN reached FNN::get.knnx()
  # and produced R's bare "there is no package called 'FNN'", with no hint that
  # use_fnn = FALSE is a working alternative.  The default is guarded by
  # requireNamespace(); an explicit request needs the same courtesy.
  if (isTRUE(use_fnn) && !requireNamespace("FNN", quietly = TRUE))
    stop("area_of_applicability(): `use_fnn = TRUE` but package 'FNN' is not ",
         "installed. Install it with install.packages(\"FNN\"), or pass ",
         "use_fnn = FALSE to use the dense fallback.", call. = FALSE)

  # A model fitted with the coordinates as predictors splits on location, so
  # the dissimilarity index has to measure location too.  Without this, a
  # prediction point far outside the training extent but with ordinary
  # covariate values reads as INSIDE the area of applicability -- exactly the
  # extrapolation the index exists to catch, and exactly the failure mode that
  # makes coordinate predictors dangerous in the first place.
  coord_model <- !is.null(model) && isTRUE(model$info$include_coords)
  coord_vars  <- if (coord_model) c("..x", "..y") else character(0)
  if (coord_model) {
    if (!inherits(train_sf, "sf") || !inherits(newdata, "sf"))
      stop("area_of_applicability(): this model uses the coordinates as ",
           "predictors, so both the training data and `newdata` must be sf ",
           "objects carrying geometry.", call. = FALSE)
    tr_crs <- sf::st_crs(train_sf)
    if (!is.na(tr_crs) && !identical(sf::st_crs(newdata), tr_crs))
      newdata <- sf::st_transform(newdata, tr_crs)
    train_sf <- .aoa_add_coords(train_sf)
    newdata  <- .aoa_add_coords(newdata)
    predictor_vars <- unique(c(predictor_vars, "..x", "..y"))
  }

  X_tr_full <- .aoa_matrix(train_sf, predictor_vars, "the training data")
  X_nw_full <- .aoa_matrix(newdata,  predictor_vars, "`newdata`")

  # Fold train/test slots from make_folds() are ..row_id VALUES, not row
  # positions; keep the IDs alongside so .aoa_fold_splits() can resolve them.
  tr_meta <- if (inherits(train_sf, "sf")) sf::st_drop_geometry(train_sf) else
    as.data.frame(train_sf)
  tr_row_ids <- if ("..row_id" %in% names(tr_meta)) tr_meta[["..row_id"]] else NULL

  # Training rows carrying NA or Inf cannot define a reference distance.
  # complete.cases() alone would let an Inf through, and it would then poison
  # that predictor's mean and sd, so .aoa_scaling() would drop the whole
  # predictor and report it as having "no variance" -- a false diagnosis that
  # silently removes a dimension the index is supposed to be measuring.
  tr_ok <- apply(X_tr_full, 1L, function(r) all(is.finite(r)))
  if (!all(tr_ok)) {
    if (sum(tr_ok) < 2L)
      stop("area_of_applicability(): fewer than two training rows with finite ",
           "values for every predictor.", call. = FALSE)
    .log_warn("area_of_applicability(): dropping %d training row(s) with missing or non-finite predictor values.",
              sum(!tr_ok))
    if (!is.null(folds))
      stop("area_of_applicability(): the training data has missing predictor ",
           "values, so `folds` row positions can no longer be matched. Remove ",
           "the incomplete rows before building folds.", call. = FALSE)
    X_tr_full <- X_tr_full[tr_ok, , drop = FALSE]
    # Unreachable when folds were supplied (the branch above stops), but keeps
    # the IDs aligned with X_tr_full for anything added later.
    if (!is.null(tr_row_ids)) tr_row_ids <- tr_row_ids[tr_ok]
  }

  sc <- .aoa_scaling(X_tr_full)
  dropped <- names(sc$keep)[!sc$keep]
  if (length(dropped) > 0L)
    .log_warn(paste0("area_of_applicability(): dropping predictor(s) with no ",
                     "variance in the training data: %s. A prediction point ",
                     "taking a different value there is extrapolation the ",
                     "dissimilarity index cannot represent."),
              paste(dropped, collapse = ", "))
  used_vars <- names(sc$keep)[sc$keep]

  Z_tr <- .aoa_apply_scaling(X_tr_full, sc)
  Z_nw <- .aoa_apply_scaling(X_nw_full, sc)

  # Weighting convention: the scaled COLUMN is multiplied by the weight, so a
  # weight w contributes w^2 to the squared distance.  This is what CAST --
  # the reference implementation, by the paper's own authors -- does, and
  # anyone cross-checking against it would otherwise find quietly different
  # numbers.  Multiplying by sqrt(w) instead (contributing w to the squared
  # distance) is the other defensible reading of "weighted Euclidean" and is
  # NOT what is used here.
  w_vec <- unname(w <- .aoa_weight_vector(weights, predictor_vars,
                                          fill_vars = coord_vars)[used_vars])
  Z_tr  <- sweep(Z_tr, 2L, w_vec, "*")
  Z_nw  <- sweep(Z_nw, 2L, w_vec, "*")

  splits <- .aoa_fold_splits(folds, nrow(Z_tr), row_ids = tr_row_ids)

  norm     <- .aoa_normalizer(Z_tr, max_n = normalizer_max_n, seed = seed)
  train_d  <- .aoa_train_dist(Z_tr, splits = splits, use_fnn = use_fnn,
                              chunk_size = chunk_size)
  train_DI <- train_d / norm$value

  thr <- if (is.null(threshold)) .aoa_threshold(train_DI) else {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold <= 0)
      stop("area_of_applicability(): `threshold` must be a single positive ",
           "number.", call. = FALSE)
    as.numeric(threshold)
  }

  # NA or Inf predictors in newdata give NA DI rather than a misleading number.
  # Same test as the training side above, deliberately: complete.cases() alone
  # lets an Inf through, and an Inf predictor then produces Inf - Inf = NaN in
  # the dense |a|^2 + |b|^2 - 2a.b path.  NaN counts as NA downstream so the
  # outcome was tolerable, but the asymmetry between the two checks was not
  # intentional and made the two paths disagree on what "usable" means.
  nw_ok <- apply(Z_nw, 1L, function(r) all(is.finite(r)))
  DI <- rep(NA_real_, nrow(Z_nw))
  if (any(nw_ok))
    DI[nw_ok] <- .aoa_min_dist(Z_nw[nw_ok, , drop = FALSE], Z_tr,
                               use_fnn = use_fnn,
                               chunk_size = chunk_size) / norm$value
  inside <- DI <= thr

  out <- newdata
  out$DI  <- DI
  out$AOA <- inside

  structure(
    list(
      aoa            = out,
      threshold      = thr,
      train_DI       = train_DI,
      normalizer     = norm$value,
      weights        = w,
      predictor_vars = used_vars,
      dropped_vars   = dropped,
      n_train        = nrow(Z_tr),
      n_new          = nrow(Z_nw),
      n_inside       = sum(inside, na.rm = TRUE),
      n_outside      = sum(!inside, na.rm = TRUE),
      n_na           = sum(is.na(DI)),
      params = list(
        folds_supplied      = !is.null(splits),
        n_folds             = if (is.null(splits)) 0L else length(splits),
        threshold_supplied  = !is.null(threshold),
        normalizer_max_n    = normalizer_max_n,
        normalizer_n_used   = norm$n_used,
        normalizer_subsampled = norm$subsampled,
        weights_supplied    = !is.null(weights),
        seed                = seed
      )
    ),
    class = "aoa"
  )
}


#' Print an area-of-applicability result
#'
#' Summarises where the model may be trusted: how many prediction locations
#' fall inside the area of applicability and how many outside, the
#' dissimilarity threshold that separated them, and the predictors the index
#' was computed over (naming any dropped for having no usable variance).  The
#' proportion outside is the headline number -- a map that extrapolates over
#' much of its extent is reporting predictions its training data cannot
#' support, whatever the cross-validation score said.
#'
#' @param x An \code{aoa} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.aoa <- function(x, ...) {
  cat("Area of applicability (Meyer & Pebesma 2021)\n\n")
  cat(sprintf("  predictors  : %d (%s)\n", length(x$predictor_vars),
              paste(x$predictor_vars, collapse = ", ")))
  if (length(x$dropped_vars) > 0L)
    cat(sprintf("  dropped     : %s (no training variance)\n",
                paste(x$dropped_vars, collapse = ", ")))
  cat(sprintf("  weighted    : %s\n",
              if (isTRUE(x$params$weights_supplied)) "yes" else
                "no (all predictors count equally)"))
  cat(sprintf("  training    : %d points\n", x$n_train))
  cat(sprintf("  reference   : %s\n",
              if (isTRUE(x$params$folds_supplied))
                sprintf("nearest point outside each of %d CV folds",
                        x$params$n_folds)
              else "nearest other training point (no folds supplied)"))
  cat(sprintf("  normaliser  : %.4f (mean pairwise distance%s)\n",
              x$normalizer,
              if (isTRUE(x$params$normalizer_subsampled))
                sprintf(", %d-point subsample", x$params$normalizer_n_used)
              else ""))
  cat(sprintf("  threshold   : %.4f%s\n", x$threshold,
              if (isTRUE(x$params$threshold_supplied)) " (supplied)" else
                " (outlier-removed max of training DI)"))
  cat("\n")

  pct <- if (x$n_new > 0L) 100 * x$n_inside / x$n_new else NA_real_
  cat(sprintf("  %d of %d prediction points inside the AOA (%.1f%%)\n",
              x$n_inside, x$n_new, pct))
  if (x$n_na > 0L)
    cat(sprintf("  %d with missing predictors (DI = NA)\n", x$n_na))
  if (x$n_outside > 0L)
    cat("\nPredictions outside the AOA are extrapolations; the ",
        "cross-validated\nperformance estimate does not cover them.\n",
        sep = "")
  invisible(x)
}
