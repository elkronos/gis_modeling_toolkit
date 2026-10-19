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
#' Categorical predictors are refused here instead of being silently
#' dummy-coded: a one-hot column has no meaningful standard deviation, so
#' scaling it by one puts an arbitrary number into a Euclidean distance.
#' \code{.aoa_design()} is where a caller who asked for
#' \code{categorical = "indicator"} or \code{"drop"} gets them handled, with
#' a scaling that is stated; this function only ever sees the numeric ones.
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
                        "no definition of its own for categorical predictors. ",
                        "Pass categorical = \"indicator\" to code each level as ",
                        "a 0/1 column (a level no training row has is then ",
                        "outside the AOA), or categorical = \"drop\" to compute ",
                        "the index on the numeric predictors alone; or recode ",
                        "them yourself if you have a scaling you can defend."),
                 what, paste(sQuote(bad), collapse = ", ")), call. = FALSE)
  m <- as.matrix(sub)
  storage.mode(m) <- "double"
  colnames(m) <- vars
  m
}


#' Is a column categorical for the purposes of the index?
#' @keywords internal
#' @noRd
.aoa_is_categorical <- function(v) is.factor(v) || is.character(v)


#' Training and prediction matrices, with categorical predictors handled
#'
#' Numeric and logical predictors go through \code{.aoa_matrix()} unchanged.
#' A factor or character predictor (which is categorical is read from the
#' TRAINING data) is treated as \code{categorical} says:
#' \describe{
#'   \item{\code{"refuse"}}{\code{.aoa_matrix()}'s error, which names the two
#'     alternatives.}
#'   \item{\code{"drop"}}{left out of both matrices.}
#'   \item{\code{"indicator"}}{one 0/1 column per level that occurs in the
#'     training rows, named \code{"<predictor>=<level>"}.  A prediction row
#'     whose level no training row has gets 0 in every one of them and is
#'     flagged in \code{unseen}, so the caller can mark it outside.}
#' }
#'
#' @return A list: \code{X_tr}, \code{X_nw} (numeric matrices with the same
#'   columns), \code{col_var} (the predictor each column came from),
#'   \code{indicator} (logical, one per column), \code{cat_vars} (the
#'   categorical predictors found), \code{vars} (the predictors still in the
#'   index) and \code{unseen} (a logical matrix, prediction rows by indicator
#'   predictors: the row's level is not among the training levels).
#' @keywords internal
#' @noRd
.aoa_design <- function(train, newdata, vars,
                        categorical = c("refuse", "indicator", "drop")) {
  categorical <- match.arg(categorical)
  tr <- if (inherits(train, "sf")) sf::st_drop_geometry(train) else
    as.data.frame(train)
  nw <- if (inherits(newdata, "sf")) sf::st_drop_geometry(newdata) else
    as.data.frame(newdata)
  for (side in list(list(tr, "the training data"), list(nw, "`newdata`"))) {
    missing_v <- setdiff(vars, names(side[[1L]]))
    if (length(missing_v) > 0L)
      stop(sprintf("area_of_applicability(): predictor(s) absent from %s: %s",
                   side[[2L]], paste(sQuote(missing_v), collapse = ", ")),
           call. = FALSE)
  }
  is_cat   <- vapply(tr[, vars, drop = FALSE], .aoa_is_categorical, logical(1))
  cat_vars <- vars[is_cat]
  num_vars <- vars[!is_cat]
  none <- matrix(logical(0), nrow = nrow(nw), ncol = 0L)
  # Numeric when the index was defined, categorical now: no `categorical`
  # setting can line the two up, so the advice .aoa_matrix() gives would
  # send the caller round in a circle.
  flipped <- num_vars[vapply(nw[, num_vars, drop = FALSE], .aoa_is_categorical,
                             logical(1))]
  flipped <- flipped[vapply(tr[, flipped, drop = FALSE], function(v)
    is.numeric(v) || is.logical(v), logical(1))]
  if (length(flipped) > 0L)
    stop(sprintf(paste0("area_of_applicability(): predictor(s) %s are numeric ",
                        "or logical in the training data but categorical in ",
                        "`newdata` (%s); convert the column(s) first."),
                 paste(sQuote(flipped), collapse = ", "),
                 paste(vapply(flipped, function(v) class(nw[[v]])[1L],
                              character(1)), collapse = ", ")), call. = FALSE)

  if (length(cat_vars) == 0L || identical(categorical, "refuse")) {
    # With nothing categorical this is exactly the two .aoa_matrix() calls
    # the function always made; with something categorical and "refuse" the
    # first of them stops, naming the predictor and the alternatives.
    X_tr <- .aoa_matrix(tr, vars, "the training data")
    X_nw <- .aoa_matrix(nw, vars, "`newdata`")
    return(list(X_tr = X_tr, X_nw = X_nw, col_var = vars,
                indicator = rep(FALSE, length(vars)),
                cat_vars = character(0), vars = vars, unseen = none))
  }

  if (identical(categorical, "drop")) {
    if (length(num_vars) == 0L)
      stop(sprintf(paste0("area_of_applicability(): categorical = \"drop\" ",
                          "leaves no predictor: every one is categorical (%s). ",
                          "Use categorical = \"indicator\"."),
                   paste(sQuote(cat_vars), collapse = ", ")), call. = FALSE)
    return(list(X_tr = .aoa_matrix(tr, num_vars, "the training data"),
                X_nw = .aoa_matrix(nw, num_vars, "`newdata`"),
                col_var = num_vars, indicator = rep(FALSE, length(num_vars)),
                cat_vars = cat_vars, vars = num_vars, unseen = none))
  }

  # "indicator": the columns keep the predictors' order.
  blocks_tr <- list(); blocks_nw <- list()
  col_var <- character(0); indicator <- logical(0)
  unseen <- matrix(FALSE, nrow = nrow(nw), ncol = length(cat_vars),
                   dimnames = list(NULL, cat_vars))
  for (v in vars) {
    if (!(v %in% cat_vars)) {
      blocks_tr[[v]] <- .aoa_matrix(tr, v, "the training data")
      blocks_nw[[v]] <- .aoa_matrix(nw, v, "`newdata`")
      col_var <- c(col_var, v); indicator <- c(indicator, FALSE)
      next
    }
    if (!.aoa_is_categorical(nw[[v]]))
      stop(sprintf(paste0("area_of_applicability(): predictor '%s' is %s in ",
                          "the training data but %s in `newdata`; its levels ",
                          "cannot be matched. Convert the column first."),
                   v, class(tr[[v]])[1L], class(nw[[v]])[1L]), call. = FALSE)
    tv <- as.character(tr[[v]])
    nv <- as.character(nw[[v]])
    # Levels that OCCUR in the training rows, not the level set: a level
    # with no training row is one the model never saw.
    lv <- sort(unique(tv[!is.na(tv)]))
    if (length(lv) == 0L)
      stop(sprintf(paste0("area_of_applicability(): categorical predictor '%s' ",
                          "is missing on every training row."), v),
           call. = FALSE)
    cn <- paste0(v, "=", lv)
    mt <- vapply(lv, function(l) as.numeric(tv == l), numeric(length(tv)))
    mn <- vapply(lv, function(l) as.numeric(nv == l), numeric(length(nv)))
    mt <- matrix(mt, nrow = length(tv), dimnames = list(NULL, cn))
    mn <- matrix(mn, nrow = length(nv), dimnames = list(NULL, cn))
    blocks_tr[[v]] <- mt
    blocks_nw[[v]] <- mn
    col_var <- c(col_var, rep(v, length(lv)))
    indicator <- c(indicator, rep(TRUE, length(lv)))
    unseen[, v] <- !is.na(nv) & !(nv %in% lv)
  }
  list(X_tr = do.call(cbind, unname(blocks_tr)),
       X_nw = do.call(cbind, unname(blocks_nw)),
       col_var = col_var, indicator = indicator, cat_vars = cat_vars,
       vars = vars, unseen = unseen)
}


#' Group the rows of a numeric matrix that are exactly equal
#'
#' Exact means \code{==} on every column: rows are sorted lexicographically
#' and each is compared with its predecessor, so nothing is rounded or
#' pasted to text on the way (\code{duplicated()} on a matrix compares 15
#' significant digits).  The first occurrence of each distinct row is its
#' representative, and with no repeated row \code{keep} is \code{1:n}.
#'
#' @return A list: \code{keep} (row numbers of the distinct rows, ascending),
#'   \code{map} (for every row, the position in \code{keep} of its
#'   representative) and \code{n_dup} (rows that repeat an earlier one).
#' @keywords internal
#' @noRd
.aoa_row_groups <- function(X) {
  n <- nrow(X)
  if (n < 2L || ncol(X) == 0L)
    return(list(keep = seq_len(n), map = seq_len(n), n_dup = 0L))
  o <- do.call(order, c(lapply(seq_len(ncol(X)), function(j) X[, j]),
                        list(method = "radix")))
  Xo <- X[o, , drop = FALSE]
  new_grp <- c(TRUE, rowSums(Xo[-1L, , drop = FALSE] !=
                               Xo[-n, , drop = FALSE]) > 0L)
  if (all(new_grp))
    return(list(keep = seq_len(n), map = seq_len(n), n_dup = 0L))
  # The radix sort is stable, so the first row of each run is the group's
  # first occurrence in the original order.
  rep_row <- o[new_grp][cumsum(new_grp)]
  first   <- integer(n)
  first[o] <- rep_row
  keep <- sort(unique(first))
  list(keep = keep, map = match(first, keep), n_dup = n - length(keep))
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
  # st_coordinates() returns one row per VERTEX, so a POLYGON or MULTIPOINT
  # layer -- this package's own create_grid_polygons() output, say -- produced
  # a longer vector than the data frame has rows and the assignment below died
  # with a bare "replacement has 45 rows, data has 9", naming nothing.  Every
  # sibling (make_folds(), estimate_sac_range(), prep_model_data()) pointizes
  # first; do the same.
  if (inherits(data_sf, "sf") &&
      !all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT"))
    data_sf <- coerce_to_points(data_sf, "auto")
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
#' The 0/1 columns of a categorical predictor (\code{categorical =
#' "indicator"}) are not scaled by their own standard deviations, which
#' would put a rare level a long way from every other row for being rare.
#' All the columns of one predictor share the scale
#' \code{sqrt(1 - sum(p^2))}, \code{p} being the levels' shares of the
#' training rows.  Two rows on different levels are then
#' \code{sqrt(2 / (1 - sum(p^2)))} apart on that predictor whichever levels
#' they are, two on the same level are 0 apart, and the mean squared
#' difference between two training rows is 2, exactly what a standardised
#' numeric predictor contributes.
#'
#' @param col_var,indicator As returned by \code{.aoa_design()}; \code{NULL}
#'   when every column is numeric.
#' @return A list with \code{center}, \code{scale} and \code{keep} (predictors
#'   with usable variance).
#' @keywords internal
#' @noRd
.aoa_scaling <- function(X_train, tol = .Machine$double.eps^0.5,
                         col_var = NULL, indicator = NULL) {
  ctr <- colMeans(X_train)
  scl <- apply(X_train, 2L, stats::sd)
  # RELATIVE tolerance.  An absolute floor on the raw standard deviation
  # (1.49e-8) drops a predictor for being recorded in the wrong UNIT rather
  # than for being constant: the same variable in metres and in gigametres has
  # the same information and the same standardised values, but only the first
  # survived -- and the second's absence flipped an obvious extrapolation
  # (DI 25.9, outside the AOA) into a point comfortably inside it (DI 0.03).
  # The manual says "predictors with zero variance are dropped"; scale the test
  # by the column's own magnitude so that is what it means.
  mag  <- apply(abs(X_train), 2L, function(z) max(z, na.rm = TRUE))
  mag[!is.finite(mag) | mag <= 0] <- 1
  keep <- is.finite(ctr) & is.finite(scl) & scl > tol * mag
  names(keep) <- colnames(X_train)
  if (!is.null(indicator) && any(indicator)) {
    for (v in unique(col_var[indicator])) {
      cols <- which(indicator & col_var == v)
      gini <- 1 - sum(ctr[cols]^2)
      # A predictor with one level among the training rows has a constant
      # column, which the variance test above has already dropped.
      if (is.finite(gini) && gini > 0) scl[cols] <- sqrt(gini)
    }
  }
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
#' The dissimilarity index is invariant to the overall scale of the weights:
#' both the nearest-neighbour distance and the normalising mean pairwise
#' distance carry the same factor.  Users need not normalise importance values
#' before passing them.  They are rescaled to mean 1 here purely so the
#' recorded values are readable.
#'
#' @param fill_vars Character vector of predictors that this package appended
#'   to \code{vars} itself (the \code{"..x"}/\code{"..y"} coordinate columns of
#'   a coordinate-using model).  The caller has never seen them, so a missing
#'   weight for one is filled in without an error, and unnamed weights may be
#'   one-per-\emph{user-visible}-predictor.
#' @keywords internal
#' @noRd
.aoa_weight_vector <- function(weights, vars, fill_vars = character(0)) {
  if (is.null(weights)) return(stats::setNames(rep(1, length(vars)), vars))
  if (!is.numeric(weights))
    stop("area_of_applicability(): `weights` must be numeric, or NULL for ",
         "equal weights, or \"auto\".", call. = FALSE)
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
  # A non-finite weight got the negative-weight advice, pmax(importance, 0),
  # which keeps NaN: a forest with no out-of-bag rows has NaN permutation
  # importance, and following the hint reproduced the same error.
  bad <- !is.finite(w)
  if (any(bad))
    stop(sprintf("area_of_applicability(): `weights` is %s for %s.%s",
                 if (all(is.nan(w[bad]))) "NaN" else "not finite",
                 paste(sQuote(names(w)[bad]), collapse = ", "),
                 if (any(is.nan(w[bad])))
                   paste0(" Permutation importance is NaN when no row is out ",
                          "of bag (a forest grown with replace = FALSE and ",
                          "sample_fraction = 1), and pmax() keeps NaN; refit ",
                          "with out-of-bag rows (replace = TRUE or ",
                          "sample_fraction < 1), or pass weights = NULL.")
                 else ""), call. = FALSE)
  if (any(w < 0))
    stop("area_of_applicability(): `weights` must be finite and non-negative. ",
         "Permutation importance is slightly negative for predictors that do ",
         "not help, so pass pmax(importance, 0).", call. = FALSE)
  w <- .aoa_equal_if_all_zero(w)
  w * (length(w) / sum(w))
}

#' Replace all-zero weights with equal ones
#'
#' The advice above, \code{pmax(importance, 0)}, is all zero whenever the
#' model found no predictor useful, which a forest with one predictor does
#' often (9 fits in 20 when that predictor carried no signal).  That used to
#' be an error, and it lost the AOA in exactly the folds where the predictor
#' was weakest.  With ONE predictor the weight cannot matter: the index is
#' invariant to the scale of the weights, so zero is just a degenerate way of
#' writing any positive value, and it is accepted silently.  With several,
#' all-zero weights say nothing about how the predictors compare, so they are
#' weighted equally, as \code{weights = NULL} would, with a warning.
#'
#' @keywords internal
#' @noRd
.aoa_equal_if_all_zero <- function(w) {
  if (length(w) == 0L || any(w != 0)) return(w)
  if (length(w) > 1L)
    .warn_and_log(paste0("area_of_applicability(): every weight is zero (%s), ",
                         "which says nothing about how the predictors compare; ",
                         "weighting them equally, as weights = NULL does. ",
                         "pmax(importance, 0) is all zero when the model found ",
                         "no predictor useful."),
                  paste(names(w), collapse = ", "))
  w[] <- 1
  w
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
  # regardless of grid size.  ||a-b||^2 = |a|^2 + |b|^2 - 2a.b is used for
  # the SEARCH only.  It loses precision for near-coincident points: for an
  # exact copy of a reference row the two terms differ by a few ulp, and the
  # result is either clamped to 0 or comes out near 1e-8.  That is negligible
  # against an ordinary threshold, but not against a threshold at or near 0,
  # where it put some exact copies of training rows outside the AOA and kept
  # a threshold that should have been 0 just above it.  So the distance
  # returned is recomputed exactly, as sqrt(sum((q - x)^2)), for every
  # reference row the search cannot tell from the nearest.
  #
  # Repeated reference rows cannot change a minimum, and a block of them is
  # a block of ties for the exact pass, so they are collapsed first.  (Not
  # under `exclude`, which indexes the rows as given; the leave-one-out
  # caller passes distinct rows.)
  if (is.null(exclude)) {
    g <- .aoa_row_groups(data)
    if (g$n_dup > 0L) {
      data <- data[g$keep, , drop = FALSE]
      nd   <- nrow(data)
    }
  }
  if (is.null(chunk_size))
    chunk_size <- max(1L, min(10000L, as.integer(floor(4e6 / max(1L, nd)))))
  # A whole number of rows per block.  A fractional chunk_size gave fractional
  # block starts, the s:e ranges below truncated, and the rows between blocks
  # kept the 0 they were initialised with: a DI of 0, inside the AOA.
  chunk_size <- max(1L, as.integer(chunk_size))
  out  <- numeric(nq)
  dsq  <- rowSums(data^2)
  dmax <- max(dsq)
  # The rounding error of the fast formula for one pair is a small multiple
  # of eps * (|q|^2 + |x|^2), growing with the number of columns summed.
  slack <- 8 * (ncol(data) + 8) * .Machine$double.eps
  for (s in seq.int(1L, nq, by = chunk_size)) {
    e   <- min(s + chunk_size - 1L, nq)
    Q   <- query[s:e, , drop = FALSE]
    qsq <- rowSums(Q^2)
    d2  <- outer(qsq, dsq, "+") - 2 * tcrossprod(Q, data)
    if (!is.null(exclude)) {
      ex <- exclude[s:e]
      ok <- !is.na(ex)
      if (any(ok)) d2[cbind(which(ok), as.integer(ex[ok]))] <- Inf
    }
    out[s:e] <- .aoa_exact_row_min(d2, Q, data, slack * (qsq + dmax))
  }
  out
}


#' Row minima of a fast squared-distance block, recomputed exactly
#'
#' @param d2 Query-by-reference squared distances from the fast formula
#'   (\code{Inf} where a pair is excluded).
#' @param Q,data The query and reference rows \code{d2} was built from.
#' @param tol One tolerance per query row: the rounding error \code{d2} can
#'   carry.  Every reference row within it of the row minimum is a candidate
#'   for nearest, and the smallest exact distance among them is returned.
#' @return One distance per query row; \code{Inf} where every pair was
#'   excluded.
#' @keywords internal
#' @noRd
.aoa_exact_row_min <- function(d2, Q, data, tol) {
  nr <- nrow(d2)
  # An exact row minimum.  max.col() is not one: it calls entries within
  # 1e-5 of the row's largest magnitude tied.
  m <- if (ncol(d2) <= nr) {
    acc <- d2[, 1L]
    for (j in seq_len(ncol(d2))[-1L]) acc <- pmin(acc, d2[, j])
    acc
  } else apply(d2, 1L, min)
  res <- rep(Inf, nr)
  res[is.nan(m)] <- NaN
  lim <- m + tol
  lim[!is.finite(m)] <- -Inf
  cand <- which(d2 <= lim, arr.ind = TRUE)
  if (nrow(cand) == 0L) return(res)
  # Bounded batches: a query row among many near-coincident reference rows
  # has as many candidates.
  step <- 5e5L
  for (b in seq.int(1L, nrow(cand), by = step)) {
    idx <- b:min(b + step - 1L, nrow(cand))
    qi  <- cand[idx, 1L]
    ex  <- sqrt(rowSums((Q[qi, , drop = FALSE] -
                           data[cand[idx, 2L], , drop = FALSE])^2))
    o   <- order(qi, ex)
    top <- !duplicated(qi[o])
    hit <- qi[o][top]
    res[hit] <- pmin(res[hit], ex[o][top])
  }
  res
}


#' Mean pairwise distance among training points (the DI normaliser)
#'
#' @keywords internal
#' @noRd
.aoa_normalizer <- function(X, max_n = 5000L, seed = 123L,
                            what = c("predictors", "location")) {
  what <- match.arg(what)
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
    .log_info("area_of_applicability(): %snormaliser estimated from a %d-point subsample of %d training points.",
              if (identical(what, "location")) "location " else "",
              as.integer(max_n), n)
  }
  m <- mean(as.numeric(stats::dist(X)))
  if (!is.finite(m) || m <= 0) {
    if (identical(what, "location"))
      stop("area_of_applicability(): every training row is at the same ",
           "location, so the mean distance between training points is zero ",
           "and the location index is undefined. Pass coords = \"no\" to ",
           "measure the predictors alone.", call. = FALSE)
    stop("area_of_applicability(): the mean pairwise training distance is ",
         "zero or undefined; the training predictors carry no spread.",
         call. = FALSE)
  }
  list(value = m, subsampled = subsampled, n_used = nrow(X))
}


#' Map coordinates of both layers, for the location index
#'
#' Works on copies: the caller's \code{newdata} is returned to it unchanged.
#' The prediction layer is brought into the training layer's CRS (or, when
#' only it has one, the training layer into its), exactly as on the
#' coordinates-as-predictors path, and non-\code{POINT} geometry is reduced
#' to representative points.  A geographic CRS is refused: the index is a
#' map distance, and a degree of longitude is not a fixed one.
#'
#' @param why A clause naming why location is measured, for the messages.
#' @return A list of two two-column matrices, \code{train} and \code{new},
#'   row for row with the layers; a row with an empty geometry is \code{NA}.
#' @keywords internal
#' @noRd
.aoa_location_coords <- function(train_sf, newdata, why) {
  if (!inherits(train_sf, "sf") || !inherits(newdata, "sf"))
    stop(sprintf(paste0("area_of_applicability(): %s, so location is measured ",
                        "and both the training data and `newdata` must be sf ",
                        "objects carrying geometry (%s). Pass coords = \"no\" ",
                        "to measure the predictors alone."),
                 why,
                 paste(c(if (!inherits(train_sf, "sf")) "the training data is not",
                         if (!inherits(newdata, "sf")) "`newdata` is not"),
                       collapse = "; ")), call. = FALSE)
  tr_crs <- sf::st_crs(train_sf)
  if (!is.na(tr_crs)) {
    if (!identical(sf::st_crs(newdata), tr_crs))
      newdata <- .transform_or_stamp(newdata, tr_crs, what = "newdata",
                                     caller = "area_of_applicability")
  } else if (!is.na(sf::st_crs(newdata))) {
    train_sf <- .transform_or_stamp(train_sf, sf::st_crs(newdata),
                                    what = "the training data",
                                    caller = "area_of_applicability")
  }
  if (isTRUE(sf::st_is_longlat(train_sf)))
    stop("area_of_applicability(): the location index is a map distance and ",
         "the training data is in a geographic CRS (longitude/latitude). ",
         "Project both layers first, with sf::st_transform() or ",
         "ensure_projected(), or pass coords = \"no\".", call. = FALSE)
  xy <- function(x) {
    if (!all(sf::st_geometry_type(x, by_geometry = TRUE) == "POINT"))
      x <- coerce_to_points(x, "auto")
    crd <- sf::st_coordinates(x)
    if (ncol(crd) < 2L || nrow(crd) != nrow(x))
      stop("area_of_applicability(): could not extract one X/Y coordinate ",
           "pair per row.", call. = FALSE)
    m <- cbind(as.numeric(crd[, 1L]), as.numeric(crd[, 2L]))
    colnames(m) <- c("x", "y")
    m
  }
  list(train = xy(train_sf), new = xy(newdata))
}


#' Raise the zero-threshold warning with a condition class
#'
#' \code{.warn_and_log()} with the class
#' \code{spatialkit_aoa_zero_threshold}, so that a script computing many
#' areas of applicability, or one inside a cross-validation loop, can catch
#' exactly this case with \code{tryCatch(spatialkit_aoa_zero_threshold = )}
#' or \code{withCallingHandlers()} without matching message text.
#'
#' @keywords internal
#' @noRd
.aoa_warn_zero_threshold <- function(fmt, ...) {
  .warn_classed("spatialkit_aoa_zero_threshold", fmt, ...)
}


#' The weights `weights = "auto"` resolves to
#'
#' Meyer and Pebesma (2021) weight the predictors by the model's variable
#' importance, so that is what a model carrying one gets:
#' \code{pmax(importance, 0)}.  Equal weights are the answer whenever the
#' importance is absent or does not mean anything, and each such case is
#' logged, not raised:
#' \itemize{
#'   \item no \code{model}, or one with no importance (\code{importance =
#'     "none"}, a GWR, a Bayesian fit, a custom learner);
#'   \item an importance that is not fully finite (it is \code{NaN} for a
#'     forest with no out-of-bag rows);
#'   \item an importance that does not name every predictor in use
#'     (\code{predictor_vars} overridden with a column the model never saw);
#'   \item no positive importance at all;
#'   \item a forest with no out-of-bag skill (\code{oob_r_squared} missing
#'     or at most 0): on a pure-noise response 2 to 7 of 10 importances are
#'     positive by chance and the largest takes 26 to 96\% of the weight.
#' }
#'
#' @param vars The predictors in use, including any in \code{fill_vars}.
#' @param fill_vars The coordinate columns of a coordinate-using forest;
#'   its importance names them, so they are weighted like the others.
#' @return A list: \code{weights} (a named numeric, or \code{NULL} for
#'   equal), \code{source} (\code{"model importance"} or \code{"equal"})
#'   and \code{note}.
#' @keywords internal
#' @noRd
.aoa_auto_weights <- function(model, vars, fill_vars = character(0)) {
  equal <- function(why, log = TRUE) {
    if (log)
      .log_info(paste0("area_of_applicability(): weights = \"auto\": %s, so ",
                       "every predictor counts equally."), why)
    list(weights = NULL, source = "equal", note = why)
  }
  if (is.null(model)) return(equal("no model was given", log = FALSE))
  imp <- model$info$importance
  if (is.null(imp) || length(imp) == 0L)
    return(equal("the model carries no variable importance"))
  if (!is.numeric(imp) || is.null(names(imp)))
    return(equal("the model's importance is not a named numeric vector"))
  need <- intersect(vars, names(imp))
  miss <- setdiff(setdiff(vars, fill_vars), names(imp))
  if (length(miss) > 0L)
    return(equal(sprintf("the model's importance has no entry for %s",
                         paste(sQuote(miss), collapse = ", "))))
  if (any(!is.finite(imp[need])))
    return(equal(paste0("the model's importance is not finite (it is NaN for ",
                        "a forest with no out-of-bag rows)")))
  if (!any(imp[need] > 0))
    return(equal("no predictor has a positive importance"))
  r2 <- .num1(model$info$oob_r_squared)
  if (!is.finite(r2))
    return(equal(paste0("the model records no out-of-bag R-squared to show ",
                        "its importances are more than noise")))
  if (r2 <= 0)
    return(equal(sprintf(paste0("the model has no out-of-bag skill (OOB ",
                                "R-squared %.3f), so its importances are noise"),
                         r2)))
  type <- model$info$importance_type
  type <- if (is.character(type) && length(type) == 1L && !is.na(type) &&
              nzchar(type)) type else "variable"
  list(weights = pmax(imp[need], 0), source = "model importance",
       note = sprintf("%s importance", type))
}


#' Normalise a fold specification into train/test position lists
#'
#' Accepts a \code{make_folds()} result, a bare list of \code{train}/\code{test}
#' splits, or a plain vector of fold labels.  The \code{train} slot is used as
#' given and never reconstructed as "everything else", so buffered and NNDM
#' folds contribute the training set that was actually available, which is the
#' point of having built them.
#'
#' @param row_ids Optional vector of \code{..row_id} values, one per training
#'   row, in training-row order.  When supplied, the \code{train}/\code{test}
#'   entries of a list-shaped \code{folds} are treated as \code{..row_id}
#'   VALUES and resolved to row positions with \code{match()}; when
#'   \code{NULL} they are treated as positions, as before.  Fold labels are
#'   always positional: they are one label per row by definition.
#' @param removed_ids Optional \code{..row_id} values of rows that
#'   \code{prep_model_data()} removed from the training data (read from its
#'   \code{"dropped"} record).  They are taken out of every train and test
#'   slot before the IDs are resolved, as \code{.remap_folds()} does for
#'   \code{cv_*()}, and a fold left with no test row is dropped.  Any other
#'   unknown ID is still an error.
#' @keywords internal
#' @noRd
.aoa_fold_splits <- function(folds, n, row_ids = NULL, removed_ids = NULL) {
  if (is.null(folds)) return(NULL)

  if (is.atomic(folds) && !is.list(folds)) {
    if (length(folds) != n)
      stop(sprintf(paste0("area_of_applicability(): `folds` has %d labels but ",
                          "the training data has %d rows."),
                   length(folds), n), call. = FALSE)
    # Labels are grouped and numbered by .folds_from_labels(), the helper
    # every cv_*() uses, so fold k here is fold k in a cv_*() result built
    # from the same labels: numbers in numeric order, a factor by its own
    # used levels, anything else in C (radix) order.  A copy of that rule
    # here drifted: two doubles that print alike (0.3 and 0.1 + 0.2) stopped
    # with R's "factor level [2] is duplicated" where cv_*() made them one
    # fold.  The helper also refuses missing labels and fewer than two
    # non-empty folds.  With ..row_id = 1..n its IDs are row positions, so
    # there is nothing to resolve.  Fall through to the shared validation
    # below, so label-built splits get the same checks as hand-built ones.
    sp <- .folds_from_labels(folds, data.frame(..row_id = seq_len(n)),
                             "area_of_applicability")
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

  # Folds built on the layer a model was fitted from name every row of it,
  # including those prep_model_data() then removed (a missing response, a
  # non-finite predictor, an empty geometry).  cv_*() drops such IDs from the
  # splits (.remap_folds()); here they were resolved against the shorter
  # training data and every one of them stopped the call with "refers to rows
  # outside 1:n" -- for the same folds cv_*() had just accepted.  Drop them,
  # and say how many.
  if (length(removed_ids) && !is.null(row_ids)) {
    named <- unique(unlist(lapply(sp, function(z) c(z$train, z$test)),
                           use.names = FALSE))
    n_hit <- sum(named %in% removed_ids)
    if (n_hit > 0L) {
      sp <- lapply(sp, function(z) {
        z$train <- z$train[!(z$train %in% removed_ids)]
        z$test  <- z$test[!(z$test %in% removed_ids)]
        z
      })
      no_test <- vapply(sp, function(z) length(z$test) == 0L, logical(1))
      .log_info(paste0("area_of_applicability(): `folds` name %d row(s) that ",
                       "prep_model_data() removed from the training data ",
                       "(missing or non-finite values, or an empty geometry); ",
                       "they were dropped from every fold%s."),
                n_hit,
                if (any(no_test))
                  sprintf(", and the %d fold(s) left with no test row with them",
                          sum(no_test))
                else "")
      sp <- sp[!no_test]
      if (length(sp) == 0L)
        stop("area_of_applicability(): no fold has a test row left once the ",
             "rows prep_model_data() removed are taken out of `folds`.",
             call. = FALSE)
    }
  }

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
                          "set, or on rows that have since been removed other ",
                          "than by prep_model_data(), whose removals are ",
                          "dropped from the folds while the training data ",
                          "carries its \"dropped\" record."),
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
                          "..row_id. (Rows prep_model_data() removed are ",
                          "dropped from the folds when the training data ",
                          "carries its \"dropped\" record, as a fit's data_sf ",
                          "does; a subset or a re-ordering of it does not.)"),
                   j, n), call. = FALSE)
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


#' Line the folds up with a training layer prep_model_data() removed rows from
#'
#' A fit's \code{data_sf} is the \code{prep_model_data()} output: rows with a
#' missing or non-finite modelling value or an empty geometry are gone, and
#' it carries no \code{..row_id} unless its input did.  Folds built on the
#' layer the model was fitted FROM -- the ones passed to \code{cv_*()} -- name
#' the input's rows.  The \code{"dropped"} record says which those were, so
#' the input's row IDs of the kept rows can be rebuilt: positions when the
#' input had no \code{..row_id} (make_folds() numbers the rows then), the
#' recorded IDs otherwise.
#'
#' When the training data carries \code{..row_id}, list-shaped folds are
#' resolved by ID already, and only the removed rows' IDs are needed.
#' Without it, the IDs of folds built on the training data itself are its
#' row positions, \code{1:n}; folds built on the layer fitted from name that
#' layer's rows, and so name a row past \code{n} whenever a row before its
#' end was removed (a trailing row's removal leaves the numbering the same).
#' So a fold list naming a row past \code{n} is read in the input's
#' numbering, and one that does not is left as positions.  A label vector
#' with one label per input row is cut down to the kept rows; one with one
#' label per training row is left alone.
#'
#' @return A list: \code{folds}, \code{row_ids} (the training rows' IDs, or
#'   the \code{row_ids} given) and \code{removed_ids} (the IDs of the removed
#'   rows, or \code{NULL}).
#' @keywords internal
#' @noRd
.aoa_rows_for_folds <- function(folds, train_sf, row_ids) {
  out <- list(folds = folds, row_ids = row_ids, removed_ids = NULL)
  rec <- if (is.data.frame(train_sf)) .get_row_record(train_sf, "dropped") else NULL
  n_rm <- suppressWarnings(as.integer(rec$n %||% 0L))
  if (is.null(rec) || length(n_rm) != 1L || is.na(n_rm) || n_rm < 1L)
    return(out)
  n    <- as.integer(nrow(train_sf))
  n_in <- n + n_rm
  kept <- setdiff(seq_len(n_in), as.integer(rec$which))
  # A record that does not add up describes other rows; ignore it.
  if (length(kept) != n) return(out)
  if (is.atomic(folds) && !is.list(folds)) {
    if (length(folds) == n_in) {
      .log_info(paste0("area_of_applicability(): `folds` has one label per row ",
                       "of the layer the model was fitted from; the %d label(s) ",
                       "of the row(s) prep_model_data() removed were dropped."),
                n_rm)
      out$folds <- folds[kept]
    }
    return(out)
  }
  if (!is.null(row_ids)) {
    if (length(rec$row_id)) out$removed_ids <- rec$row_id
    return(out)
  }
  sp <- if (is.list(folds) && is.list(folds$folds)) folds$folds else folds
  named <- if (is.list(sp))
    suppressWarnings(as.numeric(unlist(lapply(sp, function(z)
      if (is.list(z)) c(z$train, z$test)), use.names = FALSE)))
  else numeric(0)
  if (length(named) && any(named > n, na.rm = TRUE)) {
    out$row_ids     <- kept
    out$removed_ids <- as.integer(rec$which)
  }
  out
}


#' The fold provenance check of cv_*(), for area_of_applicability()
#'
#' \code{.check_fold_probe()} on the training data, with its row IDs as the
#' folds are resolved against them.  A fit's \code{data_sf} holds the points
#' \code{prep_model_data()} reduced its layer to, so folds built on the
#' polygons it was fitted from carry a probe taken on polygon centroids: the
#' two cannot be compared location by location, and the check is skipped
#' (logged) rather than refusing a documented workflow.
#' @keywords internal
#' @noRd
.aoa_check_fold_probe <- function(folds, train_sf, row_ids) {
  probe <- tryCatch(folds$params$row_probe, error = function(e) NULL)
  if (!is.list(probe) || !inherits(train_sf, "sf") || nrow(train_sf) == 0L)
    return(invisible(NULL))
  tr_points <- tryCatch(
    all(sf::st_geometry_type(train_sf, by_geometry = TRUE) == "POINT"),
    error = function(e) NA)
  if (is.logical(probe$points) && length(probe$points) == 1L &&
      !is.na(probe$points) && !identical(probe$points, tr_points)) {
    .log_info(paste0("area_of_applicability(): `folds` were built on %s ",
                     "geometry and the training data has %s geometry, so their ",
                     "row locations cannot be compared; skipping the ",
                     "provenance check."),
              if (probe$points) "POINT" else "non-POINT",
              if (isTRUE(tr_points)) "POINT" else "non-POINT")
    return(invisible(NULL))
  }
  chk <- train_sf
  chk[["..row_id"]] <- if (is.null(row_ids)) seq_len(nrow(train_sf)) else row_ids
  .check_fold_probe(folds, chk, "area_of_applicability")
}


#' Reference nearest-neighbour distances for the training data
#'
#' Without folds, each row's reference is its nearest row at a DIFFERENT
#' point of the space.  A row that repeats another's values exactly (a
#' repeat visit to a site with static covariates, covariates read off a
#' raster coarser than the sampling) is the same point, not a neighbour:
#' counted as one it has a reference distance of 0, and once about three
#' quarters of the rows have such a twin the threshold is 0 and only exact
#' copies of a training row are inside.  So the search runs on the distinct
#' rows, and every row takes its point's distance.  With no repeated row
#' that is the search on the rows as given.
#'
#' With folds the reference is the nearest row the fold left in training,
#' as defined: a twin on the other side of a split is a row the hold-out
#' really did have in its training set.
#'
#' @param groups Optional \code{.aoa_row_groups(X)}, when the caller has it.
#' @keywords internal
#' @noRd
.aoa_train_dist <- function(X, splits = NULL,
                            use_fnn = requireNamespace("FNN", quietly = TRUE),
                            chunk_size = NULL, groups = NULL) {
  n <- nrow(X)
  if (n < 2L)
    stop("area_of_applicability(): at least two training rows are needed.",
         call. = FALSE)

  if (is.null(splits)) {
    if (is.null(groups)) groups <- .aoa_row_groups(X)
    Xd <- if (groups$n_dup > 0L) X[groups$keep, , drop = FALSE] else X
    if (nrow(Xd) < 2L)
      stop("area_of_applicability(): every training row has the same values ",
           "on the predictors measured, so no row has a neighbour to take a ",
           "reference distance from.", call. = FALSE)
    d <- if (isTRUE(use_fnn))
      as.numeric(FNN::get.knn(Xd, k = 1L)$nn.dist[, 1L])
    else .aoa_min_dist(Xd, Xd, use_fnn = FALSE, chunk_size = chunk_size,
                       exclude = seq_len(nrow(Xd)))
    return(if (groups$n_dup > 0L) d[groups$map] else d)
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
#' usual rule, i.e. the largest value at or below \code{Q3 + 1.5 * IQR}, with
#' \code{stats::quantile()}'s default (type 7) quartiles.  That is the
#' "(outlier-removed) maximum" of Meyer & Pebesma (2021).  CAST, the reference
#' implementation, uses the same type-7 fence since it stopped calling
#' \code{grDevices::boxplot.stats()} (whose upper whisker is this rule on
#' Tukey's hinges), but takes the FENCE itself, capped at the largest training
#' DI, as the threshold.  The two agree whenever nothing lies above the fence;
#' otherwise CAST's threshold is the larger.  The paper's definition is kept,
#' and the help page says how to reproduce CAST's value.
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

#' How many training DI values the threshold's fence set aside
#'
#' The same \code{Q3 + 1.5 * IQR} fence as \code{.aoa_threshold()}; the count
#' above it is what the "outlier-removed" in the threshold's name refers to,
#' and it used to be discarded.
#'
#' @keywords internal
#' @noRd
.aoa_fence_outliers <- function(di) {
  d <- di[is.finite(di)]
  if (length(d) == 0L) return(NA_integer_)
  fence <- stats::quantile(d, 0.75, names = FALSE) + 1.5 * stats::IQR(d)
  as.integer(sum(d > fence))
}


# -----------------------------------------------------------------------------
# Exported
# -----------------------------------------------------------------------------

#' Area of applicability of a spatial prediction model
#'
#' Computes the dissimilarity index (DI) of Meyer & Pebesma (2021) for a set of
#' prediction locations and flags those that fall inside the model's area of
#' applicability (AOA), the region of predictor space where the model's
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
#' standard deviations, then weighted by variable importance when the model
#' carries one (see \emph{Weights}). For a
#' prediction point \eqn{p}, the DI is the distance to its nearest training
#' point in that space, divided by the mean pairwise distance among training
#' points. The same quantity is computed for the training data itself, using
#' each point's nearest neighbour \emph{among the training rows of the fold
#' that holds it out}. That means everything outside its own fold for random
#' and block folds, and the smaller training set that buffered and NNDM folds
#' actually leave (see \emph{The fold scheme changes the answer}); without
#' folds it is the nearest training row with different predictor values. The
#' threshold is then the largest
#' training DI that is not an upper outlier, i.e. not above the fence
#' \code{Q3 + 1.5 * IQR} of the training DI, with the quartiles of
#' \code{stats::quantile()}'s default type 7. Prediction points at or below
#' that threshold are inside the AOA.
#'
#' That is the paper's "outlier-removed maximum". \pkg{CAST}, the reference
#' implementation, computes the same fence with the same quartiles but uses the
#' fence itself as the threshold, capped at the largest training DI. The two
#' agree whenever no training DI lies above the fence (\code{n_outliers} is 0);
#' otherwise \pkg{CAST}'s threshold is the larger, and so is its AOA. (Earlier
#' \pkg{CAST} releases used \code{grDevices::boxplot.stats()}, which gives the
#' rule used here but with Tukey's hinges as the quartiles, so they can also
#' differ when the number of training points is even.) To apply the current
#' \pkg{CAST} rule to the same training DI, pass
#' \code{threshold = min(quantile(res$train_DI, 0.75) + 1.5 * IQR(res$train_DI),
#' max(res$train_DI))} for an earlier result \code{res}.
#'
#' The DI is invariant to the overall scale of \code{weights}: the numerator
#' and the normaliser carry the same factor. Importance values can be passed
#' as-is.
#'
#' @section Training rows with identical predictor values:
#' Repeat visits to a site with static covariates, or covariates read off a
#' raster coarser than the sampling, give training rows that are exact copies
#' of each other in predictor space. A copy is the same point of that space,
#' not a neighbour. Counted as a neighbour it gives a training DI of 0, and
#' once about three quarters of the rows have one the threshold is 0 and only
#' exact copies of a training row count as inside (60 sites visited 4 times
#' put 30 of 200 new points inside; the same sites counted once put 188).
#'
#' With \code{folds = NULL} each training row's reference is therefore its
#' nearest row with \emph{different} values: the search runs on the distinct
#' rows and every row takes its point's distance. \code{n_duplicates} counts
#' the rows that repeat an earlier one, and a layer with none is computed
#' exactly as before. The normaliser stays the mean pairwise distance over
#' all the training rows, so a prediction point's DI is the same with and
#' without folds.
#'
#' With \code{folds} supplied the definition is applied as it stands: a
#' held-out row whose copy is on the training side of its fold has a DI of 0,
#' because that hold-out did have an identical row to learn from. If the
#' threshold comes out 0 that way, a warning of class
#' \code{spatialkit_aoa_zero_threshold} gives the count, so a script or a
#' cross-validation loop can catch it with
#' \code{tryCatch(..., spatialkit_aoa_zero_threshold = )}. Folds that keep
#' the copies together (\code{make_folds(method = "leave_location_out",
#' group_var = ...)}) give the threshold its meaning back. The same applies
#' to the location index when training rows share a location.
#'
#' @section The fold scheme changes the answer, and should:
#' With \code{folds = NULL} the training reference is each point's nearest
#' neighbour anywhere in the training data, which for clustered data is very
#' close, giving a small threshold and a conservative AOA. Passing the folds
#' you actually validated with makes the reference distances larger and the AOA
#' correspondingly wider. That is not a loophole. The AOA is defined relative
#' to a performance estimate, and a spatially blocked estimate is a claim about
#' predicting further away. Pass the same \code{make_folds()} result you passed
#' to \code{\link{cv_spatial}}. Buffered and NNDM folds use the training set
#' they actually left available, not merely "everything outside the fold".
#'
#' @section Weights:
#' Meyer and Pebesma weight the predictors by the model's variable
#' importance, and with the default \code{weights = "auto"} so does this
#' function whenever \code{model} carries an importance that means
#' something: \code{pmax(model$info$importance, 0)}, which a forest from
#' \code{\link{fit_rf_model}} has. Equal weights are used, with a logged
#' note, when it does not: no \code{model}; a model with no importance
#' (\code{importance = "none"}, a GWR, a Bayesian fit); an importance that is
#' not finite (\code{NaN} for a forest with no out-of-bag rows); one that
#' names fewer predictors than are in use; no positive value; or a forest
#' with no out-of-bag skill (\code{oob_r_squared} missing or at most 0),
#' whose positive importances are noise. \code{print()} says which it was.
#' \code{weights = NULL} asks for equal weights whatever the model carries.
#'
#' What the weights buy: with 2 predictors the forest uses and 8 it does not,
#' equal weights put about 70\% of rows shifted 3 SD along a used predictor
#' (4 times the in-sample RMSE) inside the AOA and importance weights about
#' 26\%; of rows shifted 4 SD along an unused predictor (error unchanged),
#' 43\% against 93\%. What they cost: a weight multiplies the scaled column,
#' so its square enters the distance, and a weak but real predictor is
#' down-weighted by the square of its relative importance. A shift along one
#' can then go unseen (94\% inside at 1.9 times the in-sample RMSE in one
#' test, against 19\% with equal weights).
#' \code{weights = sqrt(pmax(importance, 0))} is the gentler alternative: it
#' kept most of the gain on used predictors and put 46\% of that weak-predictor
#' shift inside, at the cost of 77\% in place of 93\% inside for harmless
#' shifts.
#'
#' @section Limitations:
#' Categorical predictors (factor or character columns) have no distance of
#' their own, and by default they are refused and never silently
#' dummy-coded. \code{categorical = "indicator"} codes each level that
#' occurs in the training rows as a 0/1 column. The columns of one predictor
#' share one scale, \code{sqrt(1 - sum(p^2))} with \code{p} the levels'
#' shares of the training rows, so two rows on different levels are the same
#' distance apart whichever levels they are, two on the same level are 0
#' apart, and the predictor adds as much to the mean squared distance between
#' training rows as one standardised numeric predictor. That scaling is a
#' choice, not a result, which is why it is not the default. A prediction row
#' on a level no training row has gets \code{DI = Inf}, is outside the AOA,
#' and a warning gives the count. \code{categorical = "drop"} leaves the
#' categorical predictors out and computes the index on the others; a level
#' the training data never has then goes unseen.
#'
#' Predictors whose variance is negligible \emph{relative
#' to their own magnitude} (the test is
#' \code{sd < sqrt(.Machine$double.eps) * max(abs(x))}, so the same variable in
#' metres and in gigametres is treated identically) are dropped from the
#' distance. A prediction point taking a value there outside the training
#' range (with the same relative tolerance) is extrapolation along a direction
#' the training data never varied in: its scaled distance along it is
#' infinite, so it gets \code{DI = Inf}, is outside the AOA, and a warning
#' gives the count. A point missing that value is judged on the other
#' predictors. With equal weights every predictor counts the same, which
#' overstates dissimilarity along directions the model barely uses; see
#' \emph{Weights}.
#'
#' @section Models whose predictions depend on location:
#' For three kinds of model the prediction is a function of where the point
#' is, not only of its covariates, so the index has to measure location too.
#' Without that a prediction point far outside the training extent but with
#' ordinary covariate values reads as \emph{inside} the area of
#' applicability, which is exactly the extrapolation this index exists to
#' catch. \code{coords = "auto"} (the default) measures location for all
#' three; \code{coords = "no"} turns it off, with a logged caution;
#' \code{coords = "yes"} measures it for any model, or for none.
#'
#' \strong{A geographically weighted regression} (\code{gwr_fit}) has
#' coefficients that vary over the map, and \strong{a Bayesian spatial model}
#' (\code{bayesian_fit}) has a Gaussian-process term in the coordinates. For
#' these, and under \code{coords = "yes"}, location gets an index of its own,
#' built like the predictor index but in map space: \code{DI_location} is the
#' distance to the nearest training point over the mean distance between
#' training points; the training points' own values use the same folds; and
#' the threshold is the same outlier-removed maximum. A prediction point is
#' inside the AOA only if it is inside on both indices, and both DI columns
#' are returned so the cause can be read off. Location is kept apart from the
#' covariates because, as two more equally weighted columns, it is diluted as
#' covariates are added: with 12 of them, 89\% of points 200 to 500 units
#' beyond a 1000-unit training extent still read as inside, against none with
#' the separate index. In a test with a GWR whose slope rises eastward, points
#' one half to one and a half extents beyond the data had 13 times the
#' in-extent RMSE; 84\% of them were inside on the covariates alone and none
#' with location measured.
#'
#' Inside the extent the index marks an unsampled gap only where the gap is
#' wider than the distance its threshold allows, and that distance depends
#' on \code{folds}. With \code{folds = NULL} it is the spacing of the
#' training points (the outlier-removed largest nearest-neighbour distance).
#' With blocked folds it is the largest distance from a held-out point to
#' the nearest row its fold trained on, outliers removed, which grows with
#' the block size; a gap narrower than about twice that distance has no
#' cell beyond it and reads as inside. In a test with 400 training points
#' in two strips 400 units apart, the cut-off without folds was 46 to 54
#' units and no grid cell in the middle half of the gap, 100 units or more
#' from either strip, was inside. That held although the prediction grid
#' took its covariates from the nearest observation, so that every cell's
#' covariate DI was 0. With five block folds the cut-off was 199 to 251
#' units and 76 to 100\% of the same cells were inside (87\% on average over
#' three draws). This follows from the definition and is not a fault of the
#' index: the cross-validated error it is paired with was earned on
#' hold-outs that far from their training rows. It does mean that with
#' blocked folds the location index says nothing about a gap narrower than
#' about twice the cut-off; \code{location$threshold * location$normalizer}
#' is the cut-off in map units, and \code{print()} shows it.
#'
#' \strong{A forest fitted with \code{include_coords = TRUE}} splits on the
#' coordinates like on any other column, so there they join the predictor
#' index as \code{"..x"} and \code{"..y"} and are centred, scaled and weighted
#' with the rest.
#'
#' Either way geometry is needed on both sides, so \code{train_sf} and
#' \code{newdata} must both be \code{sf} objects; a data.frame is refused
#' rather than quietly measured without location. Non-\code{POINT}
#' \code{newdata} (grid polygons, say) is reduced to representative points
#' first, as \code{\link{coerce_to_points}} would. If exactly one side carries
#' a CRS the other is brought into it: reprojected when its coordinates look
#' like longitude/latitude, stamped otherwise, with a warning either way. This
#' is done because degrees fed into a metre-space index silently understate the
#' distances. The location index is a map distance and refuses training data
#' in a geographic CRS.
#'
#' @param newdata Prediction locations: an \code{sf} object (typically from
#'   \code{\link{predict_surface}}) or a data.frame, carrying the predictor
#'   columns. When location is measured (a \code{gwr_fit}, a
#'   \code{bayesian_fit}, a forest fitted with \code{include_coords = TRUE},
#'   or \code{coords = "yes"}) an \code{sf} object is required,
#'   non-\code{POINT} geometry is reduced to representative points, and a CRS
#'   mismatch with the training data is reconciled; see \emph{Models whose
#'   predictions depend on location}.
#' @param model A fitted \code{spatial_fit}, supplying the training data and
#'   predictor names, the importance \code{weights = "auto"} reads, and
#'   whether its predictions depend on location. Optional if \code{train_sf}
#'   and \code{predictor_vars} are given directly, which lets this be used
#'   with any model.
#' @param train_sf Training data, if not taken from \code{model}.
#' @param predictor_vars Predictor names, if not taken from \code{model}.
#' @param weights \code{"auto"} (default), \code{NULL}, or a named numeric
#'   vector of predictor importances. \code{"auto"} uses
#'   \code{pmax(model$info$importance, 0)} when \code{model} carries a finite
#'   importance with at least one positive value, naming every predictor in
#'   use, from a model with out-of-bag skill; otherwise every predictor counts
#'   equally and a logged note says why (see \emph{Weights}). \code{NULL} is
#'   the explicit request for equal weights. For a vector of your own, any
#'   positive scale works. When \code{model} was fitted with
#'   \code{include_coords = TRUE}, the coordinates \code{"..x"} and
#'   \code{"..y"} are measured as predictors too (see below), and you are not
#'   expected to supply importances for them: any you leave out default to the
#'   mean of the weights you did supply, so location counts about as much as a
#'   typical predictor. Naming them explicitly overrides that. An unnamed
#'   vector may have one value per predictor either with or without the two
#'   coordinate columns. Weights must be finite and non-negative, so pass
#'   permutation importance as \code{pmax(importance, 0)}. (A forest with no
#'   out-of-bag rows, \code{replace = FALSE} with \code{sample_fraction = 1},
#'   has \code{NaN} importance, which \code{pmax()} keeps and which is
#'   refused when passed by hand; \code{"auto"} uses equal weights for it.)
#'   \code{pmax(importance, 0)} is all zero when the model found no
#'   predictor useful, and then the weights cannot say anything: with a
#'   single predictor any weight gives the same index and zero
#'   is accepted; with several, all of them are weighted equally, as with
#'   \code{weights = NULL}, and a warning says so. The coordinate default
#'   above is the mean of the supplied weights, so a zero weight on the only
#'   covariate of a coordinate-using model zeroes the coordinates too, and
#'   that equal weighting applies.
#' @param folds Cross-validation folds: a \code{\link{make_folds}} result, a
#'   list of \code{train}/\code{test} splits, or a vector of fold labels with
#'   one entry per training row. Default \code{NULL} (plain nearest neighbour).
#'   The folds you passed to \code{cv_*()}, built on the layer \code{model}
#'   was fitted from, may name rows that \code{prep_model_data()} removed
#'   (a missing or non-finite value, an empty geometry): as in \code{cv_*()},
#'   they are dropped from the folds, and a label vector with one entry per
#'   row of that layer loses theirs. Also
#'   as in \code{cv_*()}, a \code{make_folds()} result built on other data --
#'   another layer, or these rows in another order -- is refused. That check
#'   is skipped when the folds were built on polygons and the training data
#'   are the points a fit reduced them to.
#' @param threshold Optional numeric override for the DI threshold of the
#'   predictor index. The location index, when there is one, always takes its
#'   threshold from the rule.
#' @param normalizer_max_n Subsample the training data to this many points when
#'   computing the mean pairwise distance, which is quadratic. Default 5000.
#' @param seed Seed for that subsample. Default 123.
#' @param chunk_size Query rows per distance block on the dense path. Default
#'   \code{NULL} (chosen from the training size). Otherwise a single number of
#'   at least 1; a fractional value is truncated to a whole number of rows.
#' @param use_fnn Use \pkg{FNN} for nearest-neighbour search when available.
#'   Exposed so the dense fallback can be tested. The two paths return the
#'   same distances: the dense one searches with a fast formula and then
#'   recomputes the nearest distance exactly, so an exact copy of a training
#'   row is at distance 0 on both.
#' @param coords Whether location is measured: \code{"auto"} (default)
#'   measures it for a model whose predictions depend on it (a
#'   \code{gwr_fit}, a \code{bayesian_fit}, a forest fitted with
#'   \code{include_coords = TRUE}) and not otherwise; \code{"yes"} measures
#'   it whatever the model, and with no \code{model}; \code{"no"} never
#'   does, and logs a caution for a model that depends on it. See
#'   \emph{Models whose predictions depend on location}.
#' @param categorical What to do with factor and character predictors:
#'   \code{"refuse"} (default) stops with a message naming them and the two
#'   alternatives; \code{"indicator"} codes each level as a 0/1 column and
#'   marks a prediction row on a level no training row has as outside
#'   (\code{DI = Inf}); \code{"drop"} computes the index on the other
#'   predictors and logs which were left out. See \emph{Limitations} for the
#'   scaling \code{"indicator"} uses.
#'
#' @return An object of class \code{aoa}: a list with
#'   \itemize{
#'     \item \code{aoa}: \code{newdata} with a numeric \code{DI} column
#'       and a logical \code{AOA} column added, and between them a numeric
#'       \code{DI_location} column when location has an index of its own.
#'       For a forest fitted with \code{include_coords = TRUE} it is the
#'       object the computation ran on: \code{newdata} after pointizing, CRS
#'       reconciliation and the addition of the \code{"..x"} and \code{"..y"}
#'       columns. With a location index it is \code{newdata} as passed. A
#'       row whose predictors are not all finite gets \code{NA} in \code{DI},
#'       and a row outside the training range of a predictor in
#'       \code{dropped_vars}, or on a level no training row has, gets
#'       \code{DI = Inf} and \code{AOA = FALSE} (see \emph{Limitations}).
#'       \code{AOA} is \code{DI <= threshold}, and with a location index
#'       also \code{DI_location <= location$threshold}: a row outside on
#'       either index is \code{FALSE}, and a row with one index missing and
#'       the other inside is \code{NA}.
#'     \item \code{threshold}: the DI cut-off used.
#'     \item \code{train_DI}: the training points' own DI values.
#'     \item \code{normalizer}: the mean pairwise training distance.
#'     \item \code{weights}: the weight vector actually applied, named by
#'       \code{predictor_vars}.
#'     \item \code{predictor_vars}: the predictors used, including
#'       \code{"..x"}/\code{"..y"} when a forest uses the coordinates and
#'       excluding \code{dropped_vars} and, under
#'       \code{categorical = "drop"}, \code{categorical_vars}.
#'     \item \code{dropped_vars}: predictors dropped for negligible
#'       variance.
#'     \item \code{scaling}: a list with \code{center} and \code{scale},
#'       each named by \code{predictor_vars}: the training means and standard
#'       deviations the index is computed in, so a location's DI can be
#'       traced to the predictor that put it outside. A predictor coded as
#'       indicators has one entry per level, named
#'       \code{"<predictor>=<level>"}.
#'     \item \code{n_outliers}: the number of training DI values above the
#'       \code{Q3 + 1.5 * IQR} fence, which the default threshold rule sets
#'       aside (the "outlier-removed" in its name); computed whether or not
#'       \code{threshold} was supplied.
#'     \item \code{n_duplicates}: with \code{folds = NULL}, the number of
#'       training rows that repeat an earlier row's predictor values exactly
#'       and so share its reference distance; \code{NA} when folds were
#'       supplied (see \emph{Training rows with identical predictor values}).
#'     \item \code{location}: \code{NULL}, or for a location index a list
#'       with \code{threshold}, \code{train_DI}, \code{normalizer} (the mean
#'       distance between training points, in map units, so
#'       \code{threshold * normalizer} is the cut-off as a distance),
#'       \code{n_duplicates} (training rows at an earlier row's location;
#'       \code{NA} when folds were supplied) and \code{n_outside} (prediction
#'       rows beyond the location threshold).
#'     \item \code{categorical_vars}: the factor and character predictors
#'       found, coded or left out as \code{params$categorical} says.
#'     \item \code{n_train}, \code{n_new}, \code{n_inside},
#'       \code{n_outside}, \code{n_na}: row counts. \code{n_train} counts the
#'       training rows that survived the finite-value filter; \code{n_new} is
#'       every row of \code{newdata}, so
#'       \code{n_new = n_inside + n_outside + n_na}, with \code{n_na} the
#'       rows whose \code{AOA} is \code{NA}.
#'       \code{n_outside_predictors} counts the rows beyond the predictor
#'       threshold, which is \code{n_outside} unless there is a location
#'       index.
#'     \item \code{params} records the call: \code{folds_supplied},
#'       \code{n_folds}, \code{folds_method}, \code{threshold_supplied},
#'       \code{normalizer_max_n}, \code{normalizer_n_used},
#'       \code{normalizer_subsampled}, \code{weights_supplied} (the caller
#'       passed a vector), \code{weights_source} (\code{"model importance"},
#'       \code{"supplied"} or \code{"equal"}), \code{weights_note} (the kind
#'       of importance, or why the weights are equal), \code{coords}
#'       (\code{"location index"}, \code{"predictors"} or \code{"none"}),
#'       \code{categorical} and
#'       \code{seed}.  \code{folds_method} is the \code{method} of a
#'       \code{\link{make_folds}()} result (\code{"block_kfold"},
#'       \code{"random_kfold"}, ...), \code{"labels"} for a vector of fold
#'       labels, \code{"splits"} for a bare list of train/test splits, and
#'       \code{NA} when no folds were supplied.  It is printed with the
#'       object, because the threshold's meaning depends on it: the
#'       cross-validated DI that sets it comes from the same kind of hold-out
#'       as the CV error it should be quoted beside, so an AOA built on
#'       \code{random_kfold} folds is logged as a caution and does not
#'       belong next to a blocked \code{cv_*()} result.
#'   }
#'
#' @references
#' Meyer, H. and Pebesma, E. (2021). Predicting into unknown space? Estimating
#' the area of applicability of spatial prediction models. \emph{Methods in
#' Ecology and Evolution} 12(9), 1620--1633. \doi{10.1111/2041-210X.13650}
#'
#' @seealso \code{\link{predict_surface}} to build the grid,
#'   \code{\link{make_folds}} for the fold scheme.
#' @family cross-validation
#' @family prediction
#' @examples
#' library(sf)
#' set.seed(1)
#' n <- 120
#' train <- st_as_sf(
#'   data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'              a = rnorm(n), b = rnorm(n)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' train$z <- 2 * train$a - train$b + rnorm(n, 0, 0.3)
#'
#' # Prediction points, some of them well outside the training predictor range
#' newpts <- st_as_sf(
#'   data.frame(x = 5e5 + runif(50, 0, 1000), y = 5e6 + runif(50, 0, 1000),
#'              a = c(rnorm(40), rnorm(10, 8)), b = rnorm(50)),
#'   coords = c("x", "y"), crs = 32632
#' )
#'
#' res <- area_of_applicability(newpts, train_sf = train,
#'                              predictor_vars = c("a", "b"))
#' res
#' table(res$aoa$AOA)
#' @export
area_of_applicability <- function(newdata, model = NULL, train_sf = NULL,
                                  predictor_vars = NULL, weights = "auto",
                                  folds = NULL, threshold = NULL,
                                  normalizer_max_n = 5000L, seed = 123L,
                                  chunk_size = NULL,
                                  use_fnn = requireNamespace("FNN",
                                                             quietly = TRUE),
                                  coords = c("auto", "yes", "no"),
                                  categorical = c("refuse", "indicator",
                                                  "drop")) {
  coords      <- match.arg(coords)
  categorical <- match.arg(categorical)
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

  # chunk_size reached seq.int(1, nq, by = chunk_size) unvalidated.  A computed
  # value such as 1e3 / nrow(train) = 12.5 left some rows of each block at the
  # 0 they were initialised with, so extrapolation read as inside the AOA; 0,
  # NA or a length-2 vector failed in base R without naming the argument.
  if (!is.null(chunk_size)) {
    .check_scalar(chunk_size, "chunk_size", "area_of_applicability", min = 1,
                  max = .Machine$integer.max, what = "a single positive number")
    chunk_size <- as.integer(chunk_size)
  }

  # Fold train/test slots from make_folds() are ..row_id VALUES, not row
  # positions; keep the IDs alongside so .aoa_fold_splits() can resolve them.
  tr_meta <- if (inherits(train_sf, "sf")) sf::st_drop_geometry(train_sf) else
    as.data.frame(train_sf)
  tr_row_ids <- if ("..row_id" %in% names(tr_meta)) tr_meta[["..row_id"]] else NULL
  removed_ids <- NULL
  if (!is.null(folds)) {
    fr <- .aoa_rows_for_folds(folds, train_sf, tr_row_ids)
    folds       <- fr$folds
    tr_row_ids  <- fr$row_ids
    removed_ids <- fr$removed_ids
    # The provenance check cv_*() make: fold splits are row IDs, so folds
    # built on another layer of the same size, or on these rows in another
    # order, applied silently and moved the threshold.
    .aoa_check_fold_probe(folds, train_sf, tr_row_ids)
  }

  # A model whose predictions depend on location has to have location
  # measured.  Without it, a prediction point far outside the training extent
  # but with ordinary covariate values reads as INSIDE the area of
  # applicability -- exactly the extrapolation the index exists to catch.
  # There are two kinds of such model, and location enters the index
  # differently for each:
  #   "predictors"  a forest fitted with include_coords = TRUE splits on the
  #                 coordinates like on any column, so they join the
  #                 predictor index as "..x"/"..y" and take a weight;
  #   "index"       a GWR's coefficients and a Bayesian fit's Gaussian
  #                 process are functions of location whatever the
  #                 covariates are, so location gets an index of its own --
  #                 distance to the nearest training point in map space --
  #                 and a point must be inside on both.  As two more
  #                 equally weighted columns it would be diluted: with 12
  #                 covariates, 89% of points 200 to 500 units beyond a
  #                 1000-unit extent still read as inside.
  xy_predictors  <- !is.null(model) && isTRUE(model$info$include_coords)
  location_model <- !is.null(model) &&
    inherits(model, c("gwr_fit", "bayesian_fit"))
  coord_mode <- switch(coords,
    auto = if (xy_predictors) "predictors" else if (location_model) "index"
           else "none",
    yes  = if (xy_predictors) "predictors" else "index",
    no   = "none")
  if (identical(coords, "no") && (xy_predictors || location_model))
    .log_warn(paste0("area_of_applicability(): coords = \"no\", so location is ",
                     "not measured although this model's predictions depend on ",
                     "it (%s). A point far from every training point reads as ",
                     "inside when its covariates are ordinary."),
              if (xy_predictors) "a forest fitted with include_coords = TRUE"
              else sprintf("a %s", class(model)[1L]))
  coord_model <- identical(coord_mode, "predictors")
  coord_vars  <- if (coord_model) c("..x", "..y") else character(0)
  loc_xy <- if (identical(coord_mode, "index"))
    .aoa_location_coords(train_sf, newdata,
                         why = if (location_model)
                           sprintf("this model's predictions depend on location (a %s)",
                                   class(model)[1L])
                         else "coords = \"yes\"")
  else NULL
  if (coord_model) {
    if (!inherits(train_sf, "sf") || !inherits(newdata, "sf"))
      stop("area_of_applicability(): this model uses the coordinates as ",
           "predictors, so both the training data and `newdata` must be sf ",
           "objects carrying geometry.", call. = FALSE)
    # Both directions matter.  st_transform() on CRS-less newdata aborts with
    # sf's bare "cannot transform sfc object with missing crs", and when the
    # TRAINING data has no CRS the guard short-circuited entirely, so degrees
    # were fed into a metre-space index with nothing said (DI 0.13/0.19 for the
    # same points that gave 0.24/0.34 correctly projected).
    tr_crs <- sf::st_crs(train_sf)
    if (!is.na(tr_crs)) {
      if (!identical(sf::st_crs(newdata), tr_crs))
        newdata <- .transform_or_stamp(newdata, tr_crs, what = "newdata",
                                       caller = "area_of_applicability")
    } else if (!is.na(sf::st_crs(newdata))) {
      train_sf <- .transform_or_stamp(train_sf, sf::st_crs(newdata),
                                      what = "the training data",
                                      caller = "area_of_applicability")
    }
    train_sf <- .aoa_add_coords(train_sf)
    newdata  <- .aoa_add_coords(newdata)
    predictor_vars <- unique(c(predictor_vars, "..x", "..y"))
  }

  des <- .aoa_design(train_sf, newdata, predictor_vars, categorical)
  X_tr_full <- des$X_tr
  X_nw_full <- des$X_nw
  cat_vars  <- des$cat_vars
  if (length(cat_vars) > 0L && identical(categorical, "drop"))
    .log_warn(paste0("area_of_applicability(): categorical = \"drop\": the ",
                     "index is computed on %s alone. %s %s not measured, so a ",
                     "prediction row on a level the training data never has ",
                     "reads as inside when its other predictors are ordinary."),
              paste(des$vars, collapse = ", "),
              paste(sQuote(cat_vars), collapse = ", "),
              if (length(cat_vars) == 1L) "is" else "are")
  predictor_vars <- des$vars
  col_var        <- des$col_var

  # Training rows carrying NA or Inf cannot define a reference distance.
  # complete.cases() alone would let an Inf through, and it would then poison
  # that predictor's mean and sd, so .aoa_scaling() would drop the whole
  # predictor and report it as having "no variance" -- a false diagnosis that
  # silently removes a dimension the index is supposed to be measuring.
  tr_ok <- apply(X_tr_full, 1L, function(r) all(is.finite(r)))
  # A row with no usable location has no place in a location index either.
  if (!is.null(loc_xy))
    tr_ok <- tr_ok & is.finite(loc_xy$train[, 1L]) & is.finite(loc_xy$train[, 2L])
  if (!all(tr_ok)) {
    if (sum(tr_ok) < 2L)
      stop("area_of_applicability(): fewer than two training rows with finite ",
           "values for every predictor.", call. = FALSE)
    .log_warn("area_of_applicability(): dropping %d training row(s) with missing or non-finite predictor values%s.",
              sum(!tr_ok),
              if (!is.null(loc_xy)) " or no usable location" else "")
    if (!is.null(folds))
      stop("area_of_applicability(): the training data has missing predictor ",
           "values, so `folds` row positions can no longer be matched. Remove ",
           "the incomplete rows before building folds.", call. = FALSE)
    X_tr_full <- X_tr_full[tr_ok, , drop = FALSE]
    if (!is.null(loc_xy)) loc_xy$train <- loc_xy$train[tr_ok, , drop = FALSE]
    # Unreachable when folds were supplied (the branch above stops), but keeps
    # the IDs aligned with X_tr_full for anything added later.
    if (!is.null(tr_row_ids)) tr_row_ids <- tr_row_ids[tr_ok]
  }

  sc <- .aoa_scaling(X_tr_full, col_var = col_var, indicator = des$indicator)
  # Columns and predictors are the same thing unless a categorical predictor
  # was coded as indicators: then a predictor is dropped when none of its
  # columns has variance (one level among the training rows), and a single
  # column of a predictor that is kept can be dropped too (a level whose only
  # rows were removed above).
  dropped_cols <- names(sc$keep)[!sc$keep]
  used_vars    <- unique(col_var[sc$keep])
  # Named by predictor, or by column when the predictor itself is kept.
  # Subscripted assignment, not ifelse(): with nothing dropped ifelse()
  # returns logical(0), and `dropped_vars` is documented, and was always
  # returned, as a character vector.
  dropped <- as.character(col_var[!sc$keep])
  one_col <- dropped %in% used_vars
  dropped[one_col] <- dropped_cols[one_col]
  dropped <- unique(dropped)
  if (length(dropped) > 0L)
    .log_warn(paste0("area_of_applicability(): dropping predictor(s) with no ",
                     "variance in the training data: %s. A prediction point ",
                     "taking a different value there is extrapolation along ",
                     "it, and is marked outside the AOA with DI = Inf."),
              paste(dropped, collapse = ", "))
  used_cols <- names(sc$keep)[sc$keep]

  Z_tr <- .aoa_apply_scaling(X_tr_full, sc)
  Z_nw <- .aoa_apply_scaling(X_nw_full, sc)

  # Where the weights come from.  "auto" reads the model's own variable
  # importance when it has one that means something; NULL is the explicit
  # request for equal weights; anything else is the caller's vector.
  if (identical(weights, "auto")) {
    aw       <- .aoa_auto_weights(model, predictor_vars, fill_vars = coord_vars)
    weights  <- aw$weights
    w_source <- aw$source
    w_note   <- aw$note
    w_given  <- FALSE
  } else if (is.null(weights)) {
    w_source <- "equal"
    w_note   <- "weights = NULL"
    w_given  <- FALSE
  } else {
    w_source <- "supplied"
    w_note   <- NA_character_
    w_given  <- TRUE
  }

  # Weighting convention: the scaled COLUMN is multiplied by the weight, so a
  # weight w contributes w^2 to the squared distance.  This is what CAST --
  # the reference implementation, by the paper's own authors -- does, and
  # anyone cross-checking against it would otherwise find quietly different
  # numbers.  Multiplying by sqrt(w) instead (contributing w to the squared
  # distance) is the other defensible reading of "weighted Euclidean" and is
  # NOT what is used here.
  # The weights that matter are those of the predictors kept; if the only
  # non-zero ones sat on predictors dropped above, the rest are all zero.
  w_kept <- .aoa_weight_vector(weights, predictor_vars,
                               fill_vars = coord_vars)[used_vars]
  if (!identical(w_source, "equal") && length(w_kept) > 0L &&
      all(w_kept == 0)) {
    w_note   <- sprintf("every %s weight on a predictor kept is zero",
                        if (w_given) "supplied" else "importance")
    w_source <- "equal"
  }
  w     <- .aoa_equal_if_all_zero(w_kept)
  # One weight per predictor; an indicator-coded predictor's columns share it.
  w_vec <- unname(w[col_var[sc$keep]])
  Z_tr  <- sweep(Z_tr, 2L, w_vec, "*")
  Z_nw  <- sweep(Z_nw, 2L, w_vec, "*")

  splits <- .aoa_fold_splits(folds, nrow(Z_tr), row_ids = tr_row_ids,
                             removed_ids = removed_ids)

  # What KIND of folds these are is part of what the threshold means: Meyer
  # and Pebesma define it from the cross-validated training DI, so a threshold
  # from random folds pairs with a random-CV error estimate and a threshold
  # from blocked folds with a blocked one.  make_folds() records the method;
  # a bare label vector or split list cannot say, and is recorded as such.
  folds_method <- .aoa_folds_method(folds)
  if (identical(folds_method, "random_kfold"))
    .log_warn(paste0("area_of_applicability(): `folds` are random_kfold. The DI ",
                     "threshold is then the largest dissimilarity a random ",
                     "hold-out produced, which under spatial autocorrelation ",
                     "is shorter than what a spatially blocked hold-out ",
                     "produces; this AOA pairs with a random-CV error ",
                     "estimate, not with a blocked one (Meyer & Pebesma ",
                     "2021). Build the folds with make_folds(method = ",
                     "\"block_kfold\") to match a blocked cv_*() result."))

  norm     <- .aoa_normalizer(Z_tr, max_n = normalizer_max_n, seed = seed)
  # Without folds a row that repeats another's predictor values exactly is
  # the same point of predictor space, not its neighbour; see
  # .aoa_train_dist().
  groups   <- if (is.null(splits)) .aoa_row_groups(Z_tr) else NULL
  n_dup    <- if (is.null(groups)) NA_integer_ else as.integer(groups$n_dup)
  train_d  <- .aoa_train_dist(Z_tr, splits = splits, use_fnn = use_fnn,
                              chunk_size = chunk_size, groups = groups)
  train_DI <- train_d / norm$value
  if (isTRUE(n_dup > 0L))
    .log_info(paste0("area_of_applicability(): %d of %d training rows repeat ",
                     "another row's predictor values exactly (%d distinct ",
                     "points). A repeat is the same point of predictor space, ",
                     "not a neighbour, so each row's reference is the nearest ",
                     "row with different values."),
              n_dup, nrow(Z_tr), nrow(Z_tr) - n_dup)

  thr <- if (is.null(threshold)) .aoa_threshold(train_DI) else {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold <= 0)
      stop("area_of_applicability(): `threshold` must be a single positive ",
           "number.", call. = FALSE)
    as.numeric(threshold)
  }

  # With folds, each training row's reference is the nearest row its fold
  # left in training, so a row with an exact duplicate on the other side of
  # the split -- repeat measurements at a site with static covariates,
  # covariates read off a raster coarser than the sampling -- has a training
  # DI of 0.  That is the definition, not an error: the hold-out really did
  # have an identical row to learn from.  But once about three quarters of
  # the rows have one, Q3 and the IQR are 0 and so is the threshold, and only
  # exact copies of a training row are inside.  A result that degenerate is
  # raised as a warning with a class, so a script or a cv loop can catch it.
  if (is.null(threshold) && thr == 0)
    .aoa_warn_zero_threshold(
      paste0("area_of_applicability(): the DI threshold is 0 because %d of %d ",
             "training rows have an exact duplicate in predictor space among ",
             "the rows their fold trains on (training DI = 0), so only ",
             "prediction points identical to a training row count as inside. ",
             "Repeat measurements at a site, or covariates coarser than the ",
             "sampling, do this when the folds split them. Pass folds that ",
             "keep the duplicates together (make_folds(method = ",
             "\"leave_location_out\", group_var = ...)), or folds = NULL, ",
             "which measures each row against the nearest row with different ",
             "values."),
      sum(train_DI == 0), length(train_DI))

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

  # A level no training row has: the model has nothing to say about it, and
  # no distance to a training row is defined along that predictor.
  if (ncol(des$unseen) > 0L) {
    unseen <- rowSums(des$unseen) > 0L & !is.na(DI)
    if (any(unseen)) {
      DI[unseen] <- Inf
      .warn_and_log(paste0("area_of_applicability(): %d of %d prediction ",
                           "row(s) carry a level of %s that no training row ",
                           "has; they are marked outside the AOA with DI = Inf."),
                    sum(unseen), length(DI),
                    paste(sQuote(colnames(des$unseen)[
                      colSums(des$unseen[unseen, , drop = FALSE]) > 0L]),
                      collapse = ", "))
    }
  }

  # A predictor dropped for having no training variance is not in the
  # distance, so a prediction row taking a different value there was judged
  # on the other predictors alone: a land-cover dummy that is 0 throughout the
  # training region put 37 of 40 urban rows inside the AOA, while ONE urban
  # training row would have kept the predictor and put 1 inside.  Along that
  # direction the scaled distance is (x - c) / 0, infinite, so that is the DI.
  # "Different" means outside the training range by more than the relative
  # tolerance .aoa_scaling() drops the predictor with; a missing value cannot
  # be compared and leaves the row to the other predictors.
  if (length(dropped_cols) > 0L) {
    Xd_tr <- X_tr_full[, dropped_cols, drop = FALSE]
    lo    <- apply(Xd_tr, 2L, min)
    hi    <- apply(Xd_tr, 2L, max)
    mag   <- pmax(abs(lo), abs(hi))
    mag[!is.finite(mag) | mag <= 0] <- 1
    slack <- .Machine$double.eps^0.5 * mag
    Xd_nw <- X_nw_full[, dropped_cols, drop = FALSE]
    off   <- sweep(Xd_nw, 2L, lo - slack, "<") | sweep(Xd_nw, 2L, hi + slack, ">")
    off[is.na(off)] <- FALSE
    # Rows already outside on an unseen level are not counted a second time.
    beyond <- rowSums(off) > 0L & !is.na(DI) & is.finite(DI)
    if (any(beyond)) {
      DI[beyond] <- Inf
      .warn_and_log(paste0("area_of_applicability(): %d of %d prediction ",
                           "row(s) take a value the training data never has on ",
                           "%s, dropped for having no training variance; they ",
                           "are extrapolation along it and are marked outside ",
                           "the AOA with DI = Inf."),
                    sum(beyond), length(DI),
                    paste(sQuote(dropped_cols[colSums(off[beyond, , drop = FALSE]) > 0L]),
                          collapse = ", "))
    }
  }

  inside <- DI <= thr

  # The location index: the same construction as the predictor index, in map
  # space.  Distance to the nearest training point, over the mean pairwise
  # distance between training points; the training points' own values taken
  # the same way (same folds, same rule for repeated locations); the same
  # outlier-removed maximum as the threshold.  The coordinates are centred
  # and not scaled: map distance is isotropic, and a study area longer than
  # it is wide is not a reason to count a kilometre north as more than a
  # kilometre east.
  location <- NULL
  DI_loc   <- NULL
  if (!is.null(loc_xy)) {
    xy_ctr <- colMeans(loc_xy$train)
    C_tr   <- sweep(loc_xy$train, 2L, xy_ctr, "-")
    C_nw   <- sweep(loc_xy$new,   2L, xy_ctr, "-")
    norm_loc <- .aoa_normalizer(C_tr, max_n = normalizer_max_n, seed = seed,
                                what = "location")
    loc_groups <- if (is.null(splits)) .aoa_row_groups(C_tr) else NULL
    train_DI_loc <- .aoa_train_dist(C_tr, splits = splits, use_fnn = use_fnn,
                                    chunk_size = chunk_size,
                                    groups = loc_groups) / norm_loc$value
    thr_loc <- .aoa_threshold(train_DI_loc)
    if (thr_loc == 0)
      .aoa_warn_zero_threshold(
        paste0("area_of_applicability(): the location threshold is 0 because ",
               "%d of %d training rows share their location with a row their ",
               "fold trains on, so only prediction points at a training ",
               "location count as inside. Repeat measurements at a site do ",
               "this when the folds split them. Pass folds that keep a site ",
               "together (make_folds(method = \"leave_location_out\", ",
               "group_var = ...)), or folds = NULL, which measures each site ",
               "against the nearest other site."),
        sum(train_DI_loc == 0), length(train_DI_loc))
    loc_ok <- is.finite(C_nw[, 1L]) & is.finite(C_nw[, 2L])
    DI_loc <- rep(NA_real_, nrow(C_nw))
    if (any(loc_ok))
      DI_loc[loc_ok] <- .aoa_min_dist(C_nw[loc_ok, , drop = FALSE], C_tr,
                                      use_fnn = use_fnn,
                                      chunk_size = chunk_size) / norm_loc$value
    inside_loc <- DI_loc <= thr_loc
    # Inside means inside on both.  R's three-valued `&` is the right rule
    # for a missing value: a row outside on one index is outside whatever
    # the other would have said, and a row inside on one with the other
    # missing is not classified.
    inside <- inside & inside_loc
    location <- list(
      threshold    = thr_loc,
      train_DI     = train_DI_loc,
      normalizer   = norm_loc$value,
      n_duplicates = if (is.null(loc_groups)) NA_integer_
                     else as.integer(loc_groups$n_dup),
      n_outside    = sum(!inside_loc, na.rm = TRUE)
    )
  }

  out <- newdata
  out$DI  <- DI
  if (!is.null(DI_loc)) out$DI_location <- DI_loc
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
      # The space the index is computed in: each column's training centre
      # and scale.  Without them a location's DI cannot be traced to the
      # predictor that put it outside.
      scaling        = list(center = sc$center[used_cols],
                            scale  = sc$scale[used_cols]),
      # Training DI values above the Q3 + 1.5 IQR fence: the ones the
      # default threshold rule sets aside as outliers.
      n_outliers     = .aoa_fence_outliers(train_DI),
      # Training rows that repeat another row's predictor values exactly and
      # so share its reference distance; NA when folds set the reference.
      n_duplicates   = n_dup,
      # The location index, when one was measured; NULL otherwise.
      location       = location,
      categorical_vars = cat_vars,
      n_train        = nrow(Z_tr),
      n_new          = nrow(Z_nw),
      n_inside       = sum(inside, na.rm = TRUE),
      n_outside      = sum(!inside, na.rm = TRUE),
      n_na           = sum(is.na(inside)),
      n_outside_predictors = sum(DI > thr, na.rm = TRUE),
      params = list(
        folds_supplied      = !is.null(splits),
        n_folds             = if (is.null(splits)) 0L else length(splits),
        folds_method        = folds_method,
        threshold_supplied  = !is.null(threshold),
        normalizer_max_n    = normalizer_max_n,
        normalizer_n_used   = norm$n_used,
        normalizer_subsampled = norm$subsampled,
        weights_supplied    = w_given,
        weights_source      = w_source,
        weights_note        = w_note,
        coords              = switch(coord_mode, predictors = "predictors",
                                     index = "location index", "none"),
        categorical         = categorical,
        seed                = seed
      )
    ),
    class = "aoa"
  )
}


#' The fold method behind an AOA threshold
#'
#' \code{make_folds()} results carry \code{$method}; a vector of labels
#' (\code{blockCV::cv_spatial()$folds_ids}, say) or a bare list of splits
#' cannot say how they were built.
#'
#' @return Character(1): the method, \code{"labels"}, \code{"splits"}, or
#'   \code{NA_character_} when \code{folds} is \code{NULL}.
#' @keywords internal
#' @noRd
.aoa_folds_method <- function(folds) {
  if (is.null(folds)) return(NA_character_)
  if (is.atomic(folds) && !is.list(folds)) return("labels")
  if (is.list(folds) && !is.null(folds$method) &&
      is.character(folds$method) && length(folds$method) == 1L &&
      !is.na(folds$method) && nzchar(folds$method))
    return(folds$method)
  "splits"
}

#' @keywords internal
#' @noRd
.aoa_folds_label <- function(m) {
  if (is.null(m) || is.na(m)) return("method unknown")
  switch(m,
         labels = "from fold labels; method unknown",
         splits = "from train/test splits; method unknown",
         m)
}


#' Print an area-of-applicability result
#'
#' Summarises where the model may be trusted: how many prediction locations
#' fall inside the area of applicability and how many outside, the
#' dissimilarity threshold that separated them, and the predictors the index
#' was computed over (naming any dropped for having no usable variance).  It
#' also says where the weights came from (the model's importance, a vector
#' the caller supplied, or equal weights and why), whether location was
#' measured and how, how categorical predictors were handled, and how many
#' training rows repeat another's predictor values.  With a location index
#' the rows outside are split into those outside on the predictors and those
#' outside on location.  The
#' proportion outside is the headline number. A map that extrapolates over
#' much of its extent is reporting predictions its training data cannot
#' support, whatever the cross-validation score said.
#'
#' @param x An \code{aoa} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @family print methods
#' @export
print.aoa <- function(x, ...) {
  cat("Area of applicability (Meyer & Pebesma 2021)\n\n")
  cat(sprintf("  predictors  : %d (%s)\n", length(x$predictor_vars),
              paste(x$predictor_vars, collapse = ", ")))
  if (length(x$dropped_vars) > 0L)
    cat(sprintf("  dropped     : %s (no training variance)\n",
                paste(x$dropped_vars, collapse = ", ")))
  if (length(x$categorical_vars) > 0L)
    cat(sprintf("  categorical : %s %s\n",
                paste(x$categorical_vars, collapse = ", "),
                if (identical(x$params$categorical, "drop"))
                  "left out of the index (categorical = \"drop\")"
                else "coded as one 0/1 column per level"))
  # Where the weights came from.  An object built before the source was
  # recorded says only whether the caller supplied any.
  src <- x$params$weights_source
  note <- x$params$weights_note
  note <- if (is.character(note) && length(note) == 1L && !is.na(note) &&
              nzchar(note)) note else NULL
  cat(sprintf("  weighted    : %s\n",
              if (is.null(src)) {
                if (isTRUE(x$params$weights_supplied)) "yes" else
                  "no (all predictors count equally)"
              } else switch(src,
                "model importance" = sprintf(
                  "yes, by the model's %s, floored at 0",
                  note %||% "variable importance"),
                supplied = "yes, by the weights supplied",
                sprintf("no (all predictors count equally%s)",
                        if (is.null(note)) "" else paste0(": ", note)))))
  cat(sprintf("  training    : %d points\n", x$n_train))
  n_dup <- suppressWarnings(as.integer(x$n_duplicates %||% NA_integer_))
  cat(sprintf("  reference   : %s\n",
              if (isTRUE(x$params$folds_supplied))
                sprintf("nearest point outside each of %d CV folds (%s)",
                        x$params$n_folds,
                        .aoa_folds_label(x$params$folds_method))
              else if (isTRUE(n_dup > 0L))
                sprintf(paste0("nearest training point with different values ",
                               "(no folds supplied; %d of %d rows repeat ",
                               "another's values)"), n_dup, x$n_train)
              else "nearest other training point (no folds supplied)"))
  cat(sprintf("  normaliser  : %.4f (mean pairwise distance%s)\n",
              x$normalizer,
              if (isTRUE(x$params$normalizer_subsampled))
                sprintf(", %d-point subsample", x$params$normalizer_n_used)
              else ""))
  cat(sprintf("  threshold   : %.4f%s\n", x$threshold,
              if (isTRUE(x$params$threshold_supplied)) " (supplied)" else
                " (outlier-removed max of training DI)"))
  # A zero threshold reads as "nothing is inside" with no reason given; the
  # reason is duplicated training rows, and the count says how many.
  n_zero <- sum(x$train_DI == 0, na.rm = TRUE)
  if (!isTRUE(x$params$threshold_supplied) && isTRUE(x$threshold == 0))
    cat(sprintf(paste0("                (%d of %d training DI are 0: exact ",
                       "duplicates in predictor space)\n"),
                n_zero, length(x$train_DI)))
  loc <- x$location
  if (is.list(loc)) {
    cat(sprintf(paste0("  location    : measured on its own; threshold %.4f ",
                       "(%s map units from the nearest training point)\n"),
                loc$threshold, format(signif(loc$threshold * loc$normalizer, 4))))
    if (isTRUE(loc$threshold == 0))
      cat(sprintf(paste0("                (%d of %d training rows share a ",
                         "location with a row their fold trains on)\n"),
                  sum(loc$train_DI == 0, na.rm = TRUE), length(loc$train_DI)))
  } else if (identical(x$params$coords, "predictors")) {
    cat("  location    : measured as the predictors ..x and ..y\n")
  } else if (!is.null(x$params$coords)) {
    cat("  location    : not measured\n")
  }
  cat("\n")

  pct <- if (x$n_new > 0L) 100 * x$n_inside / x$n_new else NA_real_
  cat(sprintf("  %d of %d prediction points inside the AOA (%.1f%%)\n",
              x$n_inside, x$n_new, pct))
  if (is.list(loc) && x$n_outside > 0L)
    cat(sprintf(paste0("  %d outside: %d on the predictors, %d on location ",
                       "(a point can be both)\n"),
                x$n_outside, x$n_outside_predictors %||% NA_integer_,
                loc$n_outside))
  if (x$n_na > 0L)
    cat(sprintf("  %d %s\n", x$n_na,
                if (is.list(loc))
                  "not classified (a predictor or the location is missing)"
                else "with missing predictors (DI = NA)"))
  n_inf <- sum(is.infinite(x$aoa$DI))
  if (n_inf > 0L)
    cat(sprintf("  %d outside on %s (DI = Inf)\n", n_inf,
                if (length(x$categorical_vars) > 0L &&
                    identical(x$params$categorical, "indicator"))
                  "a level no training row has, or on a dropped predictor"
                else "a dropped predictor"))
  if (x$n_outside > 0L)
    cat("\nPredictions outside the AOA are extrapolations; the ",
        "cross-validated\nperformance estimate does not cover them.\n",
        sep = "")
  invisible(x)
}
