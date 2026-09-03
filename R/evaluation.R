# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Filter a named argument list to only parameters accepted by a function
#' @keywords internal
#' @noRd
.filter_args <- function(fun, args_list) {
  fml <- try(formals(fun), silent = TRUE)
  if (inherits(fml, "try-error")) return(list())
  if ("..." %in% names(fml)) return(args_list)
  args_list[names(args_list) %in% names(fml)]
}


#' Merge caller-supplied extras onto compare_models_cv()'s base arguments
#'
#' `c(base, extra)` would produce two entries of the same name whenever an
#' extra collides with a base argument, and `do.call()` then fails with
#' "formal argument ... matched by multiple actual arguments".  A user passing
#' `rf_args = list(seed = 3)` hits this immediately, so extras replace base
#' entries by name instead of being appended alongside them.
#'
#' The four arguments that define *what is being compared* are protected: a
#' per-model override of the data, the response, the predictors or the folds
#' would silently make the models incomparable, which is the whole point of the
#' function.
#'
#' @param base Named list of arguments built by `compare_models_cv()`.
#' @param extra Named list supplied by the caller (`gwr_args` etc.).
#' @param arg_label Name of the caller-facing argument, for messages.
#' @return A named list with no duplicate names.
#' @keywords internal
#' @noRd
.merge_args <- function(base, extra, arg_label) {
  if (is.null(extra) || length(extra) == 0L) return(base)
  if (is.null(names(extra)) || any(!nzchar(names(extra))))
    stop(sprintf("compare_models_cv(): every element of `%s` must be named.",
                 arg_label), call. = FALSE)

  protected <- c("data_sf", "response_var", "predictor_vars", "folds")
  clash <- intersect(names(extra), protected)
  if (length(clash) > 0L) {
    warning(sprintf(
      paste0("compare_models_cv(): ignoring %s in `%s` -- these are set by ",
             "compare_models_cv() itself so that every model is scored on the ",
             "same data and the same folds."),
      paste(sQuote(clash), collapse = ", "), arg_label), call. = FALSE)
    extra <- extra[setdiff(names(extra), protected)]
  }
  utils::modifyList(base, extra)
}


#' Check whether a model backend is available
#' @keywords internal
#' @noRd
.model_available <- function(model_name) {
  switch(model_name,
    "GWR" = requireNamespace("GWmodel", quietly = TRUE) &&
      requireNamespace("sp", quietly = TRUE),
    "Bayesian" = requireNamespace("brms", quietly = TRUE),
    "RF" = requireNamespace("ranger", quietly = TRUE),
    FALSE
  )
}


#' Warn about argument-list entries a CV wrapper will silently discard
#'
#' \code{.filter_args()} drops anything the target function has no formal for
#' (and no \code{...} to absorb it), which is correct but silent: a mistyped or
#' misplaced entry then simply has no effect.  This names the casualties.
#'
#' @param fun The function the arguments are destined for.
#' @param args_list The user-supplied argument list.
#' @param arg_name Name of the parameter it arrived as, for the message.
#' @param fun_name Name of \code{fun}, for the message.
#' @return \code{NULL}, invisibly.  Called for the warning.
#' @keywords internal
#' @noRd
.warn_dropped_args <- function(fun, args_list, arg_name, fun_name) {
  if (!length(args_list)) return(invisible(NULL))
  fml <- try(formals(fun), silent = TRUE)
  if (inherits(fml, "try-error") || "..." %in% names(fml))
    return(invisible(NULL))
  dropped <- setdiff(names(args_list), names(fml))
  if (length(dropped) == 0L) return(invisible(NULL))
  warning(sprintf(
    paste0("compare_models_cv(): ignoring %s entr%s %s -- %s() has no such ",
           "argument. Only arguments of %s() are forwarded from `%s`; pass ",
           "anything else by calling %s() yourself."),
    arg_name, if (length(dropped) == 1L) "y" else "ies",
    paste(sQuote(dropped), collapse = ", "), fun_name, fun_name,
    arg_name, fun_name), call. = FALSE)
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Post-fit residual spatial autocorrelation check (Moran's I)
# ---------------------------------------------------------------------------

#' Index pairs of the k nearest *other* observations
#'
#' Returns the \code{(i, j)} pairs of a k-nearest-neighbour graph with every
#' self-match removed, so the weight matrix built from them has a zero
#' diagonal.  That is not cosmetic: Moran's I is defined only for
#' \eqn{w_{ii} = 0}, and both \eqn{E[I] = -1/(n-1)} and the Cliff & Ord
#' variance assume it.
#'
#' \code{FNN::get.knn()} reports a point's OWN index among its neighbours
#' whenever exact duplicate coordinates are present, which put \eqn{1/k} on the
#' diagonal and added a strictly positive \eqn{w_{ii} e_i^2} term to the
#' numerator.  Repeat observations at one site are exactly what
#' \code{make_folds(method = "leave_location_out")} is for, and
#' \code{residual_morans_i()} reads \code{fit$data_sf} without de-duplicating,
#' so this was a mainstream input rather than a corner case: with 40 sites x 4
#' repeats and a response carrying no spatial structure at all, 120 of 160 rows
#' gained a self-weight, mean I came out at +0.086 against
#' \eqn{E[I] = -0.0063}, and 77\% of samples were "significant" at
#' \eqn{p < 0.05} (nominal 5\%).  \code{spdep::knearneigh()} never returns a
#' self-match, and neither did the dense fallback, so the statistic also
#' depended silently on whether \pkg{FNN} happened to be installed.
#'
#' The kd-tree path therefore asks for \code{k + 1} neighbours, drops any
#' \code{j == i}, and keeps the \code{k} nearest of what remains.
#'
#' @param coords Numeric matrix (n x 2) of projected coordinates.
#' @param k Integer neighbour count, already clamped to \code{n - 1}.
#' @param use_fnn Logical; use \pkg{FNN}'s kd-tree rather than a dense
#'   \code{dist()} scan.
#' @return A list with integer vectors \code{i} and \code{j} of equal length.
#' @keywords internal
#' @noRd
.knn_pairs <- function(coords, k, use_fnn) {
  n <- nrow(coords)

  # At k >= n - 1 every point neighbours every other, so no lookup can add
  # anything -- and asking FNN for n - 1 neighbours of duplicated points can
  # spend one of them on a self-match, silently losing a genuine neighbour.
  if (k >= n - 1L) {
    return(list(
      i = rep(seq_len(n), each = n - 1L),
      j = as.integer(unlist(lapply(seq_len(n), function(i) seq_len(n)[-i]),
                            use.names = FALSE))
    ))
  }

  if (isTRUE(use_fnn)) {
    kq     <- k + 1L                                # k + 1 <= n - 1 here
    nn_idx <- FNN::get.knn(coords, k = kq)$nn.index  # n x kq, nearest first
    # Entries within a row are distinct point indices, so at most one of them
    # can be i itself.
    self <- nn_idx == seq_len(n)     # recycles down columns: [i, j] against i

    if (!any(self)) {
      # No slot was spent on a self-match, so all k + 1 returned neighbours are
      # genuine and the k nearest of them are a correct k-NN set.
      idx <- nn_idx[, seq_len(k), drop = FALSE]
      return(list(i = rep(seq_len(n), each = k), j = as.integer(t(idx))))
    }

    # A self-match means exact duplicate coordinates, and dropping it is NOT
    # enough: the slot it occupied displaced a genuine tied neighbour, so the
    # k that remain are not the k nearest.  Measured on 25 sites x 4 repeats,
    # k = 3: 75 of 400 retained pairs were a point at distance 121 standing in
    # for a co-located one at distance 0.  Group the duplicates and answer
    # exactly instead.
    return(.knn_pairs_dup(coords, k))
  }

  dmat <- as.matrix(stats::dist(coords))
  diag(dmat) <- Inf                                 # self is never a neighbour
  # matrix(..., nrow = n, ncol = k) forces the shape apply() will not.
  # At k = 1 the inner function returns a scalar, so apply() simplifies to a
  # length-n VECTOR and t() turns it into a 1 x n matrix -- making nn_idx[i, ]
  # fail with "subscript out of bounds" for every i > 1.
  # residual_morans_i(fit, k = 1) reached this on any machine without FNN.
  nn_idx <- matrix(t(apply(dmat, 1, function(row) order(row)[seq_len(k)])),
                   nrow = n, ncol = k)
  list(i = rep(seq_len(n), each = k), j = as.integer(t(nn_idx)))
}


#' Exact k-nearest-neighbour pairs when the coordinates contain exact duplicates
#'
#' \code{FNN::get.knn()} answers a tied query by returning \emph{some} of the
#' tied points, and the point's own index is eligible to be one of them --- so
#' with duplicates a \code{k + 1} query can come back holding self \emph{and}
#' having dropped a genuine co-located neighbour, leaving a farther point in
#' its place.  Requesting one extra neighbour therefore removes the self-weight
#' but does not restore the neighbour it displaced.
#'
#' Duplicates are not exotic here: repeat observations at one site are exactly
#' what \code{make_folds(method = "leave_location_out")} exists for, and
#' \code{residual_morans_i()} reads \code{fit$data_sf} without de-duplicating.
#'
#' The structure of the problem makes an exact answer cheap.  Points sharing a
#' coordinate are at distance 0 from each other and at an identical distance
#' from everything else, so the neighbour set is determined group-wise: take
#' the other members of the point's own group first, then fill from the nearest
#' \emph{other} groups in order.  The k-d tree runs on the group
#' representatives, which are distinct by construction and so cannot tie
#' against themselves.
#'
#' @param coords Numeric matrix (n x 2) of projected coordinates.
#' @param k Integer neighbour count, already known to be \code{< n - 1}.
#' @return A list with integer vectors \code{i} and \code{j}, \code{k} entries
#'   per row of \code{coords}.
#' @keywords internal
#' @noRd
.knn_pairs_dup <- function(coords, k) {
  n <- nrow(coords)

  # "%.17g" round-trips a double exactly, so two rows share a key iff they are
  # bit-identical.  (as.character() would stop at 15 significant digits and
  # could merge two points a few ulps apart -- harmless, but only by accident.)
  # `+ 0` normalises a negative zero: sprintf("%.17g", -0) is "-0", which
  # would put a point at (-0, y) in a different group from one at (0, y)
  # although they are the same location and st_distance() says 0.
  key <- paste(sprintf("%.17g", coords[, 1L] + 0),
               sprintf("%.17g", coords[, 2L] + 0), sep = "\r")
  gid     <- match(key, key)              # representative row index per point
  reps    <- which(!duplicated(gid))      # one row index per distinct location
  gid     <- match(gid, gid[reps])        # 1..G
  members <- split(seq_len(n), gid)
  G       <- length(reps)

  # Every point at one location: k < n - 1 already, so the first k of the other
  # members is a correct answer (they are all at distance 0).
  if (G == 1L) {
    j <- unlist(lapply(seq_len(n), function(i) seq_len(n)[-i][seq_len(k)]),
                use.names = FALSE)
    return(list(i = rep(seq_len(n), each = k), j = as.integer(j)))
  }

  # Nearest other GROUPS, in distance order.  k groups always suffice: each
  # supplies at least one member and at most k are ever needed.  Ask for one
  # extra and filter, so a self-match here could only cost an unused slot.
  kg  <- min(k + 1L, G - 1L)
  nn  <- FNN::get.knn(coords[reps, , drop = FALSE], k = kg)$nn.index
  nn  <- matrix(nn, nrow = G, ncol = kg)

  out_i <- vector("list", G)
  out_j <- vector("list", G)
  for (g in seq_len(G)) {
    mem <- members[[g]]
    m   <- length(mem)

    # Members of the nearest other groups, in order, enough to top up any point
    # of this group.  Identical for every point in the group, since they share
    # a coordinate.
    ext <- integer(0)
    if (k > m - 1L) {
      for (h in nn[g, ]) {
        if (h == g) next                  # defensive; representatives are distinct
        ext <- c(ext, members[[h]])
        if (length(ext) >= k - (m - 1L)) break
      }
      ext <- ext[seq_len(k - (m - 1L))]
    }

    js <- lapply(mem, function(i) {
      own <- mem[mem != i]
      if (length(own) >= k) own[seq_len(k)] else c(own, ext)
    })
    out_i[[g]] <- rep(mem, each = k)
    out_j[[g]] <- unlist(js, use.names = FALSE)
  }

  list(i = as.integer(unlist(out_i, use.names = FALSE)),
       j = as.integer(unlist(out_j, use.names = FALSE)))
}


#' Build a row-standardised k-nearest-neighbour sparse weight matrix
#'
#' For each observation the \code{k} closest \emph{other} observations receive
#' weight 1; all other pairs, and the diagonal, receive weight 0.
#' The resulting matrix is then row-standardised so that each row sums to 1.
#' This is the standard default in spatial statistics (Anselin, 1988) and is
#' far more robust to irregularly-spaced or clustered data than an
#' inverse-distance scheme, which gives enormous weight to very close pairs
#' and can inflate Moran's I significance.
#'
#' The neighbour lookup itself is \code{.knn_pairs()}, which is where the
#' zero-diagonal guarantee lives — see the note there on duplicate
#' coordinates.
#'
#' Uses \pkg{FNN} for O(n·k) kd-tree nearest-neighbour lookup when available,
#' avoiding the O(n²) full distance matrix.  Returns a
#' \code{Matrix::sparseMatrix} (dgCMatrix), which keeps memory proportional to
#' O(n·k) instead of O(n²).  Falls back to the dense brute-force path when
#' \pkg{FNN} or \pkg{Matrix} are not installed.
#'
#' @param coords Numeric matrix (n x 2) of projected coordinates.
#' @param k Integer number of nearest neighbours (default 8).
#' @param use_fnn Logical; use \pkg{FNN}'s kd-tree for the neighbour lookup.
#'   Defaults to whether \pkg{FNN} is installed.  Exposed so the dense
#'   brute-force fallback is reachable on a machine that has FNN.
#' @param use_matrix Logical; return a \pkg{Matrix} sparse matrix rather than a
#'   dense base matrix.  Defaults to whether \pkg{Matrix} is installed.
#'   Exposed for the same reason as \code{use_fnn}.  The fast sparse path
#'   requires both backends; with either disabled the dense path runs.
#' @return A row-standardised weight matrix W — sparse (dgCMatrix) when
#'   \pkg{Matrix} is available, otherwise a dense base matrix.
#' @keywords internal
#' @noRd
.build_knn_weights <- function(coords, k = 8L,
                              use_fnn    = requireNamespace("FNN",    quietly = TRUE),
                              use_matrix = requireNamespace("Matrix", quietly = TRUE)) {
  n <- nrow(coords)
  k <- min(as.integer(k), n - 1L)
  if (k < 1L) k <- 1L

  # Backend availability is a parameter rather than a direct requireNamespace()
  # call so the dense fallback -- and its size guard -- can be exercised on a
  # machine that has FNN installed.  Defaults preserve the previous behaviour.
  has_fnn    <- isTRUE(use_fnn)
  has_matrix <- isTRUE(use_matrix)

  # --- Fallback: dense O(n²) path when packages are missing ----
  # The guard has to test BOTH backends, not just FNN: the dense branch is
  # entered whenever either is missing and always allocates W <- matrix(0,n,n),
  # so keying it on FNN alone let an unbounded n x n allocation through whenever
  # FNN was present but Matrix was not (or use_matrix = FALSE).
  if (!(has_fnn && has_matrix) && n > 5000L)
    stop("n = ", n, " requires FNN for k-NN weights, and Matrix to hold them sparsely (the dense fallback would allocate an n*n matrix). Install both with install.packages(c(\"FNN\", \"Matrix\")).", call. = FALSE)

  # Both backends share one neighbour lookup, so the zero diagonal (and the
  # k + 1 request that FNN needs to guarantee it) cannot drift apart between
  # them.  The two paths differ only in how W is stored.
  pr    <- .knn_pairs(coords, k = k, use_fnn = has_fnn)
  row_i <- pr$i
  col_j <- pr$j
  deg   <- tabulate(row_i, nbins = n)   # neighbours actually kept per row
  deg[deg == 0L] <- 1L                  # isolate: leave a zero row, not NaN

  if (has_fnn && has_matrix) {
    # --- Fast path: O(n*k) kd-tree lookup + sparse matrix ----
    W <- Matrix::sparseMatrix(
      i = row_i, j = col_j, x = 1 / deg[row_i],             # row-standardised
      dims = c(n, n), repr = "C"
    )
  } else {
    W <- matrix(0, n, n)
    W[cbind(row_i, col_j)] <- 1 / deg[row_i]
  }

  W
}


#' Rebuild a fitted model's design matrix and test whether its residuals are
#' the OLS residuals on it
#'
#' The Cliff & Ord residual moments are exact for \eqn{e = My} with
#' \eqn{M = I - X(X'X)^{-1}X'} and nothing else, so rather than guessing from
#' the fit's class whether that holds, this checks it directly: rebuild
#' \code{X} from \code{predictor_vars} and \code{data_sf}, regress the response
#' on it, and compare the result with the residuals actually supplied.  The
#' answer is exactly the condition the formula needs, and it degrades
#' gracefully — a GWR whose bandwidth is wide enough to be global OLS passes,
#' the same GWR at a small bandwidth does not.
#'
#' @param fit A \code{spatial_fit}.
#' @param resid The residual vector already extracted and subset by \code{keep}.
#' @param keep Logical vector over the rows of \code{fit$data_sf} marking the
#'   observations that survived the finite-residual filter.
#' @return \code{NULL} when the design cannot be rebuilt, otherwise a list with
#'   \code{X} (the model matrix, rows subset by \code{keep}) and \code{is_ols}.
#' @keywords internal
#' @noRd
.morans_ols_design <- function(fit, resid, keep) {
  d  <- fit$data_sf
  rv <- fit$response_var
  pv <- fit$predictor_vars
  if (!inherits(d, "sf")) return(NULL)
  if (!is.character(rv) || length(rv) != 1L || is.na(rv)) return(NULL)
  if (is.null(pv)) pv <- character(0)
  if (!is.character(pv)) return(NULL)

  df <- tryCatch(sf::st_drop_geometry(d), error = function(e) NULL)
  if (is.null(df) || nrow(df) != length(keep)) return(NULL)
  if (!all(c(rv, pv) %in% names(df))) return(NULL)

  df <- df[keep, c(rv, pv), drop = FALSE]
  if (nrow(df) != length(resid) || anyNA(df)) return(NULL)

  y <- df[[rv]]
  if (!is.numeric(y) || !all(is.finite(y))) return(NULL)

  X <- if (length(pv) == 0L) {
    # An intercept-only model.  tr(MW) = -S0/n for any W with a zero diagonal,
    # so the residual moments reduce to E[I] = -1/(n - 1) exactly -- the
    # classical value -- which is a useful sanity check on the general path.
    matrix(1, nrow(df), 1L, dimnames = list(NULL, "(Intercept)"))
  } else {
    tryCatch(
      stats::model.matrix(~ ., data = df[, pv, drop = FALSE]),
      error = function(e) NULL
    )
  }
  if (!is.matrix(X) || nrow(X) != nrow(df) || !all(is.finite(X))) return(NULL)

  e_ols <- tryCatch(qr.resid(qr(X), y), error = function(e) NULL)
  if (is.null(e_ols) || length(e_ols) != length(resid)) return(NULL)

  scale <- max(stats::sd(y), .Machine$double.eps)
  list(X = X, is_ols = max(abs(resid - e_ols)) <= 1e-7 * scale)
}


#' Cliff & Ord (1981) sec. 8.3 moments of Moran's I for regression residuals
#'
#' For \eqn{e = M\varepsilon} with \eqn{M = I - X(X'X)^{-1}X'} and
#' \eqn{\varepsilon \sim N(0, \sigma^2 I)}:
#' \deqn{E[I] = (n/S_0)\,\mathrm{tr}(MW)/(n-p)}
#' \deqn{Var[I] = (n/S_0)^2 [\mathrm{tr}(MWMW') + \mathrm{tr}((MW)^2) +
#'   (\mathrm{tr}MW)^2]/[(n-p)(n-p+2)] - E[I]^2}
#'
#' \eqn{M} is never formed.  With \eqn{Q} an orthonormal basis of
#' \eqn{col(X)} (so \eqn{P = QQ'}), every trace reduces to a \eqn{q \times q}
#' one plus Frobenius norms of the \eqn{n \times q} products \eqn{WQ} and
#' \eqn{W'Q}, which keeps the sparse weight matrix sparse:
#' \deqn{\mathrm{tr}(MW) = \mathrm{tr}(W) - \mathrm{tr}(Q'WQ)}
#' \deqn{\mathrm{tr}((MW)^2) = \mathrm{tr}(W^2) - 2\,\mathrm{tr}(PW^2) +
#'   \mathrm{tr}((Q'WQ)^2)}
#' \deqn{\mathrm{tr}(MWMW') = \mathrm{tr}(WW') - \mathrm{tr}(PW'W) -
#'   \mathrm{tr}(PWW') + \mathrm{tr}((Q'WQ)(Q'WQ)')}
#'
#' Verified to agree with \code{spdep::lm.morantest()} and with an explicit
#' dense \eqn{M} to machine precision.
#'
#' @param W The weight matrix (base or \pkg{Matrix}).
#' @param X The design matrix.
#' @param S0 \code{sum(W)}, already computed by the caller.
#' @param is_sparse Whether \code{W} is a \pkg{Matrix} object.
#' @return \code{NULL} when the residual degrees of freedom are too small to
#'   support the formula, otherwise a list with \code{EI}, \code{VI},
#'   \code{df} and \code{p}.
#' @keywords internal
#' @noRd
.morans_residual_moments <- function(W, X, S0, is_sparse) {
  n   <- nrow(X)
  qrX <- tryCatch(qr(X), error = function(e) NULL)
  if (is.null(qrX)) return(NULL)
  q  <- qrX$rank
  df <- n - q
  # (n - p + 2) sits in the denominator and the normal approximation is
  # meaningless with a handful of residual degrees of freedom either way.
  if (!is.finite(q) || q < 1L || df < 4L) return(NULL)

  Q  <- qr.Q(qrX)[, seq_len(q), drop = FALSE]
  Wt <- if (is_sparse) Matrix::t(W) else t(W)
  WQ <- as.matrix(W  %*% Q)                     # n x q
  HQ <- as.matrix(Wt %*% Q)                     # n x q  (= W' Q)
  G  <- crossprod(Q, WQ)                        # q x q  (= Q' W Q)

  trW   <- sum(if (is_sparse) Matrix::diag(W) else diag(W))
  trMW  <- trW - sum(diag(G))
  trWWt <- sum(W * W)                           # tr(W W')
  trW2  <- sum(W * Wt)                          # tr(W^2)

  trMWMWt <- trWWt - sum(WQ * WQ) - sum(HQ * HQ) + sum(G * G)
  trMWMW  <- trW2  - 2 * sum(HQ * WQ)           + sum(G * t(G))

  EI <- (n / S0) * trMW / df
  VI <- (n / S0)^2 * (trMWMWt + trMWMW + trMW^2) / (df * (df + 2)) - EI^2
  if (!is.finite(EI) || !is.finite(VI)) return(NULL)
  list(EI = EI, VI = VI, df = df, p = q)
}


#' Compute Moran's I on the residuals of a fitted spatial model
#'
#' Given a \code{spatial_fit} object, extracts the residuals and the
#' observation coordinates, builds a spatial weight matrix, and computes
#' Moran's I together with its analytical expectation and variance under a
#' stated null.  A z-score and two-sided p-value are provided so the caller can
#' assess whether statistically significant spatial autocorrelation remains
#' after fitting.
#'
#' By default, weights are constructed as a k-nearest-neighbour (k = 8)
#' binary matrix, row-standardised.  Users may supply their own weight
#' matrix via the \code{weights} argument.
#'
#' @section Which null, and when it is approximate:
#' Two nulls are available, and the one actually used is reported back in the
#' \code{null} element of the return value.
#'
#' \code{"randomisation"} is the classical exchangeable null:
#' \eqn{E[I] = -1/(n-1)} with the Cliff & Ord randomisation variance,
#' conditioning on the observed kurtosis.  These are the moments of I for a
#' vector whose elements are equally likely in any order.
#'
#' \strong{Model residuals are not exchangeable.}  They are orthogonal to the
#' design matrix, which pushes \eqn{E[I]} materially below \eqn{-1/(n-1)}
#' whenever the covariates are spatially smooth — and pushes it further the
#' more covariates there are.  In a simulation with \eqn{n = 120}, six smooth
#' covariates and \emph{independent} errors (so the truth is "no residual
#' autocorrelation"), OLS residuals had mean \eqn{I = -0.031} against the
#' exchangeable \eqn{E[I] = -0.008}; the z-score averaged \eqn{-0.54} with
#' \eqn{sd = 0.90} instead of 0 and 1.  The cost is power, which is the point
#' of the test: at a moderate residual autocorrelation the exchangeable null
#' rejected 13\% of the time where the correct one rejected 31\%.
#'
#' \code{"residual"} therefore uses the Cliff & Ord (1981) sec. 8.3 moments for
#' regression residuals, with \eqn{M = I - X(X'X)^{-1}X'} rebuilt from
#' \code{predictor_vars} and \code{data_sf}:
#' \deqn{E[I] = (n/S_0)\,\mathrm{tr}(MW)/(n-p)}
#' \deqn{Var[I] = (n/S_0)^2[\mathrm{tr}(MWMW') + \mathrm{tr}((MW)^2) +
#'   (\mathrm{tr}MW)^2]/[(n-p)(n-p+2)] - E[I]^2}
#' These assume normal errors rather than conditioning on the observed
#' kurtosis.  On the simulation above they restored the z-score to mean
#' \eqn{-0.09}, \eqn{sd = 1.03}, and the rejection rate to 4.3\% against a
#' nominal 5\%.  They agree with \code{spdep::lm.morantest()} to machine
#' precision.
#'
#' \strong{These moments are exact for \eqn{e = My} and for nothing else}, so
#' \code{null = "auto"} does not guess from the fit's class: it rebuilds
#' \code{X}, regresses the response on it, and uses the residual moments only
#' when the supplied residuals \emph{are} those OLS residuals to numerical
#' tolerance.  A GWR wide enough to have collapsed to global OLS passes that
#' test; the same GWR at a working bandwidth does not.
#'
#' \strong{For the flexible backends neither null is exact}, and \code{"auto"}
#' leaves them on \code{"randomisation"} because forcing the OLS moments on
#' them measurably makes matters worse, not better.  Measured on null data
#' (\eqn{n = 120}, three smooth covariates, independent errors; nominal 5\%,
#' one-sided):
#' \tabular{lrr}{
#'   \strong{backend} \tab \strong{randomisation} \tab \strong{residual} \cr
#'   OLS              \tab 0.035 \tab 0.060 \cr
#'   random forest    \tab 0.128 \tab 0.200 \cr
#'   GWR              \tab 0.000 \tab 0.000
#' }
#' An in-sample random forest is anticonservative under both — its residuals
#' are shrunk and spatially heteroscedastic, so the variance is understated
#' whichever moments are used (\eqn{sd(z) \approx 1.3}) — and GWR is
#' conservative under both, because it removes far more structure than a rank-p
#' projection does.  Treat the p-value from those backends as a rough
#' indicator, and prefer cross-validated residuals or an explicit spatial
#' covariance model when the answer has to carry weight.
#'
#' A permutation null was considered and rejected: permuting the residual
#' vector destroys exactly the orthogonality that causes the bias, so its mean
#' is the exchangeable \eqn{-1/(n-1)} by construction (measured:
#' \eqn{-0.00840} against \eqn{-1/(n-1) = -0.00840}) and it reproduces the
#' randomisation null rather than correcting it.
#'
#' @param fit A \code{spatial_fit} object (from \code{fit_gwr_model},
#'   \code{fit_bayesian_spatial_model} or \code{fit_rf_model}).
#' @param alternative Character: \code{"two.sided"} (default),
#'   \code{"greater"} (positive autocorrelation), or \code{"less"}.
#' @param weights Optional user-supplied n x n weight matrix — a base
#'   matrix or a \pkg{Matrix}-package matrix (e.g. a sparse dgCMatrix).
#'   When \code{NULL} (the default), a k-nearest-neighbour binary weight matrix
#'   (k = 8, row-standardised) is built from the observation coordinates.
#'   If a non-row-standardised matrix is supplied (i.e. rows do not all sum
#'   to 1), the Cliff & Ord variance formula is still valid for general W
#'   and the computation proceeds; a note is logged (not raised as an R
#'   warning) because the magnitude of I is not directly comparable to
#'   results obtained with row-standardised weights.  A \code{weights}
#'   argument that is not an n x n matrix is an error.
#' @param k Integer number of nearest neighbours used when building the
#'   default weight matrix (ignored when \code{weights} is supplied).
#'   Default 8.
#' @param null Which null distribution the expectation, variance and p-value
#'   are computed against.  One of:
#'   \describe{
#'     \item{\code{"auto"} (default)}{\code{"residual"} when the design matrix
#'       can be rebuilt \emph{and} the fit's residuals are the OLS residuals on
#'       it, \code{"randomisation"} otherwise.}
#'     \item{\code{"randomisation"}}{Always the exchangeable moments.}
#'     \item{\code{"residual"}}{Always the Cliff & Ord residual moments.  Falls
#'       back to \code{"randomisation"} with a logged warning if the design
#'       cannot be rebuilt, and warns (but proceeds) if the residuals are not
#'       the OLS residuals on it, in which case the moments are approximate.}
#'   }
#'   See \strong{Which null, and when it is approximate} above.
#' @return A list with components:
#'   \describe{
#'     \item{observed}{Numeric scalar, Moran's I statistic.}
#'     \item{expected}{Expected I under the null named by \code{null}:
#'       \eqn{-1/(n-1)} for \code{"randomisation"},
#'       \eqn{(n/S_0)\mathrm{tr}(MW)/(n-p)} for \code{"residual"}.}
#'     \item{sd}{Standard deviation of I under that same null.}
#'     \item{z}{Standardised z-score, \eqn{(I - E[I]) / sd(I)}.}
#'     \item{p_value}{Two-sided (or one-sided) p-value from the
#'       normal approximation.}
#'     \item{n}{Number of observations used.}
#'     \item{null}{The null actually used, \code{"randomisation"} or
#'       \code{"residual"} — check this rather than assuming, since
#'       \code{"auto"} chooses per fit and \code{"residual"} can fall back.}
#'     \item{df}{Residual degrees of freedom behind the moments:
#'       \eqn{n - p} for \code{"residual"} (where \eqn{p} is the rank of the
#'       design matrix), \eqn{n - 1} for \code{"randomisation"}.}
#'   }
#'   Returns \code{NULL} with a warning if computation fails (e.g. fewer
#'   than 4 valid residuals).
#' @references Cliff, A. D. and Ord, J. K. (1981) \emph{Spatial Processes:
#'   Models and Applications}. Pion, London. Section 8.3.
#' @family model evaluation
#' @examples
#' \donttest{
#' # Works on any spatial_fit; a forest keeps the example free of the optional
#' # GWR/Stan backends.
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   # A strong east-west trend the predictor cannot explain: the residuals
#'   # should still carry spatial structure, and this is what detects it.
#'   dat$price <- 10 + 0.02 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   fit <- fit_rf_model(dat, "price", "elev", num_trees = 100, seed = 1)
#'   residual_morans_i(fit)  # z near 0 / p large = no residual structure
#' }
#' }
#' @export
residual_morans_i <- function(fit,
                              alternative = c("two.sided", "greater", "less"),
                              weights = NULL,
                              k = 8L,
                              null = c("auto", "randomisation", "residual")) {
  alternative <- match.arg(alternative)
  null        <- match.arg(null)

  if (!inherits(fit, "spatial_fit")) {
    .warn_and_log("residual_morans_i(): `fit` is not a spatial_fit object.")
    return(NULL)
  }

  # --- Extract residuals & coordinates ---
  resid <- tryCatch(residuals(fit), error = function(e) NULL)
  if (is.null(resid) || length(resid) < 4L) {
    .warn_and_log("residual_morans_i(): could not extract enough residuals (n < 4).")
    return(NULL)
  }

  coords <- tryCatch({
    pts <- ensure_projected(coerce_to_points(fit$data_sf, "auto"))
    sf::st_coordinates(pts)[, 1:2, drop = FALSE]
  }, error = function(e) NULL)

  if (is.null(coords) || nrow(coords) != length(resid)) {
    .warn_and_log("residual_morans_i(): coordinate extraction failed.")
    return(NULL)
  }

  # Drop any non-finite residuals
  ok <- is.finite(resid)
  if (sum(ok) < 4L) {
    .warn_and_log("residual_morans_i(): fewer than 4 finite residuals.")
    return(NULL)
  }
  resid  <- resid[ok]
  coords <- coords[ok, , drop = FALSE]
  n      <- length(resid)

  # --- Weight matrix ---
  if (!is.null(weights)) {
    # Accept base matrices AND Matrix-package sparse/dense matrices
    # (is.matrix() is FALSE for e.g. dgCMatrix, which previously caused a
    # misleading "wrong dimensions" warning and a silent kNN fallback).
    # A data.frame of weights is accepted as a matrix; anything else that is
    # not an n x n matrix is an ERROR.  It used to be a logged line and a
    # silent fall-back to the default kNN(8) matrix -- so a user who supplied
    # their own weights with one row too few, or as a data.frame, got the
    # statistic for a weight matrix they never chose, with no R condition
    # raised (measured: I = 0.874 returned for four malformed shapes, against
    # 0.805 for the weights actually supplied).
    if (is.data.frame(weights)) weights <- as.matrix(weights)
    is_mat_like <- is.matrix(weights) || inherits(weights, "Matrix")
    if (!is_mat_like || nrow(weights) != n || ncol(weights) != n) {
      stop(sprintf(paste0("residual_morans_i(): `weights` must be an n x n matrix ",
                          "(base matrix or Matrix), with n = %d the number of ",
                          "finite residuals; got %s of dimension %s. Omit ",
                          "`weights` to use the default k-nearest-neighbour ",
                          "matrix."),
                   n, paste(class(weights), collapse = "/"),
                   if (is_mat_like) paste(dim(weights), collapse = " x ")
                   else paste0("length ", length(weights))),
           call. = FALSE)
    } else {
      W <- weights
      # Validate row-standardisation: each row should sum to 1 (or 0 for
      # isolates).  If not, the Moran's I statistic is still computed
      # correctly for a general W via the Cliff & Ord formula, but the
      # interpretation differs from the row-standardised case.
      # (rowSums() does not S4-dispatch for Matrix classes from package
      # code when Matrix is loaded but not attached — same as the
      # variance block below.)
      rs <- if (inherits(W, "Matrix")) Matrix::rowSums(W) else rowSums(W)
      non_zero <- rs[rs != 0]
      if (length(non_zero) > 0 && !isTRUE(all.equal(non_zero, rep(1, length(non_zero)), tolerance = 1e-8))) {
        .log_warn(
          paste0("residual_morans_i(): user-supplied `weights` is not row-standardised ",
                 "(row sums range from %.4g to %.4g). The Moran's I statistic and ",
                 "variance are computed using the general Cliff & Ord formula, which ",
                 "remains valid, but the magnitude of I is not directly comparable to ",
                 "results from row-standardised weights."),
          min(non_zero), max(non_zero)
        )
      }
    }
  } else {
    W <- .build_knn_weights(coords, k = k)
  }
  S0 <- sum(W)

  # Moran's I is a function of the CENTRED residuals, so the degeneracy guard
  # has to test their sum of squares rather than the raw one.  Constant
  # non-zero residuals pass sum(resid^2) but give resid_c == 0, and the
  # variance block below then computes m2 = 0 -> b2 = NaN -> VI = NaN, so
  # `if (VI > 0)` raises "missing value where TRUE/FALSE needed" instead of
  # this function returning NULL.
  resid_c <- resid - mean(resid)
  ss_c    <- sum(resid_c^2)
  if (S0 < .Machine$double.eps || ss_c < .Machine$double.eps) {
    .warn_and_log("residual_morans_i(): degenerate weights or zero residual variance.")
    return(NULL)
  }

  # --- Moran's I ---
  # sum(resid_c * (W %*% resid_c)) rather than crossprod(): see the dispatch
  # note in the variance block below.  The two are numerically identical.
  I <- (n / S0) * sum(resid_c * (W %*% resid_c)) / ss_c

  # NOTE for everything below: t(), rowSums(), colSums() and crossprod() are
  # plain base functions that do NOT dispatch to Matrix's S4 methods when
  # called from package code with Matrix loaded-but-not-attached, so sparse W
  # needs the Matrix:: generics explicitly (or, for crossprod, the
  # primitive-only rewrite used for I above).  Only primitives -- *, %*%, and
  # sum() -- dispatch on their own; crossprod is neither a primitive nor an
  # internal generic, which is why it belongs on this list and not with them.
  is_sparse <- inherits(W, "Matrix")

  # --- Choose the null -------------------------------------------------------
  # The exchangeable moments below are the moments of I for a vector whose
  # elements are equally likely in any order.  Model residuals are not such a
  # vector: they are orthogonal to the design matrix, which drags the true
  # E[I] well below -1/(n - 1) once the covariates are spatially smooth.  See
  # ?residual_morans_i, section "Which null, and when it is approximate", for
  # the measured size and power cost and for why "auto" refuses to apply the
  # OLS residual moments to a backend whose residuals are not OLS residuals.
  des <- if (identical(null, "randomisation")) NULL else
    .morans_ols_design(fit, resid, ok)

  use_residual <- switch(null,
    randomisation = FALSE,
    auto          = !is.null(des) && isTRUE(des$is_ols),
    residual      = !is.null(des)
  )
  if (identical(null, "residual")) {
    if (is.null(des)) {
      .log_warn(paste0("residual_morans_i(): null = \"residual\" needs the design ",
                       "matrix, which cannot be rebuilt from this fit's ",
                       "`predictor_vars` and `data_sf`; using the randomisation ",
                       "null instead. The returned `null` element says which was used."))
    } else if (!isTRUE(des$is_ols)) {
      .log_warn(paste0("residual_morans_i(): null = \"residual\" was requested, but ",
                       "these residuals are not the OLS residuals of the response on ",
                       "the rebuilt design matrix, so the Cliff & Ord residual ",
                       "moments are an approximation here rather than exact."))
    }
  }

  mom <- NULL
  if (use_residual) {
    mom <- .morans_residual_moments(W = W, X = des$X, S0 = S0,
                                    is_sparse = is_sparse)
    if (is.null(mom) && !identical(null, "auto"))
      .log_warn(paste0("residual_morans_i(): the residual moments are not usable ",
                       "here (too few residual degrees of freedom); using the ",
                       "randomisation null instead."))
  }

  if (!is.null(mom)) {
    # --- Cliff & Ord (1981) sec. 8.3 residual moments ---
    null_used <- "residual"
    EI   <- mom$EI
    VI   <- mom$VI
    df_I <- mom$df
  } else {
    # --- Analytical expectation & variance (randomisation assumption) ---
    null_used <- "randomisation"
    EI   <- -1 / (n - 1)
    df_I <- n - 1

    # Cliff & Ord variance under randomisation
    # S1 = 0.5 * sum((W + t(W))^2) rewritten via the identity
    #    = sum(W^2) + sum(W * t(W))
    # so that sparse W never materialises the denser (W + t(W)) intermediate.
    Wt <- if (is_sparse) Matrix::t(W) else t(W)
    S1 <- sum(W * W) + sum(W * Wt)
    rs <- if (is_sparse) Matrix::rowSums(W) else rowSums(W)
    cs <- if (is_sparse) Matrix::colSums(W) else colSums(W)
    S2 <- sum((rs + cs)^2)
    m2 <- ss_c / n
    m4 <- sum(resid_c^4) / n
    b2 <- m4 / (m2^2)                          # kurtosis

    A  <- n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * S0^2)
    D  <- (n - 1) * (n - 2) * (n - 3) * S0^2
    C  <- (n^2 - n) * S1 - 2 * n * S2 + 6 * S0^2
    VI <- (A - b2 * C) / D - EI^2
  }

  # is.finite() as well as > 0: a degenerate kurtosis (b2) can make VI NaN,
  # and `if (NaN > 0)` is an error rather than FALSE.
  sd_I <- if (is.finite(VI) && VI > 0) sqrt(VI) else NA_real_
  z    <- if (is.finite(sd_I) && sd_I > 0) (I - EI) / sd_I else NA_real_

  p <- if (is.finite(z)) {
    switch(alternative,
      two.sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      greater   = stats::pnorm(z, lower.tail = FALSE),
      less      = stats::pnorm(z, lower.tail = TRUE)
    )
  } else {
    NA_real_
  }

  list(observed = I, expected = EI, sd = sd_I, z = z,
       p_value = p, n = n, null = null_used, df = df_I)
}


#' Compute residual Moran's I for every model in a named list of fits
#'
#' Convenience wrapper that calls \code{residual_morans_i()} on each
#' element and returns a data.frame suitable for joining to the metrics
#' table produced by \code{compare_models()}.
#'
#' @param fits Named list of \code{spatial_fit} objects.
#' @return A data.frame with columns \code{model}, \code{resid_morans_I},
#'   \code{resid_morans_z}, \code{resid_morans_p} and \code{resid_morans_null}
#'   (which null the p-value was computed against, per model).
#' @keywords internal
#' @noRd
.residual_morans_table <- function(fits) {
  # seq_along() rather than names(): fits[[nm]] returns the FIRST element of
  # that name, so duplicated names scored one fit repeatedly.  evaluate_insample()
  # now rejects duplicates outright, but positional indexing means this cannot
  # silently mis-report even if it is called directly.
  nms  <- names(fits)
  rows <- lapply(seq_along(fits), function(i) {
    mi <- residual_morans_i(fits[[i]])
    if (is.null(mi)) {
      data.frame(model = nms[i], resid_morans_I = NA_real_,
                 resid_morans_z = NA_real_, resid_morans_p = NA_real_,
                 resid_morans_null = NA_character_,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(model = nms[i], resid_morans_I = mi$observed,
                 resid_morans_z = mi$z, resid_morans_p = mi$p_value,
                 resid_morans_null = mi$null,
                 stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}


# ---------------------------------------------------------------------------
# evaluate_insample: metrics from already-fit spatial_fit objects
# ---------------------------------------------------------------------------

#' Compute in-sample (or out-of-sample) metrics for fitted spatial models
#'
#' Accepts a single \code{spatial_fit} object or a named list of them.
#' Does NOT refit — uses \code{fitted()} for in-sample and
#' \code{predict()} for new data.
#'
#' @param fits A \code{spatial_fit} object, or a named list of them
#'   (e.g. \code{list(GWR = gwr_obj, Bayesian = bayes_obj)}).  The names are
#'   used as the model labels and every element must have one; an unnamed
#'   list is an error, and so are duplicated names --- \code{model} is the key
#'   the comparison table is assembled on, so two fits sharing a name cannot be
#'   told apart in the output.
#' @param newdata Optional sf object for out-of-sample evaluation.
#'   Must contain the response variable and all predictors.
#'   If NULL, in-sample metrics are computed.
#' @param ... Extra arguments passed to predict().
#' @return A data.frame with one row per model and columns for
#'   model name and all regression metrics.
#' @family model evaluation
#' @export
evaluate_insample <- function(fits, newdata = NULL, ...) {
  # Accept a single fit

  if (inherits(fits, "spatial_fit")) {
    fits <- stats::setNames(list(fits), class(fits)[1L])
  }
  if (!is.list(fits) || length(fits) == 0L)
    stop("evaluate_insample(): `fits` must be a spatial_fit or a named list of them.")

  # The loop below is over names(fits), which is NULL for an unnamed list --
  # so lapply() would iterate zero times and this function would return NULL,
  # silently, having passed the check above.  compare_models() then turns that
  # NULL into a plain list and dies in seq_len(nrow(met_df)) with "argument
  # must be coercible to non-negative integer", nowhere near the real cause.
  nms <- names(fits)
  if (is.null(nms) || anyNA(nms) || !all(nzchar(nms)))
    stop("evaluate_insample(): `fits` must be a NAMED list of spatial_fit ",
         "objects -- the names label the models in the output. Supply them, ",
         "e.g. list(GWR = gwr_fit, Bayesian = bayes_fit).", call. = FALSE)

  # Duplicated names are silently WRONG rather than merely ambiguous.  `model`
  # is the join key compare_models() merges the metric and Moran's I tables on,
  # so two fits called "GWR" produced a 2 x 2 cross-join: four rows, every one
  # of them carrying the first fit's numbers (fits[[nm]] returns the first
  # match, so the second fit was never scored at all).
  dup <- unique(nms[duplicated(nms)])
  if (length(dup) > 0L)
    stop(sprintf(
      paste0("evaluate_insample(): `fits` has duplicated name(s) %s. The names ",
             "label the models and are the key the comparison table is built ",
             "on, so they must be unique -- give the fits distinct names, e.g. ",
             "list(GWR_bw50 = ..., GWR_bw80 = ...)."),
      paste(sQuote(dup), collapse = ", ")), call. = FALSE)

  # seq_along() rather than the names themselves: fits[[nm]] returns the first
  # element of that name, which is the other half of the duplicate-name bug.
  rows <- lapply(seq_along(fits), function(i) {
    obj <- fits[[i]]
    nm  <- nms[i]
    if (!inherits(obj, "spatial_fit")) {
      .log_warn("evaluate_insample(): '%s' is not a spatial_fit; skipping.", nm)
      return(NULL)
    }
    met <- model_metrics(obj, newdata = newdata, ...)
    cbind(data.frame(model = nm, stringsAsFactors = FALSE), met)
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}


# ---------------------------------------------------------------------------
# compare_models: pure comparison from a list of fit objects
# ---------------------------------------------------------------------------

#' Side-by-side comparison of fitted spatial models
#'
#' Takes a named list of already-fit \code{spatial_fit} objects and produces
#' a tidy comparison table including in-sample metrics and model-specific
#' information criteria (AICc, LOOIC).
#'
#' @param fits A named list of \code{spatial_fit} objects.  Names must be
#'   unique; see \code{\link{evaluate_insample}}.
#' @param newdata Optional sf for out-of-sample evaluation.
#' @param ... Extra arguments passed to predict().
#' @return A data.frame comparing all models.  Alongside the metrics it carries
#'   \code{resid_morans_I}, \code{resid_morans_z}, \code{resid_morans_p} and
#'   \code{resid_morans_null} --- the last naming which null
#'   \code{\link{residual_morans_i}} scored each model against, since that
#'   choice is per-fit and governs how much the p-value is worth.  The
#'   significant-autocorrelation warning below is driven by that p-value, so
#'   read its caveats in \code{?residual_morans_i} before treating silence as
#'   evidence of no residual structure.
#' @family model evaluation
#' @export
compare_models <- function(fits, newdata = NULL, ...) {
  # Wrap a bare fit exactly as evaluate_insample() does.  Without this, a
  # single spatial_fit passes the is.list() check below (a fit *is* a list),
  # and every downstream loop then iterates the fit's own components instead
  # of a set of models -- yielding all-NA Moran's I columns and one
  # "'fit' is not a spatial_fit object" log line per component.
  if (inherits(fits, "spatial_fit"))
    fits <- stats::setNames(list(fits), class(fits)[1L])
  if (!is.list(fits) || length(fits) == 0L)
    stop("compare_models(): `fits` must be a spatial_fit or a non-empty named list of them.")

  met_df <- evaluate_insample(fits, newdata = newdata, ...)

  # Append model-specific information criteria
  met_df$AICc  <- NA_real_
  met_df$LOOIC <- NA_real_
  met_df$bandwidth_is_fallback <- NA
  for (i in seq_len(nrow(met_df))) {
    nm <- met_df$model[i]
    # By index, not fits[[nm]]: name lookup returns the FIRST match, so with two
    # fits of one name every row read the same object.  evaluate_insample() now
    # rejects duplicate names outright, and match() keeps this loop correct
    # rather than merely lucky.  met_df can be shorter than `fits` when a
    # non-spatial_fit element was skipped, so the mapping is by name, not order.
    obj <- fits[[match(nm, names(fits))]]
    if (inherits(obj, "gwr_fit")) {
      met_df$AICc[i] <- obj$info$AICc %||% NA_real_
      is_fb <- isTRUE(obj$info$bandwidth_is_fallback)
      met_df$bandwidth_is_fallback[i] <- is_fb
      if (is_fb) {
        warning(
          sprintf(
            "compare_models(): GWR model '%s' used a fallback bandwidth (%.4f). Metrics may be unreliable; consider re-fitting with an explicit bandwidth.",
            nm, obj$info$bandwidth
          ),
          call. = FALSE
        )
      }
    }
    if (inherits(obj, "bayesian_fit"))
      met_df$LOOIC[i] <- obj$info$looic %||% NA_real_
  }

  # --- Post-fit residual spatial autocorrelation check ---
  moran_df <- .residual_morans_table(fits)
  met_df   <- merge(met_df, moran_df, by = "model", all.x = TRUE, sort = FALSE)

  # Emit warnings for models whose residuals still show significant

  # spatial autocorrelation (alpha = 0.05)
  for (i in seq_len(nrow(met_df))) {
    p_val <- met_df$resid_morans_p[i]
    I_val <- met_df$resid_morans_I[i]
    if (is.finite(p_val) && p_val < 0.05) {
      .log_warn(
        paste0("compare_models(): residuals of '%s' show significant ",
               "spatial autocorrelation (Moran's I = %.4f, p = %.4g). ",
               "The model may not fully capture the spatial structure."),
        met_df$model[i], I_val, p_val
      )
    }
  }

  met_df
}


# ---------------------------------------------------------------------------
# compare_models_cv: cross-validated comparison
# ---------------------------------------------------------------------------

#' Cross-validated comparison of spatial models
#'
#' Fits and cross-validates one or more model types, returning a unified
#' comparison table.  Unlike \code{compare_models()}, this function does
#' perform fitting (inside CV folds), because CV inherently requires
#' repeated fitting.
#'
#' @param data_sf An sf object.
#' @param response_var Response column name.
#' @param predictor_vars Predictor column names.
#' @param models Character vector: any subset of \code{c("GWR", "Bayesian",
#'   "RF")}, in any order.  Each is cross-validated on the same \code{folds}.
#'   Names outside that set raise a warning and are dropped; if nothing
#'   recognised remains, this is an error rather than a silent fallback.
#'   A recognised model whose backend package is not installed is dropped with
#'   a message so the call still returns the models that could run --- but if
#'   \emph{none} of the requested backends is installed, nothing is left to
#'   compare and the call errors with \code{"no viable models."}.  Guard with
#'   \code{requireNamespace()} when the model set is not known in advance.
#' @param k Number of folds. Default 5.
#' @param seed RNG seed. Default 123.
#' @param folds Optional precomputed fold splits.
#' @param boundary Optional polygon sf/sfc.
#' @param pointize Geometry coercion strategy.
#' @param gwr_args Extra arguments for \code{\link{cv_gwr}}.  Only names that
#'   are formal arguments of \code{cv_gwr()} are forwarded --- it has no
#'   \code{...} --- so entries meant for \code{fit_gwr_model()} alone (e.g.
#'   \code{longlat}) cannot be passed this way.  Anything dropped is named in a
#'   warning; call \code{cv_gwr()} directly if you need it.
#' @param bayes_args Extra arguments for \code{fit_bayesian_spatial_model()}.
#'   Forwarded whole as \code{cv_bayes(fit_args = )}, so an unrecognised name
#'   is an error from \code{fit_bayesian_spatial_model()} rather than a silent
#'   drop.  \code{compute_loo}, \code{boundary} and \code{pointize} are
#'   overridden by the CV internals.
#' @param rf_args Extra arguments for \code{\link{cv_rf}}, which passes
#'   anything it does not recognise on to \code{\link{fit_rf_model}} and thence
#'   to \code{ranger::ranger()}.
#' @param summary "mean" or "median" for Bayesian predictions.
#' @param quiet Logical; suppress messages.
#' @return A list with overall, by_fold, and per-model cv_results
#'   (\code{gwr_cv}, \code{bayes_cv}, \code{rf_cv} for the models that ran).
#'   Only the models that actually ran appear, so check which names are present
#'   rather than assuming one entry per requested model: a backend whose package
#'   is missing is dropped with a message.  When \strong{no} requested backend
#'   is available there is nothing to return and the function errors with
#'   \code{"no viable models."} instead of returning an empty comparison.
#' @family model evaluation
#' @examples
#' \donttest{
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 120
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000), elev = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$price <- 10 + 0.01 * st_coordinates(dat)[, 1] + 2 * dat$elev + rnorm(n)
#'   cmp <- compare_models_cv(dat, "price", "elev", models = "RF", k = 3,
#'                            rf_args = list(num_trees = 100))
#'   cmp$overall
#' }
#' }
#' @export
compare_models_cv <- function(
    data_sf, response_var, predictor_vars,
    models = c("GWR", "Bayesian"),
    k = 5, seed = 123, folds = NULL, boundary = NULL, pointize = "auto",
    gwr_args = list(), bayes_args = list(), rf_args = list(),
    summary = c("mean", "median"),
    quiet = FALSE
) {
  summary <- match.arg(summary)
  .msg <- function(...) if (!quiet) message(...)

  if (!inherits(data_sf, "sf"))
    stop("compare_models_cv(): `data_sf` must be an sf object.")

  # Unrecognised model names used to be dropped by a bare intersect(), so
  # models = "RF" silently ran GWR instead of what was asked for -- a wrong
  # answer presented as the right one.
  known    <- c("GWR", "Bayesian", "RF")
  models   <- unique(as.character(models))
  unknown  <- setdiff(models, known)
  models   <- intersect(models, known)
  if (length(unknown) > 0L)
    warning(sprintf(
      "compare_models_cv(): ignoring unknown model(s) %s. Supported: %s.",
      paste(sQuote(unknown), collapse = ", "),
      paste(sQuote(known), collapse = ", ")), call. = FALSE)
  if (length(models) == 0L)
    stop(sprintf("compare_models_cv(): no recognised model requested. Supported: %s.",
                 paste(sQuote(known), collapse = ", ")), call. = FALSE)

  for (m in models) {
    if (!.model_available(m)) {
      .msg(sprintf("compare_models_cv(): dropping %s (package/function unavailable).", m))
      models <- setdiff(models, m)
    }
  }
  if (length(models) == 0L) stop("compare_models_cv(): no viable models.")

  if (!("..row_id" %in% names(data_sf)))
    data_sf$`..row_id` <- seq_len(nrow(data_sf))

  cleanup <- .with_seed(seed)
  on.exit(cleanup(), add = TRUE)

  comparison_rows <- list(); by_fold_rows <- list(); cv_results <- list()

  if ("GWR" %in% models) {
    .msg("compare_models_cv(): running CV for GWR ...")
    .warn_dropped_args(cv_gwr, gwr_args, "gwr_args", "cv_gwr")
    base <- .merge_args(
      list(data_sf = data_sf, response_var = response_var,
           predictor_vars = predictor_vars, folds = folds,
           k = k, seed = seed, boundary = boundary, pointize = pointize),
      gwr_args, "gwr_args"
    )
    gwr_cv <- do.call(cv_gwr, .filter_args(cv_gwr, base))
    cv_results$gwr_cv <- gwr_cv
    ov <- try(as.data.frame(gwr_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, SMAPE = NA_real_,
                       R2 = NA_real_, Adj_R2 = NA_real_)
    ov$model <- "GWR"
    comparison_rows[["GWR"]] <- ov
    bf <- try(as.data.frame(gwr_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      bf$model <- "GWR"; by_fold_rows[["GWR"]] <- bf
    }
  }

  if ("Bayesian" %in% models) {
    .msg("compare_models_cv(): running CV for Bayesian ...")
    base <- list(data_sf = data_sf, response_var = response_var,
                 predictor_vars = predictor_vars, folds = folds,
                 k = k, seed = seed, boundary = boundary,
                 pointize = pointize, summary = summary,
                 fit_args = bayes_args)
    bayes_cv <- do.call(cv_bayes, .filter_args(cv_bayes, base))
    cv_results$bayes_cv <- bayes_cv
    ov <- try(as.data.frame(bayes_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, SMAPE = NA_real_,
                       R2 = NA_real_, Adj_R2 = NA_real_)
    ov$model <- "Bayesian"
    comparison_rows[["Bayesian"]] <- ov
    bf <- try(as.data.frame(bayes_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      bf$model <- "Bayesian"; by_fold_rows[["Bayesian"]] <- bf
    }
  }

  if ("RF" %in% models) {
    .msg("compare_models_cv(): running CV for RF ...")
    # cv_rf() has `...` (forwarded to fit_rf_model() and on to ranger), so
    # nothing in rf_args is filtered out; `pointize` is deliberately not
    # passed, because cv_rf() has no such formal and it would land in ranger.
    base <- .merge_args(
      list(data_sf = data_sf, response_var = response_var,
           predictor_vars = predictor_vars, folds = folds,
           k = k, seed = seed, boundary = boundary),
      rf_args, "rf_args"
    )
    rf_cv <- do.call(cv_rf, base)
    cv_results$rf_cv <- rf_cv
    ov <- try(as.data.frame(rf_cv$overall), silent = TRUE)
    if (inherits(ov, "try-error") || nrow(ov) == 0L)
      ov <- data.frame(RMSE = NA_real_, MAE = NA_real_, MAPE = NA_real_, SMAPE = NA_real_,
                       R2 = NA_real_, Adj_R2 = NA_real_)
    ov$model <- "RF"
    comparison_rows[["RF"]] <- ov
    bf <- try(as.data.frame(rf_cv$fold_metrics), silent = TRUE)
    if (!inherits(bf, "try-error") && nrow(bf)) {
      bf$model <- "RF"; by_fold_rows[["RF"]] <- bf
    }
  }

  c(list(overall = dplyr::bind_rows(comparison_rows),
         by_fold = dplyr::bind_rows(by_fold_rows)),
    cv_results)
}
