# =============================================================================
# GWR-specific model selection (wrapper for GWmodel::gwr.model.selection)
# =============================================================================
#
# Design note
# -----------
# Every call into GWmodel is confined to .gwr_ms_engine().  Everything after
# it -- criterion extraction, model-list normalisation, ranking -- is pure and
# is tested directly, and the engine itself is injectable so the assembly path
# can be exercised without GWmodel installed.  This matters because GWmodel is
# in Suggests: without the split, the entire function would be untested on any
# machine that does not have it, which is every CRAN check machine.
# =============================================================================


#' Extract the ranking criterion from GWmodel's model-selection diagnostics
#'
#' \code{GWmodel::gwr.model.selection()} documents its diagnostic table
#' (\code{GWR.df}) as "a data frame consited of four columns: bandwidth, AIC,
#' AICc, RSS" -- so \strong{AICc is column 3}, and column 2 is the uncorrected
#' AIC.  GWmodel builds the table with \code{rbind()} over unnamed vectors
#' (\code{c(bw, aic.rss[2], aic.rss[3], aic.rss[1])} in \code{Model.selection.r}),
#' so it never carries column names: the by-name branch below cannot fire on a
#' real GWmodel return, and the positional fallback is the path every real call
#' takes.  Getting it wrong is not a fallback-only risk, it is the normal case
#' -- reading column 2 ranks on AIC while labelling the result AICc, which
#' selects larger models than AICc would and can carry a pure-noise predictor.
#'
#' The by-name branch is kept for a future GWmodel that labels the table, and
#' for injected test engines.
#'
#' @param gwr_df The second element of the \code{gwr.model.selection()} return.
#' @param criterion Column name to look for, matched case-insensitively.
#' @return A list with \code{values} (numeric, non-finite coerced to
#'   \code{NA}), \code{column}, \code{column_name}, \code{by_name},
#'   \code{shape_ok} (\code{TRUE} when the table has the documented four
#'   columns, or when the column was located by name) and \code{label}.
#' @keywords internal
#' @noRd
.gwr_ms_criterion <- function(gwr_df, criterion = "AICc") {
  if (is.null(gwr_df))
    stop(".gwr_ms_criterion(): GWmodel returned no diagnostic table.",
         call. = FALSE)

  # A bare numeric vector (no dim) is unambiguous: it IS the criterion.
  if (is.atomic(gwr_df) && is.null(dim(gwr_df))) {
    vals <- suppressWarnings(as.numeric(gwr_df))
    vals[!is.finite(vals)] <- NA_real_
    return(list(values = vals, column = 1L, column_name = NA_character_,
                by_name = FALSE, shape_ok = FALSE,
                label = sprintf("%s (assumed: unlabelled vector)", criterion)))
  }

  is_df <- is.data.frame(gwr_df)
  ncols <- if (is_df) length(gwr_df) else ncol(gwr_df)
  if (is.null(ncols) || is.na(ncols) || ncols < 1L)
    stop(".gwr_ms_criterion(): diagnostic table has no columns.", call. = FALSE)

  cn      <- if (is_df) names(gwr_df) else colnames(gwr_df)
  col     <- NA_integer_
  by_name <- FALSE
  if (!is.null(cn)) {
    hit <- which(tolower(trimws(cn)) == tolower(criterion))
    if (length(hit) >= 1L) {
      col     <- as.integer(hit[[1L]])
      by_name <- TRUE
    }
  }
  # Fallback.  GWmodel's GWR.df is documented AND implemented as
  #   c(bandwidth, AIC, AICc, RSS)
  # so AICc lives in column 3.  Column 2 is the UNCORRECTED AIC, which
  # penalises the effective parameter count more weakly and therefore selects
  # larger models -- on a forward sweep it will happily keep a pure-noise
  # predictor that AICc drops, which is the failure AICc exists to prevent.
  # Anything other than the documented four columns is a shape this code has
  # never seen; read the same position but tell the caller to say so.
  shape_ok <- TRUE
  if (is.na(col)) {
    col      <- if (ncols >= 3L) 3L else ncols   # ncols is 1 or 2 here
    shape_ok <- identical(as.integer(ncols), 4L)
  }

  raw  <- if (is_df) gwr_df[[col]] else gwr_df[, col]
  vals <- suppressWarnings(as.numeric(raw))
  vals[!is.finite(vals)] <- NA_real_

  col_name <- if (!is.null(cn) && nzchar(cn[[col]])) cn[[col]] else NA_character_
  # Say which column was actually read.  Claiming "unlabelled" when the table
  # IS labelled -- just not with the name we wanted -- hides the fact that the
  # ranking may be built on a different criterion entirely.
  label <- if (by_name) col_name
    else if (!is.na(col_name))
      sprintf("%s (assumed: column %d, labelled '%s')", criterion, col, col_name)
    else sprintf("%s (assumed: column %d, unlabelled)", criterion, col)

  list(values = vals, column = col, column_name = col_name,
       by_name = by_name, shape_ok = shape_ok, label = label)
}


#' Normalise GWmodel's model list into character vectors of predictors
#'
#' Handles the three shapes the element can take: a two-element list of the
#' model's formula \emph{string} and its predictor vector -- which is what
#' GWmodel returns, \code{list(Generate.formula(DeVar, vars), vars)}, i.e.
#' \code{list("z ~ a + b", c("a", "b"))} -- a formula object, or a bare
#' character vector.
#'
#' @param model_list The first element of the \code{gwr.model.selection()}
#'   return.
#' @param response_var Response name, stripped when it appears in the element.
#' @return A list of character vectors.
#' @keywords internal
#' @noRd
.gwr_ms_varsets <- function(model_list, response_var) {
  if (!is.list(model_list) || length(model_list) == 0L)
    stop(".gwr_ms_varsets(): GWmodel returned an empty model list.",
         call. = FALSE)

  lapply(model_list, function(m) {
    if (inherits(m, "formula")) {
      return(setdiff(all.vars(m), response_var))
    }
    if (is.list(m)) {
      # GWmodel stores list("<response> ~ <v1>+<v2>", c(v1, v2)) -- the model's
      # formula as a STRING in the first slot, the predictor names in the
      # second.  Take the second slot; the setdiff() below is what keeps a
      # first-slot fallback from returning the formula text as a variable name.
      part <- if (length(m) >= 2L) m[[2L]] else m[[1L]]
      return(setdiff(as.character(unlist(part, use.names = FALSE)),
                     response_var))
    }
    setdiff(as.character(m), response_var)
  })
}


#' Assemble and rank the model-selection table
#'
#' Ties on the criterion are broken toward the smaller model, so a variable
#' that buys nothing is not carried along by an arbitrary ordering.  Models
#' whose criterion could not be evaluated sort last and never win.
#'
#' @param varsets List of character vectors, one per model.
#' @param values Numeric criterion, one per model (\code{NA} allowed).
#' @param minimise Logical; \code{TRUE} for AICc/CV, \code{FALSE} for a
#'   criterion where larger is better.
#' @return A data.frame with \code{rank}, \code{n_vars}, \code{variables},
#'   \code{criterion}.
#' @keywords internal
#' @noRd
.gwr_ms_table <- function(varsets, values, minimise = TRUE) {
  n <- length(varsets)
  if (length(values) != n)
    stop(sprintf(paste0(".gwr_ms_table(): GWmodel returned %d model(s) but ",
                        "%d criterion value(s); the two cannot be aligned."),
                 n, length(values)), call. = FALSE)

  crit <- as.numeric(values)
  crit[!is.finite(crit)] <- NA_real_

  df <- data.frame(
    rank      = NA_integer_,
    n_vars    = vapply(varsets, length, integer(1)),
    variables = vapply(varsets, function(v)
      if (length(v)) paste(v, collapse = " + ") else "<none>", character(1)),
    criterion = crit,
    stringsAsFactors = FALSE
  )

  ord <- order(if (minimise) df$criterion else -df$criterion,
               df$n_vars, na.last = TRUE)
  df  <- df[ord, , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  rownames(df) <- NULL
  attr(df, "varsets") <- varsets[ord]
  df
}


#' Run an expression, optionally swallowing GWmodel's progress output
#'
#' GWmodel reports progress with bare \code{cat()}, so neither
#' \code{suppressMessages()} nor \code{suppressWarnings()} touches it. A sweep
#' prints one "Now calibrating the model" block per candidate -- 190 of them at
#' 19 candidates -- plus the full golden-section bandwidth trace.
#'
#' @param expr Expression to evaluate (forced exactly once).
#' @param quiet Logical; discard output written to stdout.
#' @return The value of \code{expr}.
#' @keywords internal
#' @noRd
.gwr_quietly <- function(expr, quiet = TRUE) {
  if (!isTRUE(quiet)) return(expr)
  res <- NULL
  # capture.output() restores the sink on error via its own on.exit(), which a
  # hand-rolled sink() would not.
  invisible(utils::capture.output(res <- expr))
  res
}


#' Run GWmodel's model-selection sweep
#'
#' All GWmodel interaction lives here so the rest of the function is testable
#' without it.
#'
#' @return A list with \code{model_list}, \code{gwr_df}, \code{bandwidth},
#'   \code{bandwidth_source}, \code{used_dmat} and \code{raw} (GWmodel's
#'   unmodified return, which the caller passes straight through to the public
#'   result).  \code{.engine} is a documented injection point, so this list is
#'   the contract an injected engine has to satisfy -- \code{raw} included.
#' @keywords internal
#' @noRd
.gwr_ms_engine <- function(dat, response_var, candidate_vars, bandwidth,
                           adaptive, kernel, bw_approach, dmat_max_n,
                           quiet = TRUE) {
  if (!requireNamespace("GWmodel", quietly = TRUE))
    stop("gwr_model_selection(): package 'GWmodel' is required. ",
         "Install with install.packages('GWmodel').", call. = FALSE)
  if (!requireNamespace("sp", quietly = TRUE))
    stop("gwr_model_selection(): package 'sp' is required (GWmodel interop).",
         call. = FALSE)

  n_obs   <- nrow(dat)
  n_cand  <- length(candidate_vars)
  sp_dat  <- .to_sp(dat, unique(c(response_var, candidate_vars)))

  # One distance matrix reused across every candidate fit and the bandwidth
  # search.  GWmodel recomputes distances per call otherwise, and the sweep
  # makes n_cand*(n_cand+1)/2 + bandwidth-search calls.  n^2 doubles: 2000
  # points is 32 MB, which is why this is capped rather than unconditional.
  dMat <- NULL
  # `Inf` means no cap -- always precompute.  It was gated on is.finite(), so
  # Inf meant NEVER precompute: with adaptive = FALSE the sweep then died on
  # n = 100 (a fixed bandwidth needs the matrix) and blamed the caller.
  if (!is.na(dmat_max_n) && n_obs <= dmat_max_n) {
    dMat <- tryCatch(
      GWmodel::gw.dist(dp.locat = sp::coordinates(sp_dat)),
      error = function(e) {
        .log_warn("gwr_model_selection(): gw.dist() failed (%s); GWmodel will recompute distances per model.",
                  conditionMessage(e))
        NULL
      }
    )
  }

  # --- Bandwidth ------------------------------------------------------------
  # One bandwidth is used for every candidate model: that is what makes the
  # criterion comparable across them, and it is also the method's main
  # limitation (see the function's documentation).
  bandwidth_source <- "supplied"
  if (is.null(bandwidth)) {
    full_fml <- stats::reformulate(termlabels = candidate_vars,
                                   response = response_var)
    # bw.gwr() distinguishes a supplied distance matrix with missing(dMat),
    # NOT is.null(dMat).  Passing dMat = NULL explicitly therefore takes the
    # *supplied* branch, where dim(NULL)[1] != n evaluates to logical(0) and
    # the if() dies with "missing value where TRUE/FALSE needed".  dMat is NULL
    # on three ordinary paths -- dmat_max_n = 0, n_obs above dmat_max_n
    # (default 2000), and a failed gw.dist() -- so every dataset over 2000
    # points used to land in the fallback-bandwidth branch below with the
    # message blaming GWmodel.  Build the call instead and simply omit the
    # argument, which is what makes missing(dMat) true.
    bw_args <- list(formula = full_fml, data = sp_dat, approach = bw_approach,
                    kernel = kernel, adaptive = adaptive)
    if (!is.null(dMat)) bw_args$dMat <- dMat
    bw <- tryCatch(
      .gwr_quietly(suppressWarnings(do.call(GWmodel::bw.gwr, bw_args)), quiet),
      error = function(e) {
        .log_warn("gwr_model_selection(): bw.gwr() failed: %s",
                  conditionMessage(e))
        NA_real_
      }
    )
    if (!isTRUE(is.finite(bw)) || bw <= 0) {
      bw <- .fallback_bandwidth(sp_dat, adaptive)
      bandwidth_source <- "fallback"
      warning(sprintf(paste0("gwr_model_selection(): automatic bandwidth ",
                             "selection failed; using arbitrary fallback ",
                             "bandwidth = %.4f. Supply `bandwidth` explicitly ",
                             "for a defensible comparison."), bw),
              call. = FALSE)
    } else {
      bandwidth_source <- sprintf("bw.gwr(%s) on the full model", bw_approach)
    }
  } else {
    bw <- as.numeric(bandwidth)
  }

  # Clamp against the LARGEST model in the sweep (all candidates + intercept),
  # not the smallest.  A bandwidth that supports a one-variable fit can leave
  # the full model underdetermined, and GWmodel would fail partway through.
  if (isTRUE(adaptive)) {
    bw      <- as.integer(round(bw))
    min_bw  <- n_cand + 2L
    if (min_bw > n_obs)
      stop(sprintf(paste0("gwr_model_selection(): %d observations cannot ",
                          "support a local fit of %d candidate predictors ",
                          "plus an intercept."), n_obs, n_cand), call. = FALSE)
    if (bw < min_bw) {
      .log_warn("gwr_model_selection(): adaptive bandwidth %d is too small for the full %d-predictor model; clamping to %d.",
                bw, n_cand, min_bw)
      bw <- min_bw
    }
    if (bw > n_obs) bw <- as.integer(n_obs)
  }

  # A FIXED bandwidth needs the distance matrix.  Without one,
  # gwr.model.selection() sets dMat <- matrix(0, 0, 0) and then asserts
  # stopifnot(bw > min(dMat)); min() of an empty matrix is Inf, so the
  # assertion can never hold and the sweep dies with the bare
  # "(bw > min(dMat)) is not TRUE".  Say so here, where the remedies are
  # nameable, rather than letting a fixed-bandwidth sweep on n > dmat_max_n
  # stop unexplained.
  if (!isTRUE(adaptive) && is.null(dMat))
    stop(sprintf(paste0("gwr_model_selection(): a fixed bandwidth ",
                        "(adaptive = FALSE) requires the precomputed n x n ",
                        "distance matrix, and none was built for these %d ",
                        "observations%s. GWmodel asserts `bw > min(dMat)`, and ",
                        "min() of the empty matrix it substitutes is Inf, so ",
                        "the sweep cannot start. Either allow the matrix -- ",
                        "`dmat_max_n = %d` or more, costing about %.3g MB -- ",
                        "or use `adaptive = TRUE`, which needs no distance ",
                        "matrix."),
                 n_obs,
                 if (is.finite(dmat_max_n) && n_obs > dmat_max_n)
                   sprintf(" (dmat_max_n = %s)", format(dmat_max_n)) else "",
                 n_obs, (as.numeric(n_obs)^2 * 8) / 1024^2),
         call. = FALSE)

  # `approach` is deliberately NOT forwarded.  It only steers gwr.model.selection's
  # own bandwidth search, and bw is always explicit here; GWmodel's documented
  # example leaves it at its "CV" default and still reports AICc in the
  # diagnostic table, so passing it could only introduce a version-dependent
  # difference for no gain.
  #
  # `dMat` is omitted rather than passed as NULL for the same reason as in
  # bw.gwr() above: gwr.model.selection() branches on missing(dMat).
  ms_args <- list(DeVar = response_var, InDeVars = candidate_vars,
                  data = sp_dat, bw = bw, adaptive = adaptive, kernel = kernel)
  if (!is.null(dMat)) ms_args$dMat <- dMat
  res <- tryCatch(
    .gwr_quietly(do.call(GWmodel::gwr.model.selection, ms_args), quiet),
    error = function(e)
      stop(sprintf("gwr_model_selection(): gwr.model.selection() failed: %s",
                   conditionMessage(e)), call. = FALSE)
  )

  if (!is.list(res) || length(res) < 2L)
    stop("gwr_model_selection(): gwr.model.selection() returned an unexpected ",
         "shape (expected a list of a model list and a diagnostic table).",
         call. = FALSE)

  model_list <- if (!is.null(res$model.list)) res$model.list else res[[1L]]
  gwr_df     <- if (!is.null(res$GWR.df))     res$GWR.df     else res[[2L]]

  list(model_list = model_list, gwr_df = gwr_df,
       bandwidth = bw, bandwidth_source = bandwidth_source,
       used_dmat = !is.null(dMat), raw = res)
}


# -----------------------------------------------------------------------------
# Exported
# -----------------------------------------------------------------------------

#' Forward model selection for geographically weighted regression
#'
#' Wraps \code{GWmodel::gwr.model.selection()}, which grows a GWR model one
#' predictor at a time and scores every intermediate model with a corrected
#' Akaike information criterion, and returns the results as a ranked table
#' rather than the two loosely-coupled lists GWmodel produces.
#'
#' @section What this optimises, and what it does not:
#' The criterion is \strong{in-sample}. AICc penalises the effective number of
#' parameters, so it is not the same thing as maximising fit, but it is still
#' computed on the data the model was fitted to, and under spatial
#' autocorrelation an in-sample criterion is optimistic in a way that a
#' spatially blocked estimate is not. Treat this as fast screening.
#' \code{\link{select_features_forward}} performs the same forward search
#' against a spatially blocked cross-validated score; it costs far more and is
#' the one to trust when the answer matters. When the two disagree, the
#' disagreement is itself informative -- it usually means a candidate is
#' predictive only locally.
#'
#' Two further limitations are structural rather than incidental:
#'
#' \itemize{
#'   \item \strong{One bandwidth for every model.} Comparing criteria across
#'     models requires holding the smoothing fixed, but the bandwidth is itself
#'     a fitted quantity, and the value chosen for the full model is not
#'     optimal for a one-predictor model. This is how the method is defined
#'     (Lu et al. 2014); it is not an implementation shortcut. Refit the
#'     selected model with \code{bandwidth = NULL} to re-optimise once the
#'     variable set is settled.
#'   \item \strong{The null model is never evaluated.} The sweep starts from
#'     one predictor, so the result always names at least one. It cannot tell
#'     you that none of the candidates help.
#' }
#'
#' @section Cost:
#' The sweep fits \code{p * (p + 1) / 2} GWR models for \code{p} candidates --
#' 55 at p = 10, 210 at p = 20 -- each over all \code{n} locations.
#' \code{max_models} stops the call rather than letting it run for hours.
#'
#' @param data_sf An \code{sf} object with response, predictors and geometry.
#' @param response_var Response column name.
#' @param candidate_vars Character vector of at least two predictors to choose
#'   among.
#' @param bandwidth Bandwidth held fixed across all candidate models. If
#'   \code{NULL} (default) it is selected with \code{GWmodel::bw.gwr()} on the
#'   model containing every candidate. Integer neighbour count when
#'   \code{adaptive = TRUE}, distance in CRS units otherwise.
#' @param adaptive Logical; adaptive (nearest-neighbour) bandwidth. Default
#'   \code{TRUE}.
#' @param kernel Weighting kernel. One of \code{"bisquare"} (default),
#'   \code{"gaussian"}, \code{"tricube"}, \code{"boxcar"},
#'   \code{"exponential"}.
#' @param bw_approach Criterion for the bandwidth search: \code{"AICc"}
#'   (default) or \code{"CV"}. Matches \code{\link{fit_gwr_model}}'s default.
#' @param max_models Refuse the call if the sweep would exceed this many model
#'   fits. Default 200, which admits up to 19 candidates.
#' @param dmat_max_n Precompute and reuse an \code{n} x \code{n} distance
#'   matrix when \code{n} is at most this. Default 2000 (about 32 MB). Set to 0
#'   to disable.
#' @param quiet Discard GWmodel's progress output. Default \code{TRUE}, because
#'   GWmodel writes it with bare \code{cat()} -- which no
#'   \code{suppressMessages()} can silence -- and emits one block per candidate
#'   model, so it scales with the square of the candidate count. Set
#'   \code{FALSE} to watch a long sweep progress.
#' @param .engine Internal; injectable backend used for testing.
#'
#' @return An object of class \code{gwr_model_selection}, a list with:
#'   \code{best} (character vector of the selected predictors);
#'   \code{table} (ranked data.frame of every model evaluated, with columns
#'   \code{rank}, \code{n_vars}, \code{variables} and \code{criterion});
#'   \code{criterion} (label for the criterion actually read, noting when it
#'   had to be located positionally);
#'   \code{response_var} and \code{candidate_vars} (the response and the full
#'   candidate set the sweep ran over, both echoed by \code{print()});
#'   \code{bandwidth}, \code{bandwidth_source}, \code{adaptive} and
#'   \code{kernel} (the smoothing held fixed across the sweep, and where it
#'   came from);
#'   \code{n_obs}, \code{n_models}, \code{used_dmat}; and \code{raw}
#'   (GWmodel's unmodified return: the two-element list of its model list and
#'   its diagnostic table).
#'
#' @section Using $raw with GWmodel directly:
#' \code{raw} is GWmodel's own \code{list(model.list, GWR.df)}, so its two
#' elements have to be unpacked before GWmodel's own helpers will take them:
#' \code{GWmodel::gwr.model.view()} takes \code{(DeVar, InDeVars, model.list)},
#' so the call is
#' \preformatted{
#'   GWmodel::gwr.model.view(sel$response_var, sel$candidate_vars, sel$raw[[1]])
#' }
#' -- \code{sel$raw[[1]]}, not \code{sel$raw}.  The diagnostic table is
#' \code{sel$raw[[2]]}, an unlabelled numeric matrix whose columns are
#' \code{bandwidth}, \code{AIC}, \code{AICc}, \code{RSS} in that order; the
#' \code{criterion} column of \code{$table} is its third column.
#'
#' @references
#' Lu, B., Harris, P., Charlton, M. and Brunsdon, C. (2014). The GWmodel R
#' package: further topics for exploring spatial heterogeneity using
#' geographically weighted models. \emph{Geo-spatial Information Science}
#' 17(2), 85--101. \doi{10.1080/10095020.2014.917453}
#'
#' @seealso \code{\link{select_features_forward}} for the blocked
#'   cross-validated counterpart, \code{\link{fit_gwr_model}} to fit the
#'   selected model.
#' @family cross-validation
#' @examples
#' \donttest{
#' if (requireNamespace("GWmodel", quietly = TRUE) &&
#'     requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 80
#'   dat <- st_as_sf(
#'     data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
#'   sel <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30)
#'   sel$best
#'   fit <- fit_gwr_model(dat, "z", sel$best)
#' }
#' }
#' @export
gwr_model_selection <- function(data_sf, response_var, candidate_vars,
                                bandwidth = NULL,
                                adaptive = TRUE,
                                kernel = c("bisquare", "gaussian", "tricube",
                                           "boxcar", "exponential"),
                                bw_approach = c("AICc", "CV"),
                                max_models = 200L,
                                dmat_max_n = 2000L,
                                quiet = TRUE,
                                .engine = .gwr_ms_engine) {
  if (!inherits(data_sf, "sf"))
    stop("gwr_model_selection(): `data_sf` must be an sf object.", call. = FALSE)
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("gwr_model_selection(): `response_var` must be a single column name.",
         call. = FALSE)
  # match.arg() is the only kernel validation this package needs: an invalid
  # value never gets past it.  (.validate_kernel() in R/model-gwr.R exists and
  # is called by cv_gwr(), but only ever *after* that function's own
  # match.arg(), so it is unreachable there too -- see its @noRd block.)
  kernel      <- match.arg(kernel)
  bw_approach <- match.arg(bw_approach)

  # Validate `bandwidth` here, not inside .gwr_ms_engine(): the engine is an
  # injectable test seam and requires GWmodel, so a check placed there would
  # be skipped for injected engines and unreachable without GWmodel.
  # Unvalidated, NA reaches the engine's clamp as "missing value where
  # TRUE/FALSE needed", a length-2 vector as "the condition has length > 1",
  # and with adaptive = FALSE (no clamp at all) a zero or negative distance
  # goes straight into GWmodel.
  if (!is.null(bandwidth) &&
      (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
       !is.finite(bandwidth) || bandwidth <= 0))
    stop(sprintf(paste0("gwr_model_selection(): `bandwidth` must be a single ",
                        "positive finite number (or NULL to select one ",
                        "automatically). With adaptive = %s it is %s."),
                 if (isTRUE(adaptive)) "TRUE" else "FALSE",
                 if (isTRUE(adaptive))
                   "the number of nearest neighbours in each local window"
                 else
                   "a distance in the CRS units of `data_sf`"),
         call. = FALSE)

  candidate_vars <- unique(as.character(candidate_vars))
  missing_v <- setdiff(c(response_var, candidate_vars), names(data_sf))
  if (length(missing_v) > 0L)
    stop("gwr_model_selection(): column(s) absent from `data_sf`: ",
         paste(sQuote(missing_v), collapse = ", "), call. = FALSE)
  if (length(candidate_vars) < 2L)
    stop("gwr_model_selection(): at least two candidate variables are needed; ",
         "there is nothing to select among otherwise.", call. = FALSE)

  p        <- length(candidate_vars)
  # as.numeric() before the multiply: p * (p + 1L) is integer arithmetic, which
  # overflows to NA (with a warning) above p = 46341, and `NA > max_models` is
  # NA -- so the guard that exists to refuse an impossible sweep would itself
  # error on the very inputs it is meant to catch.
  n_models <- as.numeric(p) * (as.numeric(p) + 1) / 2
  if (n_models > max_models)
    stop(sprintf(paste0("gwr_model_selection(): %d candidates would fit about ",
                        "%d GWR models, above max_models = %d. Screen the ",
                        "candidates first or raise max_models deliberately."),
                 p, n_models, max_models), call. = FALSE)

  dat <- prep_model_data(data_sf = data_sf, response_var = response_var,
                         predictor_vars = candidate_vars, pointize = "auto")

  eng <- .engine(dat = dat, response_var = response_var,
                 candidate_vars = candidate_vars, bandwidth = bandwidth,
                 adaptive = adaptive, kernel = kernel,
                 bw_approach = bw_approach, dmat_max_n = dmat_max_n,
                 quiet = quiet)

  varsets <- .gwr_ms_varsets(eng$model_list, response_var)
  crit    <- .gwr_ms_criterion(eng$gwr_df, criterion = "AICc")
  tab     <- .gwr_ms_table(varsets, crit$values, minimise = TRUE)

  # GWmodel builds GWR.df with rbind() over unnamed vectors, so it never
  # carries column names and this branch is the normal path, not an edge case.
  # The documented column order is c(bandwidth, AIC, AICc, RSS): AICc is
  # column 3.  Reading column 2 would rank on the uncorrected AIC while
  # labelling the answer AICc.
  if (!crit$by_name)
    .log_info(paste0("gwr_model_selection(): GWmodel's diagnostic table carries ",
                     "no column names (it is built by rbind() over unnamed ",
                     "vectors), so AICc is read positionally from column %d of ",
                     "the documented c(bandwidth, AIC, AICc, RSS) layout. Check ",
                     "$raw if the ranking looks wrong."), crit$column)
  # A shape other than those four columns means the assumption above no longer
  # holds and the ranking may be built on a different criterion entirely.
  if (!isTRUE(crit$shape_ok) && !isTRUE(crit$by_name))
    .log_warn(paste0("gwr_model_selection(): GWmodel's diagnostic table does ",
                     "not have the documented four columns (bandwidth, AIC, ",
                     "AICc, RSS); AICc was read from column %d of a table with ",
                     "a shape this version does not recognise. Treat the ",
                     "ranking as unverified and check $raw."), crit$column)

  best_set <- attr(tab, "varsets")[[1L]]
  n_failed <- sum(is.na(tab$criterion))
  if (n_failed > 0L)
    .log_warn("gwr_model_selection(): %d of %d models produced a non-finite criterion and were ranked last.",
              n_failed, nrow(tab))
  if (n_failed == nrow(tab))
    stop("gwr_model_selection(): no model produced a finite criterion; ",
         "nothing can be selected.", call. = FALSE)

  structure(
    list(
      best             = best_set,
      table            = tab,
      criterion        = crit$label,
      response_var     = response_var,
      candidate_vars   = candidate_vars,
      bandwidth        = eng$bandwidth,
      bandwidth_source = eng$bandwidth_source,
      adaptive         = adaptive,
      kernel           = kernel,
      n_obs            = nrow(dat),
      n_models         = nrow(tab),
      used_dmat        = isTRUE(eng$used_dmat),
      raw              = eng$raw
    ),
    class = "gwr_model_selection"
  )
}


#' Print a GWR model selection result
#'
#' Shows the forward-selection trail: the response, the candidate predictors,
#' and the top-ranked models with their criterion values, so you can see both
#' which model won and by how much.  A shallow gap between the first few rows
#' means the ranking is not well identified and the choice of predictors should
#' not be treated as settled -- worth checking before reporting one model as
#' the selected one.
#'
#' @param x A \code{gwr_model_selection} object.
#' @param n Number of top-ranked models to show. Default 10.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.gwr_model_selection <- function(x, n = 10L, ...) {
  cat("Geographically weighted regression - forward model selection\n\n")
  cat(sprintf("  response    : %s\n", x$response_var))
  cat(sprintf("  candidates  : %d (%s)\n", length(x$candidate_vars),
              paste(x$candidate_vars, collapse = ", ")))
  cat(sprintf("  observations: %d\n", x$n_obs))
  cat(sprintf("  models      : %d\n", x$n_models))
  cat(sprintf("  bandwidth   : %s (%s; %s)\n", format(x$bandwidth),
              if (isTRUE(x$adaptive)) "adaptive" else "fixed",
              x$bandwidth_source))
  cat(sprintf("  kernel      : %s\n", x$kernel))
  cat(sprintf("  criterion   : %s\n\n", x$criterion))

  n_show <- max(1L, min(as.integer(n), nrow(x$table)))
  print(utils::head(x$table, n_show), row.names = FALSE)
  if (n_show < nrow(x$table))
    cat(sprintf("  ... %d more (see $table)\n", nrow(x$table) - n_show))

  cat(sprintf("\nSelected: %s ~ %s\n", x$response_var,
              paste(x$best, collapse = " + ")))
  cat("\nThe criterion is in-sample and every model shares one bandwidth.\n")
  cat("Confirm with select_features_forward() before relying on this.\n")
  invisible(x)
}
