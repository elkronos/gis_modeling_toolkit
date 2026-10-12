# =============================================================================
# GWR-specific model selection (the forward sweep of
# GWmodel::gwr.model.selection, scored by the package's local solve)
# =============================================================================
#
# Design note
# -----------
# Every model is scored inside .gwr_ms_engine().  Everything after it --
# criterion extraction, model-list normalisation, ranking -- is pure and is
# tested directly, and the engine itself is injectable so the assembly path
# can be exercised on its own.  The engine returns its sweep in the layout
# GWmodel's function returns, which is what the helpers below read.
# =============================================================================


#' Extract the ranking criterion from a model-selection diagnostic table
#'
#' The table has the layout of \code{GWmodel::gwr.model.selection()}'s
#' \code{GWR.df}, which GWmodel documents as "a data frame consited of four
#' columns: bandwidth, AIC, AICc, RSS".  So \strong{AICc is column 3}, and
#' column 2 is the uncorrected AIC.  The package's own engine names the four
#' columns and the criterion is found by name.  GWmodel builds its table with
#' \code{rbind()} over unnamed vectors
#' (\code{c(bw, aic.rss[2], aic.rss[3], aic.rss[1])} in \code{Model.selection.r}),
#' so a table that comes from it, or from an injected engine that copies it,
#' carries no column names and is read by position.
#' Reading column 2 ranks on AIC while labelling the result AICc, which selects
#' larger models than AICc would and can carry a pure-noise predictor.
#'
#' @param gwr_df The second element of the \code{gwr.model.selection()} return.
#' @param criterion Column name to look for, matched case-insensitively.
#' @return A list with \code{values} (numeric, non-finite coerced to
#'   \code{NA}), \code{raw} (the same before that coercion, so that an
#'   infinite criterion can be told from a missing one), \code{column},
#'   \code{column_name}, \code{by_name},
#'   \code{shape_ok} (\code{TRUE} when the table has the documented four
#'   columns, or when the column was located by name) and \code{label}.
#' @keywords internal
#' @noRd
.gwr_ms_criterion <- function(gwr_df, criterion = "AICc") {
  if (is.null(gwr_df))
    stop(".gwr_ms_criterion(): the engine returned no diagnostic table.",
         call. = FALSE)

  # A bare numeric vector (no dim) is unambiguous: it IS the criterion.
  if (is.atomic(gwr_df) && is.null(dim(gwr_df))) {
    vals <- suppressWarnings(as.numeric(gwr_df))
    raw_vals <- vals
    vals[!is.finite(vals)] <- NA_real_
    return(list(values = vals, raw = raw_vals, column = 1L,
                column_name = NA_character_,
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
  raw_vals <- vals
  vals[!is.finite(vals)] <- NA_real_

  col_name <- if (!is.null(cn) && nzchar(cn[[col]])) cn[[col]] else NA_character_
  # Say which column was actually read.  Claiming "unlabelled" when the table
  # IS labelled -- just not with the name we wanted -- hides the fact that the
  # ranking may be built on a different criterion entirely.
  label <- if (by_name) col_name
    else if (!is.na(col_name))
      sprintf("%s (assumed: column %d, labelled '%s')", criterion, col, col_name)
    else sprintf("%s (assumed: column %d, unlabelled)", criterion, col)

  list(values = vals, raw = raw_vals, column = col, column_name = col_name,
       by_name = by_name, shape_ok = shape_ok, label = label)
}


#' The uncorrected-AIC column of GWmodel's diagnostic table, or NULL
#'
#' Needed to tell a defined AICc from one outside its domain (see
#' \code{.gwr_aicc_undefined()}).  Found by name in a labelled table, or as
#' column 2 of the documented four-column \code{c(bandwidth, AIC, AICc, RSS)}
#' layout when AICc was read positionally from column 3.  \code{NULL} for any
#' other shape, so an unrecognised table is never guarded on a guessed column.
#'
#' @param gwr_df The diagnostic table.
#' @param crit What \code{.gwr_ms_criterion()} returned for it.
#' @return Numeric vector, one value per model, or \code{NULL}.
#' @keywords internal
#' @noRd
.gwr_ms_aic_column <- function(gwr_df, crit) {
  if (is.null(dim(gwr_df))) return(NULL)
  is_df <- is.data.frame(gwr_df)
  cn    <- if (is_df) names(gwr_df) else colnames(gwr_df)
  col   <- NA_integer_
  if (!is.null(cn)) {
    hit <- which(toupper(trimws(cn)) == "AIC")
    if (length(hit)) col <- as.integer(hit[[1L]])
  }
  if (is.na(col) && !isTRUE(crit$by_name) && isTRUE(crit$shape_ok) &&
      identical(as.integer(crit$column), 3L))
    col <- 2L
  if (is.na(col)) return(NULL)
  suppressWarnings(as.numeric(if (is_df) gwr_df[[col]] else gwr_df[, col]))
}


#' Normalise GWmodel's model list into character vectors of predictors
#'
#' Handles the three shapes the element can take: a two-element list of the
#' model's formula \emph{string} and its predictor vector, a formula object, or
#' a bare character vector.  The two-element list is what GWmodel returns,
#' \code{list(Generate.formula(DeVar, vars), vars)}, i.e.
#' \code{list("z ~ a + b", c("a", "b"))}.
#'
#' @param model_list The first element of the \code{gwr.model.selection()}
#'   return.
#' @param response_var Response name, stripped when it appears in the element.
#' @return A list of character vectors.
#' @keywords internal
#' @noRd
.gwr_ms_varsets <- function(model_list, response_var) {
  if (!is.list(model_list) || length(model_list) == 0L)
    stop(".gwr_ms_varsets(): the engine returned an empty model list.",
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
#'   \code{criterion} and \code{delta}, the criterion's distance from the
#'   best model's (0 in the first row, \code{NA} where the criterion is).  A
#'   model without predictors is labelled \code{<none>}.
#' @keywords internal
#' @noRd
.gwr_ms_table <- function(varsets, values, minimise = TRUE) {
  n <- length(varsets)
  if (length(values) != n)
    stop(sprintf(paste0(".gwr_ms_table(): the engine returned %d model(s) but ",
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
  # How far each model is from the best: the gap a reader would otherwise
  # subtract by hand before calling one model selected.
  df$delta <- if (all(is.na(df$criterion))) NA_real_
    else if (minimise) df$criterion - min(df$criterion, na.rm = TRUE)
    else max(df$criterion, na.rm = TRUE) - df$criterion
  rownames(df) <- NULL
  attr(df, "varsets") <- varsets[ord]
  df
}


#' Run an expression, optionally swallowing output written with cat()
#'
#' GWmodel's functions report progress with bare \code{cat()}, which neither
#' \code{suppressMessages()} nor \code{suppressWarnings()} touches; a call
#' to one of them (in a test that compares against GWmodel, say) is wrapped
#' in this to keep the console clean.
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


#' Run the forward model-selection sweep
#'
#' The sweep \code{GWmodel::gwr.model.selection()} defines, scored by the
#' package's own local solve.  Starting from no predictor, each step fits the
#' current model plus each remaining candidate in turn, keeps the candidate
#' with the lowest AICc, and goes on until every candidate is in:
#' \code{p (p + 1) / 2} models for \code{p} candidates, each at the one
#' bandwidth.  A model whose AICc is outside its domain is passed over when a
#' step is chosen.  The model without predictors (a local mean) is scored at
#' the same bandwidth.
#'
#' Everything that scores a model lives here, so the rest of
#' \code{gwr_model_selection()} is testable with an injected engine.
#'
#' @return A list with \code{model_list}, \code{gwr_df}, \code{bandwidth},
#'   \code{bandwidth_source}, \code{used_dmat} and \code{raw} (the sweep's
#'   model list and diagnostic table in the layout
#'   \code{GWmodel::gwr.model.selection()} returns, which the caller passes
#'   straight through to the public result).  \code{.engine} is a documented
#'   injection point, so this list is the contract an injected engine has to
#'   satisfy, \code{raw} included.  Two further elements are optional:
#'   \code{null_model} (a numeric vector named \code{bandwidth}, \code{AIC},
#'   \code{AICc} and \code{RSS} for the model without predictors; without it
#'   the table has no \code{<none>} row), \code{bandwidth_search} (what
#'   \code{.gwr_bw_search()} returned), and \code{undefined} with
#'   \code{null_undefined} (for each model of the sweep, and for the one
#'   without predictors, the reason \code{.gwr_criterion()} gives why its
#'   criterion is not usable, \code{""} when it is).
#' @keywords internal
#' @noRd
.gwr_ms_engine <- function(dat, response_var, candidate_vars, bandwidth,
                           adaptive, kernel, bw_approach, dmat_max_n,
                           quiet = TRUE) {
  # GWmodel is not needed: the sweep runs on the package's own local solve.
  if (!requireNamespace("sp", quietly = TRUE))
    stop("gwr_model_selection(): package 'sp' is required (the sweep hands ",
         "the data to the local solve as a Spatial object). Install with ",
         "install.packages('sp').", call. = FALSE)

  n_obs   <- nrow(dat)
  n_cand  <- length(candidate_vars)
  sp_dat  <- .to_sp(dat, unique(c(response_var, candidate_vars)))
  df      <- sp_dat@data
  coords  <- sp::coordinates(sp_dat)
  y       <- as.numeric(df[[response_var]])
  # The design of one model: the intercept and its predictors' columns.  The
  # candidates are numeric (gwr_model_selection() refuses anything else), so
  # this is the model matrix.
  design <- function(vars)
    cbind(Intercept = 1, as.matrix(df[, vars, drop = FALSE]))

  # --- Bandwidth ------------------------------------------------------------
  # One bandwidth is used for every candidate model: that is what makes the
  # criterion comparable across them, and it is also the method's main
  # limitation (see the function's documentation).
  bandwidth_source <- "supplied"
  bw_search <- NULL
  if (is.null(bandwidth)) {
    # The grid-then-refine search fit_gwr_model() uses, on the model with
    # every candidate; see .gwr_bw_search().
    bw_search <- tryCatch(
      .gwr_bw_search(design(candidate_vars), y, coords, kernel, adaptive,
                     approach = bw_approach),
      error = function(e) {
        .log_warn("gwr_model_selection(): the bandwidth search failed: %s",
                  conditionMessage(e))
        NULL
      }
    )
    bw <- if (is.null(bw_search)) NA_real_ else bw_search$bandwidth
    if (!isTRUE(is.finite(bw)) || bw <= 0) {
      bw <- .fallback_bandwidth(sp_dat, adaptive)
      bandwidth_source <- "fallback"
      warning(sprintf(paste0("gwr_model_selection(): automatic bandwidth ",
                             "selection failed; using arbitrary fallback ",
                             "bandwidth = %.4f. Supply `bandwidth` explicitly ",
                             "for a defensible comparison."), bw),
              call. = FALSE)
    } else if (isTRUE(bw_search$none_defined)) {
      bandwidth_source <- "widest (no bandwidth has a defined criterion)"
      .warn_and_log(paste0("gwr_model_selection(): no bandwidth from %s to ",
                           "%s%s has a defined %s for the model with every ",
                           "candidate, so none could be selected: at each of ",
                           "the %d tried %s. The widest, %s, is used. ",
                           "Supply `bandwidth`, or use fewer candidates or ",
                           "more observations."),
                    format(signif(bw_search$range[1L], 6)),
                    format(signif(bw_search$range[2L], 6)),
                    if (isTRUE(adaptive)) " neighbours" else "",
                    .gwr_bw_criterion_name(bw_search, score = TRUE),
                    bw_search$n_evaluations,
                    .gwr_bw_undefined_text(bw_search, n_obs),
                    format(signif(bw, 6)))
    } else {
      # Where locations repeat the search leaves out whole locations, for
      # either criterion asked for; see .gwr_bw_search().
      bandwidth_source <- sprintf(
        "%s search on the full model",
        if (isTRUE(bw_search$by_location)) "leave-one-location-out CV"
        else bw_approach)
      note <- .gwr_bw_search_note(bw_search, adaptive)
      if (nzchar(note))
        .log_warn(paste0("gwr_model_selection(): bandwidth %s%s selected by ",
                         "%s from %d scored between %s and %s: %s See ",
                         "$bandwidth_search."),
                  format(signif(bw, 6)),
                  if (isTRUE(adaptive)) " neighbours" else "",
                  .gwr_bw_criterion_name(bw_search),
                  bw_search$n_evaluations,
                  format(signif(bw_search$range[1L], 6)),
                  format(signif(bw_search$range[2L], 6)), note)
    }
  } else {
    bw <- as.numeric(bandwidth)
  }

  # Clamp against the LARGEST model in the sweep (all candidates + intercept),
  # not the smallest.  A bandwidth that supports a one-variable fit can leave
  # the full model underdetermined.
  if (isTRUE(adaptive)) {
    bw      <- as.integer(round(bw))
    if (n_cand + 2L > n_obs)
      stop(sprintf(paste0("gwr_model_selection(): %d observations cannot ",
                          "support a local fit of %d candidate predictors ",
                          "plus an intercept."), n_obs, n_cand), call. = FALSE)
    # The full model has n_cand + 1 parameters.  n_cand + 2 left bisquare and
    # tricube windows with n_cand + 1 points (their bw-th neighbour has weight
    # 0), an exact interpolator whose AICc of about -n^2/2 put it first; see
    # .gwr_min_adaptive_bw().  Capped at n_obs, as the stop above allows, and
    # the AICc guard in gwr_model_selection() deals with what is left.
    min_bw  <- min(.gwr_min_adaptive_bw(n_cand + 1L, kernel), n_obs)
    if (bw < min_bw) {
      # A warning: every model in the sweep is fitted at a bandwidth the
      # caller did not ask for.
      .warn_and_log(paste0("gwr_model_selection(): an adaptive bandwidth of ",
                           "%d neighbours is too small for the full ",
                           "%d-predictor model with the %s kernel; using %d.%s"),
                    bw, n_cand, kernel, min_bw,
                    if (kernel %in% c("bisquare", "tricube"))
                      paste0(" That is enough unless several neighbours tie ",
                             "at the kernel's edge (a regular grid), which ",
                             "gives them weight 0 too; then use a larger ",
                             "bandwidth.") else "")
      bw <- min_bw
    }
    # Capped at n, and said: a supplied count above n is most often a
    # distance meant for adaptive = FALSE.  See fit_gwr_model().
    if (bw > n_obs) {
      if (identical(bandwidth_source, "supplied"))
        .warn_and_log(paste0("gwr_model_selection(): an adaptive bandwidth of ",
                             "%d neighbours exceeds the %d observations; using ",
                             "%d. With adaptive = TRUE the bandwidth is a count ",
                             "of neighbours; for a distance in CRS units, set ",
                             "adaptive = FALSE."),
                      bw, n_obs, n_obs)
      bw <- as.integer(n_obs)
    }
  }

  # A supplied fixed bandwidth in the wrong units (degrees read as metres)
  # leaves every window its own observation alone; fit_gwr_model() says what
  # is wrong, and so does this.
  if (!isTRUE(adaptive) && identical(bandwidth_source, "supplied"))
    .gwr_warn_tiny_fixed_bw(dat, bw, "gwr_model_selection")

  # --- The sweep --------------------------------------------------------------
  # c(bandwidth, AIC, AICc, RSS) for one model, the row GWmodel's table has,
  # and the value a forward step is chosen on: AICc, or Inf outside its domain.
  # `why` is the reason a model that the formula gives a number for is still
  # not comparable (.gwr_criterion()): its fit reproduces the mean of the
  # observations at each location, or a term of it is estimable nowhere.
  loc <- .gwr_locations(coords)
  score <- function(vars) {
    cr <- .gwr_criterion(design(vars), y, coords, bw, kernel, adaptive,
                         approach = "AICc", loc = loc)
    list(row = c(bw, cr$AIC, cr$AICc, cr$RSS), value = cr$value,
         n_deficient = cr$n_deficient, why = cr$why)
  }
  n_sweep    <- n_cand * (n_cand + 1L) / 2L
  model_list <- vector("list", n_sweep)
  # GWmodel's table has these four columns and no names; they are named here,
  # so the criterion is found by name and not by position.
  gwr_df     <- matrix(NA_real_, n_sweep, 4L,
                       dimnames = list(NULL, c("bandwidth", "AIC", "AICc",
                                               "RSS")))
  n_def      <- integer(n_sweep)
  why        <- character(n_sweep)
  level      <- character(0)
  remaining  <- candidate_vars
  tag <- 0L
  for (i in seq_len(n_cand)) {
    step <- numeric(length(remaining))
    for (j in seq_along(remaining)) {
      vars <- c(level, remaining[j])
      tag  <- tag + 1L
      fml  <- paste0(response_var, "~", paste(vars, collapse = "+"))
      # Progress, for a long sweep watched with quiet = FALSE.
      if (!isTRUE(quiet))
        cat(sprintf("Now calibrating model %d of %d: %s\n", tag, n_sweep, fml))
      sc <- score(vars)
      model_list[[tag]] <- list(fml, vars)
      gwr_df[tag, ]     <- sc$row
      n_def[tag]        <- sc$n_deficient
      why[tag]          <- sc$why
      step[j]           <- sc$value
    }
    idx <- if (any(is.finite(step))) which.min(step) else 1L
    level     <- c(level, remaining[idx])
    remaining <- remaining[-idx]
  }
  if (!isTRUE(quiet))
    cat(sprintf("Now calibrating the model without predictors: %s~1\n",
                response_var))
  null <- score(character(0))

  # A model with windows of deficient rank is scored on its estimable fitted
  # values, which is sound, but its coefficients are NA there: say how many
  # models that concerns before one of them is selected and fitted.
  if (any(n_def > 0L, na.rm = TRUE))
    .log_warn(paste0("gwr_model_selection(): %d of %d model(s) have local ",
                     "windows of deficient rank at bandwidth %s (up to %d of ",
                     "%d windows). They are scored on their fitted values, ",
                     "which those windows do estimate; fit_gwr_model() on such ",
                     "a model returns NA for the coefficients they cannot."),
              sum(n_def > 0L, na.rm = TRUE), n_sweep, format(bw),
              max(n_def, na.rm = TRUE), n_obs)

  raw <- list(model_list, gwr_df)
  list(model_list = model_list, gwr_df = gwr_df,
       bandwidth = bw, bandwidth_source = bandwidth_source,
       used_dmat = FALSE, raw = raw,
       null_model = c(bandwidth = as.numeric(bw), AIC = null$row[2L],
                      AICc = null$row[3L], RSS = null$row[4L]),
       bandwidth_search = bw_search,
       undefined = why, null_undefined = null$why)
}


# -----------------------------------------------------------------------------
# Exported
# -----------------------------------------------------------------------------

#' Forward model selection for geographically weighted regression
#'
#' Grows a GWR model one predictor at a time and scores every intermediate
#' model with a corrected Akaike information criterion, the sweep
#' \code{GWmodel::gwr.model.selection()} defines, and returns a ranked table
#' with each model's distance from the best.  The model without predictors is
#' scored too.  Every model is fitted by the package's own local solve (see
#' \code{\link{fit_gwr_model}}, \strong{The local solve}), so a candidate
#' that is constant inside some windows, a regional indicator say, is scored
#' on what those windows do estimate and neither stops the sweep nor scores
#' as rounding error.
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
#' disagreement is itself informative. It usually means a candidate is
#' predictive only locally.
#'
#' One further limitation follows from the method itself:
#'
#' \itemize{
#'   \item \strong{One bandwidth for every model.} Comparing criteria across
#'     models requires holding the smoothing fixed, but the bandwidth is itself
#'     a fitted quantity, and the value chosen for the full model is not
#'     optimal for a one-predictor model. This is how the method is defined
#'     (Lu et al. 2014); it is not an implementation shortcut. Refit the
#'     selected model with \code{bandwidth = NULL} to re-optimise once the
#'     variable set is settled.
#' }
#'
#' @section The model without predictors, and the gaps:
#' The forward sweep starts from one predictor, so by itself it always names
#' at least one, whether or not any candidate helps.  The model without
#' predictors, a local mean at the sweep's bandwidth, is therefore scored as
#' well and ranked with the others as \code{<none>}.  When it ranks first,
#' \code{best} is \code{character(0)}: no candidate improves on it.
#'
#' The \code{delta} column of \code{table} is each model's criterion minus
#' the best one's.  Differences of AICc below about 2 do not separate two
#' models (Burnham and Anderson 2002), so \code{within_2} lists the models
#' that close to the best and \code{smallest_within_2} names the one among
#' them with the fewest predictors, which is the one to prefer when several
#' are that close.  \code{best} is still the model with the lowest criterion.
#'
#' @section Observations that share a location:
#' When rows of \code{data_sf} have the same coordinates, the bandwidth
#' search leaves out whole locations, for either \code{bw_approach} (see
#' \code{\link{fit_gwr_model}}, \strong{Observations that share a
#' location}); \code{bandwidth_source} then reads
#' \code{"leave-one-location-out CV search on the full model"}.  The models
#' are still ranked by AICc, which counts the rows at one location as
#' independent, so the gaps overstate what a predictor adds when a location
#' has an effect of its own; a warning of class
#' \code{"spatialkit_gwr_repeated_locations"} says so with the counts.  A
#' model whose fit reproduces the mean of the observations at each location
#' (every such model scores the same, whatever its predictors) has criterion
#' \code{NA} and is listed last, as is one with a term that no local window
#' can estimate; a warning names them.
#'
#' @section Cost:
#' The sweep fits \code{p * (p + 1) / 2} GWR models for \code{p} candidates
#' (55 at p = 10, 210 at p = 20), and the one without predictors, each over
#' all \code{n} locations.
#' \code{max_models} stops the call before it runs for hours.
#'
#' @param data_sf An \code{sf} object with response, predictors and geometry.
#' @param response_var Response column name.
#' @param candidate_vars Character vector naming at least two \strong{numeric}
#'   predictors to choose among.  Factor, character and logical candidates are
#'   refused: a factor is several model-matrix columns while this
#'   sweep counts it as one variable, so the criteria would not be comparable,
#'   and \code{\link{fit_gwr_model}()}, the documented next step, takes only
#'   numerics.  Encode them as numeric indicators first.
#' @param bandwidth Bandwidth held fixed across all candidate models. If
#'   \code{NULL} (default) it is searched for on the model containing every
#'   candidate, by the grid-then-refine search of
#'   \code{\link{fit_gwr_model}()} (see its \strong{Bandwidth selection});
#'   what the search scored is returned in \code{bandwidth_search}. Integer
#'   neighbour count when
#'   \code{adaptive = TRUE}; otherwise a distance in the units of the
#'   **projected** CRS the sweep runs in, which \code{prep_model_data()} may
#'   have chosen for you.  Geographic input is projected before the bandwidth
#'   is used, so a value in degrees would be read as metres; a fixed
#'   bandwidth below a ten-thousandth of the data's extent raises a warning
#'   saying so, as in \code{\link{fit_gwr_model}()}.  An adaptive
#'   count too small for the full model is raised, with a warning, to the
#'   number of candidates plus 3 for the bisquare and tricube kernels (which
#'   give the farthest neighbour in a window weight 0), plus 2 for the others.
#'   That floor is enough unless several neighbours tie at the kernel's edge
#'   (a regular grid), which leaves a window fewer weighted points; then use
#'   a larger bandwidth.  An adaptive count below 1 or above R's largest
#'   integer is refused, as in \code{\link{fit_gwr_model}()}.
#'   One above the number of observations is capped at it, with a warning.
#' @param adaptive Logical; adaptive (nearest-neighbour) bandwidth. Default
#'   \code{TRUE}.
#' @param kernel Weighting kernel. One of \code{"bisquare"} (default),
#'   \code{"gaussian"}, \code{"tricube"}, \code{"boxcar"},
#'   \code{"exponential"}.
#' @param bw_approach Criterion for the bandwidth search: \code{"AICc"}
#'   (default) or \code{"CV"}, the leave-one-out sum of squared prediction
#'   errors. Matches \code{\link{fit_gwr_model}}'s default.  When
#'   observations share a location the search leaves out whole locations
#'   whichever is asked for; see \strong{Observations that share a
#'   location}.
#' @param max_models Refuse the call if the sweep would exceed this many model
#'   fits. Default 200, which admits up to 19 candidates.
#' @param dmat_max_n Not used.  The sweep takes its distances window by
#'   window and builds no \code{n} x \code{n} distance matrix, so a fixed
#'   bandwidth works at any \code{n}.  The argument is accepted so that calls
#'   written with it run.
#' @param quiet Print no progress. Default \code{TRUE}. With \code{FALSE} one
#'   line is written (with \code{cat()}) as each model is fitted, to watch a
#'   long sweep progress.
#' @param .engine Internal; injectable backend used for testing.
#'
#' @return An object of class \code{gwr_model_selection}, a list with:
#'   \code{best} (character vector of the selected predictors, the model with
#'   the lowest criterion; \code{character(0)} when that is the model without
#'   predictors);
#'   \code{table} (ranked data.frame of every model evaluated, with columns
#'   \code{rank}, \code{n_vars}, \code{variables}, \code{criterion} and
#'   \code{delta}, the criterion minus the best model's; the model without
#'   predictors is the row \code{<none>};
#'   \code{criterion} is \code{NA}, and the model ranked last, where it
#'   could not be evaluated or where AICc is undefined because the model's
#'   effective number of parameters \eqn{tr(S)} is not below \eqn{n - 2},
#'   which raises a warning);
#'   \code{within_2} (the rows of \code{table} whose \code{delta} is at most
#'   2) and \code{smallest_within_2} (the predictors of the smallest model
#'   among them);
#'   \code{criterion} (label for the criterion actually read, noting when it
#'   had to be located positionally);
#'   \code{criterion_by_name} (logical: whether that column was found by
#'   its name rather than by the documented position),
#'   \code{criterion_column} (the column it was read from) and
#'   \code{criterion_verified} (logical: \code{FALSE} exactly when the
#'   column was read positionally from a table that did not have the four
#'   documented columns, which is the case the log calls unverified.  Gate
#'   a script on this field instead of on the label);
#'   \code{response_var} and \code{candidate_vars} (the response and the full
#'   candidate set the sweep ran over, both echoed by \code{print()});
#'   \code{bandwidth}, \code{bandwidth_source}, \code{adaptive} and
#'   \code{kernel} (the smoothing held fixed across the sweep, and where it
#'   came from: \code{"supplied"}, or the criterion's name followed by
#'   \code{"search on the full model"}, the criterion being
#'   \code{"leave-one-location-out CV"} when observations share a location);
#'   \code{bandwidth_search} (what the search scored, as
#'   \code{\link{fit_gwr_model}} keeps it; \code{NULL} when the bandwidth was
#'   supplied);
#'   \code{n_obs}, \code{n_models} (the rows of \code{table}),
#'   \code{used_dmat} (always \code{FALSE}); and \code{raw}
#'   (the forward sweep in the layout \code{GWmodel::gwr.model.selection()}
#'   returns: the two-element list of its model list and its diagnostic
#'   table, without the model that has no predictors).
#'
#' @section Using $raw with GWmodel directly:
#' \code{raw} has the layout of GWmodel's own \code{list(model.list, GWR.df)},
#' so its two elements have to be unpacked before GWmodel's helpers will take
#' them: \code{GWmodel::gwr.model.view()} takes
#' \code{(DeVar, InDeVars, model.list)}, so the call is
#' \preformatted{
#'   GWmodel::gwr.model.view(sel$response_var, sel$candidate_vars, sel$raw[[1]])
#' }
#' Pass \code{sel$raw[[1]]}, not \code{sel$raw}.  The diagnostic table is
#' \code{sel$raw[[2]]}, a numeric matrix whose columns are
#' \code{bandwidth}, \code{AIC}, \code{AICc}, \code{RSS} in that order (and
#' are named so); the \code{criterion} column of \code{$table} is its third
#' column.
#'
#' @references
#' Burnham, K. P. and Anderson, D. R. (2002). \emph{Model Selection and
#' Multimodel Inference: A Practical Information-Theoretic Approach}, 2nd
#' edition. Springer.
#'
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
#' if (requireNamespace("sp", quietly = TRUE)) {
#'   library(sf)
#'   set.seed(1)
#'   n <- 80
#'   dat <- st_as_sf(
#'     data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000),
#'                a = rnorm(n), b = rnorm(n), noise = rnorm(n)),
#'     coords = c("x", "y"), crs = 32632
#'   )
#'   dat$z <- 2 * dat$a - dat$b + rnorm(n, 0, 0.5)
#'   sel <- gwr_model_selection(dat, "z", c("a", "b", "noise"), bandwidth = 30)
#'   print(sel)              # every model tried, by AICc, with the gaps
#'   print(sel$best)         # the winning predictor set: a and b, not noise
#'   fit <- fit_gwr_model(dat, "z", sel$best)
#'   fit
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
  # value never gets past it.
  kernel      <- match.arg(kernel)
  bw_approach <- match.arg(bw_approach)

  # Validate `bandwidth` here, not inside .gwr_ms_engine(): the engine is an
  # injectable test seam, so a check placed there would be skipped for
  # injected engines.
  # Unvalidated, NA reaches the engine's clamp as "missing value where
  # TRUE/FALSE needed", a length-2 vector as "the condition has length > 1",
  # and with adaptive = FALSE (no clamp at all) a zero or negative distance
  # goes straight into the kernel.
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
  # The same range check fit_gwr_model() applies to an adaptive count.
  # Without it a count above .Machine$integer.max (a distance passed with
  # adaptive left at TRUE) became NA in the engine's as.integer() and the
  # clamp stopped with a bare "missing value where TRUE/FALSE needed", and a
  # count below 1 was rounded to "0 neighbours" and raised, where
  # fit_gwr_model() refuses it.
  if (!is.null(bandwidth) && isTRUE(adaptive))
    .check_scalar(bandwidth, "bandwidth", "gwr_model_selection", min = 1,
                  max = .Machine$integer.max,
                  what = "a single number of nearest neighbours when adaptive = TRUE")

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

  # The documented next step is fit_gwr_model(dat, response, sel$best), and
  # fit_gwr_model() refuses every non-numeric predictor -- so a factor or
  # character candidate got scored here, could be returned in $best, and then
  # errored on the very next line of the example.  A factor is also counted as
  # ONE variable while contributing several model-matrix columns, so its
  # criterion is not comparable with the others'.  Refuse it where the sweep
  # starts.
  # is.numeric ONLY, matching fit_gwr_model()'s own guard.  Accepting logicals
  # here let the sweep return a candidate in $best that the documented next
  # step -- the last line of this function's own example -- then refused.
  cand_types <- vapply(sf::st_drop_geometry(dat)[candidate_vars],
                       is.numeric, logical(1))
  if (any(!cand_types))
    stop(sprintf(paste0("gwr_model_selection(): candidate(s) %s are not numeric. ",
                        "A factor is several model-matrix columns ",
                        "while this sweep counts it as one variable, so the ",
                        "criteria would not be comparable -- and fit_gwr_model(), ",
                        "the documented next step, rejects it outright. Encode ",
                        "them as numeric indicators first (a logical column ",
                        "needs as.numeric() too -- fit_gwr_model() takes only ",
                        "numerics)."),
                 paste(sprintf("'%s'", candidate_vars[!cand_types]),
                       collapse = ", ")), call. = FALSE)

  eng <- .engine(dat = dat, response_var = response_var,
                 candidate_vars = candidate_vars, bandwidth = bandwidth,
                 adaptive = adaptive, kernel = kernel,
                 bw_approach = bw_approach, dmat_max_n = dmat_max_n,
                 quiet = quiet)

  varsets <- .gwr_ms_varsets(eng$model_list, response_var)
  crit    <- .gwr_ms_criterion(eng$gwr_df, criterion = "AICc")
  aic <- .gwr_ms_aic_column(eng$gwr_df, crit)
  # The model without predictors, when the engine scored it: a local mean at
  # the sweep's bandwidth, ranked with the others as `<none>`.  Without it the
  # sweep starts from one predictor and can only ever name at least one,
  # whether or not any candidate helps.
  null <- eng$null_model
  has_null <- is.numeric(null) && !is.null(names(null)) &&
    all(c("AIC", "AICc") %in% names(null))
  if (has_null) {
    varsets     <- c(list(character(0)), varsets)
    crit$values <- c(if (is.finite(null[["AICc"]])) null[["AICc"]] else NA_real_,
                     crit$values)
    crit$raw    <- c(null[["AICc"]], crit$raw)
    if (!is.null(aic)) aic <- c(null[["AIC"]], aic)
  }
  # A model whose AICc is outside its domain (tr S >= n - 2: its local fits
  # interpolate the data) gets a large negative AICc from the formula, minus
  # infinity when the residual sum of squares is exactly 0, and the sweep
  # ranked it first: at n = 60 a model of `a` plus three noise variables won
  # with -66689.  Set those to NA so they rank last, and say so.  Judged on
  # the criterion as the engine returned it, so that an infinite one counts.
  # Unguarded when the AIC column cannot be located.
  if (!is.null(aic) && length(aic) == length(crit$values) &&
      length(crit$raw) == length(crit$values)) {
    undef <- !is.na(crit$raw) & .gwr_aicc_undefined(aic, crit$raw)
    if (any(undef)) {
      crit$values[undef] <- NA_real_
      # No larger adaptive bandwidth exists once it is every observation.
      at_n <- isTRUE(adaptive) &&
        isTRUE(suppressWarnings(as.numeric(eng$bandwidth)) >= nrow(dat))
      .warn_and_log(paste0("gwr_model_selection(): AICc is undefined for %d ",
                           "of %d model(s) at bandwidth %s and they are ranked ",
                           "last: their effective number of parameters tr(S) ",
                           "is not below n - 2 = %d, so their local fits ",
                           "(nearly) interpolate the data. A forward step is ",
                           "chosen among the models whose AICc is defined. %s"),
                    sum(undef), length(undef), format(eng$bandwidth),
                    nrow(dat) - 2L,
                    if (at_n)
                      sprintf(paste0("Even the widest adaptive window (all %d ",
                                     "observations) leaves too few residual ",
                                     "degrees of freedom for them: use fewer ",
                                     "candidates, more observations, or a ",
                                     "gaussian or exponential kernel."),
                              nrow(dat))
                    else "Use a larger bandwidth.")
    }
  }
  # Two more kinds of model have a number from the formula and still cannot
  # be ranked, and the engine says which (`undefined`).  One reproduces the
  # mean of the observations at each location: on a layer whose locations
  # repeat, every such model scores the same, whatever its predictors.  The
  # other has a term that no window estimates, so it has no coefficient for
  # it anywhere.  Both rank last.
  why <- eng$undefined
  if (is.character(why) && length(why) == length(eng$model_list)) {
    if (has_null) why <- c(as.character(eng$null_undefined %||% "")[1L], why)
    labels <- vapply(varsets, function(v)
      if (length(v)) paste(v, collapse = " + ") else "<none>", character(1))
    sat  <- why == "df_distinct" & !is.na(crit$values)
    term <- why == "term" & !is.na(crit$values)
    if (any(sat) || any(term)) {
      crit$values[sat | term] <- NA_real_
      .warn_and_log(paste0("gwr_model_selection(): %d of %d model(s) are not ",
                           "ranked at bandwidth %s; their criterion is NA and ",
                           "they are listed last.%s%s Use a larger bandwidth."),
                    sum(sat | term), length(why), format(eng$bandwidth),
                    if (any(sat))
                      sprintf(paste0(" %d (%s) reproduce the mean of the ",
                                     "observations at each location, or ",
                                     "nearly: their effective number of ",
                                     "parameters tr(S) is not below the ",
                                     "number of distinct rows of location ",
                                     "and predictors minus 2, so AICc ",
                                     "measures only the spread among ",
                                     "observations that share a location."),
                              sum(sat), paste(labels[sat], collapse = "; "))
                    else "",
                    if (any(term))
                      sprintf(paste0(" %d (%s) have a term that no local ",
                                     "window can estimate, so no ",
                                     "coefficient for it exists anywhere."),
                              sum(term), paste(labels[term], collapse = "; "))
                    else "")
    }
  }
  # Several observations at one location: said once, with the counts.
  loc <- .gwr_locations(sf::st_coordinates(dat)[, 1:2, drop = FALSE])
  if (loc$n < nrow(dat))
    .gwr_warn(
      "spatialkit_gwr_repeated_locations",
      paste0("gwr_model_selection(): the %d observations are at %d distinct ",
             "location(s) (%d hold more than one observation, the most %d). ",
             "AICc counts the ",
             "observations at one location as independent. They are not ",
             "when a location has an effect of its own, and the gaps in the ",
             "table then overstate what a predictor adds.%s Confirm the ",
             "selection with select_features_forward(), or fit the means of ",
             "replicate observations, one row per location."),
      nrow(dat), loc$n, sum(loc$size > 1L), max(loc$size),
      if (isTRUE(eng$bandwidth_search$by_location) &&
          isFALSE(eng$bandwidth_search$none_defined))
        sprintf(paste0(" No bandwidth was supplied: it was selected by %s ",
                       "(each location predicted from the others)."),
                .gwr_bw_criterion_name(eng$bandwidth_search))
      else "")
  tab     <- .gwr_ms_table(varsets, crit$values, minimise = TRUE)

  # The diagnostic table is in GWmodel's layout, a matrix of
  # c(bandwidth, AIC, AICc, RSS): AICc is column 3, and reading column 2 would
  # rank on the uncorrected AIC while labelling the answer AICc.  An unnamed
  # table of any other shape (an injected engine's) is read at the same
  # position and the ranking is flagged as unverified.
  if (!isTRUE(crit$shape_ok) && !isTRUE(crit$by_name))
    .log_warn(paste0("gwr_model_selection(): the diagnostic table does ",
                     "not have the documented four columns (bandwidth, AIC, ",
                     "AICc, RSS); AICc was read from column %d of a table with ",
                     "a shape this version does not recognise. Treat the ",
                     "ranking as unverified and check $raw."), crit$column)

  best_set <- attr(tab, "varsets")[[1L]]
  # The models a difference of 2 does not separate from the best, and the
  # smallest of them: the usual reading of AICc differences (Burnham and
  # Anderson 2002) is that such models have about the same support, and then
  # the one with fewest predictors is the one to prefer.
  close_i <- which(!is.na(tab$delta) & tab$delta <= 2)
  small_i <- close_i[order(tab$n_vars[close_i], tab$criterion[close_i])][1L]
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
      # The rows of `table` within 2 of the best, and the predictors of the
      # smallest model among them (character(0) for the model without any).
      within_2          = tab[close_i, , drop = FALSE],
      smallest_within_2 = attr(tab, "varsets")[[small_i]],
      criterion        = crit$label,
      # How the criterion column was found: by name, or by its documented
      # position (`criterion_column`), and whether the table had the four
      # columns that position is documented for.  A script can gate on
      # `criterion_verified` where the log line above said "unverified".
      criterion_by_name  = isTRUE(crit$by_name),
      criterion_column   = crit$column,
      criterion_verified = isTRUE(crit$by_name) || isTRUE(crit$shape_ok),
      response_var     = response_var,
      candidate_vars   = candidate_vars,
      bandwidth        = eng$bandwidth,
      bandwidth_source = eng$bandwidth_source,
      adaptive         = adaptive,
      kernel           = kernel,
      n_obs            = nrow(dat),
      n_models         = nrow(tab),
      used_dmat        = isTRUE(eng$used_dmat),
      bandwidth_search = eng$bandwidth_search,
      raw              = eng$raw
    ),
    class = "gwr_model_selection"
  )
}


#' Print a GWR model selection result
#'
#' Shows the forward-selection trail: the response, the candidate predictors,
#' and the top-ranked models with their criterion values and, in the
#' \code{delta} column, how far each is from the best.  Below the table it
#' lists the models within 2 of the best, which a difference in AICc does not
#' separate (Burnham and Anderson 2002), and names the smallest of them.  When
#' several models are that close the ranking is not well identified and the
#' choice of predictors should not be treated as settled; when the model
#' without predictors (\code{<none>}) is among them, no candidate has been
#' shown to help.
#'
#' @param x A \code{gwr_model_selection} object.
#' @param n Number of top-ranked models to show. Default 10.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @references
#' Burnham, K. P. and Anderson, D. R. (2002). \emph{Model Selection and
#' Multimodel Inference: A Practical Information-Theoretic Approach}, 2nd
#' edition. Springer.
#' @family print methods
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

  rhs <- function(v) if (length(v)) paste(v, collapse = " + ") else "1"
  cat(sprintf("\nSelected: %s ~ %s\n", x$response_var, rhs(x$best)))
  if (!length(x$best))
    cat("  The model without predictors has the lowest criterion: no candidate\n",
        "  improves on a local mean at this bandwidth.\n", sep = "")

  # The gaps, worked out: which models a difference of 2 does not separate
  # from the best, and the smallest of them.  A table without `delta` (an
  # object from an earlier version) prints as it did.
  tab <- x$table
  if (is.data.frame(tab) && "delta" %in% names(tab)) {
    close <- which(!is.na(tab$delta) & tab$delta <= 2)
    if (length(close) > 1L) {
      cat(sprintf("\nWithin 2 of the best (a difference this small does not separate models):\n"))
      for (i in close)
        cat(sprintf("  %s ~ %s  (+%.2f)\n", x$response_var,
                    if (identical(tab$variables[i], "<none>")) "1"
                    else tab$variables[i], tab$delta[i]))
      small <- close[order(tab$n_vars[close], tab$criterion[close])][1L]
      cat(sprintf("Smallest of them: %s ~ %s\n", x$response_var,
                  if (identical(tab$variables[small], "<none>")) "1"
                  else tab$variables[small]))
    } else if (nrow(tab) > 1L && is.finite(tab$delta[2L])) {
      cat(sprintf("\nNo other model is within 2 of the best; the runner-up is %.2f away.\n",
                  tab$delta[2L]))
    }
  }
  cat("\nThe criterion is in-sample and every model shares one bandwidth.\n")
  cat("Confirm with select_features_forward() before relying on this.\n")
  invisible(x)
}
