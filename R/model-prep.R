#' Prepare and sanitize an sf dataset for spatial modeling
#'
#' Ensures point geometry, projected CRS, and removes rows with missing or
#' non-finite values in modeling columns, \emph{including} rows whose
#' geometry is empty or whose coordinates are not finite, which no model
#' backend can use.  All non-POINT geometries (including MULTIPOINT) are
#' coerced to representative points via \code{coerce_to_points()}, so
#' downstream coordinate extraction always aligns one row per observation.
#' Any Z or M coordinate (POINT Z from a GPS, a GeoPackage or KML) is
#' dropped, because every backend works in 2-D map distance.
#'
#' The response may not appear in \code{predictor_vars}.  Using it as its own
#' predictor is leakage no backend catches: an out-of-bag R^2 near 1 in the
#' random forest, a silently reduced design matrix in GWR, duplicated rows and
#' a phantom \code{<none>} entry in the GWR selection table.  It is refused
#' here.
#'
#' Column names must be syntactically valid R names (\code{make.names(x) == x}).
#' Every backend builds a model formula from these names.  A name R parses as
#' an expression would fit a different model from the one requested while the
#' fit object still recorded the name you asked for: \code{"B5-B4"} is
#' \code{B5 - B4}, and \code{"log(a)"} is a function call.  Rename the column
#' (for example with \code{make.names()}) before fitting.
#'
#' @param data_sf An sf object.
#' @param response_var Response variable column name.
#' @param predictor_vars Predictor column names.  May be \code{character(0)}
#'   for an intercept-only model (\code{fit_bayesian_spatial_model()} supports
#'   one; the GWR and random-forest backends do not and reject it themselves).
#' @param boundary Optional sf/sfc whose CRS the data are projected to.  It
#'   is not a clip: rows that fall outside it are kept.  Subset
#'   \code{data_sf} yourself (for example with \code{sf::st_filter()}) to
#'   restrict the rows.
#' @param pointize Strategy for non-point geometry coercion, passed to
#'   \code{\link{coerce_to_points}}.  One of \code{"auto"} (default),
#'   \code{"centroid"}, \code{"point_on_surface"}, \code{"surface"} (an alias
#'   for \code{"point_on_surface"}), \code{"line_midpoint"} or
#'   \code{"bbox_center"}.
#' @param require_response Logical; if FALSE the response column is not
#'   required to be present (useful for out-of-sample prediction where the
#'   response is unknown).  Default TRUE.
#' @return An sf object with 2-D (XY) POINT geometry, cleaned of rows
#'   carrying missing or non-finite values in the modelling columns or in the
#'   coordinates.  What
#'   was removed is recorded on the attribute \code{"dropped"}, a list with
#'   \code{n} (rows dropped), \code{n_geometry} (how many of them for an
#'   empty or non-finite geometry), \code{which} (their positions in
#'   \code{data_sf}), \code{row_id} (their \code{..row_id} values when the
#'   layer carries that column, else \code{NULL}) and \code{reason} (one per
#'   dropped row: \code{"geometry"}, \code{"missing"} or
#'   \code{"non_finite"}, in that order of precedence when several apply),
#'   \code{by_column} (for each modelling column, and for
#'   \code{"(geometry)"}, how many of the dropped rows it is missing or
#'   non-finite in; a row can count under several),
#'   \code{response_std_diff} (the mean response of the dropped rows minus
#'   that of the kept rows, in standard deviations of the response over both;
#'   \code{NA} when the response is not numeric, or is itself missing in
#'   every dropped row),
#'   plus \code{n_rows}, the number of rows returned, which the record was
#'   made for.  The same record carries two marks on the rows that were
#'   kept: \code{constant_predictors}, the names of the predictors that take
#'   one value in every returned row (\code{character(0)} when none does),
#'   and \code{n_duplicate_locations}, the number of returned rows whose
#'   coordinates repeat those of an earlier row.
#'   Every fit stores \code{n} as \code{$info$n_dropped}.  The record
#'   describes the rows this call returned and does not survive subsetting:
#'   \code{clean[i, ]}, like \code{dplyr::filter()}, \code{slice()} or
#'   \code{arrange()} of it, is a plain layer with no \code{"dropped"}
#'   attribute, and a fit given such a subset with \code{.already_prepped =
#'   TRUE} reports \code{n_dropped = 0} even when the parent layer dropped
#'   rows.  \code{sf::st_drop_geometry()} keeps the record, since the rows are
#'   the same; see \code{\link{[.spatialkit_rows}} for what binding such data
#'   frames does.  The CRS is
#'   projected whenever one can be established.  A CRS-less layer is decided by
#'   the lon/lat heuristic (see \code{\link{ensure_projected}}): if its bounding
#'   box fits the lon/lat envelope \emph{and} it either spans more than one unit
#'   on some axis or carries decimal-degree-like precision, it is read as
#'   EPSG:4326 and projected, with a warning.  A small planar survey inside that
#'   envelope is included in that, deliberately; only coordinates the heuristic
#'   declines are passed through as-is.  Set the CRS on \code{data_sf} if the
#'   data are planar.
#' @section Warnings:
#' When preparing data for a fit (\code{require_response = TRUE}) two
#' conditions are raised as R warnings, each with a class so it can be caught
#' by name with \code{tryCatch()} or \code{withCallingHandlers()}:
#' \describe{
#'   \item{\code{spatialkit_rows_dropped}}{Rows were dropped.  The message
#'     gives the number and share of rows, the count by reason and by
#'     column, and, for a numeric response, how far the mean response of
#'     the dropped rows sits from that of the kept rows.  A fit on the rows
#'     that remain describes those rows only; when the difference is large
#'     the rows are not missing at random and the fit is biased towards the
#'     complete ones.}
#'   \item{\code{spatialkit_constant_predictor}}{A predictor takes a single
#'     value in every row that remains.  Its coefficient cannot be told from
#'     the intercept: a GWR window is rank-deficient, a Bayesian slope is
#'     its prior, and a forest ignores it.  The predictor is kept; drop it
#'     from \code{predictor_vars}.}
#' }
#' Every fitter and every \code{cv_*()} function prepares its data here and
#' passes these on.  \code{\link{compare_models_cv}()},
#' \code{\link{cv_block_size_sweep}()} and
#' \code{\link{select_features_forward}()} run several cross-validations
#' of one layer, each of which prepares it, and raise a warning of either
#' class once per distinct text for the run.
#'
#' With \code{require_response = FALSE} (how every \code{predict()} method
#' prepares \code{newdata}) neither is raised: the dropped rows are written
#' to the log and to the record, and the prediction for each comes back
#' \code{NA}, so nothing is lost without a trace.  Repeated locations are
#' legitimate (repeat visits to one site) and are counted in the record and
#' the log, never warned about.
#' @family model fitting
#' @family spatial data preparation
#' @examples
#' library(sf)
#' dat <- st_as_sf(
#'   data.frame(x = 1:5, y = 5:1,
#'              resp = c(1, 2, NA, 4, 5),
#'              pred = c(1, 2, 3, 4, Inf)),
#'   coords = c("x", "y"), crs = 32632
#' )
#' # Drops rows 3 (NA) and 5 (Inf), and warns that it did, with the counts.
#' clean <- suppressWarnings(prep_model_data(dat, "resp", "pred"))
#' attr(clean, "dropped")
#' @export
prep_model_data <- function(data_sf, response_var, predictor_vars,
                            boundary = NULL,
                            pointize = c("auto", "surface", "point_on_surface",
                                         "centroid", "line_midpoint",
                                         "bbox_center"),
                            require_response = TRUE) {
  # call. = FALSE on every one of these, as elsewhere in the package: these are
  # argument-validation messages, and the call frame R would otherwise append
  # is an internal one the user did not write (prep_model_data() is reached
  # through fit_gwr_model(), the CV internals and every predict() method), so
  # printing it buries the message that names the actual problem.
  if (!inherits(data_sf, "sf"))
    stop("prep_model_data(): 'data_sf' must be an sf object.", call. = FALSE)
  pointize <- match.arg(pointize)
  if (!is.character(response_var) || length(response_var) != 1L)
    stop("prep_model_data(): 'response_var' must be a single column name.",
         call. = FALSE)
  # character(0) is allowed: an intercept-only spatial GP is a legitimate
  # model (see fit_bayesian_spatial_model()).  Backends that genuinely need a
  # predictor reject an empty set themselves, where the message can say why.
  if (!is.character(predictor_vars))
    stop("prep_model_data(): 'predictor_vars' must be a character vector.",
         call. = FALSE)

  # Every backend turns these names into a model formula -- reformulate() for
  # GWR, the model-selection sweep and the random forest's display formula,
  # sprintf() + as.formula() for the brms term.  R parses the result, so a
  # name that is not a syntactic name is not a name at all: "B5-B4" becomes
  # B5 - B4 (a DIFFERENT model, silently, with $predictor_vars still reporting
  # "B5-B4"), "log(a)" becomes a function call, and "band 4" dies inside
  # str2lang() with a parser error that names no column.  Backticking is not a
  # fix here: the sp coercion GWmodel needs runs the names through
  # make.names() anyway, so the formula and the data would disagree.  Refuse
  # early, name the columns, and say what to do about it.
  bad_nm <- setdiff(
    unique(c(if (require_response) response_var, predictor_vars)),
    make.names(unique(c(if (require_response) response_var, predictor_vars)))
  )
  if (length(bad_nm))
    stop("prep_model_data(): column name(s) ",
         paste(sprintf("'%s'", bad_nm), collapse = ", "),
         " are not syntactically valid R names, and every model backend builds",
         " a formula from them. Rename the column(s) -- e.g. names(x) <-",
         " make.names(names(x)) -- before fitting.", call. = FALSE)

  # The response among its own predictors is leakage in the random forest
  # (OOB R2 near 1, response top of the importance table), a silently reduced
  # model in GWR (model.matrix drops the duplicated column) and duplicated
  # rows plus a phantom "<none>" entry in the GWR selection table.  None of
  # the backends catches it, so catch it here.
  if (require_response && response_var %in% predictor_vars)
    stop("prep_model_data(): response '", response_var, "' is also listed in",
         " 'predictor_vars'. A model cannot use its own response as a",
         " predictor.", call. = FALSE)

  req_cols <- if (require_response) c(response_var, predictor_vars) else predictor_vars
  miss <- setdiff(req_cols, names(data_sf))
  if (length(miss))
    stop("prep_model_data(): missing required column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)

  # Coerce unless every geometry is already a plain POINT.  MULTIPOINT must
  # be coerced too: a multi-vertex MULTIPOINT survives st_coordinates() with
  # one row per *vertex*, misaligning coordinates with data rows downstream
  # (fit_bayesian_spatial_model()'s cbind, GWmodel's Spatial coercion).
  # coerce_to_points() maps MULTIPOINT to its centroid, so single-vertex
  # MULTIPOINTs (common in shapefiles) are unchanged in practice.
  if (!all(sf::st_geometry_type(data_sf, by_geometry = TRUE) == "POINT")) {
    data_sf <- coerce_to_points(data_sf, pointize)
  }

  # Drop any Z or M coordinate.  Every backend works in 2-D map distance, but
  # POINT Z is ordinary input (GPS and GeoPackage elevations, KML's altitude,
  # PointZ shapefiles) and the sf -> sp coercion GWmodel needs keeps the third
  # column.  GWmodel then refused the data ("Please input correct coordinates
  # of data points"), fitted on 3-D distances above 2500 rows, and in
  # predict() on POINT Z newdata reshaped XYZ with matrix(, ncol = 2) into the
  # wrong prediction locations, with no error.  make_folds(),
  # estimate_sac_range() and the fold-separation check already drop ZM.
  data_sf <- .drop_zm(data_sf)

  if (!is.null(boundary)) {
    bnd <- if (inherits(boundary, "sfc")) sf::st_as_sf(boundary) else boundary
    if (!inherits(bnd, "sf"))
      stop("prep_model_data(): 'boundary' must be sf/sfc when supplied.",
           call. = FALSE)
    bnd <- ensure_projected(bnd)
    data_sf <- ensure_projected(data_sf, .crs_or_null(bnd))
  } else {
    data_sf <- ensure_projected(data_sf)
  }

  # For cleaning, only check columns that are actually present and relevant.
  # When require_response = FALSE (prediction mode), exclude the response

  # variable entirely so that NAs in the response column of reused datasets
  # do not silently drop valid prediction rows.
  check_vars <- if (require_response) c(response_var, predictor_vars) else predictor_vars
  clean_cols <- intersect(check_vars, names(data_sf))
  df <- sf::st_drop_geometry(data_sf)[, clean_cols, drop = FALSE]
  ok_cc <- stats::complete.cases(df)
  num_mask <- vapply(df, is.numeric, logical(1))
  ok_fin <- if (any(num_mask)) {
    apply(as.matrix(df[, num_mask, drop = FALSE]), 1L,
          function(r) all(is.finite(r)))
  } else rep(TRUE, nrow(df))

  # The attribute columns are not the whole story: st_geometry_type() reports
  # an EMPTY POINT as "POINT", so the POINT checks above (and the one in
  # fit_bayesian_spatial_model()) pass it straight through, and nothing here
  # ever looked at the coordinates themselves.  Such a row reaches GWmodel as
  # a raw sp coercion error naming no row, ranger(include_coords = TRUE) as
  # "Missing data in columns: ..x, ..y", and -- worst -- brms at predict time,
  # where one infinite coordinate makes max(Xgp) infinite, so the GP boundary
  # is infinite and EVERY basis function on the surface is NaN.  make_folds()
  # already drops these rows and says how many; do the same here.
  xy_mat  <- sf::st_coordinates(sf::st_geometry(data_sf))
  ok_geom <- rep(TRUE, nrow(data_sf))
  if (nrow(xy_mat) == nrow(data_sf) && ncol(xy_mat) >= 2L) {
    ok_geom <- is.finite(xy_mat[, 1L]) & is.finite(xy_mat[, 2L])
  } else if (nrow(data_sf) > 0L) {
    # Empty geometries drop out of st_coordinates() entirely, so the row
    # counts disagree; fall back to the per-feature test.
    ok_geom <- !sf::st_is_empty(sf::st_geometry(data_sf))
    xy_ok <- sf::st_coordinates(sf::st_geometry(data_sf)[ok_geom])
    if (nrow(xy_ok) == sum(ok_geom) && ncol(xy_ok) >= 2L)
      ok_geom[ok_geom] <- is.finite(xy_ok[, 1L]) & is.finite(xy_ok[, 2L])
  }

  keep <- ok_cc & ok_fin & ok_geom
  dropped     <- sum(!keep)
  dropped_geo <- sum(!ok_geom)
  drop_idx    <- unname(which(!keep))
  drop_reason <- as.character(ifelse(!ok_geom[drop_idx], "geometry",
                              ifelse(!ok_cc[drop_idx], "missing", "non_finite")))

  # How many of the dropped rows each column is responsible for.  A row can
  # count under several columns, so these do not sum to `dropped`.
  col_bad <- lapply(df, function(v) {
    b <- is.na(v)
    if (is.numeric(v)) b <- b | !is.finite(v)
    b
  })
  by_column <- vapply(col_bad, function(b) sum(b & !keep), integer(1))
  if (dropped_geo > 0L) by_column <- c(by_column, `(geometry)` = dropped_geo)

  # Where the dropped rows sit on the response.  Rows that go missing on a
  # predictor for reasons tied to the outcome leave a fit that describes the
  # complete rows only: in the case that prompted this, 41% of the rows were
  # dropped and the mean response moved from 0.03 to -0.62 with nothing said.
  # Measured in SDs of the response over the rows where it is finite, so the
  # number does not depend on its units.
  response_std_diff <- NA_real_
  shift_se   <- NA_real_
  resp_means <- c(dropped = NA_real_, kept = NA_real_)
  if (require_response && dropped > 0L && sum(keep) > 0L) {
    y_all <- sf::st_drop_geometry(data_sf)[[response_var]]
    if (is.logical(y_all)) y_all <- as.numeric(y_all)
    if (is.numeric(y_all)) {
      y_drop <- y_all[!keep]; y_drop <- y_drop[is.finite(y_drop)]
      y_keep <- y_all[keep]
      s_all  <- stats::sd(c(y_drop, y_keep))
      if (length(y_drop) > 0L && is.finite(s_all) && s_all > 0) {
        resp_means <- c(dropped = mean(y_drop), kept = mean(y_keep))
        response_std_diff <- unname((resp_means[["dropped"]] -
                                       resp_means[["kept"]]) / s_all)
        # Two standard errors of a difference of means in SD units: below
        # that, a handful of dropped rows says nothing about why they are
        # missing, and the message states the difference without reading it.
        shift_se <- sqrt(1 / length(y_drop) + 1 / length(y_keep))
      }
    }
  }

  if (dropped > 0L) {
    n_miss <- sum(drop_reason == "missing")
    n_nonf <- sum(drop_reason == "non_finite")
    cols_txt <- paste(sprintf("%s %d", names(by_column), by_column)[by_column > 0L],
                      collapse = ", ")
    if (require_response) {
      shift_txt <- if (is.finite(response_std_diff))
        sprintf(paste0(" The response '%s' averages %s in the dropped rows ",
                       "that have one and %s in the kept rows, a difference ",
                       "of %.2f SD%s."),
                response_var, format(signif(resp_means[["dropped"]], 4)),
                format(signif(resp_means[["kept"]], 4)), response_std_diff,
                if (abs(response_std_diff) >= 0.2 &&
                    abs(response_std_diff) > 2 * shift_se)
                  ", more than two standard errors, so the rows are not missing at random and the fit describes the kept rows, not the layer"
                else "")
      else ""
      .warn_classed(
        "spatialkit_rows_dropped",
        paste0("prep_model_data(): dropped %d of %d row(s) (%.1f%%): %d with ",
               "a missing value, %d with a non-finite value, %d with an empty ",
               "or non-finite geometry. Rows affected by column: %s.%s The ",
               "fit uses the %d row(s) that remain; attr(, \"dropped\") lists ",
               "the rows removed."),
        dropped, nrow(data_sf), 100 * dropped / nrow(data_sf),
        n_miss, n_nonf, dropped_geo, cols_txt, shift_txt, sum(keep))
    } else {
      # Prediction mode.  Every predict() method returns NA at exactly these
      # rows, so the result itself says which they are; an R warning per call
      # would fire on every chunk of a predict_surface() grid with covariate
      # gaps.  The log line and the record stay.
      .log_warn(paste0("prep_model_data(): dropping %d row(s) with non-finite or ",
                       "missing values (%d of them with an empty or non-finite ",
                       "geometry). Rows affected by column: %s."),
                dropped, dropped_geo, cols_txt)
    }
  }

  out <- data_sf[keep, , drop = FALSE]

  # Two marks on the rows that remain.  A predictor that takes one value in
  # every row cannot be told from the intercept (it is one source of the
  # rank-deficient GWR windows), and repeated locations change what a
  # neighbour count and a distance quantile mean.  Both passed through
  # unmarked.
  out_df <- sf::st_drop_geometry(out)
  constant_predictors <- character(0)
  if (nrow(out) >= 2L && length(predictor_vars) > 0L) {
    is_const <- vapply(intersect(predictor_vars, names(out_df)), function(cn) {
      v <- out_df[[cn]]
      if (is.numeric(v)) isTRUE(max(v) == min(v)) else length(unique(v)) <= 1L
    }, logical(1))
    constant_predictors <- names(is_const)[is_const]
  }
  xy_out <- xy_mat
  n_duplicate_locations <- NA_integer_
  if (nrow(xy_out) == length(keep) && ncol(xy_out) >= 2L) {
    xy_out <- xy_out[keep, , drop = FALSE]
    # complex(): one hashable value per location, far faster than
    # duplicated() on a two-column matrix, which pastes each row.
    n_duplicate_locations <- sum(duplicated(complex(real = xy_out[, 1L],
                                                    imaginary = xy_out[, 2L])))
  } else if (nrow(out) > 0L) {
    xy_k <- sf::st_coordinates(sf::st_geometry(out))
    if (nrow(xy_k) == nrow(out) && ncol(xy_k) >= 2L)
      n_duplicate_locations <- sum(duplicated(complex(real = xy_k[, 1L],
                                                      imaginary = xy_k[, 2L])))
  }
  if (require_response) {
    if (length(constant_predictors) > 0L)
      .warn_classed(
        "spatialkit_constant_predictor",
        paste0("prep_model_data(): predictor(s) %s take a single value in all ",
               "%d row(s) kept, so their coefficients cannot be told from the ",
               "intercept (a GWR window is rank-deficient, a Bayesian slope is ",
               "its prior). They are kept in the data; remove them from ",
               "`predictor_vars`."),
        paste(sQuote(constant_predictors), collapse = ", "), nrow(out))
    if (isTRUE(n_duplicate_locations > 0L))
      .log_info(paste0("prep_model_data(): %d of %d row(s) kept share their ",
                       "coordinates with an earlier row (repeated locations); ",
                       "attr(, \"dropped\")$n_duplicate_locations records the ",
                       "count."),
                n_duplicate_locations, nrow(out))
  }
  # What was dropped, with the rows' identities.  The three masks above used
  # to be collapsed into the counts of that log line and nothing else, so a
  # fit's $n was the post-cleaning row count with no trace of how many rows
  # were lost or why -- and silently losing a third of the rows is a classic
  # cause of a suspiciously good score.  Row positions are those of `data_sf`
  # as passed (point coercion and projection keep the row order); the row IDs
  # are given as well when the layer carries `..row_id`, which is what the
  # cross-validation functions name rows by.
  # Stamped and classed, so that the record cannot outlive the rows it
  # describes: `out[i, ]` returns a plain layer without it.  The two marks on
  # the kept rows ride in the same record for that reason: as plain
  # attributes they would survive `[` and describe rows the subset no longer
  # has.
  out <- .set_row_record(out, "dropped", list(
    n          = as.integer(dropped),
    n_geometry = as.integer(dropped_geo),
    which      = drop_idx,
    row_id     = if ("..row_id" %in% names(data_sf))
      data_sf[["..row_id"]][drop_idx] else NULL,
    reason     = drop_reason,
    by_column  = by_column,
    response_std_diff     = response_std_diff,
    constant_predictors   = constant_predictors,
    n_duplicate_locations = as.integer(n_duplicate_locations)
  ))
  out
}


#' Drop the Z and M coordinates of a POINT layer, when it has any
#'
#' Tested on the coordinates themselves.  A POINT is a numeric vector of two
#' values (XY), three (XYZ or XYM) or four (XYZM), so a third coordinate on
#' any row shows as more than two values per feature.  sf's
#' \code{z_range}/\code{m_range} attributes are not enough: sf does not set
#' them on points built with \code{st_as_sf(coords = c("x", "y", "z"))}.  And
#' \code{sf::st_zm()} rebuilds every geometry with an R call each, too slow to
#' spend on every 2-D layer a \code{predict()} call prepares.
#'
#' @param x An sf object with POINT geometry.
#' @return \code{x}, with XY geometry.
#' @keywords internal
#' @noRd
.drop_zm <- function(x) {
  g <- sf::st_geometry(x)
  if (length(unlist(g, use.names = FALSE)) > 2L * length(g))
    x <- sf::st_zm(x, drop = TRUE, what = "ZM")
  x
}


#' Heuristic length-scale bounds for a squared-exponential GP
#'
#' Computes sensible prior bounds for the GP length-scale parameter \eqn{\ell}
#' of a squared-exponential (exponentiated-quadratic) kernel,
#' \eqn{k(h) = \exp(-h^2 / (2\ell^2))}{k(h) = exp(-h^2 / (2 l^2))}.
#' The "effective range" where correlation drops to ~5% is
#' \eqn{\ell \sqrt{2 \ln 20} \approx 2.45\,\ell}{l * sqrt(2 log(20)) ≈ 2.45 l}.
#'
#' Subsamples large datasets to avoid O(n^2) memory and time cost.
#'
#' These are the bounds a length-scale \emph{prior} is calibrated over, not
#' the scales a fitted model can resolve: that depends on the basis size
#' (\code{gp_k} in \code{\link{fit_bayesian_spatial_model}()}, which reports
#' it as \code{$info$gp_ell_min}, and in map units as
#' \code{$info$gp_resolution}).  Both bounds are fixed fractions of the
#' spread of pairwise distances, so they do not shrink as points are added to
#' the same area: more points help the data identify a short range, but they
#' make neither these bounds nor the derived basis finer.
#'
#' The lower bound is the 25th percentile of the pairwise distances, as an
#' effective range.  On real layouts that is 14 to 33 percent of the extent
#' along each axis (about 32 km on a uniform 100 km square), and the prior
#' built on it puts 1 percent of its mass below it.  A surface whose range
#' is shorter needs two things, not one: a basis fine enough to represent it
#' and a prior that allows it.  A larger \code{gp_k} under the prior
#' calibrated here changes almost nothing, because the prior still holds the
#' length-scale above this bound.  \code{fit_bayesian_spatial_model()}
#' therefore moves the prior's lower tail down to the length-scale the basis
#' resolves whenever \code{gp_resolution} or an explicit \code{gp_k} asks
#' for a finer basis than the one derived from these bounds.
#'
#' @param coords_xy Numeric matrix or data.frame of coordinates with at least
#'   two columns; the first two are used, and replicated rows are collapsed
#'   before the distance quantiles are taken.  \code{brms::gp()} defaults to
#'   \code{gr = TRUE} and reduces its covariates to unique rows, so a
#'   heavily-sampled station would otherwise weight the quantile by how often it
#'   was measured rather than by where the sites are.
#' @param q_small Numeric quantile for the lower bound, a single number in
#'   \code{[0, 1]}. Default 0.25.
#'   Previous versions used 0.1, but the 10th percentile can be
#'   dominated by within-cluster spacing in clustered data, producing
#'   a misleadingly small lower bound.
#' @param max_n Maximum number of \emph{distinct} points to use in the distance
#'   computation.  Default 1000. Set to \code{Inf} to use all of them.
#' @return Named numeric vector \code{c(lower, upper)} on the length-scale;
#'   \code{c(lower = 0.001, upper = 1)} when fewer than two distinct locations
#'   or no positive distances remain.
#' @family model fitting
#' @examples
#' set.seed(1)
#' xy <- cbind(runif(50), runif(50))
#' gp_lengthscale_bounds(xy)
#' @export
gp_lengthscale_bounds <- function(coords_xy, q_small = 0.25, max_n = 1000L) {
  # Without these, a vector `coords_xy` fails inside .safe_dist() with
  # "argument is of length zero" (nrow(NULL)) and an out-of-range `q_small`
  # fails inside quantile() -- neither message naming the argument at fault.
  if (!(is.matrix(coords_xy) || is.data.frame(coords_xy)))
    stop("gp_lengthscale_bounds(): 'coords_xy' must be a matrix or data.frame ",
         "of coordinates, not a ", class(coords_xy)[1L], ".", call. = FALSE)
  if (ncol(coords_xy) < 2L)
    stop("gp_lengthscale_bounds(): 'coords_xy' must have at least two columns ",
         "(x and y); it has ", ncol(coords_xy), ".", call. = FALSE)
  if (!is.numeric(q_small) || length(q_small) != 1L || !is.finite(q_small) ||
      q_small < 0 || q_small > 1)
    stop("gp_lengthscale_bounds(): 'q_small' must be a single number in [0, 1].",
         call. = FALSE)

  # "the first two are used" is a documented promise, and .safe_dist() ->
  # stats::dist() uses EVERY column.  sf::st_coordinates() on POINT Z geometry
  # returns X, Y, Z, so the documented entry point silently returned 3-D
  # bounds in mixed units (elevation metres mixed into map metres).
  coords_xy <- as.matrix(coords_xy)[, 1:2, drop = FALSE]

  # brms::gp() defaults to gr = TRUE and reduces its covariates to their
  # UNIQUE rows before doing anything else, so replicated locations (repeat
  # visits, a heavily-sampled reference station) would otherwise weight the
  # distance quantile by how often a site was measured rather than by where
  # the sites are.
  coords_xy <- coords_xy[!duplicated(coords_xy), , drop = FALSE]
  if (nrow(coords_xy) < 2L) return(c(lower = 0.001, upper = 1))

  d <- .safe_dist(coords_xy, max_n = max_n)
  d <- d[d > 0]
  if (!length(d)) return(c(lower = 0.001, upper = 1))
  dmax  <- max(d)
  dq    <- stats::quantile(d, probs = q_small, names = FALSE, type = 7)
  # For SE kernel: 5% correlation at distance d  =>  l = d / sqrt(2 ln 20)
  se_range_factor <- sqrt(2 * log(20))          # approx. 2.448
  lower <- dq   / se_range_factor
  upper <- dmax / se_range_factor
  upper <- max(upper, lower * 1.2)               # ensure separation
  c(lower = as.numeric(lower), upper = as.numeric(upper))
}


# The cap on the TOTAL number of basis functions, gp_k^2.  One constant, read
# by .gp_basis_spec() for the derived k and by fit_bayesian_spatial_model()
# for a requested gp_resolution, so the two cannot disagree.
.gp_max_basis <- 2500L

# sqrt(2 log 20): the distance, in length-scales, at which a
# squared-exponential correlation has fallen to 5 percent.  An "effective
# range" is this many length-scales.
.gp_se_range_factor <- sqrt(2 * log(20))


#' Choose GP basis count and boundary factor from the length-scale/domain ratio
#'
#' Implements the practical recommendations of Riutort-Mayol et al. (2023,
#' Statistics and Computing 33:17) for the squared-exponential kernel used by
#' \code{brms::gp()}:
#' \preformatted{
#'   c >= 3.2 * (ell/S),  c >= 1.25
#'   m >= 1.75 * c / (ell/S)
#' }
#' where \code{S} is \strong{the same domain measure \code{brms::gp(c = )}
#' multiplies}, namely the full pooled range of the column-centred coordinates.
#'
#' \strong{Why the range and not the half-range.}  Riutort-Mayol et al. state
#' their inequalities against the domain half-range, but both are really
#' constraints on the \emph{boundary} \eqn{L}: contain the longest plausible
#' range (\eqn{L \ge 3.2\,\ell_{upper}}) and resolve the shortest
#' (\eqn{m \ge 1.75\,L/\ell_{lower}}).  brms builds that boundary as
#' \preformatted{
#'   choose_L <- function(x, c) c * max(1, max(x) - min(x))
#' }
#' over the column-centred covariate matrix, pooled across dimensions
#' (\code{brms:::.data_gp()}).  That is the FULL range, about twice the per-axis
#' half-range.  Deriving \code{c} against the half-range and handing the result
#' to \code{brms::gp()} therefore built a boundary twice as wide as intended:
#' \code{k} was sized for a boundary half the real one, so the GP was
#' systematically under-resolved, and the smallest resolvable length-scale
#' (\code{gp_ell_min}) was understated by the same factor, which made the
#' post-fit adequacy diagnostic, whose whole job is to catch under-resolution,
#' twice too lenient to fire.  Working in brms's own units removes both.
#'
#' The floor is \code{1.25} where Riutort-Mayol et al. give \code{1.2}, because
#' the floor is convention-dependent too: it is brms's own default
#' (\code{c = 5/4}), and on the range convention it is the more generous of the
#' two.
#'
#' \code{m} is the count PER DIMENSION.  brms expands a full tensor grid over
#' the GP covariates, so the fitted model carries \code{m^D} basis functions
#' (\code{D = 2} for \code{gp(..x, ..y)}).  See brms/R/data-predictor.R,
#' \code{data_gp()}:  \code{out[[paste0("NBgp", pi)]] <- k ^ D}.  Because
#' \code{brms::gp()} accepts a single \code{k} shared across dimensions, the
#' per-dimension recommendations are collapsed into one value.
#'
#' The boundary factor is set from the UPPER length-scale bound, because it must
#' contain the longest plausible correlation range, while the basis count is set
#' from the LOWER bound, because it must resolve the shortest.  Deriving both
#' from the same length-scale estimate is the point: the previous rule scaled
#' the basis count with the number of observations, which made the approximation
#' cost grow as n while adding no resolution the data supported.
#'
#' @param coords_xy Numeric matrix of (already scaled) coordinates.  Only the
#'   first two columns are used, and replicated rows are collapsed, mirroring
#'   \code{brms::gp(gr = TRUE)}.
#' @param ls_bounds Named numeric \code{c(lower, upper)}, as returned by
#'   \code{gp_lengthscale_bounds()}.
#' @param k_min Integer floor on the per-dimension basis count.
#' @param max_basis Integer cap on the TOTAL basis count (\code{k^2}).  The
#'   per-dimension ceiling is derived from this as \code{floor(sqrt(max_basis))},
#'   so there is a single cap, with no second one to contradict it.
#' @param c Optional boundary factor, already validated, that \code{k} must be
#'   sized for; \code{NULL} (default) derives it here.  \code{k} grows with
#'   \code{c}, so a caller that fixes the boundary must size the basis for
#'   \emph{that} boundary: the \code{k} derived for the default \code{c} cannot
#'   resolve the lower bound inside a wider one.
#' @return A list with \code{k} (integer, per dimension), \code{c} (numeric;
#'   the \code{c} argument when one was given),
#'   \code{S} (numeric; the pooled full range of the column-centred coordinates
#'   AFTER collapsing replicated rows, i.e. exactly what \code{brms::gp(c = )}
#'   multiplies under its default \code{gr = TRUE}), \code{capped}
#'   (logical) and \code{cmeans} (the column means the coordinates were
#'   centred on, which brms stores in the fit's GP basis and centres every
#'   later \code{newdata} on).
#' @keywords internal
#' @noRd
.gp_basis_spec <- function(coords_xy, ls_bounds,
                           k_min = 10L, max_basis = .gp_max_basis, c = NULL) {
  # Reproduce brms::choose_L()'s domain measure exactly: centre each column,
  # then take the range over the POOLED matrix.  na.rm mirrors brms.
  xy <- as.matrix(coords_xy)[, 1:2, drop = FALSE]
  # brms::gp() defaults to gr = TRUE, and brms:::.data_gp() reduces Xgp to its
  # UNIQUE rows BEFORE taking cmeans and the range.  Centring and ranging over
  # all rows therefore disagreed with the boundary brms actually builds
  # whenever two observations share a location: with one heavily-sampled
  # station the boundary came out at 0.68 of the intended one, and everything
  # derived from S -- gp_c, and the gp_ell_min threshold the post-fit adequacy
  # warning compares posterior length-scale draws against -- was wrong by the
  # same factor.  Mirror the reduction.
  xy <- xy[!duplicated(xy), , drop = FALSE]
  if (!nrow(xy)) xy <- matrix(0, nrow = 1L, ncol = 2L)
  cmeans <- colMeans(xy, na.rm = TRUE)
  Xc <- sweep(xy, 2L, cmeans)
  S  <- suppressWarnings(
    max(1, max(Xc, na.rm = TRUE) - min(Xc, na.rm = TRUE)))
  if (!is.finite(S) || S <= 0) S <- 1

  # Ratios of length-scale to the domain measure brms will use.
  r_lo <- max(ls_bounds[["lower"]] / S, .Machine$double.eps)  # must resolve
  r_hi <- max(ls_bounds[["upper"]] / S, r_lo)                 # must contain

  # 1.25, not 1.2: the floor is stated on the half-range convention by
  # Riutort-Mayol et al., and 5/4 is brms's own default on this one.
  # A caller's own c replaces the derived one BEFORE k is sized from it: k was
  # once always sized for the derived c, so gp_c = 3 on a layer whose derived
  # c was 1.63 kept k = 23 where the rule gives 43, and the basis could not
  # resolve the lower length-scale bound it was supposed to.
  c_val <- if (is.null(c)) max(3.2 * r_hi, 1.25) else as.numeric(c)
  k_raw <- ceiling(1.75 * c_val / r_lo)

  k_max  <- as.integer(floor(sqrt(max_basis)))
  k_val  <- min(max(k_raw, k_min), k_max)
  capped <- k_raw > k_max

  list(k = as.integer(k_val), c = as.numeric(c_val),
       S = as.numeric(S), capped = capped, cmeans = as.numeric(cmeans))
}


# NOTE: the previously defined .matern_rho() and .safe_ll() helpers were
# removed as dead code (no callers anywhere in the package); mvtnorm was
# dropped from Suggests along with .safe_ll().
