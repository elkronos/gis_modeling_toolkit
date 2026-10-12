# tests/testthat/test-fix-resolution-moran-diagnostic.R
# ---------------------------------------------------------------------------
# moran_z is a diagnostic of the cell means, not a count rule: Moran's z of
# the cell means is near 0 at 10 to 20 cells whatever the data and grows as
# cells get finer than the range.  select_resolution() reads it only when
# asked, summary() leaves it out by default, and a profile handed to the
# tessellation functions is read by cp or not at all.
# ---------------------------------------------------------------------------

fmd_field <- function(n = 500, seed = 2, a = 100, nugget = 0.3) {
  set.seed(seed)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  D  <- as.matrix(stats::dist(xy))
  xy$w <- rnorm(n)
  xy$z <- as.numeric(t(chol(exp(-D / a) + diag(1e-8, n))) %*% rnorm(n)) +
    rnorm(n, sd = sqrt(nugget))
  sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
}
fmd_sac <- function(pts, range = 240, nugget = 0.3, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, range / 3), stringsAsFactors = FALSE)
  structure(range, class = c("sac_range", "numeric"), variogram_model = vm,
            crs = sf::st_crs(pts))
}
fmd_quiet <- function(expr) suppressWarnings(suppressMessages(expr))
fmd_bnd <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 32632))


test_that("select_resolution() refuses moran_z unless asked, and says what it measures", {
  pts <- fmd_field()
  p <- fmd_quiet(resolution_profile(pts, "z", sac = fmd_sac(pts), n_levels = 10, nstart = 3))
  ok <- is.finite(p$moran_z)
  expect_true(any(ok))

  e <- tryCatch(select_resolution(p, "moran_z"), error = function(e) e)
  expect_s3_class(e, "spatialkit_diagnostic_criterion")
  msg <- conditionMessage(e)
  expect_match(msg, "a diagnostic, not a count rule")
  expect_match(msg, "`diagnostic = TRUE`", fixed = TRUE)
  expect_match(msg, "Moran's z of the cell means over each cell's 8 nearest cells")
  expect_match(msg, sprintf("computable from %d", min(p$levels[ok])))
  expect_match(msg, "does not say how well the cells represent the response")

  # Asked for, it is read as before, and the result carries what it means.
  z <- select_resolution(p, "moran_z", diagnostic = TRUE)
  expect_identical(z$best, p$levels[ok][which.min(abs(p$moran_z[ok]))])
  expect_match(z$note, "diagnostic, not a count rule")
  expect_output(print(z), "least autocorrelated across")
  expect_error(select_resolution(p, "moran_z", diagnostic = NA), "TRUE or FALSE")

  # With predictors the statistic is of the cell-mean regression residuals,
  # and it is refused all the same.
  pw <- fmd_quiet(resolution_profile(pts, "z", predictor_vars = "w", sac = fmd_sac(pts),
                                     n_levels = 8, nstart = 3))
  ew <- tryCatch(select_resolution(pw, "moran_z"), error = function(e) e)
  expect_s3_class(ew, "spatialkit_diagnostic_criterion")
  expect_match(conditionMessage(ew), "residuals of the cell means regressed on the cell-mean predictors")

  # A profile where it is NA everywhere still says that first.
  geo <- resolution_profile(pts, n_levels = 4, nstart = 1)
  expect_error(select_resolution(geo, "moran_z"), "NA at every level")
})

test_that("summary() leaves moran_z out unless it is named, and out of the comparison", {
  pts <- fmd_field()
  p <- fmd_quiet(resolution_profile(pts, "z", sac = fmd_sac(pts), n_levels = 10, nstart = 3))
  s <- summary(p)
  expect_false("moran_z" %in% s$criterion)
  expect_true(all(c("cp", "reliability") %in% s$criterion))

  named <- summary(p, criteria = c("cp", "reliability", "moran_z"))
  expect_identical(named$criterion, c("cp", "reliability", "moran_z"))
  expect_identical(named$best[3L], select_resolution(p, "moran_z", diagnostic = TRUE)$best)
  # The levels the count criteria agree on are theirs alone.
  expect_identical(attr(named, "common"), attr(summary(p, criteria = c("cp", "reliability")),
                                               "common"))
  out <- utils::capture.output(print(named))
  expect_true(any(grepl("moran_z is a diagnostic of the cell means, not a count rule", out)))
  rng <- range(named$best[1:2])
  expect_true(any(grepl(sprintf("picks span %d to %d cells", rng[1L], rng[2L]), out)))
  expect_lte(max(nchar(out)), 80L)
})

test_that("a profile whose cp is NA is not read by another criterion in silence", {
  resolve <- spatialkit:::.resolve_cell_count
  pts <- fmd_field()
  sac <- fmd_sac(pts)

  # Geometry-only, given a variogram: cp is NA, reliability is finite.  It
  # used to be read by reliability and handed on as that many cells.
  pg <- fmd_quiet(resolution_profile(pts, sac = sac, n_levels = 8, nstart = 2))
  expect_false(any(is.finite(pg$cp)))
  expect_true(any(is.finite(pg$reliability)))
  expect_error(select_resolution(pg), "NA at every level")
  e <- tryCatch(resolve(pg, "n", "get_voronoi_seeds"), error = function(e) e)
  expect_s3_class(e, "spatialkit_no_default_criterion")
  expect_match(conditionMessage(e), "`cp`, the criterion a profile is read by, is NA at every level")
  expect_match(conditionMessage(e), "Finite on it: reliability")
  expect_match(conditionMessage(e), "select_resolution(<profile>, criterion = \"reliability\")",
               fixed = TRUE)
  expect_error(get_voronoi_seeds(fmd_bnd, method = "kmeans", n = pg, sample_points = pts),
               class = "spatialkit_no_default_criterion")
  expect_error(build_tessellation(pts, boundary = fmd_bnd, method = "hex",
                                  approx_n_cells = pg, quiet = TRUE),
               class = "spatialkit_no_default_criterion")
  # The way the message names works.
  sel <- select_resolution(pg, "reliability")
  expect_equal(fmd_quiet(resolve(sel, "n", "f"))$n, sel$best)

  # A response profile whose cp is NA at every level (a range alone, and
  # levels too fine for the fallback noise variance) used to be read by
  # moran_z.
  pf <- fmd_quiet(resolution_profile(pts, "z", sac = 240, levels = c(260, 300, 350),
                                     nstart = 1))
  expect_false(any(is.finite(pf$cp)))
  expect_true(any(is.finite(pf$moran_z)))
  ef <- tryCatch(resolve(pf, "n", "get_voronoi_seeds"), error = function(e) e)
  expect_s3_class(ef, "error")
  expect_match(conditionMessage(ef), "no criterion finite at any level")
  expect_match(conditionMessage(ef), "`moran_z` is finite on it, but it is a diagnostic")
})

test_that("a geometry-only profile with an elbow is still read by the elbow", {
  resolve <- spatialkit:::.resolve_cell_count
  set.seed(6)
  ctr <- cbind(runif(6, 100, 900), runif(6, 100, 900))
  g <- rep(1:6, each = 50)
  cl <- sf::st_as_sf(data.frame(x = ctr[g, 1] + rnorm(300, 0, 10),
                                y = ctr[g, 2] + rnorm(300, 0, 10)),
                     coords = c("x", "y"), crs = 32632)
  geo <- resolution_profile(cl, n_levels = 8, nstart = 5)
  expect_true(any(is.finite(geo$elbow)))
  r <- fmd_quiet(resolve(geo, "n", "get_voronoi_seeds"))
  expect_equal(r$n, select_resolution(geo, "elbow")$best)
  expect_identical(r$criterion, "elbow")
  expect_match(r$from, "select_resolution\\(criterion = \"elbow\"\\)")
})

test_that("determine_optimal_levels() says, with a class, when Moran's z ranked its counts", {
  fmd_dol <- function(...) {
    ws <- list()
    value <- withCallingHandlers(determine_optimal_levels(...), warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    hit <- vapply(ws, inherits, logical(1), "spatialkit_diagnostic_ranking")
    list(value = value, n = sum(hit), w = if (any(hit)) ws[[which(hit)[1L]]] else NULL)
  }
  pts <- fmd_field(seed = 701)

  # A response alone: the first-ranked count is the first the statistic can
  # be computed at, on this field as on white noise.  It came back with the
  # no-elbow warning and nothing about what had ranked it.
  r <- fmd_dol(pts, max_levels = 40, response_var = "z", criterion = "morans_i")
  d <- attr(r$value, "diagnostics")
  expect_identical(r$n, 1L)
  expect_identical(r$w$k, r$value[1L])
  expect_identical(r$w$criterion, "morans_i")
  expect_false(r$w$upgraded)
  expect_identical(r$w$first_computable, min(d$eval_ks[is.finite(d$moran_z[d$eval_ks])]))
  expect_identical(r$w$window, range(d$eval_ks))
  msg <- conditionMessage(r$w)
  expect_match(msg, "criterion = 'morans_i' ranks the candidates by |Moran's z| of the cell means",
               fixed = TRUE)
  expect_match(msg, "a diagnostic of whether neighbouring cells are independent, not a count rule")
  expect_match(msg, sprintf("k = %d is ranked first; the statistic is computable from k = %d",
                            r$value[1L], r$w$first_computable))
  expect_match(msg, "Use criterion = \"geometric\" for the elbow's ranking", fixed = TRUE)
  expect_identical(r$w$k, r$w$first_computable)
  # It can be caught by class (the no-elbow warning, raised before it on
  # these uniform points, is let through).
  expect_identical(
    withCallingHandlers(
      tryCatch(determine_optimal_levels(pts, max_levels = 40, response_var = "z",
                                        criterion = "morans_i"),
               spatialkit_diagnostic_ranking = function(w) "caught"),
      spatialkit_no_elbow = function(w) invokeRestart("muffleWarning")),
    "caught")

  # With both model variables "geometric" is read as "combined", so the way
  # back to the elbow's ranking has two parts when the criterion was named:
  # the criterion has to go and a variable has to go.  (The advice named the
  # variable alone, and with criterion = "morans_i" still in the call that
  # ranks by Moran's z all the same.)
  rb <- fmd_dol(pts, max_levels = 40, response_var = "z", predictor_vars = "w",
                criterion = "morans_i")
  expect_identical(rb$n, 1L)
  expect_false(rb$w$upgraded)
  expect_match(conditionMessage(rb$w),
               paste0("Use criterion = \"geometric\" and leave out `response_var` or ",
                      "`predictor_vars` (with both, \"geometric\" is read as \"combined\") ",
                      "for the elbow's ranking"),
               fixed = TRUE)
  # Following it gives the elbow's ranking, with no such warning.
  expect_identical(fmd_dol(pts, max_levels = 40, response_var = "z",
                           criterion = "geometric")$n, 0L)

  # "combined", asked for or taken because both model variables are given:
  # the same condition, and the message says which.
  rc <- fmd_dol(pts, max_levels = 40, response_var = "z", criterion = "combined")
  expect_identical(rc$n, 1L)
  expect_false(rc$w$upgraded)
  expect_match(conditionMessage(rc$w), "the WSS curve's sag, is flat (the curve has no elbow)",
               fixed = TRUE)
  ru <- fmd_dol(pts, max_levels = 40, response_var = "z", predictor_vars = "w")
  expect_identical(ru$n, 1L)
  expect_true(ru$w$upgraded)
  expect_identical(ru$w$criterion, "combined")
  mu <- conditionMessage(ru$w)
  expect_match(mu, "(taken because `response_var` and `predictor_vars` were both supplied)",
               fixed = TRUE)
  expect_match(mu, "the residuals of the cell means regressed on the cell-mean predictors")
  expect_match(mu, "Leave out `response_var` or `predictor_vars` for the elbow's ranking",
               fixed = TRUE)

  # No such warning where Moran's z ordered nothing: the geometric criterion,
  # a window below ten cells (the default max_levels), and an elbow below
  # ten cells, where "combined" returns the geometric ranking.
  expect_identical(fmd_dol(pts, max_levels = 40, response_var = "z")$n, 0L)
  expect_identical(fmd_dol(pts, response_var = "z", predictor_vars = "w")$n, 0L)
  set.seed(6)
  ctr <- cbind(runif(8, 0, 1000), runif(8, 0, 1000))
  g <- rep(1:8, each = 50)
  cl <- sf::st_as_sf(data.frame(x = ctr[g, 1] + rnorm(400, 0, 8), y = ctr[g, 2] + rnorm(400, 0, 8),
                                z = rnorm(400) + g, w = rnorm(400)),
                     coords = c("x", "y"), crs = 32632)
  r8 <- fmd_dol(cl, max_levels = 14, response_var = "z", predictor_vars = "w")
  expect_identical(r8$n, 0L)
  expect_identical(r8$value[1L], 8L)
})
