# tests/testthat/test-fix-cells-variogram-fallback.R
# ---------------------------------------------------------------------------
# summarize_by_cell(deff = "variogram") when the request cannot be met.
#
# It used to fall back to deff = 1: the uncorrected interval, under the same
# column names, with deff_applied FALSE.  That interval has the worst coverage
# of the grand mean (0.41 to 0.64 on real tract layers against 0.93 to 0.96
# for Kish), and the variogram's range goes unidentified on about half of
# such layers.  The fallback is "kish" now, with the classed warning kept and
# a `deff_method` column that says which correction each row carries.
# ---------------------------------------------------------------------------

# Six cells of 25 points on a correlated field with a cell effect, so the
# Kish ICC is clearly positive.
.fv_field <- function(seed = 77) {
  set.seed(seed)
  n <- 150
  x <- rep(seq(0, 500, length.out = 6), each = 25) + runif(n, 0, 60)
  y <- runif(n, 0, 60)
  z <- as.numeric(t(chol(exp(-as.matrix(stats::dist(cbind(x, y))) / 40) +
                           diag(1e-6, n))) %*% rnorm(n)) + rep(rnorm(6), each = 25)
  sf::st_as_sf(data.frame(x = x, y = y, z = z, w = z + rnorm(n),
                          poly_id = rep(1:6, each = 25)),
               coords = c("x", "y"), crs = 32632)
}

.fv_rejected <- function(reason = "fitted range exceeds the largest lag fitted") {
  structure(NA_real_, class = c("sac_range", "numeric"), rejected_reason = reason,
            variogram_model = data.frame(model = "Exp", psill = 1, range = 1e6),
            crs = sf::st_crs(32632))
}
.fv_sac <- function(model = data.frame(model = c("Nug", "Exp"), psill = c(0.2, 0.8),
                                       range = c(0, 40))) {
  structure(120, class = c("sac_range", "numeric"), variogram_model = model,
            crs = sf::st_crs(32632))
}

.fv_run <- function(expr) {
  w <- character(0); cls <- list(); conds <- list()
  val <- withCallingHandlers(
    suppressMessages(expr),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); cls[[length(cls) + 1L]] <<- class(c)
      conds[[length(conds) + 1L]] <<- c
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = w, classes = cls, conditions = conds)
}

test_that("a layer whose variogram range is refused returns Kish standard errors", {
  pts <- .fv_field()
  # The supplied fit is refused and so is the internal estimate.
  local_mocked_bindings(estimate_sac_range = function(...) .fv_rejected())
  res <- .fv_run(summarize_by_cell(pts, "z", deff = "variogram",
                                   sac = .fv_rejected(), conf_level = 0.95))
  out  <- res$value
  kish <- summarize_by_cell(pts, "z", deff = "kish", conf_level = 0.95)
  iid  <- summarize_by_cell(pts, "z", conf_level = 0.95)

  # Kish's numbers, not the uncorrected ones.
  for (col in c("..se_resp_z", "..neff_resp_z", "..df_resp_z", "..ci_lo_resp_z",
                "..ci_hi_resp_z", "cell_weight"))
    expect_equal(out[[col]], kish[[col]], tolerance = 1e-12, info = col)
  expect_true(all(out[["..se_resp_z"]] > 2 * iid[["..se_resp_z"]]))
  expect_true(all(out[["..ci_hi_resp_z"]] - out[["..ci_lo_resp_z"]] >
                    2 * (iid[["..ci_hi_resp_z"]] - iid[["..ci_lo_resp_z"]])))

  # One warning, classed, saying what was applied.
  expect_length(res$warnings, 1L)
  expect_true("spatialkit_deff_fallback" %in% res$classes[[1L]])
  expect_match(res$warnings, "requires a fitted variogram model")
  expect_match(res$warnings, "Falling back to deff = \"kish\" \\(ICC: response 'z' 0\\.[0-9]+\\)")
  expect_identical(res$conditions[[1L]]$requested, "variogram")
  expect_identical(res$conditions[[1L]]$used, "kish")

  # The rows are labelled, and the label survives binding results together.
  expect_identical(out$deff_method, rep("kish (variogram fallback)", 6))
  expect_identical(out$deff_applied, rep(TRUE, 6))
  bound <- dplyr::bind_rows(out, kish)
  expect_identical(bound$deff_method,
                   rep(c("kish (variogram fallback)", "kish"), each = 6))
  # The record names the method used, the one requested and why.
  da <- attr(out, "deff_applied")
  expect_identical(da$method, "kish")
  expect_identical(da$requested, "variogram")
  expect_match(da$fallback_reason, "none was available")
  expect_null(attr(kish, "deff_applied")$requested)
  # The ICC is recorded as for a "kish" request.
  expect_equal(attr(out, "icc")$resp, attr(kish, "icc")$resp, tolerance = 1e-12)
})

test_that("a model family the variogram route cannot read falls back to Kish too", {
  pts <- .fv_field()
  mat <- .fv_sac(data.frame(model = c("Nug", "Mat"), psill = c(0.2, 0.8),
                            range = c(0, 40)))
  res <- .fv_run(summarize_by_cell(pts, "z", "w", deff = "variogram", sac = mat))
  expect_length(res$warnings, 1L)
  expect_true("spatialkit_deff_fallback" %in% res$classes[[1L]])
  expect_match(res$warnings, "the supplied model uses 'Mat'")
  expect_match(res$warnings, "Falling back to deff = \"kish\"")
  kish <- summarize_by_cell(pts, "z", "w", deff = "kish")
  expect_equal(res$value[["..se_resp_z"]], kish[["..se_resp_z"]], tolerance = 1e-12)
  expect_equal(res$value[["..se_pred_w"]], kish[["..se_pred_w"]], tolerance = 1e-12)
  expect_identical(unique(res$value$deff_method), "kish (variogram fallback)")
})

test_that("when Kish's ICC is 0 as well, the uncorrected SE is returned and labelled", {
  # No cell effect and no spatial structure: the ICC estimate is not positive.
  set.seed(5)
  n <- 150
  flat <- sf::st_as_sf(data.frame(x = runif(n, 0, 600), y = runif(n, 0, 60),
                                  z = rep(c(-1, 1), length.out = n) + rnorm(n, 0, 0.01),
                                  poly_id = rep(1:6, each = 25)),
                       coords = c("x", "y"), crs = 32632)
  res <- .fv_run(summarize_by_cell(flat, predictor_vars = "z", deff = "variogram"))
  expect_length(res$warnings, 1L)
  expect_true("spatialkit_deff_fallback" %in% res$classes[[1L]])
  expect_match(res$warnings, "whose ICC is 0 \\(predictor 'z' 0\\.000\\), so the standard errors are the uncorrected ones")
  expect_identical(res$conditions[[1L]]$used, "none")
  iid <- summarize_by_cell(flat, predictor_vars = "z")
  expect_equal(res$value[["..se_pred_z"]], iid[["..se_pred_z"]], tolerance = 1e-12)
  expect_identical(res$value$deff_applied, rep(FALSE, 6))
  expect_identical(res$value$deff_method, rep("none (variogram fallback)", 6))
  expect_null(attr(res$value, "deff_applied"))
  expect_identical(unname(attr(res$value, "icc")$pred), 0)
})

test_that("deff_method names the correction on every path", {
  pts <- .fv_field()
  # No request: no label, the frame is as it always was.
  expect_false("deff_method" %in% names(summarize_by_cell(pts, "z")))
  expect_identical(unique(summarize_by_cell(pts, "z", deff = "variogram",
                                            sac = .fv_sac())$deff_method), "variogram")
  expect_identical(unique(summarize_by_cell(pts, "z", deff = "kish")$deff_method), "kish")
  expect_identical(unique(summarize_by_cell(pts, "z", deff = 2)$deff_method), "fixed")
  refused <- suppressWarnings(summarize_by_cell(pts, "z", deff = 0.5))
  expect_identical(unique(refused$deff_method), "none")
  expect_identical(unique(refused$deff_applied), FALSE)
  # Joined to cells, a cell no point fell in has no label, as it has no flag.
  sq <- function(x0) sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 60, 0), c(x0 + 60, 60),
                                               c(x0, 60), c(x0, 0))))
  cells <- sf::st_sf(poly_id = 1:7, geometry = sf::st_sfc(lapply(100 * 0:6, sq),
                                                          crs = 32632))
  joined <- summarize_by_cell(pts, "z", deff = "kish", cells_sf = cells)
  expect_identical(joined$deff_method, c(rep("kish", 6), NA))
})
