# tests/testthat/test-cp-nonconverged-nugget.R
# ---------------------------------------------------------------------------
# resolution_profile() gives Cp no nugget from a variogram model that did not
# converge, whatever reason its range was refused for.
#
# estimate_sac_range() gives a refused range one reason, and a falling
# variogram and a range past the largest lag come before the optimiser in it.
# A fit that stopped at gstat's iteration limit with a range past the largest
# lag is therefore "fitted range exceeds the largest lag fitted", and that is
# how most fits that do not converge end.  The profile looked for "converge"
# in the reason, so it took those fits' nuggets for Cp as fitted values while
# its documentation said a model that did not converge gives none.  Whether
# the model converged is now read from the model, where gstat's verdict rides
# (attr(variogram_model, "converged")); the reason decides only for a sac that
# carries no such flag.
# ---------------------------------------------------------------------------

# An exponential field with a known effective range on a 1000 m square (sill
# 1), simulated by Cholesky factorisation, with an optional east-west trend.
cpn_field <- function(n = 300, eff_range = 150, nugget = 0.1, slope = 0, seed = 1) {
  set.seed(seed)
  x <- runif(n, 0, 1000); y <- runif(n, 0, 1000)
  D <- as.matrix(stats::dist(cbind(x, y)))
  z <- drop(crossprod(chol(exp(-3 * D / eff_range) + diag(nugget, n)), rnorm(n)))
  sf::st_as_sf(data.frame(x = 5e5 + x, y = 5e6 + y, z = z + slope * x),
               coords = c("x", "y"), crs = 32632)
}

# A refused sac made by hand, with the attributes resolution_profile() reads
# as estimate_sac_range() sets them.  `converged` goes on the model, which is
# where the fit's verdict is kept; NULL leaves the model without one, as a
# hand-made sac or one built from REML parameters is.
cpn_refused <- function(reason, converged = NULL, nugget = 0.5, psill = 1) {
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(nugget, psill),
                   range = c(0, 3000), stringsAsFactors = FALSE)
  if (!is.null(converged)) attr(vm, "converged") <- converged
  structure(NA_real_, class = c("sac_range", "numeric"), variogram_model = vm,
            nugget = nugget, rejected_range = 9000, rejected_reason = reason,
            detrended = FALSE, crs = sf::st_crs(32632))
}

# Every R warning an expression raises, muffled.
cpn_warnings <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

# Cp took Mallows' noise variance from the finest level, and nothing from the
# variogram: no nugget, no reliability, no variogram record.
expect_cp_from_finest_level <- function(p) {
  cn <- attr(p, "cp_noise")
  expect_identical(cn$source, "finest-level residual mean square")
  expect_true(is.finite(cn$value) && cn$value > 0)
  expect_true(cn$level %in% p$levels)
  expect_true(all(is.finite(p$cp)))
  expect_true(all(is.na(p$reliability)))
  expect_null(attr(p, "variogram"))
}

# The model gstat fitted carries its verdict: one TRUE or FALSE.  Returns it.
# This is an expectation and not part of a skip condition.  A model that
# reaches the profile with no `converged` (rebuilt, subset or copied without
# its attributes somewhere on the way) sends every sac back to the reason
# text, which is the defect this file is for, and inside skip_if_not() that
# loss read as "gstat ended otherwise on this platform": the file passed with
# three tests skipped.  Only the value of the verdict, which does depend on
# the LAPACK build, is left to the skips below.  A sac with no model at all
# (both fits singular) has nothing to carry one, and gives cp no nugget
# either way.
expect_gstat_verdict <- function(sac) {
  vm <- attr(sac, "variogram_model")
  if (!is.data.frame(vm)) return(invisible(NA))
  v <- attr(vm, "converged", exact = TRUE)
  ok <- is.logical(v) && length(v) == 1L && !is.na(v)
  expect_true(ok, label = "the gstat model carrying `converged` as one TRUE or FALSE")
  invisible(if (ok) v else NA)
}

# gstat::fit.variogram() stubbed, as test-review3-S4-range-kriging.R stubs it:
# one exponential model with a nugget of 0.4 whose effective range (3 x 5000)
# is past every lag of a 1000 m square, returned with the warning gstat raises
# when its optimiser stops at the iteration limit, or without it.  Whether
# the real optimiser reaches that limit on a given draw depends on the LAPACK
# build.  What estimate_sac_range() and the profile make of the warning does
# not, and that is the part under test.
cpn_fit_stub <- function(converges) {
  function(object, model, ...) {
    m <- gstat::vgm(psill = 1, model = "Exp", range = 5000, nugget = 0.4)
    attr(m, "singular") <- FALSE
    attr(m, "SSErr") <- 1
    if (!converges)
      warning("No convergence after 200 iterations: try different initial values?")
    m
  }
}

PAST_LAGS  <- "fitted range exceeds the largest lag fitted"
DECREASING <- "empirical variogram decreases with distance"
NOT_CONV   <- "variogram model did not converge"


test_that(".vgm_converged() reads the verdict off the model, and NA when there is none", {
  vc <- spatialkit:::.vgm_converged
  vm <- data.frame(model = c("Nug", "Exp"), psill = c(0.5, 1), range = c(0, 100))
  expect_identical(vc(vm), NA)
  expect_identical(vc(NULL), NA)
  attr(vm, "converged") <- FALSE
  expect_false(vc(vm))
  attr(vm, "converged") <- TRUE
  expect_true(vc(vm))
  # Anything that is not one TRUE or FALSE says nothing.
  attr(vm, "converged") <- NA
  expect_identical(vc(vm), NA)
  attr(vm, "converged") <- "no"
  expect_identical(vc(vm), NA)
})


test_that("estimate_sac_range() hands gstat's verdict out on the model it returns", {
  skip_if_not_installed("gstat")
  # The link the rest of this file rests on, pinned where it is made.  The
  # profile reads attr(variogram_model, "converged") and no other test read
  # it, so a model that left estimate_sac_range() without the flag would
  # have brought the defect back unnoticed.
  # An accepted range.  Its fit converged, or it would have been refused.
  acc <- suppressWarnings(estimate_sac_range(cpn_field(seed = 2), "z", seed = 1))
  expect_true(is.data.frame(attr(acc, "variogram_model")))
  v <- expect_gstat_verdict(acc)
  if (is.finite(acc)) expect_true(v)
  # Refused ranges: the trend, whose fit stops at the iteration limit where
  # these tests were written, and the long range, whose fit converges there.
  # Either verdict is one TRUE or FALSE, and a refusal recorded as "did not
  # converge" cannot sit on a model that says it did.
  for (args in list(list(slope = 0.02, seed = 1),
                    list(eff_range = 5000, nugget = 0.05, seed = 1))) {
    ref <- suppressWarnings(estimate_sac_range(do.call(cpn_field, args), "z", seed = 1))
    expect_true(is.data.frame(attr(ref, "variogram_model")))
    v <- expect_gstat_verdict(ref)
    if (identical(attr(ref, "rejected_reason"), NOT_CONV)) expect_false(v)
  }
})


test_that("a model built from REML parameters carries no verdict", {
  skip_if_not_installed("gstat")
  skip_if_not_installed("nlme")
  # With detrend = "reml" the model is put together from the REML estimates
  # (gstat::vgm()), not fitted by gstat, so there is no optimiser's verdict
  # to carry.  estimate_sac_range() documents that it has none, and
  # .vgm_converged() reads NA, which leaves a refused REML range to be judged
  # by its reason.
  pts <- cpn_field(seed = 2)
  pts$w <- sf::st_coordinates(pts)[, 1] / 1000
  sac <- suppressWarnings(estimate_sac_range(pts, "z", "w", detrend = "reml", seed = 1))
  # A REML fit that fails falls back to OLS detrending and a model gstat
  # fitted, which carries its verdict like any other.
  if (!identical(attr(sac, "detrend_method"), "reml")) {
    expect_gstat_verdict(sac)
    skip("the REML fit did not succeed here")
  }
  vm <- attr(sac, "variogram_model")
  expect_true(is.data.frame(vm))
  expect_null(attr(vm, "converged", exact = TRUE))
  expect_identical(spatialkit:::.vgm_converged(vm), NA)
})


test_that("a fit gstat gave up on, past the lags, gives cp no nugget on any platform", {
  skip_if_not_installed("gstat")
  # The whole chain with nothing left to the optimiser: gstat warns that it
  # did not converge, estimate_sac_range() refuses the range as past the
  # largest lag and marks the model, and the profile reads the mark.  Nothing
  # here is skipped, so a verdict that is lost or wrong on the way fails.
  pts <- cpn_field(seed = 2)
  local_mocked_bindings(fit.variogram = cpn_fit_stub(converges = FALSE), .package = "gstat")
  sac <- suppressWarnings(estimate_sac_range(pts, "z", seed = 1))
  expect_true(is.na(sac))
  expect_identical(attr(sac, "rejected_reason"), PAST_LAGS)
  expect_identical(attr(attr(sac, "variogram_model"), "converged", exact = TRUE), FALSE)
  expect_equal(sac_nugget(sac), 0.4)

  out <- cpn_warnings(resolution_profile(pts, "z", sac = sac, n_levels = 4, nstart = 2))
  expect_cp_from_finest_level(out$value)
  expect_true(any(grepl(paste0("`sac` reports no usable range (", PAST_LAGS,
                               "), and its variogram model did not converge; ",
                               "its nugget and sill are where the optimiser stopped"),
                        out$warnings, fixed = TRUE)))
  expect_false(any(grepl("still takes its nugget", out$warnings, fixed = TRUE)))

  # The same for the variogram the profile estimates itself.
  out <- cpn_warnings(resolution_profile(pts, "z", n_levels = 4, nstart = 2))
  expect_identical(attr(attr(out$value, "sac"), "rejected_reason"), PAST_LAGS)
  expect_cp_from_finest_level(out$value)
  expect_true(any(grepl(paste0("the variogram estimated from `response_var` ",
                               "(attr(<profile>, \"sac\")) reports no usable range (",
                               PAST_LAGS, "), and its variogram model did not converge"),
                        out$warnings, fixed = TRUE)))
})


test_that("the same fit without gstat's warning still gives cp its nugget", {
  skip_if_not_installed("gstat")
  # The counterpart, and the documented behaviour that does not change: the
  # model is the one above and so is the refusal, but the optimiser did not
  # give up, so the nugget was fitted at the short lags and serves Cp.
  pts <- cpn_field(seed = 2)
  local_mocked_bindings(fit.variogram = cpn_fit_stub(converges = TRUE), .package = "gstat")
  sac <- suppressWarnings(estimate_sac_range(pts, "z", seed = 1))
  expect_true(is.na(sac))
  expect_identical(attr(sac, "rejected_reason"), PAST_LAGS)
  expect_identical(attr(attr(sac, "variogram_model"), "converged", exact = TRUE), TRUE)

  out <- cpn_warnings(resolution_profile(pts, "z", sac = sac, n_levels = 4, nstart = 2))
  p <- out$value
  expect_identical(attr(p, "cp_noise")$source, "variogram nugget")
  expect_equal(attr(p, "cp_noise")$value, 0.4)
  expect_true(all(is.na(p$reliability)))
  expect_true(any(grepl(paste0("`sac` reports no usable range (", PAST_LAGS,
                               "). `cp` still takes its nugget (0.4)"),
                        out$warnings, fixed = TRUE)))
  expect_false(any(grepl("did not converge", out$warnings)))
})


test_that("a fit that did not converge, refused as past the lags, gives cp no nugget", {
  skip_if_not_installed("gstat")
  # A strong east-west trend on an exponential field (effective range 150,
  # nugget 0.1).  The all-pairs variogram rises without a sill, gstat stops at
  # its iteration limit from every start with a range of about 29,000 against
  # a largest lag of 647, and the refusal is recorded as past the largest lag.
  # The nugget where the optimiser halted is 1.29, more than the nugget and
  # the sill of the field under the trend together (1.1).  It used to be Cp's
  # noise variance.
  pts <- cpn_field(slope = 0.02, seed = 1)
  sac <- suppressWarnings(estimate_sac_range(pts, "z", seed = 1))
  # The verdict has to be there (see expect_gstat_verdict()).  Whether gstat's
  # optimiser reaches its limit depends on the LAPACK build (see
  # test-sac-range.R), so this draw is taken only where it ends that way; the
  # stubbed fit above and the hand-made sac further down are the same case on
  # every platform.
  expect_gstat_verdict(sac)
  skip_if_not(is.na(sac) &&
                identical(attr(sac, "rejected_reason"), PAST_LAGS) &&
                isFALSE(attr(attr(sac, "variogram_model"), "converged")),
              "gstat's all-pairs fit did not end as 'not converged, past the largest lag' here")
  expect_gt(sac_nugget(sac), 0)

  out <- cpn_warnings(resolution_profile(pts, "z", sac = sac, n_levels = 6, nstart = 3))
  p <- out$value
  expect_cp_from_finest_level(p)
  # With no subsample and every row scored, Cp is Mallows' RSS / n + 2 s2 L / n
  # with s2 the residual mean square of the level it was read at.
  cn <- attr(p, "cp_noise"); n <- nrow(pts)
  j <- which(p$levels == cn$level)
  expect_equal(cn$value, p$rss[j] / (n - cn$level))
  expect_equal(p$cp, p$rss / n + cn$value * 2 * p$levels / n)
  # The warning names both facts: why the range was refused, and that the
  # fit did not converge.
  expect_true(any(grepl(paste0("`sac` reports no usable range (", PAST_LAGS,
                               "), and its variogram model did not converge; ",
                               "its nugget and sill are where the optimiser stopped"),
                        out$warnings, fixed = TRUE)))
  expect_false(any(grepl("still takes its nugget", out$warnings, fixed = TRUE)))
  expect_true(any(grepl("`sac` gives no usable variogram model, so `cp` takes its noise variance",
                        out$warnings, fixed = TRUE)))
})


test_that("the same holds for the variogram the profile estimates itself", {
  skip_if_not_installed("gstat")
  pts <- cpn_field(slope = 0.02, seed = 1)
  out <- cpn_warnings(resolution_profile(pts, "z", n_levels = 6, nstart = 3))
  p <- out$value
  sac <- attr(p, "sac")
  expect_gstat_verdict(sac)
  skip_if_not(is.na(sac) &&
                identical(attr(sac, "rejected_reason"), PAST_LAGS) &&
                isFALSE(attr(attr(sac, "variogram_model"), "converged")),
              "gstat's all-pairs fit did not end as 'not converged, past the largest lag' here")
  expect_cp_from_finest_level(p)
  expect_true(any(grepl(paste0("the variogram estimated from `response_var` ",
                               "(attr(<profile>, \"sac\")) reports no usable range (",
                               PAST_LAGS, "), and its variogram model did not converge"),
                        out$warnings, fixed = TRUE)))
})


test_that("the model's own verdict decides, whatever the reason says", {
  # The case above made by hand, so that it does not rest on how gstat's
  # optimiser ends on a given platform.
  pts <- cpn_field(seed = 2)
  for (reason in c(PAST_LAGS, DECREASING)) {
    out <- cpn_warnings(resolution_profile(pts, "z", sac = cpn_refused(reason, converged = FALSE),
                                           n_levels = 4, nstart = 2))
    expect_cp_from_finest_level(out$value)
    expect_true(any(grepl(paste0("reports no usable range (", reason,
                                 "), and its variogram model did not converge; ",
                                 "its nugget and sill are where the optimiser stopped"),
                          out$warnings, fixed = TRUE)), info = reason)
    expect_false(any(grepl("still takes its nugget", out$warnings, fixed = TRUE)), info = reason)
  }
  # A reason that already says so is not said twice.
  out <- cpn_warnings(resolution_profile(pts, "z", sac = cpn_refused(NOT_CONV, converged = FALSE),
                                         n_levels = 4, nstart = 2))
  expect_cp_from_finest_level(out$value)
  expect_true(any(grepl(paste0("reports no usable range (", NOT_CONV,
                               "); its nugget and sill are where the optimiser stopped"),
                        out$warnings, fixed = TRUE)))
  expect_length(grep("did not converge", out$warnings), 1L)
  expect_false(any(grepl("and its variogram model did not converge", out$warnings, fixed = TRUE)))
})


test_that("with no level to read a noise variance at, cp is NA and the warning allows for it", {
  # The refusal sends cp to the finest level whose cells hold at least two
  # scored rows on average.  With 200 and 250 cells on 300 points no level
  # does: cp is NA at every level and nothing is recorded as its noise
  # variance.  The warning closed on "takes its noise variance from the
  # finest level" alone, which is untrue here, and a sac refused as past the
  # lags, which used to give cp its nugget at both levels, now comes this way.
  pts <- cpn_field(seed = 2)
  sac <- cpn_refused(PAST_LAGS, converged = FALSE)
  out <- cpn_warnings(resolution_profile(pts, "z", sac = sac, levels = c(200, 250),
                                         min_cell_n = 1, nstart = 2))
  p <- out$value
  expect_true(all(is.finite(p$rss)))
  expect_true(all(is.na(p$cp)))
  expect_null(attr(p, "cp_noise"))
  expect_true(any(grepl(paste0("`cp` takes its noise variance from the finest level whose ",
                               "cells hold at least two scored rows on average, and is NA ",
                               "at every level when none does."),
                        out$warnings, fixed = TRUE)))
  expect_false(any(grepl("residual mean square of the finest level", out$warnings, fixed = TRUE)))
  # One level that does hold two rows a cell gives every level its cp.
  out <- cpn_warnings(resolution_profile(pts, "z", sac = sac, levels = c(150, 200, 250),
                                         min_cell_n = 1, nstart = 2))
  expect_identical(attr(out$value, "cp_noise")$level, 150L)
  expect_true(all(is.finite(out$value$cp)))
})


test_that("a converged fit refused as past the lags still gives cp its nugget", {
  # Unchanged: the nugget of a fit that converged was fitted at the short
  # lags, and it is the range, not the nugget, that the refusal is about.
  pts <- cpn_field(seed = 2)
  for (reason in c(PAST_LAGS, DECREASING)) {
    out <- cpn_warnings(resolution_profile(pts, "z", sac = cpn_refused(reason, converged = TRUE),
                                           n_levels = 4, nstart = 2))
    p <- out$value
    expect_identical(attr(p, "cp_noise")$source, "variogram nugget", info = reason)
    expect_equal(attr(p, "cp_noise")$value, 0.5)
    expect_equal(attr(p, "variogram")$nugget, 0.5)
    expect_equal(p$cp, p$rss / nrow(pts) + 0.5 * 2 * p$levels / nrow(pts))
    expect_true(all(is.na(p$reliability)))
    expect_true(any(grepl(paste0("reports no usable range (", reason,
                                 "). `cp` still takes its nugget (0.5)"),
                          out$warnings, fixed = TRUE)), info = reason)
    expect_false(any(grepl("did not converge", out$warnings)), info = reason)
  }
})


test_that("a converged fit of gstat's own, past the lags, still gives cp its nugget", {
  skip_if_not_installed("gstat")
  # An effective range of 5000 on a 1000 m square: the fit converges to a
  # range near 5000, past the largest lag of 647, with a nugget of 0.049
  # where the field's is 0.05.
  pts <- cpn_field(eff_range = 5000, nugget = 0.05, seed = 1)
  sac <- suppressWarnings(estimate_sac_range(pts, "z", seed = 1))
  expect_gstat_verdict(sac)
  skip_if_not(is.na(sac) &&
                identical(attr(sac, "rejected_reason"), PAST_LAGS) &&
                isTRUE(attr(attr(sac, "variogram_model"), "converged")),
              "gstat's all-pairs fit did not end as 'converged, past the largest lag' here")
  out <- cpn_warnings(resolution_profile(pts, "z", sac = sac, n_levels = 6, nstart = 3))
  p <- out$value
  expect_identical(attr(p, "cp_noise")$source, "variogram nugget")
  expect_equal(attr(p, "cp_noise")$value, sac_nugget(sac))
  expect_true(all(is.na(p$reliability)))
  expect_true(any(grepl("`cp` still takes its nugget", out$warnings, fixed = TRUE)))
  expect_false(any(grepl("did not converge", out$warnings)))
})


test_that("a sac whose model carries no verdict is judged by its reason", {
  # Made by hand, built from REML parameters, or saved before the model
  # carried the flag: the reason is all there is, and it is read as before.
  pts <- cpn_field(seed = 2)
  out <- cpn_warnings(resolution_profile(pts, "z", sac = cpn_refused(NOT_CONV),
                                         n_levels = 4, nstart = 2))
  expect_cp_from_finest_level(out$value)
  expect_true(any(grepl(paste0("reports no usable range (", NOT_CONV,
                               "); its nugget and sill are where the optimiser stopped"),
                        out$warnings, fixed = TRUE)))
  for (reason in c(PAST_LAGS, DECREASING)) {
    out <- cpn_warnings(resolution_profile(pts, "z", sac = cpn_refused(reason),
                                           n_levels = 4, nstart = 2))
    expect_identical(attr(out$value, "cp_noise")$source, "variogram nugget", info = reason)
    expect_true(any(grepl("`cp` still takes its nugget (0.5)", out$warnings, fixed = TRUE)),
                info = reason)
    expect_false(any(grepl("did not converge", out$warnings)), info = reason)
  }
  # A reason that says the fit did not converge is believed even against a
  # flag that says it did: the two cannot disagree in a sac that
  # estimate_sac_range() made, and of two accounts the cautious one is taken.
  out <- cpn_warnings(resolution_profile(pts, "z", sac = cpn_refused(NOT_CONV, converged = TRUE),
                                         n_levels = 4, nstart = 2))
  expect_cp_from_finest_level(out$value)
})
