# =============================================================================
# 10  A Bayesian GP and the coverage of its prediction intervals
# =============================================================================
#   source(system.file("scripts", "10-bayes.R", package = "spatialkit"))
#
# The slowest script of the tour. It runs four fits, one on its own and one
# for each of three cross-validation folds, and Stan compiles the model for
# each, which takes most of the time. Everything runs on a 100-point subset to
# keep the sampling short.
# =============================================================================
.tour_dir <- if (nzchar(system.file("scripts", package = "spatialkit")))
  system.file("scripts", package = "spatialkit") else "."
source(file.path(.tour_dir, "_common.R"))

if (skip_without("brms", "the whole of script 10")) {
  cat("  brms brings rstan, which compiles each model and needs a C++ toolchain.\n",
      "  cmdstanr is used instead when a CmdStan build is installed.\n", sep = "")
} else {

pts <- tour_points()
set.seed(3)
sub <- pts[sample(nrow(pts), 100), ]
cat(sprintf("  running on %d of %d points to keep the sampling tractable\n",
            nrow(sub), nrow(pts)))

step("10.1", "What length-scales the prior expects, before fitting")
# The GP is approximated by a basis expansion, and `gp_k` sets how many terms
# it gets per axis. With too few the model cannot represent short-range
# structure, whatever the data say.
# gp_lengthscale_bounds() gives the bounds the length-scale prior is calibrated
# over. The scale the basis resolves is a separate number: it depends on gp_k,
# and the fit reports it in 10.3. The lower bound is a fixed fraction of the
# pairwise distances (their 25th percentile divided by 2.45), so it does not
# move with n, and neither does the gp_k derived from it.
bounds <- gp_lengthscale_bounds(sf::st_coordinates(sub))
cat(sprintf("  length-scale prior calibrated over %.0f to %.0f m\n",
            bounds["lower"], bounds["upper"]))
cat("  A field whose true range sits below the scale the basis resolves needs a\n",
    "  bigger gp_k: neither these bounds nor the derived basis grow finer with\n",
    "  n, though denser sampling still helps the data identify a short range.\n",
    sep = "")

step("10.2", "Fit it")
t0 <- Sys.time()
bf <- fit_bayesian_spatial_model(sub, "z", "elev", chains = 2, iter = 600,
                                 gp_k = 12, seed = 1, compute_loo = TRUE)
cat(sprintf("\n  %.0f seconds including compilation\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
print(summary(bf))

step("10.3", "Read the warnings before the estimates")
# The fit reports on itself. A model that has not converged, or whose posterior
# length-scale is piled up against what the basis can resolve, is not a model
# whose intervals mean anything.
cat(sprintf("  convergence_ok : %s\n", isTRUE(bf$info$convergence_ok)))
# gp_ell_min is in the scaled coordinates brms was handed, not metres.
# Multiply back by the coordinate scaling to compare it with 10.1.
cs <- bf$info$coord_scaling
ell_crs <- bf$info$gp_ell_min * mean(c(cs$x_scale, cs$y_scale))
cat(sprintf("  basis          : gp_k %d -> %d functions\n",
            bf$info$gp_k, bf$info$gp_n_basis))
cat(sprintf("  finest scale   : %.2f in scaled units = %.0f m\n",
            bf$info$gp_ell_min, ell_crs))
# gp_k = 12 is below the 22 the rule derives for these 100 points when gp_k is
# NULL. The smaller basis keeps the fit fast, and may not reach the prior's
# lower bound.
cat(sprintf("  prior's lower  : %.0f m (10.1): %s\n", bounds["lower"],
            if (ell_crs > bounds["lower"])
              "finer than this basis resolves; a larger gp_k would reach it"
            else "within what this basis resolves"))
if (!is.null(bf$info$looic))
  cat(sprintf("  LOOIC          : %.1f\n", bf$info$looic))
cd <- bf$info$convergence_diagnostics
if (!is.null(cd) && length(cd))
  cat("  diagnostics    :", paste(utils::head(names(cd), 6), collapse = ", "), "\n")
cat("  The logged WARN line about length-scale draws below the resolvable scale\n",
    "  is the one to act on: raise gp_k and refit, or accept that the short-range\n",
    "  structure is outside this model's reach and say so.\n", sep = "")

step("10.4", "What it buys: an interval per prediction")
# Every other model in this tour returns a number. This one returns a
# distribution, which is only worth the wait if the distribution is honest.
# make_folds() warns, through cv_bayes(), that the blocks are narrower than the
# estimated range. With 100 points and three folds the blocks cannot be made
# wider and keep the folds: a sample this small cannot support both blocks as
# wide as the range and enough folds to estimate coverage.
cvb <- cv_bayes(sub, "z", "elev", k = 3, seed = 1,
                fit_args = list(chains = 2, iter = 600, gp_k = 12,
                                compute_loo = FALSE))
cat("\n  cross-validated coverage, by fold:\n")
cov_cols <- grep("^coverage_", names(cvb$fold_metrics), value = TRUE)
print(cvb$fold_metrics[, c("fold", "n_test", "RMSE", "CRPS", cov_cols)],
      row.names = FALSE, digits = 3)

step("10.5", "Is the uncertainty calibrated?")
# A 95 percent interval should contain the held-out value 95 percent of the
# time.
nominal  <- as.numeric(sub("^coverage_", "", cov_cols)) / 100
observed <- vapply(cov_cols, function(cc)
  stats::weighted.mean(cvb$fold_metrics[[cc]], cvb$fold_metrics$n_test), numeric(1))
for (j in seq_along(cov_cols))
  cat(sprintf("  nominal %.0f percent -> observed %.0f percent  (%s)\n",
              100 * nominal[j], 100 * observed[j],
              if (observed[j] < nominal[j] - 0.1) "too narrow"
              else if (observed[j] > nominal[j] + 0.1) "too wide" else "close"))
if (!skip_without("ggplot2", "the calibration plot")) {
  look_for("points on the 1:1 line. Below it the intervals are too narrow and ",
           "the model is overconfident; above it they are too wide to guide a ",
           "decision. Read the per-fold scatter as well as the average.")
  show_plot(plot_calibration(cvb), "10-calibration.png")
}

step("10.6", "Three folds of coverage is three numbers")
# Coverage is a proportion estimated from a few dozen held-out points per fold,
# so it varies from fold to fold by chance.
n_test <- mean(cvb$fold_metrics$n_test)
se95   <- sqrt(0.95 * 0.05 / n_test)
cat(sprintf("  %.0f test points per fold, so the standard error of a 95 percent\n", n_test))
cat(sprintf("  coverage estimate is about %.0f percentage points.\n", 100 * se95))
cat(sprintf("  One extra miss moves it by %.0f points on its own.\n",
            100 / n_test))
cat("  Read the calibration plot for gross miscalibration. For anything finer,\n",
    "  you need more folds and more points.\n", sep = "")

}
