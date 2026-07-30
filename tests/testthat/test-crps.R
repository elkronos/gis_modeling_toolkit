# ---------------------------------------------------------------------------
# Regression tests for .crps_energy()
#
# Guards against the cv_bayes() bug where the CRPS block called the
# non-existent matrixStats::colSort(), erroring inside try() and silently
# dropping CRPS/coverage extras from every fold whenever matrixStats was
# installed.  The computation now lives in a testable internal helper.
# ---------------------------------------------------------------------------

test_that(".crps_energy matches the naive O(m^2) energy-form CRPS", {
  set.seed(11)
  m <- 200
  n <- 7
  draws <- matrix(rnorm(m * n, mean = rep(seq_len(n), each = m)), nrow = m)
  y <- rnorm(n, mean = seq_len(n))

  crps_fast <- spatialkit:::.crps_energy(draws, y)

  # Naive reference: CRPS = E|X - y| - 0.5 * E|X - X'| over all m^2 pairs
  crps_naive <- vapply(seq_len(n), function(j) {
    x <- draws[, j]
    mean(abs(x - y[j])) - 0.5 * mean(abs(outer(x, x, "-")))
  }, numeric(1))

  expect_equal(crps_fast, crps_naive, tolerance = 1e-10)
  expect_true(all(crps_fast > 0))
})


test_that(".crps_energy handles degenerate (point-mass) predictive draws", {
  # All draws equal to a constant c: CRPS reduces to |c - y|
  draws0 <- matrix(5, nrow = 50, ncol = 3)
  expect_equal(spatialkit:::.crps_energy(draws0, c(5, 5, 5)), c(0, 0, 0))
  expect_equal(spatialkit:::.crps_energy(draws0, c(3, 5, 8)), c(2, 0, 3))
})


test_that(".crps_energy sorts within columns, not across the whole matrix", {
  # Columns with interleaved magnitudes: a whole-matrix sort (the failure
  # mode of a wrong vectorised sort) would mix values across observations.
  draws <- cbind(c(30, 10, 20), c(2, 3, 1))
  y <- c(20, 2)
  crps_naive <- vapply(1:2, function(j) {
    x <- draws[, j]
    mean(abs(x - y[j])) - 0.5 * mean(abs(outer(x, x, "-")))
  }, numeric(1))
  expect_equal(spatialkit:::.crps_energy(draws, y), crps_naive, tolerance = 1e-12)
})
