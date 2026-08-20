# ===========================================================================
# predict() must express newdata in the CRS the model was fitted in.
#
# Left alone, ensure_projected() picks a projection from newdata's own
# centroid, which for lon/lat input is whatever UTM zone the new points
# happen to fall in -- not the training zone. The backend then mixes two
# coordinate systems and returns plausible-looking wrong numbers.
# See test-predict-na-alignment.R for row alignment.
# ===========================================================================

test_that("predict.gwr_fit aligns newdata CRS to the training CRS", {
  skip_if_not_installed("sf")
  skip_if_not_installed("sp")
  skip_if_not_installed("GWmodel")

  set.seed(31)
  n <- 40
  px <- runif(n, 0, 5000)
  py <- runif(n, 0, 5000)
  df <- data.frame(
    y  = 2 + 0.001 * px + 0.5 * rnorm(n),
    x1 = rnorm(n),
    px = px, py = py
  )
  train_sf <- sf::st_as_sf(df, coords = c("px", "py"), crs = 32632)

  fit <- fit_gwr_model(train_sf, "y", "x1",
                       adaptive = TRUE, bandwidth = 20)

  newdata <- train_sf[1:10, ]
  pred_same <- predict(fit, newdata = newdata)

  # Same locations expressed in lon/lat: before the fix, ensure_projected()
  # auto-picked a UTM zone from the newdata centroid (zone 31 here, not the
  # training zone 32), silently mixing coordinate systems in gwr.predict().
  newdata_ll <- sf::st_transform(newdata, 4326)
  pred_ll <- predict(fit, newdata = newdata_ll)
  expect_equal(pred_ll, pred_same, tolerance = 1e-4)

  # Same locations in a different *projected* CRS: before the fix this was
  # left untouched (already projected), producing garbage distances.
  newdata_merc <- sf::st_transform(newdata, 3857)
  pred_merc <- predict(fit, newdata = newdata_merc)
  expect_equal(pred_merc, pred_same, tolerance = 1e-4)
})
