# cell_values_at() without seeds, a location the same distance from two
# scored locations of one lattice column (BUG-199): the scored locations are
# ordered by column (.near_columns()), then y, so the tie goes to the lower
# point in any unit or origin, as ?cell_values_at says.

.vc_crs_m    <- function() sf::st_crs(26915)
.vc_crs_usft <- function() sf::st_crs("+proj=utm +zone=15 +datum=NAD83 +units=us-ft +no_defs")
.vc_crs_m01  <- function() sf::st_crs(paste0(
  "+proj=tmerc +lat_0=0 +lon_0=-93 +k=0.9996 +x_0=500000.1 +y_0=0.1 +datum=NAD83 ",
  "+units=m +no_defs"))
.vc_unit <- c(m = 1, usft = 3937 / 1200, m01 = 1)

.vc_quiet <- function(expr) suppressWarnings(suppressMessages(
  logger::with_log_threshold(expr, threshold = logger::FATAL, namespace = "spatialkit",
                             index = 2)))

.vc_square <- function(x0, y0, s) sf::st_polygon(list(rbind(
  c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))

# A 20 x 15 lattice at 250 m in cells of 2 x 2 points, and 14 locations 200 m
# west of its first column, each midway between two of its points.
.vc_case <- function(jitter = 0) {
  g <- expand.grid(i = 0:19, j = 0:14)
  x <- 4.5e5 + g$i * 250; y <- 4.55e6 + g$j * 250
  if (jitter > 0) {
    x <- x + withr::with_seed(99, stats::runif(length(x), -jitter, jitter))
    y <- y + withr::with_seed(98, stats::runif(length(y), -jitter, jitter))
  }
  z <- withr::with_seed(1, sin((x - min(x)) / 1500) + cos((y - min(y)) / 1250) +
                          stats::rnorm(length(x), 0, 0.5))
  pts <- sf::st_as_sf(data.frame(x = x, y = y, z = z, cell_id = paste(g$i %/% 2, g$j %/% 2)),
                      coords = c("x", "y"), crs = .vc_crs_m())
  q <- sf::st_as_sf(data.frame(x = 4.5e5 - 200, y = 4.55e6 + 125 + (0:13) * 250),
                    coords = c("x", "y"), crs = .vc_crs_m())
  bl <- expand.grid(ci = 0:9, cj = 0:7)
  cells <- sf::st_sf(cell_id = paste(bl$ci, bl$cj), geometry = sf::st_sfc(lapply(
    seq_len(nrow(bl)), function(k) .vc_square(4.5e5 + 500 * bl$ci[k] - 125,
                                              4.55e6 + 500 * bl$cj[k] - 125, 500)),
    crs = .vc_crs_m()))
  list(pts = pts, q = q, cells = cells)
}

.vc_run <- function(case, crs, order = NULL) {
  p <- sf::st_transform(case$pts, crs)
  if (!is.null(order)) p <- p[order, ]
  sm <- .vc_quiet(summarize_by_cell(p, "z", id_col = "cell_id",
                                     cells_sf = sf::st_transform(case$cells, crs),
                                     estimate = "mean"))
  .vc_quiet(cell_values_at(sf::st_transform(case$q, crs), sm, p, "z", id_col = "cell_id"))
}

test_that("cell_values_at() fills a tied location from the lower point of its column in every frame", {
  case <- .vc_case()
  crs <- list(m = .vc_crs_m(), usft = .vc_crs_usft(), m01 = .vc_crs_m01())
  res <- lapply(crs, function(cr) .vc_run(case, cr))
  # Each location lies midway between rows j and j + 1 of the first column;
  # the lower row j is in the cell "0 <j %/% 2>".
  want <- paste(0, (0:13) %/% 2)
  for (f in names(crs)) {
    # Before, with the origin moved, 2 of the 14 took the upper point.
    expect_identical(res[[f]]$filled_from, want, info = f)
    expect_identical(res[[f]]$value, res$m$value, info = f)
    expect_equal(res[[f]]$dist_nearest / .vc_unit[[f]], res$m$dist_nearest, tolerance = 1e-9,
                 info = f)
    expect_equal(res[[f]]$dist_ratio, res$m$dist_ratio, tolerance = 1e-9, info = f)
  }
  # The rows of `points` in another order, and the search without FNN.
  sh <- .vc_run(case, crs$m01, order = withr::with_seed(3, sample(nrow(case$pts))))
  expect_identical(sh$filled_from, want)
  # The cell means are summed in another order, so to rounding.
  expect_equal(sh$value, res$m$value, tolerance = 1e-12)
  local_mocked_bindings(.has_backend = function(pkg, ...) FALSE)
  for (f in names(crs)) {
    nf <- .vc_run(case, crs[[f]])
    expect_identical(nf$filled_from, want, info = f)
    expect_identical(nf$value, res$m$value, info = f)
  }
})

test_that("cell_values_at() fills from the nearest scored location of a jittered lattice, as before", {
  case <- .vc_case(jitter = 60)
  res <- .vc_run(case, .vc_crs_m())
  # The nearest scored point of each location, by brute force.
  d <- sf::st_distance(case$q, case$pts)
  near <- apply(matrix(as.numeric(d), nrow(case$q)), 1L, which.min)
  expect_identical(res$filled_from, case$pts$cell_id[near])
  expect_equal(res$dist_nearest, apply(matrix(as.numeric(d), nrow(case$q)), 1L, min),
               tolerance = 1e-9)
  r2 <- .vc_run(case, .vc_crs_m01())
  expect_identical(r2$filled_from, res$filled_from)
})

test_that("?cell_values_at states the column rule for tied scored locations", {
  file <- test_path("..", "..", "man", "cell_values_at.Rd")
  rd <- if (file.exists(file)) tools::parse_Rd(file) else tools::Rd_db("spatialkit")[["cell_values_at.Rd"]]
  txt <- gsub("[[:space:]]+", " ", paste(unlist(rd), collapse = ""))
  expect_match(txt, paste0("takes the one in the lower column (x values within one part in ten ",
                           "million of the diagonal of the scored locations' bounding box, or ",
                           "within one part in 1e10 of their largest coordinate (at most a ",
                           "millionth of that diagonal, and at least four times the ",
                           "spacing of doubles at that coordinate) where that is more, counting ",
                           "as one), then the one lower in y, y values read the ",
                           "same way"), fixed = TRUE)
})
