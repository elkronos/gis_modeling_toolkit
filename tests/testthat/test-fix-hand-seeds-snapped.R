# tests/testthat/test-fix-hand-seeds-snapped.R
# ---------------------------------------------------------------------------
# Points placed by hand, used as seeds, one a hair outside a boundary that is
# not convex.
#
# build_tessellation() indexes points that are not seeds by the cell that
# holds them, and gives one that lies in no cell but within a thousandth of
# the median cell width of a cell the ID of the nearest cell.  Beside the
# notch of a square that can be another seed's cell, so two seeds share one
# ID in `index`, and cell_values_at(seeds = ) refuses them.  Its error told
# the caller to set the column from `index`, which had been done; it now says
# why `index` repeats the ID and names the route that gives each its cell.
# ---------------------------------------------------------------------------

.hs_layout <- function() {
  U <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(550, 1000), c(550, 300),
    c(450, 300), c(450, 1000), c(0, 1000), c(0, 0)))), crs = 32617))
  # The sixth point is 1e-4 above the notch's floor, a hair from the fifth
  # point's cell; its own cell is the notch's.
  p <- sf::st_as_sf(data.frame(x = c(150, 850, 150, 850, 500, 500),
                               y = c(150, 150, 650, 650, 299.99995, 300.0001)),
                    coords = c("x", "y"), crs = 32617)
  withr::local_seed(3)
  obs <- sf::st_as_sf(data.frame(x = stats::runif(300, 0, 1000), y = stats::runif(300, 0, 1000)),
                      coords = c("x", "y"), crs = 32617)
  obs <- obs[lengths(sf::st_intersects(obs, U)) > 0, ]
  obs$z <- sf::st_coordinates(obs)[, 1] / 100 + stats::rnorm(nrow(obs))
  newp <- sf::st_as_sf(data.frame(x = stats::runif(200, 0, 1000), y = stats::runif(200, 0, 1000)),
                       coords = c("x", "y"), crs = 32617)
  list(U = U, p = p, obs = obs, newp = newp)
}

.hs_call <- function(L, seeds, polygons = TRUE) {
  tess <- suppressWarnings(build_tessellation(seeds, boundary = L$U, method = "voronoi", quiet = TRUE))
  seeds$cell_id <- tess$index
  a <- suppressWarnings(assign_features_to_polygons(L$obs, tess$cells, polygon_id_col = "cell_id",
                                                    quiet = TRUE))
  v <- suppressWarnings(summarize_by_cell(a, "z", id_col = "cell_id", cells_sf = tess$cells))
  if (!polygons) v <- sf::st_drop_geometry(v)
  list(tess = tess, run = function()
    suppressWarnings(cell_values_at(v, points = a, newdata = L$newp, response_var = "z",
                                    seeds = seeds)))
}

test_that("a point placed by hand a hair outside the notch takes its neighbour's cell, and the error says so and names the route that works", {
  L <- .hs_layout()
  expect_identical(lengths(sf::st_intersects(L$p, L$U)), c(1L, 1L, 1L, 1L, 1L, 0L))
  for (polygons in c(TRUE, FALSE)) {
    plain <- .hs_call(L, L$p, polygons)
    expect_identical(plain$tess$params$snapped$which, 6L)
    expect_identical(plain$tess$index[6L], plain$tess$index[5L])
    expect_error(plain$run(), paste0(
      "cell_values_at(): 2 seeds have the cell of another in `seeds$cell_id`; each seed ",
      "stands for one cell. Set the column from the tessellation built on these seeds ",
      "(seeds$cell_id <- tess$index). Points placed by hand are indexed by the cell that ",
      "holds them, and one that lies in no cell but within a thousandth of the median cell ",
      "width of another seed's cell takes that cell in `index`; passed through ",
      "get_voronoi_seeds(method = \"provided\", seeds = ) before build_tessellation(), each ",
      "takes the cell it generates."), fixed = TRUE, info = paste("polygons:", polygons))
    # The route the error names runs, and fills no location.
    named <- .hs_call(L, get_voronoi_seeds(method = "provided", seeds = L$p), polygons)
    expect_false(anyDuplicated(named$tess$index) > 0L)
    expect_identical(named$tess$index[1:5], plain$tess$index[1:5])
    got <- named$run()
    expect_identical(sum(got$filled), 0L)
  }
})

test_that("?cell_values_at and ?create_voronoi_polygons say a point placed by hand within the snapping distance takes the nearest cell's ID", {
  # The page in man/ from the source tree, the installed page under R CMD check.
  page <- function(name) {
    file <- testthat::test_path("..", "..", "man", paste0(name, ".Rd"))
    rd <- if (file.exists(file)) tools::parse_Rd(file)
          else tools::Rd_db("spatialkit")[[paste0(name, ".Rd")]]
    flat <- function(x) paste(unlist(lapply(x, function(el)
      if (is.character(el)) el else flat(el))), collapse = "")
    gsub("[[:space:]]+", " ", flat(rd))
  }
  txt <- page("cell_values_at")
  expect_match(txt, paste0("so one outside a boundary that is not convex is NA although its ",
                           "cell is built. One that lies in no cell but within a thousandth of ",
                           "the median cell width of a cell takes the ID of the nearest cell"),
               fixed = TRUE)
  expect_match(txt, paste0("one outside the boundary is NA, or the ID of the nearest cell, ",
                           "which can be the cell of another seed, when it lies within a ",
                           "thousandth of the median cell width of it"), fixed = TRUE)
  cvp <- page("create_voronoi_polygons")
  expect_match(cvp, paste0("so a point outside the study area is NA, unless it lies within a ",
                           "thousandth of the median cell width of a cell, whose ID it then ",
                           "takes"), fixed = TRUE)
})
