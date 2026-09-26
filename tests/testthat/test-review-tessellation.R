# ===========================================================================
# Regressions from the adversarial review of the tessellation slice:
#
#   * coerce_to_points(mode = "auto") segfaulted R on an EMPTY
#     MULTILINESTRING (st_cast() gives one empty part, and st_line_sample()
#     on it crashes in sf 1.0.x), and refused an EMPTY LINESTRING outright.
#   * build_tessellation(method = "triangles") handed raw UTM-sized
#     coordinates to qhull, which lost the precision to make nearby points
#     vertices.
# ===========================================================================


# ---------------------------------------------------------------------------
# Empty lines in coerce_to_points()
# ---------------------------------------------------------------------------

# Run `code` in a fresh R process that loads this same copy of spatialkit and
# return its exit status, its output and whatever it saved as `res`.  The
# regression being guarded is a SEGFAULT: run in this process, a relapse
# would kill the whole test run instead of failing one test.
.rt_run_child <- function(code) {
  ns_path   <- getNamespaceInfo(asNamespace("spatialkit"), "path")
  installed <- file.exists(file.path(ns_path, "Meta", "package.rds"))
  load_line <- if (installed) {
    sprintf("suppressPackageStartupMessages(library(spatialkit, lib.loc = %s))",
            deparse(dirname(ns_path)))
  } else {
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE))",
            deparse(ns_path))
  }
  script <- tempfile(fileext = ".R")
  out    <- tempfile(fileext = ".rds")
  on.exit(unlink(c(script, out)), add = TRUE)
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
               load_line,
               code,
               sprintf("saveRDS(res, %s)", deparse(out))),
             script)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                     c("--vanilla", shQuote(script)),
                                     stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  list(status = if (is.null(status)) 0L else as.integer(status),
       output = output,
       res    = if (file.exists(out)) readRDS(out) else NULL)
}


test_that("coerce_to_points turns empty lines into empty POINTs instead of crashing", {
  child <- .rt_run_child(c(
    "seg <- function(x0, y0, x1, y1) rbind(c(x0, y0), c(x1, y1))",
    "mls <- function(...) sf::st_multilinestring(list(...))",
    "try_ <- function(expr) tryCatch(expr, error = function(e) conditionMessage(e))",
    "res <- list()",
    # An empty MULTILINESTRING between two real ones, projected.
    "x <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(",
    "  mls(seg(0, 0, 10, 0)), sf::st_multilinestring(), mls(seg(0, 0, 0, 20)),",
    "  crs = 32617))",
    "res$mls_proj <- try_(coerce_to_points(x, 'auto'))",
    # The same in lon/lat, projected temporarily (the default).
    "x_ll <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(",
    "  mls(seg(-80, 35, -79, 35)), sf::st_multilinestring(),",
    "  mls(seg(-80, 35, -80, 36)), crs = 4326))",
    "res$mls_ll <- try_(suppressWarnings(coerce_to_points(x_ll, 'auto')))",
    # A layer of one empty row.
    "res$mls_one <- try_(coerce_to_points(x[2, ], 'auto'))",
    # An empty PART beside a real one: the empty part used to be sampled
    # whenever it came first among parts of equal (zero) length.
    "x_part <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(",
    "  sf::st_multilinestring(list(matrix(numeric(0), 0, 2), seg(0, 0, 10, 0))),",
    "  sf::st_multilinestring(list(matrix(numeric(0), 0, 2), seg(3, 4, 3, 4))),",
    "  crs = 32617))",
    "res$mls_part <- try_(coerce_to_points(x_part, 'auto'))",
    # An empty LINESTRING, in both modes that sample lines.
    "x_ls <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(",
    "  sf::st_linestring(seg(0, 0, 10, 0)), sf::st_linestring(),",
    "  sf::st_linestring(seg(0, 0, 0, 20)), crs = 32617))",
    "res$ls_auto <- try_(coerce_to_points(x_ls, 'auto'))",
    "res$ls_mid  <- try_(coerce_to_points(x_ls, 'line_midpoint'))",
    # The documented cleaning path: prep_model_data() drops empty rows, but
    # it pointizes first, so this is where a GeoPackage null geometry crashed.
    "d <- x; d$y <- c(1, 2, 3); d$x1 <- c(0.1, 0.5, 0.9)",
    "res$prep <- try_(prep_model_data(d, 'y', 'x1'))"
  ))

  expect_equal(child$status, 0L,
               info = paste(utils::tail(child$output, 25), collapse = "\n"))
  res <- child$res
  expect_type(res, "list")
  if (!is.list(res)) return(invisible())

  xy <- function(p) unname(sf::st_coordinates(p)[, 1:2, drop = FALSE])
  for (nm in c("mls_proj", "mls_ll", "mls_one", "mls_part", "ls_auto",
               "ls_mid", "prep")) {
    expect_s3_class(res[[nm]], "sf")
  }
  if (!all(vapply(res, inherits, logical(1), "sf"))) return(invisible())

  # Row for row with the input: the empty feature is an empty POINT in its
  # own row (st_coordinates() gives it an NA row), and its neighbours get
  # their true midpoints.
  for (nm in c("mls_proj", "ls_auto", "ls_mid")) {
    got <- res[[nm]]
    expect_equal(got$id, 1:3, info = nm)
    expect_true(all(sf::st_geometry_type(got) == "POINT"), info = nm)
    expect_equal(sf::st_is_empty(got), c(FALSE, TRUE, FALSE), info = nm)
    expect_equal(xy(got)[c(1, 3), ], rbind(c(5, 0), c(0, 10)), info = nm)
  }

  ll <- res$mls_ll
  expect_equal(sf::st_is_empty(ll), c(FALSE, TRUE, FALSE))
  expect_equal(sf::st_crs(ll), sf::st_crs(4326))
  expect_equal(xy(ll)[c(1, 3), ], rbind(c(-79.5, 35), c(-80, 35.5)),
               tolerance = 1e-4)

  expect_equal(nrow(res$mls_one), 1L)
  expect_true(sf::st_is_empty(res$mls_one))

  # The real part is the one sampled; a zero-length part gives its point.
  expect_false(any(sf::st_is_empty(res$mls_part)))
  expect_equal(xy(res$mls_part), rbind(c(5, 0), c(3, 4)))

  # prep_model_data() drops the empty row and says why.
  expect_equal(nrow(res$prep), 2L)
  expect_equal(res$prep$y, c(1, 3))
  expect_equal(attr(res$prep, "dropped")$n_geometry, 1L)
  expect_equal(attr(res$prep, "dropped")$which, 2L)
})


# ---------------------------------------------------------------------------
# Delaunay triangles at UTM-sized coordinates
# ---------------------------------------------------------------------------

.rt_utm_points <- function(n, x0, y0, ext = 100, seed = 42) {
  set.seed(seed)
  sf::st_as_sf(data.frame(x = x0 + stats::runif(n, 0, ext),
                          y = y0 + stats::runif(n, 0, ext)),
               coords = c("x", "y"), crs = 32632)
}

# One key per triangle from its vertex coordinates, independent of the order
# its vertices are listed in.
.rt_tri_key <- function(cells) {
  vapply(seq_len(nrow(cells)), function(i) {
    m <- sf::st_coordinates(cells[i, ])[1:3, 1:2]
    paste(sort(sprintf("%.17g_%.17g", m[, 1], m[, 2])), collapse = "|")
  }, character(1))
}


test_that("triangles makes every point a vertex at UTM-sized coordinates", {
  skip_if_not_installed("geometry")
  # qhull lifts points onto x^2 + y^2; with raw northings near 5e6 (and 9e6
  # south of the equator) points a few metres apart were dropped as coplanar.
  # 200 points over 100 m gave a few dozen triangles instead of ~386, with
  # every point still indexed into one, so nothing looked wrong.
  n <- 200L
  for (off in list(c(5e5, 5e6), c(5e5, 9.9e6))) {
    pts <- .rt_utm_points(n, off[1], off[2])
    res <- build_tessellation(pts, method = "triangles", quiet = TRUE)
    lab <- sprintf("offset (%g, %g)", off[1], off[2])

    # Euler: a Delaunay triangulation of n points in general position with h
    # of them on the hull has 2n - h - 2 triangles.
    h <- nrow(sf::st_coordinates(
      sf::st_convex_hull(sf::st_union(sf::st_geometry(pts))))) - 1L
    expect_equal(nrow(res$cells), 2L * n - h - 2L, info = lab)

    # Every input point is a vertex, at its own coordinates exactly: the
    # centring is undone by building rings from the original coordinates,
    # not by adding the shift back.
    key <- function(m) sprintf("%.17g_%.17g", m[, 1], m[, 2])
    vtx <- key(sf::st_coordinates(res$cells)[, 1:2, drop = FALSE])
    inp <- key(sf::st_coordinates(pts))
    expect_true(all(inp %in% vtx), info = lab)
    expect_true(all(vtx %in% inp), info = lab)

    # The same points moved to the origin triangulate identically.
    shifted <- sf::st_set_geometry(pts, sf::st_geometry(pts) - off)
    sf::st_crs(shifted) <- sf::st_crs(pts)
    ref <- build_tessellation(shifted, method = "triangles", quiet = TRUE)
    expect_equal(nrow(res$cells), nrow(ref$cells), info = lab)
    expect_equal(sum(as.numeric(sf::st_area(res$cells))),
                 sum(as.numeric(sf::st_area(ref$cells))), tolerance = 1e-8,
                 info = lab)
  }
})


test_that("triangle cell_ids do not depend on the input row order", {
  skip_if_not_installed("geometry")
  # Not a regression of the old code (qhull's output order does not follow
  # the input order for points in general position), but the centring must
  # not break it: the shift is the bbox midpoint, which a permutation leaves
  # bit-identical, where a mean could differ in its last bits.
  pts  <- .rt_utm_points(150L, 5e5, 5e6, seed = 7)
  set.seed(8)
  perm <- sample(nrow(pts))
  a <- build_tessellation(pts, method = "triangles", quiet = TRUE)
  b <- build_tessellation(pts[perm, ], method = "triangles", quiet = TRUE)
  expect_identical(.rt_tri_key(a$cells), .rt_tri_key(b$cells))
  expect_identical(a$cells$cell_id, b$cells$cell_id)
  expect_identical(a$index[perm], b$index)
})
