# tests/testthat/test-core-docs.R
# -----------------------------------------------------------------------------
# What the help and the vignettes written for pooled cell values may say
# -----------------------------------------------------------------------------
# The help pages and the vignette sections that describe the pooled value,
# the interval for one observation, the size and the record of the count are
# read as text.  Five rules are held:
#
#   * none of them names existing units (counties, ZIP code areas, ZCTAs):
#     the pages do not say how the package's cells compare with any others;
#   * "calibrated" stands in them only in the sentence that denies it for
#     places away from the observations, and that sentence stands wherever
#     the intervals are described;
#   * every figure they quote is a row of dev/sim-pooled-values.csv, the
#     output of the package's own simulation, and that file is the output of
#     dev/sim-pooled-values.R as it stands;
#   * no other number stands in them but the settings of that simulation and
#     the constants of the code, each read from where it is set;
#   * a text that says how a number tends to behave quotes a figure for it,
#     or sends the reader to the section that does.
#
# The NEWS entries about pooled values are held to the first two.
#
# The word rules read whole help pages, so a sentence added anywhere on a
# page that describes pooled values is read.  The parts of a page are found
# by their titles, labels or first words, and each is held to its opening
# words, so a reader that returns another part fails.
#
# From the source tree (devtools::test()) the pages are read from man/,
# vignettes/ and NEWS.md.  Under R CMD check the tests run beside the
# installed package: the help is read from its Rd database, NEWS and the
# vignettes from the installed files; dev/ is not built, so the figures and
# the numbers are compared in the source tree only.

# ---- The readers --------------------------------------------------------------------

# The root of the source tree, or NULL when the tests do not run from one.
.docs_source <- function() {
  root <- testthat::test_path("..", "..")
  if (file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "man")))
    normalizePath(root)
  else NULL
}

# A help page as parsed Rd: the file in man/, or the installed page.
.docs_rd <- function(page) {
  src <- .docs_source()
  if (!is.null(src)) return(tools::parse_Rd(file.path(src, "man", paste0(page, ".Rd"))))
  tools::Rd_db("spatialkit")[[paste0(page, ".Rd")]]
}

# The text of a piece of parsed Rd, markup dropped.
.docs_flat <- function(x)
  paste(unlist(lapply(x, function(el) if (is.character(el)) el else .docs_flat(el))),
        collapse = "")

# Text on one line: no backticks, runs of white space as one space.
.docs_squash <- function(s) trimws(gsub("[[:space:]]+", " ", gsub("`", "", s, fixed = TRUE)))

.docs_tag <- function(node) attr(node, "Rd_tag") %||% ""

# The \section{}s of a parsed page: their text, named by title.
.docs_rd_sections <- function(rd) {
  tags <- vapply(rd, .docs_tag, character(1))
  secs <- rd[tags == "\\section"]
  stats::setNames(vapply(secs, function(s) .docs_squash(.docs_flat(s[[2L]])), character(1)),
                  vapply(secs, function(s) .docs_squash(.docs_flat(s[[1L]])), character(1)))
}

# The body of one \section{}, as parsed Rd; NULL when the page has none of
# that title.
.docs_rd_section <- function(rd, title) {
  for (node in rd)
    if (identical(.docs_tag(node), "\\section") &&
        identical(.docs_squash(.docs_flat(node[[1L]])), title)) return(node[[2L]])
  NULL
}

# The text of the top-level parts of a page that carry one tag, such as
# "\\description", "\\value" or "\\examples".
.docs_rd_part <- function(rd, tag) {
  parts <- rd[vapply(rd, .docs_tag, character(1)) == tag]
  .docs_squash(paste(vapply(parts, .docs_flat, character(1)), collapse = " "))
}

# Every \item{label}{text} inside a piece of parsed Rd, which are the entries
# of its \describe lists and of \arguments: their text, named by label.
.docs_rd_entries <- function(x) {
  out <- character(0)
  walk <- function(node) {
    if (!is.list(node)) return(invisible())
    if (identical(.docs_tag(node), "\\item") && length(node) == 2L) {
      out[[.docs_squash(.docs_flat(node[[1L]]))]] <<- .docs_squash(.docs_flat(node[[2L]]))
      return(invisible())
    }
    for (el in node) walk(el)
  }
  walk(x)
  out
}
.docs_rd_arguments <- function(rd)
  .docs_rd_entries(rd[vapply(rd, .docs_tag, character(1)) == "\\arguments"])

# The bullets of every \itemize inside a piece of parsed Rd, in order.
.docs_rd_bullets <- function(x) {
  out <- character(0)
  walk <- function(node) {
    if (!is.list(node)) return(invisible())
    if (identical(.docs_tag(node), "\\itemize")) {
      group <- cumsum(vapply(node, .docs_tag, character(1)) == "\\item")
      for (g in setdiff(unique(group), 0L))
        out[[length(out) + 1L]] <<-
          .docs_squash(paste(vapply(node[group == g], .docs_flat, character(1)), collapse = ""))
      return(invisible())
    }
    for (el in node) walk(el)
  }
  walk(x)
  out
}

# The paragraphs of a piece of parsed Rd: its text cut at blank lines.
.docs_rd_paragraphs <- function(x) {
  p <- strsplit(.docs_flat(x), "\n[ \t]*\n", perl = TRUE)[[1L]]
  p <- vapply(p, .docs_squash, character(1), USE.NAMES = FALSE)
  p[nzchar(p)]
}

# The one text of several that starts with the words given; NA when none
# does or more than one.
.docs_starting <- function(texts, first) {
  hit <- texts[startsWith(texts, first)]
  if (length(hit) == 1L) unname(hit) else NA_character_
}

# How often a text holds a phrase.
.docs_count <- function(text, phrase) {
  at <- gregexpr(phrase, text, fixed = TRUE)[[1L]]
  if (at[1L] < 0L) 0L else length(at)
}

# The lines of a vignette's source, or NULL when it is not there to read.
.docs_vignette <- function(name) {
  src  <- .docs_source()
  file <- if (!is.null(src)) file.path(src, "vignettes", paste0(name, ".Rmd"))
          else system.file("doc", paste0(name, ".Rmd"), package = "spatialkit")
  if (!nzchar(file) || !file.exists(file)) return(NULL)
  readLines(file, warn = FALSE, encoding = "UTF-8")
}

# Which markdown lines belong to a fenced code chunk, the fences included.
.docs_md_in_chunk <- function(lines) {
  fence <- grepl("^```", lines)
  (cumsum(fence) %% 2L == 1L) | fence
}

# The sections of markdown lines at one depth of heading ("## " is 2): their
# text, named by title.  A section runs to the next heading of its depth or
# above and holds its code chunks, unless `code` is FALSE.
.docs_md_sections <- function(lines, depth = 2L, code = TRUE) {
  inside <- .docs_md_in_chunk(lines)
  level  <- ifelse(!inside & grepl("^#{2,6} ", lines), nchar(sub(" .*$", "", lines)), NA_integer_)
  at     <- which(level == depth)
  if (!length(at)) return(stats::setNames(character(0), character(0)))
  stops  <- which(level <= depth)
  keep   <- if (code) rep(TRUE, length(lines)) else !inside
  text   <- vapply(at, function(a) {
    nxt <- stops[stops > a]
    to  <- if (length(nxt)) nxt[1L] - 1L else length(lines)
    idx <- seq.int(a + 1L, length.out = to - a)
    .docs_squash(paste(lines[idx][keep[idx]], collapse = " "))
  }, character(1))
  stats::setNames(text, trimws(sub("^#+ ", "", lines[at])))
}

# The front matter of a vignette: the lines between its first two "---".
.docs_md_front <- function(lines) {
  at <- which(lines == "---")
  if (length(at) < 2L) return(NA_character_)
  .docs_squash(paste(lines[seq.int(at[1L] + 1L, length.out = at[2L] - at[1L] - 1L)],
                     collapse = " "))
}

# The lines of NEWS.md, or NULL when the file is not there to read.
.docs_news <- function() {
  src  <- .docs_source()
  file <- if (!is.null(src)) file.path(src, "NEWS.md")
          else system.file("NEWS.md", package = "spatialkit")
  if (!nzchar(file) || !file.exists(file)) return(NULL)
  readLines(file, warn = FALSE, encoding = "UTF-8")
}

# The entries under the first "# " heading of NEWS lines, the development
# version: each "* " item up to the next blank line or heading.
.docs_news_entries <- function(lines) {
  top <- grep("^# ", lines)
  if (!length(top)) return(character(0))
  lines <- lines[seq.int(top[1L], if (length(top) >= 2L) top[2L] - 1L else length(lines))]
  ends  <- which(!nzchar(trimws(lines)) | grepl("^#", lines))
  vapply(grep("^\\* ", lines), function(a) {
    nxt <- ends[ends > a]
    to  <- if (length(nxt)) nxt[1L] - 1L else length(lines)
    .docs_squash(sub("^\\* ", "", paste(lines[a:to], collapse = " ")))
  }, character(1))
}

# ---- The rules ----------------------------------------------------------------------

.DOCS_DENIAL <- "Away from the observations no interval here is calibrated."

# The words of a text that name existing units.
.docs_barred <- function(text)
  unlist(regmatches(text, gregexpr("(?i)\\b(county|counties|zip|zcta)\\w*", text,
                                   perl = TRUE)))

# The places a text says "calibrated" (or any word of that stem) outside the
# sentence that denies it: each with the words round it.
.docs_calibrated <- function(text) {
  rest <- gsub(.DOCS_DENIAL, " ", text, fixed = TRUE)
  at   <- gregexpr("(?i)\\w*calibrat\\w*", rest, perl = TRUE)[[1L]]
  if (at[1L] < 0L) return(character(0))
  substring(rest, pmax(1L, at - 40L), pmin(nchar(rest), at + attr(at, "match.length") + 39L))
}

# The numbers a text states: digits with their thousands commas, decimals
# and exponent.  Digits inside a name (tau2, sigma2), exponents written with
# "^" and a version such as 2.0.0 are not numbers of the text.
.docs_numbers <- function(text)
  regmatches(text, gregexpr(
    "(?<![A-Za-z0-9_.^])[0-9]+(?:,[0-9]{3})*(?:\\.[0-9]+)?(?:e[+-]?[0-9]+)?(?![A-Za-z0-9_]|\\.[0-9])",
    text, perl = TRUE))[[1L]]

# The words by which a text says how a number tends to behave.
.docs_hedges <- function(text)
  unlist(regmatches(text, gregexpr(
    "(?i)\\b(tends?|tended|tendency|usually|typically|generally|often|mostly|in practice)\\b",
    text, perl = TRUE)))

# ---- What is read -------------------------------------------------------------------

# The help pages that describe pooled values.  The word rules are held on
# the whole text of each.
.DOCS_PAGES <- c("summarize_by_cell", "resolution_profile", "select_resolution",
                 "cell_values_at", "get_voronoi_seeds", "voronoi_seeds_kmeans",
                 "build_tessellation", "create_voronoi_polygons",
                 "assign_features_to_polygons", "spatialkit-package",
                 "summary.resolution_profile", "plot.resolution_profile")

# What those pages say with a word of the stem "calibrat" about something
# other than an interval, in its words: Moran's z in ?resolution_profile.
# It is the only exception to the rule, and it stands once.
.DOCS_OTHER_SENSE <- list(resolution_profile = "z is calibrated and flat in L")

.docs_pages <- function() {
  out <- list()
  for (page in .DOCS_PAGES) {
    rd <- .docs_rd(page)
    out[[paste0("?", page)]] <- if (is.null(rd)) NA_character_ else .docs_squash(.docs_flat(rd))
  }
  out
}

# The parts of those pages about pooled values, each found by its title, its
# label or its first words.
.docs_parts <- function() {
  out <- list()
  put <- function(name, text) out[[name]] <<- if (length(text) == 1L) unname(text) else NA_character_
  section <- function(rd, title) {
    body <- .docs_rd_section(rd, title)
    if (is.null(body)) NA_character_ else .docs_squash(.docs_flat(body))
  }
  paragraph <- function(rd, title, first) {
    body <- .docs_rd_section(rd, title)
    if (is.null(body)) NA_character_ else .docs_starting(.docs_rd_paragraphs(body), first)
  }
  entry <- function(rd, title, label) {
    body <- .docs_rd_section(rd, title)
    if (is.null(body)) NA_character_ else .docs_rd_entries(body)[label]
  }

  rd <- .docs_rd("summarize_by_cell")
  put("?summarize_by_cell: Description", .docs_rd_part(rd, "\\description"))
  for (title in c("Pooled values", "Intervals: which one answers which question",
                  "Conditions", "Standard errors and design effects", "Confidence intervals",
                  "How this was checked"))
    put(paste0("?summarize_by_cell: ", title), section(rd, title))
  put("?summarize_by_cell: Confidence intervals, the two intervals",
      paragraph(rd, "Confidence intervals", "At no setting is ..ci_* a range for one observation"))
  args <- .docs_rd_arguments(rd)
  put("?summarize_by_cell: argument estimate", args["estimate"])
  put("?summarize_by_cell: argument obs_level", args["obs_level"])
  put("?summarize_by_cell: Value", .docs_rd_part(rd, "\\value"))

  rd <- .docs_rd("resolution_profile")
  put("?resolution_profile: Observations per cell", section(rd, "Observations per cell"))
  put("?resolution_profile: The ladder and its bounds, min_cell_n",
      paragraph(rd, "The ladder and its bounds", "min_cell_n sets the ceiling of this ladder."))
  put("?resolution_profile: The ladder and its bounds, a ladder of sizes",
      paragraph(rd, "The ladder and its bounds", "A ladder that obs_per_cell sets"))
  put("?resolution_profile: The ladder and its bounds, no variogram",
      paragraph(rd, "The ladder and its bounds", "A ladder of sizes under sac = \"auto\""))
  put("?resolution_profile: The criteria, cp", entry(rd, "The criteria", "cp"))
  put("?resolution_profile: The criteria, loo", entry(rd, "The criteria", "loo"))
  put("?resolution_profile: How this was checked", section(rd, "How this was checked"))

  rd <- .docs_rd("select_resolution")
  for (label in c("cp", "loo", "size", "auto"))
    put(paste0("?select_resolution: What each criterion returns, ", label),
        entry(rd, "What each criterion returns", label))
  put("?select_resolution: A profile built with an estimate",
      section(rd, "A profile built with an estimate"))

  rd <- .docs_rd("cell_values_at")
  put("?cell_values_at: Away from the observations", section(rd, "Away from the observations"))
  put("?cell_values_at: How this was checked", section(rd, "How this was checked"))

  rd <- .docs_rd("spatialkit-package")
  put("?spatialkit: the two choices",
      paragraph(rd, "The pipeline, in order", "Two choices are the user's"))
  body <- .docs_rd_section(rd, "Defaults and their sources")
  bullets <- if (is.null(body)) character(0) else .docs_rd_bullets(body)
  found <- vapply(.DOCS_CONVENTIONS, function(first) .docs_starting(bullets, first),
                  character(1))
  put("?spatialkit: the conventions of pooled values",
      if (anyNA(found)) NA_character_ else paste(found, collapse = " "))

  put("?get_voronoi_seeds: argument partition",
      .docs_rd_arguments(.docs_rd("get_voronoi_seeds"))["partition"])
  put("?voronoi_seeds_kmeans: argument partition",
      .docs_rd_arguments(.docs_rd("voronoi_seeds_kmeans"))["partition"])
  put("?build_tessellation: argument crs",
      .docs_rd_arguments(.docs_rd("build_tessellation"))["crs"])
  out
}

# The six conventions of pooled values among the package's defaults, by
# their first words.
.DOCS_CONVENTIONS <- c("The size of", "The 8 nearest cells", "The 5 residuals",
                       "The 30 values", "The 2 percent", "The two work limits")

# The opening words of each part: what holds a part to its title.
.DOCS_OPENINGS <- c(
  "?summarize_by_cell: Description" =
    "Aggregates an sf point dataset into one row per cell.",
  "?summarize_by_cell: Pooled values" =
    "With estimate = \"pooled\" every numeric column also gets a pooled value.",
  "?summarize_by_cell: Intervals: which one answers which question" =
    "Two families of interval columns can be in the frame.",
  "?summarize_by_cell: Conditions" =
    "What the pooled columns and the count record have to say is raised as R warnings",
  "?summarize_by_cell: Standard errors and design effects" =
    "What the ..se_* columns estimate follows deff.",
  "?summarize_by_cell: Confidence intervals" =
    "With conf_level set, every numeric response and predictor column gains four more",
  "?summarize_by_cell: How this was checked" =
    "The figures on this page that describe pooled values and the two intervals come from one simulation",
  "?summarize_by_cell: Confidence intervals, the two intervals" =
    "At no setting is ..ci_* a range for one observation",
  "?summarize_by_cell: argument estimate" =
    "Which value of a cell the frame reports: \"pooled\" (default) or \"mean\".",
  "?summarize_by_cell: argument obs_level" =
    "Optional level in (0, 1), such as 0.90, of an interval for one further observation",
  "?summarize_by_cell: Value" =
    "A tibble/data.frame (or sf if cells_sf given) with per-cell summaries",
  "?resolution_profile: Observations per cell" =
    "obs_per_cell sets the ladder by a size: the mean number of rows a cell is to hold.",
  "?resolution_profile: The ladder and its bounds, min_cell_n" =
    "min_cell_n sets the ceiling of this ladder.",
  "?resolution_profile: The ladder and its bounds, a ladder of sizes" =
    "A ladder that obs_per_cell sets is a ladder of sizes",
  "?resolution_profile: The ladder and its bounds, no variogram" =
    "A ladder of sizes under sac = \"auto\" fits no variogram",
  "?resolution_profile: The criteria, cp" =
    "Mallows' C_p of the piecewise-constant approximation of the response",
  "?resolution_profile: The criteria, loo" =
    "With an estimate: the leave-one-out mean squared error of that cell value",
  "?resolution_profile: How this was checked" =
    "The sizes for pooled values were compared in a simulation of 200 layers",
  "?select_resolution: What each criterion returns, cp" =
    "The level that minimises Mallows' C_p: the count at which re-measuring the sampled sites",
  "?select_resolution: What each criterion returns, loo" =
    "The level with the lowest leave-one-out error of the estimate the profile was built for",
  "?select_resolution: What each criterion returns, size" =
    "The count a size gives, max(2, round(n / obs_per_cell))",
  "?select_resolution: What each criterion returns, auto" =
    "The rule that goes with the profile's estimate",
  "?select_resolution: A profile built with an estimate" =
    "This function's default criterion is \"auto\": \"size\" for a profile of pooled values,",
  "?cell_values_at: Away from the observations" =
    "The value and the interval describe locations among the observations.",
  "?cell_values_at: How this was checked" =
    "A simulation shows how the interval's coverage falls with dist_ratio",
  "?spatialkit: the two choices" =
    "Two choices are the user's, and each has an argument with a default.",
  "?spatialkit: the conventions of pooled values" =
    "The size of 4.5 observations per cell",
  "?get_voronoi_seeds: argument partition" =
    "For method = \"kmeans\": what the seeds minimise.",
  "?voronoi_seeds_kmeans: argument partition" =
    "What the seeds minimise.",
  "?build_tessellation: argument crs" =
    "Optional target CRS: anything sf::st_crs() accepts")

# Where the intervals are described, and the denial therefore has to stand.
.DOCS_DENIED_IN <- c("?summarize_by_cell: Intervals: which one answers which question",
                     "?summarize_by_cell: Confidence intervals, the two intervals",
                     "?summarize_by_cell: Value",
                     "?cell_values_at: Away from the observations",
                     "vignette cells: What summarize_by_cell() returns",
                     "vignette reporting: Reporting a pooled cell value")

# The entries of NEWS's first section about the pooled values of
# summarize_by_cell(), found by those two names; NULL when NEWS is not there.
.docs_news_pooled <- function() {
  lines <- .docs_news()
  if (is.null(lines)) return(NULL)
  entries <- .docs_news_entries(lines)
  entries[grepl("summarize_by_cell()", entries, fixed = TRUE) &
            grepl("pooled value", entries, ignore.case = TRUE)]
}

# The sections of the vignettes that describe pooled values, and the one
# older section a sentence about cp stands in; with them the front matter.
# `code = FALSE` leaves the code chunks out, and the openings below are
# those of the prose.
.DOCS_VIGNETTES <- list(
  resolution = list("2" = c("A size for pooled values", "Leave-one-out for plain means"),
                    "3" = "Four criteria and a diagnostic"),
  reporting  = list("2" = "Reporting a pooled cell value"),
  cells      = list("2" = c("What summarize_by_cell() returns",
                            "Pooled values and their warnings",
                            "Values at other locations")))
.docs_vignette_texts <- function(code = TRUE) {
  out <- list()
  for (name in names(.DOCS_VIGNETTES)) {
    lines <- .docs_vignette(name)
    if (is.null(lines)) return(NULL)
    out[[paste0("vignette ", name, ": front matter")]] <- .docs_md_front(lines)
    for (depth in names(.DOCS_VIGNETTES[[name]])) {
      secs <- .docs_md_sections(lines, depth = as.integer(depth), code = code)
      for (title in .DOCS_VIGNETTES[[name]][[depth]])
        out[[paste0("vignette ", name, ": ", title)]] <- unname(secs[title])
    }
  }
  out
}
.DOCS_VIGNETTE_OPENINGS <- c(
  "vignette resolution: front matter"  = "title: \"Choosing a resolution\"",
  "vignette resolution: A size for pooled values" = paste(
    "By default summarize_by_cell() reports, beside each cell's mean, a pooled value:",
    "the cell's mean pulled toward the mean of the observations in the 8 nearest cells"),
  "vignette resolution: Leave-one-out for plain means" =
    "With plain means the leave-one-out error does choose a count.",
  "vignette resolution: Four criteria and a diagnostic" =
    "elbow sees the coordinates and knows nothing about what you measured.",
  "vignette reporting: front matter"   = "title: \"Exporting and reporting results\"",
  "vignette reporting: Reporting a pooled cell value" =
    "summarize_by_cell() returns by default, beside each cell's mean, a pooled value",
  "vignette cells: front matter"       = "title: \"Building cells and summarising them\"",
  "vignette cells: What summarize_by_cell() returns" =
    "summarize_by_cell() reduces the assigned rows to one row per cell.",
  "vignette cells: Pooled values and their warnings" =
    "summarize_by_cell() reports problems with the pooled values",
  "vignette cells: Values at other locations" =
    "cell_values_at() gives a location the value of the cell it lies in")

# ---- The simulation's output ---------------------------------------------------------

# The MD5 sum of a file read as lines of text, as dev/sim-pooled-values.R
# takes it: the same whatever the line endings of a checkout.
.docs_text_md5 <- function(file) {
  tmp <- withr::local_tempfile()
  con <- file(tmp, "wb")
  writeLines(readLines(file, warn = FALSE), con, sep = "\n", useBytes = TRUE)
  close(con)
  unname(tools::md5sum(tmp))
}

# The simulation's files, read from the source tree: the script's lines and
# its settings (the constants it assigns between "# ---- Settings" and
# OUT_DEFAULT), the table of its outputs, and the record of the run.  NULL
# beside an installed package, where dev/ is not.
.docs_sim <- function() {
  src <- .docs_source()
  if (is.null(src)) return(NULL)
  files <- file.path(src, "dev", c("sim-pooled-values.R", "sim-pooled-values.csv",
                                   "sim-pooled-values-run.dcf"))
  if (!all(file.exists(files))) return(NULL)
  script <- readLines(files[1L], warn = FALSE)
  from   <- grep("^# ---- Settings", script)
  to     <- grep("^OUT_DEFAULT <- ", script)
  block  <- script[seq.int(from[1L] + 1L, to[1L] - 1L)]
  set    <- new.env(parent = baseenv())
  eval(parse(text = block[grepl("^[A-Z_]+ +<- ", block)]), set)
  list(files = files, set = as.list(set),
       out = utils::read.csv(files[2L], stringsAsFactors = FALSE),
       run = as.list(read.dcf(files[3L])[1L, ]))
}

# One row of the output, by its output, where and class (and, where those
# three leave several, a part of its quantity).
.docs_row <- function(out, output, where, class, quantity = "") {
  hit <- out[out$output == output & out$where == where & out$class == class &
               grepl(quantity, out$quantity, fixed = TRUE), ]
  if (nrow(hit) != 1L)
    stop(sprintf("%d rows of the output for: %s, %s, %s, %s", nrow(hit), output, where, class,
                 quantity))
  hit
}

# A count of layers as the pages write it.
.docs_layers <- function(k) if (k == 0L) "none" else if (k == 1L) "one" else as.character(k)

# Every statement of the run the pages make, as (the text it stands in, its
# words), built from the rows of the output.
.docs_claims <- function(sim) {
  out <- sim$out
  set <- sim$set
  row     <- function(...) .docs_row(out, ...)
  ratio   <- function(...) sprintf("%.2f", row(...)$value)
  percent <- function(...) sprintf("%.1f", 100 * row(...)$value)
  layers  <- function(...) { r <- row(...); as.integer(round(r$value * r$n)) }

  rows    <- c("1 row", "2 to 4 rows", "5 to 9 rows", "10 rows or more")
  a_cl    <- vapply(paste(rows, "of the layer in the cell"),
                    function(cl) ratio("a", "among", cl), character(1))
  a_all   <- ratio("a", "among", "all rows")
  a_block <- ratio("a", "block", "all rows")
  b_among <- percent("b", "among", "all rows")
  b_block <- percent("b", "block", "all rows")
  b_lo    <- sprintf("%.1f", 100 * row("b", "among", "all rows")$layer_min)
  b_hi    <- sprintf("%.1f", 100 * row("b", "among", "all rows")$layer_max)
  b_ratio <- vapply(c("below 1", "1 to 2", "2 to 4", "4 to 8", "above 8"), function(cl)
    percent("b", "among and block", paste("dist_ratio", cl)), character(1))
  c_among <- percent("c", "among", "rows whose cell has the interval")
  c_block <- percent("c", "block", "rows whose cell has the interval")
  d_w     <- ratio("d", "white noise", "layers", "median own weight")
  d_tau0  <- layers("d", "white noise", "layers", "tau2 = 0")
  d_n     <- row("d", "white noise", "layers", "tau2 = 0")$n
  d_cor   <- ratio("d", "white noise", "layers", "pooled value with the mean")
  d_corm  <- ratio("d", "white noise", "layers", "plain mean with the mean")
  e_among <- ratio("e", "among", "cells with at least 10 further locations")
  e_block <- sprintf("%.1f", row("e", "block", "cells with at least 10 further locations")$value)
  f_at    <- function(what, size)
    layers("f", what, sprintf("obs_per_cell = %s", size), "at this size")
  f_n     <- row("f", "pooled value", "obs_per_cell = 2", "at this size")$n
  finest  <- min(set$SIZES)

  # What the pages say in words has to hold of the rows.
  stopifnot(
    "the summary raises the warning on every white-noise layer" =
      identical(layers("d", "white noise", "layers", "summarize_by_cell()"), d_n),
    "the plain mean's lowest loo is at 18 and at 3 on the same number of layers" =
      identical(f_at("plain mean", 18), f_at("plain mean", 3)),
    "the finest size is the size 2" = finest == 2,
    "the two rows of the finest size agree" =
      identical(f_at("pooled value", finest),
                layers("f", "pooled value", sprintf("obs_per_cell = %s", finest),
                       "finest size")),
    "the pooled value's lowest loo is at 2, 3 or 4.5 on every layer" =
      f_at("pooled value", 2) + f_at("pooled value", 3) + f_at("pooled value", 4.5) == f_n,
    "the plain mean's lowest loo is at one of the six sizes named on every layer" =
      f_at("plain mean", 6) + f_at("plain mean", 4.5) + f_at("plain mean", 9) +
        f_at("plain mean", 12) + f_at("plain mean", 18) + f_at("plain mean", 3) +
        f_at("plain mean", 2) == f_n)

  hc  <- "?summarize_by_cell: How this was checked"
  op  <- "?resolution_profile: How this was checked"
  va  <- "?cell_values_at: How this was checked"
  vs  <- "vignette resolution: A size for pooled values"
  vl  <- "vignette resolution: Leave-one-out for plain means"
  vr  <- "vignette reporting: Reporting a pooled cell value"
  vw  <- "vignette cells: Pooled values and their warnings"
  claim <- function(where, ...) list(where = where, phrase = sprintf(...))
  list(
    claim(hc, "was %s of the plain mean's", a_all),
    claim(hc, paste("%s in cells that hold one of the layer's observations, %s in",
                    "cells of 2 to 4, %s in cells of 5 to 9 and %s in cells of 10 or more"),
          a_cl[1], a_cl[2], a_cl[3], a_cl[4]),
    claim(hc, "In the held-out block it was %s.", a_block),
    claim(hc, "the median own weight was %s", d_w),
    claim(hc, "tau2 was 0 on %d layers", d_tau0),
    claim(hc, "the warning of class \"spatialkit_no_structure\" was raised on all %d", d_n),
    claim(hc, "a correlation of %s with", d_cor),
    claim(hc, "where the plain means had one of %s", d_corm),
    claim(hc, "was %s of (1 - w) * tau2 among the observations and %s times it",
          e_among, e_block),
    claim(hc, paste("held %s percent of the observations held out at random among the",
                    "others (%s to %s percent in single layers) and %s percent of those in",
                    "the held-out block"), b_among, b_lo, b_hi, b_block),
    claim(hc, "held %s and %s percent of them", c_among, c_block),
    claim(op, "was at the size 2 on %d layers, at 3 on %d and at 4.5 on %d",
          f_at("pooled value", 2), f_at("pooled value", 3), f_at("pooled value", 4.5)),
    claim(op, paste("That of the plain mean was at 2 on %s: it was at 6 on %d layers, at",
                    "4.5 on %d, at 9 on %d, at 12 on %d, and at 18 and at 3 on %s each."),
          .docs_layers(f_at("plain mean", 2)), f_at("plain mean", 6), f_at("plain mean", 4.5),
          f_at("plain mean", 9), f_at("plain mean", 12), .docs_layers(f_at("plain mean", 18))),
    claim(va, paste("held %s percent of those with a dist_ratio below 1, %s percent from 1",
                    "to 2, %s percent from 2 to 4, %s percent from 4 to 8 and %s percent",
                    "above 8"), b_ratio[1], b_ratio[2], b_ratio[3], b_ratio[4], b_ratio[5]),
    claim(va, paste("It held %s percent of the observations held out at random and %s",
                    "percent of those in the block"), b_among, b_block),
    claim(vs, "the lowest loo of the pooled value was at the size 2 on %d of the %d layers (",
          f_at("pooled value", 2), f_n),
    claim(vl, "the lowest loo of the plain mean was at the finest size on %s of the %d layers.",
          .docs_layers(f_at("plain mean", finest)), f_n),
    claim(vl, "It was at 6 observations per cell on %d of them, at 4.5 on %d and at 9 on %d.",
          f_at("plain mean", 6), f_at("plain mean", 4.5), f_at("plain mean", 9)),
    claim(vr, paste("held %s percent of held-out observations with a dist_ratio below 1",
                    "and %s percent of those above 8"), b_ratio[1], b_ratio[5]),
    claim(vw, "the median own weight was %s", d_w),
    claim(vw, "a correlation of %s with", d_cor),
    claim(vw, "where the plain means had one of %s", d_corm),
    claim(vw, paste("in the simulation it was raised on all %d white-noise layers, and the",
                    "reason \"no signal between cells\" on %d of them"), d_n, d_tau0))
}

# Every text a figure can stand in: the parts of the help and the prose of
# the vignette sections.
.docs_all_texts <- function() c(.docs_parts(), .docs_vignette_texts(FALSE))

# A text with the run's statements taken out: what is left holds the
# numbers that are not figures of the simulation.
.docs_rest <- function(text, claims, where) {
  for (cl in claims) if (identical(cl$where, where)) text <- sub(cl$phrase, " ", text, fixed = TRUE)
  text
}

# The texts written for pooled values whole: the ones whose numbers are all
# accounted for, and in which a statement of how a number tends to behave
# needs a figure.  The older sections that hold a sentence or two about
# pooled values are read by the word rules only.
.DOCS_WRITTEN <- c(
  "?summarize_by_cell: Pooled values",
  "?summarize_by_cell: Intervals: which one answers which question",
  "?summarize_by_cell: Conditions",
  "?summarize_by_cell: How this was checked",
  "?summarize_by_cell: Confidence intervals, the two intervals",
  "?summarize_by_cell: argument estimate",
  "?summarize_by_cell: argument obs_level",
  "?resolution_profile: Observations per cell",
  "?resolution_profile: The ladder and its bounds, min_cell_n",
  "?resolution_profile: The ladder and its bounds, a ladder of sizes",
  "?resolution_profile: The ladder and its bounds, no variogram",
  "?resolution_profile: The criteria, loo",
  "?resolution_profile: How this was checked",
  "?select_resolution: What each criterion returns, loo",
  "?select_resolution: What each criterion returns, size",
  "?select_resolution: What each criterion returns, auto",
  "?select_resolution: A profile built with an estimate",
  "?cell_values_at: Away from the observations",
  "?cell_values_at: How this was checked",
  "?spatialkit: the two choices",
  "?spatialkit: the conventions of pooled values",
  "vignette resolution: A size for pooled values",
  "vignette resolution: Leave-one-out for plain means",
  "vignette reporting: Reporting a pooled cell value",
  "vignette cells: What summarize_by_cell() returns",
  "vignette cells: Pooled values and their warnings",
  "vignette cells: Values at other locations")

# The texts that state how a number tends to behave and send the reader to
# the section that quotes the figure for it, each with that section.
.DOCS_SEE <- c(
  "?resolution_profile: Observations per cell"           = "?resolution_profile: How this was checked",
  "?resolution_profile: The criteria, loo"               = "?resolution_profile: How this was checked",
  "?select_resolution: What each criterion returns, loo" = "?resolution_profile: How this was checked")

# ---- The readers and the rules themselves --------------------------------------------

test_that("a section is found by its title and read without its markup", {
  rd_file <- withr::local_tempfile(fileext = ".Rd")
  writeLines(c(
    "\\name{x}", "\\alias{x}", "\\title{A page}",
    "\\description{Counties are named in the \\code{description}.}",
    "\\section{Intervals: which one}{",
    "First \\code{a + b} line,", "  and a second with \\emph{stress}.",
    "\\itemize{", "\\item one", "\\item two \\verb{..ci_*}", "}", "}",
    "\\section{Other}{Nothing here.}"), rd_file)
  rd   <- tools::parse_Rd(rd_file)
  secs <- .docs_rd_sections(rd)
  expect_identical(names(secs), c("Intervals: which one", "Other"))
  expect_identical(secs[["Intervals: which one"]],
                   "First a + b line, and a second with stress. one two ..ci_*")
  expect_identical(secs[["Other"]], "Nothing here.")
  # A section's text holds no other part of the page; the page's text holds all.
  expect_false(any(grepl("Counties", secs, fixed = TRUE)))
  expect_match(.docs_squash(.docs_flat(rd)), "Counties are named in the description.",
               fixed = TRUE)
  expect_identical(.docs_rd_part(rd, "\\description"), "Counties are named in the description.")
  expect_identical(.docs_squash(.docs_flat(.docs_rd_section(rd, "Other"))), "Nothing here.")
  expect_null(.docs_rd_section(rd, "Intervals"))
  expect_identical(.docs_rd_bullets(rd), c("one", "two ..ci_*"))

  md <- c("---", "title: t", "---", "", "Intro with a county.", "", "## First  ", "",
          "Text of `one`,", "wrapped.", "", "```{r}", "x <- 1  # zip", "## not a heading",
          "```", "", "### A subsection", "stays inside.", "", "## 6. Second", "Last line.")
  ms <- .docs_md_sections(md)
  expect_identical(names(ms), c("First", "6. Second"))
  expect_identical(
    ms[["First"]],
    "Text of one, wrapped. {r} x <- 1 # zip ## not a heading ### A subsection stays inside.")
  expect_identical(ms[["6. Second"]], "Last line.")
  expect_identical(.docs_md_sections(md, code = FALSE)[["First"]],
                   "Text of one, wrapped. ### A subsection stays inside.")
  expect_identical(.docs_md_sections(md, depth = 3L),
                   c("A subsection" = "stays inside."))
  expect_identical(.docs_md_sections(c("no", "heading")),
                   stats::setNames(character(0), character(0)))
  expect_identical(.docs_md_front(md), "title: t")
  expect_identical(.docs_md_in_chunk(md)[12:15], rep(TRUE, 4L))
  expect_false(any(.docs_md_in_chunk(md)[-(12:15)]))
})

test_that("an entry, a bullet, a paragraph and a NEWS entry are found and read without their markup", {
  rd_file <- withr::local_tempfile(fileext = ".Rd")
  writeLines(c(
    "\\name{x}", "\\alias{x}", "\\title{A page}",
    "\\arguments{", "\\item{a}{The first, \\code{NULL} or a number.}", "",
    "\\item{b}{The second.}", "}",
    "\\section{The criteria}{", "Opening words of the first paragraph,", "on two lines.", "",
    "\\describe{", "\\item{\\code{cp}}{An older entry, calibrated.}",
    "\\item{\\code{loo}}{The new \\emph{entry}.}", "}", "",
    "A last paragraph", "of the section.", "}",
    "\\section{Defaults}{", "\\itemize{", "\\item The 8 nearest cells: a convention.",
    "\\item The two limits,", "  on two lines.", "}", "}"), rd_file)
  rd <- tools::parse_Rd(rd_file)
  expect_identical(.docs_rd_arguments(rd),
                   c(a = "The first, NULL or a number.", b = "The second."))
  body <- .docs_rd_section(rd, "The criteria")
  expect_identical(.docs_rd_entries(body),
                   c(cp = "An older entry, calibrated.", loo = "The new entry."))
  expect_identical(.docs_rd_paragraphs(body)[c(1L, 3L)],
                   c("Opening words of the first paragraph, on two lines.",
                     "A last paragraph of the section."))
  expect_identical(.docs_rd_bullets(.docs_rd_section(rd, "Defaults")),
                   c("The 8 nearest cells: a convention.", "The two limits, on two lines."))
  # One text of several by its first words; none or two is no answer.
  texts <- c("The two limits, on two lines.", "The 8 nearest cells.", "The 8 values.")
  expect_identical(.docs_starting(texts, "The two"), "The two limits, on two lines.")
  expect_identical(.docs_starting(texts, "The 8"), NA_character_)
  expect_identical(.docs_starting(texts, "The 5"), NA_character_)
  expect_identical(.docs_count("a figure, a figure and a number", "a figure"), 2L)
  expect_identical(.docs_count("a figure", "another"), 0L)

  news <- c("# pkg (development version)", "", "## New features", "",
            "* First `entry` of", "  two lines.", "", "* Second entry.", "  More.", "",
            "  A second paragraph.", "", "# pkg 1.0", "", "* An older entry.")
  expect_identical(.docs_news_entries(news),
                   c("First entry of two lines.", "Second entry. More."))
  expect_identical(.docs_news_entries(c("no", "heading")), character(0))
})

test_that("the rules find what they are for, and nothing else", {
  expect_identical(.docs_barred("Tracts nest in a County; counties, ZIP codes and ZCTAs do not."),
                   c("County", "counties", "ZIP", "ZCTAs"))
  expect_identical(.docs_barred("a Zip code; zip code areas; zcta5"),
                   c("Zip", "zip", "zcta5"))
  # Letters inside another word are not the word.
  expect_identical(.docs_barred("an unzipped archive and a viscounty"), character(0))
  expect_identical(.docs_barred("cells, regions and observations"), character(0))

  expect_identical(.docs_calibrated(.DOCS_DENIAL), character(0))
  expect_identical(
    .docs_calibrated(paste("An interval.", .DOCS_DENIAL, "A value.", .DOCS_DENIAL)),
    character(0))
  expect_length(.docs_calibrated("The interval is calibrated among the observations."), 1L)
  expect_length(.docs_calibrated(paste(.DOCS_DENIAL, "Its calibration was checked.")), 1L)
  expect_length(.docs_calibrated("Miscalibrated here. Uncalibrated there, far apart from the first one by many words."), 2L)
  # The denial with another word in it is not the denial.
  expect_length(.docs_calibrated("Away from the observations every interval here is calibrated."),
                1L)
  expect_identical(.docs_calibrated("No such word."), character(0))

  expect_identical(
    .docs_numbers("tau2 = max(0, x) / (1 - w); 1,500 rows at 4.5 a cell, 0.90, 2e9 and 1e10; at 18, 12 and 9."),
    c("0", "1", "1,500", "4.5", "0.90", "2e9", "1e10", "18", "12", "9"))
  expect_identical(.docs_numbers("held 89.6 percent; -0.02; the 8-nearest; about 44,700."),
                   c("89.6", "0.02", "8", "44,700"))
  # Digits inside a name, an exponent and a version are not numbers of the text.
  expect_identical(.docs_numbers("sigma2, tau2, s^2, loo_se, ..ci_lo_* and k-means++"),
                   character(0))
  expect_identical(.docs_numbers("changed from 2.0.0 (and 10.2.1), at 1.5."), "1.5")

  expect_identical(.docs_hedges("It tends to fall; in practice often, and usually so. Mostly."),
                   c("tends", "in practice", "often", "usually", "Mostly"))
  expect_identical(.docs_hedges("It extends the pretend, softened, usual cause."), character(0))
})


# ---- The help -------------------------------------------------------------------------

test_that("the parts of the help written for pooled cell values are where their titles, labels and first words say", {
  expect_length(.DOCS_PAGES, 12L)
  for (page in .DOCS_PAGES) expect_false(is.null(.docs_rd(page)), info = page)
  parts <- .docs_parts()
  expect_length(parts, 30L)
  expect_identical(names(parts), names(.DOCS_OPENINGS))
  for (what in names(parts)) {
    expect_false(is.na(parts[[what]]), info = what)
    # A part read under a title is that part: it opens with its own words.
    expect_identical(substr(parts[[what]], 1L, nchar(.DOCS_OPENINGS[[what]])),
                     .DOCS_OPENINGS[[what]], info = what)
    # A part that is there has text to read beyond them: the rules cannot
    # pass on an empty one.
    expect_gt(nchar(parts[[what]]), nchar(.DOCS_OPENINGS[[what]]) + 100L, label = what)
  }
  # No two titles are answered with one text.
  expect_identical(anyDuplicated(unlist(parts)), 0L)
  # The lead sentence, and the six conventions.
  expect_identical(.docs_count(parts[["?summarize_by_cell: Description"]], paste(
    "By default each numeric column also gets a pooled value, the cell's mean pulled toward",
    "the mean of the cells around it (estimate = \"pooled\"; estimate = \"mean\" leaves it",
    "out).")), 1L)
  expect_length(.DOCS_CONVENTIONS, 6L)
  for (first in .DOCS_CONVENTIONS)
    expect_identical(.docs_count(parts[["?spatialkit: the conventions of pooled values"]], first),
                     1L, info = first)
  # Every text held to its numbers is one of the parts or a vignette section.
  expect_true(all(.DOCS_WRITTEN %in% c(names(parts), names(.DOCS_VIGNETTE_OPENINGS))))
})

test_that("no help page that describes pooled values names existing units", {
  pages <- .docs_pages()
  expect_identical(names(pages), paste0("?", .DOCS_PAGES))
  for (what in names(pages)) {
    expect_false(is.na(pages[[what]]), info = what)
    expect_gt(nchar(pages[[what]]), 2000L, label = what)
    expect_identical(.docs_barred(pages[[what]]), character(0), info = what)
  }
  # A part is a piece of its page, so the rule reads every part too.  The
  # conventions are bullets of a longer list, read together: each is a piece.
  parts <- .docs_parts()
  conv  <- "?spatialkit: the conventions of pooled values"
  for (what in setdiff(names(parts), conv)) {
    page <- pages[[sub("spatialkit$", "spatialkit-package", sub(":.*$", "", what))]]
    expect_identical(.docs_count(page, substr(parts[[what]], 1L, 200L)), 1L, info = what)
  }
  for (first in .DOCS_CONVENTIONS)
    expect_identical(.docs_count(pages[["?spatialkit-package"]], first), 1L, info = first)
})

test_that("on those pages \"calibrated\" stands only in the sentence that denies it away from the observations, and once of Moran's z", {
  pages <- .docs_pages()
  expect_identical(names(.DOCS_OTHER_SENSE), "resolution_profile")
  expect_identical(lengths(.DOCS_OTHER_SENSE, use.names = FALSE), 1L)
  for (what in names(pages)) {
    text <- pages[[what]]
    for (phrase in .DOCS_OTHER_SENSE[[sub("^\\?", "", what)]]) {
      expect_identical(.docs_count(text, phrase), 1L, info = paste(what, phrase))
      text <- sub(phrase, " ", text, fixed = TRUE)
    }
    expect_identical(.docs_calibrated(text), character(0), info = what)
  }
})

test_that("the help says, wherever it describes the intervals, that none is calibrated away from the observations", {
  expect_length(.DOCS_DENIED_IN, 6L)
  parts <- .docs_parts()
  for (what in grep("^\\?", .DOCS_DENIED_IN, value = TRUE)) {
    expect_true(what %in% names(parts), info = what)
    expect_identical(.docs_count(parts[[what]], .DOCS_DENIAL), 1L, info = what)
  }
  expect_length(grep("^\\?", .DOCS_DENIED_IN), 4L)
})

test_that("the help says of a level above the limit of Moran's z what the code does there", {
  op <- .docs_parts()[["?resolution_profile: Observations per cell"]]
  expect_identical(.docs_count(op, sprintf(
    "moran_i and moran_z need FNN and Matrix at a level where more than %s cells hold a scored row.",
    format(.KNN_DENSE_MAX, big.mark = ","))), 1L)
  expect_identical(.docs_count(op, paste(
    "Without them a ladder of sizes has NA in both columns there, with a warning of class",
    "spatialkit_moran_skipped")), 1L)
  # Until the rule of a ladder of sizes holds on every ladder, a ladder set
  # another way stops at such a level, and the help says so.
  expect_identical(.docs_count(op, "a ladder set any other way stops there with an error"),
                   if (.moran_skip_any_ladder()) 0L else 1L)
})

test_that("the help and the vignette say that a size set the count of a selection by size, not the caller", {
  # On the default path no one names a size, and select_resolution()'s note
  # calls it the default; the count stands because a size set it and no
  # criterion, as the warning spatialkit_no_structure says.
  said <- "The count stands either way, since a size set it and no criterion"
  size <- .docs_parts()[["?select_resolution: What each criterion returns, size"]]
  expect_identical(.docs_count(size, said), 1L)
  expect_identical(.docs_count(size, "since the caller set it"), 0L)
  prose <- .docs_vignette_texts(code = FALSE)
  if (is.null(prose)) skip("the vignettes' sources are not installed")
  setting <- prose[["vignette resolution: A size for pooled values"]]
  expect_identical(.docs_count(setting, said), 1L)
  expect_identical(.docs_count(setting, "since the caller set it"), 0L)
})

test_that("the help says, in \"Pooled values\" and under spatialkit_pool_few_cells, that a value can cross the overall mean and two cells' values change places", {
  # With 2 to 9 cells the arithmetic stays, and the help and the warning say
  # what follows from it.
  parts <- .docs_parts()
  pv <- parts[["?summarize_by_cell: Pooled values"]]
  expect_identical(.docs_count(pv, paste(
    "A cell whose own weight is below its share of the rows, n / N, then gets a pooled value",
    "on the far side of the overall mean from its plain mean")), 1L)
  expect_identical(.docs_count(pv,
    "with two cells that is the other cell's mean, so the two values change places"), 1L)
  cond <- parts[["?summarize_by_cell: Conditions"]]
  expect_identical(.docs_count(cond, paste(
    "\"spatialkit_pool_few_cells\": 2 to 9 cells, so each target is the mean of all the other",
    "cells, and a cell whose own weight is below its share of the column's rows gets a pooled",
    "value on the far side of the overall mean from its plain mean (cells, and own_weight, the",
    "smallest, by column).")), 1L)
  expect_identical(.docs_count(cond, paste(
    "The text also names each column with two cells and tau2 at 0, which has cells 2 and",
    "own_weight 0: each cell's pooled value is then the other cell's mean, so the two values",
    "change places.")), 1L)
})


# ---- NEWS -----------------------------------------------------------------------------

test_that("the NEWS entries that name pooled values name no existing units and claim no calibration", {
  news <- .docs_news_pooled()
  if (is.null(news)) skip("NEWS.md is not there to read")
  # The development version describes pooled values, so the rule has
  # entries to read.
  expect_gt(length(news), 0L)
  for (entry in news) {
    expect_identical(.docs_barred(entry), character(0), info = substr(entry, 1L, 80L))
    expect_identical(.docs_calibrated(entry), character(0), info = substr(entry, 1L, 80L))
  }
})


# ---- The vignettes ---------------------------------------------------------------------

test_that("the vignette sections written for pooled cell values are where their titles say", {
  prose <- .docs_vignette_texts(code = FALSE)
  if (is.null(prose)) skip("the vignettes' sources are not installed")
  whole <- .docs_vignette_texts()
  expect_length(prose, 10L)
  expect_identical(names(prose), names(.DOCS_VIGNETTE_OPENINGS))
  expect_identical(names(whole), names(prose))
  for (what in names(prose)) {
    expect_false(is.na(prose[[what]]), info = what)
    expect_identical(substr(prose[[what]], 1L, nchar(.DOCS_VIGNETTE_OPENINGS[[what]])),
                     .DOCS_VIGNETTE_OPENINGS[[what]], info = what)
    expect_gt(nchar(prose[[what]]), nchar(.DOCS_VIGNETTE_OPENINGS[[what]]) + 100L, label = what)
    # The text with its code holds no less than the prose.
    expect_gte(nchar(whole[[what]]), nchar(prose[[what]]), label = what)
  }
  expect_identical(anyDuplicated(unlist(prose)), 0L)
  # The sentence about cp in the older section is in the text read.
  expect_identical(.docs_count(prose[["vignette resolution: Four criteria and a diagnostic"]],
                               "What cp estimates is the error of re-measuring the sampled sites"),
                   1L)
})

test_that("no vignette section written for pooled cell values names existing units", {
  texts <- .docs_vignette_texts()
  if (is.null(texts)) skip("the vignettes' sources are not installed")
  for (what in names(texts))
    expect_identical(.docs_barred(texts[[what]]), character(0), info = what)
})

test_that("in those sections \"calibrated\" stands only in the sentence that denies it, and the sentence stands where an interval is reported", {
  texts <- .docs_vignette_texts()
  if (is.null(texts)) skip("the vignettes' sources are not installed")
  for (what in names(texts))
    expect_identical(.docs_calibrated(texts[[what]]), character(0), info = what)
  for (what in grep("^vignette ", .DOCS_DENIED_IN, value = TRUE)) {
    expect_true(what %in% names(texts), info = what)
    expect_identical(.docs_count(texts[[what]], .DOCS_DENIAL), 1L, info = what)
  }
  expect_length(grep("^vignette ", .DOCS_DENIED_IN), 2L)
})


# ---- The figures ------------------------------------------------------------------------

test_that("the simulation's output is the output of the script as it stands", {
  if (is.null(.docs_source())) skip("dev/ is in the source tree only")
  sim <- .docs_sim()
  expect_false(is.null(sim), label = "the script, its output and the record of its run are in dev/")
  if (is.null(sim)) return(invisible())
  # The record of the run names the script and the table by the sums of
  # their text: an edit of either fails here until the script is run again.
  expect_identical(sim$run$script, "dev/sim-pooled-values.R")
  expect_identical(sim$run$script_md5, .docs_text_md5(sim$files[1L]))
  expect_identical(sim$run$csv_md5, .docs_text_md5(sim$files[2L]))
  expect_identical(as.integer(sim$run$layers), sim$set$LAYERS)

  out <- sim$out
  set <- sim$set
  expect_identical(names(out), c("output", "quantity", "where", "class", "value", "n",
                                 "layer_min", "layer_max"))
  expect_identical(unique(out$output), c("a", "b", "c", "d", "e", "f"))
  expect_false(anyNA(out$value))
  # How many rows and layers each figure is of.
  held <- c(among = set$N_AMONG, block = set$N_BLOCK) * set$LAYERS
  for (w in names(held)) {
    expect_identical(.docs_row(out, "a", w, "all rows")$n, held[[w]], info = w)
    expect_identical(.docs_row(out, "b", w, "all rows")$n, held[[w]], info = w)
    expect_identical(sum(out$n[out$output == "a" & out$where == w & out$class != "all rows"]),
                     held[[w]], info = w)
    expect_lte(.docs_row(out, "c", w, "rows whose cell has the interval")$n, held[[w]])
  }
  expect_identical(sum(out$n[out$output == "b" & out$where == "among and block"]), sum(held))
  expect_true(all(out$n[out$output %in% c("d", "f")] == set$LAYERS))
  # A coverage and a share lie between 0 and 1.
  share <- out$output %in% c("b", "c", "f") | (out$output == "d" & grepl("^share", out$quantity))
  expect_identical(sum(share), 28L)
  expect_true(all(out$value[share] >= 0 & out$value[share] <= 1))
  # Output f has one row per size for each estimate, of whole layers, and
  # every layer has its lowest loo at one of the sizes.
  for (what in c("pooled value", "plain mean")) {
    at <- out[out$output == "f" & out$where == what &
                grepl("at this size", out$quantity, fixed = TRUE), ]
    expect_identical(at$class, sprintf("obs_per_cell = %s", set$SIZES), info = what)
    expect_equal(at$value * at$n, round(at$value * at$n), info = what)
    expect_equal(sum(at$value), 1, info = what)
  }
  # The script profiles at the sizes and the default size the package has.
  expect_identical(set$SIZES, .OBS_PER_CELL_LADDER)
  expect_identical(set$SIZE, .OBS_PER_CELL_DEFAULT)
})

test_that("every figure the help and the vignettes quote is a row of the simulation's output", {
  sim <- .docs_sim()
  if (is.null(sim)) skip("dev/sim-pooled-values.csv is in the source tree only")
  claims <- .docs_claims(sim)
  texts  <- .docs_all_texts()
  expect_length(claims, 23L)
  for (cl in claims) {
    expect_true(cl$where %in% names(texts), info = cl$where)
    expect_identical(.docs_count(texts[[cl$where]], cl$phrase), 1L,
                     info = sprintf("%s does not say, once: %s", cl$where, cl$phrase))
  }
  # A setting the page gives in words.
  expect_identical(sim$set$NEAR, 0.8)
  expect_identical(.docs_count(texts[["?summarize_by_cell: How this was checked"]],
                               "four fifths of them round"), 1L)
})

test_that("no number stands in those texts but the simulation's figures, its settings and the code's constants", {
  sim <- .docs_sim()
  if (is.null(sim)) skip("dev/sim-pooled-values.csv is in the source tree only")
  set    <- sim$set
  claims <- .docs_claims(sim)
  texts  <- .docs_all_texts()
  num <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  sci <- function(x) sub("e\\+0*", "e", format(x, scientific = TRUE))
  lev <- function(x) sprintf("%.2f", x)

  nb     <- num(.POOL_NEIGHBOURS)                  # the neighbours of a cell
  few    <- c("2", num(.POOL_NEIGHBOURS + 1L))     # "2 to 9 cells": every other cell is one
  size   <- num(.OBS_PER_CELL_DEFAULT)
  sizes  <- vapply(.OBS_PER_CELL_LADDER, num, character(1))
  level  <- lev(set$LEVEL)                         # the one level that was simulated
  layers <- num(set$LAYERS)
  layer  <- set$N_ROWS - set$N_BLOCK - set$N_AMONG # the rows of a layer
  cells  <- num(.count_at(layer, set$SIZE))
  margin <- num(100 * .STRUCTURE_TOL)              # "within 2 percent"
  f_test <- "0.05"                                 # the level of the F test of .no_structure()
  pairs  <- c(sci(.EXACT_PAIRS_MAX), num(100 * floor(sqrt(.EXACT_PAIRS_MAX) / 100)))
  dense  <- num(.KNN_DENSE_MAX)
  # Residuals that leave .FEW_RESIDUALS beyond each end of an interval.
  needed <- num(ceiling(.FEW_RESIDUALS / round((1 - set$LEVEL) / 2, 12)))

  expected <- list(
    "?summarize_by_cell: Pooled values" = c(
      nb, "0", "1",                    # the 8 neighbours; max(0, ); w * m + (1 - w) * t
      "1", "0",                        # an own weight of 1, and of 0
      few, "0", "0",                   # 2 to 9 cells; tau2 at 0 and every own weight 0
      margin, f_test,                  # what `structured` asks
      "1"),                            # (1 - w) * tau2
    "?summarize_by_cell: Intervals: which one answers which question" = c(
      "1", "1",                        # deff = 1, and a number other than 1
      "1", "2", "1",                   # p = (1 - obs_level) / 2; quantile(e, 1 - p)
      level),                          # the level checked
    "?summarize_by_cell: Conditions" = c(
      pairs, "1", "0", "0",            # the limit; own weights of 1 and of 0; tau2 is 0
      few, "0", "2", "0",              # 2 to 9 cells; tau2 at 0, cells 2 and own_weight 0
      num(.POOL_MIN_DF), level, needed,
      sci(.POOL_XY_MAX),               # a row beyond the bound has no location
      "90"),                           # nor has one at a latitude beyond 90 degrees
    "?summarize_by_cell: How this was checked" = c(
      layers, num(set$SIDE), num(set$N_ROWS), num(set$PARENTS), num(set$SPREAD),
      num(set$RANGE), num(set$NOISE_VAR), num(set$N_BLOCK), num(set$N_AMONG), num(layer),
      cells, num(set$SIZE),
      level, level,                    # obs_level and conf_level in the run
      layers, nb),                     # the white-noise layers; the 8 nearest cells
    "?summarize_by_cell: Confidence intervals, the two intervals" = character(0),
    "?summarize_by_cell: argument estimate" = c(
      "1", nb, num(.POOL_NEIGHBOURS + 2L)),        # fewer than 10 cells: all the others
    "?summarize_by_cell: argument obs_level" = c(
      "0", "1", level,                 # a level in (0, 1), such as 0.90
      "1", "2", "1", "2",              # the (1 - obs_level) / 2 and (1 + obs_level) / 2 quantiles
      "0", "1", level),                # strictly between 0 and 1; only 0.90 checked
    "?resolution_profile: Observations per cell" = c(
      "2",                             # max(2, round(n / s))
      size, "1", "0",                  # 4.5 unless given; an own weight of 1, and of 0
      size, sizes,                     # what "auto" is
      sci(.FIT_WORK_MAX), pairs, dense, dense),
    "?resolution_profile: The ladder and its bounds, min_cell_n" =
      c(size, num(floor(.OBS_PER_CELL_DEFAULT))),  # "4.5 is read as 4"
    "?resolution_profile: The ladder and its bounds, a ladder of sizes" = "2",
    "?resolution_profile: The ladder and its bounds, no variogram" = character(0),
    "?resolution_profile: The criteria, loo" = character(0),
    "?resolution_profile: How this was checked" = c(
      layers, num(layer), num(set$RANGE), num(set$SIDE), sizes),
    "?select_resolution: What each criterion returns, loo" = f_test,
    "?select_resolution: What each criterion returns, size" = c("2", f_test),
    "?select_resolution: What each criterion returns, auto" = character(0),
    "?select_resolution: A profile built with an estimate" = character(0),
    "?cell_values_at: Away from the observations" = character(0),
    "?cell_values_at: How this was checked" = c(
      layers, num(layer), cells, num(set$N_AMONG), num(set$N_BLOCK), level),
    "?spatialkit: the two choices" = c(nb,         # the 8 nearest cells
                                       size),      # the default size
    "?spatialkit: the conventions of pooled values" = c(
      size, nb, num(.FEW_RESIDUALS), num(.POOL_MIN_DF), margin, sci(.EXACT_PAIRS_MAX),
      sci(.FIT_WORK_MAX)),
    "vignette resolution: A size for pooled values" = c(
      nb, "2", "1", "0",               # the 8 cells; max(2, ); an own weight of 1, and of 0
      size,                            # the default size
      f_test,                          # what `structured` asks
      layers, num(layer), num(max(set$SIZES)), num(min(set$SIZES))),
    "vignette resolution: Leave-one-out for plain means" = sizes,
    "vignette reporting: Reporting a pooled cell value" = c(
      nb, "600", "20",                 # the section's own layer, as its chunk draws it
      size, level),
    "vignette cells: What summarize_by_cell() returns" = c(
      "1", nb),                        # deff other than 1; the 8 nearest cells
    "vignette cells: Pooled values and their warnings" = c(
      size, margin, f_test, layers, nb,
      "1", "0", few),                  # own weights of 1 and of 0; 2 to 9 cells
    "vignette cells: Values at other locations" = character(0))

  expect_setequal(names(expected), .DOCS_WRITTEN)
  for (what in names(expected)) {
    expect_true(what %in% names(texts), info = what)
    rest <- .docs_rest(texts[[what]], claims, what)
    expect_identical(sort(.docs_numbers(rest)), sort(unname(expected[[what]])), info = what)
  }
  # A figure stands only in a text whose numbers are held.
  expect_true(all(vapply(claims, function(cl) cl$where, character(1)) %in% .DOCS_WRITTEN))
  # The layer of the reporting vignette's section is the one its chunk draws.
  code <- .docs_vignette_texts()[["vignette reporting: Reporting a pooled cell value"]]
  expect_identical(.docs_count(code, "sx <- runif(600, 0, 2e4)"), 1L)
})

test_that("a text that says how a number tends to behave quotes a figure of the simulation, or sends the reader to the section that does", {
  sim <- .docs_sim()
  if (is.null(sim)) skip("dev/sim-pooled-values.csv is in the source tree only")
  texts  <- .docs_all_texts()
  quoted <- unique(vapply(.docs_claims(sim), function(cl) cl$where, character(1)))
  hedges <- lapply(stats::setNames(nm = .DOCS_WRITTEN), function(what)
    tolower(.docs_hedges(texts[[what]])))
  hedges <- hedges[lengths(hedges) > 0L]
  # A text that sends the reader on names the section by its title, and that
  # section quotes figures.
  expect_true(all(names(.DOCS_SEE) %in% .DOCS_WRITTEN))
  expect_true(all(.DOCS_SEE %in% quoted))
  for (what in names(.DOCS_SEE))
    expect_gte(.docs_count(texts[[what]], "\"How this was checked\""), 1L, label = what)
  for (what in names(hedges))
    expect_true(what %in% quoted || what %in% names(.DOCS_SEE),
                info = sprintf("%s says '%s' and quotes no figure", what,
                               paste(hedges[[what]], collapse = "', '")))
  # The statements of that kind these texts make, by their words: the
  # pooled values' error falling to the finest size (output f), a pooled
  # value being mostly its neighbours' mean (output d), and an empty cell
  # coming more often from random starts.  One more is one more claim, to be
  # written here beside the output that shows it.
  expect_identical(hedges, list(
    "?resolution_profile: Observations per cell"           = c("tends", "tends", "often"),
    "?resolution_profile: The criteria, loo"               = "tends",
    "?select_resolution: What each criterion returns, loo" = "tends",
    "vignette resolution: A size for pooled values"        = c("tends", "tends"),
    "vignette cells: Pooled values and their warnings"     = "mostly"))
})


# ---- Sentences that stated something false about the code -------------------------------
# Each test reads one sentence that described a branch or a figure the code
# does not have, pins the sentence that replaced it, and checks the code it
# describes.

test_that("the vignette sends the edge rule for gridded locations to the help, which gives it with its condition on the rows' CRS", {
  prose <- .docs_vignette_texts(code = FALSE)
  if (is.null(prose)) skip("the vignettes' sources are not installed")
  vw <- prose[["vignette cells: Pooled values and their warnings"]]
  expect_identical(.docs_count(vw, paste0(
    "Points on a boundary between two cells, as on gridded data, are placed so that a ",
    "lattice passes the check; ?summarize_by_cell (\"Value\") lists the cases where they ",
    "are not.")), 1L)
  # The sentence without its condition is gone, and so is the one that gave
  # every row without a CRS to the polygons: rows without a CRS on cells
  # whose record names none take the profile's cell (test-core-gridded-ties.R).
  expect_identical(.docs_count(vw, "centres scored, so a lattice passes the check"), 0L)
  expect_identical(.docs_count(vw, "Rows without a CRS, rows assigned"), 0L)
  # The help states the condition.
  help <- .docs_squash(.docs_flat(.docs_rd("summarize_by_cell")))
  expect_identical(.docs_count(help, paste0(
    "for rows that have a CRS exactly when the record names one. Rows without a CRS beside ",
    "a record that names one, rows with a CRS beside a record that names none, rows assigned ",
    "to the cells after they lost the record (written to a file and read back, say) and ",
    "summarised with cells_sf set to the cells that carry it, and rows given their cells by ",
    "another join are placed by the polygons alone")), 1L)
  # In the "Value" section the vignette names.
  value <- .docs_rd_part(.docs_rd("summarize_by_cell"), "\\value")
  expect_identical(.docs_count(value, "for rows that have a CRS exactly when the record names one."),
                   1L)
})

test_that("the vignette says the inner folds of forward selection are block_kfold unless `method` names another scheme", {
  expect_identical(eval(formals(select_features_forward)$method)[1L], "block_kfold")
  lines <- .docs_vignette("models")
  if (is.null(lines)) skip("the vignettes' sources are not installed")
  sec <- .docs_md_sections(lines, depth = 2L, code = FALSE)[[
    "Choosing predictors inside the folds"]]
  expect_false(is.null(sec))
  expect_identical(.docs_count(sec, "default to the scheme you pass"), 0L)
  expect_identical(.docs_count(sec, paste0(
    "The selection's inner folds are block_kfold unless method names another scheme, ",
    "whatever folds the outer loop uses. Passing random_kfold there logs a caution:")), 1L)

  # The code: random outer folds handed to cv_spatial() leave the inner
  # folds of a selection that names no `method` at block_kfold.
  set.seed(11)
  n  <- 90
  xy <- data.frame(x = 5e5 + runif(n, 0, 1000), y = 5e6 + runif(n, 0, 1000))
  xy$a <- rnorm(n); xy$b <- rnorm(n)
  xy$z <- 1.5 * xy$a + rnorm(n, sd = 0.3)
  pts <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  lm_fit <- function(train_sf, vars) lm_spatial_fit(train_sf, "z", vars)
  inner <- character(0)
  fit_sel <- function(train_sf) {
    s <- select_features_forward(train_sf, "z", c("a", "b"), fit_fn = lm_fit, k = 3,
                                 block_size = 400, seed = 1, quiet = TRUE)
    inner <<- c(inner, s$params$method)
    lm_fit(train_sf, s$selected)
  }
  outer <- make_folds(pts, k = 3, method = "random_kfold", seed = 2)
  suppressWarnings(cv_spatial(pts, "z", c("a", "b"), fit_fn = fit_sel, folds = outer))
  expect_identical(inner, rep("block_kfold", 3L))
})

test_that("?print.spatial_fit names the LOOIC a bayesian_fit prints and does not call it in-sample", {
  rd <- .docs_rd("print.spatial_fit")
  desc <- .docs_squash(.docs_flat(.docs_rd_part(rd, "\\description")))
  expect_identical(.docs_count(desc, "Nothing printed here is an out-of-sample score."), 0L)
  expect_identical(.docs_count(desc, paste0(
    "the LOOIC when the fit computed it (compute_loo = TRUE, the default of ",
    "fit_bayesian_spatial_model()), the slope prior")), 1L)
  expect_identical(.docs_count(desc, paste0(
    "Apart from a forest's out-of-bag error and a Bayesian GP's LOOIC, an approximate ",
    "leave-one-out score of single rows, nothing printed here is an out-of-sample score, ",
    "and neither of those is a spatial one.")), 1L)

  # The code: the print shows the LOOIC when the fit holds one, and no LOOIC
  # line when it does not.
  pts <- sf::st_as_sf(data.frame(x = 5e5 + 1:10, y = 4e6 + (1:10)^2, a = 1:10,
                                 z = sin(1:10)), coords = c("x", "y"), crs = 32617)
  fit <- function(looic) new_spatial_fit("bayesian_fit", list(), z ~ a, "z", "a", pts,
                                         info = list(gp_k = 5L, gp_n_basis = 25L, looic = looic))
  expect_true(any(grepl("^  LOOIC   : 81\\.23$", capture.output(print(fit(81.23))))))
  expect_false(any(grepl("LOOIC", capture.output(print(fit(NA_real_))), fixed = TRUE)))
})

test_that("the getting-started vignette gives the distance errors nc.shp's projection log gives", {
  lines <- .docs_vignette("getting-started")
  if (is.null(lines)) skip("the vignettes' sources are not installed")
  text <- paste(lines, collapse = "\n")
  at <- regexpr("(?s)chose a local Albers equal-area projection and logged why:.*?Pass `target_crs`",
                text, perl = TRUE)
  expect_gt(at, 0L)
  para <- if (at > 0L) regmatches(text, at) else ""
  # The figures are written from the object the chunk above builds, not as
  # literals: the measured error moves with the PROJ build.
  expect_false(grepl("0.20 percent", para, fixed = TRUE))
  code <- sub("^`r ", "", sub("`$", "", regmatches(para, gregexpr("`r [^`]+`", para))[[1L]]))
  expect_length(code, 2L)

  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  env <- new.env()
  log <- capture_spatialkit_log(
    env$counties <- ensure_projected(nc[, c("BIR74", "SID74", "NWBIR74")]),
    level = logger::WARN)
  line <- grep("against the zone's", log, value = TRUE, fixed = TRUE)
  expect_length(line, 1L)
  logged <- regmatches(line, regexec(
    "worst-case distance error ([0-9.]+)% against the zone's ([0-9.]+)%", line))[[1L]]
  # The sentence names the projection chosen as a local Albers equal-area
  # one and the UTM zone as the other candidate quoted: a build that chose
  # another would leave the words false while the figures followed it.
  choice <- attr(env$counties, "crs_choice")
  expect_true(startsWith(choice$name[choice$chosen], "Albers equal-area"))
  expect_identical(sum(startsWith(choice$name, "UTM zone")), 1L)
  shown <- vapply(code, function(cd) as.character(eval(parse(text = cd), env)), "")
  # The zone first, then the projection chosen, as the sentence reads.
  expect_identical(unname(shown), logged[c(3L, 2L)])
})

test_that("?sac_nugget gives nugget_at_bound's bound as a share of the sill, and ?estimate_sac_range sends the reader there", {
  bound <- function(page) .docs_squash(.docs_flat(.docs_rd(page)))
  nug <- bound("sac_nugget")
  est <- bound("estimate_sac_range")
  # "Within 1e-4 of the sill" reads as a nugget close to the sill, or 1e-4
  # from it in the response's units; the code marks a share of the sill.
  expect_identical(.docs_count(nug, "within 1e-4"), 0L)
  expect_identical(.docs_count(est, "within 1e-4"), 0L)
  expect_identical(.docs_count(nug, paste0(
    "nugget_at_bound is TRUE when the fit put the nugget at its lower bound (0, or at most ",
    "1e-4 of the sill). That value is the bound and not an estimate")), 1L)
  expect_identical(.docs_count(est, paste0(
    "nugget (the model's nugget variance), nugget_at_bound, nugget_floor, ",
    "short_lag_semivariance, short_lag_pairs, short_lag_dist, zero_lag_semivariance and ",
    "zero_lag_pairs: the nugget and the evidence to judge it by, described in sac_nugget().")), 1L)
  # The wording ?resolution_profile uses for the same mark.
  expect_identical(.docs_count(bound("resolution_profile"),
                               "A nugget of 0, or of at most 10^{-4} of the sill"), 1L)

  # The code: the mark is nugget <= 1e-4 * sill, so it does not move with the
  # unit of the response.  A field with no noise has its nugget at 0, marked at
  # every scale; the same field with noise of sd 0.2 has a nugget near 8
  # percent of the sill, marked at none, although at a scale of 1e-3 it lies
  # within 1e-4 of the sill in the response's units.
  set.seed(5)
  n  <- 300
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  L  <- t(chol(exp(-as.matrix(stats::dist(xy)) / 120) + diag(1e-8, n)))
  z1 <- rnorm(n); z2 <- rnorm(n)
  p  <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  p$A <- drop(L %*% z2)
  p$B <- drop(L %*% z1) + 0.2 * z2
  for (v in c("A", "B")) for (k in c(1, 1e-3, 1e3)) {
    p$zk <- p[[v]] * k
    r <- NULL
    capture_spatialkit_log(r <- estimate_sac_range(p, "zk"), level = logger::WARN)
    sill <- sum(attr(r, "variogram_model")$psill)
    expect_identical(attr(r, "nugget_at_bound"), sac_nugget(r) <= 1e-4 * sill,
                     info = paste(v, k))
    expect_identical(attr(r, "nugget_at_bound"), v == "A", info = paste(v, k))
  }
})

test_that("estimate_sac_range() advises a larger `cutoff` for a range past the lags only where one is accepted", {
  # A linear trend: the variogram rises with no sill at every cutoff.
  set.seed(4)
  n  <- 300
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  xy$t <- xy$x / 50 + rnorm(n, sd = 0.3)
  p  <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  over <- function(cutoff, ...) {
    r <- NULL
    log <- capture_spatialkit_log(r <- estimate_sac_range(p, "t", cutoff = cutoff, ...),
                                  level = logger::WARN)
    expect_true(is.na(r))
    expect_identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted")
    line <- grep("exceeds the largest lag the variogram was fitted over", log,
                 value = TRUE, fixed = TRUE)
    expect_length(line, 1L)
    paste(line, collapse = "\n")
  }
  # At cutoff = 1 the lags already run to the farthest pair, and a larger
  # value is refused: the advice is to detrend or set a block size.
  at_one <- over(1)
  expect_false(grepl("Raise `cutoff`", at_one, fixed = TRUE))
  expect_true(grepl(paste0("Returning NA. The lags already run to the farthest pair of ",
                           "points (`cutoff` = 1); supply `predictor_vars` to detrend, or ",
                           "set a block size explicitly."), at_one, fixed = TRUE))
  expect_error(estimate_sac_range(p, "t", cutoff = 1.01),
               "a single number above 0 and at most 1", fixed = TRUE)
  # Detrended on a predictor that does not carry the trend, the advice to
  # try other predictors stands without the larger cutoff.
  p$noise <- rnorm(n)
  at_one_dt <- over(1, predictor_vars = "noise")
  expect_false(grepl("Raise `cutoff`", at_one_dt, fixed = TRUE))
  expect_true(grepl("(`cutoff` = 1); try predictors that carry the trend", at_one_dt,
                    fixed = TRUE))
  # Below 1 a larger cutoff is accepted, and the advice keeps it.
  for (cf in c(0.5, 0.9, 0.99))
    expect_true(grepl(paste0("Returning NA. Raise `cutoff` to fit longer lags, supply ",
                             "`predictor_vars` to detrend, or set a block size explicitly."),
                      over(cf), fixed = TRUE), info = paste("cutoff", cf))
})

test_that("the over-cutoff warning names the n_max subsample at cutoff = 1", {
  # With n above n_max the lags run to the farthest pair of the subsample,
  # not of every point, and the warning says which.
  set.seed(4)
  n  <- 300
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  xy$t <- xy$x / 50 + rnorm(n, sd = 0.3)
  p  <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  r <- NULL
  log <- capture_spatialkit_log(r <- estimate_sac_range(p, "t", cutoff = 1, n_max = 100),
                                level = logger::WARN)
  expect_true(is.na(r))
  expect_identical(attr(r, "rejected_reason"), "fitted range exceeds the largest lag fitted")
  expect_true(any(grepl(paste0("The lags already run to the farthest pair of the `n_max` ",
                                  "subsample of points (`cutoff` = 1)"), log, fixed = TRUE)))
  log_all <- capture_spatialkit_log(estimate_sac_range(p, "t", cutoff = 1), level = logger::WARN)
  expect_false(any(grepl("subsample", log_all, fixed = TRUE)))
})

test_that("a range inside the lags refused by range_frac below 1 is not said to run past them", {
  # An exponential field of scale 150 m on a 1 km square: its variogram sills
  # inside the lags at both cutoffs.  A range_frac of 0.2 refuses that range
  # although it lies within the largest lag fitted, so the warning must not
  # say it exceeds that lag or that the variogram never reached a sill, and
  # the remedy is a larger range_frac, not a larger cutoff.  The help of
  # range_frac says the bound is the longest lag only at its default.
  est <- .docs_squash(.docs_flat(.docs_rd("estimate_sac_range")))
  expect_identical(.docs_count(est, paste0(
    "range_frac * cutoff * max_dist (at the default, beyond the longest lag the ",
    "empirical variogram was fitted over) is treated as unidentified")), 1L)
  expect_identical(.docs_count(est, "(that is, beyond the longest lag"), 0L)
  set.seed(5)
  n  <- 300
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  L  <- t(chol(exp(-as.matrix(stats::dist(xy)) / 150) + diag(1e-8, n)))
  p  <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  p$A <- drop(L %*% rnorm(n))
  for (cf in c(0.5, 1)) {
    ok <- NULL
    capture_spatialkit_log(ok <- estimate_sac_range(p, "A", cutoff = cf), level = logger::WARN)
    expect_true(is.finite(ok), info = paste("cutoff", cf))
    r <- NULL
    log <- capture_spatialkit_log(r <- estimate_sac_range(p, "A", cutoff = cf, range_frac = 0.2),
                                  level = logger::WARN)
    # The value and the reason are as before: only the text changes.
    expect_true(is.na(r), info = paste("cutoff", cf))
    expect_identical(attr(r, "rejected_reason"),
                     "fitted range exceeds the largest lag fitted", info = paste("cutoff", cf))
    expect_equal(attr(r, "rejected_range"), as.numeric(ok), info = paste("cutoff", cf))
    expect_lte(attr(r, "rejected_range"), attr(r, "cutoff_dist"))
    expect_gt(attr(r, "rejected_range"), 0.2 * attr(r, "cutoff_dist"))
    line <- grep("estimate_sac_range(): fitted range", log, value = TRUE, fixed = TRUE)
    expect_length(line, 1L)
    expect_false(grepl("never reached a sill", line, fixed = TRUE), info = paste("cutoff", cf))
    expect_false(grepl("exceeds the largest lag", line, fixed = TRUE), info = paste("cutoff", cf))
    expect_false(grepl("Raise `cutoff`", line, fixed = TRUE), info = paste("cutoff", cf))
    expect_true(grepl(sprintf(paste0(
      "fitted range (%.4g) exceeds the largest range `range_frac` accepts (%.4g = 0.2 x ",
      "the largest lag fitted, %.4g), although it lies within that lag. Returning NA"),
      attr(r, "rejected_range"), 0.2 * attr(r, "cutoff_dist"), attr(r, "cutoff_dist")),
      line, fixed = TRUE), info = paste("cutoff", cf))
    expect_true(grepl("Raise `range_frac` (at 1 it accepts any range within the lags), or set a block size explicitly.",
                      line, fixed = TRUE), info = paste("cutoff", cf))
  }
  # A range past the largest lag keeps the old text with range_frac below 1.
  set.seed(4)
  xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000))
  xy$t <- xy$x / 50 + rnorm(n, sd = 0.3)
  q  <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32632)
  r <- NULL
  log <- capture_spatialkit_log(r <- estimate_sac_range(q, "t", range_frac = 0.5),
                                level = logger::WARN)
  expect_gt(attr(r, "rejected_range"), attr(r, "cutoff_dist"))
  expect_true(any(grepl("exceeds the largest lag the variogram was fitted over", log, fixed = TRUE)))
  expect_true(any(grepl("never reached a sill", log, fixed = TRUE)))
  expect_true(any(grepl("Raise `cutoff` to fit longer lags", log, fixed = TRUE)))
})

test_that("?sac_nugget says a nugget at a floor of at most 1e-4 of the sill is marked", {
  txt <- .docs_squash(.docs_flat(.docs_rd("sac_nugget")))
  expect_identical(.docs_count(txt, paste0(
    "A nugget equal to its floor (see nugget_floor) is the variance between visits, a ",
    "measurement, and nugget_at_bound stays FALSE for it unless the floor itself is at ",
    "most 1e-4 of the sill.")), 1L)
  # The code: 600 locations visited twice, on a field with no noise, so the
  # fit puts the nugget at its floor, the variance between visits.  With a
  # visit sd of 1e-3 that floor is about 1e-6 of the sill and the nugget is
  # marked; with 0.05 it is about 2.5e-3 of the sill and it is not.
  one <- function(visit_sd) {
    set.seed(5203)
    n <- 600L
    x <- runif(n, 0, 400e3); y <- runif(n, 0, 300e3)
    D <- as.matrix(stats::dist(cbind(x, y)))
    f <- drop(crossprod(chol(exp(-3 * D / 30e3) + diag(1e-8, n)), rnorm(n)))
    set.seed(1)
    d <- rbind(data.frame(x = 5e5 + x, y = 4e6 + y, z = f + rnorm(n, 0, visit_sd)),
               data.frame(x = 5e5 + x, y = 4e6 + y, z = f + rnorm(n, 0, visit_sd)))
    r <- NULL
    capture_spatialkit_log(
      r <- estimate_sac_range(sf::st_as_sf(d, coords = c("x", "y"), crs = 32615), "z"),
      level = logger::WARN)
    r
  }
  for (vs in c(1e-3, 0.05)) {
    r <- one(vs)
    sill <- sum(attr(r, "variogram_model")$psill)
    flo  <- attr(r, "nugget_floor")
    expect_gt(flo, 0)
    expect_equal(sac_nugget(r), flo, tolerance = 1e-12, info = paste("visit sd", vs))
    expect_identical(attr(r, "nugget_at_bound"), flo <= 1e-4 * sill, info = paste("visit sd", vs))
    expect_identical(attr(r, "nugget_at_bound"), vs == 1e-3, info = paste("visit sd", vs))
  }
})
