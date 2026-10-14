# tests/testthat/test-core-docs.R
# -----------------------------------------------------------------------------
# What the help, the vignettes and NEWS written for pooled cell values may say
# -----------------------------------------------------------------------------
# The help pages, the NEWS entries and the vignette sections that describe the
# pooled value, the interval for one observation, the size and the record of
# the count are read as text.  Five rules are held:
#
#   * none of them names existing units (counties, ZIP code areas, ZCTAs) or
#     the pilot study: the pages do not say how the package's regions compare
#     with any others;
#   * "calibrated" stands in them only in the sentence that denies it for
#     places away from the observations, and that sentence stands wherever
#     the intervals are described;
#   * every figure they quote is a row of dev/sim-pooled-values.csv, the
#     output of the package's own simulation, and that file is the output of
#     dev/sim-pooled-values.R as it stands;
#   * no other number stands in them but the settings of that simulation and
#     the constants of the code, each read from where it is set;
#   * a text that says how a number tends to behave quotes a figure for it.
#
# The word rules read whole help pages, so a sentence added anywhere on a
# page this work wrote into is read.  The parts of a page are found by their
# titles, labels or first words, and each is held to its opening words, so a
# reader that returns another part fails.
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

# The words of a text that name existing units or the pilot study.
.docs_barred <- function(text)
  unlist(regmatches(text, gregexpr("(?i)\\b(county|counties|zip|zcta|pilot)\\w*", text,
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
# and exponent.  Digits inside a name (tau2, sigma2) and exponents written
# with "^" are not numbers of the text.
.docs_numbers <- function(text)
  regmatches(text, gregexpr(
    "(?<![A-Za-z0-9_.^])[0-9]+(?:,[0-9]{3})*(?:\\.[0-9]+)?(?:e[+-]?[0-9]+)?(?![A-Za-z0-9_])",
    text, perl = TRUE))[[1L]]

# The words by which a text says how a number tends to behave.
.docs_hedges <- function(text)
  unlist(regmatches(text, gregexpr(
    "(?i)\\b(tends?|tended|tendency|usually|typically|generally|often|mostly|in practice)\\b",
    text, perl = TRUE)))

# ---- What is read -------------------------------------------------------------------

# The help pages this work wrote into.  The word rules are held on the whole
# text of each.
.DOCS_PAGES <- c("summarize_by_cell", "resolution_profile", "select_resolution",
                 "cell_values_at", "get_voronoi_seeds", "voronoi_seeds_kmeans",
                 "build_tessellation", "create_voronoi_polygons",
                 "assign_features_to_polygons", "spatialkit-package",
                 "summary.resolution_profile", "plot.resolution_profile")

# What those pages said with a word of the stem "calibrat" before this work,
# each in its words: the elbow rule and Moran's z in ?resolution_profile, the
# Gaussian process in ?spatialkit, and the name of a plotting function.  They
# are the only exceptions to the rule, and each stands once.
.DOCS_OLDER <- list(
  resolution_profile        = c("the rule and its calibration are in determine_optimal_levels",
                                "z is calibrated and flat in L"),
  "spatialkit-package"      = "an explicit spatial Gaussian process with calibrated uncertainty",
  plot.resolution_profile   = "plot_calibration()")

.docs_pages <- function() {
  out <- list()
  for (page in .DOCS_PAGES) {
    rd <- .docs_rd(page)
    out[[paste0("?", page)]] <- if (is.null(rd)) NA_character_ else .docs_squash(.docs_flat(rd))
  }
  out
}

# The parts of those pages this work wrote, each found by its title, its
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
                  "Conditions", "What the standard error estimates", "Confidence intervals"))
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
  put("?resolution_profile: The criteria, cp",
      entry(rd, "The criteria, and what each returns in practice", "cp"))
  put("?resolution_profile: The criteria, loo",
      entry(rd, "The criteria, and what each returns in practice", "loo"))

  rd <- .docs_rd("select_resolution")
  for (label in c("cp", "loo", "size", "auto"))
    put(paste0("?select_resolution: What each criterion returns, ", label),
        entry(rd, "What each criterion returns", label))
  put("?select_resolution: A profile built with an estimate",
      section(rd, "A profile built with an estimate"))

  rd <- .docs_rd("cell_values_at")
  put("?cell_values_at: Away from the observations", section(rd, "Away from the observations"))

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

# The six conventions this work added to "Defaults that were chosen", by
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
  "?summarize_by_cell: What the standard error estimates" =
    "What the ..se_* columns estimate follows deff.",
  "?summarize_by_cell: Confidence intervals" =
    "With conf_level set, every numeric response and predictor column gains four more",
  "?summarize_by_cell: Confidence intervals, the two intervals" =
    "At no setting is ..ci_* a range for one observation",
  "?summarize_by_cell: argument estimate" =
    "Which value of a cell the frame reports: \"mean\" (default) or \"pooled\".",
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
  "?select_resolution: What each criterion returns, cp" =
    "The level that minimises Mallows' C_p: the count at which re-measuring the sampled sites",
  "?select_resolution: What each criterion returns, loo" =
    "The level with the lowest leave-one-out error of the estimate the profile was built for",
  "?select_resolution: What each criterion returns, size" =
    "The count a size gives, max(2, round(n / obs_per_cell))",
  "?select_resolution: What each criterion returns, auto" =
    "The rule that goes with the profile's estimate",
  "?select_resolution: A profile built with an estimate" =
    "This function's default criterion is \"cp\", whatever the profile.",
  "?cell_values_at: Away from the observations" =
    "The value and the interval describe locations among the observations.",
  "?spatialkit: the two choices" =
    "Two choices are the user's, one in choosing a resolution (step 1) and one in aggregating",
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
                     "vignette reporting: Reporting a cell value")

# The NEWS entries of this work, by their first words.
.DOCS_NEWS <- c(
  "summarize_by_cell() gains estimate = \"pooled\".",
  "summarize_by_cell() gains obs_level:",
  "resolution_profile() gains estimate and obs_per_cell.",
  "resolution_profile()'s sac accepts \"auto\":",
  "resolution_profile()'s nstart, get_voronoi_seeds()'s kmeans_nstart and",
  "select_resolution() gains criterion = \"size\"",
  "A selection made by these rules carries count_rule",
  "New cell_values_at():",
  "Each warning of the pooled path of summarize_by_cell() is raised once",
  "New condition classes: spatialkit_pool_unavailable",
  "?resolution_profile describes cp as the error of re-measuring the sampled sites.",
  "?summarize_by_cell says what ..se_* and ..ci_* cover at each deff:")

# Those entries, named "NEWS: <first words>"; NULL when NEWS is not there.
.docs_news_texts <- function() {
  lines <- .docs_news()
  if (is.null(lines)) return(NULL)
  entries <- .docs_news_entries(lines)
  out <- list()
  for (first in .DOCS_NEWS) out[[paste0("NEWS: ", first)]] <- .docs_starting(entries, first)
  out
}

# The sections of the vignettes this work wrote, and the two older sections
# it added sentences to; with them the front matter, whose description line
# it edited.  `code = FALSE` leaves the code chunks out, and the openings
# below are those of the prose.
.DOCS_VIGNETTES <- list(
  resolution  = list("2" = c("Setting a size", "Leave-one-out for plain means"),
                     "3" = "Three criteria and a diagnostic"),
  reporting   = list("2" = "Reporting a cell value"),
  diagnostics = list("2" = c("Four questions a score does not answer",
                             "6. Warnings about pooled cell values")))
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
  "vignette resolution: Setting a size" =
    "The criteria above score the plain mean of a cell.",
  "vignette resolution: Leave-one-out for plain means" =
    "With plain means the leave-one-out error does choose a count.",
  "vignette resolution: Three criteria and a diagnostic" =
    "elbow sees the coordinates and knows nothing about what you measured.",
  "vignette reporting: front matter"   = "title: \"What to save, and what to report\"",
  "vignette reporting: Reporting a cell value" =
    "summarize_by_cell(estimate = \"pooled\") returns, beside each cell's mean, a pooled value",
  "vignette diagnostics: front matter" = "title: \"Diagnostics\"",
  "vignette diagnostics: Four questions a score does not answer" =
    "A cross-validated RMSE tells you how far the predictions fell from the held-out values.",
  "vignette diagnostics: 6. Warnings about pooled cell values" =
    "summarize_by_cell(estimate = \"pooled\") returns, beside each cell's mean, a pooled value")

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

# A count of layers as the pages write it, and a small count in words.
.docs_layers <- function(k) if (k == 0L) "none" else if (k == 1L) "one" else as.character(k)
.docs_word   <- function(k) c("one", "two", "three", "four", "five", "six", "seven", "eight",
                              "nine", "ten")[k]

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
  sizes   <- .docs_word(length(set$SIZES))

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

  pv  <- "?summarize_by_cell: Pooled values"
  iv  <- "?summarize_by_cell: Intervals: which one answers which question"
  ci  <- "?summarize_by_cell: Confidence intervals, the two intervals"
  op  <- "?resolution_profile: Observations per cell"
  lo  <- "?resolution_profile: The criteria, loo"
  sr  <- "?select_resolution: What each criterion returns, loo"
  va  <- "?cell_values_at: Away from the observations"
  n1  <- "NEWS: summarize_by_cell() gains estimate = \"pooled\"."
  n2  <- "NEWS: summarize_by_cell() gains obs_level:"
  vs  <- "vignette resolution: Setting a size"
  vl  <- "vignette resolution: Leave-one-out for plain means"
  vr  <- "vignette reporting: Reporting a cell value"
  vd  <- "vignette diagnostics: 6. Warnings about pooled cell values"
  claim <- function(where, ...) list(where = where, phrase = sprintf(...))
  list(
    claim(pv, "was %s of the plain mean's", a_all),
    claim(pv, paste("%s in cells that hold one of the layer's observations, %s in",
                    "cells of 2 to 4, %s in cells of 5 to 9 and %s in cells of 10 or more"),
          a_cl[1], a_cl[2], a_cl[3], a_cl[4]),
    claim(pv, "In the held-out block it was %s.", a_block),
    claim(pv, "the median own weight was %s", d_w),
    claim(pv, "tau2 was 0 on %d layers", d_tau0),
    claim(pv, "the warning was raised on all %d", d_n),
    claim(pv, "a correlation of %s with", d_cor),
    claim(pv, "where the plain means had one of %s", d_corm),
    claim(pv, "was %s of (1 - w) * tau2 among the observations and %s times it",
          e_among, e_block),
    claim(iv, paste("held %s percent of the observations held out at random among the",
                    "others (%s to %s percent in single layers) and %s percent of those in",
                    "the held-out block"), b_among, b_lo, b_hi, b_block),
    claim(iv, "held %s and %s percent of them", c_among, c_block),
    claim(ci, "held %s percent of the observations held out among the others", c_among),
    claim(ci, "at obs_level = 0.90 held %s percent", b_among),
    claim(op, "was at the size 2 on %d layers, at 3 on %d and at 4.5 on %d",
          f_at("pooled value", 2), f_at("pooled value", 3), f_at("pooled value", 4.5)),
    claim(op, paste("That of the plain mean was at 2 on %s: it was at 6 on %d layers, at",
                    "4.5 on %d, at 9 on %d, at 12 on %d, and at 18 and at 3 on %s each."),
          .docs_layers(f_at("plain mean", 2)), f_at("plain mean", 6), f_at("plain mean", 4.5),
          f_at("plain mean", 9), f_at("plain mean", 12), .docs_layers(f_at("plain mean", 18))),
    claim(lo, paste("on a ladder of the %s sizes %s to %s, its minimum was at the finest",
                    "size on %d of %d layers for pooled values and on %s for plain means"),
          sizes, max(set$SIZES), finest, f_at("pooled value", finest), f_n,
          .docs_layers(f_at("plain mean", finest))),
    claim(sr, paste("was at the finest of %s sizes on %d of %d layers, and that of the",
                    "plain mean on %s"), sizes, f_at("pooled value", finest), f_n,
          .docs_layers(f_at("plain mean", finest))),
    claim(va, paste("held %s percent of those with a dist_ratio below 1, %s percent from 1",
                    "to 2, %s percent from 2 to 4, %s percent from 4 to 8 and %s percent",
                    "above 8"), b_ratio[1], b_ratio[2], b_ratio[3], b_ratio[4], b_ratio[5]),
    claim(va, paste("It held %s percent of the observations held out at random and %s",
                    "percent of those in the block"), b_among, b_block),
    claim(n1, "among the others was %s of the plain mean's", a_all),
    claim(n2, paste("held %s percent of the observations held out among the others and %s",
                    "percent of those in a held-out block"), b_among, b_block),
    claim(n2, "read that way held %s percent of the first", c_among),
    claim(vs, "the lowest loo of the pooled value was at the size 2 on %d (",
          f_at("pooled value", 2)),
    claim(vl, paste("the lowest loo of the plain mean was at the finest size on %s of the",
                    "%d layers: it was at 6 observations per cell on %d of them, at 4.5 on",
                    "%d and at 9 on %d"), .docs_layers(f_at("plain mean", finest)), f_n,
          f_at("plain mean", 6), f_at("plain mean", 4.5), f_at("plain mean", 9)),
    claim(vr, paste("held %s percent of held-out observations with a dist_ratio below 1",
                    "and %s percent of those above 8"), b_ratio[1], b_ratio[5]),
    claim(vd, "the median own weight was %s", d_w),
    claim(vd, "a correlation of %s with", d_cor),
    claim(vd, "where the plain means had one of %s", d_corm),
    claim(vd, paste("it came on %d of the %d white-noise layers and",
                    "spatialkit_no_structure on all %d"), d_tau0, d_n, d_n))
}

# Every text a figure can stand in: the parts of the help, the NEWS entries
# and the prose of the vignette sections.
.docs_all_texts <- function() c(.docs_parts(), .docs_news_texts(), .docs_vignette_texts(FALSE))

# A text with the run's statements taken out: what is left holds the
# numbers that are not figures of the simulation.
.docs_rest <- function(text, claims, where) {
  for (cl in claims) if (identical(cl$where, where)) text <- sub(cl$phrase, " ", text, fixed = TRUE)
  text
}

# The texts this work wrote whole: the ones whose numbers are all accounted
# for, and in which a statement of how a number tends to behave needs a
# figure.  The older sections it added sentences to are read by the word
# rules only.
.DOCS_WRITTEN <- c(
  "?summarize_by_cell: Pooled values",
  "?summarize_by_cell: Intervals: which one answers which question",
  "?summarize_by_cell: Conditions",
  "?summarize_by_cell: Confidence intervals, the two intervals",
  "?summarize_by_cell: argument estimate",
  "?summarize_by_cell: argument obs_level",
  "?resolution_profile: Observations per cell",
  "?resolution_profile: The ladder and its bounds, min_cell_n",
  "?resolution_profile: The ladder and its bounds, a ladder of sizes",
  "?resolution_profile: The ladder and its bounds, no variogram",
  "?resolution_profile: The criteria, loo",
  "?select_resolution: What each criterion returns, loo",
  "?select_resolution: What each criterion returns, size",
  "?select_resolution: What each criterion returns, auto",
  "?select_resolution: A profile built with an estimate",
  "?cell_values_at: Away from the observations",
  "?spatialkit: the two choices",
  "?spatialkit: the conventions of pooled values",
  paste0("NEWS: ", .DOCS_NEWS),
  "vignette resolution: Setting a size",
  "vignette resolution: Leave-one-out for plain means",
  "vignette reporting: Reporting a cell value",
  "vignette diagnostics: 6. Warnings about pooled cell values")


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
  expect_identical(.docs_barred("the pilot study; a Pilot; zip code areas; zcta5"),
                   c("pilot", "Pilot", "zip", "zcta5"))
  # Letters inside another word are not the word.
  expect_identical(.docs_barred("an unzipped archive, a copilot and a viscounty"),
                   character(0))
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
  # Digits inside a name, and an exponent, are not numbers of the text.
  expect_identical(.docs_numbers("sigma2, tau2, s^2, loo_se, ..ci_lo_* and k-means++"),
                   character(0))

  expect_identical(.docs_hedges("It tends to fall; in practice often, and usually so. Mostly."),
                   c("tends", "in practice", "often", "usually", "Mostly"))
  expect_identical(.docs_hedges("It extends the pretend, softened, usual cause."), character(0))
})


# ---- The help -------------------------------------------------------------------------

test_that("the parts of the help written for pooled cell values are where their titles, labels and first words say", {
  expect_length(.DOCS_PAGES, 12L)
  for (page in .DOCS_PAGES) expect_false(is.null(.docs_rd(page)), info = page)
  parts <- .docs_parts()
  expect_length(parts, 27L)
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
    "It can also return a pooled value for each cell, its mean pulled toward the mean of",
    "the cells around it, and a range for one further observation.")), 1L)
  expect_length(.DOCS_CONVENTIONS, 6L)
  for (first in .DOCS_CONVENTIONS)
    expect_identical(.docs_count(parts[["?spatialkit: the conventions of pooled values"]], first),
                     1L, info = first)
  # Every text held to its numbers is one of the parts, an entry of NEWS or
  # a vignette section.
  expect_true(all(.DOCS_WRITTEN %in% c(names(parts), paste0("NEWS: ", .DOCS_NEWS),
                                       names(.DOCS_VIGNETTE_OPENINGS))))
})

test_that("no help page this work wrote into names existing units or the pilot study", {
  pages <- .docs_pages()
  expect_identical(names(pages), paste0("?", .DOCS_PAGES))
  for (what in names(pages)) {
    expect_false(is.na(pages[[what]]), info = what)
    expect_gt(nchar(pages[[what]]), 2000L, label = what)
    expect_identical(.docs_barred(pages[[what]]), character(0), info = what)
  }
  # A part is a piece of its page, so the rule reads every part too.
  parts <- .docs_parts()
  for (what in names(parts)) {
    page <- pages[[sub("spatialkit$", "spatialkit-package", sub(":.*$", "", what))]]
    expect_identical(.docs_count(page, substr(parts[[what]], 1L, 200L)), 1L, info = what)
  }
})

test_that("on those pages \"calibrated\" stands only in the sentence that denies it away from the observations, and in what they said of other things before this work", {
  pages <- .docs_pages()
  expect_identical(names(.DOCS_OLDER),
                   c("resolution_profile", "spatialkit-package", "plot.resolution_profile"))
  expect_identical(lengths(.DOCS_OLDER, use.names = FALSE), c(2L, 1L, 1L))
  for (what in names(pages)) {
    text <- pages[[what]]
    for (phrase in .DOCS_OLDER[[sub("^\\?", "", what)]]) {
      expect_identical(.docs_count(text, phrase), 1L, info = paste(what, phrase))
      text <- sub(phrase, " ", text, fixed = TRUE)
    }
    expect_identical(.docs_calibrated(text), character(0), info = what)
  }
})

test_that("the help says, wherever it describes the intervals, that none is calibrated away from the observations", {
  expect_length(.DOCS_DENIED_IN, 5L)
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


# ---- NEWS -----------------------------------------------------------------------------

test_that("the NEWS entries of this work are found, name no existing units and claim no calibration", {
  news <- .docs_news_texts()
  if (is.null(news)) skip("NEWS.md is not there to read")
  expect_length(.DOCS_NEWS, 12L)
  expect_identical(names(news), paste0("NEWS: ", .DOCS_NEWS))
  for (what in names(news)) {
    expect_false(is.na(news[[what]]), info = what)
    expect_identical(.docs_barred(news[[what]]), character(0), info = what)
    expect_identical(.docs_calibrated(news[[what]]), character(0), info = what)
  }
  expect_identical(anyDuplicated(unlist(news)), 0L)
})


# ---- The vignettes ---------------------------------------------------------------------

test_that("the vignette sections written for pooled cell values are where their titles say", {
  prose <- .docs_vignette_texts(code = FALSE)
  if (is.null(prose)) skip("the vignettes' sources are not installed")
  whole <- .docs_vignette_texts()
  expect_length(prose, 9L)
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
  # The sentences added to the two older sections are in the texts read.
  expect_identical(.docs_count(prose[["vignette resolution: Three criteria and a diagnostic"]],
                               "What cp estimates is the error of re-measuring the sampled sites"),
                   1L)
  expect_identical(.docs_count(prose[["vignette diagnostics: Four questions a score does not answer"]],
                               "a sixth reads the warnings summarize_by_cell() raises about pooled cell values"),
                   1L)
})

test_that("no vignette section written for pooled cell values names existing units or the pilot study", {
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
  expect_length(grep("^vignette ", .DOCS_DENIED_IN), 1L)
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

test_that("every figure the help, the vignettes and NEWS quote is a row of the simulation's output", {
  sim <- .docs_sim()
  if (is.null(sim)) skip("dev/sim-pooled-values.csv is in the source tree only")
  claims <- .docs_claims(sim)
  texts  <- .docs_all_texts()
  expect_length(claims, 29L)
  for (cl in claims) {
    expect_true(cl$where %in% names(texts), info = cl$where)
    expect_identical(.docs_count(texts[[cl$where]], cl$phrase), 1L,
                     info = sprintf("%s does not say, once: %s", cl$where, cl$phrase))
  }
  # A setting the page gives in words.
  expect_identical(sim$set$NEAR, 0.8)
  expect_identical(.docs_count(texts[["?summarize_by_cell: Pooled values"]],
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
      layers, num(set$SIDE), num(set$N_ROWS), num(set$PARENTS), num(set$SPREAD),
      num(set$RANGE), num(set$NOISE_VAR), num(set$N_BLOCK), num(set$N_AMONG), num(layer),
      cells, num(set$SIZE),
      margin, f_test,                  # what `structured` asks
      layers, nb, "1"),                # the white-noise layers; (1 - w) * tau2
    "?summarize_by_cell: Intervals: which one answers which question" = c(
      "1", "1",                        # deff = 1, and a number other than 1
      "1", "2", "1",                   # p = (1 - obs_level) / 2; quantile(e, 1 - p)
      level, level, level),            # the level checked; obs_level and conf_level in the run
    "?summarize_by_cell: Conditions" = c(
      pairs, "1", "0", "0",            # the limit; own weights of 1 and of 0; tau2 is 0
      few, num(.POOL_MIN_DF), level, needed),
    "?summarize_by_cell: Confidence intervals, the two intervals" = level,
    "?summarize_by_cell: argument estimate" = c(
      "1", nb, num(.POOL_NEIGHBOURS + 2L)),        # fewer than 10 cells: all the others
    "?summarize_by_cell: argument obs_level" = c(
      "0", "1", level,                 # a level in (0, 1), such as 0.90
      "1", "2", "1", "2",              # the (1 - obs_level) / 2 and (1 + obs_level) / 2 quantiles
      "0", "1", level),                # strictly between 0 and 1; only 0.90 checked
    "?resolution_profile: Observations per cell" = c(
      "2",                             # max(2, round(n / s))
      layers, num(layer), num(set$RANGE), num(set$SIDE), sizes,
      size, "1", "0",                  # 4.5 unless given; an own weight of 1, and of 0
      size, sizes,                     # what "auto" is
      sci(.FIT_WORK_MAX), pairs, dense, dense),
    "?resolution_profile: The ladder and its bounds, min_cell_n" =
      c(size, num(floor(.OBS_PER_CELL_DEFAULT))),  # "4.5 is read as 4"
    "?resolution_profile: The ladder and its bounds, a ladder of sizes" = "2",
    "?resolution_profile: The ladder and its bounds, no variogram" = character(0),
    "?resolution_profile: The criteria, loo" = character(0),
    "?select_resolution: What each criterion returns, loo" = f_test,
    "?select_resolution: What each criterion returns, size" = c("2", f_test),
    "?select_resolution: What each criterion returns, auto" = character(0),
    "?select_resolution: A profile built with an estimate" = character(0),
    "?cell_values_at: Away from the observations" = c(
      layers, num(layer), cells, num(set$N_AMONG), num(set$N_BLOCK), level),
    "?spatialkit: the two choices" = c("1", "4"),  # the steps of the pipeline
    "?spatialkit: the conventions of pooled values" = c(
      size, nb, num(.FEW_RESIDUALS), num(.POOL_MIN_DF), margin, sci(.EXACT_PAIRS_MAX),
      sci(.FIT_WORK_MAX)),
    "vignette resolution: Setting a size" = c(
      nb, "2", "1", "0",               # the 8 cells; max(2, ); an own weight of 1, and of 0
      layers, num(layer), num(max(set$SIZES)), num(min(set$SIZES))),
    "vignette resolution: Leave-one-out for plain means" = sizes,
    "vignette reporting: Reporting a cell value" = c(
      nb, "600", "20",                 # the section's own layer, as its chunk draws it
      size, level),
    "vignette diagnostics: 6. Warnings about pooled cell values" = c(
      size, margin, f_test, layers, nb,
      "1", "0", few))                  # own weights of 1 and of 0; 2 to 9 cells
  news <- stats::setNames(rep(list(character(0)), length(.DOCS_NEWS)),
                          paste0("NEWS: ", .DOCS_NEWS))
  news[[1L]]  <- c(nb, layers, num(layer), num(set$SIZE))
  news[[2L]]  <- c(level, level)
  news[[5L]]  <- num(.NSTART_DENSITY)              # the density kind's own starts
  news[[6L]]  <- size
  news[[12L]] <- "1"                               # deff = 1
  expected <- c(expected, news)

  expect_setequal(names(expected), .DOCS_WRITTEN)
  for (what in names(expected)) {
    expect_true(what %in% names(texts), info = what)
    rest <- .docs_rest(texts[[what]], claims, what)
    expect_identical(sort(.docs_numbers(rest)), sort(unname(expected[[what]])), info = what)
  }
  # A figure stands only in a text whose numbers are held.
  expect_true(all(vapply(claims, function(cl) cl$where, character(1)) %in% .DOCS_WRITTEN))
  # The layer of the reporting vignette's section is the one its chunk draws.
  code <- .docs_vignette_texts()[["vignette reporting: Reporting a cell value"]]
  expect_identical(.docs_count(code, "sx <- runif(600, 0, 2e4)"), 1L)
})

test_that("a text that says how a number tends to behave quotes a figure of the simulation", {
  sim <- .docs_sim()
  if (is.null(sim)) skip("dev/sim-pooled-values.csv is in the source tree only")
  texts  <- .docs_all_texts()
  quoted <- unique(vapply(.docs_claims(sim), function(cl) cl$where, character(1)))
  hedges <- lapply(stats::setNames(nm = .DOCS_WRITTEN), function(what)
    tolower(.docs_hedges(texts[[what]])))
  hedges <- hedges[lengths(hedges) > 0L]
  for (what in names(hedges))
    expect_true(what %in% quoted,
                info = sprintf("%s says '%s' and quotes no figure", what,
                               paste(hedges[[what]], collapse = "', '")))
  # The statements of that kind these texts make, by their words: the
  # pooled values' error falling to the finest size (output f), a pooled
  # value being mostly its neighbours' mean (output d), and an empty cell
  # coming more often from random starts.  One more is one more claim, to be
  # written here beside the output that shows it.
  expect_identical(hedges, list(
    "?summarize_by_cell: Pooled values"                          = "mostly",
    "?resolution_profile: Observations per cell"                 = c("tends", "tends", "often"),
    "?select_resolution: What each criterion returns, loo"       = "tends",
    "vignette resolution: Setting a size"                        = c("tends", "tends"),
    "vignette diagnostics: 6. Warnings about pooled cell values" = "mostly"))
})
