#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Does GWmodel's gwr.model.selection() really report AICc in column 3?
#
# The whole third-pass GWR fix rests on this. GWmodel's own documentation and
# its Model.selection.r source say GWR.df is
#
#     c(bandwidth, AIC, AICc, RSS)
#
# built by rbind() over UNNAMED vectors -- so the table never carries column
# names and the positional read is the path every real call takes. Column 2 is
# the UNCORRECTED AIC. Reading it ranks on AIC while labelling the answer AICc,
# which selects larger models and can keep a pure-noise predictor.
#
# That was verified against a faithful API stub, not the real package. This
# script settles it on a real GWmodel install. It does NOT trust the column
# order: it recomputes AICc independently for every candidate model with
# gwr.basic() and asks which column of GWR.df those numbers actually match.
#
#   Rscript dev/verify-gwr-aicc.R
# ---------------------------------------------------------------------------

for (p in c("sf", "sp", "GWmodel")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("this script needs '", p, "' installed.", call. = FALSE)
}
suppressPackageStartupMessages({library(sf); library(sp); library(GWmodel)})

cat("GWmodel version:", as.character(utils::packageVersion("GWmodel")), "\n\n")

set.seed(7)
n  <- 200
xy <- data.frame(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                 a = rnorm(n), b = rnorm(n), noise = rnorm(n))
xy$y_resp <- 5 + 2 * xy$a - 1.5 * xy$b + rnorm(n, 0, 1)   # `noise` is pure noise
spd <- SpatialPointsDataFrame(coords = as.matrix(xy[, c("x", "y")]),
                              data   = xy[, c("a", "b", "noise", "y_resp")])

bw  <- 60
sel <- GWmodel::gwr.model.selection(DeVar = "y_resp",
                                    InDeVars = c("a", "b", "noise"),
                                    data = spd, bw = bw,
                                    adaptive = TRUE, kernel = "bisquare")
model_list <- sel[[1]]
gwr_df     <- sel[[2]]

cat("GWR.df dim:      ", paste(dim(gwr_df), collapse = " x "), "\n")
cat("GWR.df colnames: ",
    if (is.null(colnames(gwr_df))) "<none -- positional read is the only option>"
    else paste(colnames(gwr_df), collapse = ", "), "\n\n")

# Independently recompute AIC and AICc for each candidate with gwr.basic().
vars_of <- function(m) as.character(unlist(m[[2]], use.names = FALSE))
ref <- t(vapply(model_list, function(m) {
  v   <- vars_of(m)
  fml <- stats::reformulate(termlabels = v, response = "y_resp")
  g   <- GWmodel::gwr.basic(fml, data = spd, bw = bw,
                            adaptive = TRUE, kernel = "bisquare")
  d <- g$GW.diagnostic
  c(AIC = as.numeric(d$AIC), AICc = as.numeric(d$AICc))
}, c(AIC = 0, AICc = 0)))

tab <- data.frame(
  model = vapply(model_list, function(m) paste(vars_of(m), collapse = " + "),
                 character(1)),
  col2  = as.numeric(gwr_df[, 2]),
  col3  = as.numeric(gwr_df[, 3]),
  ref_AIC = ref[, "AIC"], ref_AICc = ref[, "AICc"]
)
print(tab, row.names = FALSE, digits = 7)

close_to <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-6))
cat("\ncolumn 2 == independently computed AIC :", close_to(tab$col2, tab$ref_AIC), "\n")
cat("column 3 == independently computed AICc:", close_to(tab$col3, tab$ref_AICc), "\n")

pick <- function(v) tab$model[which.min(v)]
cat("\nselected by column 2 (AIC) :", pick(tab$col2), "\n")
cat("selected by column 3 (AICc):", pick(tab$col3), "\n")

ok <- close_to(tab$col3, tab$ref_AICc) && !close_to(tab$col2, tab$ref_AICc)
cat("\n", if (ok)
  "PASS: AICc is column 3, and column 2 is not AICc. The fix reads the right column."
  else
  "FAIL: this GWmodel does NOT lay GWR.df out as c(bandwidth, AIC, AICc, RSS). Do not submit; tell Claude.",
  "\n", sep = "")

if (grepl("noise", pick(tab$col2)) && !grepl("noise", pick(tab$col3)))
  cat("Confirmed on this data: AIC keeps the pure-noise predictor, AICc drops it.\n")
