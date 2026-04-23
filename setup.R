# setup.R — install all R packages required for Cond2ST experiments.
#
# Usage (from package root):
#   Rscript setup.R
#
# For exact-version reproducibility, prefer:
#   R> renv::restore()
# which uses ./renv.lock to pin every transitive dependency.

cat("=== Cond2ST: package setup ===\n\n")

# ---- CRAN -------------------------------------------------------------------
cran_pkgs <- c(
  # simulation
  "MASS", "pbapply", "data.table", "tmvtnorm",
  # ML / regression
  "ranger", "xgboost", "glmnet", "e1071", "kernlab",
  "caret", "mgcv", "nnet", "recipes", "SuperLearner",
  # parallelisation
  "foreach", "doSNOW", "doRNG",
  # auxiliary
  "progress", "prettyunits", "plyr", "CVST", "ANN2", "stabs", "devtools",
  "densratio", "Rcpp",
  # CIT / DRT baselines (current CRAN versions)
  "CondIndTests", "vimp",
  # visualisation
  "ggplot2", "dplyr", "tidyr", "stringr", "scales",
  "grid", "gridExtra", "latex2exp", "ggrepel", "forcats", "purrr", "readr",
  # utilities
  "withr", "R.utils", "microbenchmark"
)

to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  cat("Installing CRAN packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  cat("All CRAN packages already installed.\n")
}

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

# ---- GitHub-only ------------------------------------------------------------
github_pkgs <- list(
  list(pkg = "RCIT",  repo = "ericstrobl/RCIT"),
  list(pkg = "KDist", repo = "zhangxiany-tamu/KDist")
)

for (p in github_pkgs) {
  if (!requireNamespace(p$pkg, quietly = TRUE)) {
    cat("Installing from GitHub:", p$repo, "\n")
    remotes::install_github(p$repo)
  } else {
    cat("Already installed:", p$pkg, "\n")
  }
}

# ---- CRAN-archived (no longer on the active CRAN index) --------------------
# GeneralisedCovarianceMeasure 0.2.0 was archived from CRAN; pin to that version
# from the CRAN archive. The GitHub mirror (runesen/gcm) is no longer reachable.
if (!requireNamespace("GeneralisedCovarianceMeasure", quietly = TRUE)) {
  archive_url <- paste0(
    "https://cran.r-project.org/src/contrib/Archive/",
    "GeneralisedCovarianceMeasure/GeneralisedCovarianceMeasure_0.2.0.tar.gz"
  )
  cat("Installing from CRAN archive: GeneralisedCovarianceMeasure 0.2.0\n")
  install.packages(archive_url, repos = NULL, type = "source")
} else {
  cat("Already installed: GeneralisedCovarianceMeasure\n")
}

cat("\nAll packages ready.\n")
