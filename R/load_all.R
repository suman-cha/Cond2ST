# R/load_all.R
# Single entry point for all experiments.
# Usage: source("R/load_all.R")

# ── Locate R/ directory ───────────────────────────────────────────────────────
# All scripts must be run from the package root as the working directory.
.lib_dir <- "R"
if (!file.exists(file.path(.lib_dir, "regression.R"))) {
    stop("Cannot find R/regression.R. ",
         "Please set working directory to the package root before running.")
}

# ── Load required packages (suppress startup messages) ────────────────────────
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(ranger)
    library(xgboost)
    library(glmnet)
    library(kernlab)
    library(CVST)
    library(densratio)
    library(RCIT)
    library(CondIndTests)
    library(vimp)
})

# ── Source library modules in dependency order ────────────────────────────────
source(file.path(.lib_dir, "regression.R"))
source(file.path(.lib_dir, "density_ratio.R"))
source(file.path(.lib_dir, "cp.R"))
source(file.path(.lib_dir, "mmd.R"))
source(file.path(.lib_dir, "subsampling.R"))
source(file.path(.lib_dir, "cit.R"))
source(file.path(.lib_dir, "hypothesis_tests.R"))
source(file.path(.lib_dir, "classifier.R"))
source(file.path(.lib_dir, "scenarios.R"))

# ── MMD backend check ─────────────────────────────────────────────────────────
# `.mmd_use_cpp` is set in R/mmd.R at source time. The C++ kernels are
# >10x faster than the R fallback at n>=300, so surfacing which path is live
# makes performance regressions (e.g., missing src/mmd_kernels.cpp or
# Rcpp compile failure) visible rather than silent.
if (isTRUE(.mmd_use_cpp)) {
    message("[mmd] C++ backend active (src/mmd_kernels.cpp)")
} else {
    warning("[mmd] Using pure-R MMD fallback; experiments will be >10x slower. ",
            "Check that src/mmd_kernels.cpp is present and Rcpp compiled.",
            call. = FALSE, immediate. = TRUE)
}

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Create a directory recursively if it does not exist.
#'
#' @param ... Path components passed to \code{file.path}.
#' @return The (possibly created) path, returned invisibly.
ensure_dir <- function(...) {
    path <- file.path(...)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
}

RESULTS_DIR <- "results"

# ── Parallel backend (opt-in) ─────────────────────────────────────────────────
# Opt in by setting `options(cdtst.cores = N)` before source-ing this file, or
# via the `CDTST_CORES` env var. Each rep body already calls
#   gen_pair(..., seed = seed_base + sim)
# and every test starts with `if (!is.null(seed)) set.seed(seed)`, so
# per-worker RNG state is fully reset before any RNG consumption that affects
# output — parallelism is expected to be bit-exact with the serial path.
.cdtst_cluster <- NULL

#' Register a SOCK cluster for parallel Monte-Carlo replication.
#'
#' Reads the worker count from \code{getOption("cdtst.cores")} or the
#' \code{CDTST_CORES} environment variable. Each worker re-runs
#' \code{source("R/load_all.R")} so that Rcpp-compiled MMD kernels are
#' linked into the worker process. Idempotent: re-calling returns the
#' existing cluster.
#'
#' @return The registered cluster handle (invisibly), or \code{NULL} if no
#'   parallelism was requested.
cdtst_parallel_setup <- function() {
    if (nzchar(Sys.getenv("CDTST_PARALLEL_CHILD", ""))) return(invisible(NULL))
    cores <- getOption("cdtst.cores",
                       as.integer(Sys.getenv("CDTST_CORES", "1")))
    if (is.null(cores) || is.na(cores) || cores <= 1L) return(invisible(NULL))
    if (inherits(.cdtst_cluster, "cluster")) return(invisible(.cdtst_cluster))
    suppressPackageStartupMessages({
        requireNamespace("parallel")
        requireNamespace("doSNOW")
        requireNamespace("foreach")
    })
    cl <- parallel::makeCluster(cores, type = "SOCK")
    wd <- getwd()
    parallel::clusterExport(cl, "wd", envir = environment())
    parallel::clusterEvalQ(cl, {
        Sys.setenv(CDTST_PARALLEL_CHILD = "1")
        setwd(wd)
        source("R/load_all.R")
    })
    doSNOW::registerDoSNOW(cl)
    .cdtst_cluster <<- cl
    message(sprintf("[parallel] registered cluster with %d cores", cores))
    invisible(cl)
}

# Drop-in replacement for `sapply(seq_len(n_rep), fn)`. Uses foreach when a
# cluster has been registered; otherwise falls back to serial sapply so the
# serial run produces bit-identical output.
#
# foreach's auto-export can't always traverse through the `fn` argument's
# closure when the runner is launched via `Rscript` (closure env == globalenv)
# — explicitly export names visible in fn's environment so workers have what
# they need (seed_base, sc, tc, eps, alpha, etc.).
#' Drop-in parallel \code{sapply} for Monte-Carlo replication.
#'
#' Behaves like \code{sapply(seq_len(n_rep), fn)}: returns a vector
#' (typically of 0/1 rejections) of length \code{n_rep}. If a SOCK cluster
#' has been registered via \code{cdtst_parallel_setup}, distributes the
#' replications via \code{foreach::\%dopar\%}, otherwise falls back to
#' serial \code{sapply} so the serial path is bit-identical.
#'
#' @param n_rep Number of replications.
#' @param fn A function \code{function(sim) ...} where \code{sim} is the
#'   1-based replication index. Closure variables visible to \code{fn} are
#'   auto-exported to workers.
#' @return Numeric / integer vector of length \code{n_rep}.
cdtst_sapply <- function(n_rep, fn) {
    if (inherits(.cdtst_cluster, "cluster")) {
        fn_env <- environment(fn)
        # Export every binding visible to fn's defining scope directly to
        # workers. `foreach`'s .export looks up names in the caller of
        # `%dopar%` (here: cdtst_sapply's frame), which doesn't contain the
        # runner-local vars — so we call clusterExport ourselves with the
        # correct envir. This also covers transitive closures like gen_pair
        # referencing n_rep.
        if (is.environment(fn_env)) {
            # all.names = TRUE includes dot-prefixed vars (runners sometimes use
            # them for per-cell state), at the cost of picking up `.Random.seed`
            # etc. — filtered out below.
            export_names <- ls(fn_env, all.names = TRUE)
            # Remove load_all.R-defined helpers — workers already have their
            # own (correctly-linked for Rcpp) copies — plus R-internal dot-vars.
            if (exists(".cdtst_shared_globals", envir = globalenv(), inherits = FALSE)) {
                export_names <- setdiff(export_names,
                                        get(".cdtst_shared_globals", envir = globalenv()))
            }
            export_names <- setdiff(export_names,
                                    c(".Random.seed", ".cdtst_cluster"))
            if (length(export_names) > 0) {
                parallel::clusterExport(.cdtst_cluster, export_names, envir = fn_env)
            }
        }
        foreach::foreach(sim = seq_len(n_rep), .combine = c,
                         .packages = c("data.table")) %dopar% fn(sim)
    } else {
        sapply(seq_len(n_rep), fn)
    }
}

library(foreach)   # makes `%dopar%` available to runners

# Snapshot the set of names defined by load_all.R itself. When the runner is
# launched via `Rscript` directly, fn_env == globalenv(), so `ls(fn_env)` would
# include these library helpers — shipping them to workers can BREAK things
# (e.g. Rcpp-compiled wrappers like mmd_linear_cpp carry DLL-specific symbol
# pointers that are invalid on a different process). Workers already have
# correct copies via clusterEvalQ(source("R/load_all.R")), so we filter
# these out of the export set in cdtst_sapply.
.cdtst_shared_globals <- ls(globalenv(), all.names = FALSE)

cdtst_parallel_setup()

rm(.lib_dir)
