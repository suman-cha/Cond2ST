# =============================================================================
# scripts/main/real_data.R -- Main-track real data experiments
#
# Runs the 6 DRT tests that appear in the paper's real-data figure:
#   CV_LinearMMD, CV_BlockMMD, CV_QuadraticMMD, CV_CLF (learned+ranger),
#   CP, DCP
# across Diamonds and Superconductivity, LL + KLR estimators.
#
# Appendix-only tests (plain LinearMMD/BlockMMD/QuadraticMMD/CLF plus CIT
# and CK) live in scripts/appendix/real_data.R.
#
# Sampling:
#   Group 1: uniform subsample of half the dataset
#   Group 2: biased subsample weighted by dnorm(X[,1], 0, 1)
#   Null: Y uniformly sampled from Y_subset
#   Alt : Y weighted by dunif(Y,0,1) for group 1, exp(-Y) for group 2
#
# n in {200, 400, 800, 1200, 1600, 2000}, n_rep=500, alpha=0.05
#
# Output: results/main/real_data/simulation_results_{dataset}_drt.csv
#   schema: dataset, test_name, estimator, n, hypothesis, alpha, n_rep,
#           rejection_rate
#
# Env overrides: CDTST_N_REP, CDTST_N_VALUES, CDTST_DATASETS, CDTST_CORES
#
# Usage (from package root):
#   CDTST_CORES=8 Rscript scripts/main/real_data.R
# =============================================================================

source("R/load_all.R")
set.seed(1203)

# ---- Parameters -------------------------------------------------------------
n_values  <- c(200, 400, 800, 1200, 1600, 2000)
n_rep     <- as.integer(Sys.getenv("CDTST_N_REP", "500"))
alpha     <- 0.05
seed_base <- 1203
dataset_names <- c("diamonds", "superconductivity")

.env_n_values <- Sys.getenv("CDTST_N_VALUES", "")
if (nzchar(.env_n_values)) {
    n_values <- as.integer(strsplit(.env_n_values, ",", fixed = TRUE)[[1]])
}
.env_datasets <- Sys.getenv("CDTST_DATASETS", "")
if (nzchar(.env_datasets)) {
    dataset_names <- strsplit(.env_datasets, ",", fixed = TRUE)[[1]]
}
rm(.env_n_values, .env_datasets)

# ---- Data loading -----------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
data("diamonds")
diamonds_X_raw <- as.matrix(diamonds[, c("carat", "depth", "table", "x", "y", "z")])
colnames(diamonds_X_raw) <- paste0("V", 1:6)
normalize <- function(v) (v - min(v)) / (max(v) - min(v))
diamonds_X <- apply(diamonds_X_raw, 2, normalize)
diamonds_Y <- normalize(diamonds$price)

sc_candidates <- c(
    file.path(RESULTS_DIR, "..", "data", "superconductivity.csv"),
    "data/superconductivity.csv",
    "../data/superconductivity.csv"
)
sc_path <- NULL
for (cand in sc_candidates) if (file.exists(cand)) { sc_path <- cand; break }

has_superconductivity <- FALSE
if (!is.null(sc_path)) {
    sc_dat   <- fread(sc_path)
    Y_col    <- ncol(sc_dat)
    sc_X_raw <- as.matrix(sc_dat[, seq_len(Y_col - 1), with = FALSE])
    sc_X     <- apply(sc_X_raw, 2, normalize)
    sc_Y     <- normalize(sc_dat[[Y_col]])
    sc_X     <- sc_X[, apply(sc_X, 2, var) > 1e-10, drop = FALSE]
    has_superconductivity <- TRUE
    cat(sprintf("[Data] Superconductivity: %d obs x %d features\n",
                nrow(sc_X), ncol(sc_X)))
} else {
    warning("Superconductivity CSV not found; skipping that dataset.")
}

# ---- Sampling (shared for both datasets) ------------------------------------
sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
    if (is_x1) {
        X_idx    <- sample(seq_len(nrow(X)), nrow(X) %/% 2, replace = FALSE)
        X_sub    <- X[X_idx, , drop = FALSE]
        Y_subset <- Y[X_idx]
    } else {
        prob  <- dnorm(X[, 1], 0, 1); prob <- prob / sum(prob)
        X_idx <- sample(seq_len(nrow(X)), nrow(X) %/% 2, replace = FALSE, prob = prob)
        X_sub    <- X[X_idx, , drop = FALSE]
        Y_subset <- Y[X_idx]
    }
    x_idx <- sample(seq_len(nrow(X_sub)), n, replace = FALSE)
    x     <- X_sub[x_idx, , drop = FALSE]
    if (is_null) {
        y <- sample(Y_subset, size = n, replace = FALSE)
    } else {
        u <- if (is_x1) dunif(Y_subset, 0, 1) else exp(-Y_subset)
        u <- u / sum(u)
        y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
    }
    list(x = x, y = y)
}

# ---- Main-track DRT tests ---------------------------------------------------
drt_tests <- list(
    CV_LinearMMD_test    = CV_LinearMMD_test,
    CV_BlockMMD_test     = CV_BlockMMD_test,
    CV_QuadraticMMD_test = CV_QuadraticMMD_test,
    CV_CLF_test          = CV_CLF_test,          # learned + ranger
    CP_test              = CP_test,
    DCP_test             = DCP_test
)
drt_estimators <- c("LL", "KLR")

# ---- Runner -----------------------------------------------------------------
run_dataset <- function(X_norm, Y_norm, dataset_name) {
    seed_base   <- seed_base
    n_rep       <- n_rep
    alpha       <- alpha
    sample_data <- sample_data
    drt_tests   <- drt_tests

    cat(sprintf("\n########## Dataset: %s ##########\n", dataset_name))
    out_dir <- file.path(RESULTS_DIR, "main", "real_data")
    ensure_dir(out_dir)
    results <- list()

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (test_name in names(drt_tests)) {
                for (est in drt_estimators) {
                    cat(sprintf("[%s][DRT] %s|%s | n=%d | %s\n",
                                dataset_name, test_name, est, n, h_label))

                    res <- cdtst_sapply(n_rep, function(sim) {
                        seed <- seed_base + sim
                        set.seed(seed)
                        d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
                        set.seed(seed + n_rep)
                        d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
                        tm <- system.time({
                            val <- tryCatch(
                                drt_tests[[test_name]](d1$x, d2$x, d1$y, d2$y,
                                                       est.method = est,
                                                       alpha = alpha, seed = seed),
                                error = function(e) NA_real_
                            )
                        })
                        c(val, unname(tm["elapsed"]))
                    })

                    if (is.matrix(res)) {
                        vals <- res[1, ]; times <- res[2, ]
                    } else {
                        vals <- res[seq(1, length(res), by = 2)]
                        times <- res[seq(2, length(res), by = 2)]
                    }
                    rr <- mean(vals, na.rm = TRUE)
                    cat(sprintf("  RR=%.3f  time=%.2f±%.2fs\n",
                                rr, mean(times, na.rm = TRUE),
                                sd(times, na.rm = TRUE)))

                    results[[length(results) + 1]] <- data.table(
                        dataset = dataset_name, test_name = test_name,
                        estimator = est, n = n, hypothesis = h_label,
                        alpha = alpha, n_rep = n_rep, rejection_rate = rr,
                        time_mean = mean(times, na.rm = TRUE),
                        time_sd = sd(times, na.rm = TRUE)
                    )
                }
            }
        }
    }

    dt   <- rbindlist(results)
    file <- file.path(out_dir,
                       sprintf("simulation_results_%s_drt.csv", dataset_name))
    fwrite(dt, file)
    cat(sprintf("Saved: %s (%d rows)\n", file, nrow(dt)))
}

# ---- Execute ----------------------------------------------------------------
if ("diamonds" %in% dataset_names) {
    run_dataset(diamonds_X, diamonds_Y, "diamonds")
}
if ("superconductivity" %in% dataset_names && has_superconductivity) {
    run_dataset(sc_X, sc_Y, "superconductivity")
}

cat("\n========== Main Real Data (DRT) Complete ==========\n")
