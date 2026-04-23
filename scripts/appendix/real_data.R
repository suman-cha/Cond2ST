# =============================================================================
# scripts/appendix/real_data.R -- Appendix-track real data experiments
#
# Runs the 10 tests that go into the appendix tables:
#   DRT (non-CV, 4): LinearMMD, BlockMMD, QuadraticMMD, CLF (learned+ranger)
#   CIT         (4): GCM (ranger), PCM (ranger), RCIT, WGSC (xgb) — Alg1, eps
#   CK          (2): CMMD, CGED
# across Diamonds and Superconductivity.
#
# Main-track counterpart is scripts/main/real_data.R (6 CV DRT tests only).
#
# Sampling, normalization, and hypothesis structure are identical to the
# main-track runner; this file is a standalone copy to keep the two tracks
# independently executable.
#
# Output: results/appendix/real_data/simulation_results_{dataset}_{drt,cit,ck}.csv
#
# Env overrides: CDTST_N_REP, CDTST_N_VALUES, CDTST_DATASETS, CDTST_CORES
#
# Usage (from package root):
#   CDTST_CORES=8 Rscript scripts/appendix/real_data.R
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

skip_drt <- nzchar(Sys.getenv("CDTST_SKIP_DRT", ""))
skip_cit <- nzchar(Sys.getenv("CDTST_SKIP_CIT", ""))
skip_ck  <- nzchar(Sys.getenv("CDTST_SKIP_CK",  ""))

# ---- Data loading (identical to main runner) --------------------------------
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

# ---- Sampling (identical to main runner) ------------------------------------
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

# ---- Appendix test configurations -------------------------------------------
drt_tests <- list(
    LinearMMD_test    = LinearMMD_test,
    BlockMMD_test     = BlockMMD_test,
    QuadraticMMD_test = QuadraticMMD_test,
    CLF_test          = CLF_test          # learned + ranger
)
drt_estimators <- c("LL", "KLR")

cit_configs <- list(
    list(name = "GCM_test",  fn = GCM_test,  reg_name = "rf",
         reg = ranger_reg_method, breg = ranger_reg_method_binary),
    list(name = "PCM_test",  fn = PCM_test,  reg_name = "rf",
         reg = ranger_reg_method, breg = ranger_reg_method_binary),
    list(name = "RCIT_test", fn = RCIT_test, reg_name = NA_character_,
         reg = NULL, breg = NULL),
    list(name = "WGSC_test", fn = WGSC_test, reg_name = "xgb",
         reg = xgboost_reg_method, breg = xgboost_reg_method_binary)
)

ck_tests <- list(CMMD_test = CMMD_test, CGED_test = CGED_test)

# ---- Runner -----------------------------------------------------------------
run_dataset <- function(X_norm, Y_norm, dataset_name) {
    seed_base   <- seed_base
    n_rep       <- n_rep
    alpha       <- alpha
    sample_data <- sample_data
    drt_tests   <- drt_tests
    cit_configs <- cit_configs
    ck_tests    <- ck_tests

    cat(sprintf("\n########## Dataset: %s ##########\n", dataset_name))
    out_dir <- file.path(RESULTS_DIR, "appendix", "real_data")
    ensure_dir(out_dir)

    # --- DRT ---
    if (!skip_drt) {
    cat(sprintf("\n--- [%s] DRT (appendix) ---\n", dataset_name))
    drt_results <- list()
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
                    drt_results[[length(drt_results) + 1]] <- data.table(
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
    drt_dt <- rbindlist(drt_results)
    drt_file <- file.path(out_dir,
                           sprintf("simulation_results_%s_drt.csv", dataset_name))
    fwrite(drt_dt, drt_file)
    cat(sprintf("Saved: %s (%d rows)\n", drt_file, nrow(drt_dt)))
    } else {
        cat(sprintf("\n--- [%s] DRT SKIPPED ---\n", dataset_name))
    }

    # --- CIT ---
    if (!skip_cit) {
    cat(sprintf("\n--- [%s] CIT (appendix) ---\n", dataset_name))
    cit_results <- list()
    for (n in n_values) {
        eps <- 1 / sqrt(log(n))
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (tc in cit_configs) {
                cat(sprintf("[%s][CIT] %s(%s) | n=%d | %s\n",
                            dataset_name, tc$name, tc$reg_name, n, h_label))
                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    set.seed(seed)
                    d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
                    set.seed(seed + n_rep)
                    d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
                    tm <- system.time({
                        val <- tryCatch({
                            args <- list(d1$x, d2$x, d1$y, d2$y,
                                         alg1 = TRUE, epsilon = eps,
                                         alpha = alpha, seed = seed)
                            if (!is.null(tc$reg)) {
                                args$regr.method <- tc$reg
                                args$binary.regr.method <- tc$breg
                            }
                            do.call(tc$fn, args)
                        }, error = function(e) NA_real_)
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
                cit_results[[length(cit_results) + 1]] <- data.table(
                    dataset = dataset_name, test_name = tc$name,
                    estimator = tc$reg_name, n = n, hypothesis = h_label,
                    alpha = alpha, n_rep = n_rep, rejection_rate = rr,
                    time_mean = mean(times, na.rm = TRUE),
                    time_sd = sd(times, na.rm = TRUE)
                )
            }
        }
    }
    cit_dt <- rbindlist(cit_results)
    cit_file <- file.path(out_dir,
                           sprintf("simulation_results_%s_cit.csv", dataset_name))
    fwrite(cit_dt, cit_file)
    cat(sprintf("Saved: %s (%d rows)\n", cit_file, nrow(cit_dt)))
    } else {
        cat(sprintf("\n--- [%s] CIT SKIPPED ---\n", dataset_name))
    }

    # --- CK ---
    if (!skip_ck) {
    cat(sprintf("\n--- [%s] CK (appendix) ---\n", dataset_name))
    ck_results <- list()
    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (test_name in names(ck_tests)) {
                cat(sprintf("[%s][CK] %s | n=%d | %s\n",
                            dataset_name, test_name, n, h_label))
                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    set.seed(seed)
                    d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
                    set.seed(seed + n_rep)
                    d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
                    tm <- system.time({
                        val <- tryCatch(
                            ck_tests[[test_name]](d1$x, d2$x, d1$y, d2$y,
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
                ck_results[[length(ck_results) + 1]] <- data.table(
                    dataset = dataset_name, test_name = test_name,
                    estimator = NA_character_, n = n, hypothesis = h_label,
                    alpha = alpha, n_rep = n_rep, rejection_rate = rr,
                    time_mean = mean(times, na.rm = TRUE),
                    time_sd = sd(times, na.rm = TRUE)
                )
            }
        }
    }
    ck_dt <- rbindlist(ck_results)
    ck_file <- file.path(out_dir,
                          sprintf("simulation_results_%s_ck.csv", dataset_name))
    fwrite(ck_dt, ck_file)
    cat(sprintf("Saved: %s (%d rows)\n", ck_file, nrow(ck_dt)))
    } else {
        cat(sprintf("\n--- [%s] CK SKIPPED ---\n", dataset_name))
    }
}

# ---- Execute ----------------------------------------------------------------
if ("diamonds" %in% dataset_names) {
    run_dataset(diamonds_X, diamonds_Y, "diamonds")
}
if ("superconductivity" %in% dataset_names && has_superconductivity) {
    run_dataset(sc_X, sc_Y, "superconductivity")
}

cat("\n========== Appendix Real Data Complete ==========\n")
