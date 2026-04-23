# =============================================================================
# scripts/appendix/density_ratio.R -- Density ratio estimation MSE
#
# Estimates marginal and conditional density ratios using LL and KLR on
# Diamonds (low-dim, 6 features) and Superconductivity (high-dim, ~80 features)
# datasets, then computes MSE relative to oracle ratios.
#
# Oracle ratios:
#   Marginal: g(x) = P_uniform(x) / P_biased(x), where biased ~ dnorm(x1,0,1)
#   Conditional: v(y) = P_null(y) / P_alt(y), where alt uses exp(-y) weighting
#
# Estimators: LL (linear logistic), KLR (kernel logistic regression)
# n in {200, 500, 1000, 2000}, n_rep=500
#
# Output: results/appendix/density_ratio/density_ratio_mse.csv
# Columns: dataset, estimator, n, hypothesis, mse_g12, mse_g22, mse_v12, mse_v22
#
# Usage (from package root):
#   Rscript scripts/appendix/density_ratio.R
# =============================================================================

source("R/load_all.R")
set.seed(1203)

# =============================================================================
# 1. Parameters
# =============================================================================

n_values   <- c(200, 500, 1000, 2000)
n_rep      <- 500
seed_base  <- 1203
estimators <- c("LL", "KLR")

# =============================================================================
# 2. Data Loading
# =============================================================================

# ---- Diamonds ---------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))
data("diamonds")
diamonds_X_raw <- as.matrix(diamonds[, c("carat", "depth", "table", "x", "y", "z")])
colnames(diamonds_X_raw) <- paste0("V", 1:6)
diamonds_Y_raw <- diamonds$price

normalize <- function(v) (v - min(v)) / (max(v) - min(v))
diamonds_X <- apply(diamonds_X_raw, 2, normalize)
diamonds_Y <- normalize(diamonds_Y_raw)

# ---- Superconductivity -----------------------------------------------------
sc_candidates <- c(
    file.path(RESULTS_DIR, "..", "data", "superconductivity.csv"),
    "data/superconductivity.csv",
    "../data/superconductivity.csv"
)
sc_path <- NULL
for (cand in sc_candidates) {
    if (file.exists(cand)) { sc_path <- cand; break }
}

has_superconductivity <- FALSE
if (!is.null(sc_path)) {
    sc_dat <- fread(sc_path)
    Y_col  <- ncol(sc_dat)
    sc_X_raw <- as.matrix(sc_dat[, seq_len(Y_col - 1), with = FALSE])
    sc_Y_raw <- sc_dat[[Y_col]]
    sc_X <- apply(sc_X_raw, 2, normalize)
    sc_Y <- normalize(sc_Y_raw)
    # Remove zero-variance columns
    col_var <- apply(sc_X, 2, var)
    sc_X <- sc_X[, col_var > 1e-10, drop = FALSE]
    has_superconductivity <- TRUE
    cat(sprintf("[Data] Superconductivity loaded: %d obs x %d features\n",
                nrow(sc_X), ncol(sc_X)))
} else {
    warning("Superconductivity CSV not found. Skipping.")
}

# =============================================================================
# 3. Oracle Density Ratios
# =============================================================================

# MSE helper
mse <- function(true, est) mean((true - est)^2)

# True marginal density ratio: P_uniform / P_biased
# P_uniform is constant, P_biased proportional to dnorm(x1, 0, 1)
true_marginal_density_ratio <- function(X1, X2, x_subset, is_x1 = TRUE) {
    prob_X1 <- rep(1 / nrow(X1), nrow(X1))
    feature_to_bias <- X2[, 1]
    prob_X2 <- dnorm(feature_to_bias, mean = 0, sd = 1)
    prob_X2 <- prob_X2 / sum(prob_X2)

    if (is_x1) {
        indices <- match(x_subset[, 1], X1[, 1])
        ratio <- prob_X1[indices] / prob_X2[indices]
    } else {
        indices <- match(x_subset[, 1], X2[, 1])
        ratio <- prob_X1[indices] / prob_X2[indices]
    }
    ratio
}

# True conditional density ratio: P_null(y) / P_alt(y)
true_conditional_density_ratio <- function(Y_subset, is_null = TRUE, is_x1 = TRUE) {
    if (is_null) {
        return(rep(1, length(Y_subset)))
    }
    prob_Y_given_X1 <- dunif(Y_subset, min = 0, max = 1)
    prob_Y_given_X1 <- prob_Y_given_X1 / sum(prob_Y_given_X1)
    prob_Y_given_X2 <- exp(-Y_subset)
    prob_Y_given_X2 <- prob_Y_given_X2 / sum(prob_Y_given_X2)
    prob_Y_given_X1 / prob_Y_given_X2
}

# =============================================================================
# 4. Sampling Function (returns extra info for oracle computation)
# =============================================================================

sample_data_oracle <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
    # Split entire dataset in half for the two groups
    X_idx <- sample(seq_len(nrow(X)), nrow(X) %/% 2, replace = FALSE)
    X1 <- X[X_idx, , drop = FALSE]
    Y1 <- Y[X_idx]

    feature_to_bias <- X[, 1]
    prob <- dnorm(feature_to_bias, mean = 0, sd = 1)
    prob <- prob / sum(prob)
    X2_idx <- sample(seq_len(nrow(X)), nrow(X) %/% 2, replace = FALSE, prob = prob)
    X2 <- X[X2_idx, , drop = FALSE]
    Y2 <- Y[X2_idx]

    if (is_x1) {
        x_idx <- sample(seq_len(nrow(X1)), n, replace = FALSE)
        x <- X1[x_idx, , drop = FALSE]
        Y_subset <- Y1[x_idx]
    } else {
        x_idx <- sample(seq_len(nrow(X2)), n, replace = FALSE)
        x <- X2[x_idx, , drop = FALSE]
        Y_subset <- Y2[x_idx]
    }

    if (is_null) {
        y <- sample(Y_subset, size = n, replace = FALSE)
    } else {
        u <- if (is_x1) dunif(Y_subset, min = 0, max = 1) else exp(-Y_subset)
        u <- u / sum(u)
        y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
    }
    list(x = x, y = y, X1 = X1, X2 = X2, y_subset = Y_subset)
}

# =============================================================================
# 5. Single Simulation Run
# =============================================================================

run_simulation <- function(X, Y, n, is_null, estimator, seed) {
    set.seed(seed)
    d1 <- sample_data_oracle(X, Y, n, is_null, TRUE)
    set.seed(seed + 500)
    d2 <- sample_data_oracle(X, Y, n, is_null, FALSE)

    # Split each sample in half for estimation / evaluation
    split_idx <- sample(seq_len(n), n %/% 2)
    x11 <- d1$x[split_idx, , drop = FALSE];  y11 <- d1$y[split_idx]
    x12 <- d1$x[-split_idx, , drop = FALSE]; y12 <- d1$y[-split_idx]
    x21 <- d2$x[split_idx, , drop = FALSE];  y21 <- d2$y[split_idx]
    x22 <- d2$x[-split_idx, , drop = FALSE]; y22 <- d2$y[-split_idx]

    # Oracle ratios
    true_g12 <- true_marginal_density_ratio(d1$X1, d2$X2, x12, is_x1 = TRUE)
    true_g22 <- true_marginal_density_ratio(d1$X1, d2$X2, x22, is_x1 = FALSE)
    true_v12 <- true_conditional_density_ratio(d1$y_subset, is_null, is_x1 = TRUE)
    true_v22 <- true_conditional_density_ratio(d2$y_subset, is_null, is_x1 = FALSE)

    # Estimated ratios
    est <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22,
                      est.method = estimator, seed = seed)

    c(mse_g12 = mse(true_g12, est$g12.est),
      mse_g22 = mse(true_g22, est$g22.est),
      mse_v12 = mse(true_v12, est$v12.est),
      mse_v22 = mse(true_v22, est$v22.est))
}

# =============================================================================
# 6. Main Experiment
# =============================================================================

run_mse_experiment <- function(X_norm, Y_norm, dataset_name) {
    cat(sprintf("\n########## Dataset: %s ##########\n", dataset_name))
    local_results <- list()

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (est in estimators) {
                cat(sprintf("[%s] %s | n=%d | %s ...", dataset_name, est, n, h_label))

                # sapply (not cdtst_sapply) because we need a 4 x n_rep matrix
                # and cdtst_sapply's foreach path uses .combine=c which flattens.
                result <- sapply(seq_len(n_rep), function(sim) {
                    seed <- seed_base + sim
                    tryCatch(
                        run_simulation(X_norm, Y_norm, n, is_null, est, seed),
                        error = function(e) rep(NA_real_, 4)
                    )
                })

                # Median across replications (robust to outliers)
                median_result <- apply(result, 1, median, na.rm = TRUE)

                local_results[[length(local_results) + 1]] <- data.table(
                    dataset   = dataset_name,
                    estimator = est,
                    n         = n,
                    hypothesis = h_label,
                    mse_g12   = median_result["mse_g12"],
                    mse_g22   = median_result["mse_g22"],
                    mse_v12   = median_result["mse_v12"],
                    mse_v22   = median_result["mse_v22"]
                )
                cat(sprintf(" g12=%.4f g22=%.4f v12=%.4f v22=%.4f\n",
                            median_result["mse_g12"], median_result["mse_g22"],
                            median_result["mse_v12"], median_result["mse_v22"]))
            }
        }
    }
    rbindlist(local_results)
}

# =============================================================================
# 7. Execute and Save
# =============================================================================

all_results <- list()

all_results[[1]] <- run_mse_experiment(diamonds_X, diamonds_Y, "diamonds")

if (has_superconductivity) {
    all_results[[2]] <- run_mse_experiment(sc_X, sc_Y, "superconductivity")
}

results_dt <- rbindlist(all_results)
out_file <- file.path(RESULTS_DIR, "appendix/density_ratio/density_ratio_mse.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary ---
cat("\n========== Density Ratio MSE Summary ==========\n")
print(results_dt)
