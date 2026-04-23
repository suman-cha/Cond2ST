# =============================================================================
# scripts/appendix/oracle_cit.R
#
# Oracle CIT vs C2ST (Algorithm 1) comparison.
#
# Purpose: Quantify the power loss induced by Algorithm 1's subsampling step
#          by comparing CIT performance in two settings:
#          1. C2ST setting: Two independent samples of size n each (total 2n),
#             apply Algorithm 1 (subsampling) before running CIT
#          2. Oracle CIT setting: 2n i.i.d. triples from P_{XYZ},
#             run CIT directly on all samples
#
# CIT Methods: GCM (ranger), PCM (xgboost), RCIT, WGSC (xgboost)
#   -- each with all 3 regression methods (lm, ranger, xgboost)
# Scenarios: S1U, S2U, S3U (unbounded density ratio)
# n_per_group: {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
#
# Output: results/appendix/oracle/oracle_comparison.csv
#
# Usage (from package root):
#   Rscript scripts/appendix/oracle_cit.R
# =============================================================================

source("R/load_all.R")
set.seed(1203)

# =============================================================================
# 1. Experiment Parameters
# =============================================================================

n_values   <- c(200, 500, 1000, 2000)
n_rep      <- 500
alpha      <- 0.05
seed_base  <- 1203

scenario_names <- c("S1U", "S2U", "S3U")

# =============================================================================
# 2. Oracle Data Generator
#    Generate n_total i.i.d. triples (X_i, Y_i, Z_i) from P_{XYZ}
#    where Z_i ~ Bernoulli(0.5), (X_i, Y_i) | Z_i ~ P_{XY}^{(Z_i)}
# =============================================================================

generate_oracle_data <- function(sc, n_total, is_null, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    Z <- rbinom(n_total, 1, 0.5)    # Z in {0, 1}
    idx1 <- which(Z == 0)
    idx2 <- which(Z == 1)

    gen_data <- sc$gen_data
    gen_y    <- sc$gen_y

    # Determine dimensionality from a small probe
    probe <- gen_data(2L, 1L)
    p <- ncol(probe)

    X <- matrix(0, nrow = n_total, ncol = p)
    Y <- numeric(n_total)

    if (length(idx1) > 0) {
        X[idx1, ] <- gen_data(length(idx1), 1L)
        Y[idx1]   <- gen_y(X[idx1, , drop = FALSE], is_null = TRUE)
    }

    if (length(idx2) > 0) {
        X[idx2, ] <- gen_data(length(idx2), 2L)
        Y[idx2]   <- gen_y(X[idx2, , drop = FALSE], is_null = is_null)
    }

    list(X = X, Y = Y, Z = Z, n1 = length(idx1), n2 = length(idx2))
}

# =============================================================================
# 3. Oracle CIT Test Runner
#    Apply CIT test directly on (X, Y, Z) without Algorithm 1
# =============================================================================

run_oracle_cit <- function(X, Y, Z, cit_method, reg_method, binary_reg_method,
                           alpha = 0.05, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # CIT mapping: X_cit = Y (response), Y_cit = Z (binary), Z_cit = X (covariate)
    p_value <- tryCatch({
        if (cit_method == "GCM") {
            gcm_test_binary(X = Y, Y = Z, Z = X,
                            reg_method = reg_method,
                            binary_reg_method = binary_reg_method, seed = seed)
        } else if (cit_method == "PCM") {
            pcm_test_binary(Y = Z, X = Y, Z = X,
                            reg_method = reg_method,
                            binary_reg_method = binary_reg_method, seed = seed)
        } else if (cit_method == "RCIT") {
            RCIT::RCIT(x = Y, y = Z, z = X,
                       approx = "lpd4", num_f = 100, num_f2 = 5, seed = seed)$p
        } else if (cit_method == "WGSC") {
            wgsc_binary(Y = Z, X = Y, Z = X,
                        reg_method = reg_method,
                        binary_reg_method = binary_reg_method, seed = seed)
        } else {
            stop("Unknown CIT method: ", cit_method)
        }
    }, error = function(e) NA)

    as.integer(!is.na(p_value) && p_value < alpha)
}

# =============================================================================
# 4. Test Configurations
#    Each CIT method x each regression method = full grid
# =============================================================================

regression_methods <- list(
    lm      = list(reg = lm_reg_method,      breg = lm_reg_method_binary),
    ranger  = list(reg = ranger_reg_method,   breg = ranger_reg_method_binary),
    xgboost = list(reg = xgboost_reg_method,  breg = xgboost_reg_method_binary)
)

cit_methods <- c("GCM", "PCM", "RCIT", "WGSC")

# Build complete config list
test_configs <- list()
for (cit_m in cit_methods) {
    if (cit_m == "RCIT") {
        # RCIT has no regression method
        test_configs[[length(test_configs) + 1]] <- list(
            method = "RCIT", reg_name = NA_character_,
            reg = NULL, binary_reg = NULL
        )
    } else {
        for (rn in names(regression_methods)) {
            rm_obj <- regression_methods[[rn]]
            test_configs[[length(test_configs) + 1]] <- list(
                method = cit_m, reg_name = rn,
                reg = rm_obj$reg, binary_reg = rm_obj$breg
            )
        }
    }
}

# =============================================================================
# 5. Main Simulation
# =============================================================================

results_list <- list()

for (sc_name in scenario_names) {
    sc <- SCENARIOS[[sc_name]]
    cat(sprintf("\n========== Scenario: %s ==========\n", sc_name))

    for (n in n_values) {
        n_total <- 2 * n   # total sample size = n (group1) + n (group2)
        eps <- 1 / sqrt(log(n))

        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"

            for (tc in test_configs) {
                label <- if (is.na(tc$reg_name)) tc$method
                         else paste0(tc$method, "_", tc$reg_name)
                cat(sprintf("[%s] %s | n_per_group=%d (total=%d) | %s\n",
                            sc_name, label, n, n_total, h_label))

                # --- C2ST (with Algorithm 1) ---
                cat("  C2ST (Alg1)...")
                res_c2st <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    set.seed(seed)
                    x1 <- sc$gen_data(n, 1L)
                    y1 <- sc$gen_y(x1, is_null = TRUE)
                    set.seed(seed + n_rep)
                    x2 <- sc$gen_data(n, 2L)
                    y2 <- sc$gen_y(x2, is_null = is_null)

                    tryCatch({
                        if (tc$method == "GCM") {
                            GCM_test(x1, x2, y1, y2, alpha = alpha,
                                     epsilon = eps, regr.method = tc$reg,
                                     binary.regr.method = tc$binary_reg,
                                     alg1 = TRUE, seed = seed)
                        } else if (tc$method == "PCM") {
                            PCM_test(x1, x2, y1, y2, alpha = alpha,
                                     epsilon = eps, regr.method = tc$reg,
                                     binary.regr.method = tc$binary_reg,
                                     alg1 = TRUE, seed = seed)
                        } else if (tc$method == "RCIT") {
                            RCIT_test(x1, x2, y1, y2, alpha = alpha,
                                      epsilon = eps, alg1 = TRUE, seed = seed)
                        } else if (tc$method == "WGSC") {
                            WGSC_test(x1, x2, y1, y2, alpha = alpha,
                                      epsilon = eps, regr.method = tc$reg,
                                      binary.regr.method = tc$binary_reg,
                                      alg1 = TRUE, seed = seed)
                        }
                    }, error = function(e) NA_integer_)
                })
                rr_c2st <- mean(res_c2st, na.rm = TRUE)
                cat(sprintf(" RR=%.3f\n", rr_c2st))

                results_list[[length(results_list) + 1]] <- data.table(
                    scenario = sc_name, test_name = tc$method,
                    regression_method = tc$reg_name,
                    n_per_group = n, n_total = n_total,
                    h_label = h_label, setting = "C2ST_Alg1",
                    rejection_rate = rr_c2st
                )

                # --- Oracle CIT (direct i.i.d. samples, no subsampling) ---
                cat("  Oracle CIT...")
                res_oracle <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim + 50000  # different seed from C2ST
                    oracle_data <- generate_oracle_data(
                        sc, n_total = n_total, is_null = is_null, seed = seed
                    )

                    # Skip if either group has < 5 samples (extremely unlikely)
                    if (oracle_data$n1 < 5 || oracle_data$n2 < 5) return(NA_integer_)

                    tryCatch(
                        run_oracle_cit(
                            oracle_data$X, oracle_data$Y, oracle_data$Z,
                            cit_method = tc$method,
                            reg_method = tc$reg,
                            binary_reg_method = tc$binary_reg,
                            alpha = alpha, seed = seed
                        ),
                        error = function(e) NA_integer_
                    )
                })
                rr_oracle <- mean(res_oracle, na.rm = TRUE)
                cat(sprintf(" RR=%.3f\n", rr_oracle))

                results_list[[length(results_list) + 1]] <- data.table(
                    scenario = sc_name, test_name = tc$method,
                    regression_method = tc$reg_name,
                    n_per_group = n, n_total = n_total,
                    h_label = h_label, setting = "Oracle_CIT",
                    rejection_rate = rr_oracle
                )

                cat(sprintf("  >> C2ST=%.3f, Oracle=%.3f, Diff=%.3f\n",
                            rr_c2st, rr_oracle, rr_c2st - rr_oracle))
            }
        }
    }
}

# =============================================================================
# 6. Save Results
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/oracle/oracle_comparison.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# Print summary
cat("\n========== Summary ==========\n")
summary_dt <- results_dt[, .(rejection_rate = mean(rejection_rate)),
                          by = .(scenario, test_name, regression_method,
                                 n_per_group, h_label, setting)]
summary_wide <- dcast(summary_dt,
                       scenario + test_name + regression_method + n_per_group + h_label ~ setting,
                       value.var = "rejection_rate")
if ("C2ST_Alg1" %in% names(summary_wide) && "Oracle_CIT" %in% names(summary_wide)) {
    summary_wide[, power_diff := C2ST_Alg1 - Oracle_CIT]
}
print(summary_wide)
