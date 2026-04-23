# =============================================================================
# scripts/appendix/kstar.R -- Refined k* bound comparison
#
# Compares Algorithm 1 with two k* formulations:
#   (1) "chernoff" : k via Chernoff bound (apply_alg1)
#   (2) "refined"  : k via exact Binomial tail (apply_alg1_refined)
#
# CIT Methods: GCM (ranger), PCM (ranger), RCIT, WGSC (xgb)
# Scenarios: S1U, S1B, S2U, S2B, S3U, S3B (all 6)
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
#
# Output: results/appendix/kstar/kstar_comparison.csv
# Columns: scenario, test_name, k_method, n, hypothesis, rejection_rate
#
# Usage (from package root):
#   Rscript scripts/appendix/kstar.R
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

scenario_names <- c("S1U", "S1B", "S2U", "S2B", "S3U", "S3B")
k_methods      <- c("chernoff", "refined")

# =============================================================================
# 2. CIT Test Configurations
# =============================================================================

cit_configs <- list(
    list(name = "GCM_test",  fn = GCM_test,  reg_name = "ranger",
         reg = ranger_reg_method, breg = ranger_reg_method_binary),
    list(name = "PCM_test",  fn = PCM_test,  reg_name = "ranger",
         reg = ranger_reg_method, breg = ranger_reg_method_binary),
    list(name = "RCIT_test", fn = RCIT_test, reg_name = NA_character_,
         reg = NULL, breg = NULL),
    list(name = "WGSC_test", fn = WGSC_test, reg_name = "xgb",
         reg = xgboost_reg_method, breg = xgboost_reg_method_binary)
)

# =============================================================================
# 3. apply_alg1 swap mechanism
#    Temporarily replace apply_alg1 in the global environment so that
#    CIT wrappers (GCM_test, etc.) call the selected variant.
# =============================================================================

apply_alg1_chernoff <- apply_alg1   # save original

swap_alg1 <- function(k_method) {
    if (k_method == "refined") {
        assign("apply_alg1", apply_alg1_refined, envir = .GlobalEnv)
    } else {
        assign("apply_alg1", apply_alg1_chernoff, envir = .GlobalEnv)
    }
}

restore_alg1 <- function() {
    assign("apply_alg1", apply_alg1_chernoff, envir = .GlobalEnv)
}

# =============================================================================
# 4. Reference table: display k* values
# =============================================================================

cat("\n===== k* vs k*_refined reference table =====\n")
cat(sprintf("%-6s %-10s %-12s %-12s\n", "n", "epsilon", "k_chernoff", "k_refined"))
cat(strrep("-", 50), "\n")
for (n in n_values) {
    n1 <- n; n2 <- n; n_total <- n1 + n2
    eps <- 1 / sqrt(log(n))
    # Chernoff
    n_min <- min(n1, n2)
    c_val <- (3 * log(eps)) / (2 * n_min)
    k_c <- 1 - c_val - sqrt((1 - c_val)^2 - 1)
    # Refined (binary search)
    find_k_single <- function(nj) {
        k_low <- 0.501; k_high <- 0.999; tol <- 1e-6
        while (k_high - k_low > tol) {
            k_mid <- (k_low + k_high) / 2
            tilde_n <- ceiling(k_mid * n_total)
            tail_p <- pbinom(nj, size = tilde_n, prob = nj / n_total,
                             lower.tail = FALSE)
            if (tail_p <= eps) k_low <- k_mid else k_high <- k_mid
        }
        k_low
    }
    k_r <- min(find_k_single(n1), find_k_single(n2))
    cat(sprintf("%-6d %-10.4f %-12.6f %-12.6f\n", n, eps, k_c, k_r))
}
cat("\n")

# =============================================================================
# 5. Main Simulation
# =============================================================================

results_list <- list()

for (sc_name in scenario_names) {
    sc <- SCENARIOS[[sc_name]]
    cat(sprintf("\n========== Scenario: %s ==========\n", sc_name))

    for (n in n_values) {
        eps <- 1 / sqrt(log(n))

        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"

            for (k_method in k_methods) {
                swap_alg1(k_method)

                for (tc in cit_configs) {
                    cat(sprintf("[%s] %s | k=%s | n=%d | %s\n",
                                sc_name, tc$name, k_method, n, h_label))

                    rej_vec <- cdtst_sapply(n_rep, function(sim) {
                        seed <- seed_base + sim
                        set.seed(seed)
                        x1 <- sc$gen_data(n, 1L)
                        y1 <- sc$gen_y(x1, is_null = TRUE)
                        set.seed(seed + n_rep)
                        x2 <- sc$gen_data(n, 2L)
                        y2 <- sc$gen_y(x2, is_null = is_null)

                        tryCatch({
                            args <- list(x1, x2, y1, y2,
                                         alg1 = TRUE, epsilon = eps,
                                         alpha = alpha, seed = seed)
                            if (!is.null(tc$reg)) {
                                args$regr.method <- tc$reg
                                args$binary.regr.method <- tc$breg
                            }
                            do.call(tc$fn, args)
                        }, error = function(e) NA_integer_)
                    })

                    rr <- mean(rej_vec, na.rm = TRUE)
                    cat(sprintf("  RR=%.3f\n", rr))

                    results_list[[length(results_list) + 1]] <- data.table(
                        scenario       = sc_name,
                        test_name      = tc$name,
                        k_method       = k_method,
                        n              = n,
                        hypothesis     = h_label,
                        rejection_rate = rr
                    )
                }
            }
        }
    }
}

restore_alg1()

# =============================================================================
# 6. Save and Summarize
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/kstar/kstar_comparison.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary: Chernoff vs Refined ---
cat("\n========== Chernoff vs Refined: Rejection Rate Difference ==========\n")
dt_c <- results_dt[k_method == "chernoff",
                    .(rr_chernoff = rejection_rate),
                    by = .(scenario, test_name, n, hypothesis)]
dt_r <- results_dt[k_method == "refined",
                    .(rr_refined = rejection_rate),
                    by = .(scenario, test_name, n, hypothesis)]
dt_cmp <- merge(dt_c, dt_r, by = c("scenario", "test_name", "n", "hypothesis"))
dt_cmp[, delta := rr_refined - rr_chernoff]
setorder(dt_cmp, scenario, test_name, n, hypothesis)
print(dt_cmp)
