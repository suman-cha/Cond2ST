# =============================================================================
# scripts/appendix/epsilon.R
#
# Epsilon sensitivity in Algorithm 1.
#
# Purpose: Evaluate how different choices of epsilon in Algorithm 1 affect
#          CIT test performance. Larger epsilon discards fewer samples
#          (higher power) but risks violating the approximate i.i.d. condition.
#
# CIT Methods: GCM (ranger), PCM (ranger), RCIT, WGSC (xgb)
# Epsilon values: 1/n, 1/sqrt(log(n)), 1/log(n), 1/sqrt(n)
# Scenarios: S1U, S2U, S3U
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
#
# Output: results/appendix/epsilon/epsilon_sensitivity.csv
#
# Usage (from package root):
#   Rscript scripts/appendix/epsilon.R
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
epsilon_types  <- c("1/n", "1/sqrt(log(n))", "1/log(n)", "1/sqrt(n)")

# =============================================================================
# 2. Test Configurations
# =============================================================================

test_configs <- list(
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
# 3. Helper: compute epsilon value from type string and n
# =============================================================================

compute_epsilon <- function(epsilon_type, n) {
    switch(epsilon_type,
        "1/n"            = 1 / n,
        "1/sqrt(log(n))" = 1 / sqrt(log(n)),
        "1/log(n)"       = 1 / log(n),
        "1/sqrt(n)"      = 1 / sqrt(n),
        stop("Unknown epsilon_type: ", epsilon_type)
    )
}

# =============================================================================
# 4. Main Simulation
# =============================================================================

results_list <- list()

for (sc_name in scenario_names) {
    sc <- SCENARIOS[[sc_name]]
    cat(sprintf("\n========== Scenario: %s ==========\n", sc_name))

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"

            for (epsilon_type in epsilon_types) {
                eps_val <- compute_epsilon(epsilon_type, n)

                for (tc in test_configs) {
                    cat(sprintf("[%s] %s | n=%d | eps=%s (%.4f) | %s\n",
                                sc_name, tc$name, n, epsilon_type, eps_val, h_label))

                    res <- cdtst_sapply(n_rep, function(sim) {
                        seed <- seed_base + sim
                        set.seed(seed)
                        x1 <- sc$gen_data(n, 1L)
                        y1 <- sc$gen_y(x1, is_null = TRUE)
                        set.seed(seed + n_rep)
                        x2 <- sc$gen_data(n, 2L)
                        y2 <- sc$gen_y(x2, is_null = is_null)

                        tryCatch({
                            args <- list(x1, x2, y1, y2,
                                         seed = seed, alg1 = TRUE,
                                         epsilon = eps_val, alpha = alpha)
                            if (!is.null(tc$reg)) {
                                args$regr.method <- tc$reg
                                args$binary.regr.method <- tc$breg
                            }
                            do.call(tc$fn, args)
                        }, error = function(e) NA_integer_)
                    })

                    rr <- mean(res, na.rm = TRUE)
                    cat(sprintf("  RR=%.3f\n", rr))

                    results_list[[length(results_list) + 1]] <- data.table(
                        scenario     = sc_name,
                        test_name    = tc$name,
                        n            = n,
                        epsilon_type = epsilon_type,
                        epsilon_val  = eps_val,
                        h_label      = h_label,
                        rejection_rate = rr
                    )
                }
            }
        }
    }
}

# =============================================================================
# 5. Save and Summarize
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/epsilon/epsilon_sensitivity.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary by epsilon type ---
cat("\n========== Rejection Rates by Epsilon ==========\n")
for (sc_name in scenario_names) {
    dt_wide <- dcast(results_dt[scenario == sc_name],
                      test_name + n + h_label ~ epsilon_type,
                      value.var = "rejection_rate")
    setorder(dt_wide, test_name, n, h_label)

    for (h in c("Null", "Alternative")) {
        cat(sprintf("\n--- %s / %s ---\n", sc_name, h))
        print(dt_wide[h_label == h, -"h_label"])
    }
}
