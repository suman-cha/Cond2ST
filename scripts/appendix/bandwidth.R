# =============================================================================
# scripts/appendix/bandwidth.R -- Bandwidth sensitivity for CV_LinearMMD
#
# Tests sensitivity of CV_LinearMMD (and LinearMMD for comparison) to
# bandwidth choice on bounded scenarios S1B, S2B, S3B.
#
# Bandwidths: {0.1, 0.5, 1.0, "median", 5, 10}
#   "median" -> compute median heuristic per sample pair, then pass numeric value
#
# Scenarios: S1B, S2B, S3B
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05, est.method="LL"
#
# Output: results/appendix/bandwidth/bandwidth_sensitivity.csv
# Columns: scenario, test_name, bandwidth, n, hypothesis, rejection_rate
#
# Usage (from package root):
#   Rscript scripts/appendix/bandwidth.R
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
est_method <- "LL"

scenario_names <- c("S1B", "S2B", "S3B")
bandwidths     <- c(0.1, 0.5, 1.0, "median", 5, 10)

# =============================================================================
# 2. Test Registry
# =============================================================================

test_functions <- list(
    LinearMMD_test    = LinearMMD_test,
    CV_LinearMMD_test = CV_LinearMMD_test
)

# =============================================================================
# 3. Main Simulation
# =============================================================================

results_list <- list()

for (sc_name in scenario_names) {
    sc <- SCENARIOS[[sc_name]]
    cat(sprintf("\n========== Scenario: %s ==========\n", sc_name))

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"

            for (test_name in names(test_functions)) {
                cat(sprintf("[%s] %s | n=%d | %s\n",
                            sc_name, test_name, n, h_label))

                result_matrix <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    set.seed(seed)
                    x1 <- sc$gen_data(n, 1L)
                    y1 <- sc$gen_y(x1, is_null = TRUE)
                    set.seed(seed + n_rep)
                    x2 <- sc$gen_data(n, 2L)
                    y2 <- sc$gen_y(x2, is_null = is_null)

                    tryCatch(
                        test_functions[[test_name]](
                            x1, x2, y1, y2,
                            bandwidths   = bandwidths,
                            est.method   = est_method,
                            alpha        = alpha,
                            seed         = seed
                        ),
                        error = function(e) rep(NA_integer_, length(bandwidths))
                    )
                })

                # cdtst_sapply returns a 6xN matrix in serial mode but a flat
                # vector of length 6*N in parallel mode (foreach .combine=c).
                # Reshape to unified rows=bandwidths, cols=reps layout.
                if (!is.matrix(result_matrix)) {
                    result_matrix <- matrix(result_matrix,
                                            nrow = length(bandwidths),
                                            ncol = n_rep)
                }

                # result_matrix: rows = bandwidths, cols = replications
                for (bw_idx in seq_along(bandwidths)) {
                    bw_label <- as.character(bandwidths[bw_idx])
                    rr <- mean(result_matrix[bw_idx, ], na.rm = TRUE)

                    results_list[[length(results_list) + 1]] <- data.table(
                        scenario       = sc_name,
                        test_name      = test_name,
                        bandwidth      = bw_label,
                        n              = n,
                        hypothesis     = h_label,
                        rejection_rate = rr
                    )
                }

                cat(sprintf("  Bandwidths done. RR range: [%.3f, %.3f]\n",
                            min(rowMeans(result_matrix, na.rm = TRUE)),
                            max(rowMeans(result_matrix, na.rm = TRUE))))
            }
        }
    }
}

# =============================================================================
# 4. Save and Summarize
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/bandwidth/bandwidth_sensitivity.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary ---
cat("\n========== Bandwidth Sensitivity Summary ==========\n")
for (sc_name in scenario_names) {
    cat(sprintf("\n--- Scenario: %s ---\n", sc_name))
    for (h in c("Null", "Alternative")) {
        cat(sprintf("  [%s]\n", h))
        sub <- results_dt[scenario == sc_name & hypothesis == h]
        dt_wide <- dcast(sub, test_name + n ~ bandwidth,
                         value.var = "rejection_rate")
        setorder(dt_wide, test_name, n)
        print(dt_wide)
    }
}
