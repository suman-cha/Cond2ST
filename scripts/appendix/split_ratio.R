# =============================================================================
# scripts/appendix/split_ratio.R
#
# Split ratio sensitivity analysis.
#
# Purpose: Demonstrate the trade-off between density ratio estimation quality
#          and effective test sample size when varying the split proportion.
#
# Tests (8, every DRT with a meaningful split-proportion knob):
#   LinearMMD, CV_LinearMMD, BlockMMD, QuadraticMMD, CV_QuadraticMMD,
#   CLF, CV_CLF, CP
# (CV_BlockMMD and DCP are excluded because they don't expose a split-ratio parameter.)
# Split proportions: prop in {0.2, 0.5, 0.8}
#   - For MMD / CP: `prop` = fraction allocated to density ratio estimation
#   - For CV_MMD:   `prop` = same, combined with K=2 folds
#   - For CLF:      `est.prop` = fraction for ratio est. (`split.prop = 0.5` fixed)
#   - For CV_CLF:   `est.prop` = fraction within each fold (K=2)
# Scenarios: S1U, S1B, S2U, S2B, S3U, S3B
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
# est.method = "LL"
#
# Output: results/appendix/split_ratio/split_ratio_sensitivity.csv
#
# Usage (from package root):
#   Rscript scripts/appendix/split_ratio.R
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

scenario_names <- c("S1U", "S1B", "S2U", "S2B", "S3U", "S3B")
prop_values    <- c(0.2, 0.5, 0.8)

# =============================================================================
# 2. Test Configurations
#
# Each test function is called with its native `prop` parameter.
# CLF_test uses `est.prop` (with `split.prop = 0.5` fixed).
# =============================================================================

test_configs <- list(
    list(
        name = "LinearMMD",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            LinearMMD_test(x1, x2, y1, y2,
                           prop = prop, alpha = alpha,
                           est.method = est, seed = seed)
        }
    ),
    list(
        name = "CV_LinearMMD",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            CV_LinearMMD_test(x1, x2, y1, y2,
                              prop = prop, K = 2, alpha = alpha,
                              est.method = est, seed = seed)
        }
    ),
    list(
        name = "BlockMMD",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            BlockMMD_test(x1, x2, y1, y2,
                          prop = prop, alpha = alpha,
                          est.method = est, seed = seed)
        }
    ),
    list(
        name = "QuadraticMMD",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            QuadraticMMD_test(x1, x2, y1, y2,
                              prop = prop, alpha = alpha,
                              est.method = est, seed = seed)
        }
    ),
    list(
        name = "CV_QuadraticMMD",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            CV_QuadraticMMD_test(x1, x2, y1, y2,
                                 prop = prop, K = 2, alpha = alpha,
                                 est.method = est, seed = seed)
        }
    ),
    list(
        name = "CLF",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            CLF_test(x1, x2, y1, y2,
                     clf_split = 0.5, est.prop = prop,
                     alpha = alpha, est.method = est, seed = seed)
        }
    ),
    list(
        name = "CV_CLF",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            CV_CLF_test(x1, x2, y1, y2,
                        clf_ratio_prop = prop, K = 2, alpha = alpha,
                        est.method = est, seed = seed)
        }
    ),
    list(
        name = "CP",
        call_fn = function(x1, x2, y1, y2, prop, est, alpha, seed) {
            CP_test(x1, x2, y1, y2,
                    prop = prop, alpha = alpha,
                    est.method = est, seed = seed)
        }
    )
)

# =============================================================================
# 3. Main Simulation
# =============================================================================

results_list <- list()

for (sc_name in scenario_names) {
    sc <- SCENARIOS[[sc_name]]
    cat(sprintf("\n========== Scenario: %s ==========\n", sc$tag))

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"

            for (prop in prop_values) {
                for (tc in test_configs) {
                    # CV_QuadraticMMD with K=2 requires test-fold size m = (1-prop)*n
                    # to satisfy K <= floor(n/m). At prop=0.2, m=0.8n, floor(n/m)=1 < K=2.
                    # Structurally infeasible; skip rather than emit NA rows.
                    if (tc$name == "CV_QuadraticMMD" && prop == 0.2) {
                        cat(sprintf("[%s] %s | n=%d | prop=%.1f | %s -- SKIP (K=2 infeasible)\n",
                                    sc$tag, tc$name, n, prop, h_label))
                        next
                    }
                    cat(sprintf("[%s] %s | n=%d | prop=%.1f | %s\n",
                                sc$tag, tc$name, n, prop, h_label))

                    rej_vec <- cdtst_sapply(n_rep, function(sim) {
                        seed <- seed_base + sim
                        set.seed(seed)
                        x1 <- sc$gen_data(n, 1L)
                        y1 <- sc$gen_y(x1, is_null = TRUE)
                        set.seed(seed + n_rep)
                        x2 <- sc$gen_data(n, 2L)
                        y2 <- sc$gen_y(x2, is_null = is_null)

                        tryCatch(
                            tc$call_fn(x1, x2, y1, y2,
                                       prop = prop, est = est_method,
                                       alpha = alpha, seed = seed),
                            error = function(e) NA_integer_
                        )
                    })

                    rr      <- mean(rej_vec, na.rm = TRUE)
                    na_rate <- mean(is.na(rej_vec))
                    cat(sprintf("  RR=%.3f (NA=%.3f)\n", rr, na_rate))

                    results_list[[length(results_list) + 1]] <- data.table(
                        scenario       = sc$tag,
                        test_name      = tc$name,
                        prop           = prop,
                        n              = n,
                        hypothesis     = h_label,
                        rejection_rate = rr,
                        na_rate        = na_rate
                    )
                }
            }
        }
    }
}

# =============================================================================
# 4. Save and Summarize
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/split_ratio/split_ratio_sensitivity.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary: Type I error by prop ---
cat("\n========== Type I Error by Split Proportion ==========\n")
null_dt <- results_dt[hypothesis == "Null"]
null_wide <- dcast(null_dt, scenario + test_name + n ~ prop,
                   value.var = "rejection_rate")
setorder(null_wide, scenario, test_name, n)
print(null_wide)

# --- Summary: Power by prop ---
cat("\n========== Power by Split Proportion ==========\n")
alt_dt <- results_dt[hypothesis == "Alternative"]
alt_wide <- dcast(alt_dt, scenario + test_name + n ~ prop,
                  value.var = "rejection_rate")
setorder(alt_wide, scenario, test_name, n)
print(alt_wide)

# --- Trade-off table ---
cat("\n========== Trade-off Summary ==========\n")
tradeoff <- merge(
    null_dt[, .(typeI = rejection_rate), by = .(scenario, test_name, n, prop)],
    alt_dt[, .(power = rejection_rate), by = .(scenario, test_name, n, prop)],
    by = c("scenario", "test_name", "n", "prop"),
    all = TRUE
)
tradeoff[, typeI_inflation := abs(typeI - alpha)]
setorder(tradeoff, scenario, test_name, n, prop)
print(tradeoff)
