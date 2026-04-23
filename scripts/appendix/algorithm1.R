# =============================================================================
# scripts/appendix/algorithm1.R
#
# Algorithm 1 sensitivity: with/without, using different regression methods.
#
# Purpose: Show that Algorithm 1 subsampling does not significantly degrade
#          CIT performance, and compare across regression backends.
#
# CIT Methods: GCM, PCM, WGSC -- each with regression methods: lm, ranger, xgboost
#              RCIT -- no regression method
# alg1 in {TRUE, FALSE}
# Scenarios: S1U, S2U, S3U
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
#
# Output: results/appendix/algorithm1/alg1_sensitivity.csv
#
# Usage (from package root):
#   Rscript scripts/appendix/algorithm1.R
# =============================================================================

source("R/load_all.R")
set.seed(1203)

# =============================================================================
# 1. Experiment Parameters
# =============================================================================

n_values   <- c(200, 500, 1000, 2000)
n_rep      <- as.integer(Sys.getenv("CDTST_N_REP", "500"))
alpha      <- 0.05
seed_base  <- 1203

scenario_names <- c("S1U", "S2U", "S3U")
alg1_list      <- c(TRUE, FALSE)

.env_n_values <- Sys.getenv("CDTST_N_VALUES", "")
if (nzchar(.env_n_values)) {
    n_values <- as.integer(strsplit(.env_n_values, ",", fixed = TRUE)[[1]])
}
.env_scenarios <- Sys.getenv("CDTST_SCENARIOS", "")
if (nzchar(.env_scenarios)) {
    scenario_names <- strsplit(.env_scenarios, ",", fixed = TRUE)[[1]]
}
rm(.env_n_values, .env_scenarios)

# =============================================================================
# 2. Test & Regression Configurations
# =============================================================================

# CIT test functions (all accept regr.method, binary.regr.method, alg1)
cit_test_functions <- list(
    GCM_test  = GCM_test,
    PCM_test  = PCM_test,
    WGSC_test = WGSC_test
)

regression_methods <- list(
    lm      = list(reg = lm_reg_method,      breg = lm_reg_method_binary),
    ranger  = list(reg = ranger_reg_method,   breg = ranger_reg_method_binary),
    xgboost = list(reg = xgboost_reg_method,  breg = xgboost_reg_method_binary)
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

            # --- CIT methods with regression ---
            for (test_name in names(cit_test_functions)) {
                for (regressor in names(regression_methods)) {
                    for (alg1 in alg1_list) {
                        cat(sprintf("[%s] %s | regr=%s | alg1=%s | n=%d | %s\n",
                                    sc_name, test_name, regressor,
                                    as.character(alg1), n, h_label))

                        rm_obj <- regression_methods[[regressor]]

                        res <- cdtst_sapply(n_rep, function(sim) {
                            seed <- seed_base + sim
                            set.seed(seed)
                            x1 <- sc$gen_data(n, 1L)
                            y1 <- sc$gen_y(x1, is_null = TRUE)
                            set.seed(seed + n_rep)
                            x2 <- sc$gen_data(n, 2L)
                            y2 <- sc$gen_y(x2, is_null = is_null)

                            tryCatch({
                                test_args <- list(
                                    x1, x2, y1, y2,
                                    seed = seed, alg1 = alg1, alpha = alpha,
                                    regr.method = rm_obj$reg,
                                    binary.regr.method = rm_obj$breg
                                )
                                do.call(cit_test_functions[[test_name]], test_args)
                            }, error = function(e) NA_integer_)
                        })

                        rr <- mean(res, na.rm = TRUE)
                        cat(sprintf("  RR=%.3f\n", rr))

                        results_list[[length(results_list) + 1]] <- data.table(
                            scenario          = sc_name,
                            test_name         = test_name,
                            regression_method = regressor,
                            n                 = n,
                            h_label           = h_label,
                            alg1              = alg1,
                            rejection_rate    = rr
                        )
                    }
                }
            }

            # --- RCIT (no regression method) ---
            for (alg1 in alg1_list) {
                cat(sprintf("[%s] RCIT_test | alg1=%s | n=%d | %s\n",
                            sc_name, as.character(alg1), n, h_label))

                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    set.seed(seed)
                    x1 <- sc$gen_data(n, 1L)
                    y1 <- sc$gen_y(x1, is_null = TRUE)
                    set.seed(seed + n_rep)
                    x2 <- sc$gen_data(n, 2L)
                    y2 <- sc$gen_y(x2, is_null = is_null)

                    tryCatch(
                        RCIT_test(x1, x2, y1, y2,
                                  seed = seed, alg1 = alg1, alpha = alpha),
                        error = function(e) NA_integer_
                    )
                })

                rr <- mean(res, na.rm = TRUE)
                cat(sprintf("  RR=%.3f\n", rr))

                results_list[[length(results_list) + 1]] <- data.table(
                    scenario          = sc_name,
                    test_name         = "RCIT_test",
                    regression_method = NA_character_,
                    n                 = n,
                    h_label           = h_label,
                    alg1              = alg1,
                    rejection_rate    = rr
                )
            }
        }
    }
}

# =============================================================================
# 4. Save and Summarize
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/algorithm1/alg1_sensitivity.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary: Wide format with alg1 TRUE/FALSE columns ---
cat("\n========== Rejection Rates: With vs Without Algorithm 1 ==========\n")
dt_wide <- dcast(results_dt,
                  scenario + test_name + regression_method + n + h_label ~ paste0("Alg1_", alg1),
                  value.var = "rejection_rate")
setorder(dt_wide, scenario, test_name, regression_method, n, h_label)

for (sc_name in scenario_names) {
    cat(sprintf("\n--- Scenario: %s ---\n", sc_name))
    for (h in c("Null", "Alternative")) {
        cat(sprintf("  [%s]\n", h))
        sub <- dt_wide[scenario == sc_name & h_label == h]
        print(sub[, -c("scenario", "h_label")])
    }
}
