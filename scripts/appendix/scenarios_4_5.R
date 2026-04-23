# =============================================================================
# scripts/appendix/scenarios_4_5.R -- Scenarios 4-5 (U and B variants)
#
# Runs all 3 test groups for Scenarios 4 and 5:
#   CIT : GCM (ranger), PCM (ranger), RCIT, WGSC (xgb) -- with Algorithm 1
#   DRT : CV_LinearMMD, CV_BlockMMD, CV_QuadraticMMD, CV_CLF,
#         CP, DCP -- est.method="LL"
#   CK  : CMMD, CGED
#
# IMPORTANT: S4U/S4B generate_y takes an extra `group` parameter.
#   Under the alternative, y2 must be generated with group=2 to get the
#   location shift delta=0.5.
#
# Scenarios: S4U, S4B, S5U, S5B
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
# eps = 1/sqrt(log(n)) for CIT methods
#
# Output:
#   results/appendix/scenario_4_5/simulation_results_{tag}_{cit/drt/ck}.csv
#
# Usage (from package root):
#   Rscript scripts/appendix/scenarios_4_5.R
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

scenario_names <- c("S4U", "S4B", "S5U", "S5B")

.env_n_values <- Sys.getenv("CDTST_N_VALUES", "")
if (nzchar(.env_n_values)) {
    n_values <- as.integer(strsplit(.env_n_values, ",", fixed = TRUE)[[1]])
}
.env_scenarios <- Sys.getenv("CDTST_SCENARIOS", "")
if (nzchar(.env_scenarios)) {
    scenario_names <- strsplit(.env_scenarios, ",", fixed = TRUE)[[1]]
}
rm(.env_n_values, .env_scenarios)

# Scenarios that need group= parameter in gen_y
needs_group <- c("S4U", "S4B")

# =============================================================================
# 2. Test Configurations
# =============================================================================

# --- CIT ---
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

# --- DRT ---
drt_tests <- list(
    CV_LinearMMD_test    = CV_LinearMMD_test,
    CV_BlockMMD_test     = CV_BlockMMD_test,
    CV_QuadraticMMD_test = CV_QuadraticMMD_test,
    CV_CLF_test          = CV_CLF_test,
    CP_test              = CP_test,
    DCP_test             = DCP_test
)

# --- CK ---
ck_tests <- list(
    CMMD_test = CMMD_test,
    CGED_test = CGED_test
)

# =============================================================================
# 3. Helper: generate one pair of samples (handles group= for S4)
# =============================================================================

gen_pair <- function(sc, sc_name, n, is_null, seed) {
    set.seed(seed)
    x1 <- sc$gen_data(n, 1L)
    if (sc_name %in% needs_group) {
        y1 <- sc$gen_y(x1, is_null = TRUE, group = 1L)
    } else {
        y1 <- sc$gen_y(x1, is_null = TRUE)
    }
    set.seed(seed + n_rep)
    x2 <- sc$gen_data(n, 2L)
    if (sc_name %in% needs_group) {
        y2 <- sc$gen_y(x2, is_null = is_null, group = 2L)
    } else {
        y2 <- sc$gen_y(x2, is_null = is_null)
    }
    list(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
}

# =============================================================================
# 4. Main Simulation Loop
# =============================================================================

out_dir <- file.path(RESULTS_DIR, "appendix/scenario_4_5")
ensure_dir(out_dir)

for (sc_name in scenario_names) {
    sc  <- SCENARIOS[[sc_name]]
    tag <- sc$tag
    cat(sprintf("\n########## Scenario: %s ##########\n", tag))

    # ------------------------------------------------------------------
    # 4a. CIT methods
    # ------------------------------------------------------------------
    cat(sprintf("\n--- [%s] CIT methods ---\n", tag))
    cit_results <- list()

    for (n in n_values) {
        eps <- 1 / sqrt(log(n))
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (tc in cit_configs) {
                cat(sprintf("[%s][CIT] %s(%s) | n=%d | %s\n",
                            tag, tc$name, tc$reg_name, n, h_label))

                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    d <- gen_pair(sc, sc_name, n, is_null, seed)

                    tryCatch({
                        args <- list(d$x1, d$x2, d$y1, d$y2,
                                     alg1 = TRUE, epsilon = eps,
                                     alpha = alpha, seed = seed)
                        if (!is.null(tc$reg)) {
                            args$regr.method <- tc$reg
                            args$binary.regr.method <- tc$breg
                        }
                        do.call(tc$fn, args)
                    }, error = function(e) NA_integer_)
                })

                rr <- mean(res, na.rm = TRUE)
                cat(sprintf("  RR=%.3f\n", rr))

                cit_results[[length(cit_results) + 1]] <- data.table(
                    scenario = tag, test_type = "CIT", test_name = tc$name,
                    estimator = tc$reg_name, n = n, hypothesis = h_label,
                    alpha = alpha, rejection_rate = rr
                )
            }
        }
    }

    cit_dt   <- rbindlist(cit_results)
    cit_file <- file.path(out_dir, sprintf("simulation_results_%s_cit.csv", tag))
    fwrite(cit_dt, cit_file)
    cat(sprintf("Saved: %s\n", cit_file))

    # ------------------------------------------------------------------
    # 4b. DRT methods
    # ------------------------------------------------------------------
    cat(sprintf("\n--- [%s] DRT methods ---\n", tag))
    drt_results <- list()

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (test_name in names(drt_tests)) {
                cat(sprintf("[%s][DRT] %s | n=%d | %s\n",
                            tag, test_name, n, h_label))

                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    d <- gen_pair(sc, sc_name, n, is_null, seed)

                    tryCatch(
                        drt_tests[[test_name]](d$x1, d$x2, d$y1, d$y2,
                                               est.method = "LL",
                                               alpha = alpha, seed = seed),
                        error = function(e) NA_integer_
                    )
                })

                rr <- mean(res, na.rm = TRUE)
                cat(sprintf("  RR=%.3f\n", rr))

                drt_results[[length(drt_results) + 1]] <- data.table(
                    scenario = tag, test_type = "DRT", test_name = test_name,
                    estimator = "LL", n = n, hypothesis = h_label,
                    alpha = alpha, rejection_rate = rr
                )
            }
        }
    }

    drt_dt   <- rbindlist(drt_results)
    drt_file <- file.path(out_dir, sprintf("simulation_results_%s_drt.csv", tag))
    fwrite(drt_dt, drt_file)
    cat(sprintf("Saved: %s\n", drt_file))

    # ------------------------------------------------------------------
    # 4c. CK (competitor) methods
    # ------------------------------------------------------------------
    if (nzchar(Sys.getenv("CDTST_SKIP_CK", ""))) {
        cat(sprintf("\n--- [%s] CK methods SKIPPED (CDTST_SKIP_CK) ---\n", tag))
        next
    }
    cat(sprintf("\n--- [%s] CK methods ---\n", tag))
    ck_results <- list()

    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (test_name in names(ck_tests)) {
                cat(sprintf("[%s][CK] %s | n=%d | %s\n",
                            tag, test_name, n, h_label))

                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    d <- gen_pair(sc, sc_name, n, is_null, seed)

                    tryCatch(
                        ck_tests[[test_name]](d$x1, d$x2, d$y1, d$y2,
                                              alpha = alpha, seed = seed),
                        error = function(e) NA_integer_
                    )
                })

                rr <- mean(res, na.rm = TRUE)
                cat(sprintf("  RR=%.3f\n", rr))

                ck_results[[length(ck_results) + 1]] <- data.table(
                    scenario = tag, test_type = "Competitor", test_name = test_name,
                    estimator = NA_character_, n = n, hypothesis = h_label,
                    alpha = alpha, rejection_rate = rr
                )
            }
        }
    }

    ck_dt   <- rbindlist(ck_results)
    ck_file <- file.path(out_dir, sprintf("simulation_results_%s_ck.csv", tag))
    fwrite(ck_dt, ck_file)
    cat(sprintf("Saved: %s\n", ck_file))
}

cat("\n========== All Scenarios 4-5 Complete ==========\n")
