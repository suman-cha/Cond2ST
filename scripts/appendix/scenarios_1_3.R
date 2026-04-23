# =============================================================================
# scripts/appendix/scenarios_1_3.R
#
# Appendix-track scenarios 1-3. Runs the 4 non-CV DRT tests that are reported
# in the appendix tables but not on the main-paper figures:
#   LinearMMD, BlockMMD, QuadraticMMD, CLF (learned+ranger)
#
# Main-track (CV_LinearMMD, CV_BlockMMD, CV_QuadraticMMD, CV_CLF_test,
# CP, DCP, RCIT, GCM, PCM, WGSC, CMMD, CGED) lives in
# scripts/main/scenarios_1_3.R.
#
# Scenarios: S1U, S1B, S2U, S2B, S3U, S3B
# n in {200, 500, 1000, 2000}, n_rep=500, alpha=0.05
# est.method = "LL" (single estimator; KLR is only swept in real-data runner)
#
# Output: results/appendix/scenarios_1_3/simulation_results_{tag}_drt.csv
#
# Env overrides: CDTST_N_REP, CDTST_N_VALUES, CDTST_SCENARIOS, CDTST_CORES
#
# Usage (from package root):
#   CDTST_CORES=8 Rscript scripts/appendix/scenarios_1_3.R
# =============================================================================

source("R/load_all.R")
set.seed(1203)

# ---- Parameters -------------------------------------------------------------
n_values  <- c(200, 500, 1000, 2000)
n_rep     <- as.integer(Sys.getenv("CDTST_N_REP", "500"))
alpha     <- 0.05
seed_base <- 1203
scenario_names <- c("S1U", "S1B", "S2U", "S2B", "S3U", "S3B")

.env_n_values <- Sys.getenv("CDTST_N_VALUES", "")
if (nzchar(.env_n_values)) {
    n_values <- as.integer(strsplit(.env_n_values, ",", fixed = TRUE)[[1]])
}
.env_scenarios <- Sys.getenv("CDTST_SCENARIOS", "")
if (nzchar(.env_scenarios)) {
    scenario_names <- strsplit(.env_scenarios, ",", fixed = TRUE)[[1]]
}
rm(.env_n_values, .env_scenarios)

# ---- Appendix DRT test list -------------------------------------------------
drt_tests <- list(
    LinearMMD_test    = LinearMMD_test,
    BlockMMD_test     = BlockMMD_test,
    QuadraticMMD_test = QuadraticMMD_test,
    CLF_test          = CLF_test          # learned + ranger
)

# ---- Data helper (matches scripts/main/scenarios_1_3.R) -------------------
gen_pair <- function(sc, n, is_null, seed) {
    set.seed(seed)
    x1 <- sc$gen_data(n, 1L); y1 <- sc$gen_y(x1, is_null = TRUE)
    set.seed(seed + n_rep)
    x2 <- sc$gen_data(n, 2L); y2 <- sc$gen_y(x2, is_null = is_null)
    list(x1 = x1, x2 = x2, y1 = y1, y2 = y2)
}

# ---- Main loop --------------------------------------------------------------
out_dir <- file.path(RESULTS_DIR, "appendix", "scenarios_1_3")
ensure_dir(out_dir)

for (sc_name in scenario_names) {
    sc  <- SCENARIOS[[sc_name]]
    tag <- sc$tag
    cat(sprintf("\n########## Scenario: %s (appendix DRT) ##########\n", tag))

    results <- list()
    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            for (test_name in names(drt_tests)) {
                cat(sprintf("[%s][DRT] %s | n=%d | %s\n",
                            tag, test_name, n, h_label))

                res <- cdtst_sapply(n_rep, function(sim) {
                    seed <- seed_base + sim
                    d <- gen_pair(sc, n, is_null, seed)
                    tryCatch(
                        drt_tests[[test_name]](d$x1, d$x2, d$y1, d$y2,
                                               est.method = "LL",
                                               alpha = alpha, seed = seed),
                        error = function(e) NA_integer_
                    )
                })

                rr <- mean(res, na.rm = TRUE)
                cat(sprintf("  RR=%.3f\n", rr))

                results[[length(results) + 1]] <- data.table(
                    scenario = tag, test_type = "DRT", test_name = test_name,
                    estimator = "LL", n = n, hypothesis = h_label,
                    alpha = alpha, rejection_rate = rr
                )
            }
        }
    }

    dt   <- rbindlist(results)
    file <- file.path(out_dir, sprintf("simulation_results_%s_drt.csv", tag))
    fwrite(dt, file)
    cat(sprintf("Saved: %s (%d rows)\n", file, nrow(dt)))
}

cat("\n========== Appendix Scenarios 1-3 Complete ==========\n")
