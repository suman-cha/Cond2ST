# =============================================================================
# scripts/appendix/computation_cost.R -- Computational cost benchmarks
#
# Measures elapsed time (system.time) for all methods on Scenario S1U.
# Each method/variant is timed over n_timing replications per sample size; the
# mean and sd of elapsed seconds are recorded.
#
# Variants:
#   DRT : est.method in {"LL", "KLR"}
#   CIT : regr.method in {lm, ranger, xgboost} (RCIT has no regression knob)
#   CK  : no variant
#
# NOTE: no per-call timeout cap -- each call is measured in full.
#
# Scenario: S1U only
# n in {200, 500, 1000, 2000}, n_timing reps per timing measurement
#
# Output: results/appendix/computation_cost/computation_cost.csv
# Columns: test_name, variant, n, mean_time_sec, sd_time_sec
#
# Usage (from package root):
#   Rscript scripts/appendix/computation_cost.R
# =============================================================================

source("R/load_all.R")
set.seed(1203)

# =============================================================================
# 1. Experiment Parameters
# =============================================================================

n_values   <- c(200, 500, 1000, 2000)
n_timing   <- 10        # replications per timing measurement
alpha      <- 0.05
seed_base  <- 1203

# Use S1U scenario (alternative hypothesis for meaningful computation)
sc <- SCENARIOS[["S1U"]]

# =============================================================================
# 2. Test Registry
# =============================================================================

# DRT methods: compare LL vs KLR density-ratio estimators
drt_tests <- list(
    LinearMMD_test       = LinearMMD_test,
    CV_LinearMMD_test    = CV_LinearMMD_test,
    BlockMMD_test        = BlockMMD_test,
    CV_BlockMMD_test     = CV_BlockMMD_test,
    QuadraticMMD_test    = QuadraticMMD_test,
    CV_QuadraticMMD_test = CV_QuadraticMMD_test,
    CP_test              = CP_test,
    DCP_test             = DCP_test,
    CLF_test             = CLF_test,
    CV_CLF_test          = CV_CLF_test
)
drt_est_methods <- c("LL", "KLR")

# CIT methods: compare lm / ranger / xgboost regressors (RCIT = regression-free)
cit_fns <- list(
    GCM_test  = GCM_test,
    PCM_test  = PCM_test,
    WGSC_test = WGSC_test
)
cit_regressors <- list(
    lm      = list(reg = lm_reg_method,      breg = lm_reg_method_binary),
    ranger  = list(reg = ranger_reg_method,  breg = ranger_reg_method_binary),
    xgboost = list(reg = xgboost_reg_method, breg = xgboost_reg_method_binary)
)

# CK methods: no variant knob
ck_tests <- list(
    CMMD_test = CMMD_test,
    CGED_test = CGED_test
)

# =============================================================================
# 3. Helper: time a single test call (no timeout cap)
# =============================================================================

time_one_call <- function(test_fn, x1, x2, y1, y2, extra_args, seed) {
    t0 <- proc.time()
    tryCatch(
        do.call(test_fn, c(list(x1, x2, y1, y2), extra_args, list(seed = seed))),
        error = function(e) NA_integer_
    )
    t1 <- proc.time()
    (t1 - t0)[["elapsed"]]
}

gen_sample <- function(sim, n) {
    seed <- seed_base + sim
    set.seed(seed)
    x1 <- sc$gen_data(n, 1L)
    y1 <- sc$gen_y(x1, is_null = FALSE)
    set.seed(seed + n_timing)
    x2 <- sc$gen_data(n, 2L)
    y2 <- sc$gen_y(x2, is_null = FALSE)
    list(x1 = x1, x2 = x2, y1 = y1, y2 = y2, seed = seed)
}

# =============================================================================
# 4. Main Timing Loop
# =============================================================================

results_list <- list()

for (n in n_values) {
    eps <- 1 / sqrt(log(n))
    cat(sprintf("\n========== n = %d ==========\n", n))

    # --- DRT methods: LL + KLR ---
    for (test_name in names(drt_tests)) {
        for (est in drt_est_methods) {
            cat(sprintf("[DRT] %s | est=%s | n=%d ...", test_name, est, n))

            times <- sapply(seq_len(n_timing), function(sim) {
                s <- gen_sample(sim, n)
                time_one_call(drt_tests[[test_name]],
                              s$x1, s$x2, s$y1, s$y2,
                              list(est.method = est, alpha = alpha), s$seed)
            })

            results_list[[length(results_list) + 1]] <- data.table(
                test_name     = test_name,
                variant       = est,
                n             = n,
                mean_time_sec = mean(times),
                sd_time_sec   = sd(times)
            )
            cat(sprintf(" mean=%.3fs sd=%.3fs\n", mean(times), sd(times)))
        }
    }

    # --- CIT methods: lm / ranger / xgboost (plus RCIT, no regressor) ---
    for (test_name in names(cit_fns)) {
        for (reg_name in names(cit_regressors)) {
            cat(sprintf("[CIT] %s | reg=%s | n=%d ...", test_name, reg_name, n))

            rg <- cit_regressors[[reg_name]]
            times <- sapply(seq_len(n_timing), function(sim) {
                s <- gen_sample(sim, n)
                extra <- list(alg1 = TRUE, epsilon = eps, alpha = alpha,
                              regr.method = rg$reg,
                              binary.regr.method = rg$breg)
                time_one_call(cit_fns[[test_name]],
                              s$x1, s$x2, s$y1, s$y2, extra, s$seed)
            })

            results_list[[length(results_list) + 1]] <- data.table(
                test_name     = test_name,
                variant       = reg_name,
                n             = n,
                mean_time_sec = mean(times),
                sd_time_sec   = sd(times)
            )
            cat(sprintf(" mean=%.3fs sd=%.3fs\n", mean(times), sd(times)))
        }
    }

    # RCIT: regression-free
    cat(sprintf("[CIT] RCIT_test | variant=none | n=%d ...", n))
    times <- sapply(seq_len(n_timing), function(sim) {
        s <- gen_sample(sim, n)
        time_one_call(RCIT_test, s$x1, s$x2, s$y1, s$y2,
                      list(alg1 = TRUE, epsilon = eps, alpha = alpha), s$seed)
    })
    results_list[[length(results_list) + 1]] <- data.table(
        test_name     = "RCIT_test",
        variant       = "none",
        n             = n,
        mean_time_sec = mean(times),
        sd_time_sec   = sd(times)
    )
    cat(sprintf(" mean=%.3fs sd=%.3fs\n", mean(times), sd(times)))

    # --- CK methods ---
    for (test_name in names(ck_tests)) {
        cat(sprintf("[CK]  %s | n=%d ...", test_name, n))

        times <- sapply(seq_len(n_timing), function(sim) {
            s <- gen_sample(sim, n)
            time_one_call(ck_tests[[test_name]],
                          s$x1, s$x2, s$y1, s$y2,
                          list(alpha = alpha), s$seed)
        })

        results_list[[length(results_list) + 1]] <- data.table(
            test_name     = test_name,
            variant       = "none",
            n             = n,
            mean_time_sec = mean(times),
            sd_time_sec   = sd(times)
        )
        cat(sprintf(" mean=%.3fs sd=%.3fs\n", mean(times), sd(times)))
    }
}

# =============================================================================
# 5. Save and Summarize
# =============================================================================

results_dt <- rbindlist(results_list)
out_file <- file.path(RESULTS_DIR, "appendix/computation_cost/computation_cost.csv")
ensure_dir(dirname(out_file))
fwrite(results_dt, out_file)

cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(results_dt)))

# --- Summary ---
cat("\n========== Computation Cost Summary (mean sec) ==========\n")
dt_wide <- dcast(results_dt, test_name + variant ~ n,
                 value.var = "mean_time_sec")
print(dt_wide)
