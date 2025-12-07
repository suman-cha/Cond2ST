#=============================================================================
# Scenario 5: Density-Ratio Stress Test (U: Unbounded, B: Bounded by truncation)
#   - X: product-Beta on (0,1)^p
#   - 5U: group1 Beta(0.5,2), group2 Beta(2,2)  -> r(x) unbounded near 0
#   - 5B: SAME Betas but truncate support to [eps, 1 - eps] -> r(x) bounded
#=============================================================================
rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(pbapply)
    library(data.table)
})

# ---- Choose scenario: "U" for unbounded, "B" for bounded-by-truncation -----
scenario_flag <- "U"    # set to "B" for truncation-bounded
tag <- if (scenario_flag == "U") "S5U" else "S5B"

source("./experiments/all_tests.R")

# ---- Beta product sampler ---------------------------------------------------
rprod_beta <- function(n, p, a, b) {
    a <- rep(a, length.out = p)
    b <- rep(b, length.out = p)
    out <- matrix(NA_real_, n, p)
    for (j in seq_len(p)) out[, j] <- rbeta(n, shape1 = a[j], shape2 = b[j])
    out
}

# ---- Truncated Beta product sampler on [lower, upper] ----------------------
# Rejection sampling: draw from Beta on (0,1), keep if in [lower, upper]
rprod_beta_trunc_interval <- function(n, p, a, b, lower = 0.02, upper = 0.98) {
    stopifnot(0 < lower, lower < upper, upper < 1)
    a <- rep(a, length.out = p)
    b <- rep(b, length.out = p)
    acc <- NULL
    while (is.null(acc) || nrow(acc) < n) {
        m <- max(2*n, 2000L)  # batch size for efficiency
        X <- rprod_beta(m, p, a, b)
        keep <- apply(X >= lower & X <= upper, 1L, all)
        if (any(keep)) {
            Xk <- X[keep, , drop = FALSE]
            acc <- if (is.null(acc)) Xk else rbind(acc, Xk)
        }
    }
    acc[seq_len(n), , drop = FALSE]
}

# ---- Group-specific X laws (same Beta params; truncation only for 5B) -------
# Intentional: use the SAME (0.5,2) vs (2,2) in both 5U and 5B.
#   - 5U: untruncated -> r(x) diverges near 0
#   - 5B: truncated to [eps, 1-eps] -> r(x) bounded on compact set
get_beta_params <- function(group) {
    if (group == 1L) list(a = 0.5, b = 2) else list(a = 2, b = 2)
}

# truncation edge (can tune: 0.01~0.05; 0.02 works well)
eps_trunc <- 0.1

generate_data <- function(n, p, group, scenario = scenario_flag) {
    pars <- get_beta_params(group)
    if (scenario == "U") {
        rprod_beta(n, p, pars$a, pars$b)
    } else {
        rprod_beta_trunc_interval(n, p, pars$a, pars$b, lower = eps_trunc, upper = 1 - eps_trunc)
    }
}

# ---- Y|X mechanism (same across groups under H0; shifted under H1) ----------
generate_y <- function(x, is_null = TRUE) {
    n <- nrow(x)
    # smooth & bounded on (0,1)^p (and still fine on truncated support)
    y_mean <- sin(2*pi*x[,1]) + 0.5 * log(1 + 10*x[,2]) + 0.3 * (x[,3] - 0.5)^2
    epsilon <- rnorm(n, 0, 1)
    mean_shift <- ifelse(is_null, 0, 0.5)  # group 2 shift under H1
    y_mean + epsilon + mean_shift
}

# ---- Test registries --------------------------------------------------------
drt_test_functions <- list(
    LinearMMD_test = LinearMMD_test,
    CV_LinearMMD_test = CV_LinearMMD_test,
    CLF_test = CLF_test,
    CV_CLF_test = CV_CLF_test,
    CP_test = CP_test,
    debiased_test = debiased_test,
    BlockMMD_test = BlockMMD_test,
    CV_BlockMMD_test = CV_BlockMMD_test,
    bootstrap_MMD_test = bootstrap_MMD_test
)

cit_test_functions <- list(
    RCIT_test = RCIT_test,
    PCM_test  = PCM_test,
    GCM_test  = GCM_test,
    WGSC_test = WGSC_test,
    KCI_test  = KCI_test
)

# ---- Simulation settings ----------------------------------------------------
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
d <- 5
alpha <- 0.05
results_list <- list()

# ---- Simulation loop --------------------------------------------------------
for (n in n_values) {
    for (is_null in c(TRUE, FALSE)) {
        h_label <- if (is_null) "Null" else "Alternative"
        for (test_type in c("DRT", "CIT")) {
            test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
            for (test_name in names(test_functions)) {
                result <- pbapply::pbsapply(1:n_sims, function(sim) {
                    seed <- 1203 + sim
                    set.seed(seed)
                    x1 <- generate_data(n, d, group = 1L, scenario = scenario_flag)
                    y1 <- generate_y(x1, is_null = TRUE)
                    set.seed(seed + n_sims)
                    x2 <- generate_data(n, d, group = 2L, scenario = scenario_flag)
                    y2 <- generate_y(x2, is_null = is_null)
                    test_args <- list(x1, x2, y1, y2, seed = seed, alpha = alpha)
                    do.call(test_functions[[test_name]], test_args)
                }, simplify = "array")
                mean_result <- mean(result)
                results_list[[length(results_list) + 1]] <- data.table(
                    scenario = tag,
                    test_type = test_type,
                    test_name = test_name,
                    n = n,
                    h_label = h_label,
                    rejection_rate = mean_result
                )
                cat(sprintf("[%s][%s] %s | n: %d | %s | Rejection Rate: %.4f\n",
                            tag, test_type, test_name, n, h_label, mean_result))
                cat(strrep("-", 80), "\n")
            }
        }
    }
}

# ---- Save results -----------------------------------------------------------
results_dt <- rbindlist(results_list)
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("\n", strrep("=", 80), "\n")
cat("Results saved to", filename, "\n")
cat(strrep("=", 80), "\n")
