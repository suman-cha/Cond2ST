#=============================================================================
# Scenario 5: Density-Ratio Stress Test (U: Unbounded, B: Bounded)
#   - X is non-Gaussian: product-Beta on (0,1)^p
#   - 5U: r(x)=f1(x)/f2(x) is unbounded (diverges near 0)
#   - 5B: r(x)=f1(x)/f2(x) is bounded on compact support
#=============================================================================
rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(pbapply)
    library(data.table)
})

# ---- Choose scenario: "U" for unbounded, "B" for bounded -------------------
scenario_flag <- "B"    # set to "B" for bounded
tag <- if (scenario_flag == "U") "S5U" else "S5B"

source("./experiments/all_tests.R")

# ---- Beta product sampler ---------------------------------------------------
# Draw an n x p matrix with independent Beta marginals per coordinate.
rprod_beta <- function(n, p, a, b) {
    # a, b can be scalars or length-p vectors
    a <- rep(a, length.out = p)
    b <- rep(b, length.out = p)
    out <- matrix(NA_real_, n, p)
    for (j in seq_len(p)) out[, j] <- rbeta(n, shape1 = a[j], shape2 = b[j])
    out
}

# ---- Group-specific X laws ensuring desired density-ratio geometry ----------
# Scenario 5U (Unbounded): group1 Beta(0.5,2), group2 Beta(2,2)  => r unbounded
# Scenario 5B (Bounded):   group1 Beta(2,4),   group2 Beta(3,5)  => r bounded
get_beta_params <- function(group, scenario = c("U","B"), p) {
    scenario <- match.arg(scenario)
    if (scenario == "U") {
        if (group == 1L) list(a = 4, b = 2) else list(a = 2, b = 2)
    } else {
        if (group == 1L) list(a = 4,   b = 4) else list(a = 5, b = 5)
    }
}

generate_data <- function(n, p, group, scenario = scenario_flag) {
    pars <- get_beta_params(group, scenario, p)
    rprod_beta(n, p, pars$a, pars$b)
}

# ---- Y|X mechanism (same across groups under H0; shifted under H1) ----------
# Nonlinear additive noise model; X is non-Gaussian (Beta product).
generate_y <- function(x, is_null = TRUE) {
    n <- nrow(x)
    # Nonlinear, smooth mean; bounded on (0,1)^p
    y_mean <- sin(2*pi*x[,1]) + 0.5 * log(1 + 10*x[,2]) + 0.3 * (x[,3] - 0.5)^2
    epsilon <- rnorm(n, 0, 1)  # homoscedastic noise
    mean_shift <- ifelse(is_null, 0, 0.5)  # shift only for group 2 under H1
    y_mean + epsilon + mean_shift
}

# ---- Test registries (as in your script) -----------------------------------
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
    GCM_test  = GCM_test
    # WGSC_test = WGSC_test,
    # KCI_test  = KCI_test
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
                    
                    # Group 1
                    x1 <- generate_data(n, d, group = 1L, scenario = scenario_flag)
                    y1 <- generate_y(x1, is_null = TRUE)
                    
                    # Group 2
                    set.seed(seed + n_sims)
                    x2 <- generate_data(n, d, group = 2L, scenario = scenario_flag)
                    y2 <- generate_y(x2, is_null = is_null)  # shift only under H1
                    
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
