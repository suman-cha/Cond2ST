#=============================================================================
# Scenario 4: Linear Model with Heavy-tailed Covariates (t_ν)
#   - 4U: Unbounded density ratio via different tail indices (df_j)
#   - 4B: Bounded density ratio by truncating support to a compact cube
#   - Y|X: same as Scenario 1 (linear model with Gaussian noise)
#=============================================================================
rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(pbapply)
    library(data.table)
})

# ---------------------------- Scenario switch -------------------------------
scenario_flag <- "U"         # "U" -> S4U (unbounded),  "B" -> S4B (bounded)
tag <- if (scenario_flag == "U") "S4U" else "S4B"

source("./experiments/all_tests.R")

# ------------------------------ Parameters ----------------------------------
# Dimension d >= 4 for the specified beta and mu^(1); extra entries are zeros.
d <- 10

# Group-specific location shifts (mu) per your description:
mu1 <- c(1, 1, -1, -1, rep(0, max(0, d - 4)))
mu2 <- rep(0, d)

# Regression vector beta per your description:
beta <- c(1, -1, 1, -1, rep(0, max(0, d - 4)))

# For Scenario 4U (unbounded r): choose df with different tail indices.
#   df1 smaller (heavier tail) than df2  => r(x) = f1/f2 grows unbounded in the tails.
df_U_g1 <- 5  # >=3 as requested
df_U_g2 <- 4 # lighter-tailed

# For Scenario 4B (bounded r): truncate to a compact set; df can be any >=3.
df_B_g1 <- 5
df_B_g2 <- 4
# Truncation radius for each coordinate (compact cube [-R, R]^d):
R_trunc <- 3

# -------------------------- Samplers for X ----------------------------------
# Unbounded case: independent Student t per coordinate, different dfs per group.
rprod_t <- function(n, p, df) {
    # df can be scalar or length-p vector
    df <- rep(df, length.out = p)
    out <- matrix(NA_real_, n, p)
    for (j in seq_len(p)) out[, j] <- rt(n, df = df[j])
    out
}

# Bounded case: truncated independent t by rejection to the cube [-R_trunc, R_trunc]^p.
rprod_t_trunc <- function(n, p, df, R) {
    acc <- NULL
    # sample in batches to be efficient
    while (is.null(acc) || nrow(acc) < n) {
        m <- max(2*n, 1000L)
        Xb <- rprod_t(m, p, df)
        keep <- apply(abs(Xb) <= R, 1L, all)
        if (any(keep)) {
            X_keep <- Xb[keep, , drop = FALSE]
            if (is.null(acc)) acc <- X_keep else acc <- rbind(acc, X_keep)
        }
    }
    acc[seq_len(n), , drop = FALSE]
}

# Wrapper to generate X^(j) for group j with scenario switch and location mu^(j)
generate_data <- function(n, p, group, scenario = c("U","B")) {
    scenario <- match.arg(scenario)
    if (scenario == "U") {
        df <- if (group == 1L) df_U_g1 else df_U_g2
        base <- rprod_t(n, p, df)
    } else {
        df <- if (group == 1L) df_B_g1 else df_B_g2
        base <- rprod_t_trunc(n, p, df, R_trunc)
    }
    mu <- if (group == 1L) mu1 else mu2
    sweep(base, 2L, mu, FUN = "+")
}

# ---------------------------- Y | X (Scenario 1) ----------------------------
# y^{(j)} | x^{(j)} = δ^{(j)} + x^{(j)⊤} β + ε^{(j)},   ε^{(j)} ~ N(0,1)
# H0: δ^{(1)} = δ^{(2)} = 0
# H1: δ^{(1)} = 0, δ^{(2)} = 0.5
generate_y <- function(x, is_null = TRUE, group = 1L) {
    n <- nrow(x)
    delta <- if (is_null) 0 else if (group == 2L) 0.5 else 0
    drop(delta + x %*% beta + rt(n, 2))
}

# ---------------------------- Test registries -------------------------------
drt_test_functions <- list(
    LinearMMD_test      = LinearMMD_test,
    CV_LinearMMD_test   = CV_LinearMMD_test,
    CLF_test            = CLF_test,
    CV_CLF_test         = CV_CLF_test,
    CP_test             = CP_test,
    debiased_test       = debiased_test,
    BlockMMD_test       = BlockMMD_test,
    CV_BlockMMD_test    = CV_BlockMMD_test,
    bootstrap_MMD_test  = bootstrap_MMD_test
)

cit_test_functions <- list(
    RCIT_test = RCIT_test,
    PCM_test  = PCM_test,
    GCM_test  = GCM_test,
    WGSC_test = WGSC_test,
    KCI_test  = KCI_test
)

# ----------------------------- Simulation loop ------------------------------
n_values <- c(200, 500, 1000, 2000)
n_sims   <- 500
alpha    <- 0.05
results_list <- list()

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
                    y1 <- generate_y(x1, is_null = TRUE,  group = 1L)
                    
                    # Group 2
                    set.seed(seed + n_sims)
                    x2 <- generate_data(n, d, group = 2L, scenario = scenario_flag)
                    y2 <- generate_y(x2, is_null = is_null, group = 2L)
                    
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

# ------------------------------ Save results --------------------------------
results_dt <- rbindlist(results_list)
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("\n", strrep("=", 80), "\n")
cat("Results saved to", filename, "\n")
cat(strrep("=", 80), "\n")
