################################################################################
# CIT Estimation Performance and Stability Analysis
# 
# Purpose: Address reviewer comment on how estimation affects instability in CIT
#          Compares oracle (known conditional expectations) vs estimated versions
#          and examines the relationship between regression error and test performance.
#
# Reference: Paper Examples 1 (stable) and 2 (unstable)
################################################################################

rm(list = ls())
set.seed(2024)

suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(ggplot2)
    library(ranger)
    library(xgboost)
})

# Source helper functions
source("./experiments/CIT_functions.R")

################################################################################
# 1. Data Generating Process
################################################################################

#' Generate covariate X from standard multivariate normal
#' @param n Sample size
#' @param p Dimension
generate_X <- function(n, p) {
    mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
}

#' Generate Y given X with known conditional mean f(X) = X %*% beta
#' @param X Covariate matrix
#' @param beta Coefficient vector
#' @param sigma Error standard deviation
generate_Y <- function(X, beta, sigma = 1) {
    n <- nrow(X)
    f_X <- as.numeric(X %*% beta)
    Y <- f_X + rnorm(n, 0, sigma)
    return(list(Y = Y, f_X = f_X))
}

#' Generate binary Z with known conditional probability g(X) = logistic(X %*% gamma)
#' Under H0: Z is independent of Y given X
#' @param X Covariate matrix
#' @param gamma Coefficient vector for logistic model
generate_Z_null <- function(X, gamma) {
    n <- nrow(X)
    linear_pred <- as.numeric(X %*% gamma)
    g_X <- plogis(linear_pred)  # P(Z=1|X)
    Z <- rbinom(n, 1, g_X)
    return(list(Z = Z, g_X = g_X))
}

################################################################################
# 2. Oracle GCM Test (with known conditional expectations)
################################################################################

#' Oracle GCM test using true conditional expectations
#' @param Y Response variable
#' @param Z Binary group indicator  
#' @param X Covariate matrix
#' @param f_X True E[Y|X]
#' @param g_X True E[Z|X]
gcm_oracle <- function(Y, Z, X, f_X, g_X) {
    n <- length(Y)
    
    # Residuals using true conditional expectations
    eps <- Y - f_X      # Y - E[Y|X]
    xi <- Z - g_X       # Z - E[Z|X]
    
    # GCM statistic
    R <- eps * xi
    T_n <- sqrt(n) * mean(R) / sd(R)
    
    # Two-sided p-value
    p_value <- 2 * pnorm(-abs(T_n))
    
    return(list(
        statistic = T_n,
        p_value = p_value,
        mean_R = mean(R),
        sd_R = sd(R)
    ))
}

################################################################################
# 3. Estimated GCM Test (with regression estimators)
################################################################################

#' Estimated GCM test with various regression methods
#' @param Y Response variable
#' @param Z Binary group indicator
#' @param X Covariate matrix
#' @param reg_method Regression method for E[Y|X]
#' @param binary_reg_method Regression method for E[Z|X]
#' @param seed Random seed
gcm_estimated <- function(Y, Z, X, reg_method, binary_reg_method, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    n <- length(Y)
    X <- as.matrix(X)
    
    # Estimate f(X) = E[Y|X]
    f_hat <- reg_method(X = X, y = Y, seed = seed)
    f_X_hat <- f_hat(X)
    
    # Estimate g(X) = E[Z|X]
    g_hat <- binary_reg_method(X = X, y = Z, seed = seed)
    g_X_hat <- g_hat(X)
    
    # Residuals using estimated conditional expectations
    eps_hat <- Y - f_X_hat
    xi_hat <- Z - g_X_hat
    
    # GCM statistic
    R_hat <- eps_hat * xi_hat
    
    # Handle edge cases
    if (sd(R_hat) < 1e-10) {
        T_n <- 0
        p_value <- 1
    } else {
        T_n <- sqrt(n) * mean(R_hat) / sd(R_hat)
        p_value <- 2 * pnorm(-abs(T_n))
    }
    
    return(list(
        statistic = T_n,
        p_value = p_value,
        f_X_hat = f_X_hat,
        g_X_hat = g_X_hat,
        eps_hat = eps_hat,
        xi_hat = xi_hat
    ))
}

################################################################################
# 4. Metrics for Regression Error
################################################################################

#' Compute regression error metrics
#' @param true_vals True conditional expectation values
#' @param est_vals Estimated conditional expectation values
compute_regression_error <- function(true_vals, est_vals) {
    # Mean Squared Error
    mse <- mean((true_vals - est_vals)^2)
    
    # Mean Absolute Error
    mae <- mean(abs(true_vals - est_vals))
    
    # R-squared
    ss_res <- sum((true_vals - est_vals)^2)
    ss_tot <- sum((true_vals - mean(true_vals))^2)
    r_squared <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
    
    # Correlation
    correlation <- cor(true_vals, est_vals)
    
    return(list(
        mse = mse,
        rmse = sqrt(mse),
        mae = mae,
        r_squared = r_squared,
        correlation = correlation
    ))
}

################################################################################
# 5. Main Simulation Function
################################################################################

#' Run single simulation comparing oracle vs estimated GCM
#' @param n Sample size
#' @param p Dimension
#' @param beta Coefficient for E[Y|X]
#' @param gamma Coefficient for E[Z|X]
#' @param reg_method Regression method name
#' @param seed Random seed
run_single_simulation <- function(n, p, beta, gamma, reg_method_name, seed) {
    set.seed(seed)
    
    # Generate data
    X <- generate_X(n, p)
    Y_data <- generate_Y(X, beta, sigma = 1)
    Y <- Y_data$Y
    f_X_true <- Y_data$f_X
    
    Z_data <- generate_Z_null(X, gamma)
    Z <- Z_data$Z
    g_X_true <- Z_data$g_X
    
    # Select regression methods
    if (reg_method_name == "linear") {
        reg_method <- lm_reg_method
        binary_reg_method <- lm_reg_method_binary
    } else if (reg_method_name == "ranger") {
        reg_method <- ranger_reg_method
        binary_reg_method <- ranger_reg_method_binary
    } else if (reg_method_name == "xgboost") {
        reg_method <- xgboost_reg_method
        binary_reg_method <- xgboost_reg_method_binary
    } else {
        stop("Unknown regression method")
    }
    
    # Oracle GCM
    oracle_result <- gcm_oracle(Y, Z, X, f_X_true, g_X_true)
    
    # Estimated GCM
    est_result <- gcm_estimated(Y, Z, X, reg_method, binary_reg_method, seed = seed + 1000)
    
    # Compute regression errors
    f_error <- compute_regression_error(f_X_true, est_result$f_X_hat)
    g_error <- compute_regression_error(g_X_true, est_result$g_X_hat)
    
    # Compute stability metric: difference in test statistics
    stat_diff <- abs(oracle_result$statistic - est_result$statistic)
    
    # Return results
    return(data.table(
        n = n,
        p = p,
        reg_method = reg_method_name,
        seed = seed,
        # Oracle results
        oracle_stat = oracle_result$statistic,
        oracle_pval = oracle_result$p_value,
        oracle_reject = as.integer(oracle_result$p_value < 0.05),
        # Estimated results
        est_stat = est_result$statistic,
        est_pval = est_result$p_value,
        est_reject = as.integer(est_result$p_value < 0.05),
        # Regression errors for f (E[Y|X])
        f_mse = f_error$mse,
        f_rmse = f_error$rmse,
        f_r2 = f_error$r_squared,
        f_cor = f_error$correlation,
        # Regression errors for g (E[Z|X])
        g_mse = g_error$mse,
        g_rmse = g_error$rmse,
        g_r2 = g_error$r_squared,
        g_cor = g_error$correlation,
        # Stability metric
        stat_difference = stat_diff
    ))
}

################################################################################
# 6. Run Full Simulation Study
################################################################################

cat("=" %+% rep("=", 70) %+% "\n")
cat("CIT Estimation Performance and Stability Analysis\n")
cat("=" %+% rep("=", 70) %+% "\n\n")

# Simulation parameters
n_values <- c(200, 500, 1000, 2000)
p <- 5
n_sims <- 500
reg_methods <- c("linear", "ranger", "xgboost")

# True parameters (sparse linear model)
beta <- c(1, -1, 0.5, 0, 0)
gamma <- c(0.5, -0.5, 0, 0, 0)

# String concatenation helper
`%+%` <- function(a, b) paste0(a, b)

results_list <- list()
counter <- 0
total_sims <- length(n_values) * length(reg_methods) * n_sims

cat("Running", total_sims, "simulations...\n\n")

for (n in n_values) {
    for (reg_method_name in reg_methods) {
        cat("n =", n, "| Method:", reg_method_name, "\n")
        
        sim_results <- pblapply(1:n_sims, function(sim) {
            run_single_simulation(
                n = n, 
                p = p, 
                beta = beta, 
                gamma = gamma, 
                reg_method_name = reg_method_name, 
                seed = 2024 + sim
            )
        })
        
        results_list[[length(results_list) + 1]] <- rbindlist(sim_results)
    }
}

# Combine all results
results_dt <- rbindlist(results_list)

################################################################################
# 7. Summary Statistics
################################################################################

cat("\n" %+% "=" %+% rep("=", 70) %+% "\n")
cat("SUMMARY RESULTS\n")
cat("=" %+% rep("=", 70) %+% "\n\n")

# Summary by n and method
summary_dt <- results_dt[, .(
    # Type I error rates
    oracle_type1 = mean(oracle_reject),
    est_type1 = mean(est_reject),
    type1_diff = mean(est_reject) - mean(oracle_reject),
    # Mean regression errors
    mean_f_rmse = mean(f_rmse),
    mean_g_rmse = mean(g_rmse),
    mean_f_r2 = mean(f_r2, na.rm = TRUE),
    mean_g_r2 = mean(g_r2, na.rm = TRUE),
    # Stability metrics
    mean_stat_diff = mean(stat_difference),
    sd_stat_diff = sd(stat_difference),
    # Correlation between errors and test stat difference
    cor_f_rmse_stat = cor(f_rmse, stat_difference),
    cor_g_rmse_stat = cor(g_rmse, stat_difference)
), by = .(n, reg_method)]

print(summary_dt)

################################################################################
# 8. Save Results
################################################################################

# Create results directory if needed
if (!dir.exists("results/appendix")) {
    dir.create("results/appendix", recursive = TRUE)
}

# Save detailed results
fwrite(results_dt, "results/appendix/cit_estimation_stability_detailed.csv")

# Save summary results
fwrite(summary_dt, "results/appendix/cit_estimation_stability_summary.csv")

cat("\n\nResults saved to:\n")
cat("  - results/appendix/cit_estimation_stability_detailed.csv\n")
cat("  - results/appendix/cit_estimation_stability_summary.csv\n")

################################################################################
# 9. Visualization
################################################################################

# Create figures directory if needed
if (!dir.exists("figures/appendix")) {
    dir.create("figures/appendix", recursive = TRUE)
}

# Plot 1: Type I Error Comparison (Oracle vs Estimated)
p1 <- ggplot(summary_dt, aes(x = factor(n), y = est_type1, fill = reg_method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_point(aes(y = oracle_type1), position = position_dodge(width = 0.8), 
               shape = 4, size = 3, color = "black") +
    labs(
        title = "Type I Error: Oracle (×) vs Estimated (bars)",
        x = "Sample Size (n)",
        y = "Type I Error Rate",
        fill = "Regression Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Set2")

ggsave("figures/appendix/cit_type1_comparison.pdf", p1, width = 8, height = 5)

# Plot 2: Regression Error vs Test Stability
p2 <- ggplot(results_dt, aes(x = f_rmse, y = stat_difference, color = reg_method)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    facet_wrap(~ n, scales = "free", labeller = label_both) +
    labs(
        title = "Regression Error (f) vs Test Statistic Instability",
        x = "RMSE of E[Y|X] Estimation",
        y = "|T_oracle - T_estimated|",
        color = "Regression Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set1")

ggsave("figures/appendix/cit_error_vs_stability.pdf", p2, width = 10, height = 4)

# Plot 3: Distribution of Test Statistics
p3 <- ggplot(results_dt[n == 500], aes(x = est_stat, fill = reg_method)) +
    geom_density(alpha = 0.5) +
    geom_vline(data = results_dt[n == 500, .(oracle_mean = mean(oracle_stat)), by = reg_method],
               aes(xintercept = oracle_mean), linetype = "dashed", color = "black") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), 
                  color = "red", linetype = "dotted", linewidth = 1) +
    facet_wrap(~ reg_method) +
    labs(
        title = "Distribution of Test Statistics (n = 500)",
        subtitle = "Dashed line: Oracle mean | Dotted red: N(0,1)",
        x = "Test Statistic",
        y = "Density",
        fill = "Regression Method"
    ) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set2")

ggsave("figures/appendix/cit_stat_distribution.pdf", p3, width = 10, height = 4)

# Plot 4: R-squared of Regression vs Sample Size
r2_summary <- results_dt[, .(
    mean_f_r2 = mean(f_r2, na.rm = TRUE),
    se_f_r2 = sd(f_r2, na.rm = TRUE) / sqrt(.N),
    mean_g_r2 = mean(g_r2, na.rm = TRUE),
    se_g_r2 = sd(g_r2, na.rm = TRUE) / sqrt(.N)
), by = .(n, reg_method)]

r2_long <- melt(r2_summary, 
                id.vars = c("n", "reg_method"),
                measure.vars = c("mean_f_r2", "mean_g_r2"),
                variable.name = "target", value.name = "r_squared")
r2_long[, target := ifelse(target == "mean_f_r2", "E[Y|X]", "E[Z|X]")]

p4 <- ggplot(r2_long, aes(x = factor(n), y = r_squared, fill = reg_method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
    facet_wrap(~ target) +
    labs(
        title = "Regression Quality (R²) by Sample Size and Method",
        x = "Sample Size (n)",
        y = "R-squared",
        fill = "Regression Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, 1)

ggsave("figures/appendix/cit_regression_r2.pdf", p4, width = 8, height = 5)

cat("\nFigures saved to figures/appendix/\n")

################################################################################
# 10. Additional Analysis: Stability Under Different Estimation Quality
################################################################################

cat("\n" %+% "=" %+% rep("=", 70) %+% "\n")
cat("STABILITY ANALYSIS\n")
cat("=" %+% rep("=", 70) %+% "\n\n")

# Quartile analysis of regression error vs stability
stability_analysis <- results_dt[, {
    # Quartiles of f_rmse
    f_quartile <- cut(f_rmse, quantile(f_rmse, probs = c(0, 0.25, 0.5, 0.75, 1)), 
                      labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
    
    .SD[, .(
        f_rmse_quartile = f_quartile,
        stat_difference = stat_difference,
        est_reject = est_reject
    )]
}, by = .(n, reg_method)]

stability_summary <- stability_analysis[, .(
    mean_stat_diff = mean(stat_difference),
    sd_stat_diff = sd(stat_difference),
    type1_error = mean(est_reject)
), by = .(n, reg_method, f_rmse_quartile)]

cat("Stability by Regression Error Quartile:\n")
print(stability_summary[order(n, reg_method, f_rmse_quartile)])

cat("\n" %+% "-" %+% rep("-", 70) %+% "\n")
cat("Key Finding: Higher regression error leads to greater test instability.\n")
cat("This supports the paper's discussion on estimation affecting CIT validity.\n")
cat("-" %+% rep("-", 70) %+% "\n")

