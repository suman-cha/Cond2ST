################################################################################
# Comment C.1 Analysis: Oracle vs. Subsampling Efficiency Comparison
# Quantifying finite-sample efficiency loss due to O(√n log(1/ε)) sample discarding
################################################################################

rm(list = ls())
set.seed(1203)

suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(tmvtnorm)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
})

source("./experiments/all_tests.R")

# ============================================================================
# 1. DATA GENERATION MECHANISMS (Adapted from Scenarios 1, 2, 3)
# ============================================================================

# ------------------------------
# Scenario 1: Mean Shift
# ------------------------------
generate_scenario1_data <- function(n, p, group) {
    mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
    sigma <- diag(1, p)
    X <- mvrnorm(n, mu = mu, Sigma = sigma)
    return(X)
}

generate_scenario1_y <- function(X, is_null = TRUE, sigma = 2) {
    n <- nrow(X)
    epsilon <- rt(n, df = sigma)
    f0 <- X %*% c(1, -1, -1, 1, rep(0, ncol(X) - 4))
    mean_shift <- if (is_null) 0 else 0.5
    Y <- f0 + epsilon + mean_shift
    return(Y)
}

# ------------------------------
# Scenario 2: Variance Structure
# ------------------------------
g_function <- function(X, rho) {
    X_adjusted <- X
    diag(X_adjusted) <- diag(X_adjusted) - 0.5
    norm_diff <- norm(X_adjusted, "F")
    return(10 + rho * exp(-norm_diff / 64))
}

generate_scenario2_data <- function(n, p, group) {
    mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
    sigma <- diag(1, p)
    X <- mvrnorm(n, mu = mu, Sigma = sigma)
    return(X)
}

generate_scenario2_y <- function(X, rho = 10, is_null = TRUE) {
    n <- nrow(X)
    p <- ncol(X)
    
    if (is_null) {
        beta <- rep(1, p)
    } else {
        beta <- c(rep(1, p - 1), 0)
    }
    
    beta <- matrix(beta, ncol = 1)
    mean_X <- X %*% beta
    var_X <- g_function(X, rho)
    
    if (is_null) {
        Y <- rnorm(n, mean = mean_X, sd = 10)
    } else {
        Y <- rnorm(n, mean = mean_X, sd = sqrt(var_X))
    }
    
    return(Y)
}

# ------------------------------
# Scenario 3: Nonlinear Transformation
# ------------------------------
generate_scenario3_data <- function(n, p, group) {
    mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
    sigma <- diag(1, p)
    X <- mvrnorm(n, mu = mu, Sigma = sigma)
    return(X)
}

generate_scenario3_y <- function(X, is_null) {
    n <- nrow(X)
    
    transformations <- list(
        function(z) z,
        function(z) z^2,
        function(z) z^3,
        function(z) sin(z),
        function(z) tanh(z)
    )
    
    random_transformation <- sample(transformations, 1)[[1]]
    
    if (is_null) {
        Y <- cos(rowSums(X) + 2 * rnorm(n))
    } else {
        Y <- random_transformation(rowSums(X) + 2 * rnorm(n))
    }
    
    return(Y)
}

# ============================================================================
# 2. TWO DATA GENERATION PARADIGMS
# ============================================================================

# ------------------------------
# Paradigm A: C2ST with Fixed Group Sizes (Algorithm 1 with subsampling)
# ------------------------------
generate_c2st_data <- function(n1, n2, p, scenario, is_null, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    # Generate fixed-size groups
    if (scenario == "S1") {
        X1 <- generate_scenario1_data(n1, p, group = 1)
        Y1 <- generate_scenario1_y(X1, is_null = TRUE)
        X2 <- generate_scenario2_data(n2, p, group = 2)
        Y2 <- generate_scenario1_y(X2, is_null = is_null)
    } else if (scenario == "S2") {
        X1 <- generate_scenario2_data(n1, p, group = 1)
        Y1 <- generate_scenario2_y(X1, is_null = TRUE)
        X2 <- generate_scenario2_data(n2, p, group = 2)
        Y2 <- generate_scenario2_y(X2, is_null = is_null)
    } else if (scenario == "S3") {
        X1 <- generate_scenario3_data(n1, p, group = 1)
        Y1 <- generate_scenario3_y(X1, is_null = TRUE)
        X2 <- generate_scenario3_data(n2, p, group = 2)
        Y2 <- generate_scenario3_y(X2, is_null = is_null)
    }
    
    return(list(X1 = X1, Y1 = Y1, X2 = X2, Y2 = Y2))
}

# ------------------------------
# Paradigm B: Oracle i.i.d. P_XYZ (No subsampling needed)
# ------------------------------
generate_oracle_data <- function(n, p, scenario, is_null, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    # Generate i.i.d. (X, Y, Z) ~ P_XYZ
    # Z ~ Bernoulli(0.5) for balanced groups
    Z <- rbinom(n, 1, 0.5) + 1  # Z ∈ {1, 2}
    
    X <- matrix(NA, nrow = n, ncol = p)
    Y <- numeric(n)
    
    for (i in 1:n) {
        if (scenario == "S1") {
            X[i, ] <- generate_scenario1_data(1, p, group = Z[i])
            if (Z[i] == 1) {
                Y[i] <- generate_scenario1_y(matrix(X[i, ], nrow = 1), is_null = TRUE)
            } else {
                Y[i] <- generate_scenario1_y(matrix(X[i, ], nrow = 1), is_null = is_null)
            }
        } else if (scenario == "S2") {
            X[i, ] <- generate_scenario2_data(1, p, group = Z[i])
            if (Z[i] == 1) {
                Y[i] <- generate_scenario2_y(matrix(X[i, ], nrow = 1), is_null = TRUE)
            } else {
                Y[i] <- generate_scenario2_y(matrix(X[i, ], nrow = 1), is_null = is_null)
            }
        } else if (scenario == "S3") {
            X[i, ] <- generate_scenario3_data(1, p, group = Z[i])
            if (Z[i] == 1) {
                Y[i] <- generate_scenario3_y(matrix(X[i, ], nrow = 1), is_null = TRUE)
            } else {
                Y[i] <- generate_scenario3_y(matrix(X[i, ], nrow = 1), is_null = is_null)
            }
        }
    }
    
    return(list(X = X, Y = Y, Z = Z))
}

# ============================================================================
# 3. TESTING FRAMEWORK
# ============================================================================

# CIT tests to compare
cit_test_functions <- list(
    RCIT_test = RCIT_test,
    PCM_test = PCM_test,
    GCM_test = GCM_test, 
    WGSC_test = WGSC_test
)

# Wrapper for C2ST setting (with subsampling)
run_c2st_test <- function(n1, n2, p, scenario, is_null, test_name, seed) {
    data <- generate_c2st_data(n1, n2, p, scenario, is_null, seed)
    
    test_func <- cit_test_functions[[test_name]]
    result <- test_func(data$X1, data$X2, data$Y1, data$Y2, seed = seed)
    
    return(result)
}

# Wrapper for Oracle setting (no subsampling)
run_oracle_test <- function(n, p, scenario, is_null, test_name, seed) {
    data <- generate_oracle_data(n, p, scenario, is_null, seed)
    
    # Split into groups based on Z
    idx1 <- which(data$Z == 1)
    idx2 <- which(data$Z == 2)
    
    X1 <- data$X[idx1, , drop = FALSE]
    Y1 <- data$Y[idx1]
    X2 <- data$X[idx2, , drop = FALSE]
    Y2 <- data$Y[idx2]
    
    test_func <- cit_test_functions[[test_name]]
    result <- test_func(X1, X2, Y1, Y2, seed = seed)
    
    return(result)
}

# ============================================================================
# 4. MAIN SIMULATION LOOP
# ============================================================================

run_comment_c1_experiment <- function() {
    
    # Parameters
    n_values <- c(200, 500, 1000, 2000)
    n_sims <- 500
    alpha <- 0.05
    p <- 10
    scenarios <- c("S1", "S2", "S3")
    
    results_list <- list()
    
    for (scenario in scenarios) {
        cat("\n========================================\n")
        cat("SCENARIO:", scenario, "\n")
        cat("========================================\n")
        
        for (n in n_values) {
            cat("\n--- Sample size n =", n, "---\n")
            
            for (is_null in c(TRUE, FALSE)) {
                h_label <- if (is_null) "Null" else "Alternative"
                
                for (test_name in names(cit_test_functions)) {
                    
                    # ==============================
                    # C2ST Setting (with subsampling)
                    # ==============================
                    cat("  [C2ST] Testing", test_name, "-", h_label, "...\n")
                    
                    n1 <- n2 <- n  # Balanced groups
                    
                    c2st_results <- pbsapply(1:n_sims, function(sim) {
                        seed <- 1203 + sim
                        tryCatch({
                            run_c2st_test(n1, n2, p, scenario, is_null, test_name, seed)
                        }, error = function(e) NA)
                    })
                    
                    c2st_rejection_rate <- mean(c2st_results, na.rm = TRUE)
                    
                    results_list[[length(results_list) + 1]] <- data.table(
                        scenario = scenario,
                        paradigm = "C2ST (Subsampling)",
                        test_name = test_name,
                        n = n,
                        h_label = h_label,
                        rejection_rate = c2st_rejection_rate,
                        n_na = sum(is.na(c2st_results))
                    )
                    
                    # ==============================
                    # Oracle Setting (no subsampling)
                    # ==============================
                    cat("  [Oracle] Testing", test_name, "-", h_label, "...\n")
                    
                    oracle_results <- pbsapply(1:n_sims, function(sim) {
                        seed <- 1203 + sim
                        tryCatch({
                            run_oracle_test(n, p, scenario, is_null, test_name, seed)
                        }, error = function(e) NA)
                    })
                    
                    oracle_rejection_rate <- mean(oracle_results, na.rm = TRUE)
                    
                    results_list[[length(results_list) + 1]] <- data.table(
                        scenario = scenario,
                        paradigm = "Oracle (i.i.d. P_XYZ)",
                        test_name = test_name,
                        n = n,
                        h_label = h_label,
                        rejection_rate = oracle_rejection_rate,
                        n_na = sum(is.na(oracle_results))
                    )
                    
                    cat("    C2ST Rejection:", round(c2st_rejection_rate, 4), 
                        "| Oracle Rejection:", round(oracle_rejection_rate, 4), "\n")
                }
            }
        }
    }
    
    return(rbindlist(results_list))
}

# ============================================================================
# 5. RUN EXPERIMENT
# ============================================================================

results <- run_comment_c1_experiment()

# Save results
fwrite(results, "results/comment_c1_oracle_vs_subsampling.csv")
cat("\n✓ Results saved to results/comment_c1_oracle_vs_subsampling.csv\n")

# ============================================================================
# 6. VISUALIZATION: Type I Error & Power Comparison
# ============================================================================

create_c1_plots <- function(results) {
    
    # Separate null and alternative
    null_results <- results[h_label == "Null"]
    alt_results <- results[h_label == "Alternative"]
    
    # ==============================
    # Plot 1: Type I Error Control
    # ==============================
    plot_type1 <- ggplot(null_results, aes(x = n, y = rejection_rate, 
                                           color = paradigm, shape = test_name)) +
        geom_point(size = 3.5, alpha = 0.8) +
        geom_line(aes(linetype = paradigm), size = 1) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 0.8) +
        facet_wrap(~scenario, nrow = 1) +
        scale_color_manual(
            values = c("C2ST (Subsampling)" = "#E41A1C", "Oracle (i.i.d. P_XYZ)" = "#377EB8")
        ) +
        scale_shape_manual(values = c(16, 17, 15)) +
        labs(
            x = "Sample Size (n)",
            y = "Type I Error Rate",
            title = "Comment C.1: Type I Error - Oracle vs. C2ST with Subsampling",
            subtitle = "Dashed red line indicates nominal level α = 0.05",
            color = "Data Generation",
            shape = "Test Method",
            linetype = "Data Generation"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 13, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 11, hjust = 0, color = "gray40"),
            panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 10, face = "bold"),
            strip.text = element_text(size = 10, face = "bold"),
            panel.grid.minor = element_blank()
        )
    
    # ==============================
    # Plot 2: Statistical Power
    # ==============================
    plot_power <- ggplot(alt_results, aes(x = n, y = rejection_rate, 
                                          color = paradigm, shape = test_name)) +
        geom_point(size = 3.5, alpha = 0.8) +
        geom_line(aes(linetype = paradigm), size = 1) +
        facet_wrap(~scenario, nrow = 1) +
        scale_color_manual(
            values = c("C2ST (Subsampling)" = "#E41A1C", "Oracle (i.i.d. P_XYZ)" = "#377EB8")
        ) +
        scale_shape_manual(values = c(16, 17, 15)) +
        labs(
            x = "Sample Size (n)",
            y = "Power (Rejection Rate)",
            title = "Comment C.1: Power Comparison - Oracle vs. C2ST with Subsampling",
            subtitle = "Quantifying finite-sample efficiency loss due to O(√n log(1/ε)) discarded samples",
            color = "Data Generation",
            shape = "Test Method",
            linetype = "Data Generation"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 13, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 11, hjust = 0, color = "gray40"),
            panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 10, face = "bold"),
            strip.text = element_text(size = 10, face = "bold"),
            panel.grid.minor = element_blank()
        )
    
    # ==============================
    # Plot 3: Efficiency Loss Quantification
    # ==============================
    # Calculate power difference
    power_diff <- merge(
        alt_results[paradigm == "Oracle (i.i.d. P_XYZ)", 
                    .(scenario, test_name, n, oracle_power = rejection_rate)],
        alt_results[paradigm == "C2ST (Subsampling)", 
                    .(scenario, test_name, n, c2st_power = rejection_rate)],
        by = c("scenario", "test_name", "n")
    )
    power_diff[, efficiency_loss := oracle_power - c2st_power]
    
    plot_loss <- ggplot(power_diff, aes(x = n, y = efficiency_loss, 
                                        color = test_name, shape = scenario)) +
        geom_point(size = 3.5, alpha = 0.8) +
        geom_line(aes(linetype = test_name), size = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
        facet_wrap(~scenario, nrow = 1) +
        scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
        scale_shape_manual(values = c(16, 17, 15)) +
        labs(
            x = "Sample Size (n)",
            y = "Power Loss (Oracle - C2ST)",
            title = "Comment C.1: Finite-Sample Efficiency Loss Due to Subsampling",
            subtitle = "Positive values indicate oracle advantage; values near 0 indicate negligible loss",
            color = "Test Method",
            shape = "Scenario",
            linetype = "Test Method"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 13, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 11, hjust = 0, color = "gray40"),
            panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 10, face = "bold"),
            strip.text = element_text(size = 10, face = "bold"),
            panel.grid.minor = element_blank()
        )
    
    return(list(plot_type1 = plot_type1, plot_power = plot_power, plot_loss = plot_loss))
}

# Generate plots
plots <- create_c1_plots(results)

# Save plots
ggsave("figures/c1_type1_error.pdf", plots$plot_type1, 
       width = 14, height = 5, dpi = 300)
ggsave("figures/c1_power_comparison.pdf", plots$plot_power, 
       width = 14, height = 5, dpi = 300)
ggsave("figures/c1_efficiency_loss.pdf", plots$plot_loss, 
       width = 14, height = 5, dpi = 300)

cat("✓ Publication-quality figures saved to figures/\n")

# ============================================================================
# 7. SUMMARY STATISTICS
# ============================================================================

cat("\n=== SUMMARY STATISTICS ===\n")

cat("\nType I Error Summary:\n")
print(results[h_label == "Null", .(
    mean_type1 = mean(rejection_rate, na.rm = TRUE),
    sd_type1 = sd(rejection_rate, na.rm = TRUE)
), by = .(scenario, paradigm, test_name)])

cat("\nPower Summary:\n")
print(results[h_label == "Alternative", .(
    mean_power = mean(rejection_rate, na.rm = TRUE),
    sd_power = sd(rejection_rate, na.rm = TRUE)
), by = .(scenario, paradigm, test_name)])

cat("\n✓ Comment C.1 analysis complete!\n")
