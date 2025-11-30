#' @title Comprehensive Analysis of BlockMMD Statistics Distribution
#' @description Analyzes BlockMMD and LinearMMD statistics distributions across 
#'              different splitting ratios and block sizes for comparative evaluation.

rm(list=ls())
library(tmvtnorm)
library(CVST)
library(data.table)
library(ggplot2)
library(gridExtra)
source("./experiments/utils.R") # contains MMDb, MMDl, estimate_r

# Function to generate covariates
generate_data <- function(n, p, group) {
    mu <- c(1, 1, -1, -1, rep(0, p - 4))
    sigma <- diag(1, p)
    lb <- rep(-.5, p)
    ub <- rep(.5, p)
    if (group == 1) {
        x <- rtmvnorm(n, mean = mu, sigma = sigma, lower = lb, upper = ub, algorithm = "gibbs")
    } else {
        x <- rtmvnorm(n, mean = rep(0, p), sigma = sigma, lower = lb, upper = ub, algorithm = "gibbs")
    }
    return(x)
}

# Function to generate response
generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, 1, 1, rep(0, dim(x)[2] - 4))
    mean_shift <- if (is_null) 0 else 0.5
    y <- f0 + epsilon + mean_shift
    return(y)
}

simulate_comprehensive <- function(n = 2000, n_sim = 500, seed = 1203, 
                                   hypothesis = c("null", "alternative"),
                                   prop_vals = c(0.3, 0.5, 0.7), 
                                   block_powers = c(0.3, 0.4, 0.5, 0.6),
                                   est_method = "LL") {
    
    hypothesis <- match.arg(hypothesis)
    d <- 10
    is_null_hyp <- (hypothesis == "null")
    title_main <- ifelse(is_null_hyp, "Null", "Alternative")
    
    # Initialize results storage
    all_results <- data.table()
    
    cat("Starting comprehensive simulation for", title_main, "hypothesis\n")
    cat("Total combinations:", length(prop_vals) * length(block_powers), "\n")
    
    # Loop through all combinations
    for (prop in prop_vals) {
        for (block_power in block_powers) {
            
            n_train <- ceiling(n * prop)
            n_test <- n - n_train
            B_size <- max(2, floor(n_test^block_power))
            
            # Check constraint: B_size should be reasonable for n_test
            if (B_size >= n_test / 2) {
                cat("Skipping prop =", prop, ", power =", block_power, 
                    ": B_size too large (", B_size, ">=", n_test/2, ")\n")
                next
            }
            
            cat(sprintf("Running: prop=%.1f, power=%.1f -> n_train=%d, n_test=%d, B_size=%d\n", 
                        prop, block_power, n_train, n_test, B_size))
            
            # Storage for this combination
            stats_block <- numeric(n_sim)
            stats_linear <- numeric(n_sim)
            
            # Run simulations
            for (i in seq_len(n_sim)) {
                set.seed(seed + i)
                x1 <- generate_data(n, d, group = 1)
                y1 <- generate_y(x1, is_null = TRUE)
                set.seed(seed + i + n_sim)
                x2 <- generate_data(n, d, group = 2)
                y2 <- generate_y(x2, is_null = is_null_hyp)
                
                # Split data
                test_idx <- (n_train + 1):n
                x11 <- x1[1:n_train, , drop=FALSE]; x12 <- x1[test_idx, , drop=FALSE]
                y11 <- y1[1:n_train]; y12 <- y1[test_idx]
                x21 <- x2[1:n_train, , drop=FALSE]; x22 <- x2[test_idx, , drop=FALSE]
                y21 <- y2[1:n_train]; y22 <- y2[test_idx]
                
                # Estimate density ratio
                ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est_method, seed+i)
                r_X <- 1 / ratios$g22.est
                
                h_x <- h_y <- 1
                
                # Calculate statistics
                tryCatch({
                    stats_block[i] <- MMDb(x12, x22, y12, y22, B_size, h_x, h_y, r_X, seed+i)
                    stats_linear[i] <- MMDl(x12, x22, y12, y22, h_x, h_y, r_X, seed+i)
                }, error = function(e) {
                    stats_block[i] <<- NA
                    stats_linear[i] <<- NA
                })
            }
            
            # Remove NAs and infinite values
            valid_block <- is.finite(stats_block)
            valid_linear <- is.finite(stats_linear)
            
            # Store results
            block_dt <- data.table(
                statistic = stats_block[valid_block],
                method = "BlockMMD",
                splitting_ratio = prop,
                block_power = block_power,
                block_size = B_size,
                n_train = n_train,
                n_test = n_test,
                hypothesis = title_main,
                sim_id = which(valid_block)
            )
            
            linear_dt <- data.table(
                statistic = stats_linear[valid_linear],
                method = "LinearMMD",
                splitting_ratio = prop,
                block_power = block_power,
                block_size = B_size,
                n_train = n_train,
                n_test = n_test,
                hypothesis = title_main,
                sim_id = which(valid_linear)
            )
            
            all_results <- rbind(all_results, block_dt, linear_dt)
            
            # Print summary statistics for this combination
            cat(sprintf("  BlockMMD: mean=%.3f, sd=%.3f, valid=%d/%d\n",
                        mean(stats_block[valid_block]), sd(stats_block[valid_block]), 
                        sum(valid_block), n_sim))
            cat(sprintf("  LinearMMD: mean=%.3f, sd=%.3f, valid=%d/%d\n",
                        mean(stats_linear[valid_linear]), sd(stats_linear[valid_linear]), 
                        sum(valid_linear), n_sim))
        }
    }
    
    return(all_results)
}

# Function to create comprehensive plots
create_distribution_plots <- function(results_dt, save_plots = TRUE) {
    
    # Plot 1: Density plots by splitting ratio and block size (BlockMMD only)
    p1 <- ggplot(results_dt[method == "BlockMMD"], 
                 aes(x = statistic, color = factor(block_size))) +
        geom_density(alpha = 0.7) +
        facet_grid(hypothesis ~ splitting_ratio, 
                   labeller = labeller(splitting_ratio = label_both)) +
        labs(title = "BlockMMD Statistics Distribution by Splitting Ratio and Block Size",
             x = "Statistic Value", y = "Density", 
             color = "Block Size") +
        theme_minimal() +
        theme(strip.text = element_text(size = 10))
    
    # Plot 2: Box plots comparing methods
    p2 <- ggplot(results_dt, aes(x = factor(splitting_ratio), y = statistic, 
                                 fill = method)) +
        geom_boxplot(alpha = 0.7) +
        facet_grid(hypothesis ~ block_size, 
                   labeller = labeller(block_size = label_both)) +
        labs(title = "Statistics Distribution: BlockMMD vs LinearMMD",
             x = "Splitting Ratio", y = "Statistic Value", 
             fill = "Method") +
        theme_minimal()
    
    # Plot 3: Q-Q plots for normality check
    p3 <- ggplot(results_dt[method == "BlockMMD"], aes(sample = statistic)) +
        stat_qq() + stat_qq_line(color = "red") +
        facet_grid(splitting_ratio ~ block_size, 
                   labeller = labeller(.default = label_both)) +
        labs(title = "Q-Q Plots: BlockMMD Statistics vs Normal Distribution") +
        theme_minimal()
    
    # Plot 4: Summary statistics heatmap
    summary_stats <- results_dt[, .(
        mean_stat = mean(statistic, na.rm = TRUE),
        sd_stat = sd(statistic, na.rm = TRUE),
        count = .N
    ), by = .(method, splitting_ratio, block_size, hypothesis)]
    
    p4 <- ggplot(summary_stats[method == "BlockMMD"], 
                 aes(x = factor(splitting_ratio), y = factor(block_size), 
                     fill = mean_stat)) +
        geom_tile() +
        geom_text(aes(label = round(mean_stat, 2)), color = "white", size = 3) +
        facet_wrap(~ hypothesis) +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0) +
        labs(title = "Mean BlockMMD Statistics Heatmap",
             x = "Splitting Ratio", y = "Block Size", 
             fill = "Mean Statistic") +
        theme_minimal()
    
    # Display plots
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    
    if (save_plots) {
        ggsave("blockMMD_density_by_params.png", p1, width = 12, height = 8)
        ggsave("blockMMD_vs_linearMMD_boxplot.png", p2, width = 12, height = 8)
        ggsave("blockMMD_qq_plots.png", p3, width = 12, height = 8)
        ggsave("blockMMD_mean_heatmap.png", p4, width = 10, height = 6)
        cat("Plots saved as PNG files\n")
    }
    
    return(list(density_plot = p1, boxplot = p2, qq_plot = p3, heatmap = p4))
}

# Function to print summary statistics
print_summary_analysis <- function(results_dt) {
    cat("\n" , rep("=", 80), "\n")
    cat("COMPREHENSIVE SUMMARY ANALYSIS\n")
    cat(rep("=", 80), "\n")
    
    # Overall summary by method and hypothesis
    overall_summary <- results_dt[, .(
        mean_stat = mean(statistic, na.rm = TRUE),
        sd_stat = sd(statistic, na.rm = TRUE),
        median_stat = median(statistic, na.rm = TRUE),
        q25 = quantile(statistic, 0.25, na.rm = TRUE),
        q75 = quantile(statistic, 0.75, na.rm = TRUE),
        count = .N
    ), by = .(method, hypothesis)]
    
    cat("Overall Statistics by Method and Hypothesis:\n")
    print(overall_summary)
    
    # Effect of splitting ratio on BlockMMD
    split_effect <- results_dt[method == "BlockMMD", .(
        mean_stat = mean(statistic, na.rm = TRUE),
        sd_stat = sd(statistic, na.rm = TRUE),
        count = .N
    ), by = .(splitting_ratio, hypothesis)]
    
    cat("\nEffect of Splitting Ratio on BlockMMD:\n")
    print(split_effect)
    
    # Effect of block size on BlockMMD
    block_effect <- results_dt[method == "BlockMMD", .(
        mean_stat = mean(statistic, na.rm = TRUE),
        sd_stat = sd(statistic, na.rm = TRUE),
        count = .N
    ), by = .(block_size, hypothesis)]
    
    cat("\nEffect of Block Size on BlockMMD:\n")
    print(block_effect)
}

# Main execution
cat("Starting comprehensive BlockMMD distribution analysis...\n")

# Run simulations for both hypotheses
results_null <- simulate_comprehensive(hypothesis = "null")
results_alt <- simulate_comprehensive(hypothesis = "alternative")

# Combine results
all_results <- rbind(results_null, results_alt)

# Create and display plots
plots <- create_distribution_plots(all_results, save_plots = TRUE)

# Print comprehensive summary
print_summary_analysis(all_results)

# Save results
fwrite(all_results, "blockMMD_distribution_analysis_results.csv")
cat("\nResults saved to: blockMMD_distribution_analysis_results.csv\n")

# Additional analysis: Type I error rates under null hypothesis
if (nrow(results_null) > 0) {
    type1_rates <- results_null[method == "BlockMMD", .(
        type1_rate = mean(statistic > qnorm(0.95), na.rm = TRUE),  # approximate 5% threshold
        type1_rate_01 = mean(statistic > qnorm(0.99), na.rm = TRUE)   # 1% threshold
    ), by = .(splitting_ratio, block_size)]
    
    cat("\nType I Error Rates (BlockMMD under null hypothesis):\n")
    cat("5% threshold (should be ≈ 0.05):\n")
    print(type1_rates[, .(splitting_ratio, block_size, type1_rate)])
}

cat("\nAnalysis complete!\n")
