rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(tmvtnorm)
})

tag <- "S1B_MMDb_blocksize_splitting_ablation"
source("./experiments/all_tests.R")

# Function to generate covariates using truncated multivariate normal distributions
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

# Function to generate response variable
generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, 1, 1, rep(0, dim(x)[2] - 4))
    mean_shift <- if (is_null) 0 else .5
    y <- f0 + epsilon + mean_shift
    return(y)
}

# Modified test functions to handle custom splitting ratios
BlockMMD_test_custom <- function(x1, x2, y1, y2, B_size = NULL, prop_ratio = 0.5, 
                                 alpha = 0.05, bandwidths = c(1), est.method = 'LL', seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    stopifnot(length(y1) == length(y2))
    total_n <- length(y1)
    
    # Custom splitting: N samples for ratio estimation, (n-N) for statistics
    N <- ceiling(total_n * prop_ratio)  # ratio estimation samples
    n_stat <- total_n - N              # statistic computation samples
    
    # Split data
    x11 <- x1[1:N, , drop = FALSE]; x12 <- x1[(N+1):total_n, , drop = FALSE]
    x21 <- x2[1:N, , drop = FALSE]; x22 <- x2[(N+1):total_n, , drop = FALSE]
    y11 <- y1[1:N]; y12 <- y1[(N+1):total_n]
    y21 <- y2[1:N]; y22 <- y2[(N+1):total_n]
    
    # Estimate density ratio
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    r_X2 <- 1/ratios$g22.est
    
    if (is.null(B_size)) B_size <- floor(n_stat^(0.5))
    
    num_blocks <- floor(n_stat / B_size)
    if (num_blocks < 2) {
        warning("Test sample is too small for the given B_size, resulting in less than 2 blocks.")
        return(rep(NA, length(bandwidths)))
    }
    
    rejections <- numeric(length(bandwidths))
    for (i in seq_along(bandwidths)) {
        bw <- bandwidths[i]
        
        if (bw == "median") {
            h_x <- median.bandwidth(x12, x22)
            h_y <- median.bandwidth(matrix(y12), matrix(y22))
        } else {
            h_x <- as.numeric(bw)
            h_y <- as.numeric(bw)
        }
        
        test_stat <- MMDb(x12, x22, y12, y22, B_size,
                          h_x = h_x, h_y = h_y, r_X = r_X2, seed = seed)
        
        pval <- 1 - pnorm(test_stat)
        rejections[i] <- as.integer(pval < alpha)
    }
    return(rejections)
}

# Test functions
drt_test_functions <- list(
    BlockMMD_test_custom = BlockMMD_test_custom
)
cit_test_functions <- list()

# Parameters
n_values <- c(500, 1000, 2000, 4000)
splitting_ratios <- c(0.6, 0.7, 0.8, 0.9)  # N/n ratios (more samples for ratio estimation)
block_size_powers <- c(0.3, 0.4, 0.5, 0.6)  # Powers for (n-N)^power
n_sims <- 500
alpha <- 0.05
d <- 10
results_list <- list()

# Simulation loop
for (n in n_values) {
    for (prop_ratio in splitting_ratios) {
        N <- ceiling(n * prop_ratio)
        n_stat <- n - N
        
        for (block_power in block_size_powers) {
            B_size <- max(2, floor(n_stat^block_power))
            
            # Check constraint: B_size^(1+2*block_power) << N
            constraint_value <- B_size^(1 + 2 * block_power)
            constraint_ratio <- constraint_value / N
            
            # Skip if constraint is not satisfied (constraint_ratio should be << 1, say < 0.1)
            if (constraint_ratio >= 0.1) {
                cat("Skipping n=", n, ", prop_ratio=", prop_ratio, ", block_power=", block_power, 
                    ": constraint not satisfied (", round(constraint_ratio, 3), ")\n")
                next
            }
            
            for (is_null in c(TRUE, FALSE)) {
                h_label <- if (is_null) "Null" else "Alternative"
                
                for (test_type in c("DRT")) {  # Only DRT for this ablation
                    test_functions <- drt_test_functions
                    
                    for (test_name in names(test_functions)) {
                        
                        # Run the simulations
                        result <- pbapply::pbsapply(1:n_sims, function(sim) {
                            seed <- 1203 + sim
                            set.seed(seed)
                            
                            x1 <- generate_data(n, d, group = 1)
                            y1 <- generate_y(x1, is_null = TRUE)
                            set.seed(seed + n_sims)
                            x2 <- generate_data(n, d, group = 2)
                            y2 <- generate_y(x2, is_null = is_null)
                            
                            test_args <- list(x1, x2, y1, y2, B_size = B_size, 
                                              prop_ratio = prop_ratio, seed = seed)
                            do.call(test_functions[[test_name]], test_args)
                        }, simplify = "array")
                        
                        mean_result <- mean(result, na.rm = TRUE)
                        results_list[[length(results_list) + 1]] <- data.table(
                            test_type = test_type,
                            test_name = test_name,
                            n = n,
                            N_ratio_est = N,
                            n_stat = n_stat,
                            splitting_ratio = prop_ratio,
                            block_size = B_size,
                            block_power = block_power,
                            constraint_ratio = constraint_ratio,
                            h_label = h_label,
                            rejection_rate = mean_result
                        )
                        
                        cat("[Test]", test_name, "| n:", n, "| N:", N, "| n_stat:", n_stat, 
                            "| B:", B_size, "| power:", block_power, "| constraint:", 
                            round(constraint_ratio, 3), "|", h_label,
                            "| Rejection Rate:", round(mean_result, 3), "\n")
                    }
                }
            }
        }
        cat(strrep("-", 100), "\n")
    }
}

# Combine and save results
results_dt <- rbindlist(results_list)
dir.create("results/ablations/MMDb", showWarnings = FALSE, recursive = TRUE)
filename <- paste0("results/ablations/MMDb/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")

# Print summary of experimental conditions
cat("\nExperimental Summary:\n")
summary_dt <- results_dt[, .(
    constraint_ratio = mean(constraint_ratio),
    rejection_rate_null = mean(rejection_rate[h_label == "Null"]),
    rejection_rate_alt = mean(rejection_rate[h_label == "Alternative"])
), by = .(n, splitting_ratio, block_power, block_size)]

print(summary_dt)
