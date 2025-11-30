rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(tmvtnorm)
})
tag <- "S2B_bandwidth_comparison"
source("./experiments/all_tests.R")

# Define the function g(Z)
g <- function(Z, rho) {
    Z_adjusted <- Z
    diag(Z_adjusted) <- diag(Z_adjusted) - 0.5
    norm_diff <- norm(Z_adjusted, "F")
    return(10 + rho * exp(-norm_diff/64))
}

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


generate_y <- function(x, rho=10, is_null=TRUE) {
    n <- nrow(x)
    p <- ncol(x)
    
    if (is_null) {
        beta <- c(rep(1, p))
    } else {
        beta <- c(rep(1, p-1),0)
    }
    
    beta <- matrix(beta, ncol = 1)  
    mean_X <- x %*% beta  
    var_X <- g(x, rho)
    
    if (is_null) {
        y <- rnorm(n, mean = mean_X, sd = 10)
    } else {
        y <- rnorm(n, mean = mean_X, sd = sqrt(var_X))
    }
    
    return(y)
}


# Parameters 
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d <- 10
bandwidths <- c(0.1, 0.5, 1.0, "median", 5, 10)
results_list <- list()

test_functions <- list(
    LinearMMD_test = LinearMMD_test,
    CV_LinearMMD_test = CV_LinearMMD_test
)

# Simulation loop
for(n in n_values) {
    for(is_null in c(TRUE, FALSE)) {
        h_label <- if(is_null) "Null" else "Alternative"
        
        for(test_name in names(test_functions)) {
            result_matrix <- pbapply::pbsapply(1:n_sims, function(sim) {
                seed <- 1203 + sim
                set.seed(seed)
                
                x1 <- generate_data(n, d, group = 1)
                y1 <- generate_y(x1, is_null = TRUE)
                set.seed(seed + n_sims)
                x2 <- generate_data(n, d, group = 2)
                y2 <- generate_y(x2, is_null = is_null)
                
                test_args <- list(
                    x1, x2, y1, y2, 
                    seed = seed,
                    bandwidths = bandwidths
                )
                do.call(test_functions[[test_name]], test_args)
            })
            
            bw_results <- data.table(
                bandwidth = sapply(bandwidths, as.character),
                rejection_rate = apply(result_matrix, 1, mean)
            )
            
            for(bw_idx in seq_along(bandwidths)) {
                bw <- bandwidths[bw_idx]
                bw_label <- if(is.character(bw)) bw else as.character(bw)
                mean_result <- mean(result_matrix[bw_idx, ])
                
                results_list[[length(results_list) + 1]] <- data.table(
                    test_type = "DRT",
                    test_name = test_name,
                    n = n,
                    h_label = h_label,
                    bandwidth = bw_label,
                    rejection_rate = mean_result
                )
            }
            
            cat(sprintf("[Test] %s | n: %d | %s\n", test_name, n, h_label))
            print(bw_results)
            cat(strrep("-", 80), "\n")
        }
    }
}


# Combine results and save
results_dt <- rbindlist(results_list)
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")