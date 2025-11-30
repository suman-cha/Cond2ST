rm(list=ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(data.table)     
    library(tmvtnorm)
})
source("./experiments/all_tests.R")
tag <- "MMDb_S2B"

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

# Test functions
drt_test_functions <- list(
    BlockMMD_test = BlockMMD_test,
    CV_BlockMMD_test = CV_BlockMMD_test
)

cit_test_functions <- list()

# Parameters
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d <- 10
results_list <- list()


# Simulation loop
for (n in n_values) {
    for (is_null in c(TRUE, FALSE)) {
        h_label <- if (is_null) "Null" else "Alternative"
        
        for (test_type in c("DRT", "CIT")) {
            test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
            
            for (test_name in names(test_functions)) {
                
                # Run the simulations
                result <- pbapply::pbsapply(1:n_sims, function(sim) {
                    seed <- 1203 + sim
                    set.seed(seed)
                    
                    # Generate data for Group 1 and 2
                    x1 <- generate_data(n, d, group = 1)
                    y1 <- generate_y(x1, is_null = TRUE)
                    set.seed(seed + n_sims)
                    x2 <- generate_data(n, d, group = 2)
                    y2 <- generate_y(x2, is_null = is_null)
                    
                    test_args <- list(x1, x2, y1, y2, seed = seed)
                    
                    do.call(test_functions[[test_name]], test_args)
                }, simplify = "array")
                
                mean_result <- mean(result)
                results_list[[length(results_list) + 1]] <- data.table(
                    test_type = test_type,
                    test_name = test_name,
                    n = n,
                    h_label = h_label,
                    rejection_rate = mean_result
                )
                
                # Print results
                cat("[Test]", test_name, "| n:", n, "|", h_label, "| Rejection Rate:", mean_result, "\n", strrep("-", 80), "\n")
            }
        }
    }
}

# Combine all results into a single data.table
results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/ablations/MMDb/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")