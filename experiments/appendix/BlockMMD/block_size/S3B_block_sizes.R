rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(tmvtnorm)
})
tag <- "S3B_MMDb_blocksize_ablation"
source("./experiments/all_tests.R")

# Bounded
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

generate_y <- function(x, is_null) {
    n <- nrow(x)
    # List of possible transformations
    transformations <- list(
        function(z) z,
        function(z) z^2,
        function(z) z^3,
        function(z) sin(z),
        function(z) tanh(z)
    )
    # Randomly choose a transformation
    random_transformation <- sample(transformations, 1)[[1]]
    
    if (is_null) {
        y <- cos(rowSums(x) + 2 * rnorm(n))
    } else {
        y<- random_transformation(rowSums(x) + 2 * rnorm(n))
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
block_size_multipliers <- c(0.5, 1, 1.5, 2)  # as multiples of sqrt(n)
n_sims <- 500
alpha <- 0.05
d <- 10
results_list <- list()

# Simulation loop
for (n in n_values) {
    B_defaults <- floor(sqrt(n))
    
    for (bmult in block_size_multipliers) {
        B_size <- max(2, floor(bmult * B_defaults))
        
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            
            for (test_type in c("DRT", "CIT")) {
                test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
                
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
                        
                        test_args <- list(x1, x2, y1, y2, B_size = B_size, seed = seed)
                        do.call(test_functions[[test_name]], test_args)
                    }, simplify = "array")
                    
                    mean_result <- mean(result)
                    results_list[[length(results_list) + 1]] <- data.table(
                        test_type = test_type,
                        test_name = test_name,
                        n = n,
                        block_size = B_size,
                        h_label = h_label,
                        rejection_rate = mean_result
                    )
                    
                    cat("[Test]", test_name, "| n:", n, "| B:", B_size, "|", h_label,
                        "| Rejection Rate:", mean_result, "\n", strrep("-", 80), "\n")
                }
            }
        }
    }
}

# Combine and save results
results_dt <- rbindlist(results_list)
dir.create("results/ablations/MMDb", showWarnings = FALSE, recursive = TRUE)
filename <- paste0("results/ablations/MMDb/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")