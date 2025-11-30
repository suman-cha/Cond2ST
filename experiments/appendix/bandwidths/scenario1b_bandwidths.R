rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(tmvtnorm)
})

tag <- "S1B_bandwidth_comparison"
source("./experiments/all_tests.R")

generate_data <- function(n, p, group) {
    mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
    sigma <- diag(1, p)
    mvrnorm(n, mu = mu, Sigma = sigma)
}

generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, -1, 1, rep(0, dim(x)[2] - 4))
    mean_shift <- if (is_null) 0 else .5
    f0 + epsilon + mean_shift
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