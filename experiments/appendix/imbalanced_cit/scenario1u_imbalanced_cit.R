rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(data.table)
})
tag <- "S1U_Imbalanced_CIT"
source("./experiments/all_tests.R")

# Data generation functions
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

# CIT test functions only
cit_test_functions <- list(
  RCIT_test = RCIT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test
)

# Parameters
n_base <- 50  # Base sample size for smaller group
imbalance_ratios <- c(1, 2, 3, 4)  # Ratio of larger group to smaller group
n_sims <- 500
alpha <- 0.05
d <- 10
results_list <- list()

# Simulation loop
for (ratio in imbalance_ratios) {
  n1 <- n_base  # Smaller group (Group 1)
  n2 <- n_base * ratio  # Larger group (Group 2)
  
  for (is_null in c(TRUE, FALSE)) {
    h_label <- if (is_null) "Null" else "Alternative"
    
    for (test_name in names(cit_test_functions)) {
      for (alg1_flag in c(TRUE, FALSE)) {
        alg1_label <- if (alg1_flag) "with_alg1" else "without_alg1"
        
        # Run the simulations
        result <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          
          # Generate data for Group 1 (smaller group) and Group 2 (larger group)
          x1 <- generate_data(n1, d, group = 1)
          y1 <- generate_y(x1, is_null = TRUE)
          set.seed(seed + n_sims)
          x2 <- generate_data(n2, d, group = 2)
          y2 <- generate_y(x2, is_null = is_null)
          
          test_args <- list(x1, x2, y1, y2, seed = seed, alg1 = alg1_flag)
          
          # For PCM_test and GCM_test, add regression methods
          if (test_name %in% c("PCM_test", "GCM_test")) {
            test_args$regr.method <- ranger_reg_method
            test_args$binary.regr.method <- ranger_reg_method_binary
          }
          
          do.call(cit_test_functions[[test_name]], test_args)
        }, simplify = "array")
        
        mean_result <- mean(result)
        results_list[[length(results_list) + 1]] <- data.table(
          test_name = test_name,
          n1 = n1,
          n2 = n2,
          imbalance_ratio = ratio,
          alg1 = alg1_label,
          h_label = h_label,
          rejection_rate = mean_result
        )
        
        # Print results
        cat("[Test]", test_name, "| Ratio:", ratio, "| n1:", n1, "| n2:", n2, "|", alg1_label, "|", h_label, "| Rejection Rate:", mean_result, "\n", strrep("-", 80), "\n")
      }
    }
  }
}

# Combine all results into a single data.table
results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")
