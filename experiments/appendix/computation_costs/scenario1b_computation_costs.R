rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(data.table)
  library(tmvtnorm)
})
tag <- "scenario1b"
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
  f0 <- x %*% c(1, -1, -1, 1, rep(0, dim(x)[2] - 4))
  mean_shift <- if (is_null) 0 else .5
  y <- f0 + epsilon + mean_shift
  return(y)
}

# Test functions
drt_test_functions <- list(
  LinearMMD_test = LinearMMD_test,
  CV_LinearMMD_test = CV_LinearMMD_test,
  BlockMMD_test = BlockMMD_test,
  bootstrap_MMD_test = bootstrap_MMD_test,
  CLF_test = CLF_test,
  CV_CLF_test = CV_CLF_test,
  CP_test = CP_test,
  debiased_test = debiased_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test,
  KCI_test = KCI_test
)

# Regression methods for CIT tests
regression_methods <- list(
  lm = list(regr.method = lm_reg_method, binary.regr.method = lm_reg_method_binary),
  ranger = list(regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary),
  xgboost = list(regr.method = xgboost_reg_method, binary.regr.method = xgboost_reg_method_binary)
)

# Parameters
n_values <- c(200, 500, 1000, 2000)
n_sims <- 50  # Reduced number of simulations for computational cost measurement
d <- 10
results_list <- list()

# Only test under null hypothesis
is_null <- TRUE

# Simulation loop
for (n in n_values) {
  for (test_type in c("DRT", "CIT")) {
    test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
    
    for (test_name in names(test_functions)) {
      
      if (test_type == "DRT") {
        # DRT tests: no regression method variation
        cat("[Measuring]", test_name, "| n:", n, "\n")
        
        time_results <- numeric(n_sims)
        for (sim in 1:n_sims) {
          seed <- 1203 + sim
          set.seed(seed)
          
          # Generate data for Group 1 and 2
          x1 <- generate_data(n, d, group = 1)
          y1 <- generate_y(x1, is_null = TRUE)
          set.seed(seed + n_sims)
          x2 <- generate_data(n, d, group = 2)
          y2 <- generate_y(x2, is_null = is_null)
          
          test_args <- list(x1, x2, y1, y2, seed = seed)
          
          # Measure execution time
          time_result <- system.time({
            result <- do.call(test_functions[[test_name]], test_args)
          })
          time_results[sim] <- time_result["elapsed"]
        }
        
        avg_time <- mean(time_results)
        std_time <- sd(time_results)
        
        results_list[[length(results_list) + 1]] <- data.table(
          test_type = test_type,
          test_name = test_name,
          n = n,
          regression_method = NA_character_,
          estimator = NA_character_,
          avg_time_seconds = avg_time,
          std_time_seconds = std_time
        )
        
        cat("  Avg time:", avg_time, "seconds | SD:", std_time, "\n", strrep("-", 80), "\n")
        
      } else {
        # CIT tests
        if (test_name == "RCIT_test") {
          # RCIT: no regression method variation
          cat("[Measuring]", test_name, "| n:", n, "\n")
          
          time_results <- numeric(n_sims)
          for (sim in 1:n_sims) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # Generate data for Group 1 and 2
            x1 <- generate_data(n, d, group = 1)
            y1 <- generate_y(x1, is_null = TRUE)
            set.seed(seed + n_sims)
            x2 <- generate_data(n, d, group = 2)
            y2 <- generate_y(x2, is_null = is_null)
            
            test_args <- list(x1, x2, y1, y2, seed = seed)
            
            # Measure execution time
            time_result <- system.time({
              result <- do.call(test_functions[[test_name]], test_args)
            })
            time_results[sim] <- time_result["elapsed"]
          }
          
          avg_time <- mean(time_results)
          std_time <- sd(time_results)
          
          results_list[[length(results_list) + 1]] <- data.table(
            test_type = test_type,
            test_name = test_name,
            n = n,
            regression_method = NA_character_,
            estimator = NA_character_,
            avg_time_seconds = avg_time,
            std_time_seconds = std_time
          )
          
          cat("  Avg time:", avg_time, "seconds | SD:", std_time, "\n", strrep("-", 80), "\n")
          
        } else {
          # Other CIT tests: test with different regression methods
          for (reg_name in names(regression_methods)) {
            cat("[Measuring]", test_name, "| n:", n, "| Regression:", reg_name, "\n")
            
            time_results <- numeric(n_sims)
            for (sim in 1:n_sims) {
              seed <- 1203 + sim
              set.seed(seed)
              
              # Generate data for Group 1 and 2
              x1 <- generate_data(n, d, group = 1)
              y1 <- generate_y(x1, is_null = TRUE)
              set.seed(seed + n_sims)
              x2 <- generate_data(n, d, group = 2)
              y2 <- generate_y(x2, is_null = is_null)
              
              test_args <- c(
                list(x1, x2, y1, y2, seed = seed),
                regression_methods[[reg_name]]
              )
              
              # Measure execution time
              time_result <- system.time({
                result <- do.call(test_functions[[test_name]], test_args)
              })
              time_results[sim] <- time_result["elapsed"]
            }
            
            avg_time <- mean(time_results)
            std_time <- sd(time_results)
            
            results_list[[length(results_list) + 1]] <- data.table(
              test_type = test_type,
              test_name = test_name,
              n = n,
              regression_method = reg_name,
              estimator = NA_character_,
              avg_time_seconds = avg_time,
              std_time_seconds = std_time
            )
            
            cat("  Avg time:", avg_time, "seconds | SD:", std_time, "\n", strrep("-", 80), "\n")
          }
        }
      }
    }
  }
}

# Combine all results into a single data.table
results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/computation_costs_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")

