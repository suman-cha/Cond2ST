rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(pbapply)
  library(MASS)
  library(dplyr)
  library(data.table)
  library(parallel)
})
source("./experiments/all_tests.R")
tag <- "real_high_dim"

cur_wd <- getwd()
file_path <- file.path(cur_wd, "experiments", "real_examples", "data", "superconductivity.csv")
data <- read.csv(file_path)

X <- as.matrix(data[, !names(data) %in% c("critical_temp")])
Y <- data[, "critical_temp"]

# Remove missing data
X <- na.omit(X)
Y <- Y[complete.cases(Y)]

normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

X_norm <- apply(X, 2, normalize)
Y_norm <- normalize(Y)

sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
  if (is_x1) {
    # uniform sampling for X1
    X_idx <- sample(1:nrow(X), nrow(X)%/%2, replace = FALSE)
    X1 <- X[X_idx, , drop = FALSE]
    Y_subset <- Y[X_idx]
    
    # subsample from X1 to construct x1
    x_idx <- sample(1:nrow(X1), n, replace=FALSE)
    x <- X1[x_idx,,drop=FALSE]
  } else {
    # biased sampling for X2 based on normal distribution
    feature_to_bias <- X[, 1]  
    prob <- dnorm(feature_to_bias, 0, 1)
    prob <- prob / sum(prob)  
    X_idx <- sample(1:nrow(X), nrow(X)%/%2, replace = FALSE, prob = prob)
    X2 <- X[X_idx, , drop = FALSE]
    Y_subset <- Y[X_idx]  
    
    # subsample from X2 to construct x2
    x_idx <- sample(1:nrow(X2), n, replace=FALSE)
    x <- X2[x_idx,,drop=FALSE]
  }
  
  if (is_null) {
    # Null hypothesis: uniform sampling from Y values
    y <- sample(Y_subset, size=n, replace=FALSE)
  } else {
    # Alternative hypothesis: introduce bias in Y1 and Y2
    if (is_x1) {
      u <- dunif(Y_subset, 0, 1)
    } else {
      u <- exp(-Y_subset)
    }
    u <- u / sum(u)
    y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
  }
  
  return(list(x = x, y = y))
}

# Define test functions
drt_test_functions <- list(
  LinearMMD_test = LinearMMD_test,
  CLF_test = CLF_test,
  CP_test = CP_test,
  CV_LinearMMD_test = CV_LinearMMD_test,
  CV_CLF_test = CV_CLF_test,
  debiased_test = debiased_test,
  BlockMMD_test = BlockMMD_test,
  bootstrap_MMD_test = bootstrap_MMD_test
)

cit_test_functions <- list(
  RCIT_test = RCIT_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test,
  PCM_test = PCM_test,
  KCI_test = KCI_test
)

# Regression methods for CIT tests
regression_methods <- list(
  lm = list(regr.method = lm_reg_method, binary.regr.method = lm_reg_method_binary),
  ranger = list(regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary),
  xgboost = list(regr.method = xgboost_reg_method, binary.regr.method = xgboost_reg_method_binary)
)

n_values <- c(200, 400, 800, 1200, 1600, 2000)
n_sims <- 50  # Reduced number of simulations for computational cost measurement
estimators <- c("LL", "KLR")
results_list <- list()

# Only test under null hypothesis
is_null <- TRUE

for (n in n_values) {
  for (test_type in c("DRT", "CIT")) {
    test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
    
    for (test_name in names(test_functions)) {
      if (test_type == "DRT") {
        # DRT tests: test with different estimators
        for (est in estimators) {
          cat("[Measuring]", test_name, "| n:", n, "| Estimator:", est, "\n")
          
          time_results <- numeric(n_sims)
          for (sim in 1:n_sims) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # generate data
            d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
            set.seed(seed + n_sims)
            d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
            
            test_args <- list(d1$x, d2$x, d1$y, d2$y, est.method = est, seed = seed)
            
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
            estimator = est,
            avg_time_seconds = avg_time,
            std_time_seconds = std_time
          )
          
          cat("  Avg time:", avg_time, "seconds | SD:", std_time, "\n", strrep("-", 80), "\n")
        }
      } else {
        # CIT tests
        if (test_name == "RCIT_test") {
          # RCIT: no regression method variation
          cat("[Measuring]", test_name, "| n:", n, "\n")
          
          time_results <- numeric(n_sims)
          for (sim in 1:n_sims) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # generate data
            d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
            set.seed(seed + n_sims)
            d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
            
            epsilon <- 1 / sqrt(log(n))
            test_args <- list(d1$x, d2$x, d1$y, d2$y, alg1 = TRUE, epsilon = epsilon, seed = seed)
            
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
              
              # generate data
              d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
              set.seed(seed + n_sims)
              d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
              
              epsilon <- 1 / sqrt(log(n))
              test_args <- c(
                list(d1$x, d2$x, d1$y, d2$y, alg1 = TRUE, epsilon = epsilon, seed = seed),
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

results_dt <- rbindlist(results_list)

# Save the results
filename <- paste0("results/computation_costs_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")

