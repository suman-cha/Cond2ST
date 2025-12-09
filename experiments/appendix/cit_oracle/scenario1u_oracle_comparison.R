rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(data.table)
})
tag <- "scenario1u"
source("./experiments/all_tests.R")

# Data generation functions (from original scenario1u)
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

# Oracle data generation: Z ~ Bernoulli(0.5), then (X,Y) | Z
generate_oracle_data <- function(n, p, is_null_overall = TRUE) {
  # Generate Z ~ Bernoulli(0.5)
  Z <- rbinom(n, 1, 0.5)
  
  # Initialize X and Y
  X <- matrix(0, nrow = n, ncol = p)
  Y <- numeric(n)
  
  # Generate (X, Y) | Z for each sample
  for (i in 1:n) {
    group <- Z[i] + 1  # Convert 0/1 to 1/2
    X[i, ] <- generate_data(1, p, group)
    
    # For null: both groups generate Y the same way (is_null = TRUE)
    # For alternative: group 1 (Z=0) uses is_null=TRUE, group 2 (Z=1) uses is_null=FALSE
    if (is_null_overall) {
      Y[i] <- generate_y(matrix(X[i, ], nrow = 1), is_null = TRUE)
    } else {
      Y[i] <- generate_y(matrix(X[i, ], nrow = 1), is_null = (Z[i] == 0))
    }
  }
  
  return(list(X = X, Y = Y, Z = Z))
}

# CIT test functions
cit_test_functions <- list(
  RCIT_test = RCIT_test,
  PCM_test = PCM_test,
  GCM_test = GCM_test,
  WGSC_test = WGSC_test,
  KCI_test = KCI_test
)

# Regression methods for CIT tests (except RCIT)
regression_methods <- list(
  lm = list(regr.method = lm_reg_method, binary.regr.method = lm_reg_method_binary),
  ranger = list(regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary),
  xgboost = list(regr.method = xgboost_reg_method, binary.regr.method = xgboost_reg_method_binary)
)

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
    
    for (test_name in names(cit_test_functions)) {
      
      if (test_name == "RCIT_test") {
        # RCIT: no regression method variation
        
        # ===== Oracle Setting =====
        cat("[Oracle]", test_name, "| n:", n, "|", h_label, "\n")
        result_oracle <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          
          # Generate oracle data
          oracle_data <- generate_oracle_data(n, d, is_null_overall = is_null)
          
          # Test Y ⊥ Z | X directly (no subsampling needed)
          # Call RCIT with alg1 = FALSE
          test_args <- list(
            x1 = oracle_data$X[oracle_data$Z == 0, , drop = FALSE],
            x2 = oracle_data$X[oracle_data$Z == 1, , drop = FALSE],
            y1 = oracle_data$Y[oracle_data$Z == 0],
            y2 = oracle_data$Y[oracle_data$Z == 1],
            alg1 = FALSE,
            seed = seed
          )
          
          do.call(cit_test_functions[[test_name]], test_args)
        }, simplify = "array")
        
        mean_result_oracle <- mean(result_oracle)
        results_list[[length(results_list) + 1]] <- data.table(
          test_name = test_name,
          n = n,
          h_label = h_label,
          regression_method = NA_character_,
          setting = "Oracle",
          rejection_rate = mean_result_oracle,
          avg_sample_size = n
        )
        cat("  Rejection Rate:", mean_result_oracle, "\n")
        
        # ===== Subsampling Setting =====
        cat("[Subsampling]", test_name, "| n:", n, "|", h_label, "\n")
        sample_sizes <- numeric(n_sims)
        result_subsampling <- pbapply::pbsapply(1:n_sims, function(sim) {
          seed <- 1203 + sim
          set.seed(seed)
          
          # Generate data for Group 1 and 2 (standard approach)
          x1 <- generate_data(n, d, group = 1)
          y1 <- generate_y(x1, is_null = TRUE)
          set.seed(seed + n_sims)
          x2 <- generate_data(n, d, group = 2)
          y2 <- generate_y(x2, is_null = is_null)
          
          epsilon <- 1 / sqrt(log(n))
          test_args <- list(x1, x2, y1, y2, alg1 = TRUE, epsilon = epsilon, seed = seed)
          
          # Track effective sample size (would need to modify test function to return this)
          # For now, we just run the test
          do.call(cit_test_functions[[test_name]], test_args)
        }, simplify = "array")
        
        mean_result_subsampling <- mean(result_subsampling)
        results_list[[length(results_list) + 1]] <- data.table(
          test_name = test_name,
          n = n,
          h_label = h_label,
          regression_method = NA_character_,
          setting = "Subsampling",
          rejection_rate = mean_result_subsampling,
          avg_sample_size = NA_real_  # Would need to track from alg1
        )
        cat("  Rejection Rate:", mean_result_subsampling, "\n", strrep("-", 80), "\n")
        
      } else {
        # Other CIT tests: test with different regression methods
        for (reg_name in names(regression_methods)) {
          
          # ===== Oracle Setting =====
          cat("[Oracle]", test_name, "| n:", n, "| Regression:", reg_name, "|", h_label, "\n")
          result_oracle <- pbapply::pbsapply(1:n_sims, function(sim) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # Generate oracle data
            oracle_data <- generate_oracle_data(n, d, is_null_overall = is_null)
            
            # Test Y ⊥ Z | X directly (no subsampling needed)
            test_args <- c(
              list(
                x1 = oracle_data$X[oracle_data$Z == 0, , drop = FALSE],
                x2 = oracle_data$X[oracle_data$Z == 1, , drop = FALSE],
                y1 = oracle_data$Y[oracle_data$Z == 0],
                y2 = oracle_data$Y[oracle_data$Z == 1],
                alg1 = FALSE,
                seed = seed
              ),
              regression_methods[[reg_name]]
            )
            
            do.call(cit_test_functions[[test_name]], test_args)
          }, simplify = "array")
          
          mean_result_oracle <- mean(result_oracle)
          results_list[[length(results_list) + 1]] <- data.table(
            test_name = test_name,
            n = n,
            h_label = h_label,
            regression_method = reg_name,
            setting = "Oracle",
            rejection_rate = mean_result_oracle,
            avg_sample_size = n
          )
          cat("  Rejection Rate:", mean_result_oracle, "\n")
          
          # ===== Subsampling Setting =====
          cat("[Subsampling]", test_name, "| n:", n, "| Regression:", reg_name, "|", h_label, "\n")
          result_subsampling <- pbapply::pbsapply(1:n_sims, function(sim) {
            seed <- 1203 + sim
            set.seed(seed)
            
            # Generate data for Group 1 and 2 (standard approach)
            x1 <- generate_data(n, d, group = 1)
            y1 <- generate_y(x1, is_null = TRUE)
            set.seed(seed + n_sims)
            x2 <- generate_data(n, d, group = 2)
            y2 <- generate_y(x2, is_null = is_null)
            
            epsilon <- 1 / sqrt(log(n))
            test_args <- c(
              list(x1, x2, y1, y2, alg1 = TRUE, epsilon = epsilon, seed = seed),
              regression_methods[[reg_name]]
            )
            
            do.call(cit_test_functions[[test_name]], test_args)
          }, simplify = "array")
          
          mean_result_subsampling <- mean(result_subsampling)
          results_list[[length(results_list) + 1]] <- data.table(
            test_name = test_name,
            n = n,
            h_label = h_label,
            regression_method = reg_name,
            setting = "Subsampling",
            rejection_rate = mean_result_subsampling,
            avg_sample_size = NA_real_  # Would need to track from alg1
          )
          cat("  Rejection Rate:", mean_result_subsampling, "\n", strrep("-", 80), "\n")
        }
      }
    }
  }
}

# Combine all results into a single data.table
results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/oracle_comparison_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")

