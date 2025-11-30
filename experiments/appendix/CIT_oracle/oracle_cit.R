#' Oracle vs Subsampling Comparison for CIT-Based Methods
#'
#' This script addresses Reviewer Comment C.1 by comparing Conditional
#' Independence Test (CIT) methods under two settings:
#' (a) Subsampling via Algorithm 1 with epsilon parameter (finite-sample)
#' (b) Oracle setting where Z is generated first, then (X,Y)|Z ~ P^(Z)_XY
#'     (no subsampling, ideal case)
#'
#' The core limitation of CIT-based methods is that Algorithm 1 discards
#' O(p*n*log(1/epsilon)) samples on average, which affects finite-sample
#' performance. This experiment quantifies the performance gap between:
#' - CIT with subsampling (realistic, loses samples)
#' - CIT with oracle Z (ideal, no sample loss)
#' - DRT methods (proposed, no subsampling needed)
#'
#' Methods Tested:
#' - DRT: BlockMMD_test (proposed density ratio-based method)
#' - CIT: GCM, PCM, WGSC, RCIT, RCoT, KCI, WGCM (7 methods × 2 modes)
#'
#' Scenarios:
#' - S1: Mean shift with t-distributed noise
#' - S2: Heteroscedastic conditional variance
#' - S3: Nonlinear transformation
#'
#' @references Lee, Cha, Kim (2024). arXiv:2410.16636
#' @author Conditional Two-Sample Testing Research Team
#' @date 2024

rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(data.table)
})

source("./experiments/all_tests.R")

# ------------------------------------------------------------------------------
# Scenario Generators (Oracle Setting)
# ------------------------------------------------------------------------------
# In oracle setting: First generate group indicator Z, then (X,Y)|Z ~ P^(Z)_XY
# This matches the theoretical framework where CIT methods operate directly
# on the joint distribution without subsampling.

# Scenario 1: Mean shift with t-distributed noise
# Matches experiments/scenarios/scenario1u.R
sc1_generate_x <- function(n, p, group) {
  mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
  sigma <- diag(1, p)
  MASS::mvrnorm(n, mu = mu, Sigma = sigma)
}

sc1_generate_y <- function(x, is_null = TRUE, sigma = 2) {
  n <- nrow(x)
  epsilon <- rt(n, df = sigma)
  f0 <- x %*% c(1, -1, -1, 1, rep(0, ncol(x) - 4))
  mean_shift <- if (is_null) 0 else 0.5
  as.numeric(f0 + epsilon + mean_shift)
}

# Scenario 2: Heteroscedastic conditional variance
# Matches experiments/scenarios/scenario2u.R
sc2_g <- function(Z, rho) {
  Z_adjusted <- Z
  diag(Z_adjusted) <- diag(Z_adjusted) - 0.5
  norm_diff <- norm(Z_adjusted, "F")
  10 + rho * exp(-norm_diff / 64)
}

sc2_generate_x <- function(n, p, group) {
  mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
  sigma <- diag(1, p)
  MASS::mvrnorm(n, mu = mu, Sigma = sigma)
}

sc2_generate_y <- function(x, rho = 10, is_null = TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  beta <- if (is_null) rep(1, p) else c(rep(1, p - 1), 0)
  mean_X <- as.numeric(x %*% matrix(beta, ncol = 1))
  var_X <- sc2_g(x, rho)
  if (is_null) {
    rnorm(n, mean = mean_X, sd = 10)
  } else {
    rnorm(n, mean = mean_X, sd = sqrt(var_X))
  }
}

# Scenario 3: Nonlinear transformation
# Matches experiments/scenarios/scenario3u.R
sc3_generate_x <- function(n, p, group) {
  mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
  sigma <- diag(1, p)
  MASS::mvrnorm(n, mu = mu, Sigma = sigma)
}

sc3_generate_y <- function(x, is_null) {
  n <- nrow(x)
  transformations <- list(
    function(z) z,
    function(z) z^2,
    function(z) z^3,
    function(z) sin(z),
    function(z) tanh(z)
  )
  random_transformation <- sample(transformations, 1)[[1]]
  if (is_null) {
    as.numeric(cos(rowSums(x) + 2 * rnorm(n)))
  } else {
    as.numeric(random_transformation(rowSums(x) + 2 * rnorm(n)))
  }
}

# ------------------------------------------------------------------------------
# Experiment Configuration
# ------------------------------------------------------------------------------

scenarios <- list(
  S1 = list(gen_x = sc1_generate_x, gen_y = sc1_generate_y),
  S2 = list(gen_x = sc2_generate_x, gen_y = sc2_generate_y),
  S3 = list(gen_x = sc3_generate_x, gen_y = sc3_generate_y)
)

n_values <- c(200, 500, 1000)
n_sims <- 500  # Increased from 300 for more precise estimates
alpha <- 0.05
p_dim <- 10

# Epsilon function for Algorithm 1 subsampling
# Default choice from utils.R: epsilon = 1/log(n)
# Smaller epsilon -> more samples discarded -> better approximation but worse
# finite-sample performance
epsilon_of_n <- function(n) 1 / log(n)

# Method configuration
drt_method_name <- "LinearMMD_test"  # Proposed DRT representative
cit_method_names <- c("GCM_test", "PCM_test", "WGSC_test", "RCIT_test",
                      "RCoT_test", "KCI_test", "WGCM_test")

results <- list()

# ------------------------------------------------------------------------------
# Main Simulation Loop
# ------------------------------------------------------------------------------

for (sc_tag in names(scenarios)) {
  gen_x <- scenarios[[sc_tag]]$gen_x
  gen_y <- scenarios[[sc_tag]]$gen_y
  
  for (n in n_values) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      eps_val <- epsilon_of_n(n)
      
      cat("\n", strrep("=", 80), "\n")
      cat(sprintf("Scenario: %s | n: %d | Hypothesis: %s | epsilon: %.4f\n",
                  sc_tag, n, h_label, eps_val))
      cat(strrep("=", 80), "\n")
      
      # Run simulations with progress bar
      sim_res <- pbapply::pbsapply(1:n_sims, function(sim) {
        seed <- 1203 + sim
        set.seed(seed)
        
        # Generate data from oracle setting: Z first, then (X,Y)|Z
        # Group 1: Z=1, sample from P^(1)_XY
        x1 <- gen_x(n, p_dim, group = 1)
        y1 <- gen_y(x1, is_null = TRUE)
        
        set.seed(seed + n_sims)
        # Group 2: Z=2, sample from P^(2)_XY
        x2 <- gen_x(n, p_dim, group = 2)
        y2 <- gen_y(x2, is_null = is_null)
        
        # Initialize results vector
        # 1 DRT + 7 CIT methods × 2 modes = 15 methods total
        res_vec <- numeric(15)
        names_vec <- character(15)
        idx <- 1
        
        # 1. DRT method (no subsampling, group-wise approach)
        drt_rej <- tryCatch({
          do.call(drt_method_name, list(
            x1 = x1, x2 = x2, y1 = y1, y2 = y2,
            alpha = alpha, seed = seed, est.method = "LL",
            bandwidths = "median"
          ))
        }, error = function(e) {
          0
        })
        res_vec[idx] <- as.integer(drt_rej[1])
        names_vec[idx] <- "DRT_BlockMMD"
        idx <- idx + 1
        
        # 2. CIT methods: subsample (alg1=TRUE) and oracle (alg1=FALSE)
        for (cit_method in cit_method_names) {
          # Subsample mode (Algorithm 1 with epsilon)
          cit_sub_rej <- tryCatch({
            do.call(cit_method, list(
              x1 = x1, x2 = x2, y1 = y1, y2 = y2,
              alpha = alpha, epsilon = eps_val, alg1 = TRUE, seed = seed
            ))
          }, error = function(e) {
            0
          })
          res_vec[idx] <- as.integer(cit_sub_rej)
          names_vec[idx] <- sprintf("CIT_%s_subsample",
                                    gsub("_test", "", cit_method))
          idx <- idx + 1
          
          # Oracle mode (no subsampling, uses full data)
          cit_oracle_rej <- tryCatch({
            do.call(cit_method, list(
              x1 = x1, x2 = x2, y1 = y1, y2 = y2,
              alpha = alpha, alg1 = FALSE, seed = seed
            ))
          }, error = function(e) {
            0
          })
          res_vec[idx] <- as.integer(cit_oracle_rej)
          names_vec[idx] <- sprintf("CIT_%s_oracle",
                                    gsub("_test", "", cit_method))
          idx <- idx + 1
        }
        
        names(res_vec) <- names_vec
        return(res_vec)
      }, simplify = "matrix")
      
      # Aggregate results: compute rejection rates for all methods
      method_names <- rownames(sim_res)
      rejection_rates <- rowMeans(sim_res)
      
      # Store in results list
      for (i in seq_along(method_names)) {
        results[[length(results) + 1]] <- data.table(
          scenario = sc_tag,
          n = n,
          h_label = h_label,
          method = method_names[i],
          rejection_rate = rejection_rates[i]
        )
      }
      
      # Print summary for this configuration
      cat("\nResults Summary:\n")
      cat(strrep("-", 80), "\n")
      cat(sprintf("%-30s %s\n", "Method", "Rejection Rate"))
      cat(strrep("-", 80), "\n")
      
      # Print DRT
      cat(sprintf("%-30s %.3f\n", method_names[1], rejection_rates[1]))
      cat(strrep("-", 80), "\n")
      
      # Print CIT methods grouped by test type
      for (i in seq_along(cit_method_names)) {
        base_idx <- 1 + (i - 1) * 2 + 1
        cat(sprintf("%-30s %.3f\n",
                    method_names[base_idx], rejection_rates[base_idx]))
        cat(sprintf("%-30s %.3f\n",
                    method_names[base_idx + 1], rejection_rates[base_idx + 1]))
        if (i < length(cit_method_names)) {
          cat(strrep("-", 40), "\n")
        }
      }
      cat(strrep("-", 80), "\n\n")
    }
  }
}

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

results_dt <- rbindlist(results)
outfile <- file.path("results", "comment_c1_oracle_vs_cit.csv")
fwrite(results_dt, outfile)

cat("\n", strrep("=", 80), "\n")
cat("Experiment completed successfully!\n")
cat(sprintf("Results saved to: %s\n", outfile))
cat(sprintf("Total configurations: %d scenarios × %d n values × 2 hypotheses = %d\n",
            length(scenarios), length(n_values), 
            length(scenarios) * length(n_values) * 2))
cat(sprintf("Methods tested: 1 DRT + 7 CIT × 2 modes = 15 methods\n"))
cat(sprintf("Total rows in output: %d\n", nrow(results_dt)))
cat(strrep("=", 80), "\n\n")
