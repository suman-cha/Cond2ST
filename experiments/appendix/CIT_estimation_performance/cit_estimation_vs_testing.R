#' CIT Estimation Performance vs Testing Performance Analysis
#'
#' This script addresses the reviewer's comment: "How estimation affects 
#' instability and as a consequence the alternative approach"
#'
#' Compares 4 CIT methods (GCM, PCM, WGSC, RCIT) by measuring:
#' 1. Estimation quality (regression prediction errors)
#' 2. Testing performance (rejection rates, power)
#' 3. Test stability (variance of test statistics)
#' 4. Correlation between estimation errors and testing instability
#'
#' Methods:
#' - GCM: Regression-based (X_on_Z and mhat)
#' - PCM: Regression-based (ghat, mtilde, mhat)
#' - WGSC: Regression-based (full and reduced models)
#' - RCIT: Kernel-based (measures test statistic stability)
#'
#' @references Lee, Cha, Kim (2024). arXiv:2410.16636
#' @author Conditional Two-Sample Testing Research Team
#' @date 2024-11-17

rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(vimp)
  library(ranger)
  library(RCIT)
})

source("./experiments/utils.R")
source("./experiments/CIT_functions.R")

# ------------------------------------------------------------------------------
# Modified CIT Functions with Estimation Metrics
# ------------------------------------------------------------------------------

#' Modified GCM Test with Estimation Metrics
#'
#' Returns both test decision and estimation errors for regression components
#'
#' @return List with rejection, mse_X_on_Z, mse_mhat, test_statistic, p_value
gcm_test_with_metrics <- function(Y, X, Z, reg_method, binary_reg_method,
                                   reg_params = list(), seed = NULL,
                                   alpha = 0.05) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(Y)
  
  # Fit X_on_Z regression
  X_on_Z <- do.call(reg_method, c(
    list(X = Z, y = X, seed = seed),
    reg_params[["X_on_Z"]]
  ))
  X_on_Z_pred <- X_on_Z(Z)
  eps <- X - X_on_Z_pred
  mse_X_on_Z <- mean((X - X_on_Z_pred)^2)
  
  # Fit mhat regression
  mhat <- do.call(binary_reg_method, c(
    list(X = Z, y = Y, seed = seed),
    reg_params[["mhat"]]
  ))
  mhat_pred <- mhat(Z)
  xi <- Y - mhat_pred
  mse_mhat <- mean((Y - mhat_pred)^2)
  
  # Compute test statistic
  R <- eps * xi
  test_statistic <- mean(R) / stats::sd(R) * sqrt(n)
  p_value <- 2 * pnorm(-abs(test_statistic))
  rejection <- as.integer(p_value < alpha)
  
  return(list(
    rejection = rejection,
    mse_X_on_Z = mse_X_on_Z,
    mse_mhat = mse_mhat,
    mse_total = (mse_X_on_Z + mse_mhat) / 2,
    test_statistic = test_statistic,
    p_value = p_value
  ))
}

#' Modified PCM Test with Estimation Metrics
#'
#' Returns both test decision and estimation errors for multiple components
#'
#' @return List with rejection, mse_ghat, mse_mtilde, mse_mhat, test_statistic
pcm_test_with_metrics <- function(Y, X, Z, reg_method, binary_reg_method,
                                   var_min = 0.01, reg_params = list(),
                                   seed = NULL, alpha = 0.05) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  Y <- as.numeric(Y)
  n <- length(Y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  
  # Sample splitting
  n_dir <- floor(n / 2)
  n_main <- n - n_dir
  direction_indices <- sample(1:n, n_dir)
  main_indices <- setdiff(1:n, direction_indices)
  
  X_dir <- X[direction_indices, , drop = FALSE]
  Z_dir <- Z[direction_indices, , drop = FALSE]
  Y_dir <- Y[direction_indices]
  
  # Fit ghat (full model)
  if (!is.null(seed)) set.seed(seed + 1)
  ghat <- do.call(binary_reg_method, c(
    list(X = cbind(X_dir, Z_dir), y = Y_dir, seed = seed + 1),
    reg_params[["ghat"]]
  ))
  gtilde <- function(X, Z) ghat(cbind(X, Z))
  ghat_dir <- ghat(cbind(X_dir, Z_dir))
  gtilde_dir <- gtilde(X_dir, Z_dir)
  mse_ghat <- mean((Y_dir - ghat_dir)^2)
  
  # Fit mtilde (reduced model on ghat predictions)
  if (!is.null(seed)) set.seed(seed + 2)
  mtilde <- do.call(reg_method, c(
    list(X = Z_dir, y = ghat_dir, seed = seed + 2),
    reg_params[["mtilde"]]
  ))
  mtilde_dir <- mtilde(Z_dir)
  mse_mtilde <- mean((ghat_dir - mtilde_dir)^2)
  
  htilde_dir <- gtilde_dir - mtilde_dir
  rho <- mean((Y_dir - ghat_dir + gtilde_dir - mtilde_dir) * htilde_dir)
  sgn <- ifelse((rho < 0), -1, 1)
  
  # Variance estimation
  sqr_resid_dir <- (Y_dir - ghat_dir)^2
  if (!is.null(seed)) set.seed(seed + 3)
  vtilde_fit <- do.call(reg_method, c(
    list(X = cbind(X_dir, Z_dir), y = sqr_resid_dir, seed = seed + 3),
    reg_params[["vtilde"]]
  ))
  vtilde_dir <- pmax(vtilde_fit(cbind(X_dir, Z_dir)), var_min)
  a <- function(c) mean(sqr_resid_dir / (pmax(vtilde_dir, 0) + c),
                        na.rm = TRUE)
  
  if (a(0) <= 1) {
    chat <- 0
  } else {
    chat <- tryCatch({
      uniroot(function(c) a(c) - 1, c(0, 10), extendInt = "yes")$root
    }, error = function(e) 0)
  }
  
  vhat <- function(X, Z) {
    ghat_eval <- ghat(cbind(X, Z))
    return(pmax(ghat_eval * (1 - ghat_eval) + chat, var_min))
  }
  
  Z_main <- Z[main_indices, , drop = FALSE]
  X_main <- X[main_indices, , drop = FALSE]
  Y_main <- Y[main_indices]
  
  fhat_main <- sgn * (gtilde(X_main, Z_main) - mtilde(Z_main)) /
    vhat(X_main, Z_main)
  
  # Fit mhat
  if (!is.null(seed)) set.seed(seed + 4)
  mhat <- do.call(binary_reg_method, c(
    list(X = Z_main, y = Y_main, seed = seed + 4),
    reg_params[["mhat"]]
  ))
  mhat_pred <- mhat(Z_main)
  eps <- Y_main - mhat_pred
  mse_mhat <- mean((Y_main - mhat_pred)^2)
  
  # Fit mhat_fhat
  if (!is.null(seed)) set.seed(seed + 5)
  mhat_fhat <- do.call(reg_method, c(
    list(X = Z_main, y = fhat_main, seed = seed + 5),
    reg_params[["mhat_fhat"]]
  ))
  xi <- fhat_main - mhat_fhat(Z_main)
  
  R <- xi * eps
  test_statistic <- sqrt(length(R)) * mean(R) / stats::sd(R)
  p_value <- 1 - pnorm(test_statistic)
  rejection <- as.integer(p_value < alpha)
  
  return(list(
    rejection = rejection,
    mse_ghat = mse_ghat,
    mse_mtilde = mse_mtilde,
    mse_mhat = mse_mhat,
    mse_total = (mse_ghat + mse_mtilde + mse_mhat) / 3,
    test_statistic = test_statistic,
    p_value = p_value
  ))
}

#' Modified WGSC Test with Estimation Metrics
#'
#' Returns both test decision and estimation errors for full/reduced models
#'
#' @return List with rejection, mse_full, mse_reduced, test_statistic, p_value
wgsc_test_with_metrics <- function(Y, X, Z, reg_method, binary_reg_method,
                                    reg_params = list(), seed = NULL,
                                    alpha = 0.05) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(Y)
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  
  full_fitted <- numeric(n)
  reduced_fitted <- numeric(n)
  data_folds <- sample(rep(1:10, length.out = n))
  sample_splitting_folds <- vimp::make_folds(unique(data_folds), V = 2)
  
  # Cross-fitting for full and reduced models
  for (j in 1:10) {
    full_fit <- do.call(
      binary_reg_method,
      c(
        list(
          X = cbind(X, Z)[data_folds != j, ],
          y = Y[data_folds != j],
          seed = seed
        ),
        reg_params[["ghat"]]
      )
    )
    full_fitted[data_folds == j] <- full_fit(cbind(X, Z)[data_folds == j, ])
    
    reduced_fit <- do.call(
      reg_method,
      c(
        list(
          X = Z[data_folds != j, ],
          y = full_fit(cbind(X, Z)[data_folds != j, ]),
          seed = seed
        ),
        reg_params[["mtilde"]]
      )
    )
    reduced_fitted[data_folds == j] <- reduced_fit(Z[data_folds == j, ])
  }
  
  # Compute MSE for full and reduced models
  mse_full <- mean((Y - full_fitted)^2)
  mse_reduced <- mean((full_fitted - reduced_fitted)^2)
  
  # Compute test statistic using vimp
  suppressWarnings(
    est <- vimp::cv_vim(
      Y = Y, cross_fitted_f1 = full_fitted,
      cross_fitted_f2 = reduced_fitted, V = 2, type = "r_squared",
      cross_fitting_folds = data_folds,
      sample_splitting_folds = sample_splitting_folds,
      run_regression = FALSE, alpha = alpha
    )
  )
  
  p_value <- est$p_value
  rejection <- as.integer(p_value < alpha)
  test_statistic <- est$test_statistic
  
  return(list(
    rejection = rejection,
    mse_full = mse_full,
    mse_reduced = mse_reduced,
    mse_total = (mse_full + mse_reduced) / 2,
    test_statistic = test_statistic,
    p_value = p_value
  ))
}

#' Modified RCIT Test with Metrics
#'
#' Returns test decision and test statistic for stability analysis
#'
#' @return List with rejection, test_statistic, p_value
rcit_test_with_metrics <- function(Y, X, Z, seed = NULL, alpha = 0.05) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  result <- tryCatch({
    RCIT(x = X, y = Y, z = Z, approx = "lpd4", num_f = 100,
         num_f2 = 5, seed = seed)
  }, error = function(e) {
    list(Sta = NA, p = 1)
  })
  
  test_statistic <- result$Sta
  p_value <- result$p
  rejection <- as.integer(p_value < alpha)
  
  return(list(
    rejection = rejection,
    test_statistic = test_statistic,
    p_value = p_value
  ))
}

# ------------------------------------------------------------------------------
# Scenario Generators
# ------------------------------------------------------------------------------

# Scenario 1: Mean shift with t-distributed noise (from scenario1u.R)
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

# Scenario 2: Heteroscedastic conditional variance (from scenario2u.R)
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

# Scenario 3: Nonlinear transformation (from scenario3u.R)
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
# Main Simulation Function
# ------------------------------------------------------------------------------

run_cit_estimation_experiment <- function(n, scenario_name, gen_x, gen_y,
                                           p_dim, is_null, seed) {
  set.seed(seed)
  
  # Generate data for two groups
  x1 <- gen_x(n, p_dim, group = 1)
  y1 <- gen_y(x1, is_null = TRUE)
  
  set.seed(seed + 10000)
  x2 <- gen_x(n, p_dim, group = 2)
  y2 <- gen_y(x2, is_null = is_null)
  
  # Merge data for CIT format: test Y indep G | X
  X_merged <- rbind(x1, x2)
  Y_merged <- c(y1, y2)
  Z_merged <- c(rep(0, n), rep(1, n))
  
  # Subsample for computational efficiency (use 80% of data)
  sample_size <- floor(0.8 * length(Y_merged))
  sample_idx <- sample(seq_along(Y_merged), sample_size)
  X_sub <- X_merged[sample_idx, , drop = FALSE]
  Y_sub <- Y_merged[sample_idx]
  Z_sub <- Z_merged[sample_idx]
  
  # Run CIT methods with metrics
  results <- list()
  
  # GCM
  results$gcm <- tryCatch({
    gcm_test_with_metrics(
      Y = Z_sub, X = Y_sub, Z = X_sub,
      reg_method = ranger_reg_method,
      binary_reg_method = ranger_reg_method_binary,
      seed = seed, alpha = 0.05
    )
  }, error = function(e) {
    list(rejection = 0, mse_total = NA, test_statistic = NA, p_value = 1,
         mse_X_on_Z = NA, mse_mhat = NA)
  })
  
  # PCM
  results$pcm <- tryCatch({
    pcm_test_with_metrics(
      Y = Z_sub, X = Y_sub, Z = X_sub,
      reg_method = ranger_reg_method,
      binary_reg_method = ranger_reg_method_binary,
      seed = seed, alpha = 0.05
    )
  }, error = function(e) {
    list(rejection = 0, mse_total = NA, test_statistic = NA, p_value = 1,
         mse_ghat = NA, mse_mtilde = NA, mse_mhat = NA)
  })
  
  # WGSC
  results$wgsc <- tryCatch({
    wgsc_test_with_metrics(
      Y = Z_sub, X = Y_sub, Z = X_sub,
      reg_method = ranger_reg_method,
      binary_reg_method = ranger_reg_method_binary,
      seed = seed, alpha = 0.05
    )
  }, error = function(e) {
    list(rejection = 0, mse_total = NA, test_statistic = NA, p_value = 1,
         mse_full = NA, mse_reduced = NA)
  })
  
  # RCIT
  results$rcit <- tryCatch({
    rcit_test_with_metrics(
      Y = Z_sub, X = Y_sub, Z = X_sub,
      seed = seed, alpha = 0.05
    )
  }, error = function(e) {
    list(rejection = 0, test_statistic = NA, p_value = 1)
  })
  
  return(results)
}

# ------------------------------------------------------------------------------
# Experiment Configuration
# ------------------------------------------------------------------------------

scenarios <- list(
  S1 = list(
    name = "Mean Shift (t-noise)",
    gen_x = sc1_generate_x,
    gen_y = sc1_generate_y
  ),
  S2 = list(
    name = "Heteroscedastic Variance",
    gen_x = sc2_generate_x,
    gen_y = sc2_generate_y
  ),
  S3 = list(
    name = "Nonlinear Transformation",
    gen_x = sc3_generate_x,
    gen_y = sc3_generate_y
  )
)

n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
p_dim <- 10
results_list <- list()

# ------------------------------------------------------------------------------
# Main Simulation Loop
# ------------------------------------------------------------------------------

cat("\n", strrep("=", 80), "\n")
cat("CIT Estimation Performance vs Testing Performance Analysis\n")
cat(strrep("=", 80), "\n\n")

for (sc_tag in names(scenarios)) {
  sc <- scenarios[[sc_tag]]
  
  for (n in n_values) {
    for (is_null in c(TRUE, FALSE)) {
      h_label <- if (is_null) "Null" else "Alternative"
      
      cat(sprintf("\n[Scenario: %s | n: %d | Hypothesis: %s]\n",
                  sc_tag, n, h_label))
      cat(strrep("-", 80), "\n")
      
      # Run simulations with progress bar
      sim_results <- pbapply::pblapply(seq_len(n_sims), function(sim) {
        seed <- 1203 + sim
        run_cit_estimation_experiment(
          n = n,
          scenario_name = sc_tag,
          gen_x = sc$gen_x,
          gen_y = sc$gen_y,
          p_dim = p_dim,
          is_null = is_null,
          seed = seed
        )
      })
      
      # Aggregate results for each method
      methods <- c("gcm", "pcm", "wgsc", "rcit")
      
      for (method in methods) {
        # Extract metrics across simulations
        rejections <- sapply(sim_results, function(x) x[[method]]$rejection)
        test_stats <- sapply(sim_results, function(x) x[[method]]$test_statistic)
        p_values <- sapply(sim_results, function(x) x[[method]]$p_value)
        
        if (method == "gcm") {
          mse_total <- sapply(sim_results, function(x) x[[method]]$mse_total)
          mse_1 <- sapply(sim_results, function(x) x[[method]]$mse_X_on_Z)
          mse_2 <- sapply(sim_results, function(x) x[[method]]$mse_mhat)
        } else if (method == "pcm") {
          mse_total <- sapply(sim_results, function(x) x[[method]]$mse_total)
          mse_1 <- sapply(sim_results, function(x) x[[method]]$mse_ghat)
          mse_2 <- sapply(sim_results, function(x) x[[method]]$mse_mtilde)
          mse_3 <- sapply(sim_results, function(x) x[[method]]$mse_mhat)
        } else if (method == "wgsc") {
          mse_total <- sapply(sim_results, function(x) x[[method]]$mse_total)
          mse_1 <- sapply(sim_results, function(x) x[[method]]$mse_full)
          mse_2 <- sapply(sim_results, function(x) x[[method]]$mse_reduced)
        } else {  # RCIT
          mse_total <- rep(NA, n_sims)
          mse_1 <- rep(NA, n_sims)
          mse_2 <- rep(NA, n_sims)
        }
        
        # Store aggregated results
        results_list[[length(results_list) + 1]] <- data.table(
          scenario = sc_tag,
          n = n,
          h_label = h_label,
          method = toupper(method),
          rejection_rate = mean(rejections, na.rm = TRUE),
          mean_mse_total = mean(mse_total, na.rm = TRUE),
          sd_mse_total = stats::sd(mse_total, na.rm = TRUE),
          mean_test_stat = mean(test_stats, na.rm = TRUE),
          sd_test_stat = stats::sd(test_stats, na.rm = TRUE),
          cv_test_stat = stats::sd(test_stats, na.rm = TRUE) /
            abs(mean(test_stats, na.rm = TRUE)),
          mean_p_value = mean(p_values, na.rm = TRUE),
          sd_p_value = stats::sd(p_values, na.rm = TRUE)
        )
        
        # Store individual simulation results for correlation analysis
        for (sim in seq_len(n_sims)) {
          results_list[[length(results_list) + 1]] <- data.table(
            scenario = sc_tag,
            n = n,
            h_label = h_label,
            method = toupper(method),
            sim_id = sim,
            rejection = rejections[sim],
            mse_total = mse_total[sim],
            test_statistic = test_stats[sim],
            p_value = p_values[sim],
            type = "individual"
          )
        }
      }
      
      # Print summary
      cat(sprintf("  GCM: Rejection Rate = %.3f\n",
                  mean(sapply(sim_results, function(x) x$gcm$rejection))))
      cat(sprintf("  PCM: Rejection Rate = %.3f\n",
                  mean(sapply(sim_results, function(x) x$pcm$rejection))))
      cat(sprintf("  WGSC: Rejection Rate = %.3f\n",
                  mean(sapply(sim_results, function(x) x$wgsc$rejection))))
      cat(sprintf("  RCIT: Rejection Rate = %.3f\n",
                  mean(sapply(sim_results, function(x) x$rcit$rejection))))
    }
  }
}

# ------------------------------------------------------------------------------
# Save Results
# ------------------------------------------------------------------------------

results_dt <- rbindlist(results_list, fill = TRUE)

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Save main results
outfile <- file.path("results", "cit_estimation_vs_testing.csv")
fwrite(results_dt, outfile)

cat("\n", strrep("=", 80), "\n")
cat("Simulation completed successfully!\n")
cat(sprintf("Results saved to: %s\n", outfile))
cat(sprintf("Total configurations: %d scenarios × %d n values × 2 hypotheses\n",
            length(scenarios), length(n_values)))
cat(sprintf("Methods tested: %d (GCM, PCM, WGSC, RCIT)\n", 4))
cat(sprintf("Simulations per configuration: %d\n", n_sims))
cat(strrep("=", 80), "\n\n")

cat("Next step: Run visualization script\n")
cat("  Source: experiments/appendix/CIT_estimation_performance/",
    "create_visualizations.R\n")

