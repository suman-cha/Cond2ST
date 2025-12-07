# Install the required packages
required_pkgs <- c("MASS",
                   "splines",
                   "progress",
                   "prettyunits",
                   "pbapply",
                   "glmnet",
                   "SuperLearner",
                   "plyr",
                   "kernlab",
                   "CVST",
                   "ANN2",
                   "doRNG",
                   "doSNOW",
                   "foreach",
                   "stabs",
                   "devtools",
                   "GeneralisedCovarianceMeasure")


install_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)){
    install.packages(new_pkgs, repos="http://cran.us.r-project.org")
  }
}

install_pkgs(required_pkgs)
source("./experiments/utils.R")
source('./experiments/CP_FunFiles.R')
source("./experiments/CIT_functions.R")

# Load libraries
for (pkg in required_pkgs){
  suppressPackageStartupMessages({
    library(pkg, character.only = TRUE)
    library(RCIT)
  })
}

# install_github("ericstrobl/RCIT")
# library(RCIT)

# Density Ratio-based Test (DRT)
# Debiased conditional two-sample test (Chen and Lei 2024)
debiased_test <- function(x1, x2, y1, y2, alpha=0.05, est.method="LL", S=2, seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0
  
  n1 <- length(y1)
  n2 <- length(y2)
  d <- if (is.null(dim(x1))) 1 else dim(x1)[2]
  
  data0 <- as.data.frame(x1)
  data0$y <- y1
  data0$class <- rep(0, nrow(data0))
  
  data1 <- as.data.frame(x2)
  data1$y <- y2
  data1$class <- rep(1, nrow(data1))
  
  # Ensure consistent column names between data0 and data1
  colnames(data1) <- colnames(data0)
  
  # Sample split
  split0_ind <- sample(seq_len(nrow(data0)), size = floor(0.5 * nrow(data0)))
  split1_ind <- sample(seq_len(nrow(data1)), size = floor(0.5 * nrow(data1)))
  split0 <- data0[split0_ind, ]
  split1 <- data1[split1_ind, ]
  colnames(split1) <- colnames(split0)
  split <- rbind(split0, split1)
  
  test0 <- data0[-split0_ind, ]
  test1 <- data1[-split1_ind, ]
  
  # Construct functions for marginal and joint ratios
  marg_ratio <- estimate_marginal_ratio(split, nrow(split0), nrow(split1), est.method)
  joint_ratio <- estimate_joint_ratio(split, est.method)
  
  a_calc <- function(r0, r1) {
    if (r0 < r1) {
      return(1)
    } else {
      return(0)
    }
  }
  
  a_outer <- function(point0, point1) {
    cr0 <- joint_ratio(point0) * marg_ratio(point0)
    cr1 <- joint_ratio(point1) * marg_ratio(point1)
    return(outer(cr0, cr1, Vectorize(a_calc)))
  }
  
  # Form cross-fit folds
  ind0 <- sample(seq(nrow(test0)))
  ind1 <- sample(seq(nrow(test1)))
  fold0 <- split(ind0, ceiling(seq_along(ind0) / (nrow(test0) / S)))
  fold1 <- split(ind1, ceiling(seq_along(ind1) / (nrow(test1) / S)))
  
  store_values <- matrix(nrow = nrow(test0), ncol = nrow(test1))
  
  for (j in 1:S) {
    for (k in 1:S) {
      # Split data into estimate and nuisance
      est_data0 <- test0[fold0[[j]], ]
      est_data1 <- test1[fold1[[k]], ]
      est_data <- rbind(est_data0, est_data1)
      nuisance0 <- test0[-fold0[[j]], ]
      nuisance1 <- test1[-fold1[[k]], ]
      nuisance <- rbind(nuisance0, nuisance1)
      
      # Estimate nuisance parameters
      alpha_data <- as.data.frame(nuisance0[, 1:d])
      as <- a_outer(nuisance0, nuisance1)
      alpha_data$successes <- rowSums(as)
      alpha_data$failures <- nrow(nuisance1) - rowSums(as)

      # Use linear model for alpha model
      colnames(alpha_data) <- c(paste0("V", 1:d), "successes", "failures")

      alpha_model <- glm(cbind(successes, failures) ~ ., data = alpha_data, family = binomial())
      gamma <- estimate_marginal_ratio(nuisance, nrow(nuisance0), nrow(nuisance1), est.method)
      
      # Estimate
      as <- a_outer(est_data0, est_data1)
      gamma0 <- gamma(est_data0)
      colnames(est_data0) <- c(paste0("V", 1:d))
      colnames(est_data1) <- c(paste0("V", 1:d))
      
      alpha0 <- predict(alpha_model, est_data0, type = "response")
      alpha1 <- predict(alpha_model, est_data1, type = "response")
      for (u in 1:length(fold0[[j]])) {
        for (v in 1:length(fold1[[k]])) {
          store_values[fold0[[j]][u], fold1[[k]][v]] <- gamma0[u] * as[u, v] + alpha1[v] - gamma0[u] * alpha0[u]
        }
      }
    }
  }
  
  value <- sum(store_values) / (nrow(test0) * nrow(test1))
  
  # Variance estimate
  variance0 <- rep(0, nrow(test0))
  variance1 <- rep(0, nrow(test1))
  
  for (u in 1:nrow(test0)) {
    variance0[u] <- mean(store_values[u, ]) - 0.5
  }
  
  for (u in 1:nrow(test1)) {
    variance1[u] <- mean(store_values[, u]) - 0.5
  }
  
  variance0 <- sum(variance0^2) / (nrow(test0) - 1)
  variance1 <- sum(variance1^2) / (nrow(test1) - 1)
  variance <- (2 * variance0 + 2 * variance1) / (nrow(test0) + nrow(test1))
  
  t <- (0.5 - value) / sqrt(variance)
  rejection <- ifelse(t >= qnorm(1 - alpha), 1, 0)
  
  return(rejection)
}

# Hu and Lei (2024)
CP_test <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, est.method="LL", seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0
  
  n1 <- length(y1); n2 <- length(y2)
  n12 <- ceiling(n1 * prop); n22 <- ceiling(n2 * prop)
  n11 <- n1 - n12; n21 <- n2 - n22
  x11 <- x1[1:n11, , drop=F]; x12 <- x1[-(1:n11),,drop=F]
  y11 <- y1[1:n11]; y12 <- y1[-(1:n11)]
  x21 <- x2[1:n21, , drop=F]; x22 <- x2[-(1:n21), ,drop=F]
  y21 <- y2[1:n21]; y22 <- y2[-(1:n21)]
  ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
  temp <- getFinalStat(ratios$g12.est, ratios$g22.est, ratios$v12.est, ratios$v22.est)
  pvalue <- pnorm(median(temp$z.hm))
  if (pvalue < alpha){
    rejection <- 1
  }
  return(rejection)
}

LinearMMD_test <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, bandwidths=c(1), 
                           est.method="LL", seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    rejections <- numeric(length(bandwidths))
    
    # Assume equal sample size
    stopifnot(length(y1) == length(y2))
    total_sample_size <- length(y1)
    n <- ceiling(total_sample_size * prop)
    x11 <- x1[1:n, , drop=F]; x12 <- x1[-(1:n),,drop=F]
    y11 <- y1[1:n]; y12 <- y1[-(1:n)]
    x21 <- x2[1:n, , drop=F]; x22 <- x2[-(1:n),,drop=F]
    y21 <- y2[1:n]; y22 <- y2[-(1:n)]
    
    # Estimate density ratio once
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    r_X2 <- 1/ratios$g22.est
    
    # Calculate for all bandwidths
    for(i in seq_along(bandwidths)) {
        bw <- bandwidths[i]
        if(bw == "median") {
            h_x <- median.bandwidth(x12, x22)
            h_y <- median.bandwidth(matrix(y12), matrix(y22))
        } else {
            h_x <- as.numeric(bw)
            h_y <- as.numeric(bw)
        }
        
        MMDl_hat2 <- MMDl(x12, x22, y12, y22, h_x=h_x, h_y=h_y, r_X2, seed)
        pvalue <- 1 - pnorm(MMDl_hat2)
        rejections[i] <- as.integer(pvalue < alpha)
    }
    
    return(rejections)
}

CV_LinearMMD_test <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, bandwidths=c(1), K=2, 
                              est.method="LL", seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    rejections <- numeric(length(bandwidths))
    
    # Assume equal sample size 
    stopifnot(length(y1) == length(y2))
    total_sample_size <- length(y1)
    fold_size <- total_sample_size %/% K
    n <- floor(total_sample_size * prop)
    S_bar_values <- matrix(0, nrow = K, ncol = length(bandwidths))
    
    for(j in 1:K) {
        start_idx <- (j - 1) * fold_size + 1
        end_idx <- j * fold_size
        
        test_idx <- start_idx:end_idx
        est_idx <- setdiff(1:total_sample_size, test_idx)
        
        x11 <- x1[est_idx, , drop = FALSE]; y11 <- y1[est_idx]
        x12 <- x1[test_idx, , drop = FALSE]; y12 <- y1[test_idx]
        x21 <- x2[est_idx, , drop = FALSE]; y21 <- y2[est_idx]
        x22 <- x2[test_idx, , drop = FALSE]; y22 <- y2[test_idx]
        
        # Estimate density ratio once per fold
        ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
        r_X2 <- 1/ratios$g22.est
        
        # Calculate for all bandwidths
        for(bw_idx in seq_along(bandwidths)) {
            bw <- bandwidths[bw_idx]
            if(bw == "median") {
                h_x <- median.bandwidth(x12, x22)
                h_y <- median.bandwidth(matrix(y12), matrix(y22))
            } else {
                h_x <- as.numeric(bw)
                h_y <- as.numeric(bw)
            }
            
            S_bar_values[j, bw_idx] <- MMDl(x12, x22, y12, y22, h_x=h_x, h_y=h_y, r_X2, seed) * 
                sqrt(n) / sqrt(floor(length(y12)/2))
        }
    }
    
    # Compute final results for each bandwidth
    for(bw_idx in seq_along(bandwidths)) {
        MMD_hat_cv <- sum(S_bar_values[, bw_idx]) / K
        pvalue <- 1 - pnorm(MMD_hat_cv)
        rejections[bw_idx] <- as.integer(pvalue < alpha)
    }
    
    return(rejections)
}

BlockMMD_test <- function(x1, x2, y1, y2, B_size = NULL, prop=0.5, alpha=0.05, 
                           bandwidths=c(1), est.method='LL', seed=NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Assume equal sample size    
    stopifnot(length(y1) == length(y2))
    total_n <- length(y1)
    
    n <- ceiling(total_n * prop)
    x11 <- x1[1:n,,drop=FALSE]; x12 <- x1[-(1:n),,drop=FALSE]
    x21 <- x2[1:n,,drop=FALSE]; x22 <- x2[-(1:n),,drop=FALSE]
    y11 <- y1[1:n]; y12 <- y1[-(1:n)]
    y21 <- y2[1:n]; y22 <- y2[-(1:n)]
    
    # Estimate density ratio
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    r_X2 <- 1/ratios$g22.est
    
    n_test <- length(y12)
    if (is.null(B_size)) B_size <- floor(n_test^(0.5))
    
    num_blocks <- floor(n_test / B_size)
    if (num_blocks < 2) {
        warning("Test sample is too small for the given B_size, resulting in less than 2 blocks.")
        return(rep(NA, length(bandwidths)))
    }
    
    rejections <- numeric(length(bandwidths))
    for (i in seq_along(bandwidths)) {
        bw <- bandwidths[i]
        
        if (bw == "median") {
            h_x <- median.bandwidth(x12, x22)
            h_y <- median.bandwidth(matrix(y12), matrix(y22))
        } else {
            h_x <- as.numeric(bw)
            h_y <- as.numeric(bw)
        }
        
        # Compute the test statistic
        test_stat <- MMDb(x12, x22, y12, y22, B_size,
                          h_x=h_x, h_y=h_y, r_X=r_X2, seed=seed)
        
        pval <- 1 - pnorm(test_stat)
        rejections[i] <- as.integer(pval < alpha)
    }
    return(rejections)
}

# CV_BlockMMD_test <- function(x1, x2, y1, y2, B_size = NULL, K = 2, alpha = 0.05,
#                              bandwidths = c(1), est.method = "LL", seed = NULL) {
#     if (!is.null(seed)) set.seed(seed)
#     
#     stopifnot(nrow(x1) == length(y1), nrow(x2) == length(y2))
#     total_n <- length(y1)
#     folds <- sample(rep(1:K, length.out = total_n))
#     
#     rejections <- numeric(length(bandwidths))
#     fold_size <- floor(total_n / K)
#     if (is.null(B_size)) B_size <- floor(fold_size^(0.5))
#     
#     for (i in seq_along(bandwidths)) {
#         bw <- bandwidths[i]
#         fold_stats <- numeric(K)
#         
#         for (j in 1:K) {
#             test_idx <- which(folds == j)
#             est_idx <- which(folds != j)
#             
#             x1_est <- x1[est_idx, , drop = FALSE]; y1_est <- y1[est_idx]
#             x1_test <- x1[test_idx, , drop = FALSE]; y1_test <- y1[test_idx]
#             x2_est <- x2[est_idx, , drop = FALSE]; y2_est <- y2[est_idx]
#             x2_test <- x2[test_idx, , drop = FALSE]; y2_test <- y2[test_idx]
#             
#             n_test <- length(y1_test)
#             S <- floor(n_test / B_size)
#             if (S < 2) {
#                 warning(paste("Fold", j, "has insufficient test data for >=2 blocks. Skipping."))
#                 fold_stats[j] <- NA
#                 next
#             }
#             
#             ratios <- estimate_r(x1_est, x1_test, x2_est, x2_test,
#                                  y1_est, y1_test, y2_est, y2_test,
#                                  est.method = est.method, seed = seed + j)
#             r_X <- 1 / ratios$g22.est
# 
#             if (bw == "median") {
#                 h_x <- median.bandwidth(x1_test, x2_test)
#                 h_y <- median.bandwidth(matrix(y1_test), matrix(y2_test))
#             } else {
#                 h_x <- as.numeric(bw)
#                 h_y <- as.numeric(bw)
#             }
#             
#             stat_j <- MMDb(x12 = x1_test, x22 = x2_test,
#                            y12 = y1_test, y22 = y2_test,
#                            B_size = B_size, h_x = h_x, h_y = h_y,
#                            r_X = r_X, seed = seed + j)
#             fold_stats[j] <- stat_j
#         }
#         
#         final_stat <- mean(fold_stats, na.rm = TRUE)
#         if (is.nan(final_stat)) final_stat <- 0
#         pval <- 1 - pnorm(final_stat)
#         rejections[i] <- as.integer(pval < alpha)
#     }
#     
#     return(rejections)
# }

CV_BlockMMD_test <- function(x1, x2, y1, y2, B_size = NULL, K = 2, alpha = 0.05,
                             bandwidths = c(1), est.method = "LL", seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    stopifnot(nrow(x1) == length(y1), nrow(x2) == length(y2))
    total_n <- length(y1)
    folds <- sample(rep(1:K, length.out = total_n))
    
    rejections <- numeric(length(bandwidths))
    fold_size <- floor(total_n / K)
    if (is.null(B_size)) B_size <- floor(fold_size^(0.5))
    
    for (i in seq_along(bandwidths)) {
        bw <- bandwidths[i]
        all_block_values <- list()  # Store all block MMD values across folds
        
        for (j in 1:K) {
            test_idx <- which(folds == j)
            est_idx <- which(folds != j)
            
            x1_est <- x1[est_idx, , drop = FALSE]; y1_est <- y1[est_idx]
            x1_test <- x1[test_idx, , drop = FALSE]; y1_test <- y1[test_idx]
            x2_est <- x2[est_idx, , drop = FALSE]; y2_est <- y2[est_idx]
            x2_test <- x2[test_idx, , drop = FALSE]; y2_test <- y2[test_idx]
            
            n_test <- length(y1_test)
            S_fold <- floor(n_test / B_size)
            if (S_fold < 1) {
                warning(paste("Fold", j, "has insufficient test data for blocks. Skipping."))
                next
            }
            
            ratios <- estimate_r(x1_est, x1_test, x2_est, x2_test,
                                 y1_est, y1_test, y2_est, y2_test,
                                 est.method = est.method, seed = seed + j)
            r_X <- 1 / ratios$g22.est
            
            if (bw == "median") {
                h_x <- median.bandwidth(x1_test, x2_test)
                h_y <- median.bandwidth(matrix(y1_test), matrix(y2_test))
            } else {
                h_x <- as.numeric(bw)
                h_y <- as.numeric(bw)
            }
            
            # Compute individual block MMD values for this fold
            fold_block_values <- numeric(S_fold)
            for (b in 1:S_fold) {
                idx <- ((b-1)*B_size+1):(b*B_size)
                if (max(idx) > n_test) {
                    idx <- idx[idx <= n_test]  # Handle last block if incomplete
                    if (length(idx) < 2) next  # Skip if too few samples
                }
                
                xb1 <- x1_test[idx,,drop=FALSE]
                xb2 <- x2_test[idx,,drop=FALSE]
                yb1 <- y1_test[idx]
                yb2 <- y2_test[idx]
                rb <- r_X[idx]
                
                sum_k <- 0
                B_actual <- length(idx)
                for (ii in 1:(B_actual-1)) {
                    for (jj in (ii+1):B_actual) {
                        k_zz <- gaussian.kernel(xb1[ii, ], xb1[jj, ], h_x) * 
                            gaussian.kernel(yb1[ii], yb1[jj], h_y)
                        k_ww <- gaussian.kernel(xb2[ii, ], xb2[jj, ], h_x) * 
                            gaussian.kernel(yb2[ii], yb2[jj], h_y)
                        k_wz <- gaussian.kernel(xb2[ii, ], xb1[jj, ], h_x) * 
                            gaussian.kernel(yb2[ii], yb1[jj], h_y)
                        k_zw <- gaussian.kernel(xb1[ii, ], xb2[jj, ], h_x) * 
                            gaussian.kernel(yb1[ii], yb2[jj], h_y)
                        
                        sum_k <- sum_k + k_zz + rb[ii]*rb[jj]*k_ww - rb[ii]*k_wz - rb[jj]*k_zw
                    }
                }
                n_pairs <- B_actual * (B_actual - 1) / 2
                if (n_pairs > 0) {
                    fold_block_values[b] <- sum_k / n_pairs
                }
            }
            
            all_block_values[[j]] <- fold_block_values[fold_block_values != 0]
        }
        
        # Combine all block MMD values from all folds
        combined_blocks <- unlist(all_block_values)
        if (length(combined_blocks) < 2) {
            rejections[i] <- 0
            next
        }
        
        S_total <- length(combined_blocks)
        S_bar <- mean(combined_blocks)
        sigma_hat <- sqrt(var(combined_blocks))
        final_stat <- ifelse(sigma_hat > 0, sqrt(S_total) * S_bar / sigma_hat, 0)
        
        pval <- 1 - pnorm(final_stat)
        rejections[i] <- as.integer(pval < alpha)
    }
    
    return(rejections)
}

CLF_test <- function(x1, x2, y1, y2, split.prop=0.5, est.prop=0.8, alpha=0.05, est.method="LL", seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0

  # Training samples
  n1 <- length(y1); n2 <- length(y2)
  stopifnot(length(y1) == length(y2)) # Assume equal sample size
  n12 <- ceiling(n1 * split.prop); n22 <- ceiling(n2 * split.prop)
  n11 <- n1 - n12; n21 <- n2 - n22
  x11 <- x1[1:n11, , drop=F]; x12 <- x1[-(1:n11),,drop=F]
  y11 <- y1[1:n11]; y12 <- y1[-(1:n11)]
  y21 <- y2[1:n21]; y22 <- y2[-(1:n21)]
  
  # Estimation / Testing samples
  n1_rat <- ceiling(n12 * est.prop) ; n2_rat <- ceiling(n22 * est.prop)
  n1_te <- n12 - n1_rat ; n2_te <- n22 - n2_rat
  
  x21 <- x2[1:n21, , drop=F]; x22 <- x2[-(1:n21), ,drop=F]
  x12_rat <- x12[1:n1_rat, , drop=F] ; x12_te <- x12[-(1:n1_rat), , drop=F]
  x22_rat <- x22[1:n2_rat, , drop=F] ; x22_te <- x22[-(1:n2_rat), , drop=F]
  
  y12_rat <- y12[1:n1_rat] ; y12_te <- y12[-(1:n1_rat)]
  y22_rat <- y22[1:n2_rat] ; y22_te <- y22[-(1:n2_rat)]
  
  te_ratios <- estimate_r(x11, x12_te, x21, x22_te, y11, y12_te, y21, y22_te, est.method, seed)
  
  g12.est <- 1/te_ratios$g12.est
  g22.est <- 1/te_ratios$g22.est
  
  J12.est <- te_ratios$v12.est * g12.est
  J22.est <- te_ratios$v22.est * g22.est
  
  pred_V1 <- ifelse((g12.est/(J12.est+g12.est)) > 0.5, 1, 0)
  pred_V2 <- ifelse((g22.est/(J22.est+g22.est)) > 0.5, 1, 0)
  
  Acc_ratios <- estimate_r(x12_rat, x12_te, x22_rat, x22_te, y12_rat, y12_te, y22_rat, y22_te, est.method, seed)
  ratio <- 1/Acc_ratios$g22.est
  
  A1_hat <- as.integer(pred_V1 == 0)
  A2_hat <- ratio * as.integer(pred_V2 == 1)
  
  bar_A1 <- mean(A1_hat)
  bar_A2 <- mean(A2_hat)
  sigma1_squared <- var(A1_hat)
  sigma2_squared <- var(A2_hat)
  sig <- sigma1_squared + sigma2_squared
  
  if (sig > 0) {
    Acc_hat <- sqrt(n1_te) * (bar_A1 + bar_A2 - 1) / sqrt(sig)
  } else {
    Acc_hat <- 0.0
  }
  p_value <- 1 - pnorm(Acc_hat)
  
  if (p_value < alpha){
    rejection <- 1
  }
  return(rejection)
}
CV_CLF_test <- function(x1, x2, y1, y2, est.prop=0.8, alpha=0.05, K=2, est.method="LL", seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0
  
  # Assume equal sample size 
  stopifnot(length(y1) == length(y2))
  total_sample_size <- length(y1)
  fold_size <- total_sample_size %/% K
  Acc_values <- numeric(K)
  
  m <- fold_size
  
  for (j in 1:K) {
    start_idx <- (j - 1) * fold_size + 1
    end_idx <- j * fold_size
    
    test_idx <- start_idx:end_idx
    est_idx <- setdiff(1:total_sample_size, test_idx)
    n12 <- length(test_idx)
    n22 <- length(test_idx)
    
    x11 <- x1[est_idx, , drop = FALSE]; y11 <- y1[est_idx]
    x12 <- x1[test_idx, , drop = FALSE]; y12 <- y1[test_idx]
    x21 <- x2[est_idx, , drop = FALSE]; y21 <- y2[est_idx]
    x22 <- x2[test_idx, , drop = FALSE]; y22 <- y2[test_idx]
    
    n1_rat <- ceiling(n12 * est.prop) ; n2_rat <- ceiling(n22 * est.prop)
    n1_te <- n12 - n1_rat ; n2_te <- n22 - n2_rat
    
    x12_rat <- x12[1:n1_rat, , drop=F] ; x12_te <- x12[-(1:n1_rat), , drop=F]
    x22_rat <- x22[1:n2_rat, , drop=F] ; x22_te <- x22[-(1:n2_rat), , drop=F]
    
    y12_rat <- y12[1:n1_rat] ; y12_te <- y12[-(1:n1_rat)]
    y22_rat <- y22[1:n2_rat] ; y22_te <- y22[-(1:n2_rat)]
    
    # For compute test statistic
    te_ratios <- estimate_r(x11, x12_te, x21, x22_te, y11, y12_te, y21, y22_te, est.method, seed)
    
    g12.est <- 1/te_ratios$g12.est
    g22.est <- 1/te_ratios$g22.est
    
    J12.est <- te_ratios$v12.est * g12.est
    J22.est <- te_ratios$v22.est * g22.est
    
    pred_V1 <- ifelse((g12.est/(J12.est+g12.est)) > 0.5, 1, 0)
    pred_V2 <- ifelse((g22.est/(J22.est+g22.est)) > 0.5, 1, 0)
    
    Acc_ratios <- estimate_r(x12_rat, x12_te, x22_rat, x22_te, y12_rat, y12_te, y22_rat, y22_te, est.method, seed)
    ratio <- 1/Acc_ratios$g22.est
    
    A1_hat <- as.integer(pred_V1 == 0)
    A2_hat <- ratio * as.integer(pred_V2 == 1)
    
    bar_A1 <- mean(A1_hat)
    bar_A2 <- mean(A2_hat)
    sigma1_squared <- var(A1_hat)
    sigma2_squared <- var(A2_hat)
    
    sig <- sigma1_squared + sigma2_squared
    
    if (sig > 0) {
      Acc_hat <- sqrt(n1_te) * (bar_A1 + bar_A2 - 1) / sqrt(sig)
    } else {
      Acc_hat <- 0.0
    }
    Acc_values[j] <- Acc_hat
  }
  
  Acc_values_cv <- sum(Acc_values) / sqrt(K)
  p_value <- 1 - pnorm(Acc_values_cv)
  rejection <- if (p_value < alpha) 1 else 0
  
  return(rejection)
}
CLF_test2 <- function(x1, x2, y1, y2, prop = 0.5, alpha = 0.05, est.method = "LL", seed = NULL, K = 5) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n1 <- length(y1); n2 <- length(y2)
  n <- min(n1, n2)
  m <- floor(n/2)
  
  # Ensure x1 and x2 are data frames
  x1 <- as.data.frame(x1[1:n,])
  x2 <- as.data.frame(x2[1:n,])
  
  # Combine and shuffle the data
  combined_data <- rbind(
    cbind(x1, Y = y1[1:n], Group = 0),
    cbind(x2, Y = y2[1:n], Group = 1)
  )
  combined_data <- as.data.frame(combined_data)
  shuffled_indices <- sample(1:(2*n))
  combined_data <- combined_data[shuffled_indices,]
  
  # Convert Group to factor
  combined_data$Group <- as.factor(combined_data$Group)
  
  # Split into D_a and D_b
  D_a <- combined_data[1:n,]
  D_b <- combined_data[(n+1):(2*n),]
  
  # Train classifier on D_b
  classifier <- glm(Group ~ ., data = D_b[, -which(names(D_b) == "Y")], family = binomial())
  
  # K-fold cross-validation on D_a
  folds <- cut(seq(1, nrow(D_a)), breaks = K, labels = FALSE)
  Acc_cv <- numeric(K)
  
  for (i in 1:K) {
    # Split D_a into D_a* and D_a**
    test_indices <- which(folds == i)
    D_a_star <- D_a[test_indices,]
    D_a_star_star <- D_a[-test_indices,]
    
    # Estimate density ratio on D_a**
    ratios <- estimate_r(
      as.matrix(D_a_star_star[D_a_star_star$Group == 0, -c(which(names(D_a_star_star) %in% c("Y", "Group")))]), 
      as.matrix(D_a_star[D_a_star$Group == 0, -c(which(names(D_a_star) %in% c("Y", "Group")))]), 
      as.matrix(D_a_star_star[D_a_star_star$Group == 1, -c(which(names(D_a_star_star) %in% c("Y", "Group")))]), 
      as.matrix(D_a_star[D_a_star$Group == 1, -c(which(names(D_a_star) %in% c("Y", "Group")))]), 
      D_a_star_star$Y[D_a_star_star$Group == 0], 
      D_a_star$Y[D_a_star$Group == 0], 
      D_a_star_star$Y[D_a_star_star$Group == 1], 
      D_a_star$Y[D_a_star$Group == 1], 
      est.method, seed
    )
    
    r_X <- 1 / ratios$g22.est
    
    # Predict on D_a*
    predictions <- predict(classifier, newdata = D_a_star[, -which(names(D_a_star) == "Y")], type = "response")
    
    # Calculate A_1 and A_2
    A_1 <- mean(ifelse(predictions[D_a_star$Group == 0] > 0.5, 1, 0))
    A_2 <- mean(r_X * ifelse(predictions[D_a_star$Group == 1] <= 0.5, 1, 0))
    
    # Calculate variances
    var_1 <- var(ifelse(predictions[D_a_star$Group == 0] > 0.5, 1, 0))
    var_2 <- var(r_X * ifelse(predictions[D_a_star$Group == 1] <= 0.5, 1, 0))
    
    # Calculate Acc for this fold
    Acc_cv[i] <- sqrt(length(test_indices)) * (A_1 + A_2 - 1) / sqrt(var_1 + var_2)
  }
  
  # Calculate final test statistic
  Acc_final <- sum(Acc_cv) / sqrt(K)
  
  # Calculate p-value and make decision
  p_value <- 1 - pnorm(Acc_final)
  rejection <- as.integer(p_value < alpha)
  
  return(rejection)
}

bootstrap_MMD_test <- function(x1, x2, y1, y2, B=1000, prop=0.8, alpha=0.05, bandwidths=c(1), est.method='LL', seed=NULL){
    if(!is.null(seed)){
        set.seed(seed)
    }
    
    stopifnot(length(y1) == length(y2))
    total_n <- length(y1)
    
    # Split the data 
    n_est <- ceiling(total_n * prop) 
    est_indices <- 1:n_est
    test_indices <- (n_est+1):total_n
    
    x1_est <- x1[est_indices,,drop=FALSE]; x1_test <- x1[test_indices,,drop=FALSE]
    x2_est <- x2[est_indices,,drop=FALSE]; x2_test <- x2[test_indices,,drop=FALSE]
    y1_est <- y1[est_indices]; y1_test <- y1[test_indices]
    y2_est <- y2[est_indices]; y2_test <- y2[test_indices]
    
    ratios <- estimate_r(x1_est, x1_test, x2_est, x2_test,
                         y1_est, y1_test, y2_est, y2_test,
                         est.method, seed)
    r_X2 <- 1 / ratios$g22.est
    
    rejections <- numeric(length(bandwidths))
    for (i in seq_along(bandwidths)){
        bw <- bandwidths[i]
        if (bw == 'median'){
            h_x <- median.bandwidth(x1_test, x2_test)
            h_y <- median.bandwidth(matrix(y1_test), matrix(y2_test))
        } else {
            h_x <- as.numeric(bw)
            h_y <- as.numeric(bw)
        }
        
        boot_results <- bootstrap_MMD(x1_test, x2_test, y1_test, y2_test,
                                      h_x=h_x, h_y=h_y, r_X=r_X2, B=B, seed=seed) 
        
        obs_stat <- boot_results$obs_stat
        bootstrap_stats <- boot_results$bootstrap_stats
        
        # Calculate the p-values as the proportion of bootstrap statistics
        # that are greater than or equal to the observed statistic
        p_value <- mean(bootstrap_stats >= obs_stat)
        
        rejections[i] <- as.integer(p_value < alpha) 
    }
    return(rejections)
}


# Conditional Independence Test (CIT)
GCM_test <- function(x1, x2, y1, y2, alpha=0.05, epsilon=NULL, regr.method=ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1=TRUE, seed=NULL){
    if (!is.null(seed)) {
        set.seed(seed)
    }
    rejection <- 0
    n1 <- length(y1); n2 <- length(y2)
    
    if (alg1) {
        alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
        tilde_n1 <- alg1_res$tilde_n1
        tilde_n2 <- alg1_res$tilde_n2
        tilde_n <- tilde_n1 + tilde_n2
        
        if (tilde_n1 > n1 || tilde_n2 > n2) {
            print("tilde_n1 > n1 || tilde_n2 > n2")
            rejection <- 0
            return(rejection)
        } else {
            x1_sample <- x1[sample(1:nrow(x1), tilde_n1, replace = FALSE), , drop = FALSE]
            x2_sample <- x2[sample(1:nrow(x2), tilde_n2, replace = FALSE), , drop = FALSE]
            y1_sample <- y1[sample(1:length(y1), tilde_n1, replace = FALSE)]
            y2_sample <- y2[sample(1:length(y2), tilde_n2, replace = FALSE)]
            X_merged <- rbind(x1_sample, x2_sample)
            Y_merged <- c(y1_sample, y2_sample)
            Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
            
            merged_indices <- sample(1:tilde_n, size = tilde_n, replace = FALSE)
            X_merged <- X_merged[merged_indices, , drop = FALSE]
            Y_merged <- Y_merged[merged_indices]
            Z_merged <- Z_merged[merged_indices]
            gcm.pvalue <- gcm_test_binary(X=Y_merged, Y=Z_merged, Z=X_merged, reg_method=regr.method, binary_reg_method = binary.regr.method, seed=seed)
        }
    } else{
        X <- rbind(x1, x2)
        Y <- c(y1, y2)
        Z <- c(rep(0, n1), rep(1, n2))
        gcm.pvalue <- gcm_test_binary(X=Y, Y=Z, Z=X, reg_method=regr.method, binary_reg_method = binary.regr.method, seed=seed)
    }
    
    rejection <- ifelse(gcm.pvalue < alpha, 1, 0)
    return(rejection)
}
PCM_test <- function(x1, x2, y1, y2, alpha = 0.05, epsilon = NULL, regr.method = ranger_reg_method, binary.regr.method = ranger_reg_method_binary, alg1 = TRUE, seed=NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    rejection <- 0
    n1 <- length(y1); n2 <- length(y2)
    
    if (alg1) {
        alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
        tilde_n1 <- alg1_res$tilde_n1
        tilde_n2 <- alg1_res$tilde_n2
        # tilde_n1 <- min(tilde_n1, tilde_n2)
        # tilde_n2 <- tilde_n1
        tilde_n <- tilde_n1 + tilde_n2
        
        if (tilde_n1 > n1 || tilde_n2 > n2) {
            rejection <- 0
            print("tilde_n1 > n1 || tilde_n2 > n2")
            return(rejection)
        } else {
            idx1 <- sample.int(nrow(x1), tilde_n1)              
            x1_sample <- x1[idx1, , drop = FALSE]
            y1_sample <- y1[idx1]
            idx2 <- sample.int(nrow(x2), tilde_n2)
            x2_sample <- x2[idx2, , drop = FALSE]
            y2_sample <- y2[idx2]
            X_merged <- rbind(x1_sample, x2_sample)
            Y_merged <- c(y1_sample, y2_sample)
            Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
            
            merged_indices <- sample(1:tilde_n, size = tilde_n, replace = FALSE)
            X_merged <- X_merged[merged_indices, , drop = FALSE]
            Y_merged <- Y_merged[merged_indices]
            Z_merged <- Z_merged[merged_indices]
            pcm.pvalue <- pcm_test_binary(Y = Z_merged, X = Y_merged, Z = X_merged, reg_method = regr.method, binary_reg_method = binary.regr.method, seed=seed)
        }
    } else {
        X <- rbind(x1, x2)
        Y <- c(y1, y2)
        Z <- c(rep(0, n1), rep(1, n2))
        pcm.pvalue <- pcm_test_binary(X = Y, Y = Z, Z = X, reg_method = regr.method, binary_reg_method = binary.regr.method, seed=seed)
    }
    
    rejection <- ifelse(pcm.pvalue < alpha, 1, 0)
    return(rejection)
}
WGSC_test <- function(x1, x2, y1, y2, alpha=0.05, epsilon=NULL, regr.method=xgboost_reg_method, binary.regr.method = xgboost_reg_method_binary, alg1=TRUE, seed=NULL){
    if (!is.null(seed)) {
        set.seed(seed)
    }
    rejection <- 0
    n1 <- length(y1); n2 <- length(y2)
    
    if (alg1) {
        alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
        tilde_n1 <- alg1_res$tilde_n1
        tilde_n2 <- alg1_res$tilde_n2
        tilde_n <- tilde_n1 + tilde_n2
        
        if (tilde_n1 > n1 || tilde_n2 > n2) {
            print("tilde_n1 > n1 || tilde_n2 > n2")
            rejection <- 0
            return(rejection)
        } else {
            x1_sample <- x1[sample(1:nrow(x1), tilde_n1, replace = FALSE), , drop = FALSE]
            x2_sample <- x2[sample(1:nrow(x2), tilde_n2, replace = FALSE), , drop = FALSE]
            y1_sample <- y1[sample(1:length(y1), tilde_n1, replace = FALSE)]
            y2_sample <- y2[sample(1:length(y2), tilde_n2, replace = FALSE)]
            X_merged <- rbind(x1_sample, x2_sample)
            Y_merged <- c(y1_sample, y2_sample)
            Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
            
            merged_indices <- sample(1:tilde_n, size = tilde_n, replace = FALSE)
            X_merged <- X_merged[merged_indices, , drop = FALSE]
            Y_merged <- Y_merged[merged_indices]
            Z_merged <- Z_merged[merged_indices]
            wgsc.pvalue <- wgsc_binary(Y=Z_merged, X=Y_merged, Z=X_merged, reg_method=regr.method, binary_reg_method = binary.regr.method, seed=seed)
        }
    } else{
        X <- rbind(x1, x2)
        Y <- c(y1, y2)
        Z <- c(rep(0, n1), rep(1, n2))
        wgsc.pvalue <- wgsc_binary(Y=Z, X=Y, Z=X, reg_method=regr.method, binary_reg_method = binary.regr.method, seed=seed)
    }
    if (wgsc.pvalue < alpha){
        rejection <- 1
    }
    return(rejection)
}
RCIT_test <- function(x1, x2, y1, y2, alpha=0.05, epsilon=NULL, alg1=TRUE, seed=NULL, ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    rejection <- 0
    n1 <- length(y1); n2 <- length(y2)
    
    if (alg1 == TRUE) {
        alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
        tilde_n1 <- alg1_res$tilde_n1
        tilde_n2 <- alg1_res$tilde_n2
        # tilde_n1 <- min(tilde_n1, tilde_n2)
        # tilde_n2 <- tilde_n1
        tilde_n <- tilde_n1 + tilde_n2
        if (tilde_n1 > n1 || tilde_n2 > n2) {
            print("tilde_n1 > n1 || tilde_n2 > n2")
            rejection <- 0
            return(rejection)
        } else {
            x1_sample <- x1[sample(1:nrow(x1), tilde_n1, replace = FALSE), , drop = FALSE]
            x2_sample <- x2[sample(1:nrow(x2), tilde_n2, replace = FALSE), , drop = FALSE]
            y1_sample <- y1[sample(1:length(y1), tilde_n1, replace = FALSE)]
            y2_sample <- y2[sample(1:length(y2), tilde_n2, replace = FALSE)]
            X_merged <- rbind(x1_sample, x2_sample)
            Y_merged <- c(y1_sample, y2_sample)
            Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
            
            merged_indices <- sample(1:tilde_n, size=tilde_n, replace = FALSE)
            X_merged <- X_merged[merged_indices, , drop = FALSE]
            Y_merged <- Y_merged[merged_indices]
            Z_merged <- Z_merged[merged_indices]
            rcit.pvalue <- RCIT(x=Y_merged, y=Z_merged, z=X_merged, approx = "lpd4", num_f = 100, num_f2 = 5, seed=seed)$p
        }
    } else {  
        X <- rbind(x1, x2)
        Y <- c(y1, y2)
        Z <- c(rep(0, n1), rep(1, n2))
        rcit.pvalue <- RCIT(x=Y, y=Z, z=X, approx = "lpd4", num_f = 100, num_f2 = 5, seed=seed)$p
    }
    
    rejection <- ifelse(rcit.pvalue < alpha, 1, 0)
    return(rejection)
}



RCoT_test <- function(x1, x2, y1, y2, alpha=0.05, epsilon=NULL, alg1=TRUE, seed=NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0
  n1 <- length(y1); n2 <- length(y2)
  
  if (alg1) {
    alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
    tilde_n1 <- alg1_res$tilde_n1
    tilde_n2 <- alg1_res$tilde_n2
    tilde_n <- tilde_n1 + tilde_n2
    
    if (tilde_n1 > n1 || tilde_n2 > n2) {
      print("tilde_n1 > n1 || tilde_n2 > n2")
      rejection <- 0
      return(rejection)
    } else {
      x1_sample <- x1[sample(1:nrow(x1), tilde_n1, replace = FALSE), , drop = FALSE]
      x2_sample <- x2[sample(1:nrow(x2), tilde_n2, replace = FALSE), , drop = FALSE]
      y1_sample <- y1[sample(1:length(y1), tilde_n1, replace = FALSE)]
      y2_sample <- y2[sample(1:length(y2), tilde_n2, replace = FALSE)]
      X_merged <- rbind(x1_sample, x2_sample)
      Y_merged <- c(y1_sample, y2_sample)
      Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
      
      merged_indices <- sample(1:tilde_n, size = tilde_n, replace = FALSE)
      X_merged <- X_merged[merged_indices, , drop = FALSE]
      Y_merged <- Y_merged[merged_indices]
      Z_merged <- Z_merged[merged_indices]
      rcot.pvalue <- RCoT(x=Y_merged, y=Z_merged, z=X_merged, approx = "lpd4", num_f = 100, num_f2 = 5, seed = seed)$p
    }
  } else {
    X <- rbind(x1, x2)
    Y <- c(y1, y2)
    Z <- c(rep(0, n1), rep(1, n2))
    rcot.pvalue <- RCoT(x=Y, y=Z, z=X, approx = "lpd4", num_f = 100, num_f2 = 5, seed = seed)$p
  }
  
  rejection <- ifelse(rcot.pvalue < alpha, 1, 0)
  return(rejection)
}
KCI_test <- function(x1, x2, y1, y2, alpha=0.05, epsilon=NULL, GP=TRUE, alg1 = TRUE, seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0
  n1 <- length(y1); n2 <- length(y2)
  
  if (alg1) {
    alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
    tilde_n1 <- alg1_res$tilde_n1
    tilde_n2 <- alg1_res$tilde_n2
    tilde_n <- tilde_n1 + tilde_n2
    
    if (tilde_n1 > n1 || tilde_n2 > n2) {
      print("tilde_n1 > n1 || tilde_n2 > n2")
      rejection <- 0
      return(rejection)
    } else {
      x1_sample <- x1[sample(1:nrow(x1), tilde_n1, replace = FALSE), , drop = FALSE]
      x2_sample <- x2[sample(1:nrow(x2), tilde_n2, replace = FALSE), , drop = FALSE]
      y1_sample <- y1[sample(1:length(y1), tilde_n1, replace = FALSE)]
      y2_sample <- y2[sample(1:length(y2), tilde_n2, replace = FALSE)]
      X_merged <- rbind(x1_sample, x2_sample)
      Y_merged <- c(y1_sample, y2_sample)
      Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
      
      merged_indices <- sample(1:tilde_n, size = tilde_n, replace = FALSE)
      X_merged <- X_merged[merged_indices, , drop = FALSE]
      Y_merged <- Y_merged[merged_indices]
      Z_merged <- Z_merged[merged_indices]
      kci.pvalue <- kci_test(X = Y_merged, Y = Z_merged, Z = X_merged, GP=GP)
    }
  } else {
    X <- rbind(x1, x2)
    Y <- c(y1, y2)
    Z <- c(rep(0, n1), rep(1, n2))
    kci.pvalue <- kci_test(X = Y, Y = Z, Z = X, GP=GP)
  }
  rejection <- ifelse(kci.pvalue < alpha, 1, 0)
  return(rejection)
}
WGCM_test <- function(x1, x2, y1, y2, alpha=0.05, epsilon=NULL, regr.method=lm_reg_method, binary.regr.method = lm_reg_method_binary, alg1=TRUE, seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  rejection <- 0
  n1 <- length(y1); n2 <- length(y2)
  
  if (alg1) {
    alg1_res <- apply_alg1(x1, x2, y1, y2, seed, epsilon)
    tilde_n1 <- alg1_res$tilde_n1
    tilde_n2 <- alg1_res$tilde_n2
    tilde_n <- tilde_n1 + tilde_n2
    
    if (tilde_n1 > n1 || tilde_n2 > n2) {
      print("tilde_n1 > n1 || tilde_n2 > n2")
      rejection <- 0
      return(rejection)
    } else {
      x1_sample <- x1[sample(1:nrow(x1), tilde_n1, replace = FALSE), , drop = FALSE]
      x2_sample <- x2[sample(1:nrow(x2), tilde_n2, replace = FALSE), , drop = FALSE]
      y1_sample <- y1[sample(1:length(y1), tilde_n1, replace = FALSE)]
      y2_sample <- y2[sample(1:length(y2), tilde_n2, replace = FALSE)]
      X_merged <- rbind(x1_sample, x2_sample)
      Y_merged <- c(y1_sample, y2_sample)
      Z_merged <- c(rep(0, tilde_n1), rep(1, tilde_n2))
      
      merged_indices <- sample(1:tilde_n, size = tilde_n, replace = FALSE)
      X_merged <- X_merged[merged_indices, , drop = FALSE]
      Y_merged <- Y_merged[merged_indices]
      Z_merged <- Z_merged[merged_indices]
      wgcm.pvalue <- wGCM_fix_binary(X=Y_merged, Y=Z_merged, Z=X_merged, reg_method=regr.method, binary_reg_method = binary.regr.method, seed=seed)
    }
  } else{
    X <- rbind(x1, x2)
    Y <- c(y1, y2)
    Z <- c(rep(0, n1), rep(1, n2))
    wgcm.pvalue <- wGCM_fix_binary(X=Y, Y=Z, Z=X, reg_method=regr.method, binary_reg_method = binary.regr.method, seed=seed)
  }
  
  rejection <- ifelse(wgcm.pvalue < alpha, 1, 0)
  return(rejection)
}