#' Quadratic-time MMD test with DRPT star-sampler permutation
#' 
#' @param x1, x2, y1, y2 Numeric matrices/vectors for covariates and responses of two samples
#' @param prop Proportion for splitting (as in LinearMMD_test)
#' @param alpha Significance level
#' @param bandwidths Vector of bandwidths for the Gaussian kernel
#' @param est.method Method for density ratio estimation
#' @param H Number of permutations (default 199)
#' @param S Number of MCMC steps for star-sampler (default 50)
#' @param kernel Kernel function (default: Gaussian)
#' @param seed Random seed (optional)
#' @return Vector of rejections for each bandwidth

rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
    library(tmvtnorm)
    library(DRPT)
})
source("./experiments/utils.R")

MMDq <- function(V1, V2, r_X2, h_x, h_y, kernel) {
    n1 <- nrow(V1); n2 <- nrow(V2)
    term1 <- sum(outer(1:n1, 1:n1, Vectorize(function(i, j)
        if (i != j) kernel_joint(V1[i, ], V1[j, ], h_x, h_y, kernel) else 0))) / (n1 * (n1 - 1))
    term2 <- sum(outer(1:n2, 1:n2, Vectorize(function(i, j)
        if (i != j) r_X2[i] * r_X2[j] * kernel_joint(V2[i, ], V2[j, ], h_x, h_y, kernel) else 0))) / (n2 * (n2 - 1))
    term3 <- sum(outer(1:n1, 1:n2, Vectorize(function(i, j)
        r_X2[j] * kernel_joint(V1[i, ], V2[j, ], h_x, h_y, kernel)))) / (n1 * n2)
    
    return(term1 + term2 - 2 * term3)
}

kernel_joint <- function(a, b, h_x, h_y, kernel) {
    stopifnot(length(a) >= 2, length(b) >= 2)
    x_a <- a[-length(a)]
    x_b <- b[-length(b)]
    y_a <- a[length(a)]
    y_b <- b[length(b)]
    return(kernel(x_a, x_b, h_x) * kernel(y_a, y_b, h_y))
}



MMDq_test <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, bandwidths=c(1), 
                      est.method="LL", H=199, S=50, kernel=NULL, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    stopifnot(length(y1) == length(y2))
    total_sample_size <- length(y1)
    n <- ceiling(total_sample_size * prop)
    x11 <- x1[1:n, , drop=F]; x12 <- x1[-(1:n),,drop=F]
    y11 <- y1[1:n]; y12 <- y1[-(1:n)]
    x21 <- x2[1:n, , drop=F]; x22 <- x2[-(1:n),,drop=F]
    y21 <- y2[1:n]; y22 <- y2[-(1:n)]
    
    # Estimate density ratio
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    r_X2 <- 1/ratios$g22.est
    
    # Combine data for permutation
    V1 <- cbind(x12, as.numeric(y12))
    V2 <- cbind(x22, as.numeric(y22))
    n1 <- nrow(V1)
    n2 <- nrow(V2)
    Z <- rbind(V1, V2)
    r_full <- c(rep(1, n1), r_X2) # r=1 for group 1, r_X2 for group 2
    
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
        
        obs_stat <- MMDq(V1, V2, r_X2, h_x, h_y, kernel=kernel)
        
        # Permutations via star-sampler from DRPT package
        r_fun <- function(z) {
            # z는 starSampler 내부에서 combined data의 행렬 인덱스(1~n1+n2)로 전달됨
            r_full[z]
        }
        perms <- DRPT::starSampler(V1, V2, r = r_fun, H = H, S = S)
        perm_stats <- numeric(H)
        for (h in 1:H) {
            idx <- perms[[h+1]] # first is original, next H are permutations
            idx1 <- idx[1:n1]
            idx2 <- idx[(n1+1):(n1+n2)]
            V1_perm <- Z[idx1, , drop=F]
            V2_perm <- Z[idx2, , drop=F]
            r_X2_perm <- r_full[idx2]
            perm_stats[h] <- MMDq(V1_perm, V2_perm, r_X2_perm, h_x, h_y, kernel=kernel)
        }
        pval <- mean(c(obs_stat, perm_stats) >= obs_stat)
        rejections[i] <- as.integer(pval < alpha)
    }
    return(rejections)
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

# Function to generate response variable
generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, 1, 1, rep(0, dim(x)[2] - 4))
    mean_shift <- if (is_null) 0 else .5
    y <- f0 + epsilon + mean_shift
    return(y)
}

# Test functions
drt_test_functions <- list(MMDq_test = MMDq_test)
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
                    
                    test_args <- list(x1, x2, y1, y2, seed = seed, kernel=gaussian.kernel)
                    if (!is.function(test_args$kernel)) stop("kernel must be a function")
                    
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
