#' @title Comparative Simulation of BlockMMD and LinearMMD
#' @description Simulates and plots BlockMMD and LinearMMD statistics for conditional two-sample testing
#'              under various block sizes and hypotheses.

rm(list=ls())
library(tmvtnorm)
library(CVST)
source("./experiments/utils.R") # contains MMDb, MMDl, estimate_r

# Function to generate covariates
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

# Function to generate response
generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, 1, 1, rep(0, dim(x)[2] - 4))
    mean_shift <- if (is_null) 0 else 0.5
    y <- f0 + epsilon + mean_shift
    return(y)
}

simulate_and_plot <- function(n = 2000, n_sim = 1000, seed = 1203, 
                              hypothesis = c("null", "alternative"),
                              prop = 0.5, est_method = "LL") {
    hypothesis <- match.arg(hypothesis)
    d <- 10
    is_null_hyp <- (hypothesis == "null")
    title_main <- ifelse(is_null_hyp, "Null", "Alternative")
    
    block_powers <- c(5/16, 7/16)
    
    for (p in block_powers) {
        B_size <- max(2, floor(n^p))
        cat("Running block size power:", p, "-> Block size =", B_size, "(hypothesis:", hypothesis, ")\n")
        
        # initialize statistics
        stats_block <- numeric(n_sim)
        stats_linear <- numeric(n_sim)
        
        for (i in seq_len(n_sim)) {
            set.seed(seed + i)
            x1 <- generate_data(n, d, group = 1)
            y1 <- generate_y(x1, is_null = TRUE)
            set.seed(seed + i + n_sim)
            x2 <- generate_data(n, d, group = 2)
            y2 <- generate_y(x2, is_null = is_null_hyp)
            
            # split
            n_train <- ceiling(n * prop)
            test_idx <- (n_train + 1):n
            
            x11 <- x1[1:n_train,,drop=FALSE]; x12 <- x1[test_idx,,drop=FALSE]
            y11 <- y1[1:n_train]; y12 <- y1[test_idx]
            x21 <- x2[1:n_train,,drop=FALSE]; x22 <- x2[test_idx,,drop=FALSE]
            y21 <- y2[1:n_train]; y22 <- y2[test_idx]
            
            ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est_method, seed+i)
            r_X <- 1 / ratios$g22.est
            
            h_x <- h_y <- 1
            
            stats_block[i] <- MMDb(x12, x22, y12, y22, B_size, h_x, h_y, r_X, seed+i)
            stats_linear[i] <- MMDl(x12, x22, y12, y22, h_x, h_y, r_X, seed+i)
        }
        
        stats_block <- stats_block[is.finite(stats_block)]
        stats_linear <- stats_linear[is.finite(stats_linear)]
        
        # Plot
        par(mfrow = c(2,2), mar = c(4,4,2,1), oma = c(0,0,3,0))
        
        hist(stats_block, freq = FALSE, breaks = 30, main = "BlockMMD", xlab = "value", col = "lightblue")
        if (length(stats_block) > 1) curve(dnorm(x, mean(stats_block), sd(stats_block)), add=TRUE, col="red", lwd=2)
        
        qqnorm(stats_block, main = "BlockMMD Q-Q Plot"); qqline(stats_block, col = "red", lwd = 2)
        
        hist(stats_linear, freq = FALSE, breaks = 30, main = "LinearMMD", xlab = "value", col = "lightgreen")
        if (length(stats_linear) > 1) curve(dnorm(x, mean(stats_linear), sd(stats_linear)), add=TRUE, col="blue", lwd=2)
        
        qqnorm(stats_linear, main = "LinearMMD Q-Q Plot"); qqline(stats_linear, col = "blue", lwd = 2)
        
        mtext(paste0(title_main, " (block size = ", B_size, ")"), outer = TRUE, cex = 1.5, font = 2)
    }
}

# Example: run both cases
simulate_and_plot(hypothesis = "null")
simulate_and_plot(hypothesis = "alternative")