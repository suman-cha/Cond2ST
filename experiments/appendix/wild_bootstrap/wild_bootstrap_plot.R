# -----------------------------------------------
# Wild-bootstrap MMD: normality visualization
# -----------------------------------------------
rm(list=ls())
suppressPackageStartupMessages({
    library(tmvtnorm)   # rtmvnorm, dtmvnorm
})
# 'estimate_r' 등 유틸은 기존과 동일하게 사용
#  - 반드시 이 파일 안에 estimate_r가 있어야 합니다.
source("./experiments/utils.R") # contains estimate_r, etc.

# ---------- Gaussian kernel & pairwise helpers ----------
gaussian.kernel <- function(x, y = NULL, h = 1) {
    if (is.null(y)) {
        res <- (exp(-0.5 * (x/h)^2) / (h * sqrt(2 * pi)))
    } else {
        dist <- sum((x - y)^2)
        res <- exp(-dist / (2 * h^2))
    }
    return(res)
}

sq_dist_mat <- function(A, B = NULL) {
    # returns ||a_i - b_j||^2 matrix
    if (is.null(B)) B <- A
    A <- as.matrix(A); B <- as.matrix(B)
    An <- rowSums(A^2); Bn <- rowSums(B^2)
    # broadcasting
    D2 <- outer(An, Bn, "+") - 2 * (A %*% t(B))
    # numeric guard
    D2[D2 < 0] <- 0
    return(D2)
}
gaussian_gram_X <- function(A, B = NULL, h = 1) {
    D2 <- sq_dist_mat(A, B)
    exp(- D2 / (2 * h^2))
}
gaussian_gram_y <- function(a, b = NULL, h = 1) {
    if (is.null(b)) b <- a
    D2 <- (outer(as.numeric(a), as.numeric(b), "-"))^2
    exp(- D2 / (2 * h^2))
}

# ---------- Synthetic data generators (same as your code) ----------
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
generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, 1, 1, rep(0, ncol(x) - 4))
    mean_shift <- if (is_null) 0 else 0.5
    y <- f0 + epsilon + mean_shift
    return(as.numeric(y))
}

# ---------- True r_X(x) under the truncated MVN model ----------
true_rX <- function(x, mu1, mu2, lb, ub, sigma = NULL) {
    if (is.null(sigma)) sigma <- diag(1, ncol(x))
    f1 <- dtmvnorm(x, mean = mu1, sigma = sigma, lower = lb, upper = ub)
    f2 <- dtmvnorm(x, mean = mu2, sigma = sigma, lower = lb, upper = ub)
    # guard against division-by-zero
    f2[pmax(f2, 0) == 0] <- .Machine$double.eps
    r <- f1 / f2
    return(as.numeric(r))
}

# ---------- Build H-matrix for quadratic-time MMD kernel ----------
# Z_i = (V_i^(1), V_i^(2)), where V_i^(1) = (x1_i, y1_i), V_i^(2) = (x2_i, y2_i)
build_H_matrix <- function(x1, y1, x2, y2, r_vec, h_x = 1, h_y = 1) {
    # all are length n after split; r_vec corresponds to x2's rows
    n <- nrow(x1)
    stopifnot(nrow(x2) == n, length(y1) == n, length(y2) == n, length(r_vec) == n)
    
    # Gram matrices for X and Y
    Kx_11 <- gaussian_gram_X(x1, x1, h_x)
    Kx_22 <- gaussian_gram_X(x2, x2, h_x)
    Kx_21 <- gaussian_gram_X(x2, x1, h_x)  # rows: x2, cols: x1
    Kx_12 <- t(Kx_21)
    
    Ky_11 <- gaussian_gram_y(y1, y1, h_y)
    Ky_22 <- gaussian_gram_y(y2, y2, h_y)
    Ky_21 <- gaussian_gram_y(y2, y1, h_y)
    Ky_12 <- t(Ky_21)
    
    K_11 <- Kx_11 * Ky_11
    K_22 <- Kx_22 * Ky_22
    K_21 <- Kx_21 * Ky_21
    K_12 <- Kx_12 * Ky_12
    
    R_outer <- outer(r_vec, r_vec, "*")   # r_i * r_j
    R_row   <- outer(r_vec, rep(1, n), "*")  # r_i * 1_j
    R_col   <- t(R_row)                      # 1_i * r_j
    
    H <- K_11 + R_outer * K_22 - R_row * K_21 - R_col * K_12
    diag(H) <- 0
    return(H)
}

# ---------- Single wild-bootstrap draw given H ----------
wild_boot_stat <- function(H) {
    n <- nrow(H)
    W <- sample(c(-1, 1), size = n, replace = TRUE)
    as.numeric(t(W) %*% (H %*% W)) / (n * (n - 1))
}

# ---------- Main simulation & plotting ----------
simulate_and_plot_wild <- function(
        n = 1000, n_sim = 1000, seed = 1203,
        hypothesis = c("null", "alternative"),
        prop = 0.5, est_method = "LL",
        h_x = 1, h_y = 1, d = 10
) {
    hypothesis <- match.arg(hypothesis)
    is_null_hyp <- (hypothesis == "null")
    title_main <- ifelse(is_null_hyp, "Null", "Alternative")
    
    # model parameters used in data generator
    mu1 <- c(1, 1, -1, -1, rep(0, d - 4))
    mu2 <- rep(0, d)
    lb  <- rep(-.5, d)
    ub  <- rep(.5,  d)
    
    stats_oracle <- numeric(n_sim)
    stats_est    <- numeric(n_sim)
    
    for (i in seq_len(n_sim)) {
        set.seed(seed + i)
        # sample group 1 and 2
        x1 <- generate_data(n, d, group = 1)
        y1 <- generate_y(x1, is_null = TRUE)                    # Y1 follows null structure
        set.seed(seed + i + n_sim)
        x2 <- generate_data(n, d, group = 2)
        y2 <- generate_y(x2, is_null = is_null_hyp)             # Y2: null/alt controlled
        
        # split train/test for density-ratio estimation (like your code)
        n_train <- ceiling(n * prop)
        test_idx <- (n_train + 1):n
        
        x11 <- x1[1:n_train, , drop=FALSE];  x12 <- x1[test_idx, , drop=FALSE]
        y11 <- y1[1:n_train];                 y12 <- y1[test_idx]
        x21 <- x2[1:n_train, , drop=FALSE];  x22 <- x2[test_idx, , drop=FALSE]
        y21 <- y2[1:n_train];                 y22 <- y2[test_idx]
        
        # (A) Oracle r_X at test X2
        rX_true <- true_rX(x22, mu1 = mu1, mu2 = mu2, lb = lb, ub = ub)
        
        # (B) Estimated r_X at test X2 via estimate_r (same convention as your code)
        ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est_method, seed + i)
        # In your code: r_X <- 1 / ratios$g22.est
        rX_est <- 1 / as.numeric(ratios$g22.est)
        
        # Build H and compute one wild-bootstrap draw per dataset
        H_orc <- build_H_matrix(x12, y12, x22, y22, rX_true, h_x, h_y)
        H_est <- build_H_matrix(x12, y12, x22, y22, rX_est,  h_x, h_y)
        
        stats_oracle[i] <- wild_boot_stat(H_orc)
        stats_est[i]    <- wild_boot_stat(H_est)
    }
    
    # remove non-finite (safety)
    stats_oracle <- stats_oracle[is.finite(stats_oracle)]
    stats_est    <- stats_est[is.finite(stats_est)]
    
    # standardized (z-scores) for visual comparison to N(0,1)
    z_orc <- (stats_oracle - mean(stats_oracle)) / sd(stats_oracle)
    z_est <- (stats_est    - mean(stats_est))    / sd(stats_est)
    
    # Plot
    op <- par(mfrow = c(2,2), mar = c(4,4,2,1), oma = c(0,0,3,0))
    on.exit(par(op))
    
    # Oracle r_X
    hist(z_orc, freq = FALSE, breaks = 30, main = "Oracle r: Histogram", xlab = "z",
         col = "lightblue", border = "white")
    curve(dnorm(x, 0, 1), add = TRUE, lwd = 2)
    qqnorm(z_orc, main = "Oracle r: Q-Q Plot"); qqline(z_orc, lwd = 2)
    
    # Estimated r_X
    hist(z_est, freq = FALSE, breaks = 30, main = "Estimated r: Histogram", xlab = "z",
         col = "lightgreen", border = "white")
    curve(dnorm(x, 0, 1), add = TRUE, lwd = 2)
    qqnorm(z_est, main = "Estimated r: Q-Q Plot"); qqline(z_est, lwd = 2)
    
    mtext(sprintf("%s – Wild-bootstrap MMD", title_main),
          outer = TRUE, cex = 1.4, font = 2)
    
    invisible(list(
        raw = list(oracle = stats_oracle, estimated = stats_est),
        z   = list(oracle = z_orc, estimated = z_est)
    ))
}

# ---------------- Example runs ----------------
# Null: normality 시각화
out_null <- simulate_and_plot_wild(
    n = 1000, n_sim = 1000, seed = 1203,
    hypothesis = "null", prop = 0.5, est_method = "LL",
    h_x = 1, h_y = 1, d = 10
)

# Alternative도 보고 싶으면 다음 주석 해제
# out_alt <- simulate_and_plot_wild(
#   n = 1000, n_sim = 1000, seed = 1203,
#   hypothesis = "alternative", prop = 0.5, est_method = "LL",
#   h_x = 1, h_y = 1, d = 10
# )
