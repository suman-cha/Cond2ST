# =============================================================================
# R/mmd.R -- MMD kernel functions and test statistics.
#
# Implements three estimators of the squared MMD between (X, Y) ~ P_1 and
# (X, Y) ~ P_2 after weighting the second sample by an estimated density
# ratio r_X = p_1(X) / p_2(X):
#   * MMDl           : linear-time U-statistic (low variance, fast)
#   * MMDb           : block U-statistic (variance/cost tradeoff)
#   * quadratic_MMD  : quadratic U-statistic + parametric bootstrap
#
# All three accept the same (x12, x22, y12, y22, r_X, h_x, h_y) signature.
# When src/mmd_kernels.cpp is available the Rcpp backend is used; otherwise
# the pure-R implementation is invoked. The two paths agree to within
# IEEE-754 round-off.
# =============================================================================

# ── Rcpp accelerated kernels ──────────────────────────────────────────────────
# Try to load C++ acceleration if available. Falls back to pure R.
.mmd_use_cpp <- FALSE
for (.mmd_cpp_dir in c("src", file.path("..", "src"))) {
    .mmd_cpp_file <- file.path(.mmd_cpp_dir, "mmd_kernels.cpp")
    if (file.exists(.mmd_cpp_file)) {
        if (!exists("mmd_linear_cpp", mode = "function"))
            tryCatch(Rcpp::sourceCpp(.mmd_cpp_file), error = function(e) NULL)
        .mmd_use_cpp <- exists("mmd_linear_cpp", mode = "function")
        break
    }
}
rm(.mmd_cpp_dir, .mmd_cpp_file)
# ─────────────────────────────────────────────────────────────────────────────

# Gaussian kernel
gaussian.kernel <- function(x, y = NULL, h=1) {
  if (is.null(y)) {
    res <- (exp(-0.5 * (x/h)^2) / (h * sqrt(2 * pi)))
  } else {
    dist <- sum((x - y)^2)
    res <- exp(-dist / (2 * h^2))
  }
  return(res)
}

median.bandwidth <- function(x, y){
    dists <- as.vector(dist(rbind(x, y)))
    bw <- median(dists[dists > 0])
    return(bw)
}

#' Linear-time MMD U-statistic
#'
#' Pairs adjacent observations within each split to form a linear-time
#' (\eqn{O(n)}) unbiased estimator of the squared MMD, weighted by the
#' supplied density ratio \code{r_X}. Returns the studentized statistic.
#'
#' @param x12,x22 Hold-out covariate matrices for groups 1 and 2.
#' @param y12,y22 Hold-out response vectors (same length as the rows of
#'   \code{x12}, \code{x22}).
#' @param h_x,h_y Gaussian-kernel bandwidths for the X and Y blocks.
#' @param r_X Estimated marginal density ratio evaluated on the second
#'   group's hold-out covariates.
#' @param seed Optional integer seed (only matters when the Rcpp backend
#'   draws auxiliary randomness; the linear estimator itself is
#'   deterministic given the data).
#' @return Studentized MMD statistic (rejects when large).
MMDl <- function(x12, x22, y12, y22, h_x=1, h_y=1, r_X, seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (.mmd_use_cpp) {
    return(mmd_linear_cpp(x12, x22, y12, y22, r_X, h_x, h_y))
  }
  # x12, x22, y12, y22 are splitted data 2n -> n
  stopifnot(length(y12) == length(y22))
  n <- length(y12)
  m <- floor(n/2)
  S_hat_values <- numeric(m)

  for(i in 1:m){
    k_zz <- gaussian.kernel(x12[i, ], x12[i+m, ], h_x) * gaussian.kernel(y12[i], y12[i+m], h_y)
    k_ww <- gaussian.kernel(x22[i, ], x22[i+m, ], h_x) * gaussian.kernel(y22[i], y22[i+m], h_y)
    k_wz <- gaussian.kernel(x22[i, ], x12[i+m, ], h_x) * gaussian.kernel(y22[i], y12[i+m], h_y)
    k_zw <- gaussian.kernel(x12[i, ], x22[i+m, ], h_x) * gaussian.kernel(y12[i], y22[i+m], h_y)

    S_hat_values[i] <- k_zz + r_X[i]*r_X[i+m]*k_ww - r_X[i]*k_wz - r_X[i+m]*k_zw
  }
  S_bar <- mean(S_hat_values)
  sigma_hat <- sum((S_hat_values - S_bar)^2) / (m - 1)
  MMDl_hat2 <- sqrt(m) * S_bar / sqrt(sigma_hat)
  return(MMDl_hat2)
}

#' Block MMD U-statistic
#'
#' Splits the hold-out sample into \code{floor(n / B_size)} disjoint blocks
#' of size \code{B_size}, computes a U-statistic per block, and returns the
#' studentized average. Trades variance for cost between \code{MMDl}
#' (\code{B_size = 2}) and \code{quadratic_MMD} (single block).
#'
#' @inheritParams MMDl
#' @param B_size Block size; must be at least 2 and at most \eqn{n / 2}.
#' @return Studentized MMD statistic.
MMDb <- function(x12, x22, y12, y22, B_size, h_x=1, h_y=1, r_X, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (.mmd_use_cpp) {
        return(mmd_block_cpp(x12, x22, y12, y22, r_X, B_size, h_x, h_y))
    }

    n <- length(y12)
    S <- floor(n / B_size)

    if (S < 2) stop("Block size too large for given sample size")

    MMD_values <- numeric(S)

    for (b in 1:S){
        idx <- ((b-1)*B_size+1):(b*B_size)
        xb1 <- x12[idx,,drop=FALSE]
        xb2 <- x22[idx,,drop=FALSE]
        yb1 <- y12[idx]
        yb2 <- y22[idx]
        rb <- r_X[idx]

        sum_k <- 0
        # (i, j) where i < j
        for (i in 1:(B_size-1)) {
            for (j in (i+1):B_size) {
                k_zz <- gaussian.kernel(xb1[i, ], xb1[j, ], h_x) * gaussian.kernel(yb1[i], yb1[j], h_y)
                k_ww <- gaussian.kernel(xb2[i, ], xb2[j, ], h_x) * gaussian.kernel(yb2[i], yb2[j], h_y)
                k_wz <- gaussian.kernel(xb2[i, ], xb1[j, ], h_x) * gaussian.kernel(yb2[i], yb1[j], h_y)
                k_zw <- gaussian.kernel(xb1[i, ], xb2[j, ], h_x) * gaussian.kernel(yb1[i], yb2[j], h_y)

                sum_k <- sum_k + k_zz + rb[i]*rb[j]*k_ww - rb[i]*k_wz - rb[j]*k_zw

            }
        }
        n_pairs <- B_size*(B_size - 1) / 2
        # n_pairs <- (B_size - 1) / 2
        MMD_values[b] <- sum_k / n_pairs
    }

    S_bar <- mean(MMD_values)
    sigma_hat <- sqrt(var(MMD_values))
    stat <- ifelse(sigma_hat >0, sqrt(S)*S_bar/sigma_hat, 0)
    return(stat)
}

# Regime II (gamma = 1) wrapper: rescales the raw U-stat and its bootstrap
# draws by sqrt(n(n-1)) as required by Theorem 3.2 of the paper. The
# rejection decision is unaffected (monotone transform) but the returned
# statistics are now on the regime-II scale.
quadratic_MMD_regime2 <- function(x12, x22, y12, y22, h_x = 1, h_y = 1, r_X, B,
                                  seed = NULL) {
    n <- length(y12)
    stopifnot(length(y12) == length(y22))
    if (n < 2) {
        stop("Quadratic MMD requires at least two test observations per group")
    }

    raw <- quadratic_MMD(
        x12 = x12, x22 = x22, y12 = y12, y22 = y22,
        h_x = h_x, h_y = h_y, r_X = r_X, B = B, seed = seed
    )

    scale <- sqrt(n * (n - 1))
    list(
        obs_stat        = scale * raw$obs_stat,
        bootstrap_stats = scale * raw$bootstrap_stats,
        obs_u_stat      = raw$obs_stat
    )
}

#' Quadratic-time MMD with parametric bootstrap calibration
#'
#' Computes the quadratic U-statistic estimator of the squared MMD
#' (\eqn{O(n^2)} cost, lowest variance) and calibrates rejection by
#' generating \code{B} parametric-bootstrap replicates under the null.
#'
#' @inheritParams MMDl
#' @param B Bootstrap replication count for null calibration (e.g. 199 or
#'   299).
#' @return A list with the observed statistic and bootstrap p-value.
quadratic_MMD <- function(x12, x22, y12, y22, h_x=1, h_y=1, r_X, B, seed=NULL){
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (.mmd_use_cpp) {
        return(bootstrap_mmd_cpp(x12, x22, y12, y22, r_X, h_x, h_y, B))
    }

    n <- length(y12)
    stopifnot(length(y12) == length(y22))

    H_hat <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n){
        for (j in 1:n) {
            if (i != j) {
                k_zz <- gaussian.kernel(x12[i,,drop=FALSE], x12[j,,drop=FALSE],
                                        h_x) * gaussian.kernel(y12[i], y12[j], h_y)
                k_ww <- gaussian.kernel(x22[i,,drop=FALSE], x22[j,,drop=FALSE],
                                        h_x) * gaussian.kernel(y22[i], y22[j], h_y)
                k_wz <- gaussian.kernel(x22[i,,drop=FALSE], x12[j,,drop=FALSE],
                                        h_x) * gaussian.kernel(y22[i], y12[j], h_y)
                k_zw <- gaussian.kernel(x12[i,,drop=FALSE], x22[j,,drop=FALSE],
                                        h_x) * gaussian.kernel(y12[i], y22[j], h_y)
                H_hat[i, j] <- k_zz + r_X[i]*r_X[j]*k_ww - r_X[i]*k_wz - r_X[j]*k_zw
            }
        }
    }
    obs_stat <- sum(H_hat) / (n*(n-1))
    bootstrap_stats <- numeric(B)

    for (b in 1:B){
        # Generate n i.i.d. Gaussian random variables
        W <- rnorm(n)

        # Calculate the wild boostrap statistic
        wild_bootstrap_sum <- sum(outer(W,W) * H_hat)
        bootstrap_stats[b] <- wild_bootstrap_sum / (n*(n-1))
    }
    return(list(obs_stat=obs_stat, bootstrap_stats=bootstrap_stats))
}
