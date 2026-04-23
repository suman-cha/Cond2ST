# =============================================================================
# scenarios.R -- All 10 data generating processes for the cdtst paper
#
# Each scenario provides:
#   generate_data_S{N}{U/B}(n, group) -> n x d matrix of covariates X
#   generate_y_S{N}{U/B}(x, ...)     -> length-n vector of responses Y
#
# Scenario overview:
#   S1: Location shift, d=10
#   S2: Heteroscedastic variance, d=10
#   S3: Post-nonlinear transform, d=10
#   S4: Student-t covariates, d=10
#   S5: Beta covariates, d=5
#   U = unbounded, B = bounded (truncated)
#
# Dependencies: MASS (S1-S3), tmvtnorm (S1B-S3B)
# =============================================================================


# =============================================================================
# Shared constants & helpers for S1-S3
# =============================================================================

.D_S13   <- 10L
.MU1_S13 <- c(1, 1, -1, -1, rep(0, .D_S13 - 4))
.BETA_S1 <- c(1, -1, -1, 1, rep(0, .D_S13 - 4))

.gen_x_normal <- function(n, group) {
    mu <- if (group == 1L) .MU1_S13 else rep(0, .D_S13)
    MASS::mvrnorm(n, mu = mu, Sigma = diag(1, .D_S13))
}

.gen_x_truncnorm <- function(n, group) {
    lb <- rep(-0.5, .D_S13); ub <- rep(0.5, .D_S13)
    mu <- if (group == 1L) .MU1_S13 else rep(0, .D_S13)
    tmvtnorm::rtmvnorm(n, mean = mu, sigma = diag(1, .D_S13),
                        lower = lb, upper = ub, algorithm = "gibbs")
}

.g_heteroscedastic <- function(x, rho = 10) {
    x_centered <- sweep(x, 2, 0.5)
    norm_sq <- rowSums(x_centered^2)
    rho * (1 + exp(-norm_sq / 64))
}


# =============================================================================
# S1 -- Location shift
# Y = X*beta + t(2) + delta; delta=0 (null), 0.5 (alt)
# =============================================================================

generate_data_S1U <- .gen_x_normal
generate_data_S1B <- .gen_x_truncnorm

generate_y_S1 <- function(x, is_null = TRUE) {
    n  <- nrow(x)
    f0 <- as.numeric(x %*% .BETA_S1)
    f0 + rt(n, df = 2) + if (is_null) 0 else 0.5
}
generate_y_S1U <- generate_y_S1
generate_y_S1B <- generate_y_S1


# =============================================================================
# S2 -- Heteroscedastic variance
# Null: Y ~ N(X*1_p, 10^2); Alt: Y ~ N(X*beta_alt, g(X))
# g(x) = rho * (1 + exp(-||x - 0.5||^2 / 64)), rho=10
# =============================================================================

generate_data_S2U <- .gen_x_normal
generate_data_S2B <- .gen_x_truncnorm

generate_y_S2 <- function(x, is_null = TRUE) {
    n <- nrow(x); p <- ncol(x)
    beta <- if (is_null) rep(1, p) else c(rep(1, p - 1), 0)
    mean_X <- as.numeric(x %*% beta)
    if (is_null) {
        rnorm(n, mean = mean_X, sd = 10)
    } else {
        rnorm(n, mean = mean_X, sd = sqrt(.g_heteroscedastic(x)))
    }
}
generate_y_S2U <- generate_y_S2
generate_y_S2B <- generate_y_S2


# =============================================================================
# S3 -- Post-nonlinear transform
# Y = f(X*1_p + 2*eps); null: f=cos; alt: f in {x, x^2, x^3, sin, tanh}
# =============================================================================

generate_data_S3U <- .gen_x_normal
generate_data_S3B <- .gen_x_truncnorm

generate_y_S3 <- function(x, is_null = TRUE) {
    n <- nrow(x); d <- ncol(x)
    eps <- rnorm(n)
    linear_part <- as.numeric(x %*% rep(1, d)) + 2 * eps
    if (is_null) {
        cos(linear_part)
    } else {
        transforms <- list(function(z) z, function(z) z^2, function(z) z^3,
                           function(z) sin(z), function(z) tanh(z))
        tf <- sample(transforms, 1)[[1]]
        tf(linear_part)
    }
}
generate_y_S3U <- generate_y_S3
generate_y_S3B <- generate_y_S3


# =============================================================================
# S4U -- Student-t covariates, unbounded
# X: t(df) + mu, df=5 (group1) or 4 (group2)
# Y = delta + X*beta + t(2); delta=0 under null, 0.5 for group2 under alt
# =============================================================================

generate_data_S4U <- function(n, group) {
    d   <- 10
    mu1 <- c(1, 1, -1, -1, rep(0, max(0, d - 4)))
    df  <- if (group == 1L) 5 else 4
    mu  <- if (group == 1L) mu1 else rep(0, d)
    X   <- matrix(rt(n * d, df = df), nrow = n, ncol = d)
    sweep(X, 2, mu, "+")
}


# =============================================================================
# S4B -- Student-t covariates, bounded (rejection sampling to [-R, R]^d)
# =============================================================================

.rprod_t_trunc <- function(n, p, df, R) {
    acc <- NULL
    while (is.null(acc) || nrow(acc) < n) {
        m   <- max(2 * n, 1000L)
        Xb  <- matrix(rt(m * p, df = df), nrow = m, ncol = p)
        keep <- apply(abs(Xb) <= R, 1L, all)
        if (any(keep)) {
            Xk  <- Xb[keep, , drop = FALSE]
            acc <- if (is.null(acc)) Xk else rbind(acc, Xk)
        }
    }
    acc[seq_len(n), , drop = FALSE]
}

generate_data_S4B <- function(n, group) {
    d       <- 10
    R_trunc <- 3
    mu1     <- c(1, 1, -1, -1, rep(0, max(0, d - 4)))
    df      <- if (group == 1L) 5 else 4
    mu      <- if (group == 1L) mu1 else rep(0, d)
    X       <- .rprod_t_trunc(n, d, df, R_trunc)
    sweep(X, 2, mu, "+")
}

generate_y_S4 <- function(x, is_null = TRUE, group = 1L) {
    d     <- 10
    beta  <- c(1, -1, 1, -1, rep(0, max(0, d - 4)))
    n     <- nrow(x)
    delta <- if (is_null) 0 else if (group == 2L) 0.5 else 0
    drop(delta + x %*% beta + rt(n, 2))
}
generate_y_S4U <- generate_y_S4
generate_y_S4B <- generate_y_S4


# =============================================================================
# S5U -- Beta covariates, unbounded (full [0,1] support)
# X_ij ~ Beta(a, 2); a=0.5 (group1) or 2 (group2)
# Y = sin(2*pi*x1) + 0.5*log(1+10*x2) + 0.3*(x3-0.5)^2 + N(0,1) + delta
# =============================================================================

generate_data_S5U <- function(n, group) {
    d <- 5
    a <- if (group == 1L) 0.5 else 2
    b <- 2
    matrix(rbeta(n * d, shape1 = a, shape2 = b), nrow = n, ncol = d)
}


# =============================================================================
# S5B -- Beta covariates, bounded (truncated to [0.1, 0.9]^d)
# =============================================================================

.rprod_beta_trunc <- function(n, p, a, b, lower = 0.1, upper = 0.9) {
    acc <- NULL
    while (is.null(acc) || nrow(acc) < n) {
        m <- max(2 * n, 2000L)
        X <- matrix(rbeta(m * p, shape1 = a, shape2 = b), nrow = m, ncol = p)
        keep <- apply(X >= lower & X <= upper, 1L, all)
        if (any(keep)) {
            Xk  <- X[keep, , drop = FALSE]
            acc <- if (is.null(acc)) Xk else rbind(acc, Xk)
        }
    }
    acc[seq_len(n), , drop = FALSE]
}

generate_data_S5B <- function(n, group) {
    d         <- 5
    eps_trunc <- 0.1
    a <- if (group == 1L) 0.5 else 2
    b <- 2
    .rprod_beta_trunc(n, d, a, b, lower = eps_trunc, upper = 1 - eps_trunc)
}

generate_y_S5 <- function(x, is_null = TRUE) {
    n      <- nrow(x)
    y_mean <- sin(2 * pi * x[, 1]) + 0.5 * log(1 + 10 * x[, 2]) +
              0.3 * (x[, 3] - 0.5)^2
    y_mean + rnorm(n, 0, 1) + if (is_null) 0 else 0.5
}
generate_y_S5U <- generate_y_S5
generate_y_S5B <- generate_y_S5


# =============================================================================
# Convenience lookup table
# =============================================================================

#' Scenario lookup table
#'
#' Maps a scenario tag (\code{"S1U"}, ..., \code{"S5B"}) to its data
#' generators. Each entry is a list with fields:
#' \itemize{
#'   \item \code{tag}: the scenario name.
#'   \item \code{gen_data(n, group)}: returns an \eqn{n \times d} covariate
#'     matrix from the marginal of group \eqn{j \in \{1, 2\}}.
#'   \item \code{gen_y(x, is_null = TRUE, ...)}: returns the response vector
#'     under either the null \eqn{H_0} (\code{is_null = TRUE}) or the
#'     alternative \eqn{H_1}.
#' }
#' Use \code{SCENARIO_META} for descriptive metadata.
SCENARIOS <- list(
    S1U = list(tag = "S1U", gen_data = generate_data_S1U, gen_y = generate_y_S1U),
    S1B = list(tag = "S1B", gen_data = generate_data_S1B, gen_y = generate_y_S1B),
    S2U = list(tag = "S2U", gen_data = generate_data_S2U, gen_y = generate_y_S2U),
    S2B = list(tag = "S2B", gen_data = generate_data_S2B, gen_y = generate_y_S2B),
    S3U = list(tag = "S3U", gen_data = generate_data_S3U, gen_y = generate_y_S3U),
    S3B = list(tag = "S3B", gen_data = generate_data_S3B, gen_y = generate_y_S3B),
    S4U = list(tag = "S4U", gen_data = generate_data_S4U, gen_y = generate_y_S4U),
    S4B = list(tag = "S4B", gen_data = generate_data_S4B, gen_y = generate_y_S4B),
    S5U = list(tag = "S5U", gen_data = generate_data_S5U, gen_y = generate_y_S5U),
    S5B = list(tag = "S5B", gen_data = generate_data_S5B, gen_y = generate_y_S5B)
)

SCENARIO_META <- list(
    S1U = list(d = 10, x_dist = "Normal",            y_type = "location_shift"),
    S1B = list(d = 10, x_dist = "Truncated Normal",   y_type = "location_shift"),
    S2U = list(d = 10, x_dist = "Normal",            y_type = "heteroscedastic"),
    S2B = list(d = 10, x_dist = "Truncated Normal",   y_type = "heteroscedastic"),
    S3U = list(d = 10, x_dist = "Normal",            y_type = "nonlinear_transform"),
    S3B = list(d = 10, x_dist = "Truncated Normal",   y_type = "nonlinear_transform"),
    S4U = list(d = 10, x_dist = "Student-t",          y_type = "location_shift"),
    S4B = list(d = 10, x_dist = "Truncated Student-t", y_type = "location_shift"),
    S5U = list(d = 5,  x_dist = "Beta",              y_type = "nonlinear_additive"),
    S5B = list(d = 5,  x_dist = "Truncated Beta",     y_type = "nonlinear_additive")
)

cat("[scenarios.R] Loaded 10 scenarios + SCENARIOS lookup\n")
