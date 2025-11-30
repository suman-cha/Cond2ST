rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(data.table)
})

# This script addresses Comment C.1 by comparing a CIT method under
# (a) subsampling (Alg.1 with epsilon) and (b) an oracle Z setting (no subsampling),
# and also reporting the proposed DRT's performance, across Scenarios 1–3.

source("./experiments/all_tests.R")

# ----------------------
# Scenario generators
# ----------------------

# Scenario 1 (matches experiments/scenarios/scenario1u.R)
sc1_generate_x <- function(n, p, group) {
    mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
    sigma <- diag(1, p)
    MASS::mvrnorm(n, mu = mu, Sigma = sigma)
}
sc1_generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, -1, 1, rep(0, ncol(x) - 4))
    mean_shift <- if (is_null) 0 else .5
    as.numeric(f0 + epsilon + mean_shift)
}

# Scenario 2 (matches experiments/scenarios/scenario2u.R)
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
    if (is_null) rnorm(n, mean = mean_X, sd = 10) else rnorm(n, mean = mean_X, sd = sqrt(var_X))
}

# Scenario 3 (matches experiments/scenarios/scenario3u.R)
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

# ----------------------
# Experiment configuration
# ----------------------

scenarios <- list(
    S1 = list(gen_x = sc1_generate_x, gen_y = sc1_generate_y),
    S2 = list(gen_x = sc2_generate_x, gen_y = sc2_generate_y),
    S3 = list(gen_x = sc3_generate_x, gen_y = sc3_generate_y)
)

n_values <- c(200, 500, 1000)
n_sims <- 300
alpha <- 0.05
p_dim <- 10

# epsilon chooser (default 1 / log(n) as used in utils)
epsilon_of_n <- function(n) 1 / log(n)

# Choose one DRT representative and one CIT representative
drt_method_name <- "BlockMMD_test"    # proposed DRT representative
cit_method_name <- "PCM_test"         # CIT representative

results <- list()

for (sc_tag in names(scenarios)) {
    gen_x <- scenarios[[sc_tag]]$gen_x
    gen_y <- scenarios[[sc_tag]]$gen_y
    
    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            eps_val <- epsilon_of_n(n)
            
            # Run simulations
            sim_res <- pbapply::pbsapply(1:n_sims, function(sim) {
                seed <- 1203 + sim
                set.seed(seed)
                
                # Generate groups P^(1), P^(2)
                x1 <- gen_x(n, p_dim, group = 1)
                y1 <- gen_y(x1, is_null = TRUE)
                
                set.seed(seed + n_sims)
                x2 <- gen_x(n, p_dim, group = 2)
                y2 <- gen_y(x2, is_null = is_null)
                
                # DRT (group-wise; no notion of subsampling)
                drt_rej <- do.call(drt_method_name, list(
                    x1 = x1, x2 = x2, y1 = y1, y2 = y2,
                    alpha = alpha, seed = seed, est.method = "LL",
                    bandwidths = "median"
                ))
                drt_rej <- as.integer(drt_rej[1])
                
                # CIT with subsampling (Alg.1, uses epsilon)
                cit_sub_rej <- do.call(cit_method_name, list(
                    x1 = x1, x2 = x2, y1 = y1, y2 = y2,
                    alpha = alpha, epsilon = eps_val, alg1 = TRUE, seed = seed
                ))
                
                # CIT with oracle Z (no subsampling)
                cit_oracle_rej <- do.call(cit_method_name, list(
                    x1 = x1, x2 = x2, y1 = y1, y2 = y2,
                    alpha = alpha, alg1 = FALSE, seed = seed
                ))
                
                c(drt = drt_rej, cit_sub = cit_sub_rej, cit_oracle = cit_oracle_rej)
            }, simplify = "matrix")
            
            # Aggregate
            drt_rate <- mean(sim_res["drt", ])
            cit_sub_rate <- mean(sim_res["cit_sub", ])
            cit_oracle_rate <- mean(sim_res["cit_oracle", ])
            
            results[[length(results) + 1]] <- data.table(
                scenario = sc_tag,
                n = n,
                h_label = h_label,
                method = c("DRT_BlockMMD", "CIT_PCM_subsample", "CIT_PCM_oracle"),
                rejection_rate = c(drt_rate, cit_sub_rate, cit_oracle_rate)
            )
            
            cat("[", sc_tag, "] n:", n, "|", h_label,
                "| DRT:", sprintf("%.3f", drt_rate),
                "| CIT_sub:", sprintf("%.3f", cit_sub_rate),
                "| CIT_oracle:", sprintf("%.3f", cit_oracle_rate),
                "\n",
                strrep("-", 80), "\n")
        }
    }
}

results_dt <- rbindlist(results)
outfile <- file.path("results", "comment_c1_oracle_vs_cit.csv")
fwrite(results_dt, outfile)
cat("Results saved to", outfile, "\n")


