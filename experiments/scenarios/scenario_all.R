# -----------------------------------------------------------------------------
# Heavy-tailed / non-Gaussian X experiments for conditional two-sample testing
# BOUNDED vs UNBOUNDED density-ratio regimes
# -----------------------------------------------------------------------------
rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(mvtnorm)   # rmvt
    library(pbapply)
    library(data.table)
})

# Tag used in filenames
tag <- "S_heavyX"

# Your test wrappers (DRT and CIT) must be sourced here
# and should define: LinearMMD_test, CV_LinearMMD_test, CLF_test, CV_CLF_test,
# CP_test, debiased_test, BlockMMD_test, CV_BlockMMD_test,
# RCIT_test, PCM_test, GCM_test, WGSC_test
source("./experiments/all_tests.R")

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# --------------------------
# Utilities
# --------------------------

# i.i.d. Laplace sampler (independent across coordinates)
rlaplace_mat <- function(n, p, mu = rep(0, p), scale = 1) {
    # Laplace(0, b): density = (1/(2b)) exp(-|x|/b)
    z <- rexp(n * p, rate = 1/scale)
    s <- sample(c(-1, 1), n * p, replace = TRUE)
    mat <- matrix(s * z, nrow = n, ncol = p)
    sweep(mat, 2, mu, FUN = "+")
}

# Generic accept-reject truncation to a hyper-rectangle [lower, upper]
sample_truncated <- function(draw_fun, n, lower, upper, batch = max(1000, 5*n), ...) {
    p <- length(lower)
    out <- matrix(NA_real_, nrow = n, ncol = p)
    filled <- 0L
    while (filled < n) {
        cand <- draw_fun(batch, p = p, ...)
        keep <- rep(TRUE, nrow(cand))
        for (j in 1:p) {
            keep <- keep & (cand[, j] >= lower[j] & cand[, j] <= upper[j])
        }
        if (any(keep)) {
            add <- cand[keep, , drop = FALSE]
            k <- min(n - filled, nrow(add))
            out[(filled + 1):(filled + k), ] <- add[1:k, , drop = FALSE]
            filled <- filled + k
        }
    }
    out
}

# Linear signal for Y
make_beta <- function(p) {
    # same pattern as your Simulation-1 (first four active)
    b <- rep(0, p); b[1:4] <- c(-1, -1, 1, 1)
    b
}

# Y generator: H0: same conditional, H1: constant shift (keeps focus on X design)
generate_y <- function(x, is_null = TRUE, sigma_eps = 1, delta = 0.5) {
    n <- nrow(x)
    f0 <- as.numeric(x %*% make_beta(ncol(x)))
    eps <- rnorm(n, 0, sigma_eps)
    f0 + eps + if (is_null) 0 else delta
}

# --------------------------
# X-scenarios
# --------------------------
# Each returns list(x1=..., x2=...) with dimension n x p.
# Bounded-ratio cases use truncation to a common compact set.
# Unbounded-ratio cases differ in tail-decay (df/scale), making sup_x r_X(x) infinite.

generate_X <- function(n, p, scenario, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    mu_shift <- c(1, 1, -1, -1, rep(0, max(0, p - 4)))  # same as baseline
    I_p <- diag(1, p)
    B <- 1  # truncation bound (wide to keep accept rate high)
    lower <- rep(-B, p); upper <- rep(B, p)
    
    switch(
        scenario,
        
        # ------------------ BOUNDED density ratio ------------------
        # Both groups truncated to same compact support; densities continuous & bounded away
        # from 0 on the set -> density ratio bounded.
        
        "B_t_trunc" = {
            # Heavy-tailed (Student t) but truncated to a common box
            draw_t1 <- function(m, p) rmvt(m, sigma = I_p, df = 3, delta = mu_shift)
            draw_t2 <- function(m, p) rmvt(m, sigma = I_p, df = 3, delta = rep(0, p))
            x1 <- sample_truncated(draw_t1, n, lower, upper)
            x2 <- sample_truncated(draw_t2, n, lower, upper)
            list(x1 = x1, x2 = x2)
        },
        
        "B_laplace_trunc" = {
            # Non-Gaussian (Laplace) with mean shift; truncated to a common box
            draw_L1 <- function(m, p) rlaplace_mat(m, p, mu = mu_shift, scale = 1)
            draw_L2 <- function(m, p) rlaplace_mat(m, p, mu = rep(0, p), scale = 1)
            x1 <- sample_truncated(draw_L1, n, lower, upper)
            x2 <- sample_truncated(draw_L2, n, lower, upper)
            list(x1 = x1, x2 = x2)
        },
        
        # ------------------ UNBOUNDED density ratio ------------------
        # Different tail parameters (df/scale) -> likelihood ratio diverges in the tails.
        
        "U_t_df" = {
            # Same location/scale, different df (heavier vs lighter tails)
            # r_X(x) = f_t,df=3(x) / f_t,df=10(x) is unbounded as ||x|| -> ∞.
            x1 <- rmvt(n, sigma = I_p, df = 3,  delta = rep(0, p))
            x2 <- rmvt(n, sigma = I_p, df = 10, delta = rep(0, p))
            list(x1 = x1, x2 = x2)
        },
        
        "U_laplace_scale" = {
            # Different Laplace scales: scale-1 vs scale-0.5 => unbounded ratio in the tails
            x1 <- rlaplace_mat(n, p, mu = rep(0, p), scale = 1.0)
            x2 <- rlaplace_mat(n, p, mu = rep(0, p), scale = 0.5)
            list(x1 = x1, x2 = x2)
        },
        
        # Fallback
        {
            stop("Unknown scenario: ", scenario)
        }
    )
}

# --------------------------
# Test function registries
# --------------------------
drt_test_functions <- list(
    LinearMMD_test   = LinearMMD_test,
    CV_LinearMMD_test = CV_LinearMMD_test,
    CLF_test         = CLF_test,
    CV_CLF_test      = CV_CLF_test,
    CP_test          = CP_test,
    debiased_test    = debiased_test,
    BlockMMD_test    = BlockMMD_test,
    CV_BlockMMD_test = CV_BlockMMD_test
)

cit_test_functions <- list(
    RCIT_test = RCIT_test,
    PCM_test  = PCM_test,
    GCM_test  = GCM_test,
    WGSC_test = WGSC_test
)

# --------------------------
# Experiment grid
# --------------------------
n_values   <- c(200, 500, 1000, 2000)
n_sims     <- 500
alpha      <- 0.05
p          <- 10
estimators <- c("LL", "KLR")     # for DRT’s classification-based ratio estimation
scenarios  <- c("B_t_trunc", "B_laplace_trunc", "U_t_df", "U_laplace_scale")

results_all <- list()

# --------------------------
# Main simulation loop
# --------------------------
for (sc in scenarios) {
    cat("\n========== Scenario:", sc, "==========\n")
    results_list <- list()
    
    for (n in n_values) {
        for (is_null in c(TRUE, FALSE)) {
            h_label <- if (is_null) "Null" else "Alternative"
            
            for (test_type in c("DRT", "CIT")) {
                test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
                
                for (test_name in names(test_functions)) {
                    
                    if (test_type == "DRT") {
                        # Loop over density-ratio estimators (LL/KLR)
                        for (est in estimators) {
                            res <- pbapply::pbsapply(1:n_sims, function(sim) {
                                seed <- 1203 + sim
                                set.seed(seed)
                                
                                Xpair <- generate_X(n, p, scenario = sc, seed = seed)
                                x1 <- Xpair$x1; x2 <- Xpair$x2
                                
                                y1 <- generate_y(x1, is_null = TRUE)
                                y2 <- generate_y(x2, is_null = is_null)
                                
                                test_args <- list(x1, x2, y1, y2, est.method = est, seed = seed)
                                do.call(test_functions[[test_name]], test_args)
                            }, simplify = "array")
                            
                            rr <- mean(res)
                            results_list[[length(results_list) + 1]] <- data.table(
                                scenario = sc,
                                test_type = test_type,
                                test_name = test_name,
                                estimator = est,
                                n = n,
                                h_label = h_label,
                                alpha = alpha,
                                rejection_rate = rr
                            )
                            cat("[", sc, "]", test_name, "| n:", n, "| Estimator:", est, "|", h_label,
                                "| Rejection Rate:", sprintf("%.3f", rr), "\n", strrep("-", 80), "\n")
                        }
                        
                    } else {
                        # CIT family: Algorithm 1 (stable construction) + epsilon as in your real-data code
                        res <- pbapply::pbsapply(1:n_sims, function(sim) {
                            seed <- 1203 + sim
                            set.seed(seed)
                            
                            Xpair <- generate_X(n, p, scenario = sc, seed = seed)
                            x1 <- Xpair$x1; x2 <- Xpair$x2
                            
                            y1 <- generate_y(x1, is_null = TRUE)
                            y2 <- generate_y(x2, is_null = is_null)
                            
                            epsilon <- 1 / sqrt(log(n))
                            test_args <- list(x1, x2, y1, y2, alg1 = TRUE, epsilon = epsilon, seed = seed)
                            do.call(test_functions[[test_name]], test_args)
                        }, simplify = "array")
                        
                        rr <- mean(res)
                        results_list[[length(results_list) + 1]] <- data.table(
                            scenario = sc,
                            test_type = test_type,
                            test_name = test_name,
                            estimator = NA_character_,
                            n = n,
                            h_label = h_label,
                            alpha = alpha,
                            rejection_rate = rr
                        )
                        cat("[", sc, "]", test_name, "| n:", n, "|", h_label,
                            "| Rejection Rate:", sprintf("%.3f", rr), "\n", strrep("-", 80), "\n")
                    }
                }
            }
        }
    }
    
    results_dt <- rbindlist(results_list)
    out_file <- sprintf("results/simulation_results_%s_%s.csv", tag, sc)
    fwrite(results_dt, out_file)
    cat(">>> Saved:", out_file, "\n")
    results_all[[sc]] <- results_dt
}

# Optionally, bind all scenarios together:
all_dt <- rbindlist(results_all, use.names = TRUE, fill = TRUE)
fwrite(all_dt, sprintf("results/simulation_results_%s_ALL.csv", tag))
cat(">>> Saved:", sprintf("results/simulation_results_%s_ALL.csv", tag), "\n")
