rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(pbapply)
    library(MASS)
    library(dplyr)
    library(data.table)
    library(R.utils)
    library(microbenchmark)
})

tag <- "S1U_computationl_costs"
source("./experiments/all_tests.R")

max_time <- 180

# Data generation
generate_data <- function(n, p, group) {
    mu <- if (group == 1) c(1, 1, -1, -1, rep(0, p - 4)) else rep(0, p)
    sigma <- diag(1, p)
    mvrnorm(n, mu = mu, Sigma = sigma)
}

generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- rt(n, df = sigma)
    f0 <- x %*% c(1, -1, -1, 1, rep(0, dim(x)[2] - 4))
    mean_shift <- if (is_null) 0 else .5
    f0 + epsilon + mean_shift
}

# Tests
drt_test_functions <- list(
    LinearMMD_test      = LinearMMD_test,
    CV_LinearMMD_test   = CV_LinearMMD_test,
    CLF_test            = CLF_test,
    CV_CLF_test         = CV_CLF_test,
    CP_test             = CP_test,
    debiased_test       = debiased_test,
    BlockMMD_test       = BlockMMD_test,
    CV_BlockMMD_test    = CV_BlockMMD_test,
    bootstrap_MMD_test  = bootstrap_MMD_test
)

cit_test_functions <- list(
    RCIT_test = RCIT_test,
    PCM_test  = PCM_test,
    GCM_test  = GCM_test,
    WGSC_test = WGSC_test
)

# Settings
n_values <- c(200, 500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d <- 10
results_list <- list()

# Main simulation loop
for (n in n_values) {
    for (is_null in c(TRUE, FALSE)) {
        h_label <- if (is_null) "Null" else "Alternative"
        
        for (test_type in c("DRT", "CIT")) {
            test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
            
            for (test_name in names(test_functions)) {
                cat("[Start] Test:", test_name, "| n:", n, "|", h_label, "\n")
                
                # Simulation + timing
                result_matrix <- pbapply::pbsapply(1:n_sims, function(sim) {
                    seed <- 1203 + sim
                    set.seed(seed)
                    
                    x1 <- generate_data(n, d, group = 1)
                    y1 <- generate_y(x1, is_null = TRUE)
                    set.seed(seed + n_sims)
                    x2 <- generate_data(n, d, group = 2)
                    y2 <- generate_y(x2, is_null = is_null)
                    
                    test_args <- list(x1, x2, y1, y2, seed = seed)
                    
                    t0 <- proc.time()
                    out <- do.call(test_functions[[test_name]], test_args)
                    t1 <- proc.time()
                    time_elapsed <- (t1 - t0)[["elapsed"]]
                    time_elapsed <- pmin(time_elapsed, max_time)
                    
                    c(result = out, time = time_elapsed)
                }, simplify = "array")
                
                mean_result <- mean(result_matrix["result", ])
                mean_time <- mean(result_matrix["time", ])
                
                results_list[[length(results_list) + 1]] <- data.table(
                    test_type = test_type,
                    test_name = test_name,
                    n = n,
                    h_label = h_label,
                    rejection_rate = mean_result,
                    avg_time_sec = mean_time
                )
                
                cat("[Test]", test_name, "| n:", n, "|", h_label, "| Rejection Rate:", mean_result,
                    "| Avg Time (sec):", round(mean_time, 4), "\n", strrep("-", 80), "\n")
            }
        }
    }
}

# Save all results
results_dt <- rbindlist(results_list)
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")