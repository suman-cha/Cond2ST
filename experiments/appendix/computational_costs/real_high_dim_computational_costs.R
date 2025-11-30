rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(pbapply)
    library(MASS)
    library(dplyr)
    library(data.table)
    library(R.utils)
})

MAX_TIME <- 180
ERROR_LOG <- "error_log.txt"

source("./experiments/all_tests.R")
tag <- "real_high_dim_computational_costs"

cur_wd <- getwd()
file_path <- file.path(cur_wd, "experiments", "real_examples", "data", "superconductivity.csv")
data <- read.csv(file_path)

X <- as.matrix(data[, !names(data) %in% c("critical_temp")])
Y <- data[, "critical_temp"]

X <- na.omit(X)
Y <- Y[complete.cases(Y)]

normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}
X_norm <- apply(X, 2, normalize)
Y_norm <- normalize(Y)

sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
    if (is_x1) {
        X_idx <- sample(1:nrow(X), nrow(X) %/% 2, replace = FALSE)
        X1 <- X[X_idx, , drop = FALSE]
        Y_subset <- Y[X_idx]
        x_idx <- sample(1:nrow(X1), n, replace = FALSE)
        x <- X1[x_idx, , drop = FALSE]
    } else {
        feature_to_bias <- X[, 1]
        prob <- dnorm(feature_to_bias, 0, 1)
        prob <- prob / sum(prob)
        X_idx <- sample(1:nrow(X), nrow(X) %/% 2, replace = FALSE, prob = prob)
        X2 <- X[X_idx, , drop = FALSE]
        Y_subset <- Y[X_idx]
        x_idx <- sample(1:nrow(X2), n, replace = FALSE)
        x <- X2[x_idx, , drop = FALSE]
    }
    
    if (is_null) {
        y <- sample(Y_subset, size = n, replace = FALSE)
    } else {
        u <- if (is_x1) dunif(Y_subset, 0, 1) else exp(-Y_subset)
        u <- u / sum(u)
        y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
    }
    return(list(x = x, y = y))
}

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

n_values    <- c(200, 400, 800, 1200, 1600, 2000)
n_sims      <- 100
estimators  <- c("LL", "KLR")

reg_methods <- list(
    "lm"     = list(regr.method = "lm_reg_method",     binary.regr.method = "lm_reg_method_binary"),
    "ranger" = list(regr.method = "ranger_reg_method", binary.regr.method = "ranger_reg_method_binary"),
    "xgb"    = list(regr.method = "xgboost_reg_method",    binary.regr.method = "xgboost_reg_method_binary")
)

results_list <- list()
cat("", file = ERROR_LOG)

for (n in n_values) {
    for (is_null in c(TRUE, FALSE)) {
        h_label <- if (is_null) "Null" else "Alternative"
        
        for (test_type in c("DRT", "CIT")) {
            test_functions <- if (test_type == "DRT") drt_test_functions else cit_test_functions
            
            for (test_name in names(test_functions)) {
                if (test_type == "DRT") {
                    for (est in estimators) {
                        result_matrix <- pbapply::pbsapply(1:n_sims, function(sim) {
                            seed <- 1203 + sim
                            set.seed(seed)
                            d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
                            set.seed(seed + n_sims)
                            d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
                            
                            test_args <- list(
                                d1$x, d2$x, d1$y, d2$y,
                                est.method = est,
                                seed       = seed
                            )
                            
                            res <- tryCatch({
                                t0 <- Sys.time()
                                out <- withTimeout({
                                    do.call(test_functions[[test_name]], test_args)
                                }, timeout = MAX_TIME)
                                t1 <- Sys.time()
                                elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
                                
                                c(result = out, time = elapsed)
                            }, TimeoutException = function(ex) {
                                c(result = NA, time = MAX_TIME)
                            }, error = function(e) {
                                msg <- paste(Sys.time(), "[Error in DRT]", test_name, "|", e$message)
                                cat(msg, "\n", file = ERROR_LOG, append = TRUE)
                                c(result = NA, time = NA)
                            })
                            
                            return(res)
                        })
                        
                        valid_results <- result_matrix["result", !is.na(result_matrix["result", ])]
                        mean_result   <- ifelse(length(valid_results) > 0, mean(valid_results), NA)
                        mean_time     <- mean(as.numeric(result_matrix["time", ]), na.rm = TRUE)
                        
                        results_list[[length(results_list) + 1]] <- data.table(
                            test_type      = test_type,
                            test_name      = test_name,
                            n              = n,
                            h_label        = h_label,
                            estimator      = est,
                            rejection_rate = mean_result,
                            avg_time_sec   = mean_time
                        )
                        cat(sprintf("[Test] %s | n: %d | Estimator: %s | %s | Rejection Rate: %.4f | Avg Time (sec): %.6f\n",
                                    test_name, n, est, h_label, mean_result, mean_time))
                    }
                    
                } else {  
                    for (reg_name in names(reg_methods)) {
                        reg_config <- reg_methods[[reg_name]]
                        
                        result_matrix <- pbapply::pbsapply(1:n_sims, function(sim) {
                            seed <- 1203 + sim
                            set.seed(seed)
                            d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
                            set.seed(seed + n_sims)
                            d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
                            
                            epsilon <- 1 / sqrt(log(n))
                            test_args <- list(
                                d1$x, d2$x, d1$y, d2$y,
                                alg1                = TRUE,
                                epsilon             = epsilon,
                                regr.method         = reg_config$`regr.method`,
                                binary.regr.method  = reg_config$`binary.regr.method`,
                                seed                = seed
                            )
                            
                            res <- tryCatch({
                                t0 <- Sys.time()
                                out <- withTimeout({
                                    do.call(test_functions[[test_name]], test_args)
                                }, timeout = MAX_TIME)
                                t1 <- Sys.time()
                                elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
                                
                                c(result = out, time = elapsed)
                            }, TimeoutException = function(ex) {
                                c(result = NA, time = MAX_TIME)
                            }, error = function(e) {
                                msg <- paste(Sys.time(), "[Error in CIT]", test_name,
                                             "| Reg Method:", reg_name, "|", e$message)
                                cat(msg, "\n", file = ERROR_LOG, append = TRUE)
                                c(result = NA, time = NA)
                            })
                            
                            return(res)
                        })
                        
                        valid_results <- result_matrix["result", !is.na(result_matrix["result", ])]
                        mean_result   <- ifelse(length(valid_results) > 0, mean(valid_results), NA)
                        mean_time     <- mean(as.numeric(result_matrix["time", ]), na.rm = TRUE)
                        
                        results_list[[length(results_list) + 1]] <- data.table(
                            test_type      = test_type,
                            test_name      = test_name,
                            n              = n,
                            h_label        = h_label,
                            estimator      = reg_name,
                            rejection_rate = mean_result,
                            avg_time_sec   = mean_time
                        )
                        cat(sprintf("[Test] %s | n: %d | Reg Method: %s | %s | Rejection Rate: %.4f | Avg Time (sec): %.6f\n",
                                    test_name, n, reg_name, h_label, mean_result, mean_time))
                    }
                }
            }
        }
    }
}

results_dt <- rbindlist(results_list)
filename <- paste0("results/ablations/xgb_simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")