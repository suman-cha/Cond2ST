rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
    library(MASS)
    library(pbapply)
    library(dplyr)
    library(ggplot2)
    library(data.table)
})
source("./experiments/all_tests.R")

tag <- "real_high_dim_bandwidths_comparion"

cur_wd <- getwd()
file_path <- file.path(cur_wd, "experiments", "real_examples", "data", "superconductivity.csv")
data <- read.csv(file_path)

X <- as.matrix(data[, !names(data) %in% c("critical_temp")])
Y <- data[, "critical_temp"]

# Remove missing data
X <- na.omit(X)
Y <- Y[complete.cases(Y)]

normalize <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}

X_norm <- apply(X, 2, normalize)
Y_norm <- normalize(Y)

sample_data <- function(X, Y, n, is_null = TRUE, is_x1 = TRUE) {
    if (is_x1) {
        X_idx <- sample(1:nrow(X), nrow(X)%/%2, replace = FALSE)
        X1 <- X[X_idx, , drop = FALSE]
        Y_subset <- Y[X_idx]
        x_idx <- sample(1:nrow(X1), n, replace=FALSE)
        x <- X1[x_idx,,drop=FALSE]
    } else {
        feature_to_bias <- X[, 1]  
        prob <- dnorm(feature_to_bias, 0, 1)
        prob <- prob / sum(prob)  
        X_idx <- sample(1:nrow(X), nrow(X)%/%2, replace = FALSE, prob = prob)
        X2 <- X[X_idx, , drop = FALSE]
        Y_subset <- Y[X_idx]  
        x_idx <- sample(1:nrow(X2), n, replace=FALSE)
        x <- X2[x_idx,,drop=FALSE]
    }
    
    if (is_null) {
        y <- sample(Y_subset, size=n, replace=FALSE)
    } else {
        u <- if (is_x1) dunif(Y_subset, 0, 1) else exp(-Y_subset)
        u <- u / sum(u)
        y <- sample(Y_subset, size = n, prob = u, replace = FALSE)
    }
    
    return(list(x = x, y = y))
}

# Parameters
n_values <- c(200,400,800,1200,1600,2000)
n_sims <- 500
bandwidths <- c(0.1, 0.5, 1.0, "median", 5, 10)
estimators <- c("LL", "KLR")
results_list <- list()

test_functions <- list(
    LinearMMD_test = LinearMMD_test,
    CV_LinearMMD_test = CV_LinearMMD_test
)

# Simulation loop
for(n in n_values) {
    for(is_null in c(FALSE, TRUE)) {
        h_label <- if(is_null) "Null" else "Alternative"
        
        for(test_name in names(test_functions)) {
            for(est in estimators) {
                if(test_name == "LinearMMD_test" || test_name == "CV_LinearMMD_test") {
                    result_matrix <- pbapply::pbsapply(1:n_sims, function(sim) {
                        seed <- 1203 + sim
                        set.seed(seed)
                        
                        d1 <- sample_data(X_norm, Y_norm, n, is_null, TRUE)
                        set.seed(seed + n_sims)
                        d2 <- sample_data(X_norm, Y_norm, n, is_null, FALSE)
                        
                        test_args <- list(
                            x1 = d1$x,
                            x2 = d2$x,
                            y1 = d1$y,
                            y2 = d2$y,
                            est.method = est,
                            seed = seed,
                            bandwidths = bandwidths
                        )
                        do.call(test_functions[[test_name]], test_args)
                    })
                    
                    # Process results for each bandwidth
                    for(bw_idx in seq_along(bandwidths)) {
                        bw <- bandwidths[bw_idx]
                        mean_result <- mean(result_matrix[bw_idx,])
                        
                        results_list[[length(results_list)+1]] <- data.table(
                            test_name = test_name,
                            n = n,
                            h_label = h_label,
                            estimator = est,
                            bandwidth = if(is.character(bw)) bw else as.character(bw),
                            rejection_rate = mean_result
                        )
                    }
                    
                    cat("[Test]", test_name, "| n:", n, "| Estimator:", est, "|", h_label,
                        "| Rejection Rates:", paste(round(rowMeans(result_matrix),3), collapse=", "), "\n",
                        strrep("-",80), "\n")
                }
            }
        }
    }
}

# Save results
results_dt <- rbindlist(results_list)
filename <- paste0("results/simulation_results_", tag, ".csv")
fwrite(results_dt, filename, row.names=FALSE)
cat("Results saved to", filename, "\n")
