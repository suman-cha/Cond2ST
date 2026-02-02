rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(MASS)
  library(pbapply)
  library(data.table)
})
tag <- "S1U_split_ratio"
source("./experiments/all_tests.R")

# Data generation functions
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

# Test functions
drt_test_functions <- list(
  LinearMMD_test = LinearMMD_test,
  BlockMMD_test = BlockMMD_test,
  bootstrap_MMD_test = bootstrap_MMD_test,
  CLF_test = CLF_test
)

# Parameters
n_values <- c(500, 1000, 2000)
n_sims <- 500
alpha <- 0.05
d <- 10
split_ratios <- c(0.3, 0.5, 0.8)
results_list <- list()

# Simulation loop
for (n in n_values) {
  for (is_null in c(FALSE, TRUE)) {  # Start with Alternative (FALSE), then Null (TRUE)
    h_label <- if (is_null) "Null" else "Alternative"
    
    for (split_prop in split_ratios) {
      for (test_name in names(drt_test_functions)) {
        
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
          
          # Prepare test arguments based on test function
          if (test_name %in% c("LinearMMD_test", "BlockMMD_test", "bootstrap_MMD_test")) {
            test_args <- list(x1, x2, y1, y2, prop = split_prop, seed = seed)
          } else if (test_name == "CLF_test") {
            test_args <- list(x1, x2, y1, y2, split.prop = split_prop, seed = seed)
          }
          
          tryCatch({
            do.call(drt_test_functions[[test_name]], test_args)
          }, error = function(e) {
            cat("Error in", test_name, "with prop", split_prop, ":", e$message, "\n")
            return(NA)
          })
        }, simplify = "array")
        
        # Remove NA values
        result <- result[!is.na(result)]
        if (length(result) == 0) {
          mean_result <- NA
        } else {
          mean_result <- mean(result)
        }
        
        results_list[[length(results_list) + 1]] <- data.table(
          test_type = "DRT",
          test_name = test_name,
          n = n,
          h_label = h_label,
          split_prop = split_prop,
          rejection_rate = mean_result,
          n_successful = length(result)
        )
        
        # Print results
        cat("[Test]", test_name, "| n:", n, "| prop:", split_prop, "|", h_label, 
            "| Rejection Rate:", mean_result, "| Successful:", length(result), "/", n_sims, "\n", strrep("-", 80), "\n")
      }
    }
  }
}

# Combine all results into a single data.table
results_dt <- rbindlist(results_list)

# Define filename and save to CSV
filename <- paste0("results/simulation_results_", tag, ".csv")
if (!dir.exists("results")) {
  dir.create("results")
}
fwrite(results_dt, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")

