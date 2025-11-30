#' Create Publication-Quality Visualizations
#'
#' This script analyzes the simulation results and creates visualizations
#' showing the relationship between estimation errors and testing performance
#' for CIT methods (GCM, PCM, WGSC, RCIT)
#'
#' Visualizations:
#' 1. Estimation error vs sample size (4 panels for 4 methods)
#' 2. Test statistic variance vs sample size
#' 3. Scatter plots: Estimation error vs rejection rate
#' 4. Correlation heatmap: Estimation error vs test instability
#'
#' @references Lee, Cha, Kim (2024). arXiv:2410.16636
#' @author Conditional Two-Sample Testing Research Team
#' @date 2024-11-17

rm(list = ls())
set.seed(1203)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Load Results
# ------------------------------------------------------------------------------

results_file <- "results/cit_estimation_vs_testing.csv"

if (!file.exists(results_file)) {
  stop("Results file not found. Please run cit_estimation_vs_testing.R first.")
}

cat("Loading results from:", results_file, "\n")
results_dt <- fread(results_file)

# Separate aggregated and individual results
results_agg <- results_dt[is.na(sim_id) | type != "individual"]
results_ind <- results_dt[type == "individual"]

cat("Loaded", nrow(results_agg), "aggregated records and",
    nrow(results_ind), "individual simulation records\n\n")

# ------------------------------------------------------------------------------
# Compute Stability and Correlation Metrics
# ------------------------------------------------------------------------------

cat("Computing stability and correlation metrics...\n")

# Compute correlation between MSE and rejection for each configuration
correlation_results <- results_ind[!is.na(mse_total)] %>%
  group_by(scenario, n, h_label, method) %>%
  summarize(
    cor_mse_rejection = cor(mse_total, rejection, use = "complete.obs"),
    cor_mse_pvalue = cor(mse_total, p_value, use = "complete.obs"),
    mean_mse = mean(mse_total, na.rm = TRUE),
    sd_mse = sd(mse_total, na.rm = TRUE),
    cv_mse = sd(mse_total, na.rm = TRUE) / mean(mse_total, na.rm = TRUE),
    mean_rejection = mean(rejection, na.rm = TRUE),
    sd_rejection = sd(rejection, na.rm = TRUE),
    n_obs = .N,
    .groups = "drop"
  ) %>%
  as.data.table()

cat("Correlation analysis completed for",
    nrow(correlation_results), "configurations\n\n")

# Save correlation results
cor_file <- "results/cit_estimation_correlation_analysis.csv"
fwrite(correlation_results, cor_file)
cat("Correlation results saved to:", cor_file, "\n\n")

# ------------------------------------------------------------------------------
# Create Publication-Quality Figures Directory
# ------------------------------------------------------------------------------

fig_dir <- "figures/cit_estimation_performance"
if (!dir.exists("figures")) {
  dir.create("figures")
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Figure 1: Estimation Error vs Sample Size
# ------------------------------------------------------------------------------

cat("Creating Figure 1: Estimation Error vs Sample Size...\n")

# Filter for alternative hypothesis and non-NA MSE values
plot_data <- results_agg[h_label == "Alternative" & !is.na(mean_mse_total)]

# Create color palette
method_colors <- c("GCM" = "#E41A1C", "PCM" = "#377EB8",
                   "WGSC" = "#4DAF4A", "RCIT" = "#984EA3")

pdf(file.path(fig_dir, "fig1_mse_vs_sample_size.pdf"),
    width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.2, cex.axis = 1.1)

for (method_name in c("GCM", "PCM", "WGSC")) {
  method_data <- plot_data[method == method_name]
  
  # Plot for each scenario
  plot(NULL, xlim = c(min(plot_data$n), max(plot_data$n)),
       ylim = c(0, max(method_data$mean_mse_total, na.rm = TRUE) * 1.2),
       xlab = "Sample Size (n)", ylab = "Mean Squared Error",
       main = sprintf("%s Estimation Error", method_name),
       log = "x")
  
  scenarios <- unique(method_data$scenario)
  scenario_symbols <- c(S1 = 19, S2 = 17, S3 = 15)
  scenario_labels <- c(S1 = "Mean Shift", S2 = "Heteroscedastic",
                       S3 = "Nonlinear")
  
  for (sc in scenarios) {
    sc_data <- method_data[scenario == sc]
    sc_data <- sc_data[order(sc_data$n)]
    
    lines(sc_data$n, sc_data$mean_mse_total, col = method_colors[method_name],
          lwd = 2, type = "b", pch = scenario_symbols[sc], cex = 1.5)
    
    # Add error bars (±1 SE)
    arrows(sc_data$n,
           sc_data$mean_mse_total - sc_data$sd_mse_total / sqrt(500),
           sc_data$n,
           sc_data$mean_mse_total + sc_data$sd_mse_total / sqrt(500),
           length = 0.05, angle = 90, code = 3,
           col = method_colors[method_name], lwd = 1.5)
  }
  
  legend("topright", legend = scenario_labels[scenarios],
         pch = scenario_symbols[scenarios], col = method_colors[method_name],
         lty = 1, lwd = 2, cex = 0.9, bg = "white")
}

# Fourth panel: Summary comparison across methods
plot(NULL, xlim = c(min(plot_data$n), max(plot_data$n)),
     ylim = c(0, max(plot_data$mean_mse_total, na.rm = TRUE) * 1.2),
     xlab = "Sample Size (n)", ylab = "Mean Squared Error",
     main = "All Methods Comparison (S1)",
     log = "x")

s1_data <- plot_data[scenario == "S1"]
for (method_name in c("GCM", "PCM", "WGSC")) {
  method_data <- s1_data[method == method_name]
  method_data <- method_data[order(method_data$n)]
  lines(method_data$n, method_data$mean_mse_total,
        col = method_colors[method_name], lwd = 2, type = "b", pch = 19,
        cex = 1.5)
}

legend("topright", legend = c("GCM", "PCM", "WGSC"),
       col = method_colors[c("GCM", "PCM", "WGSC")],
       lty = 1, lwd = 2, pch = 19, cex = 0.9, bg = "white")

dev.off()
cat("Figure 1 saved to:", file.path(fig_dir, "fig1_mse_vs_sample_size.pdf"),
    "\n\n")

# ------------------------------------------------------------------------------
# Figure 2: Test Statistic Stability (CV) vs Sample Size
# ------------------------------------------------------------------------------

cat("Creating Figure 2: Test Statistic Stability...\n")

pdf(file.path(fig_dir, "fig2_test_stat_stability.pdf"),
    width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.2, cex.axis = 1.1)

# Panel 1: Coefficient of Variation of Test Statistics
plot_data_cv <- results_agg[h_label == "Alternative" & !is.na(cv_test_stat)]

plot(NULL, xlim = c(min(plot_data_cv$n), max(plot_data_cv$n)),
     ylim = c(0, max(plot_data_cv$cv_test_stat, na.rm = TRUE) * 1.2),
     xlab = "Sample Size (n)",
     ylab = "CV of Test Statistic",
     main = "Test Statistic Stability (S1)",
     log = "x")

s1_cv <- plot_data_cv[scenario == "S1"]
for (method_name in c("GCM", "PCM", "WGSC", "RCIT")) {
  method_data <- s1_cv[method == method_name]
  method_data <- method_data[order(method_data$n)]
  lines(method_data$n, method_data$cv_test_stat,
        col = method_colors[method_name], lwd = 2, type = "b", pch = 19,
        cex = 1.5)
}

legend("topright", legend = c("GCM", "PCM", "WGSC", "RCIT"),
       col = method_colors, lty = 1, lwd = 2, pch = 19, cex = 0.9,
       bg = "white")

# Panel 2: SD of Test Statistics
plot(NULL, xlim = c(min(plot_data_cv$n), max(plot_data_cv$n)),
     ylim = c(0, max(plot_data_cv$sd_test_stat, na.rm = TRUE) * 1.2),
     xlab = "Sample Size (n)",
     ylab = "SD of Test Statistic",
     main = "Test Statistic Variability (S1)",
     log = "x")

for (method_name in c("GCM", "PCM", "WGSC", "RCIT")) {
  method_data <- s1_cv[method == method_name]
  method_data <- method_data[order(method_data$n)]
  lines(method_data$n, method_data$sd_test_stat,
        col = method_colors[method_name], lwd = 2, type = "b", pch = 19,
        cex = 1.5)
}

legend("topright", legend = c("GCM", "PCM", "WGSC", "RCIT"),
       col = method_colors, lty = 1, lwd = 2, pch = 19, cex = 0.9,
       bg = "white")

dev.off()
cat("Figure 2 saved to:",
    file.path(fig_dir, "fig2_test_stat_stability.pdf"), "\n\n")

# ------------------------------------------------------------------------------
# Figure 3: Scatter Plots - MSE vs Rejection Rate
# ------------------------------------------------------------------------------

cat("Creating Figure 3: MSE vs Rejection Rate Scatter Plots...\n")

pdf(file.path(fig_dir, "fig3_mse_vs_rejection.pdf"),
    width = 12, height = 10)
par(mfrow = c(3, 3), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.2, cex.axis = 1.1)

# Create scatter plots for each scenario and method
for (sc in c("S1", "S2", "S3")) {
  for (method_name in c("GCM", "PCM", "WGSC")) {
    plot_data_scatter <- results_ind[
      scenario == sc & method == method_name & h_label == "Alternative" &
        !is.na(mse_total)
    ]
    
    if (nrow(plot_data_scatter) > 0) {
      # Aggregate by binning MSE values
      n_bins <- 20
      plot_data_scatter$mse_bin <- cut(plot_data_scatter$mse_total,
                                       breaks = n_bins)
      
      agg_data <- plot_data_scatter %>%
        group_by(mse_bin) %>%
        summarize(
          mean_mse = mean(mse_total, na.rm = TRUE),
          mean_rejection = mean(rejection, na.rm = TRUE),
          n_obs = n(),
          .groups = "drop"
        ) %>%
        filter(n_obs >= 10)  # Only show bins with sufficient data
      
      plot(agg_data$mean_mse, agg_data$mean_rejection,
           xlab = "Mean Squared Error",
           ylab = "Rejection Rate",
           main = sprintf("%s - %s", method_name,
                          c(S1 = "Mean Shift", S2 = "Heteroscedastic",
                            S3 = "Nonlinear")[sc]),
           pch = 19, col = method_colors[method_name], cex = 1.5)
      
      # Add smoothing line
      if (nrow(agg_data) >= 3) {
        loess_fit <- loess(mean_rejection ~ mean_mse, data = agg_data,
                           span = 0.75)
        x_seq <- seq(min(agg_data$mean_mse), max(agg_data$mean_mse),
                     length.out = 100)
        y_pred <- predict(loess_fit, newdata = data.frame(mean_mse = x_seq))
        lines(x_seq, y_pred, col = method_colors[method_name], lwd = 2)
      }
      
      # Add correlation coefficient
      cor_val <- cor(plot_data_scatter$mse_total,
                     plot_data_scatter$rejection,
                     use = "complete.obs")
      text(max(agg_data$mean_mse) * 0.7, max(agg_data$mean_rejection) * 0.9,
           sprintf("r = %.3f", cor_val), cex = 1.1)
    }
  }
}

dev.off()
cat("Figure 3 saved to:", file.path(fig_dir, "fig3_mse_vs_rejection.pdf"),
    "\n\n")

# ------------------------------------------------------------------------------
# Figure 4: Correlation Heatmap
# ------------------------------------------------------------------------------

cat("Creating Figure 4: Correlation Heatmap...\n")

# Prepare data for heatmap
heatmap_data <- correlation_results[h_label == "Alternative" &
                                      !is.na(cor_mse_rejection)]

# Create matrix for heatmap
cor_matrix <- dcast(heatmap_data, scenario + n ~ method,
                    value.var = "cor_mse_rejection")

pdf(file.path(fig_dir, "fig4_correlation_heatmap.pdf"),
    width = 10, height = 8)

# Create heatmap using base R
layout(matrix(c(1, 1, 1, 2), nrow = 1), widths = c(3, 3, 3, 0.5))
par(mar = c(8, 10, 4, 1))

# Prepare matrix data
scenarios <- unique(heatmap_data$scenario)
n_vals <- unique(heatmap_data$n)
methods <- c("GCM", "PCM", "WGSC")

# Create matrix with scenario×n as rows, methods as columns
mat_data <- matrix(NA, nrow = length(scenarios) * length(n_vals),
                   ncol = length(methods))
rownames(mat_data) <- paste0(
  rep(scenarios, each = length(n_vals)), "_n",
  rep(n_vals, length(scenarios))
)
colnames(mat_data) <- methods

for (i in seq_len(nrow(heatmap_data))) {
  row <- heatmap_data[i, ]
  if (row$method %in% methods) {
    row_name <- paste0(row$scenario, "_n", row$n)
    mat_data[row_name, row$method] <- row$cor_mse_rejection
  }
}

# Plot heatmap
image(1:ncol(mat_data), 1:nrow(mat_data), t(mat_data),
      col = colorRampPalette(c("blue", "white", "red"))(100),
      xlab = "", ylab = "", main = "Correlation: MSE vs Rejection Rate",
      axes = FALSE)

# Add axes
axis(1, at = 1:ncol(mat_data), labels = colnames(mat_data), las = 2,
     cex.axis = 1.1)
axis(2, at = 1:nrow(mat_data), labels = rownames(mat_data), las = 2,
     cex.axis = 0.8)

# Add grid
abline(h = seq(0.5, nrow(mat_data) + 0.5, 1), col = "gray", lwd = 0.5)
abline(v = seq(0.5, ncol(mat_data) + 0.5, 1), col = "gray", lwd = 0.5)

# Add text annotations
for (i in 1:nrow(mat_data)) {
  for (j in 1:ncol(mat_data)) {
    if (!is.na(mat_data[i, j])) {
      text(j, i, sprintf("%.2f", mat_data[i, j]), cex = 0.9)
    }
  }
}

# Add color bar
par(mar = c(8, 1, 4, 3))
color_bar <- seq(-1, 1, length.out = 100)
image(1, color_bar, matrix(color_bar, nrow = 1),
      col = colorRampPalette(c("blue", "white", "red"))(100),
      xlab = "", ylab = "", axes = FALSE)
axis(4, at = seq(-1, 1, 0.5), las = 2, cex.axis = 1.0)

dev.off()
cat("Figure 4 saved to:", file.path(fig_dir, "fig4_correlation_heatmap.pdf"),
    "\n\n")

# ------------------------------------------------------------------------------
# Figure 5: Power vs MSE for All Methods
# ------------------------------------------------------------------------------

cat("Creating Figure 5: Power vs MSE Comparison...\n")

pdf(file.path(fig_dir, "fig5_power_vs_mse_comparison.pdf"),
    width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.2, cex.axis = 1.1)

# Panel 1: Rejection rate vs mean MSE (n=1000)
plot_n1000 <- correlation_results[n == 1000 & h_label == "Alternative" &
                                    method != "RCIT"]

plot(plot_n1000$mean_mse, plot_n1000$mean_rejection,
     xlab = "Mean MSE", ylab = "Rejection Rate (Power)",
     main = "Power vs MSE (n=1000)",
     xlim = c(0, max(plot_n1000$mean_mse) * 1.1),
     ylim = c(0, 1),
     pch = 19, col = "white", cex = 1.5)

for (method_name in c("GCM", "PCM", "WGSC")) {
  method_data <- plot_n1000[method == method_name]
  points(method_data$mean_mse, method_data$mean_rejection,
         pch = 19, col = method_colors[method_name], cex = 2)
  text(method_data$mean_mse, method_data$mean_rejection,
       labels = method_data$scenario,
       pos = 3, cex = 0.8, col = method_colors[method_name])
}

legend("bottomleft", legend = c("GCM", "PCM", "WGSC"),
       col = method_colors[c("GCM", "PCM", "WGSC")],
       pch = 19, cex = 1.0, bg = "white")

# Panel 2: CV of test statistic vs CV of MSE
plot_cv <- correlation_results[h_label == "Alternative" & method != "RCIT" &
                                 !is.na(cv_mse)]

plot(plot_cv$cv_mse, plot_cv$sd_rejection,
     xlab = "CV of MSE", ylab = "SD of Rejection",
     main = "Instability: MSE vs Test Decision",
     pch = 19, col = "white", cex = 1.5)

for (method_name in c("GCM", "PCM", "WGSC")) {
  method_data <- plot_cv[method == method_name]
  points(method_data$cv_mse, method_data$sd_rejection,
         pch = 19, col = method_colors[method_name], cex = 1.5)
}

# Add trend line
if (nrow(plot_cv[!is.na(cv_mse) & !is.na(sd_rejection)]) >= 3) {
  fit <- lm(sd_rejection ~ cv_mse, data = plot_cv)
  abline(fit, col = "black", lwd = 2, lty = 2)
  
  r2 <- summary(fit)$r.squared
  text(max(plot_cv$cv_mse, na.rm = TRUE) * 0.5,
       max(plot_cv$sd_rejection, na.rm = TRUE) * 0.9,
       sprintf("R² = %.3f", r2), cex = 1.1)
}

legend("topright", legend = c("GCM", "PCM", "WGSC"),
       col = method_colors[c("GCM", "PCM", "WGSC")],
       pch = 19, cex = 1.0, bg = "white")

dev.off()
cat("Figure 5 saved to:",
    file.path(fig_dir, "fig5_power_vs_mse_comparison.pdf"), "\n\n")

# ------------------------------------------------------------------------------
# Summary Statistics Table
# ------------------------------------------------------------------------------

cat("Creating summary statistics table...\n")

summary_table <- correlation_results[h_label == "Alternative"] %>%
  group_by(method, scenario) %>%
  summarize(
    mean_cor_mse_rejection = mean(cor_mse_rejection, na.rm = TRUE),
    sd_cor_mse_rejection = sd(cor_mse_rejection, na.rm = TRUE),
    mean_mse = mean(mean_mse, na.rm = TRUE),
    mean_power = mean(mean_rejection, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.table()

summary_file <- "results/cit_estimation_summary_table.csv"
fwrite(summary_table, summary_file)
cat("Summary table saved to:", summary_file, "\n\n")

# Print summary to console
cat("\n", strrep("=", 80), "\n")
cat("SUMMARY: Correlation between MSE and Rejection Rate\n")
cat(strrep("=", 80), "\n")
print(summary_table)

# ------------------------------------------------------------------------------
# Completion Message
# ------------------------------------------------------------------------------

cat("\n", strrep("=", 80), "\n")
cat("Visualization and analysis completed successfully!\n")
cat(strrep("=", 80), "\n\n")

cat("Generated figures:\n")
cat("  1. fig1_mse_vs_sample_size.pdf - Estimation errors vs n\n")
cat("  2. fig2_test_stat_stability.pdf - Test statistic stability\n")
cat("  3. fig3_mse_vs_rejection.pdf - MSE vs rejection scatter plots\n")
cat("  4. fig4_correlation_heatmap.pdf - Correlation heatmap\n")
cat("  5. fig5_power_vs_mse_comparison.pdf - Power vs MSE comparison\n\n")

cat("Generated data files:\n")
cat("  1. cit_estimation_correlation_analysis.csv - Correlation metrics\n")
cat("  2. cit_estimation_summary_table.csv - Summary statistics\n\n")

cat("Key findings to report:\n")
cat("  - Negative correlations indicate: Higher MSE → Lower power\n")
cat("  - CV of MSE predicts test decision instability\n")
cat("  - Methods with better regression → more stable tests\n\n")

