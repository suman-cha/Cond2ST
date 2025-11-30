# # =============================================================================
# # Computational Cost Visualization for Top-Tier Journal
# # =============================================================================
# 
# rm(list = ls())
# suppressPackageStartupMessages({
#     library(tidyverse)
#     library(scales)
#     library(ggrepel)
#     library(gridExtra)
#     library(kableExtra)
# })
# 
# # =============================================================================
# # CONFIGURATION
# # =============================================================================
# 
# CONFIG <- list(
#     files = list(
#         S1U = "~/WORKSPACE/Cond2st/results/ablations/compute_costs/simulation_results_S1U_computationl_costs.csv",
#         real_low = "~/WORKSPACE/Cond2st/results/ablations/compute_costs/simulation_results_real_low_dim_computational_costs.csv",
#         real_high = "~/WORKSPACE/Cond2st/results/ablations/compute_costs/simulation_results_real_high_dim_computational_costs.csv",
#         cit_real_high = "~/WORKSPACE/Cond2st/results/ablations/compute_costs/CIT_simulation_results_real_high_dim_computational_costs.csv"
#     ),
#     output_dir = "Figures/computational_costs",
#     plot_width = 12,
#     plot_height = 7,
#     dpi = 300
# )
# 
# # Create output directory
# dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)
# 
# # =============================================================================
# # DATA LOADING
# # =============================================================================
# 
# data_s1u <- read_csv(CONFIG$files$S1U, show_col_types = FALSE) %>%
#     mutate(scenario = "Scenario 1U")
# 
# data_real_low <- read_csv(CONFIG$files$real_low, show_col_types = FALSE) %>%
#     filter(estimator == "LL") %>%
#     mutate(scenario = "Real Low Dim")
# 
# data_real_high <- read_csv(CONFIG$files$real_high, show_col_types = FALSE) %>%
#     filter(estimator == "LL") %>%
#     mutate(scenario = "Real High Dim")
# 
# data_cit_real_high <- read_csv(CONFIG$files$cit_real_high, show_col_types = FALSE) %>%
#     filter(estimator == "lm") %>%
#     mutate(scenario = "Real High Dim")
# 
# dat <- bind_rows(
#     data_s1u,
#     data_real_low,
#     data_real_high,
#     data_cit_real_high
# ) %>%
#     mutate(
#         test_name_clean = str_remove(test_name, "_test$"),
#         test_name_clean = str_replace_all(test_name_clean, "_", "-"),
#         scenario = factor(scenario, levels = c("Scenario 1U", "Real Low Dim", "Real High Dim"))
#     )
# 
# # =============================================================================
# # COLOR PALETTE
# # =============================================================================
# 
# method_colors <- c(
#     "LinearMMD" = "#E41A1C",
#     "CV-LinearMMD" = "#377EB8",
#     "CLF" = "#4DAF4A",
#     "CV-CLF" = "#984EA3",
#     "CP" = "#FF7F00",
#     "debiased" = "#A65628",
#     "BlockMMD" = "#F781BF",
#     "bootstrap-MMD" = "#999999",
#     "RCIT" = "#1B9E77",
#     "GCM" = "#D95F02",
#     "PCM" = "#7570B3",
#     "WGSC" = "#E7298A"
# )
# 
# method_linetypes <- c(
#     "LinearMMD" = "solid",
#     "CV-LinearMMD" = "solid",
#     "CLF" = "solid",
#     "CV-CLF" = "solid",
#     "CP" = "solid",
#     "debiased" = "solid",
#     "BlockMMD" = "solid",
#     "bootstrap-MMD" = "solid",
#     "RCIT" = "dashed",
#     "GCM" = "dashed",
#     "PCM" = "dashed",
#     "WGSC" = "dashed"
# )
# 
# # =============================================================================
# # THEME
# # =============================================================================
# 
# theme_pub <- theme_minimal(base_size = 11, base_family = "serif") +
#     theme(
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
#         panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
#         panel.background = element_rect(fill = "white", color = NA),
#         
#         legend.position = "bottom",
#         legend.title = element_text(size = 10, face = "bold"),
#         legend.text = element_text(size = 8),
#         legend.key.width = unit(1.0, "cm"),
#         legend.key.height = unit(0.4, "cm"),
#         
#         axis.title = element_text(size = 11, face = "bold"),
#         axis.text = element_text(size = 9, color = "grey20"),
#         axis.ticks = element_line(color = "grey60", linewidth = 0.4),
#         
#         strip.text = element_text(size = 10, face = "bold", color = "grey10"),
#         strip.background = element_rect(fill = "grey90", color = "grey60", linewidth = 0.5),
#         
#         plot.margin = margin(10, 10, 10, 10)
#     )
# 
# # =============================================================================
# # FIGURE 1: COMPUTATIONAL TIME (All scenarios, Null)
# # =============================================================================
# 
# dat_null <- dat %>% filter(h_label == "Null")
# 
# p1 <- dat_null %>%
#     ggplot(aes(x = factor(n), y = avg_time_sec, 
#                color = test_name_clean, 
#                linetype = test_name_clean,
#                group = test_name_clean)) +
#     geom_line(aes(group = test_name_clean), linewidth = 0.7, alpha = 0.9) +
#     geom_point(size = 2, alpha = 0.95) +
#     scale_y_log10(
#         breaks = c(0.001, 0.01, 0.1, 1, 10),
#         labels = c("0.001", "0.01", "0.1", "1", "10")
#     ) +
#     scale_color_manual(
#         name = "Method",
#         values = method_colors
#     ) +
#     scale_linetype_manual(
#         name = "Method",
#         values = method_linetypes
#     ) +
#     facet_wrap(~ scenario, nrow = 1, scales = "free_x") +
#     labs(
#         x = "Sample Size per Group",
#         y = "Average Computation Time (seconds, log scale)"
#     ) +
#     theme_pub +
#     guides(
#         color = guide_legend(nrow = 2, override.aes = list(linewidth = 1.2, size = 2.5)),
#         linetype = guide_legend(nrow = 2, override.aes = list(linewidth = 1.2))
#     )
# 
# # =============================================================================
# # FIGURE 2A: TRADE-OFF - Cost vs Power (S1U only, has Alternative)
# # =============================================================================
# 
# dat_power_s1u <- dat %>%
#     filter(scenario == "Scenario 1U", h_label == "Alternative") %>%
#     group_by(test_name_clean) %>%
#     summarise(
#         avg_time = mean(avg_time_sec, na.rm = TRUE),
#         power = mean(rejection_rate, na.rm = TRUE),
#         .groups = "drop"
#     ) %>%
#     mutate(scenario = "Scenario 1U")
# 
# p2a <- ggplot(dat_power_s1u, 
#               aes(x = avg_time, y = power, color = test_name_clean)) +
#     geom_point(size = 4, alpha = 0.9) +
#     geom_text_repel(
#         aes(label = test_name_clean), 
#         size = 3, 
#         family = "serif",
#         max.overlaps = 20, 
#         box.padding = 0.5,
#         segment.color = "grey50",
#         segment.size = 0.3
#     ) +
#     scale_x_log10(
#         breaks = c(0.01, 0.1, 1, 10),
#         labels = c("0.01", "0.1", "1", "10")
#     ) +
#     scale_y_continuous(
#         limits = c(0, 1),
#         breaks = seq(0, 1, by = 0.2),
#         labels = number_format(accuracy = 0.01)
#     ) +
#     scale_color_manual(values = method_colors) +
#     labs(
#         x = "Average Computation Time (seconds, log scale)",
#         y = "Power",
#         title = "Scenario 1U: Computational Cost vs Power"
#     ) +
#     theme_pub +
#     theme(legend.position = "none")
# 
# # =============================================================================
# # FIGURE 2B: TRADE-OFF - Cost vs Type I Error (All scenarios, Null)
# # Use largest sample size for comparison
# # =============================================================================
# 
# # Get largest n for each scenario
# dat_tradeoff_null <- dat_null %>%
#     group_by(scenario) %>%
#     filter(n == max(n)) %>%
#     ungroup() %>%
#     group_by(scenario, test_name_clean) %>%
#     summarise(
#         avg_time = mean(avg_time_sec, na.rm = TRUE),
#         type1_error = mean(rejection_rate, na.rm = TRUE),
#         n_size = first(n),
#         .groups = "drop"
#     )
# 
# p2b <- ggplot(dat_tradeoff_null, 
#               aes(x = avg_time, y = type1_error, 
#                   color = test_name_clean)) +
#     geom_hline(yintercept = 0.05, linetype = "dashed", 
#                color = "red", linewidth = 0.4, alpha = 0.7) +
#     geom_point(size = 3.5, alpha = 0.9) +
#     geom_text_repel(
#         aes(label = test_name_clean), 
#         size = 2.5, 
#         family = "serif",
#         max.overlaps = 20, 
#         box.padding = 0.3,
#         segment.color = "grey50",
#         segment.size = 0.2
#     ) +
#     scale_x_log10(
#         breaks = c(0.01, 0.1, 1, 10),
#         labels = c("0.01", "0.1", "1", "10")
#     ) +
#     scale_y_continuous(
#         limits = c(0, 0.5),
#         breaks = seq(0, 0.5, by = 0.1),
#         labels = number_format(accuracy = 0.01)
#     ) +
#     scale_color_manual(values = method_colors) +
#     facet_wrap(~ scenario, nrow = 1) +
#     labs(
#         x = "Average Computation Time (seconds, log scale)",
#         y = "Type I Error Rate"
#     ) +
#     theme_pub +
#     theme(legend.position = "none")
# 
# # =============================================================================
# # FIGURE 3: COMBINED TRADE-OFF (Power for S1U, Type I Error for Real Data)
# # =============================================================================
# 
# # Combine S1U power and Real data type I error
# dat_tradeoff_combined <- bind_rows(
#     dat_power_s1u %>% 
#         mutate(metric = "Power", value = power) %>%
#         select(scenario, test_name_clean, avg_time, metric, value),
#     
#     dat_tradeoff_null %>% 
#         filter(scenario %in% c("Real Low Dim", "Real High Dim")) %>%
#         mutate(metric = "Type I Error", value = type1_error) %>%
#         select(scenario, test_name_clean, avg_time, metric, value)
# )
# 
# p3 <- ggplot(dat_tradeoff_combined, 
#              aes(x = avg_time, y = value, color = test_name_clean)) +
#     geom_hline(data = filter(dat_tradeoff_combined, metric == "Type I Error"),
#                yintercept = 0.05, linetype = "dashed", 
#                color = "red", linewidth = 0.4, alpha = 0.7) +
#     geom_point(size = 3.5, alpha = 0.9) +
#     geom_text_repel(
#         aes(label = test_name_clean), 
#         size = 2.5, 
#         family = "serif",
#         max.overlaps = 20, 
#         box.padding = 0.3,
#         segment.color = "grey50",
#         segment.size = 0.2
#     ) +
#     scale_x_log10(
#         breaks = c(0.01, 0.1, 1, 10),
#         labels = c("0.01", "0.1", "1", "10")
#     ) +
#     scale_y_continuous(
#         limits = c(0, 1),
#         breaks = seq(0, 1, by = 0.2),
#         labels = number_format(accuracy = 0.01)
#     ) +
#     scale_color_manual(values = method_colors) +
#     facet_grid(metric ~ scenario, scales = "free_y") +
#     labs(
#         x = "Average Computation Time (seconds, log scale)",
#         y = "Performance Metric"
#     ) +
#     theme_pub +
#     theme(legend.position = "none")
# 
# # =============================================================================
# # SAVE PLOTS
# # =============================================================================
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_time_all_methods.pdf"),
#     p1, width = CONFIG$plot_width, height = CONFIG$plot_height
# )
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_time_all_methods.png"),
#     p1, width = CONFIG$plot_width, height = CONFIG$plot_height, dpi = CONFIG$dpi
# )
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_cost_vs_power_S1U.pdf"),
#     p2a, width = 7, height = 6
# )
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_cost_vs_type1error.pdf"),
#     p2b, width = CONFIG$plot_width, height = 5
# )
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_cost_vs_type1error.png"),
#     p2b, width = CONFIG$plot_width, height = 5, dpi = CONFIG$dpi
# )
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_tradeoff_combined.pdf"),
#     p3, width = CONFIG$plot_width, height = 7
# )
# 
# ggsave(
#     file.path(CONFIG$output_dir, "computational_tradeoff_combined.png"),
#     p3, width = CONFIG$plot_width, height = 7, dpi = CONFIG$dpi
# )
# 
# # =============================================================================
# # SUMMARY STATISTICS
# # =============================================================================
# 
# cat("\n=================================================================\n")
# cat("         COMPUTATIONAL COST & PERFORMANCE SUMMARY\n")
# cat("=================================================================\n\n")
# 
# cat("POWER (Scenario 1U only):\n")
# cat(strrep("-", 70), "\n")
# print(dat_power_s1u %>% arrange(desc(power)), n = Inf)
# 
# cat("\n\nTYPE I ERROR CONTROL (All scenarios, largest n):\n")
# cat(strrep("-", 70), "\n")
# print(dat_tradeoff_null %>% arrange(scenario, avg_time), n = Inf)
# 
# cat("\n=================================================================\n")
# cat("All plots saved to:", CONFIG$output_dir, "\n")
# cat("Files:\n")
# cat("  - computational_time_all_methods.pdf (main figure)\n")
# cat("  - computational_cost_vs_power_S1U.pdf (S1U power trade-off)\n")
# cat("  - computational_cost_vs_type1error.pdf (type I error trade-off)\n")
# cat("  - computational_tradeoff_combined.pdf (combined view)\n")
# cat("=================================================================\n")

# =============================================================================
# Computational Cost vs Performance - Scenario 1U (CORRECTED)
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
    library(tidyverse)
    library(scales)
    library(ggrepel)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG <- list(
    file = "~/WORKSPACE/Cond2st/results/ablations/compute_costs/simulation_results_S1U_computationl_costs.csv",
    output_dir = "~/WORKSPACE/Cond2st/results/ablations/compute_costs/figures",
    plot_width = 12,
    plot_height = 10,
    dpi = 300
)

CONFIG$file <- path.expand(CONFIG$file)
CONFIG$output_dir <- path.expand(CONFIG$output_dir)
dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# DATA LOADING
# =============================================================================

if (!file.exists(CONFIG$file)) {
    stop("File not found: ", CONFIG$file)
}

dat <- read_csv(CONFIG$file, show_col_types = FALSE) %>%
    mutate(
        test_name_clean = str_remove(test_name, "_test$"),
        test_name_clean = str_replace_all(test_name_clean, "_", "-")
    )

cat("\n=================================================================\n")
cat("Data loaded:", nrow(dat), "observations\n")
cat("Methods:", length(unique(dat$test_name_clean)), "\n")
cat("Sample sizes:", paste(sort(unique(dat$n)), collapse = ", "), "\n")
cat("=================================================================\n\n")

# =============================================================================
# COLOR PALETTE
# =============================================================================

method_colors <- c(
    "LinearMMD" = "#E41A1C",
    "CV-LinearMMD" = "#377EB8",
    "CLF" = "#4DAF4A",
    "CV-CLF" = "#984EA3",
    "CP" = "#FF7F00",
    "debiased" = "#A65628",
    "BlockMMD" = "#F781BF",
    "CV-BlockMMD" = "#66C2A5",
    "bootstrap-MMD" = "#999999",
    "RCIT" = "#1B9E77",
    "GCM" = "#D95F02",
    "PCM" = "#7570B3",
    "WGSC" = "#E7298A"
)

# =============================================================================
# THEME
# =============================================================================

theme_pub <- theme_minimal(base_size = 11, base_family = "serif") +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
        panel.background = element_rect(fill = "white", color = NA),
        
        legend.position = "none",
        
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9, color = "grey20"),
        axis.ticks = element_line(color = "grey60", linewidth = 0.4),
        
        strip.text = element_text(size = 10, face = "bold", color = "grey10"),
        strip.background = element_rect(fill = "grey90", color = "grey60", linewidth = 0.5),
        
        plot.margin = margin(8, 8, 8, 8)
    )

# =============================================================================
# FIGURE 1: 4x2 GRID
# TOP ROW: Type I Error (n=200, 500, 1000, 2000)
# BOTTOM ROW: Power (n=200, 500, 1000, 2000)
# =============================================================================

dat_by_n <- dat %>%
    mutate(
        metric = ifelse(h_label == "Null", "Type I Error", "Power"),
        # Type I Error가 위에 오도록 factor level 설정
        metric = factor(metric, levels = c("Type I Error", "Power")),
        # Sample size도 순서대로 factor 설정
        n_label = factor(paste0("n = ", n), 
                         levels = c("n = 200", "n = 500", "n = 1000", "n = 2000"))
    ) %>%
    select(test_name_clean, n, n_label, metric, avg_time_sec, rejection_rate)

# Type I Error 패널용 reference line 데이터
dat_ref_line <- dat_by_n %>%
    filter(metric == "Type I Error") %>%
    distinct(metric, n_label)

p_figure1 <- dat_by_n %>%
    ggplot(aes(x = avg_time_sec, y = rejection_rate, color = test_name_clean)) +
    # Type I Error 패널에만 reference line
    geom_hline(
        data = dat_ref_line,
        yintercept = 0.05,
        linetype = "dashed",
        color = "red",
        linewidth = 0.4,
        alpha = 0.7
    ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
        aes(label = test_name_clean),
        size = 2.2,
        family = "serif",
        max.overlaps = 30,
        box.padding = 0.25,
        point.padding = 0.2,
        segment.color = "grey50",
        segment.size = 0.15,
        min.segment.length = 0
    ) +
    scale_x_log10(
        breaks = c(0.01, 0.1, 1, 10),
        labels = c("0.01", "0.1", "1", "10")
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = number_format(accuracy = 0.01)
    ) +
    scale_color_manual(values = method_colors) +
    facet_grid(metric ~ n_label) +
    labs(
        x = "Average Computation Time (seconds, log scale)",
        y = ""
    ) +
    theme_pub +
    theme(
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        panel.spacing = unit(0.5, "lines")
    )

# =============================================================================
# FIGURE 2: n=2000 ONLY
# LEFT: Type I Error (with red line)
# RIGHT: Power (NO red line)
# =============================================================================

dat_n2000 <- dat %>%
    filter(n == 2000) %>%
    mutate(
        metric = ifelse(h_label == "Null", "Type I Error", "Power"),
        metric = factor(metric, levels = c("Type I Error", "Power"))
    )

# Type I Error만 선택 (Power 제외!)
dat_n2000_type1_only <- dat_n2000 %>%
    filter(metric == "Type I Error")

p_figure2 <- dat_n2000 %>%
    ggplot(aes(x = avg_time_sec, y = rejection_rate, color = test_name_clean)) +
    # Type I Error 데이터만 사용 (Power 패널에는 line이 나타나지 않음)
    geom_hline(
        data = dat_n2000_type1_only,
        aes(yintercept = 0.05),
        linetype = "dashed",
        color = "red",
        linewidth = 0.5,
        alpha = 0.8,
        inherit.aes = FALSE  # 중요: aes 상속 차단
    ) +
    geom_point(size = 4.5, alpha = 0.9) +
    geom_text_repel(
        aes(label = test_name_clean),
        size = 3,
        family = "serif",
        max.overlaps = 30,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.3,
        min.segment.length = 0
    ) +
    scale_x_log10(
        breaks = c(0.01, 0.1, 1, 10),
        labels = c("0.01", "0.1", "1", "10")
    ) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = number_format(accuracy = 0.01)
    ) +
    scale_color_manual(values = method_colors) +
    facet_wrap(~ metric, ncol = 2) +
    labs(
        x = "Average Computation Time (seconds, log scale)",
        y = ""
    ) +
    theme_pub +
    theme(
        strip.text = element_text(size = 11, face = "bold")
    )

# =============================================================================
# SAVE FIGURES
# =============================================================================

ggsave(
    file.path(CONFIG$output_dir, "Figure1_cost_by_samplesize.pdf"),
    p_figure1,
    width = 14,
    height = 10
)

ggsave(
    file.path(CONFIG$output_dir, "Figure1_cost_by_samplesize.png"),
    p_figure1,
    width = 14,
    height = 10,
    dpi = CONFIG$dpi
)

ggsave(
    file.path(CONFIG$output_dir, "Figure2_cost_n2000.pdf"),
    p_figure2,
    width = 12,
    height = 6
)

ggsave(
    file.path(CONFIG$output_dir, "Figure2_cost_n2000.png"),
    p_figure2,
    width = 12,
    height = 6,
    dpi = CONFIG$dpi
)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n=================================================================\n")
cat("Figures saved successfully!\n")
cat("=================================================================\n\n")

cat("Figure 1 Layout:\n")
cat("  TOP ROW: Type I Error (n=200, 500, 1000, 2000)\n")
cat("  BOTTOM ROW: Power (n=200, 500, 1000, 2000)\n")
cat("  Red line (y=0.05) ONLY in Type I Error panels\n\n")

cat("Figure 2 Layout:\n")
cat("  LEFT: Type I Error (with red line at y=0.05)\n")
cat("  RIGHT: Power (NO red line)\n")
cat("  Both panels: y-axis [0, 1]\n\n")

cat("Output directory:", CONFIG$output_dir, "\n")
cat("=================================================================\n")

# Quick verification
cat("\nVerification:\n")
cat("Sample size order in data:\n")
print(levels(dat_by_n$n_label))
cat("\nMetric order in data:\n")
print(levels(dat_by_n$metric))
