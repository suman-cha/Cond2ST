#=============================================================================
# Main Plotting Script for Scenarios 1-3
# Publication-quality heatmaps for top theoretical statistics journals
#=============================================================================
rm(list = ls())

suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(latex2exp)
    library(stringr)
    library(grid)
    library(gridExtra)
    library(scales)
})

# Source plotting utilities
source("./experiments/plot_utils.R")

#------------------------------------------------------------------------------
# Configuration for Scenarios 1-3
#------------------------------------------------------------------------------

# Methods for Scenarios 1-3 (older naming convention)
c2st_methods <- c("debiased_test", "CP_test", "CVCLF_test", "CLF_test",
                  "CVLinearMMD_test", "LinearMMD_test")

cit_methods <- c("WGSC_test", "GCM_test", "PCM_test", "RCoT_test", "RCIT_test")

method_labels <- c(
    "CP_test"           = "CP",
    "debiased_test"     = "DCP",
    "CVLinearMMD_test"  = "$^{\\dagger}$MMD-$\\ell$",
    "LinearMMD_test"    = "MMD-$\\ell$",
    "CVCLF_test"        = "$^{\\dagger}$CLF",
    "CLF_test"          = "CLF",
    "GCM_test"          = "GCM",
    "PCM_test"          = "PCM",
    "WGSC_test"         = "WGSC",
    "RCIT_test"         = "RCIT",
    "RCoT_test"         = "RCoT"
)

#------------------------------------------------------------------------------
# Data Loading and Preparation
#------------------------------------------------------------------------------

load_and_prepare_data <- function(file_path, scenario) {
    data <- fread(file_path)

    scenario_name <- paste0(
        "Scenario ", str_extract(scenario, "\\d"),
        ifelse(str_detect(scenario, "U"), "(U)", "(B)")
    )

    data$Scenario <- scenario_name
    data$test_name <- factor(data$test_name, levels = unique(data$test_name))
    data$extra_param <- factor(data$extra_param, levels = unique(data$extra_param))
    data$h_label <- factor(data$h_label, levels = c("Null", "Alternative"))

    return(data)
}

prepare_data <- function(data) {
    data_c2st <- data %>%
        filter(test_type == "C2ST" & extra_param == "LL" & test_name %in% c2st_methods)
    data_cit <- data %>%
        filter(test_type == "CIT" & extra_param == "TRUE" & test_name %in% cit_methods)

    data_filtered <- rbind(data_c2st, data_cit)
    data_filtered$Method <- data_filtered$test_name
    data_filtered$Method <- factor(data_filtered$Method, levels = c(c2st_methods, cit_methods))

    # Rename test_type for display
    data_filtered$test_type <- ifelse(data_filtered$test_type == "C2ST", "DRT", "CIT")

    return(data_filtered)
}

#------------------------------------------------------------------------------
# Publication-Quality Heatmap Generator
#------------------------------------------------------------------------------

generate_plot_combined <- function(data, scenario_name) {
    # Create Scenario_Hypothesis interaction
    data$Scenario_Hypothesis <- with(data, interaction(
        Scenario,
        factor(h_label, levels = c("Null", "Alternative")),
        sep = " | ", lex.order = TRUE
    ))

    data$Scenario_Hypothesis <- factor(data$Scenario_Hypothesis,
        levels = c(
            paste0("Scenario ", scenario_name, "(U) | Null"),
            paste0("Scenario ", scenario_name, "(U) | Alternative"),
            paste0("Scenario ", scenario_name, "(B) | Null"),
            paste0("Scenario ", scenario_name, "(B) | Alternative")
        )
    )

    data$Method <- factor(data$Method, levels = c(c2st_methods, cit_methods))
    data$test_type <- factor(data$test_type, levels = c("DRT", "CIT"))

    # Professional heatmap
    p <- ggplot(data, aes(x = factor(n), y = Method, fill = mean_result)) +
        geom_tile(color = "white", linewidth = 0.4) +
        scale_fill_gradientn(
            colors = c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404"),
            values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
            limits = c(0, 1),
            name = "Rejection\nRate",
            breaks = c(0, 0.25, 0.5, 0.75, 1),
            labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
        ) +
        geom_text(
            aes(label = sprintf("%.2f", mean_result)),
            color = "black",
            size = 3.5,
            family = "serif"
        ) +
        facet_grid(test_type ~ Scenario_Hypothesis, scales = "free_y", space = "free") +
        labs(
            title = paste("Scenario", scenario_name),
            x = "Sample Size (n)",
            y = "Method"
        ) +
        scale_y_discrete(labels = function(x) TeX(method_labels[x])) +
        theme_minimal(base_size = 11, base_family = "serif") +
        theme(
            # Title
            plot.title = element_text(
                size = 16,
                face = "bold",
                hjust = 0.5,
                margin = margin(b = 15)
            ),

            # Axis
            axis.title = element_text(size = 12, face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10)),
            axis.title.y = element_text(margin = margin(r = 5)),
            axis.text.x = element_text(size = 10, margin = margin(t = 5)),
            axis.text.y = element_text(size = 10, margin = margin(r = 5)),
            axis.ticks = element_blank(),

            # Panel
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
            panel.spacing = unit(0.5, "lines"),

            # Strip
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),

            # Legend
            legend.position = "right",
            legend.title = element_text(size = 10, face = "bold", margin = margin(b = 8)),
            legend.text = element_text(size = 9),
            legend.key.height = unit(2, "cm"),
            legend.key.width = unit(0.4, "cm"),

            # Margins
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
        )

    return(p)
}

#------------------------------------------------------------------------------
# Main: Generate Plots for Scenarios 1-3
#------------------------------------------------------------------------------

save_dir <- "Figures"
if (!dir.exists(save_dir)) {
    dir.create(save_dir)
}

for (scenario_number in 1:3) {
    file_path_U <- paste0("results/simulation_results_S", scenario_number, "U.csv")
    file_path_B <- paste0("results/simulation_results_S", scenario_number, "B.csv")

    # Check if files exist
    if (!file.exists(file_path_U) || !file.exists(file_path_B)) {
        cat("Skipping Scenario", scenario_number, "- result files not found\n")
        next
    }

    data_B <- load_and_prepare_data(file_path_B, paste0("S", scenario_number, "B"))
    data_U <- load_and_prepare_data(file_path_U, paste0("S", scenario_number, "U"))
    data_combined <- rbind(prepare_data(data_B), prepare_data(data_U))

    data_combined$Scenario <- factor(data_combined$Scenario,
        levels = c(
            paste0("Scenario ", scenario_number, "(U)"),
            paste0("Scenario ", scenario_number, "(B)")
        )
    )

    plot <- generate_plot_combined(data_combined, scenario_number)

    # Save as publication-quality PDF
    filename <- file.path(save_dir, paste0("Scenario", scenario_number, ".pdf"))
    save_publication_pdf(plot, filename, width = 10, height = 6)
}

cat("\nAll scenario plots generated successfully.\n")
