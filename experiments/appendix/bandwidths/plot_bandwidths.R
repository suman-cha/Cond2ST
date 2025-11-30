# =============================================================================
# Bandwidth sweep plots for Linear-time MMD (and cross-fitted version)
# - Produces journal-ready figures with no confidence intervals.
# - Assumes CSVs include: scenario, test_type, test_name, n, h_label,
#   rejection_rate, and a bandwidth column (named one of: bandwidth, bw, sigma).
# =============================================================================
rm(list = ls())
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(scales)
    library(purrr)
    library(forcats)
})

# -------------------------- file list (edit as needed) -----------------------
files <- c(
    # "simulation_results_real_high_dim_bandwidth_comparison.csv",
    # "simulation_results_real_low_dim_bandwidth_comparison.csv",
    "simulation_results_S1B_bandwidth_comparison.csv",
    "simulation_results_S2B_bandwidth_comparison.csv",
    "simulation_results_S3B_bandwidth_comparison.csv"
)

# If running from project root and files live in /mnt/data/, prepend that path:
files <- file.path("results", "ablations", "bandwidths", files)  # adjust to your working dir

stopifnot(all(file.exists(files)))

# -------------------------- Data Loading -------------------------------------
read_and_process <- function(filepath) {
    df <- read_csv(filepath, show_col_types = FALSE)
    
    scenario <- tools::file_path_sans_ext(basename(filepath))
    scenario <- str_remove(scenario, "^simulation_results_")
    scenario <- str_remove(scenario, "_bandwidth_comparison$")
    df$scenario <- scenario
    
    df <- df %>%
        mutate(
            test_name = case_when(
                str_detect(test_name, "CV.*Linear") ~ "CV LinearMMD",
                str_detect(test_name, "Linear") ~ "LinearMMD",
                TRUE ~ test_name
            )
        )
    
    return(df)
}

dat <- map_dfr(files, read_and_process)

# Clean and prepare data - KEEP median
dat <- dat %>%
    filter(
        !is.na(rejection_rate),
        test_name %in% c("LinearMMD", "CV LinearMMD")
    ) %>%
    mutate(
        scenario = factor(scenario, levels = c("S1B", "S2B", "S3B")),
        test_name = factor(test_name, levels = c("LinearMMD", "CV LinearMMD")),
        n = as.numeric(n),
        # Handle bandwidth: keep median as string, convert others to numeric
        bandwidth_clean = ifelse(bandwidth == "median", "median", bandwidth),
        bandwidth_numeric = suppressWarnings(as.numeric(bandwidth)),
        bandwidth_label = factor(
            bandwidth_clean,
            levels = c("0.1", "0.5", "1", "median", "5", "10")
        ),
        n_plot = case_when(
            n == 200 ~ 1,
            n == 500 ~ 2,
            n == 1000 ~ 3,
            n == 2000 ~ 4,
            TRUE ~ NA_real_
        )
    )

scenario_labels <- c(
    "S1B" = "Scenario 1B",
    "S2B" = "Scenario 2B",
    "S3B" = "Scenario 3B"
)

dat <- dat %>%
    mutate(scenario_label = factor(scenario_labels[as.character(scenario)],
                                   levels = scenario_labels))

# -------------------------- Plot Theme ---------------------------------------
theme_pub <- theme_minimal(base_size = 11, base_family = "serif") +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.6),
        panel.background = element_rect(fill = "white", color = NA),
        
        legend.position = "bottom",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.4, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.box.spacing = unit(0.4, "cm"),
        legend.background = element_rect(fill = "white", color = NA),
        
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9, color = "grey20"),
        axis.ticks = element_line(color = "grey60", linewidth = 0.4),
        axis.ticks.length = unit(0.15, "cm"),
        
        strip.text = element_text(size = 10, face = "bold", color = "grey10"),
        strip.background = element_rect(fill = "grey90", color = "grey60", linewidth = 0.5),
        
        plot.margin = margin(4, 8, 4, 8)
    )

# 6 colors for 6 bandwidths
bandwidth_colors <- c(
    "0.1" = "#E41A1C",
    "0.5" = "#377EB8",
    "1" = "#4DAF4A",
    "median" = "#FF7F00",
    "5" = "#984EA3",
    "10" = "#A65628"
)

bandwidth_linetypes <- c(
    "0.1" = "solid",
    "0.5" = "dashed",
    "1" = "solid",
    "median" = "dotdash",
    "5" = "longdash",
    "10" = "dotted"
)

bandwidth_shapes <- c(
    "0.1" = 16,
    "0.5" = 17,
    "1" = 15,
    "median" = 18,
    "5" = 16,
    "10" = 17
)

# -------------------------- Plotting Functions -------------------------------
make_bandwidth_plot <- function(data, hyp_label, test_method, show_alpha = TRUE, 
                                show_x_label = TRUE, show_legend = TRUE) {
    df <- data %>% 
        filter(h_label == hyp_label, test_name == test_method)
    
    if (nrow(df) == 0) {
        message("No data for h_label = ", hyp_label, ", test_name = ", test_method)
        return(NULL)
    }
    
    p <- ggplot(df, aes(x = n_plot, y = rejection_rate,
                        color = bandwidth_label, 
                        linetype = bandwidth_label,
                        shape = bandwidth_label,
                        group = bandwidth_label)) +
        geom_line(linewidth = 0.6, alpha = 0.9) +
        geom_point(size = 1.8, alpha = 0.95, stroke = 0.5) +
        scale_x_continuous(
            breaks = c(1, 2, 3, 4),
            labels = c("200", "500", "1000", "2000")
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2),
            labels = number_format(accuracy = 0.01),
            expand = expansion(mult = c(0.02, 0.02))
        ) +
        scale_color_manual(
            name = "Bandwidth",
            values = bandwidth_colors
        ) +
        scale_linetype_manual(
            name = "Bandwidth",
            values = bandwidth_linetypes
        ) +
        scale_shape_manual(
            name = "Bandwidth",
            values = bandwidth_shapes
        ) +
        labs(
            x = if (show_x_label) "Sample Size" else "",
            y = if (hyp_label == "Null") "Type I Error Rate" else "Power"
        ) +
        facet_wrap(~ scenario_label, nrow = 1, ncol = 3) +
        theme_pub
    
    if (!show_x_label) {
        p <- p + theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()
        )
    }
    
    if (!show_legend) {
        p <- p + theme(legend.position = "none")
    } else {
        p <- p + guides(
            color = guide_legend(
                nrow = 1,
                override.aes = list(size = 3, linewidth = 1)
            ),
            linetype = guide_legend(
                nrow = 1,
                override.aes = list(size = 3, linewidth = 1)
            ),
            shape = guide_legend(
                nrow = 1,
                override.aes = list(size = 3, linewidth = 1)
            )
        )
    }
    
    if (show_alpha && hyp_label == "Null") {
        p <- p + geom_hline(
            yintercept = 0.05, 
            color = "#D62728",
            linetype = "dotted", 
            linewidth = 0.35,
            alpha = 0.8
        )
    }
    
    return(p)
}

# -------------------------- Generate Plots -----------------------------------
p_null_linear <- make_bandwidth_plot(dat, "Null", "LinearMMD", 
                                     show_alpha = TRUE, 
                                     show_x_label = FALSE,
                                     show_legend = FALSE)

p_alt_linear <- make_bandwidth_plot(dat, "Alternative", "LinearMMD", 
                                    show_alpha = FALSE, 
                                    show_x_label = TRUE,
                                    show_legend = TRUE)

p_null_cv <- make_bandwidth_plot(dat, "Null", "CV LinearMMD", 
                                 show_alpha = TRUE, 
                                 show_x_label = FALSE,
                                 show_legend = FALSE)

p_alt_cv <- make_bandwidth_plot(dat, "Alternative", "CV LinearMMD", 
                                show_alpha = FALSE, 
                                show_x_label = TRUE,
                                show_legend = TRUE)

# -------------------------- Save Combined Plots ------------------------------
out_dir <- "Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_combined_plot <- function(p_top, p_bottom, filename, width = 12, height = 7) {
    if (is.null(p_top) || is.null(p_bottom)) {
        warning("One or both plots are NULL, skipping: ", filename)
        return(invisible(NULL))
    }
    
    filepath <- file.path(out_dir, filename)
    
    pdf(file = filepath, width = width, height = height, 
        family = "serif", useDingbats = FALSE)
    grid.arrange(p_top, p_bottom, nrow = 2, heights = c(1, 1.12))
    dev.off()
    
    message("✓ Saved: ", filepath)
    
    png_path <- str_replace(filepath, "\\.pdf$", ".png")
    png(filename = png_path, width = width * 120, height = height * 120, 
        res = 300, type = "cairo")
    grid.arrange(p_top, p_bottom, nrow = 2, heights = c(1, 1.12))
    dev.off()
    
    message("✓ Saved: ", png_path, " (preview)")
}

save_combined_plot(p_null_linear, p_alt_linear, 
                   "bandwidth_comparison_LinearMMD.pdf", 
                   width = 12, height = 7)

save_combined_plot(p_null_cv, p_alt_cv, 
                   "bandwidth_comparison_CV_LinearMMD.pdf", 
                   width = 12, height = 7)

message("\n========================================")
message("Combined figures saved!")
message("Layout: Type I Error (top) + Power (bottom)")
message("Bandwidths: 0.1, 0.5, 1, median, 5, 10")
message("Files:")
message("  - bandwidth_comparison_LinearMMD.pdf")
message("  - bandwidth_comparison_CV_LinearMMD.pdf")
message("Output directory: ", normalizePath(out_dir))
message("========================================")