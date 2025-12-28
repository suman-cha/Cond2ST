rm(list=ls())
library(reshape2)
library(ggplot2)
library(colorspace)
library(extrafont)
library(gridExtra)
library(dplyr)
library(latex2exp)
library(stringr)

# Function to load and prepare data from a CSV file
load_and_prepare_data <- function(file_path, scenario) {
  data <- read.csv(file_path)
  
  scenario_name <- paste0("Scenario ", str_extract(scenario, "\\d"), 
                          ifelse(str_detect(scenario, "U"), "(U)", "(B)"))
  
  data$Scenario <- scenario_name
  data$test_name <- factor(data$test_name, levels = unique(data$test_name))
  data$h_label <- factor(data$h_label, levels = c("Null", "Alternative"))
  
  return(data)
}

c2st_methods <- c("debiased_test", "CP_test", "CV_CLF_test", "CLF_test", "bootstrap_MMD_test", "CV_BlockMMD_test", "BlockMMD_test", "CV_LinearMMD_test", "LinearMMD_test")
cit_methods <- c("WGSC_test", "PCM_test", "GCM_test", "RCIT_test")
method_labels <- c(
  "LinearMMD_test" = "MMD-\u2113",
  "CV_LinearMMD_test" = "$^\\dagger$MMD-\u2113",
  "BlockMMD_test" = "MMD-b",
  "CV_BlockMMD_test" = "$^\\dagger$MMD-b",
  "bootstrap_MMD_test" = "MMD-wb",
  "CLF_test" = "CLF",
  "CV_CLF_test" = "$^\\dagger$CLF",
  "CP_test" = "CP",
  "debiased_test" = "DCP",
  "RCIT_test" = "RCIT",
  "GCM_test" = "GCM",
  "PCM_test" = "PCM",
  "WGSC_test" = "WGSC"
)

# Function to filter and prepare data
prepare_data <- function(data) {
  data_drt <- data %>%
    filter(test_type == "DRT" & test_name %in% c2st_methods)
  data_cit <- data %>%
    filter(test_type == "CIT" & test_name %in% cit_methods)
  
  data_filtered <- rbind(data_drt, data_cit)
  data_filtered$Method <- data_filtered$test_name
  data_filtered$Method <- factor(data_filtered$Method, levels = c(c2st_methods, cit_methods))
  data_filtered$test_type <- factor(data_filtered$test_type, levels = c("DRT", "CIT"))
  
  return(data_filtered)
}

generate_plot_combined <- function(data, scenario_name) {
  data$Scenario_Hypothesis <- with(data, interaction(Scenario, 
                                                     factor(h_label, levels = c("Null", "Alternative")), 
                                                     sep = " | ", lex.order = TRUE))
  
  data$Scenario_Hypothesis <- factor(data$Scenario_Hypothesis, 
                                     levels = c(paste0("Scenario ", scenario_name, "(U) | Null"),
                                                paste0("Scenario ", scenario_name, "(U) | Alternative"),
                                                paste0("Scenario ", scenario_name, "(B) | Null"),
                                                paste0("Scenario ", scenario_name, "(B) | Alternative")))
  
  data$Method <- factor(data$Method, levels = c(c2st_methods, cit_methods))
  
  ggplot(data, aes(x = factor(n), y = Method, fill = rejection_rate)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradientn(colors = c("lightyellow", "orange", "red"), 
                         values = scales::rescale(c(0, 0.5, 1)),
                         limits = c(0, 1),
                         name = "Rejection \n     rate") +
    geom_text(aes(label = sprintf("%.3f", rejection_rate)), color = "black", size = 12.6, family = "serif") +  
    facet_grid(test_type ~ Scenario_Hypothesis, scales = "free_y", space = "free") +
    labs(
      title = paste("Rejection Rate for Scenario", scenario_name),
      x = "Sample Size (n)", 
      y = "Method"
    ) +
    theme_minimal(base_size = 19.2, base_family = "serif") +  
    scale_y_discrete(labels = function(x) TeX(method_labels[x])) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 36, margin = margin(t = 18)), 
      axis.text.y = element_text(size = 36, margin = margin(r = 6)),  
      strip.text = element_text(size = 38.4, face = "bold"),  
      strip.background = element_rect(fill = "gray95", color = "black", size = 2.4), 
      plot.title = element_text(hjust = 0.5, size = 60, face = "bold", margin = margin(b = 36)),  
      panel.spacing = unit(0.9, "lines"), 
      panel.border = element_rect(color = "black", fill = NA, size = 2.4), 
      plot.margin = margin(t = 24, r = 24, b = 24, l = 24),  
      legend.position = "right",
      legend.key.height = unit(7.2, "cm"), 
      legend.title = element_text(size = 36, margin = margin(b = 30)),  
      legend.text = element_text(size = 32.4), 
      axis.title.x = element_text(size = 48, face = "bold", margin = margin(t = 36, b = 12)),  
      axis.title.y = element_text(size = 48, face = "bold", angle = 90, margin = margin(r = 3.6)),  
      axis.title.y.right = element_text(margin = margin(l = 24)),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank()
    )
}


for (scenario_number in 1:5) {
  file_path_U <- paste0("results/simulation/simulation_results_S", scenario_number, "U.csv")
  file_path_B <- paste0("results/simulation/simulation_results_S", scenario_number, "B.csv")
  save_dir <- "figures"
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  data_B <- load_and_prepare_data(file_path_B, paste0("S", scenario_number, "B"))
  data_U <- load_and_prepare_data(file_path_U, paste0("S", scenario_number, "U"))
  data_combined <- rbind(prepare_data(data_B), prepare_data(data_U))
  
  data_combined$Scenario <- factor(data_combined$Scenario, 
                                   levels = c(paste0("Scenario ", scenario_number, "(U)"), 
                                              paste0("Scenario ", scenario_number, "(B)")))
  
  plot <- generate_plot_combined(data_combined, scenario_number)
  
  save_plot_as_pdf <- function(plot, filename, aspect_ratio = 2.0, height = 20) {
    width <- height * aspect_ratio
    cairo_pdf(filename, width = width, height = height, family = "serif")
    print(plot)
    dev.off()
  }
  
  save_plot_as_pdf(plot, file.path(save_dir, paste0("Scenario", scenario_number, ".pdf")))
}