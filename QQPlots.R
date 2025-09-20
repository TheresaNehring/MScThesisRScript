# QQ PLOTS VISUALIZATION SCRIPT
# Purpose: Create quantile-quantile plots to assess normality of data distributions

library(tidyverse)
library(ggplot2)
library(patchwork)

turtles_labeled <- stragglergreenturtles %>%
  mutate(
    Group = factor(Group, levels = c("Experimental", "Validation")),
    Sampling.timepoint = factor(Sampling.timepoint, levels = c("Intake", "Release"))
  )

qq_variable_units <- setNames(ordered_variables, ordered_variables)
group_labels <- c("Experimental" = "Experimental", "Validation" = "Validation")
timepoint_labels <- c("Intake" = "Intake", "Release" = "Release")

output_dir <- "Data Quality and Characterization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

qqplots_dir <- file.path(output_dir, "QQ_Plots")
if (!dir.exists(qqplots_dir)) dir.create(qqplots_dir)

create_qq_plot <- function(var) {
  if (all(is.na(turtles_labeled[[var]]))) return(NULL)
  
  ref_vals <- reference_intervals %>% filter(Variable == var)
  
  p <- ggplot(turtles_labeled, aes(sample = .data[[var]])) +
    stat_qq() +
    stat_qq_line() +
    {if (nrow(ref_vals) > 0 && !is.na(ref_vals$Lower) && !is.na(ref_vals$Upper)) {
      list(
        geom_hline(yintercept = ref_vals$Lower, linetype = "dashed", 
                   color = "red", linewidth = 0.7, alpha = 0.6),
        geom_hline(yintercept = ref_vals$Upper, linetype = "dashed", 
                   color = "red", linewidth = 0.7, alpha = 0.6)
      )
    }} +
    facet_grid(rows = vars(Group), cols = vars(Sampling.timepoint), 
               labeller = labeller(Group = group_labels, 
                                   Sampling.timepoint = timepoint_labels)) +
    labs(x = "Theoretical Quantiles", y = qq_variable_units[[var]]) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.text.y = element_text(face = "bold", size = 9, angle = 270),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# Generate individual QQ plots
for (var in ordered_variables) {
  p <- create_qq_plot(var)
  if (!is.null(p)) {
    clean_name <- get_clean_filename(var)
    ggsave(filename = file.path(qqplots_dir, paste0("QQ_", clean_name, ".pdf")), 
           plot = p, width = 8, height = 6)
  }
}

# Create morphometric QQ plots grid
morpho_qq_list <- list()
for (var in morphometrics_variables) {
  p <- create_qq_plot(var)
  if (!is.null(p)) morpho_qq_list[[var]] <- p
}

if (length(morpho_qq_list) > 0) {
  combined_morpho_qq <- wrap_plots(morpho_qq_list, ncol = 3) +
    plot_annotation(title = "Morphometrics - QQ Plots", 
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  ggsave(filename = file.path(qqplots_dir, "morphometric_qq_grid.pdf"), 
         plot = combined_morpho_qq, width = 15, height = 10)
}

# Create blood chemistry QQ plots grid
blood_chem_qq_list <- list()
for (var in blood_chemistry_variables) {
  p <- create_qq_plot(var)
  if (!is.null(p)) blood_chem_qq_list[[var]] <- p
}

if (length(blood_chem_qq_list) > 0) {
  combined_blood_chem_qq <- wrap_plots(blood_chem_qq_list, ncol = 3) +
    plot_annotation(title = "Blood Chemistry - QQ Plots", 
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  ggsave(filename = file.path(qqplots_dir, "blood_chemistry_qq_grid.pdf"), 
         plot = combined_blood_chem_qq, width = 15, height = 10)
}

# Create hematology QQ plots grid
hematology_qq_list <- list()
for (var in hematology_variables) {
  p <- create_qq_plot(var)
  if (!is.null(p)) hematology_qq_list[[var]] <- p
}

if (length(hematology_qq_list) > 0) {
  combined_hematology_qq <- wrap_plots(hematology_qq_list, ncol = 3) +
    plot_annotation(title = "Hematology - QQ Plots", 
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  ggsave(filename = file.path(qqplots_dir, "hematology_qq_grid.pdf"), 
         plot = combined_hematology_qq, width = 15, height = 14)
}

# Create comprehensive QQ plots grid
all_qq_list <- list()
for (var in ordered_variables) {
  p <- create_qq_plot(var)
  if (!is.null(p)) all_qq_list[[var]] <- p
}

if (length(all_qq_list) > 0) {
  comprehensive_qq_grid <- wrap_plots(all_qq_list, ncol = 3) +
    plot_annotation(title = "All Variables - QQ Plots", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))
  
  ggsave(filename = file.path(qqplots_dir, "comprehensive_qq_grid.pdf"), 
         plot = comprehensive_qq_grid, width = 18, height = 24)
}