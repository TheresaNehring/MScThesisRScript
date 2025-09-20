# HISTOGRAMS VISUALIZATION SCRIPT
# Purpose: Create comprehensive histograms showing data distributions

library(tidyverse)
library(ggplot2)
library(patchwork)

color_palette <- list(
  experimental_group = "#8B67C4",
  validation_group = "#6BAF92", 
  intake = "#E6C84E",
  release = "#5B8AC7",
  morphometrics = "#A07855",
  blood_chemistry = "#E6954E", 
  hematology = "#E67C74"
)

turtles_labeled <- stragglergreenturtles %>%
  mutate(
    Group = factor(Group, levels = c("Experimental", "Validation")),
    Sampling.timepoint = factor(Sampling.timepoint, levels = c("Intake", "Release"))
  )

histogram_variable_units <- setNames(ordered_variables, ordered_variables)
group_labels <- c("Experimental" = "Experimental", "Validation" = "Validation")
timepoint_labels <- c("Intake" = "Intake", "Release" = "Release")

output_dir <- "Data Quality and Characterization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

histograms_dir <- file.path(output_dir, "Histograms")
if (!dir.exists(histograms_dir)) dir.create(histograms_dir)

create_histogram <- function(var) {
  if (all(is.na(turtles_labeled[[var]]))) return(NULL)
  
  ref_vals <- reference_intervals %>% filter(Variable == var)
  
  p <- ggplot(turtles_labeled, aes(x = .data[[var]])) +
    # Histogram with grey fill and black outline
    geom_histogram(aes(y = ..density..), 
                   fill = "grey80", color = "black", bins = 15, alpha = 0.7) +
    # Density line
    geom_density(color = "blue", linewidth = 1) +
    # Reference intervals (if available)
    {if (nrow(ref_vals) > 0 && !is.na(ref_vals$Lower) && !is.na(ref_vals$Upper)) {
      list(
        geom_vline(xintercept = ref_vals$Lower, linetype = "dashed", 
                   color = "red", linewidth = 0.7, alpha = 0.6),
        geom_vline(xintercept = ref_vals$Upper, linetype = "dashed", 
                   color = "red", linewidth = 0.7, alpha = 0.6)
      )
    }} +
    facet_grid(rows = vars(Group), cols = vars(Sampling.timepoint), 
               labeller = labeller(Group = group_labels, 
                                   Sampling.timepoint = timepoint_labels)) +
    labs(x = histogram_variable_units[[var]], y = "Density") +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text.x = element_text(face = "bold", size = 9),
      strip.text.y = element_text(face = "bold", size = 9, angle = 270),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  return(p)
}

# Generate individual histograms
for (var in ordered_variables) {
  p <- create_histogram(var)
  if (!is.null(p)) {
    clean_name <- get_clean_filename(var)
    ggsave(filename = file.path(histograms_dir, paste0("Histogram_", clean_name, ".pdf")), 
           plot = p, width = 8, height = 6)
  }
}

# Create morphometric histograms grid
morpho_hist_list <- list()
for (var in morphometrics_variables) {
  p <- create_histogram(var)
  if (!is.null(p)) morpho_hist_list[[var]] <- p
}

if (length(morpho_hist_list) > 0) {
  combined_morpho_hist <- wrap_plots(morpho_hist_list, ncol = 3) +
    plot_annotation(
      title = "Morphometrics - Histograms", 
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  ggsave(filename = file.path(histograms_dir, "morphometric_histograms_grid.pdf"), 
         plot = combined_morpho_hist, width = 15, height = 10)
}

# Create blood chemistry histograms grid
blood_chem_hist_list <- list()
for (var in blood_chemistry_variables) {
  p <- create_histogram(var)
  if (!is.null(p)) blood_chem_hist_list[[var]] <- p
}

if (length(blood_chem_hist_list) > 0) {
  combined_blood_chem_hist <- wrap_plots(blood_chem_hist_list, ncol = 3) +
    plot_annotation(
      title = "Blood Chemistry - Histograms", 
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  ggsave(filename = file.path(histograms_dir, "blood_chemistry_histograms_grid.pdf"), 
         plot = combined_blood_chem_hist, width = 15, height = 10)
}

# Create hematology histograms grid
hematology_hist_list <- list()
for (var in hematology_variables) {
  p <- create_histogram(var)
  if (!is.null(p)) hematology_hist_list[[var]] <- p
}

if (length(hematology_hist_list) > 0) {
  combined_hematology_hist <- wrap_plots(hematology_hist_list, ncol = 3) +
    plot_annotation(
      title = "Hematology - Histograms", 
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  ggsave(filename = file.path(histograms_dir, "hematology_histograms_grid.pdf"), 
         plot = combined_hematology_hist, width = 15, height = 14)
}

# Create comprehensive histograms grid
all_hist_list <- list()
for (var in ordered_variables) {
  p <- create_histogram(var)
  if (!is.null(p)) all_hist_list[[var]] <- p
}

if (length(all_hist_list) > 0) {
  comprehensive_hist_grid <- wrap_plots(all_hist_list, ncol = 3) +
    plot_annotation(title = "All Variables - Histograms", 
                    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = "black")))
  
  ggsave(filename = file.path(histograms_dir, "comprehensive_histograms_grid.pdf"), 
         plot = comprehensive_hist_grid, width = 18, height = 24)
}