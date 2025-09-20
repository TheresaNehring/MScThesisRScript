# BOXPLOTS VISUALIZATION SCRIPT

library(ggplot2)
library(patchwork)
library(readxl)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)

# Load statistical results
exp_holding_results <- read_excel("Holding Effect/holding_effect_analysis.xlsx", sheet = "Initial_Results")
val_morpho_holding_results <- read_excel("Validation Holding Effect/validation_holding_effect_morphometrics_analysis.xlsx", sheet = "Initial_Results")

# Find clinical grade column
clinical_grade_column <- "Clinical.condition.grade..1.Healthy..2.Compromised..3.Critical."

# Verify the column exists
if(!clinical_grade_column %in% colnames(stragglergreenturtles)) {
  stop("Clinical grade column not found in dataset")
}

# Prepare data with outlier detection
data_long <- stragglergreenturtles %>%
  pivot_longer(cols = all_of(ordered_variables), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Group = case_when(Group == 1 ~ "Experimental", Group == 2 ~ "Validation", TRUE ~ as.character(Group)),
    Sampling.timepoint = case_when(Sampling.timepoint == 0 ~ "Intake", Sampling.timepoint == 1 ~ "Release", TRUE ~ as.character(Sampling.timepoint))
  ) %>%
  mutate(
    Sampling.timepoint = factor(Sampling.timepoint, levels = c("Intake", "Release")),
    Variable = factor(Variable, levels = ordered_variables),
    Variable_with_units = Variable
  ) %>%
  group_by(Group, Sampling.timepoint, Variable) %>%
  mutate(
    Q1 = quantile(Value, 0.25, type = 6, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, type = 6, na.rm = TRUE),
    IQR_val = Q3 - Q1,
    lower_extreme = Q1 - 3 * IQR_val,
    upper_extreme = Q3 + 3 * IQR_val,
    is_severe_outlier = Value < lower_extreme | Value > upper_extreme
  ) %>%
  ungroup()

# Add clinical grade coloring if column exists
if (!is.null(clinical_grade_column)) {
  id_col_in_long <- grep("Hatchling", colnames(data_long), value = TRUE)[1]
  
  data_long <- data_long %>%
    left_join(
      stragglergreenturtles %>% 
        select(Hatchling.ID, Clinical_Grade = all_of(clinical_grade_column)) %>%
        distinct(),
      by = setNames("Hatchling.ID", id_col_in_long)
    ) %>%
    distinct() %>%
    mutate(
      Point_Color_Category = case_when(
        Clinical_Grade == 1 ~ "Healthy",
        Clinical_Grade == 2 ~ "Compromised", 
        Clinical_Grade == 3 ~ "Critical",
        TRUE ~ "Healthy"
      ),
      Point_Shape_Category = ifelse(is_severe_outlier, "Severe Outlier", "Normal")
    )
  
  clinical_colors <- c(
    "Healthy" = "#2E8B57",
    "Compromised" = "#FF8C00",
    "Critical" = "#8B0000"
  )
} else {
  data_long <- data_long %>%
    mutate(
      Point_Color_Category = "Normal",
      Point_Shape_Category = ifelse(is_severe_outlier, "Severe Outlier", "Normal")
    )
  
  clinical_colors <- c("Normal" = "#000000")
}

outlier_shapes <- c("Normal" = 16, "Severe Outlier" = 17)

# Create output directories
output_dir <- "Data Quality and Characterization/Boxplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

color_palette <- list(
  intake = "#E6C84E",
  release = "#5B8AC7",
  morphometrics = "#A07855",
  blood_chemistry = "#E6954E", 
  hematology = "#E67C74"
)

create_legend_grob <- function() {
  clinical_legend_data <- data.frame(
    x = c(1, 2, 3),
    y = rep(1, 3),
    color = c("Healthy", "Compromised", "Critical"),
    shape = rep("Normal", 3)
  )
  
  outlier_legend_data <- data.frame(
    x = c(5, 6),
    y = rep(1, 2),
    color = rep("Normal", 2),
    shape = c("Normal", "Severe Outlier")
  )
  
  legend_plot <- ggplot() +
    geom_point(data = clinical_legend_data, 
               aes(x = x, y = y, color = color, shape = shape), 
               size = 4) +
    geom_point(data = outlier_legend_data, 
               aes(x = x, y = y, color = color, shape = shape), 
               size = 4) +
    annotate("text", x = 1, y = 0.7, label = "Healthy", size = 5, hjust = 0.5, fontface = "bold") +
    annotate("text", x = 2, y = 0.7, label = "Compromised", size = 5, hjust = 0.5, fontface = "bold") +
    annotate("text", x = 3, y = 0.7, label = "Critical", size = 5, hjust = 0.5, fontface = "bold") +
    annotate("text", x = 5, y = 0.7, label = "Normal", size = 5, hjust = 0.5, fontface = "bold") +
    annotate("text", x = 6, y = 0.7, label = "Outlier (>3×IQR)", size = 5, hjust = 0.5, fontface = "bold") +
    annotate("text", x = 2, y = 1.4, label = "Clinical Grade", size = 6, hjust = 0.5, fontface = "bold") +
    annotate("text", x = 5.5, y = 1.4, label = "Outlier Status", size = 6, hjust = 0.5, fontface = "bold") +
    scale_color_manual(values = c(
      "Healthy" = "#2E8B57",
      "Compromised" = "#FF8C00", 
      "Critical" = "#8B0000",
      "Normal" = "#000000"
    ), guide = "none") +
    scale_shape_manual(values = c("Normal" = 16, "Severe Outlier" = 17), guide = "none") +
    coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0.5, 1.6)) +
    theme_void() +
    theme(
      plot.margin = margin(10, 10, 10, 10),
      panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(legend_plot)
}

format_pvalue <- function(pval_string) {
  if (is.na(pval_string) || is.null(pval_string) || pval_string == "-") return(NA)
  
  if (grepl("<0.001", pval_string)) {
    return("P < 0.001***")
  } else if (grepl("\\*", pval_string)) {
    clean_pval <- gsub("\\*+", "", pval_string)
    star_count <- nchar(pval_string) - nchar(clean_pval)
    stars <- paste(rep("*", star_count), collapse = "")
    return(paste0("P = ", clean_pval, stars))
  } else if (pval_string %in% c("No variation", "Insufficient data")) {
    return(NA)
  } else {
    numeric_p <- as.numeric(pval_string)
    if (!is.na(numeric_p)) {
      stars <- if (numeric_p < 0.001) "***" 
      else if (numeric_p < 0.01) "**" 
      else if (numeric_p < 0.05) "*" 
      else ""
      return(paste0("P = ", pval_string, stars))
    } else {
      return(paste0("P = ", pval_string))
    }
  }
}

get_stats_for_variable <- function(var) {
  exp_holding_stat <- exp_holding_results %>% filter(Variable == var)
  val_morpho_stat <- val_morpho_holding_results %>% filter(Variable == var)
  
  results <- list()
  
  if (nrow(exp_holding_stat) > 0 && !is.na(exp_holding_stat$`p-value`) && exp_holding_stat$`p-value` != "-") {
    pval_formatted <- format_pvalue(exp_holding_stat$`p-value`)
    if (!is.na(pval_formatted)) {
      results$experimental_intake_vs_release <- pval_formatted
    }
  }
  
  if (nrow(val_morpho_stat) > 0 && !is.na(val_morpho_stat$`p-value`) && val_morpho_stat$`p-value` != "-") {
    pval_formatted <- format_pvalue(val_morpho_stat$`p-value`)
    if (!is.na(pval_formatted)) {
      results$validation_intake_vs_release <- pval_formatted
    }
  }
  
  return(results)
}

create_boxplot <- function(var, data_subset) {
  ref_vals <- reference_intervals %>% filter(Variable == var)
  var_with_units <- unique(data_subset$Variable_with_units)[1]
  stats <- get_stats_for_variable(var)
  timepoint_colors <- c("Intake" = color_palette$intake, "Release" = color_palette$release)
  
  exp_data <- data_subset %>% filter(Group == "Experimental")
  val_data <- data_subset %>% filter(Group == "Validation")
  
  create_group_plot <- function(group_data, group_name) {
    p <- ggplot(group_data, aes(x = Sampling.timepoint, y = Value)) +
      geom_boxplot(aes(fill = Sampling.timepoint), outlier.shape = NA, alpha = 0.8, 
                   position = position_dodge(width = 0.75), linewidth = 0.8,  # Increased line width
                   width = 0.6, color = "black") +
      geom_point(aes(color = Point_Color_Category, shape = Point_Shape_Category), alpha = 0.8, size = 3.5,  # Increased point size
                 position = position_jitter(width = 0.15)) +
      scale_color_manual(values = clinical_colors, guide = "none") +
      scale_shape_manual(values = outlier_shapes, guide = "none") +
      scale_fill_manual(values = timepoint_colors, guide = "none") +
      theme_classic(base_size = 20) +  # Increased from 16 to 20
      labs(y = if(group_name == "Experimental") var_with_units else "", 
           x = "Sampling timepoint") +
      ggtitle(paste("Group:", group_name)) +
      theme(
        axis.text.x = element_text(size = 18, color = "black"),  # Removed bold
        axis.text.y = element_text(size = 17, color = "black", face = "bold"),
        axis.title.y = element_text(size = 19, color = "black", face = "bold"),
        axis.title.x = element_text(size = 18, color = "black", face = "bold"),
        plot.margin = margin(30, 20, 20, 20),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.ticks.length = unit(0.3, "cm")
      )
    
    if (nrow(ref_vals) > 0 && !is.na(ref_vals$Lower) && !is.na(ref_vals$Upper)) {
      p <- p +
        geom_hline(yintercept = ref_vals$Lower, linetype = "dashed", color = "red", alpha = 0.8, linewidth = 0.9) +  # Increased line width and alpha
        geom_hline(yintercept = ref_vals$Upper, linetype = "dashed", color = "red", alpha = 0.8, linewidth = 0.9)   # Increased line width and alpha
    }
    
    return(p)
  }
  
  p_exp <- create_group_plot(exp_data, "Experimental")
  p_val <- create_group_plot(val_data, "Validation")
  
  y_max <- max(data_subset$Value, na.rm = TRUE)
  y_min <- min(data_subset$Value, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (!is.null(stats$experimental_intake_vs_release) && !is.na(stats$experimental_intake_vs_release)) {
    is_significant <- grepl("\\*", stats$experimental_intake_vs_release)
    exp_label <- ifelse(is_significant, stats$experimental_intake_vs_release, "NS")
    bracket_y1 <- y_max + 0.12 * y_range
    
    p_exp <- p_exp + 
      annotate("segment", x = 1, xend = 2, y = bracket_y1, yend = bracket_y1, color = "black", linewidth = 0.8) +  # Increased line width
      annotate("segment", x = 1, xend = 1, y = bracket_y1 - 0.02 * y_range, yend = bracket_y1, color = "black", linewidth = 0.8) +
      annotate("segment", x = 2, xend = 2, y = bracket_y1 - 0.02 * y_range, yend = bracket_y1, color = "black", linewidth = 0.8) +
      annotate("text", x = 1.5, y = bracket_y1 + 0.04 * y_range, label = exp_label, size = 6, hjust = 0.5, color = "black", fontface = "bold")  # Increased from 4.5 to 6
  }
  
  if (!is.null(stats$validation_intake_vs_release) && !is.na(stats$validation_intake_vs_release)) {
    is_significant <- grepl("\\*", stats$validation_intake_vs_release)
    val_label <- ifelse(is_significant, stats$validation_intake_vs_release, "NS")
    bracket_y1 <- y_max + 0.12 * y_range
    
    p_val <- p_val + 
      annotate("segment", x = 1, xend = 2, y = bracket_y1, yend = bracket_y1, color = "black", linewidth = 0.8) +  # Increased line width
      annotate("segment", x = 1, xend = 1, y = bracket_y1 - 0.02 * y_range, yend = bracket_y1, color = "black", linewidth = 0.8) +
      annotate("segment", x = 2, xend = 2, y = bracket_y1 - 0.02 * y_range, yend = bracket_y1, color = "black", linewidth = 0.8) +
      annotate("text", x = 1.5, y = bracket_y1 + 0.04 * y_range, label = val_label, size = 6, hjust = 0.5, color = "black", fontface = "bold")  # Increased from 4.5 to 6
  }
  
  legend_grob <- create_legend_grob()
  
  combined_plot <- (p_exp + p_val) / legend_grob + 
    plot_layout(heights = c(4, 1)) +
    plot_annotation(
      title = var_with_units,
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))  # Reduced from 22 to 18
    )
  
  return(combined_plot)
}

# Create grid boxplots with legend
create_boxplot_grid <- function(var, data_subset) {
  ref_vals <- reference_intervals %>% filter(Variable == var)
  var_with_units <- unique(data_subset$Variable_with_units)[1]
  stats <- get_stats_for_variable(var)
  timepoint_colors <- c("Intake" = color_palette$intake, "Release" = color_palette$release)
  
  p <- ggplot(data_subset, aes(x = Sampling.timepoint, y = Value)) +
    geom_boxplot(aes(fill = Sampling.timepoint), outlier.shape = NA, alpha = 0.8, 
                 position = position_dodge(width = 0.75), linewidth = 0.6, 
                 width = 0.6, color = "black") +
    geom_point(aes(color = Point_Color_Category, shape = Point_Shape_Category), alpha = 0.8, size = 2.5, 
               position = position_jitter(width = 0.15)) +
    scale_color_manual(values = c(
      "Healthy" = "#2E8B57",
      "Compromised" = "#FF8C00", 
      "Critical" = "#8B0000",
      "Normal" = "#000000"
    ), guide = "none") +
    scale_shape_manual(values = c("Normal" = 16, "Severe Outlier" = 17), guide = "none") +
    scale_fill_manual(values = timepoint_colors, guide = "none") +
    theme_classic(base_size = 14) +
    labs(y = var_with_units, x = "Sampling timepoint") +
    facet_grid(. ~ Group, labeller = labeller(Group = function(x) paste("Group:", x))) +
    theme(
      axis.text.x = element_text(size = 11, color = "black"),  # Reduced from 12 to 11
      axis.text.y = element_text(size = 10, color = "black"),  # Reduced from 11 to 10
      axis.title.y = element_text(size = 12, color = "black", face = "bold"),  # Reduced from 13 to 12
      axis.title.x = element_text(size = 11, color = "black", face = "bold"),  # Reduced from 12 to 11
      plot.margin = margin(25, 20, 35, 20),  # Increased bottom margin from 30 to 35
      strip.text.x = element_text(size = 12, face = "bold", color = "black"),  # Reduced from 13 to 12
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5)
    )
  
  if (nrow(ref_vals) > 0 && !is.na(ref_vals$Lower) && !is.na(ref_vals$Upper)) {
    p <- p +
      geom_hline(yintercept = ref_vals$Lower, linetype = "dashed", color = "red", alpha = 0.8, linewidth = 0.6) +
      geom_hline(yintercept = ref_vals$Upper, linetype = "dashed", color = "red", alpha = 0.8, linewidth = 0.6)
  }
  
  y_max <- max(data_subset$Value, na.rm = TRUE)
  y_min <- min(data_subset$Value, na.rm = TRUE)
  y_range <- y_max - y_min
  bracket_y1 <- y_max + 0.12 * y_range  # Increased from 0.08 back to 0.12 for more space
  
  annotation_data <- data.frame()
  
  if (!is.null(stats$experimental_intake_vs_release) && !is.na(stats$experimental_intake_vs_release)) {
    is_significant <- grepl("\\*", stats$experimental_intake_vs_release)
    exp_label <- ifelse(is_significant, stats$experimental_intake_vs_release, "NS")
    
    exp_annotation <- data.frame(
      Group = "Experimental", x1 = 1, x2 = 2, y = bracket_y1, label = exp_label
    )
    annotation_data <- rbind(annotation_data, exp_annotation)
  }
  
  if (!is.null(stats$validation_intake_vs_release) && !is.na(stats$validation_intake_vs_release)) {
    is_significant <- grepl("\\*", stats$validation_intake_vs_release)
    val_label <- ifelse(is_significant, stats$validation_intake_vs_release, "NS")
    
    val_annotation <- data.frame(
      Group = "Validation", x1 = 1, x2 = 2, y = bracket_y1, label = val_label
    )
    annotation_data <- rbind(annotation_data, val_annotation)
  }
  
  if (nrow(annotation_data) > 0) {
    p <- p + 
      geom_segment(data = annotation_data, aes(x = x1, xend = x2, y = y, yend = y), color = "black", linewidth = 0.5, inherit.aes = FALSE) +
      geom_segment(data = annotation_data, aes(x = x1, xend = x1, y = y - 0.012 * y_range, yend = y), color = "black", linewidth = 0.5, inherit.aes = FALSE) +
      geom_segment(data = annotation_data, aes(x = x2, xend = x2, y = y - 0.012 * y_range, yend = y), color = "black", linewidth = 0.5, inherit.aes = FALSE) +
      geom_text(data = annotation_data, aes(x = 1.5, y = y + 0.04 * y_range, label = label), size = 3.5, hjust = 0.5, color = "black", fontface = "bold", inherit.aes = FALSE)  # Increased spacing from 0.025 to 0.04
  }
  
  return(p)
}

# Function to create grid with legend
create_grid_with_legend <- function(plot_list, title, color) {
  if (length(plot_list) == 0) return(NULL)
  
  main_grid <- wrap_plots(plot_list, ncol = 3)
  legend_grob <- create_legend_grob()
  
  combined_plot <- main_grid / legend_grob + 
    plot_layout(heights = c(10, 1)) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, color = color),  # Reduced from 22 to 18
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  return(combined_plot)
}
for (var in ordered_variables) {
  data_var <- data_long %>% filter(Variable == var)
  if (nrow(data_var) == 0 || all(is.na(data_var$Value))) next
  
  p <- create_boxplot(var, data_var)
  clean_var_name <- get_clean_filename(var)
  pdf_filename <- file.path(output_dir, paste0("boxplot_", clean_var_name, ".pdf"))
# Generate individual boxplots
for (var in ordered_variables) {
  data_var <- data_long %>% filter(Variable == var)
  if (nrow(data_var) == 0 || all(is.na(data_var$Value))) next
  
  p <- create_boxplot(var, data_var)
  clean_var_name <- get_clean_filename(var)
  pdf_filename <- file.path(output_dir, paste0("boxplot_", clean_var_name, ".pdf"))
  ggsave(filename = pdf_filename, plot = p, width = 14, height = 10, device = "pdf")
}

# Generate morphometrics grid
morphometrics_data <- data_long %>% filter(Variable %in% morphometrics_variables)
if (nrow(morphometrics_data) > 0) {
  plot_list_morpho <- list()
  for (var in morphometrics_variables) {
    data_var <- morphometrics_data %>% filter(Variable == var)
    if (nrow(data_var) > 0 && !all(is.na(data_var$Value))) {
      plot_list_morpho[[var]] <- create_boxplot_grid(var, data_var)
    }
  }
  
  combined_morpho <- create_grid_with_legend(
    plot_list_morpho, 
    "Morphometrics - Boxplots", 
    "black"
  )
  
  if (!is.null(combined_morpho)) {
    ggsave(filename = file.path(output_dir, "morphometrics_grid.pdf"), 
           plot = combined_morpho, width = 18, height = 15, device = "pdf")
  }
}

# Generate blood chemistry grid
blood_chemistry_data <- data_long %>% filter(Variable %in% blood_chemistry_variables)
if (nrow(blood_chemistry_data) > 0) {
  plot_list_blood <- list()
  for (var in blood_chemistry_variables) {
    data_var <- blood_chemistry_data %>% filter(Variable == var)
    if (nrow(data_var) > 0 && !all(is.na(data_var$Value))) {
      plot_list_blood[[var]] <- create_boxplot_grid(var, data_var)
    }
  }
  
  combined_blood <- create_grid_with_legend(
    plot_list_blood, 
    "Blood Chemistry - Boxplots", 
    "black"
  )
  
  if (!is.null(combined_blood)) {
    ggsave(filename = file.path(output_dir, "blood_chemistry_grid.pdf"), 
           plot = combined_blood, width = 18, height = 15, device = "pdf")
  }
}

# Generate hematology grid
hematology_data <- data_long %>% filter(Variable %in% hematology_variables)
if (nrow(hematology_data) > 0) {
  plot_list_hematologic <- list()
  for (var in hematology_variables) {
    data_var <- hematology_data %>% filter(Variable == var)
    if (nrow(data_var) > 0 && !all(is.na(data_var$Value))) {
      plot_list_hematologic[[var]] <- create_boxplot_grid(var, data_var)
    }
  }
  
  combined_hematologic <- create_grid_with_legend(
    plot_list_hematologic, 
    "Hematology - Boxplots", 
    "black"
  )
  
  if (!is.null(combined_hematologic)) {
    ggsave(filename = file.path(output_dir, "hematology_grid.pdf"), 
           plot = combined_hematologic, width = 20, height = 22, device = "pdf")  # Increased dimensions
  }
}

# Generate comprehensive grid
all_plot_list <- list()
for (var in ordered_variables) {
  data_var <- data_long %>% filter(Variable == var)
  if (nrow(data_var) > 0 && !all(is.na(data_var$Value))) {
    all_plot_list[[var]] <- create_boxplot_grid(var, data_var)
  }
}

if (length(all_plot_list) > 0) {
  comprehensive_grid <- create_grid_with_legend(
    all_plot_list, 
    "All Variables - Comprehensive Boxplot Analysis", 
    "black"
  )
  
  if (!is.null(comprehensive_grid)) {
    ggsave(filename = file.path(output_dir, "comprehensive_measurements_grid.pdf"), 
           plot = comprehensive_grid, width = 18, height = 40, device = "pdf")
  }
}
}
