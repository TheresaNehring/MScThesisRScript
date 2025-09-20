# HOLDING EFFECT ANALYSIS WITH OUTLIER SCREENING
# Statistical tests comparing intake and release measurements
# Experimental group: all variables | Validation group: morphometrics only
# Note: Run DatasetPrep script first!

library(tidyverse)
library(writexl)

# Define variable types for each variable
variable_types <- data.frame(
  Variable = c(
    "Mass [g]", "SCL [mm]", "SCW [mm]", "BD [mm]", "BCI",
    "Glucose [mg/dL]", "PCV [%]", "Lactate [mmol/L]", "BHB [mmol/L]", 
    "Total protein [g/dL]", "Corticosterone [ng/mL]",
    "WBC estimate [× 10³/μL]", "Mature heterophils [× 10³/μL]", "Immature heterophils [× 10³/μL]",
    "Total heterophils [× 10³/μL]", "Lymphocytes [× 10³/μL]", "Heterophil:lymphocyte ratio",
    "Monocytes [× 10³/μL]", "Eosinophils [× 10³/μL]", "Basophils [× 10³/μL]", 
    "Immature RBC/100 mature RBC"
  ),
  Variable_Type = c(
    rep("Morphometrics", 5),
    rep("Blood Chemistry", 6),
    rep("Hematology", 10)
  )
)

# Define morphometric variables for validation group
morphometric_variables <- c("Mass [g]", "SCL [mm]", "SCW [mm]", "BD [mm]", "BCI")

# Helper functions
get_significance <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}

format_p_with_sig <- function(p_value) {
  if (is.na(p_value)) return("-")
  formatted <- if (p_value < 0.001) "<0.001" else sprintf("%.3f", p_value)
  sig <- get_significance(p_value)
  return(paste0(formatted, sig))
}

# Outlier detection using IQR method (3×IQR) for differences
detect_outliers_in_differences <- function(differences, hatchling_ids, variable_name) {
  if (length(differences) < 3) return(data.frame())
  
  # Calculate IQR outliers
  q1 <- quantile(differences, 0.25, na.rm = TRUE)
  q3 <- quantile(differences, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 3 * iqr
  upper_bound <- q3 + 3 * iqr
  
  # Identify outliers
  is_outlier <- differences < lower_bound | differences > upper_bound
  
  if (any(is_outlier)) {
    outlier_data <- data.frame(
      Variable = variable_name,
      `Hatchling ID` = hatchling_ids[is_outlier],
      Difference = round(differences[is_outlier], 3),
      `Lower Bound` = round(lower_bound, 3),
      `Upper Bound` = round(upper_bound, 3),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    return(outlier_data)
  } else {
    return(data.frame())
  }
}

# Main analysis function
perform_holding_analysis <- function(data, variables_to_analyze, group_name) {
  
  # Create wide format for paired analysis
  wide_data <- data %>%
    dplyr::select(Hatchling.ID, Sampling.timepoint, all_of(variables_to_analyze)) %>%
    pivot_wider(id_cols = Hatchling.ID, 
                names_from = Sampling.timepoint,
                values_from = all_of(variables_to_analyze))
  
  initial_results <- list()
  all_outliers <- list()
  differences_data <- list()
  
  for (var in variables_to_analyze) {
    
    # Create column names for intake and release timepoints
    col_intake <- paste0(var, "_Intake")
    col_release <- paste0(var, "_Release")
    
    # Check if both columns exist in the data
    if (col_intake %in% names(wide_data) && col_release %in% names(wide_data)) {
      
      # Filter for complete pairs (both intake and release measurements)
      complete_pairs <- wide_data %>%
        filter(!is.na(.data[[col_intake]]) & !is.na(.data[[col_release]]))
      
      n_pairs <- nrow(complete_pairs)
      
      if (n_pairs >= 3) {
        # Get paired values
        intake_values <- complete_pairs[[col_intake]]
        release_values <- complete_pairs[[col_release]]
        differences <- release_values - intake_values
        
        # Skip if all differences are zero or identical
        if (all(differences == 0, na.rm = TRUE) || 
            length(unique(differences[!is.na(differences)])) <= 1) {
          
          initial_results[[var]] <- data.frame(
            Variable = var,
            Test = "No variation",
            `p-value` = "-",
            `Effect size` = 0,
            `Effect direction` = "No Change",
            `n (I/R)` = paste(n_pairs, n_pairs, sep = "/"),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
          next
        }
        
        # Detect outliers in differences
        outliers <- detect_outliers_in_differences(differences, complete_pairs$Hatchling.ID, var)
        if (nrow(outliers) > 0) {
          all_outliers[[var]] <- outliers
        }
        
        # Store differences data for later use
        differences_data[[var]] <- list(
          data = complete_pairs,
          differences = differences,
          intake_values = intake_values,
          release_values = release_values
        )
        
        # Test normality of differences
        normality_test <- tryCatch({
          if (length(differences) >= 3 && var(differences, na.rm = TRUE) > 0) {
            shapiro.test(differences)
          } else {
            list(p.value = NA)
          }
        }, error = function(e) list(p.value = NA))
        
        is_normal <- !is.na(normality_test$p.value) && normality_test$p.value > 0.05
        
        # Perform appropriate statistical test
        if (is_normal) {
          test_result <- t.test(intake_values, release_values, paired = TRUE)
          sd_diff <- sd(differences, na.rm = TRUE)
          effect_size <- if (sd_diff > 0) mean(differences, na.rm = TRUE) / sd_diff else 0
          test_name <- "Paired t-test"
        } else {
          test_result <- wilcox.test(intake_values, release_values, paired = TRUE, exact = FALSE)
          ranks <- rank(abs(differences), ties.method = "average", na.last = "keep")
          pos_ranks <- sum(ranks[differences > 0], na.rm = TRUE)
          neg_ranks <- sum(ranks[differences < 0], na.rm = TRUE)
          total_ranks <- pos_ranks + neg_ranks
          effect_size <- if (total_ranks > 0) (pos_ranks - neg_ranks) / total_ranks else 0
          test_name <- "Wilcoxon signed-rank"
        }
        
        # Determine effect direction
        mean_diff <- mean(differences, na.rm = TRUE)
        effect_direction <- if (mean_diff > 0) {
          if (test_result$p.value < 0.05) "I < R" else "I < R (ns)"
        } else if (mean_diff < 0) {
          if (test_result$p.value < 0.05) "I > R" else "I > R (ns)"
        } else {
          "No Change"
        }
        
        initial_results[[var]] <- data.frame(
          Variable = var,
          Test = test_name,
          `p-value` = format_p_with_sig(test_result$p.value),
          `Effect size` = round(effect_size, 3),
          `Effect direction` = effect_direction,
          `n (I/R)` = paste(n_pairs, n_pairs, sep = "/"),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        
      } else {
        initial_results[[var]] <- data.frame(
          Variable = var,
          Test = "Insufficient data",
          `p-value` = "-",
          `Effect size` = NA,
          `Effect direction` = "Cannot determine",
          `n (I/R)` = paste(n_pairs, n_pairs, sep = "/"),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Analysis After Outlier Removal
  combined_outliers <- bind_rows(all_outliers)
  final_results <- list()
  
  for (var in names(differences_data)) {
    if (is.null(differences_data[[var]])) {
      next
    }
    
    original_data <- differences_data[[var]]$data
    
    # Identify outliers for this variable
    var_outliers <- if(nrow(combined_outliers) > 0) {
      combined_outliers %>%
        filter(Variable == var) %>%
        pull(`Hatchling ID`)
    } else {
      character(0)
    }
    
    # Remove outliers if any
    if (length(var_outliers) > 0) {
      cleaned_data <- original_data %>%
        filter(!Hatchling.ID %in% var_outliers)
    } else {
      cleaned_data <- original_data
    }
    
    # Check if we still have enough data
    if (nrow(cleaned_data) < 3) {
      next
    }
    
    # Recalculate with cleaned data
    col_intake <- paste0(var, "_Intake")
    col_release <- paste0(var, "_Release")
    
    intake_clean <- cleaned_data[[col_intake]]
    release_clean <- cleaned_data[[col_release]]
    differences_clean <- release_clean - intake_clean
    
    # Test normality of cleaned differences
    normality_test_clean <- tryCatch({
      if (length(differences_clean) >= 3 && var(differences_clean, na.rm = TRUE) > 0) {
        shapiro.test(differences_clean)
      } else {
        list(p.value = NA)
      }
    }, error = function(e) list(p.value = NA))
    
    is_normal_clean <- !is.na(normality_test_clean$p.value) && normality_test_clean$p.value > 0.05
    
    # Perform appropriate statistical test on cleaned data
    if (is_normal_clean) {
      final_test_result <- t.test(intake_clean, release_clean, paired = TRUE)
      sd_diff_clean <- sd(differences_clean, na.rm = TRUE)
      final_effect_size <- if (sd_diff_clean > 0) mean(differences_clean, na.rm = TRUE) / sd_diff_clean else 0
      test_name_final <- "Paired t-test"
    } else {
      final_test_result <- wilcox.test(intake_clean, release_clean, paired = TRUE, exact = FALSE)
      ranks_clean <- rank(abs(differences_clean), ties.method = "average", na.last = "keep")
      pos_ranks_clean <- sum(ranks_clean[differences_clean > 0], na.rm = TRUE)
      neg_ranks_clean <- sum(ranks_clean[differences_clean < 0], na.rm = TRUE)
      total_ranks_clean <- pos_ranks_clean + neg_ranks_clean
      final_effect_size <- if (total_ranks_clean > 0) (pos_ranks_clean - neg_ranks_clean) / total_ranks_clean else 0
      test_name_final <- "Wilcoxon signed-rank"
    }
    
    # Determine final effect direction
    mean_diff_final <- mean(differences_clean, na.rm = TRUE)
    final_effect_direction <- if (mean_diff_final > 0) {
      if (final_test_result$p.value < 0.05) "I < R" else "I < R (ns)"
    } else if (mean_diff_final < 0) {
      if (final_test_result$p.value < 0.05) "I > R" else "I > R (ns)"
    } else {
      "No Change"
    }
    
    final_results[[var]] <- data.frame(
      Variable = var,
      Test = test_name_final,
      `p-value` = format_p_with_sig(final_test_result$p.value),
      `Effect size` = round(final_effect_size, 3),
      `Effect direction` = final_effect_direction,
      `n (I/R)` = paste(nrow(cleaned_data), nrow(cleaned_data), sep = "/"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results and add variable types
  initial_table <- bind_rows(initial_results) %>%
    left_join(variable_types, by = "Variable") %>%
    arrange(factor(Variable, levels = variables_to_analyze)) %>%
    dplyr::select(
      `Variable type` = Variable_Type,
      Variable,
      Test,
      `p-value`,
      `Effect size`,
      `Effect direction`,
      `n (I/R)`
    )
  
  final_results_table <- bind_rows(final_results) %>%
    left_join(variable_types, by = "Variable") %>%
    arrange(factor(Variable, levels = variables_to_analyze)) %>%
    dplyr::select(
      `Variable type` = Variable_Type,
      Variable,
      Test,
      `p-value`,
      `Effect size`,
      `Effect direction`,
      `n (I/R)`
    )
  
  # Summary of significant effects
  significant_effects_initial <- initial_table %>%
    filter(grepl("\\*", `p-value`)) %>%
    filter(`Effect direction` %in% c("I < R", "I > R"))
  
  significant_effects_final <- final_results_table %>%
    filter(grepl("\\*", `p-value`)) %>%
    filter(`Effect direction` %in% c("I < R", "I > R"))
  
  return(list(
    initial_results = initial_table,
    outliers = combined_outliers,
    final_results = final_results_table,
    significant_initial = significant_effects_initial,
    significant_final = significant_effects_final
  ))
}

# MAIN ANALYSIS

# Filter data for each group
experimental_data <- stragglergreenturtles %>%
  filter(Group == "Experimental")

validation_data <- stragglergreenturtles %>%
  filter(Group == "Validation")

# Perform analysis for experimental group (all variables)
exp_results <- perform_holding_analysis(experimental_data, ordered_variables, "Experimental")

# Perform analysis for validation group (morphometrics only)
val_results <- perform_holding_analysis(validation_data, morphometric_variables, "Validation")

# EXPORT RESULTS

# Create output directories
exp_output_dir <- "Holding Effect"
val_output_dir <- "Validation Holding Effect"

if (!dir.exists(exp_output_dir)) {
  dir.create(exp_output_dir, recursive = TRUE)
}
if (!dir.exists(val_output_dir)) {
  dir.create(val_output_dir, recursive = TRUE)
}

# Export experimental group results
exp_export_list <- list(
  "Initial_Results" = exp_results$initial_results,
  "Outliers_Identified" = if(nrow(exp_results$outliers) > 0) exp_results$outliers else data.frame(Note = "No outliers identified"),
  "Final_Results_After_Outlier_Removal" = if(nrow(exp_results$final_results) > 0) exp_results$final_results else data.frame(Note = "No final results available"),
  "Significant_Effects_Initial" = if(nrow(exp_results$significant_initial) > 0) exp_results$significant_initial else data.frame(Note = "No significant effects in initial analysis"),
  "Significant_Effects_Final" = if(nrow(exp_results$significant_final) > 0) exp_results$significant_final else data.frame(Note = "No significant effects after outlier removal")
)

write_xlsx(exp_export_list, file.path(exp_output_dir, "holding_effect_analysis.xlsx"))

# Export validation group results
val_export_list <- list(
  "Initial_Results" = val_results$initial_results,
  "Outliers_Identified" = if(nrow(val_results$outliers) > 0) val_results$outliers else data.frame(Note = "No outliers identified"),
  "Final_Results_After_Outlier_Removal" = if(nrow(val_results$final_results) > 0) val_results$final_results else data.frame(Note = "No final results available"),
  "Significant_Effects_Initial" = if(nrow(val_results$significant_initial) > 0) val_results$significant_initial else data.frame(Note = "No significant effects in initial analysis"),
  "Significant_Effects_Final" = if(nrow(val_results$significant_final) > 0) val_results$significant_final else data.frame(Note = "No significant effects after outlier removal")
)

write_xlsx(val_export_list, file.path(val_output_dir, "validation_holding_effect_morphometrics_analysis.xlsx"))