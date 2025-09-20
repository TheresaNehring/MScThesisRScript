# DATA DISTRIBUTION - NORMALITY TESTING
# Purpose: Test normality of raw data and differences (Release - Intake) using Shapiro-Wilk tests

library(tidyverse)
library(writexl)

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

morphometric_variables <- c("Mass [g]", "SCL [mm]", "SCW [mm]", "BD [mm]", "BCI")

format_p_value <- function(p_value) {
  case_when(
    is.na(p_value) ~ "-",
    p_value < 0.001 ~ "< 0.001***",
    p_value < 0.01 ~ paste0(sprintf("%.3f", p_value), "**"),
    p_value < 0.05 ~ paste0(sprintf("%.3f", p_value), "*"),
    TRUE ~ sprintf("%.3f", p_value)
  )
}

safe_shapiro_test <- function(values) {
  values <- na.omit(values)
  
  if (length(values) < 3) return(NA)
  if (length(unique(values)) == 1) return(NA)
  
  tryCatch({
    shapiro.test(values)$p.value
  }, error = function(e) NA)
}


raw_normality_results <- stragglergreenturtles %>%
  pivot_longer(cols = all_of(ordered_variables), names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = ordered_variables)) %>%
  
  group_by(Group, Sampling.timepoint, Variable) %>%
  
  summarise(
    n = sum(!is.na(Value)),
    shapiro_p = safe_shapiro_test(Value),
    .groups = "drop"
  ) %>%
  
  mutate(p_value_formatted = format_p_value(shapiro_p)) %>%
  
  left_join(variable_types, by = "Variable") %>%
  
  arrange(Group, Sampling.timepoint, factor(Variable, levels = ordered_variables)) %>%
  
  select(
    Group,
    `Sampling timepoint` = Sampling.timepoint,
    `Variable type` = Variable_Type,
    Variable,
    `p-value` = p_value_formatted,
    n
  )

experimental_raw_normality <- raw_normality_results %>% 
  filter(Group == "Experimental") %>%
  select(-Group)

validation_raw_normality <- raw_normality_results %>% 
  filter(Group == "Validation") %>%
  select(-Group)

non_normal_raw <- stragglergreenturtles %>%
  pivot_longer(cols = all_of(ordered_variables), names_to = "Variable", values_to = "Value") %>%
  group_by(Group, Sampling.timepoint, Variable) %>%
  summarise(
    n = sum(!is.na(Value)),
    shapiro_p = safe_shapiro_test(Value),
    .groups = "drop"
  ) %>%
  filter(!is.na(shapiro_p) & shapiro_p < 0.05) %>%
  left_join(variable_types, by = "Variable") %>%
  mutate(p_value_formatted = format_p_value(shapiro_p)) %>%
  select(Group, `Sampling timepoint` = Sampling.timepoint, 
         `Variable type` = Variable_Type, Variable, 
         `p-value` = p_value_formatted, n) %>%
  arrange(Group, `Sampling timepoint`, Variable)


test_difference_normality <- function(data, group_name, variables_to_test) {
  
  wide_data <- data %>%
    filter(Group == group_name) %>%
    select(Hatchling.ID, Sampling.timepoint, all_of(variables_to_test)) %>%
    pivot_wider(id_cols = Hatchling.ID, 
                names_from = Sampling.timepoint,
                values_from = all_of(variables_to_test))
  
  results_list <- list()
  
  for (var in variables_to_test) {
    col_intake <- paste0(var, "_Intake")
    col_release <- paste0(var, "_Release")
    
    if (col_intake %in% names(wide_data) && col_release %in% names(wide_data)) {
      
      n_intake <- sum(!is.na(wide_data[[col_intake]]))
      n_release <- sum(!is.na(wide_data[[col_release]]))
      
      complete_pairs <- wide_data %>%
        filter(!is.na(.data[[col_intake]]) & !is.na(.data[[col_release]]))
      
      n_pairs <- nrow(complete_pairs)
      
      if (n_pairs >= 3) {
        differences <- complete_pairs[[col_release]] - complete_pairs[[col_intake]]
        
        if (length(unique(differences)) == 1) {
          shapiro_p <- NA
          test_recommended <- "All differences identical"
        } else {
          shapiro_p <- safe_shapiro_test(differences)
          is_normal <- ifelse(is.na(shapiro_p), NA, shapiro_p > 0.05)
          test_recommended <- ifelse(is_normal %in% TRUE, "Paired t-test", "Wilcoxon signed rank")
        }
      } else {
        shapiro_p <- NA
        test_recommended <- "Insufficient data"
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        Variable = var,
        `p-value` = format_p_value(shapiro_p),
        `n (I/R)` = paste0(n_intake, "/", n_release),
        `Recommended test` = test_recommended,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  
  return(do.call(rbind, results_list))
}

exp_difference_results <- test_difference_normality(stragglergreenturtles, "Experimental", ordered_variables) %>%
  left_join(variable_types, by = "Variable") %>%
  arrange(factor(Variable, levels = ordered_variables)) %>%
  select(`Variable type` = Variable_Type, Variable, `p-value`, `n (I/R)`, `Recommended test`)

val_difference_results <- test_difference_normality(stragglergreenturtles, "Validation", morphometric_variables) %>%
  left_join(variable_types, by = "Variable") %>%
  arrange(factor(Variable, levels = morphometric_variables)) %>%
  select(`Variable type` = Variable_Type, Variable, `p-value`, `n (I/R)`, `Recommended test`)

exp_non_normal_diff <- exp_difference_results %>%
  filter(grepl("\\*", `p-value`)) %>%
  filter(`Recommended test` == "Wilcoxon signed rank")

val_non_normal_diff <- val_difference_results %>%
  filter(grepl("\\*", `p-value`)) %>%
  filter(`Recommended test` == "Wilcoxon signed rank")


output_dir <- "Data Quality and Characterization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

export_list <- list(
  "Raw_Experimental" = experimental_raw_normality,
  "Raw_Validation" = validation_raw_normality,
  "Raw_All_Groups" = raw_normality_results,
  "Raw_Non_Normal" = if(nrow(non_normal_raw) > 0) non_normal_raw else 
    data.frame(Note = "All raw variables follow normal distribution"),
  
  "Diff_Experimental" = exp_difference_results,
  "Diff_Validation" = val_difference_results,
  "Diff_Non_Normal_Exp" = if(nrow(exp_non_normal_diff) > 0) exp_non_normal_diff else
    data.frame(Note = "All experimental differences follow normal distribution"),
  "Diff_Non_Normal_Val" = if(nrow(val_non_normal_diff) > 0) val_non_normal_diff else
    data.frame(Note = "All validation differences follow normal distribution")
)

write_xlsx(export_list, file.path(output_dir, "datadistribution_normalitytesting.xlsx"))