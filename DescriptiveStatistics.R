# DESCRIPTIVE STATISTICS
# Purpose: Generate descriptive statistics for all variables by group and timepoint

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

descriptive_stats <- stragglergreenturtles %>%
  pivot_longer(cols = all_of(ordered_variables), 
               names_to = "Variable", 
               values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = ordered_variables)) %>%
  
  group_by(Group, Sampling.timepoint, Variable) %>%
  
  summarise(
    Mean = ifelse(sum(!is.na(Value)) > 0, mean(Value, na.rm = TRUE), NA),
    SD = ifelse(sum(!is.na(Value)) > 1, sd(Value, na.rm = TRUE), NA),
    Median = ifelse(sum(!is.na(Value)) > 0, median(Value, na.rm = TRUE), NA),
    Min = ifelse(sum(!is.na(Value)) > 0, min(Value, na.rm = TRUE), NA),
    Max = ifelse(sum(!is.na(Value)) > 0, max(Value, na.rm = TRUE), NA),
    n = sum(!is.na(Value)),
    .groups = "drop"
  ) %>%
  
  mutate(
    Mean_SD = ifelse(is.na(Mean) | is.na(SD), "NA", 
                     paste0(round(Mean, 2), " ± ", round(SD, 2))),
    Median = ifelse(is.na(Median), NA, round(Median, 2)),
    Min_Max = ifelse(is.na(Min) | is.na(Max), "NA", 
                     paste0(round(Min, 2), " - ", round(Max, 2)))
  ) %>%
  
  left_join(variable_types, by = "Variable") %>%
  
  arrange(Group, Sampling.timepoint, factor(Variable, levels = ordered_variables)) %>%
  
  select(
    Group,
    `Sampling timepoint` = Sampling.timepoint, 
    `Variable type` = Variable_Type,
    Variable, 
    `Mean / SD` = Mean_SD, 
    Median, 
    `Min-Max` = Min_Max, 
    n
  )

experimental_descriptives <- descriptive_stats %>% 
  filter(Group == "Experimental") %>%
  select(-Group)

validation_descriptives <- descriptive_stats %>% 
  filter(Group == "Validation") %>%
  select(-Group)

# Create output directory
output_dir <- "Data Quality and Characterization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Export results
write_xlsx(list(
  "Experimental_Group" = experimental_descriptives,
  "Validation_Group" = validation_descriptives,
  "All_Stats" = descriptive_stats
), file.path(output_dir, "descriptive_statistics.xlsx"))