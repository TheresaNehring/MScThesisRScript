# OUTLIER SCREENING
# Purpose: Identify outliers >3×IQR across all variables

library(tidyverse)
library(openxlsx)

stragglergreenturtles <- stragglergreenturtles %>%
  mutate(
    Group_Label = Group,
    Timepoint_Label = Sampling.timepoint,
    Group_Time = paste0(Group, "_", Sampling.timepoint)
  )

# Outlier detection function (3×IQR only)
detect_outliers <- function(x) {
  q1 <- quantile(x, 0.25, type = 6, na.rm = TRUE)
  q3 <- quantile(x, 0.75, type = 6, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower_bound <- q1 - 3 * iqr
  upper_bound <- q3 + 3 * iqr
  
  is_outlier <- x < lower_bound | x > upper_bound
  
  return(is_outlier)
}

# Initialize results
individual_outliers <- data.frame()

# Process each turtle
for (i in 1:nrow(stragglergreenturtles)) {
  turtle_outlier_vars <- c()
  
  for (var in ordered_variables) {
    if (var %in% names(stragglergreenturtles) && !is.na(stragglergreenturtles[i, var])) {
      
      group_time <- stragglergreenturtles$Group_Time[i]
      subset_data <- stragglergreenturtles %>% filter(Group_Time == group_time)
      
      outlier_flags <- detect_outliers(subset_data[[var]])
      
      turtle_id <- stragglergreenturtles$Hatchling.ID[i]
      turtle_position <- which(subset_data$Hatchling.ID == turtle_id)
      
      if (length(turtle_position) > 0 && outlier_flags[turtle_position]) {
        turtle_outlier_vars <- c(turtle_outlier_vars, var)
      }
    }
  }
  
  individual_outliers <- rbind(individual_outliers, data.frame(
    Hatchling_ID = stragglergreenturtles$Hatchling.ID[i],
    Group = stragglergreenturtles$Group_Label[i],
    Sampling_timepoint = stragglergreenturtles$Timepoint_Label[i],
    N_Outlier_Variables = length(turtle_outlier_vars),
    Outlier_Variables = ifelse(length(turtle_outlier_vars) > 0,
                               paste(turtle_outlier_vars, collapse = "; "),
                               "None"),
    stringsAsFactors = FALSE
  ))
}

# Filter for only individuals with outliers
outliers_only <- individual_outliers %>%
  filter(N_Outlier_Variables > 0) %>%
  select(Group, `Sampling timepoint` = Sampling_timepoint, 
         `Hatchling ID` = Hatchling_ID, 
         `Outlier Variables` = Outlier_Variables) %>%
  arrange(Group, `Sampling timepoint`, `Hatchling ID`)

# Export results
output_dir <- "Data Quality and Characterization"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

wb <- createWorkbook()
addWorksheet(wb, "Outliers_3xIQR")
writeData(wb, "Outliers_3xIQR", outliers_only)
setColWidths(wb, "Outliers_3xIQR", cols = 1:4, widths = "auto")

addWorksheet(wb, "Complete_Analysis")
writeData(wb, "Complete_Analysis", individual_outliers)
setColWidths(wb, "Complete_Analysis", cols = 1:ncol(individual_outliers), widths = "auto")

saveWorkbook(wb, file.path(output_dir, "outlier_screening.xlsx"), overwrite = TRUE)