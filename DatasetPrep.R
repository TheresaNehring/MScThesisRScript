# DATASET PREPARATION 
# Purpose: Load and prepare the sea turtle dataset for analysis

# Set working directory
setwd("~/Documents/Uni/MSc Molecular and Microbial Biology UALG/Masters thesis/Data/Data + Results")

# Load required libraries
library(tidyverse)

# Load dataset
stragglergreenturtles <- read.csv("GreenTurtlesDataR.csv")
names(stragglergreenturtles) <- trimws(names(stragglergreenturtles))

# Variable mapping and standardization
variables <- c(
  "Mass..g." = "Mass [g]", 
  "SCL..mm." = "SCL [mm]", 
  "SCW..mm." = "SCW [mm]", 
  "BD..mm." = "BD [mm]",
  "BCI" = "BCI", 
  "Glucose..mg.dL." = "Glucose [mg/dL]", 
  "PCV...." = "PCV [%]",
  "Lactate..mmol.L." = "Lactate [mmol/L]", 
  "BHB..mmol.L." = "BHB [mmol/L]", 
  "Total.protein..g.dL." = "Total protein [g/dL]",
  "Corticosterone..ng.mL." = "Corticosterone [ng/mL]",
  "WBC.estimate...10..μL." = "WBC estimate [× 10³/μL]", 
  "Mature.heterophils...10..μL." = "Mature heterophils [× 10³/μL]", 
  "Immature.heterophils...10..μL." = "Immature heterophils [× 10³/μL]",
  "Total.heterophils...10..μL." = "Total heterophils [× 10³/μL]",
  "Lymphocytes...10..μL." = "Lymphocytes [× 10³/μL]", 
  "Heterophil.lymphocyte.ratio" = "Heterophil:lymphocyte ratio", 
  "Monocytes...10..μL." = "Monocytes [× 10³/μL]", 
  "Eosinophils...10..μL." = "Eosinophils [× 10³/μL]", 
  "Basophils...10..μL." = "Basophils [× 10³/μL]",
  "Immature.RBC.100.mature.RBC" = "Immature RBC/100 mature RBC"
)

# Apply variable name mapping
names(stragglergreenturtles) <- sapply(names(stragglergreenturtles), function(x) {
  if(x %in% names(variables)) variables[[x]] else x
})

# Convert categorical variables to factors
stragglergreenturtles <- stragglergreenturtles %>%
  mutate(
    Group = factor(Group),
    Sampling.timepoint = factor(Sampling.timepoint)
  )

# Define variable groups
morphometrics_variables <- c("Mass [g]", "SCL [mm]", "SCW [mm]", "BD [mm]", "BCI")

blood_chemistry_variables <- c(
  "Glucose [mg/dL]", "PCV [%]", "Lactate [mmol/L]", 
  "BHB [mmol/L]", "Total protein [g/dL]", "Corticosterone [ng/mL]"
)

hematology_variables <- c(
  "WBC estimate [× 10³/μL]", "Mature heterophils [× 10³/μL]", 
  "Immature heterophils [× 10³/μL]", "Total heterophils [× 10³/μL]", 
  "Lymphocytes [× 10³/μL]", "Heterophil:lymphocyte ratio", 
  "Monocytes [× 10³/μL]", "Eosinophils [× 10³/μL]", 
  "Basophils [× 10³/μL]", "Immature RBC/100 mature RBC"
)

ordered_variables <- c(morphometrics_variables, blood_chemistry_variables, hematology_variables)

# Reference intervals
reference_intervals <- data.frame(
  Variable = c(
    "Mass [g]", "SCL [mm]", "SCW [mm]", "BD [mm]", "BCI",
    "Glucose [mg/dL]", "PCV [%]", "Lactate [mmol/L]", "BHB [mmol/L]", 
    "Total protein [g/dL]", "Corticosterone [ng/mL]",
    "WBC estimate [× 10³/μL]", "Mature heterophils [× 10³/μL]", "Immature heterophils [× 10³/μL]",
    "Total heterophils [× 10³/μL]", "Lymphocytes [× 10³/μL]", "Heterophil:lymphocyte ratio",
    "Monocytes [× 10³/μL]", "Eosinophils [× 10³/μL]", "Basophils [× 10³/μL]", 
    "Immature RBC/100 mature RBC"
  ),
  Lower = c(20.9, 47.8, 35.2, 18.1, 1.67, 39.63, 25, NA, NA, NA, NA,
            5.9, 2.8, 0.07, 3.1, 1.3, 1.19, 0.4, 0, 0.1, 3),
  Upper = c(27.1, 52.7, 39.5, 21.5, 2.15, 140.4, 43, NA, NA, NA, NA,
            17.4, 9.9, 1.7, 10.8, 4.3, 4.26, 2.4, 0.1, 1.4, 24)
)

# Utility function for clean file naming
variable_name_mapping <- list(
  "Mass [g]" = "Mass",
  "SCL [mm]" = "SCL", 
  "SCW [mm]" = "SCW",
  "BD [mm]" = "BD",
  "BCI" = "BCI",
  "Glucose [mg/dL]" = "Glucose",
  "PCV [%]" = "PCV",
  "Lactate [mmol/L]" = "Lactate",
  "BHB [mmol/L]" = "BHB",
  "Total protein [g/dL]" = "Total_Protein",
  "Corticosterone [ng/mL]" = "Corticosterone",
  "WBC estimate [× 10³/μL]" = "WBC_Estimate",
  "Mature heterophils [× 10³/μL]" = "Mature_Heterophils",
  "Immature heterophils [× 10³/μL]" = "Immature_Heterophils", 
  "Total heterophils [× 10³/μL]" = "Total_Heterophils",
  "Lymphocytes [× 10³/μL]" = "Lymphocytes",
  "Heterophil:lymphocyte ratio" = "Heterophil_Lymphocyte_Ratio",
  "Monocytes [× 10³/μL]" = "Monocytes",
  "Eosinophils [× 10³/μL]" = "Eosinophils",
  "Basophils [× 10³/μL]" = "Basophils",
  "Immature RBC/100 mature RBC" = "Immature_RBC"
)

get_clean_filename <- function(variable_name) {
  return(variable_name_mapping[[variable_name]])
}
}