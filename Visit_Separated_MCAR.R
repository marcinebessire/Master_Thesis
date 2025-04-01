#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)

# -------------------------------------------------
# Title: Separated based on Visit
# -------------------------------------------------

#load data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/FAO_data.csv", check.names = FALSE) #34 metabolites

FAO_data$Time_min <- as.numeric(FAO_data$Time_min)

# ------------------------------
# Part 1: Patient 4 add row 120
# ------------------------------

FAO_data_full <- FAO_data

# Ensure the row is added only once
if (!any(FAO_data_full$Patient == "P4" & FAO_data_full$Visit == "Visit 1" & FAO_data_full$Time_min == 120)) {
  
  # Get all rows for Patient 4, Visit 1
  patient4_visit1 <- FAO_data_full %>% filter(Patient == "P4", Visit == "Visit 1")
  
  if (nrow(patient4_visit1) > 0) {
    # Use the first row as a template for metadata
    new_row <- patient4_visit1[1, ]
    
    # Set Time_min to 120
    new_row$Time_min <- 120
    
    # Set all measurement columns (from 6 onward) to NA
    new_row[, 6:ncol(FAO_data_full)] <- NA
    
    # Add the new row to the dataset
    FAO_data_full <- bind_rows(FAO_data_full, new_row)
  }
}


#sort by numeric patient number and then Date
FAO_data_full <- FAO_data_full %>%
  mutate(Patient_num = as.numeric(gsub("P", "", Patient))) %>%
  arrange(Patient_num, Date, Time_min) %>%
  select(-Patient_num)  #remove after sorting

# ------------------------------
# Part 2: MCAR simulation
# ------------------------------

#function to generate missing values per visit dataframe
MCAR_manipulation_visit <- function(data, missing_percentage){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #get the Visit 1 and Visit 2 indices
  visit1_indices <- which(data_copy$Visit == "Visit 1")
  visit2_indices <- which(data_copy$Visit == "Visit 2")
  
  #iterate through each numeric column
  for (col in 6:ncol(data_copy)) {
    #get number of missing values to introduce (rounded down to ensure pairing)
    num_mv_total <- round(nrow(data_copy) * missing_percentage)
    num_mv <- floor(num_mv_total / 2)  #even split between Visit 1 and Visit 2
    
    if (num_mv > 0 && length(visit1_indices) > 0 && length(visit2_indices) > 0) {
      #randomly choose indices for missing values in each visit group
      missing_indices_v1 <- sample(visit1_indices, num_mv, replace = FALSE)
      missing_indices_v2 <- sample(visit2_indices, num_mv, replace = FALSE)
      
      #set selected values to NA
      data_copy[missing_indices_v1, col] <- NA
      data_copy[missing_indices_v2, col] <- NA
    }
  }
  
  return(data_copy)
}
#call function
FAO_10pct_mcar <- MCAR_manipulation_visit(FAO_data_full, 0.1)
FAO_20pct_mcar <- MCAR_manipulation_visit(FAO_data_full, 0.2)
FAO_25pct_mcar <- MCAR_manipulation_visit(FAO_data_full, 0.25)
FAO_30pct_mcar <- MCAR_manipulation_visit(FAO_data_full, 0.3)
FAO_35pct_mcar <- MCAR_manipulation_visit(FAO_data_full, 0.35)
FAO_40pct_mcar <- MCAR_manipulation_visit(FAO_data_full, 0.4)

#split into two dataframe visit 1 and visit 2
#10%
FAO_v1_10pct_mcar <- FAO_10pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_10pct_mcar <- FAO_10pct_mcar %>% filter(Visit == "Visit 2") #60
#20%
FAO_v1_20pct_mcar <- FAO_20pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_20pct_mcar <- FAO_20pct_mcar %>% filter(Visit == "Visit 2") #60
#25%
FAO_v1_25pct_mcar <- FAO_25pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_25pct_mcar <- FAO_25pct_mcar %>% filter(Visit == "Visit 2") #60
#30%
FAO_v1_30pct_mcar <- FAO_30pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_30pct_mcar <- FAO_30pct_mcar %>% filter(Visit == "Visit 2") #60
#35%
FAO_v1_35pct_mcar <- FAO_35pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_35pct_mcar <- FAO_35pct_mcar %>% filter(Visit == "Visit 2") #60
#40%
FAO_v1_40pct_mcar <- FAO_40pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_40pct_mcar <- FAO_40pct_mcar %>% filter(Visit == "Visit 2") #60

