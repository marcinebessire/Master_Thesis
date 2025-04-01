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
# Part 2: MNAR simulation
# ------------------------------

#function to simulate MNAR
MNAR_manipulation_visit <- function(data, missing_percentage){
  data_copy <- data
  
  #get Visit indices
  visit1_indices <- which(data_copy$Visit == "Visit 1")
  visit2_indices <- which(data_copy$Visit == "Visit 2")
  
  #for each measurement column (assuming starts at col 6)
  for (col in 6:ncol(data_copy)) {
    
    #visit 1
    if (length(visit1_indices) > 0) {
      values_v1 <- data_copy[visit1_indices, col]
      num_mv_v1 <- floor(length(values_v1) * missing_percentage)
      
      #get indices of the lowest values in Visit 1
      if (num_mv_v1 > 0) {
        lowest_v1 <- order(values_v1, na.last = NA)[1:num_mv_v1]
        rows_to_na_v1 <- visit1_indices[lowest_v1]
        data_copy[rows_to_na_v1, col] <- NA
      }
    }
    
    #visit 2
    if (length(visit2_indices) > 0) {
      values_v2 <- data_copy[visit2_indices, col]
      num_mv_v2 <- floor(length(values_v2) * missing_percentage)
      
      #get indices of the lowest values in Visit 2
      if (num_mv_v2 > 0) {
        lowest_v2 <- order(values_v2, na.last = NA)[1:num_mv_v2]
        rows_to_na_v2 <- visit2_indices[lowest_v2]
        data_copy[rows_to_na_v2, col] <- NA
      }
    }
  }
  
  return(data_copy)
}


#call function
FAO_10pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.10)
FAO_20pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.20)
FAO_25pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.25)
FAO_30pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.30)
FAO_35pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.35)
FAO_40pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.40)


#split into two dataframe visit 1 and visit 2
#10%
FAO_v1_10pct_mnar <- FAO_10pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_10pct_mnar <- FAO_10pct_mnar %>% filter(Visit == "Visit 2") #60
#20%
FAO_v1_20pct_mnar <- FAO_20pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_20pct_mnar <- FAO_20pct_mnar %>% filter(Visit == "Visit 2") #60
#25%
FAO_v1_25pct_mnar <- FAO_25pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_25pct_mnar <- FAO_25pct_mnar %>% filter(Visit == "Visit 2") #60
#30%
FAO_v1_30pct_mnar <- FAO_30pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_30pct_mnar <- FAO_30pct_mnar %>% filter(Visit == "Visit 2") #60
#35%
FAO_v1_35pct_mnar <- FAO_35pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_35pct_mnar <- FAO_35pct_mnar %>% filter(Visit == "Visit 2") #60
#40%
FAO_v1_40pct_mnar <- FAO_40pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_40pct_mnar <- FAO_40pct_mnar %>% filter(Visit == "Visit 2") #60

