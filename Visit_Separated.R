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

# -------------------------------------------------
# Part 1: MCAR simulation
# -------------------------------------------------

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
FAO_10pct_mcar <- MCAR_manipulation_visit(FAO_data, 0.1)
FAO_20pct_mcar <- MCAR_manipulation_visit(FAO_data, 0.2)
FAO_25pct_mcar <- MCAR_manipulation_visit(FAO_data, 0.25)
FAO_30pct_mcar <- MCAR_manipulation_visit(FAO_data, 0.3)
FAO_35pct_mcar <- MCAR_manipulation_visit(FAO_data, 0.35)
FAO_40pct_mcar <- MCAR_manipulation_visit(FAO_data, 0.4)

#split into two dataframe visit 1 and visit 2
#10%
FAO_v1_10pct_mcar <- FAO_10pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_10pct_mcar <- FAO_10pct_mcar %>% filter(Visit == "Visit 2") #59
#20%
FAO_v1_20pct_mcar <- FAO_20pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_20pct_mcar <- FAO_20pct_mcar %>% filter(Visit == "Visit 2") #59
#25%
FAO_v1_25pct_mcar <- FAO_25pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_25pct_mcar <- FAO_25pct_mcar %>% filter(Visit == "Visit 2") #59
#30%
FAO_v1_30pct_mcar <- FAO_30pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_30pct_mcar <- FAO_30pct_mcar %>% filter(Visit == "Visit 2") #59
#35%
FAO_v1_35pct_mcar <- FAO_35pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_35pct_mcar <- FAO_35pct_mcar %>% filter(Visit == "Visit 2") #59
#40%
FAO_v1_40pct_mcar <- FAO_40pct_mcar %>% filter(Visit == "Visit 1") #60
FAO_v2_40pct_mcar <- FAO_40pct_mcar %>% filter(Visit == "Visit 2") #59

# -------------------------------------------------
# Part 2.1: MNAR simulation
# -------------------------------------------------

#function to generate missing values per visit dataframe
MNAR_manipulation_visit <- function(data, missing_percentage){
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
    
    if (num_mv > 0 && length(visit1_indices) >= num_mv && length(visit2_indices) >= num_mv) {
      #column values for each visit
      col_v1 <- data_copy[visit1_indices, col]
      col_v2 <- data_copy[visit2_indices, col]
      
      #indices with lowest value between each group
      lowest_v1_indices <- visit1_indices[order(col_v1, na.last = NA)][1:num_mv]
      lowest_v2_indices <- visit2_indices[order(col_v2, na.last = NA)][1:num_mv]
      
      #set to NA
      data_copy[lowest_v1_indices, col] <- NA
      data_copy[lowest_v2_indices, col] <- NA
     }
  }
  
  return(data_copy)
}

#call function
FAO_10pct_mnar <- MNAR_manipulation_visit(FAO_data, 0.1)
FAO_20pct_mnar <- MNAR_manipulation_visit(FAO_data, 0.2)
FAO_25pct_mnar <- MNAR_manipulation_visit(FAO_data, 0.25)
FAO_30pct_mnar <- MNAR_manipulation_visit(FAO_data, 0.3)
FAO_35pct_mnar <- MNAR_manipulation_visit(FAO_data, 0.35)
FAO_40pct_mnar <- MNAR_manipulation_visit(FAO_data, 0.4)

#split into two dataframe visit 1 and visit 2
#10%
FAO_v1_10pct_mnar <- FAO_10pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_10pct_mnar <- FAO_10pct_mnar %>% filter(Visit == "Visit 2") #59
#20%
FAO_v1_20pct_mnar <- FAO_20pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_20pct_mnar <- FAO_20pct_mnar %>% filter(Visit == "Visit 2") #59
#25%
FAO_v1_25pct_mnar <- FAO_25pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_25pct_mnar <- FAO_25pct_mnar %>% filter(Visit == "Visit 2") #59
#30%
FAO_v1_30pct_mnar <- FAO_30pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_30pct_mnar <- FAO_30pct_mnar %>% filter(Visit == "Visit 2") #59
#35%
FAO_v1_35pct_mnar <- FAO_35pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_35pct_mnar <- FAO_35pct_mnar %>% filter(Visit == "Visit 2") #59
#40%
FAO_v1_40pct_mnar <- FAO_40pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_40pct_mnar <- FAO_40pct_mnar %>% filter(Visit == "Visit 2") #59
