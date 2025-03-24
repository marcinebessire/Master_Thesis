#load necessary libraries 
library(readxl)
library(stringr)
library(dplyr)
library(lubridate) #to handle dates better
library(tidyverse)

#load excel sheet with data 
whole_data <- "/Users/marcinebessire/Desktop/Master_Thesis/clinical_data.xlsx"

#Get each sheet
#BAS data
BAS_data <- read_excel(whole_data, sheet = 2, .name_repair = "minimal")
#FAO data
FAO_data <- read_excel(whole_data, sheet = 3, .name_repair = "minimal")
#intact lipids data
intact_lipids_data <- read_excel(whole_data, sheet = 4, .name_repair = "minimal")

#--------------------------------------------
# Part 1: Clean data set for Sheet 2 and 3 
# -------------------------------------------

#function for first step of cleaning data (remove column and make column names etc)
clean_data <- function(data){
  #remove columns before 3rd columns and remove first and second row
  data_cleaned <- data[, 3:ncol(data)]
  data_cleaned <- data_cleaned[3:nrow(data_cleaned),]
  
  #make second row the column name 
  colnames(data_cleaned) <- as.character(data_cleaned[1, ])  #assign second row as column names
  
  data_cleaned <- data_cleaned[-1, ]  #and remove the row that is now column names
  
  #remove sample ID row
  data_cleaned <- data_cleaned[2:nrow(data_cleaned),]
  
  #remove Zero blank row
  data_cleaned <- data_cleaned[-1,]
  
  #rename column 1 to ID for further processing
  names(data_cleaned)[1]<- paste("ID")
  
  return(data_cleaned)
  
}

#call function to clean data
BAS_data_cleaned <- clean_data(BAS_data)
FAO_data_cleaned <- clean_data(FAO_data)


#clean and make new columns
process_patient_data <- function(df) {
  df_clean <- df %>%
    #remove rows where first column starts with "NIST_"
    filter(!str_detect(.[[1]], "^.*NIST_")) %>%
  
    #extract patient info from first column
    mutate(
      ID = .[[1]],  #keep original identifier
      
      #extract patient number e.g. "P10"
      Patient = str_extract(ID, "P\\d+"),
      
      #extract numeric part of patient
      Patient_num = as.numeric(str_remove(Patient, "P")),
      
      #extract date in format DDMMYY (e.g. 250222)
      Date_raw = str_extract(ID, "\\d{6}"),
      Date = dmy(Date_raw),  #convert
      
      #extract time in minutes (after E)
      Time_min = as.numeric(str_extract(ID, "(?<=E)\\d+")),
      
      #format for output
      Date_formatted = format(Date, "%m/%d/%Y")
    ) %>%
    
    #sort by patient and full date
    arrange(Patient_num, Date, Time_min) %>%
    select(ID, Patient, Date, Time_min, everything(), -Date_raw, -Date_formatted, -Patient_num)
  
  return(df_clean)
}

#call function 
BAS_data_extended <- process_patient_data(BAS_data_cleaned)
FAO_data_extended <- process_patient_data(FAO_data_cleaned)

#function to convert data to numeric and remove whole columns with all same vlaue
convert_columns_to_numeric <- function(data) {
  data[, 4:ncol(data)] <- lapply(data[, 4:ncol(data)], function(x) {
    suppressWarnings(as.numeric(x)) #convert to numeric and replace non-numeric values with NA
  })
  
  #remove columns with dame values 
  data <- data[, sapply(data, function(col) length(unique(na.omit(col))) > 1)]
  return(data)
}


#apply function to convert to numeric
BAS_data_final <- convert_columns_to_numeric(BAS_data_extended)
FAO_data_final <- convert_columns_to_numeric(FAO_data_extended)

#save to csv file 
write_csv(BAS_data_final, "/Users/marcinebessire/Desktop/Master_Thesis/BAS_data.csv")
write_csv(FAO_data_final, "/Users/marcinebessire/Desktop/Master_Thesis/FAO_data.csv")


#--------------------------------------------
# Part 2: Clean dataset for Sheet 4
# -------------------------------------------

#make copy 
intact_lipids_data_copy <- intact_lipids_data

#define the values to ignore
ignore_values <- c("Species", "NA", "Trial NA", "MS1") 
new_colnames <- names(intact_lipids_data_copy)

#loop through each column index
for (i in seq_along(names(intact_lipids_data_copy))) {
  first_row_val <- as.character(intact_lipids_data_copy[1, i])  #get first row value
  column_name <- names(intact_lipids_data_copy)[i]  #get column name
  
  #if the value is not in the ignore list, merge it with column name
  if (!(is.na(first_row_val) || first_row_val %in% ignore_values)) {
    new_colnames[i] <- paste0(column_name, "_", first_row_val)
  }
}

#assign the new names
colnames(intact_lipids_data_copy) <- new_colnames

#remove first row after merging
intact_lipids_data_copy <- intact_lipids_data_copy[-1,]

#rename column 1 to ID for further processing
names(intact_lipids_data_copy)[1]<- paste("ID")

#add unique names 
colnames(intact_lipids_data_copy) <- make.names(colnames(intact_lipids_data_copy), unique = TRUE)

#call function
intact_lipids_cleaned <- process_patient_data(intact_lipids_data_copy)

intact_lipids_ordered <- intact_lipids_cleaned %>%
  group_by(Patient) %>%
  mutate(Visit = ifelse(row_number() == 1, "Visit 1", "Visit 2")) %>%
  select(ID, Patient, Date, Time_min, Visit, everything())

#now convert to n umeric and renomve columns with same values
convert_columns_to_numeric_lipids <- function(data) {
  #take numeric data
  numeric_data <- data[, 6:ncol(data)]
  
  #keep only columns that more than 1 different value
  cols_to_keep <- sapply(numeric_data, function(col) length(unique(na.omit(col))) > 1)
  numeric_data <- numeric_data[, cols_to_keep]
  
  #convert columns to numeric
  numeric_data <- lapply(numeric_data, function(x) suppressWarnings(as.numeric(x)))
  numeric_data <- as.data.frame(numeric_data)

  
  #combine with metadata
  data_cleaned <- cbind(data[, 1:5], numeric_data)
  
  return(data_cleaned)
}

#call function
intact_lipids_final <- convert_columns_to_numeric_lipids(intact_lipids_ordered)


