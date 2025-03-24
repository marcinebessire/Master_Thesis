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
  
  return(data_cleaned)
  
}

#call function to clean data
BAS_data_cleaned <- clean_data(BAS_data)
FAO_data_cleaned <- clean_data(FAO_data)


#function to extract participant number "PX", date, time and visit 1 or visit 2
make_new_columns <- function(data){
  data <- data %>%
    #filter out rows starting with NIST_
    filter(
      !(str_detect(as.character(data[[1]]), "^NIST_"))
    ) %>%
    mutate(
      Participant = str_extract(as.character(data[[1]]), "P[0-9]+"), #extract "P" followed by digits
      Date = str_extract(data[[1]], "[0-9]{6}"), #6 digit date
      Time = str_extract(data[[1]], "E[0-9]+") %>% str_replace("E", "") #extract number after E to get time
    ) %>%
    mutate(
      Day = str_sub(Date,1,2), #extract first two digits for day
      Month = str_sub(Date,3,4), #extract next two digits for month
      Year = paste0("20", str_sub(Date,5,6)), #extract year in format e.g. 2020
      Full_Date = as.Date(paste0(Year, "-", Month, "-", Day), format = "%Y-%m-%d"),
      MonthDay = format(Full_Date, "%m/%d") #format for visualization
    ) %>%
    mutate(
      MonthDaycalc = as.numeric(format(Full_Date, "%m%d")), #numeric sorting
      Participant_Number = as.numeric(str_extract(Participant, "[0-9]+")) #extract numeric part of Participant for sorting
    ) %>%
    arrange(Participant_Number, Full_Date) %>% #sort by date and participant
    group_by(Participant) %>%
    select(ID,Participant, MonthDay, Year, Time, Visit, everything()) %>%
    select(-Date, -Day, -Month, -MonthDaycalc, -Participant_Number, -Full_Date) %>%
    
    return(data)
}

#call function 
BAS_data_extended <- make_new_columns(BAS_data_cleaned)
FAO_data_extended <- make_new_columns(FAO_data_cleaned)

#function to convert data to numeric and remove whole columns with all same vlaue
convert_columns_to_numeric <- function(data) {
  data[, 7:ncol(data)] <- lapply(data[, 7:ncol(data)], function(x) {
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

intact_lipids_data_cleaned <- intact_lipids_data_copy %>%
  mutate( #mutate to create or edit existing columns in a dataframe
    #extract date from anywhere in the Name column (because not the same)
    Whole_Date = str_extract(Name, "\\d{8}"), #\\d for matching any digits
    Whole_Date = as.Date(Whole_Date, format="%Y%m%d"), 
    Year = year(Whole_Date), 
    MonthDay = format(Whole_Date, "%m-%d"),
    
    #extract only the meaningful part of ID 
    ID = str_remove(Name, ".*?_NIST_?"),  #remove everything before and including "NIST_"
    ID = str_remove(ID, "^DI_"),  #remove "DI_" prefix if present, ^ for the beginning 
    ID = str_remove(ID, "\\d{8}"),  #remove date if it appears again (for the wrongly ordered)
    ID = str_replace(ID, "^_|_$", ""),  #remove underscores
    
    #assign trial numbers based on MonthDay
    Trial_number = dense_rank(MonthDay),  
    Trial = paste0("Trial ", Trial_number)  
  ) %>%
  select(Name, ID, Year, MonthDay, Trial, everything(),-Whole_Date, -Trial_number)


# Part 4 ------
# Convert data columns (keeping metadata unchanged)
info_cols <- c("Name", "ID", "Year", "MonthDay", "Trial") 

df_data_cleaned24 <- df_data_cleaned24 %>%
  mutate(across(-all_of(info_cols), ~ {
    numeric_col <- suppressWarnings(as.numeric(.))
    ifelse(is.infinite(numeric_col), NA, numeric_col) #put all Inf values to NA 
  }))

print(df_data_cleaned24)

df_data_cleaned24 <- df_data_cleaned24 %>%
  arrange(MonthDay) # Orders the data by Trial number

#write new CSV
write_csv(df_data_cleaned24, "/Users/marcinebessire/Desktop/project/cleaned_data_2024.csv")

# Part 5 ----
# *optional* remove all the "new" ones (new chip used) started from 09-10 (aka starting from Trial 6)
cut_df_data_cleaned24 <- df_data_cleaned24 %>%
  filter(MonthDay < "09-10") 

#save the cut dataset
write_csv(cut_df_data_cleaned24, "/Users/marcinebessire/Desktop/project/cut_cleaned_data_2024.csv")
