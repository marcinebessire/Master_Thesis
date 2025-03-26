#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(openxlsx)
library(readxl)
library(zoo) #for interpolation
library(purrr)


#load whole data
mcar_data <- "/Users/marcinebessire/Desktop/Master_Thesis/FAO_MCAR_patient_visit_sep.xlsx"
original_data <- "/Users/marcinebessire/Desktop/Master_Thesis/FAO_original_patient_visit_sep.xlsx"

# #open each sheet in excel
# open_excel_sheets <- function(file_path){
#   #get sheet name
#   sheet_names <- excel_sheets(file_path)
# 
#   for (sheet in sheet_names){
#     df <- read_excel(file_path, sheet = sheet)
# 
#     #create R variable name
#     var_name <- make.names(sheet)
# 
#     #assing the dataframe to variable
#     assign(var_name, df, envir = .GlobalEnv)
#   }
# }
# 
# #call function to open all sheets
# open_excel_sheets(mcar_data)
# open_excel_sheets(original_data)

# ----------------------------------------------------
# Part 1: Area under the curve (AUC) and interpolation
# ----------------------------------------------------

#interpolation function
interpolate_mcar <- function(df){
  time <- df$Time_min
  
  #metaa data
  metadata <- c("ID", "Patient", "Date", "Time_min", "Visit")
  
  metabolite_cols <- setdiff(names(df), metadata)
  
  #interpolate each metabolite column
  df[metabolite_cols] <- lapply(df[metabolite_cols], function(col) {
    na.approx(col, x = time, na.rm = FALSE) #linear interpolation, values filled based on time spacing
  })
  
  return(df)
}

#trapezoidal AUC 
compute_auc <- function(df){
  time <- df$Time_min
  metadata <- c("ID", "Patient", "Date", "Time_min", "Visit")
  
  metabolite_cols <- setdiff(names(df), metadata)
  
  auc <- sapply(metabolite_cols, function(metabolite) {
    conc <- df[[metabolite]]
    valid <- !is.na(time) & !is.na(conc)
    
    if (sum(valid) < 2) return(NA) #if not enough data
    
    sum(diff(time[valid]) * (head(conc[valid], -1) + tail(conc[valid], -1)) / 2)
  })
  
  return(auc)
}

#function to process all sheets
process_all_sheets <- function(dataset, apply_interpolation = TRUE) {
  sheet_names <- excel_sheets(dataset)
  
  #empty list
  auc_list <- list()
  
  for (sheet in sheet_names){
    df <- read_excel(dataset, sheet = sheet)
    
    #only interpolate if there are missing values 
    if (apply_interpolation) {
      df <- interpolate_mcar(df)
    }
    
    auc_values <- compute_auc(df)
    auc_list[[sheet]] <- auc_values
  }
  
  return(auc_list)
} 

#call function to process sheets
auc_mcar <- process_all_sheets(mcar_data, apply_interpolation = TRUE)
auc_original <- process_all_sheets(original_data, apply_interpolation = FALSE)

#compare now the auc results
compare_auc <- function(auc_list1, auc_list2){
  sheet_names <- intersect(names(auc_list1), names(auc_list2))
  
  comparison_df <- bind_rows(
    lapply(sheet_names, function(sheet) {
      auc1 <- auc_list1[[sheet]]
      auc2 <- auc_list2[[sheet]]
      
      common_metabolites <- intersect(names(auc1), names(auc2))
      
      auc_mcar_vals <- unname(auc1[common_metabolites])
      auc_orig_vals <- unname(auc2[common_metabolites])
      
      #compute difference
      diff_vals <- auc_mcar_vals - auc_orig_vals
      
      #compute relative difference s(afely )handle division by 0
      rel_diff <- ifelse(
        is.na(auc_orig_vals) | auc_orig_vals == 0,
        NA,
        (diff_vals / auc_orig_vals) * 100
      )
    
      
      data.frame(
        Sheet = sheet,
        Metabolite = common_metabolites,
        AUC_MCAR = auc_mcar_vals,
        AUC_Original = auc_orig_vals,
        Absolute_Diff = diff_vals,
        Relative_Diff_Percent = rel_diff
      )
    })
  )
  
  return(comparison_df)
}

#call comparison function
auc_comparison <- compare_auc(auc_mcar, auc_original)

#Make density plot
prepare_density_data <- function(file_path_mcar, file_path_orig) {
  sheets <- intersect(excel_sheets(file_path_mcar), excel_sheets(file_path_orig))
  
  all_data <- purrr::map_dfr(sheets, function(sheet) {
    df_mcar <- read_excel(file_path_mcar, sheet = sheet)
    df_orig <- read_excel(file_path_orig, sheet = sheet)
    
    #interpolate MCAR only if it has NAs
    if (anyNA(df_mcar)) df_mcar <- interpolate_mcar(df_mcar)
    
    #define metadata
    metadata <- c("ID", "Patient", "Date", "Time_min", "Visit")
    metabolite_cols <- setdiff(names(df_mcar), metadata)
    
    #add sheet-based Patient and Visit info
    patient <- str_extract(sheet, "p\\d+")
    visit <- str_extract(sheet, "v\\d+")
    
    df_mcar_long <- df_mcar %>%
      select(Time_min, all_of(metabolite_cols)) %>%
      pivot_longer(-Time_min, names_to = "Metabolite", values_to = "Value") %>%
      mutate(Source = "MCAR", Patient = patient, Visit = visit)
    
    df_orig_long <- df_orig %>%
      select(Time_min, all_of(metabolite_cols)) %>%
      pivot_longer(-Time_min, names_to = "Metabolite", values_to = "Value") %>%
      mutate(Source = "Original", Patient = patient, Visit = visit)
    
    bind_rows(df_mcar_long, df_orig_long)
  })
  
  return(all_data)
}

#call function to prepare for plotting
density_data <- prepare_density_data(mcar_data, original_data)

#make density plot (whole FAO data for each patient)
ggplot(density_data, aes(x = Value, fill = Source, color = Source)) +
  geom_density(alpha = 0.4, adjust = 1.2, na.rm = TRUE) +
  facet_grid(Patient ~ Visit, scales = "free") +
  labs(
    title = "Density Distribution of Metabolite Values: MCAR vs Original",
    x = "Metabolite Value",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  ) +
  xlim(-50,200)

#plot Density for each metabolite and each patient and visit separately 

#output PDF path
pdf_path <- "/Users/marcinebessire/Desktop/Master_Thesis/Density_original_vs_mcar_all.pdf"

#open pdf 
pdf(pdf_path, width = 12, height = 6)

#get unique pateint visit 
combo_list <- density_data %>%
  distinct(Patient, Visit)

#loop over each Patient + Visit combo and plot
for (i in 1:nrow(combo_list)) {
  patient <- combo_list$Patient[i]
  visit <- combo_list$Visit[i]
  
  df_subset <- density_data %>%
    filter(Patient == patient, Visit == visit)
  
  if (nrow(df_subset) == 0) next  # skip if empty
  
  p <- ggplot(df_subset, aes(x = Value, fill = Source, color = Source)) +
    geom_density(alpha = 0.4, adjust = 1.2, na.rm = TRUE) +
    facet_wrap(~ Metabolite, scales = "free") +
    labs(
      title = paste("Patient:", patient, "| Visit:", visit),
      x = "Value",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      strip.text = element_text(size = 9)
    )
  
  print(p)  
}

#close the PDF device
dev.off()
