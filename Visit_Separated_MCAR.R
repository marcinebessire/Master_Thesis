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
#original
FAO_original_v1 <- FAO_data %>% filter(Visit == "Visit 1") #59
FAO_original_v2 <- FAO_data %>% filter(Visit == "Visit 2") #60
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

# --------------------------------------
# TITLE: IMPUTATION METHODS
# --------------------------------------
# --------------------------------------
# Part 1: Linear interpolation
# --------------------------------------

#interpolate missing data 
interpolate_missing <- function(data){
  data_copy <- data
  
  #apply interpolation to each numeric column from 6th to last
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[, 6:ncol(data_copy)], function(col){
    if (is.numeric(col)) {
      return(na.approx(col, na.rm = FALSE, rule = 2)) #uses the first non-NA value as a flat extrapolation (same for)
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}

#call interpolation function on mcar data
#Visit1
v1_10pct_mcar_interpolation <- interpolate_missing(FAO_v1_10pct_mcar)
v1_20pct_mcar_interpolation <- interpolate_missing(FAO_v1_20pct_mcar)
v1_25pct_mcar_interpolation <- interpolate_missing(FAO_v1_25pct_mcar)
v1_30pct_mcar_interpolation <- interpolate_missing(FAO_v1_30pct_mcar)
v1_35pct_mcar_interpolation <- interpolate_missing(FAO_v1_35pct_mcar)
v1_40pct_mcar_interpolation <- interpolate_missing(FAO_v1_40pct_mcar)
#Visit2
v2_10pct_mcar_interpolation <- interpolate_missing(FAO_v2_10pct_mcar)
v2_20pct_mcar_interpolation <- interpolate_missing(FAO_v2_20pct_mcar)
v2_25pct_mcar_interpolation <- interpolate_missing(FAO_v2_25pct_mcar)
v2_30pct_mcar_interpolation <- interpolate_missing(FAO_v2_30pct_mcar)
v2_35pct_mcar_interpolation <- interpolate_missing(FAO_v2_35pct_mcar)
v2_40pct_mcar_interpolation <- interpolate_missing(FAO_v2_40pct_mcar)

# --------------------------------------
# Part 2: Kalman Smoothing
# --------------------------------------

#create function to impute MV with Kalman Smoothing
kalman_imputation <- function(data){
  data_copy <- data
  
  #apply kalman smoothing to each numeric column
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[,6:ncol(data_copy)], function(col){
    if (is.numeric(col)){
      return(na_kalman(col, model = "StructTS", smooth = TRUE))
    } else {
      return(col)
    }
  })
  return(data_copy)
}

#call kalman function
#Visit 1
v1_10pct_mcar_kalman <- kalman_imputation(FAO_v1_10pct_mcar)
v1_20pct_mcar_kalman <- kalman_imputation(FAO_v1_20pct_mcar)
v1_25pct_mcar_kalman <- kalman_imputation(FAO_v1_25pct_mcar)
v1_30pct_mcar_kalman <- kalman_imputation(FAO_v1_30pct_mcar)
v1_35pct_mcar_kalman <- kalman_imputation(FAO_v1_35pct_mcar)
v1_40pct_mcar_kalman <- kalman_imputation(FAO_v1_40pct_mcar)
#Visit 2
v2_10pct_mcar_kalman <- kalman_imputation(FAO_v2_10pct_mcar)
v2_20pct_mcar_kalman <- kalman_imputation(FAO_v2_20pct_mcar)
v2_25pct_mcar_kalman <- kalman_imputation(FAO_v2_25pct_mcar)
v2_30pct_mcar_kalman <- kalman_imputation(FAO_v2_30pct_mcar)
v2_35pct_mcar_kalman <- kalman_imputation(FAO_v2_35pct_mcar)
v2_40pct_mcar_kalman <- kalman_imputation(FAO_v2_40pct_mcar)


# -----------------------------------------------------------
# Part 3: Linear Weighted Moving Average (LWMA)
# ----------------------------------------------------------

#function to make weighted moving average imputation
weighted_mov_average <- function(data, window = 3){ 
  data_copy <- data
  
  data_copy[,6:ncol(data_copy)] <- lapply(data_copy[,6:ncol(data_copy)], function(col){
    if (is.numeric(col)){
      return(na_ma(col, k = window, weighting = "linear")) #give more value to the closer points and linearly decreasing
    } else {
      return(col)
    }
  })
  return(data_copy)
}

#call LWMA
#Visit 1
v1_10pct_mcar_lwma <- weighted_mov_average(FAO_v1_10pct_mcar)
v1_20pct_mcar_lwma <- weighted_mov_average(FAO_v1_20pct_mcar)
v1_25pct_mcar_lwma <- weighted_mov_average(FAO_v1_25pct_mcar)
v1_30pct_mcar_lwma <- weighted_mov_average(FAO_v1_30pct_mcar)
v1_35pct_mcar_lwma <- weighted_mov_average(FAO_v1_35pct_mcar)
v1_40pct_mcar_lwma <- weighted_mov_average(FAO_v1_40pct_mcar)
#Visit 2
v2_10pct_mcar_lwma <- weighted_mov_average(FAO_v2_10pct_mcar)
v2_20pct_mcar_lwma <- weighted_mov_average(FAO_v2_20pct_mcar)
v2_25pct_mcar_lwma <- weighted_mov_average(FAO_v2_25pct_mcar)
v2_30pct_mcar_lwma <- weighted_mov_average(FAO_v2_30pct_mcar)
v2_35pct_mcar_lwma <- weighted_mov_average(FAO_v2_35pct_mcar)
v2_40pct_mcar_lwma <- weighted_mov_average(FAO_v2_40pct_mcar)


# --------------------------------------
# TITLE: KINETICS PLOT BEFORE AND AFTER
# --------------------------------------
# -------------------------------
# Plot before and After 
# -------------------------------
  
#plot dots and line

plot_imputed_vs_original <- function(data, visit, percent, file_path){
  # Ensure correct type
  data <- data %>%
    mutate(
      Visit = as.factor(Visit),
      Time_min = as.numeric(Time_min),
      Method = as.factor(Method)
    )
  
  # Reshape to long format
  long_df <- data %>%
    pivot_longer(cols = 6:(ncol(data) - 2), # last two cols: Method, MissingPct
                 names_to = "Metabolite",
                 values_to = "Concentration")
  
  # Unique patients
  patients <- unique(long_df$Patient)
  
  # Create PDF
  pdf(file_path, width = 16, height = 10)
  
  for (patient in patients) {
    plot_data <- long_df %>% filter(Patient == patient)
    
    p <- ggplot(plot_data, aes(x = Time_min, y = Concentration, color = Method)) +
      geom_point(size = 2) +
      geom_line(aes(group = Method), linewidth = 1) +
      facet_wrap(~ Metabolite, scales = "free_y", ncol = 6) +
      labs(
        title = paste(patient, ":", "Kinetics of all Metabolites (", visit, ", ", percent, "%)"),
        subtitle = "Original vs Imputed Methods",
        x = "Time [min]",
        y = "Concentration [µM]"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14)
      )
    
    print(p)
  }
  
  dev.off()
}


# ----------------------------
# Part 1: Combine Visit 1
# ----------------------------

#10pct missingness
v1_10pct_all <- bind_rows(
  FAO_v1_10pct_mcar %>% mutate(MissingPct = 10, Method = "Original"),
  v1_10pct_mcar_interpolation %>% mutate(MissingPct = 10, Method = "Interpolation"),
  v1_10pct_mcar_kalman %>% mutate(MissingPct = 10, Method = "Kalman"),
  v1_10pct_mcar_lwma %>% mutate(MissingPct = 10, Method = "LWMA")
)

#20pct missingness
v1_20pct_all <- bind_rows(
  FAO_v1_20pct_mcar %>% mutate(MissingPct = 20, Method = "Original"),
  v1_20pct_mcar_interpolation %>% mutate(MissingPct = 20, Method = "Interpolation"),
  v1_20pct_mcar_kalman %>% mutate(MissingPct = 20, Method = "Kalman"),
  v1_20pct_mcar_lwma %>% mutate(MissingPct = 20, Method = "LWMA")
)

#25pct missingness
v1_25pct_all <- bind_rows(
  FAO_v1_25pct_mcar %>% mutate(MissingPct = 25, Method = "Original"),
  v1_25pct_mcar_interpolation %>% mutate(MissingPct = 25, Method = "Interpolation"),
  v1_25pct_mcar_kalman %>% mutate(MissingPct = 25, Method = "Kalman"),
  v1_25pct_mcar_lwma %>% mutate(MissingPct = 25, Method = "LWMA")
)

#30pct missingness
v1_30pct_all <- bind_rows(
  FAO_v1_30pct_mcar %>% mutate(MissingPct = 30, Method = "Original"),
  v1_30pct_mcar_interpolation %>% mutate(MissingPct = 30, Method = "Interpolation"),
  v1_30pct_mcar_kalman %>% mutate(MissingPct = 30, Method = "Kalman"),
  v1_30pct_mcar_lwma %>% mutate(MissingPct = 30, Method = "LWMA")
)

#35pct missingness
v1_35pct_all <- bind_rows(
  FAO_v1_35pct_mcar %>% mutate(MissingPct = 35, Method = "Original"),
  v1_35pct_mcar_interpolation %>% mutate(MissingPct = 35, Method = "Interpolation"),
  v1_35pct_mcar_kalman %>% mutate(MissingPct = 35, Method = "Kalman"),
  v1_35pct_mcar_lwma %>% mutate(MissingPct = 35, Method = "LWMA")
)

#40pct missingness
v1_40pct_all <- bind_rows(
  FAO_v1_40pct_mcar %>% mutate(MissingPct = 40, Method = "Original"),
  v1_40pct_mcar_interpolation %>% mutate(MissingPct = 40, Method = "Interpolation"),
  v1_40pct_mcar_kalman %>% mutate(MissingPct = 40, Method = "Kalman"),
  v1_40pct_mcar_lwma %>% mutate(MissingPct = 40, Method = "LWMA")
)


#named list of all datasets
all_datasets_v1 <- list(
  "10" = v1_10pct_all,
  "20" = v1_20pct_all,
  "25" = v1_25pct_all,
  "30" = v1_30pct_all,
  "35" = v1_35pct_all,
  "40" = v1_40pct_all
)

#loop through and plot each
for (pct in names(all_datasets_v1)) {
  file_path <- paste0("/Users/marcinebessire/Desktop/Master_Thesis/Visit_Separated/Visit1_", pct, "pct_imputation_plots.pdf")
  plot_imputed_vs_original(
    data = all_datasets_v1[[pct]],
    visit = "Visit1",
    percent = pct,
    file_path = file_path
  )
}

# ----------------------------
# Part 2: Combine Visit 2
# ----------------------------

#10pct missingness
v2_10pct_all <- bind_rows(
  FAO_v2_10pct_mcar %>% mutate(MissingPct = 10, Method = "Original"),
  v2_10pct_mcar_interpolation %>% mutate(MissingPct = 10, Method = "Interpolation"),
  v2_10pct_mcar_kalman %>% mutate(MissingPct = 10, Method = "Kalman"),
  v2_10pct_mcar_lwma %>% mutate(MissingPct = 10, Method = "LWMA")
)

#20pct missingness
v2_20pct_all <- bind_rows(
  FAO_v2_20pct_mcar %>% mutate(MissingPct = 20, Method = "Original"),
  v2_20pct_mcar_interpolation %>% mutate(MissingPct = 20, Method = "Interpolation"),
  v2_20pct_mcar_kalman %>% mutate(MissingPct = 20, Method = "Kalman"),
  v2_20pct_mcar_lwma %>% mutate(MissingPct = 20, Method = "LWMA")
)

#25pct missingness
v2_25pct_all <- bind_rows(
  FAO_v2_25pct_mcar %>% mutate(MissingPct = 25, Method = "Original"),
  v2_25pct_mcar_interpolation %>% mutate(MissingPct = 25, Method = "Interpolation"),
  v2_25pct_mcar_kalman %>% mutate(MissingPct = 25, Method = "Kalman"),
  v2_25pct_mcar_lwma %>% mutate(MissingPct = 25, Method = "LWMA")
)

#30pct missingness
v2_30pct_all <- bind_rows(
  FAO_v2_30pct_mcar %>% mutate(MissingPct = 30, Method = "Original"),
  v2_30pct_mcar_interpolation %>% mutate(MissingPct = 30, Method = "Interpolation"),
  v2_30pct_mcar_kalman %>% mutate(MissingPct = 30, Method = "Kalman"),
  v2_30pct_mcar_lwma %>% mutate(MissingPct = 30, Method = "LWMA")
)

#35pct missingness
v2_35pct_all <- bind_rows(
  FAO_v2_35pct_mcar %>% mutate(MissingPct = 35, Method = "Original"),
  v2_35pct_mcar_interpolation %>% mutate(MissingPct = 35, Method = "Interpolation"),
  v2_35pct_mcar_kalman %>% mutate(MissingPct = 35, Method = "Kalman"),
  v2_35pct_mcar_lwma %>% mutate(MissingPct = 35, Method = "LWMA")
)

#40pct missingness
v2_40pct_all <- bind_rows(
  FAO_v2_40pct_mcar %>% mutate(MissingPct = 40, Method = "Original"),
  v2_40pct_mcar_interpolation %>% mutate(MissingPct = 40, Method = "Interpolation"),
  v2_40pct_mcar_kalman %>% mutate(MissingPct = 40, Method = "Kalman"),
  v2_40pct_mcar_lwma %>% mutate(MissingPct = 40, Method = "LWMA")
)


#named list of all datasets
all_datasets_v2 <- list(
  "10" = v2_10pct_all,
  "20" = v2_20pct_all,
  "25" = v2_25pct_all,
  "30" = v2_30pct_all,
  "35" = v2_35pct_all,
  "40" = v2_40pct_all
)

#loop through and plot each
for (pct in names(all_datasets_v2)) {
  file_path <- paste0("/Users/marcinebessire/Desktop/Master_Thesis/Visit_Separated/Visit2_", pct, "pct_imputation_plots.pdf")
  plot_imputed_vs_original(
    data = all_datasets_v2[[pct]],
    visit = "Visit1",
    percent = pct,
    file_path = file_path
  )
}



