#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(zoo) #for interpolation
library(imputeTS) #for imputation methods
library(pracma) #for AUC calculation

# --------------------------------------
# TITLE: PATIENT AND VISIT SEPARATED
# --------------------------------------

#load data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/FAO_data.csv", check.names = FALSE) #34 metabolites

# --------------------------------------
# TITLE: MISSING VALUE SIMULATION
# --------------------------------------
# ---------------------------------------------------
# Part 1: Split Dataframe according to Patient and Visit
# ----------------------------------------------------

p1_visit1 <- FAO_data[1:6,]
p1_visit2 <- FAO_data[7:12,]
p2_visit1 <- FAO_data[13:18,]
p2_visit2 <- FAO_data[19:24,]
p3_visit1 <- FAO_data[25:30,]
p3_visit2 <- FAO_data[31:36,]
p4_visit1 <- FAO_data[37:41,] #misses 120 min
p4_visit2 <- FAO_data[42:47,]
p5_visit1 <- FAO_data[48:53,]
p5_visit2 <- FAO_data[54:59,]
p6_visit1 <- FAO_data[60:65,]
p6_visit2 <- FAO_data[66:71,]
p7_visit1 <- FAO_data[72:77,]
p7_visit2 <- FAO_data[78:83,]
p8_visit1 <- FAO_data[84:89,]
p8_visit2 <- FAO_data[90:95,]
p9_visit1 <- FAO_data[96:101,]
p9_visit2 <- FAO_data[102:107,]
p10_visit1 <- FAO_data[108:113,]
p10_visit2 <- FAO_data[114:119,]

# --------------------------------------
# Part 2: MCAR Simulation
# --------------------------------------

#function to impute row T = 120 for patient 4
add_row_120 <- function(data){
  data_copy <- data

  #new row using metadata from existing row (e.g., row with Time_min == 60)
  template_row <- data_copy[which.min(abs(data_copy$Time_min - 120)), ]  # Closest row (e.g., 60)
  
  new_row <- template_row
  new_row$Time_min <- 120
  new_row[6:ncol(new_row)] <- NA  #set all metabolite values to NA for this row
  
  #bind and re-sort
  data_copy <- bind_rows(data_copy, new_row) %>%
    arrange(Time_min)
  
  return(data_copy)
}

#call function for patient 4 visit 1 (NA introduced)
p4_visit1_full <- add_row_120(p4_visit1)


# #function to introduce 1 MCAR per dataframe randomly
# #middle values only (one missing value)
# MCAR_manipulation_middle <- function(data){
#   #copy dataset to avoid modifying the original
#   data_copy <- data
#   
#   #filter eligible rows (30, 60 or 120)
#   middle_rows <- which(data_copy$Time_min %in% c(30, 60, 120))
#   
#   for (col in colnames(data_copy[6:ncol(data_copy)])) {
#     rand_row <- sample(middle_rows, 1)
#     data_copy[rand_row, col] <- NA
#   }
#   
#   return(data_copy)
# }

#5th (120 min) values only (one missing value) and generate the column for patient which has no 120 min
MCAR_manipulation_middle <- function(data){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #find rows with T = 120
  rows_120 <- which(data_copy$Time_min == 120)
  
  #set them to NA if they exist
  if (length(rows_120) > 0){
    data_copy[rows_120, 6:ncol(data_copy)] <- NA
  }
  
  return(data_copy)
}


#call function
#p1
p1_v1_mcar <- MCAR_manipulation_middle(p1_visit1)
p1_v2_mcar <- MCAR_manipulation_middle(p1_visit2)
#p2
p2_v1_mcar <- MCAR_manipulation_middle(p2_visit1)
p2_v2_mcar <- MCAR_manipulation_middle(p2_visit2)
#p3
p3_v1_mcar <- MCAR_manipulation_middle(p3_visit1)
p3_v2_mcar <- MCAR_manipulation_middle(p3_visit2)
#p4
p4_v1_mcar <- MCAR_manipulation_middle(p4_visit1_full)
p4_v2_mcar <- MCAR_manipulation_middle(p4_visit2)
#p5
p5_v1_mcar <- MCAR_manipulation_middle(p5_visit1)
p5_v2_mcar <- MCAR_manipulation_middle(p5_visit2)
#p6
p6_v1_mcar <- MCAR_manipulation_middle(p6_visit1)
p6_v2_mcar <- MCAR_manipulation_middle(p6_visit2)
#p7
p7_v1_mcar <- MCAR_manipulation_middle(p7_visit1)
p7_v2_mcar <- MCAR_manipulation_middle(p7_visit2)
#p8
p8_v1_mcar <- MCAR_manipulation_middle(p8_visit1)
p8_v2_mcar <- MCAR_manipulation_middle(p8_visit2)
#p9
p9_v1_mcar <- MCAR_manipulation_middle(p9_visit1)
p9_v2_mcar <- MCAR_manipulation_middle(p9_visit2)
#p10
p10_v1_mcar <- MCAR_manipulation_middle(p10_visit1)
p10_v2_mcar <- MCAR_manipulation_middle(p10_visit2)

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
      return(na.approx(col, na.rm = FALSE, rule = 2)) #keep NA if interpolation not possible
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}

#call interpolation function on mcar data
#p1
p1_v1_mcar_interpolation <- interpolate_missing(p1_v1_mcar)
p1_v2_mcar_interpolation <- interpolate_missing(p1_v2_mcar)
#p2
p2_v1_mcar_interpolation <- interpolate_missing(p2_v1_mcar)
p2_v2_mcar_interpolation <- interpolate_missing(p2_v2_mcar)
#p3
p3_v1_mcar_interpolation <- interpolate_missing(p3_v1_mcar)
p3_v2_mcar_interpolation <- interpolate_missing(p3_v2_mcar)
#p4
p4_v1_mcar_interpolation <- interpolate_missing(p4_v1_mcar)
p4_v2_mcar_interpolation <- interpolate_missing(p4_v2_mcar)
#p5
p5_v1_mcar_interpolation <- interpolate_missing(p5_v1_mcar)
p5_v2_mcar_interpolation <- interpolate_missing(p5_v2_mcar)
#p6
p6_v1_mcar_interpolation <- interpolate_missing(p6_v1_mcar)
p6_v2_mcar_interpolation <- interpolate_missing(p6_v2_mcar)
#p7
p7_v1_mcar_interpolation <- interpolate_missing(p7_v1_mcar)
p7_v2_mcar_interpolation <- interpolate_missing(p7_v2_mcar)
#p8
p8_v1_mcar_interpolation <- interpolate_missing(p8_v1_mcar)
p8_v2_mcar_interpolation <- interpolate_missing(p8_v2_mcar)
#p9
p9_v1_mcar_interpolation <- interpolate_missing(p9_v1_mcar)
p9_v2_mcar_interpolation <- interpolate_missing(p9_v2_mcar)
#p10
p10_v1_mcar_interpolation <- interpolate_missing(p10_v1_mcar)
p10_v2_mcar_interpolation <- interpolate_missing(p10_v2_mcar)

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

#call function for kalman smoothing

#MCAR
#p1
p1_v1_mcar_kalman <- kalman_imputation(p1_v1_mcar)
p1_v2_mcar_kalman <- kalman_imputation(p1_v2_mcar)
#2
p2_v1_mcar_kalman <- kalman_imputation(p2_v1_mcar)
p2_v2_mcar_kalman <- kalman_imputation(p2_v2_mcar)
#3
p3_v1_mcar_kalman <- kalman_imputation(p3_v1_mcar)
p3_v2_mcar_kalman <- kalman_imputation(p3_v2_mcar)
#4
p4_v1_mcar_kalman <- kalman_imputation(p4_v1_mcar)
p4_v2_mcar_kalman <- kalman_imputation(p4_v2_mcar)
#5
p5_v1_mcar_kalman <- kalman_imputation(p5_v1_mcar)
p5_v2_mcar_kalman <- kalman_imputation(p5_v2_mcar)
#6
p6_v1_mcar_kalman <- kalman_imputation(p6_v1_mcar)
p6_v2_mcar_kalman <- kalman_imputation(p6_v2_mcar)
#7
p7_v1_mcar_kalman <- kalman_imputation(p7_v1_mcar)
p7_v2_mcar_kalman <- kalman_imputation(p7_v2_mcar)
#8
p8_v1_mcar_kalman <- kalman_imputation(p8_v1_mcar)
p8_v2_mcar_kalman <- kalman_imputation(p8_v2_mcar)
#9
p9_v1_mcar_kalman <- kalman_imputation(p9_v1_mcar)
p9_v2_mcar_kalman <- kalman_imputation(p9_v2_mcar)
#10
p10_v1_mcar_kalman <- kalman_imputation(p10_v1_mcar)
p10_v2_mcar_kalman <- kalman_imputation(p10_v2_mcar)

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

#call function for LWMA

#MCAR
#p1
p1_v1_mcar_lwma <- weighted_mov_average(p1_v1_mcar)
p1_v2_mcar_lwma <- weighted_mov_average(p1_v2_mcar)
#2
p2_v1_mcar_lwma <- weighted_mov_average(p2_v1_mcar)
p2_v2_mcar_lwma <- weighted_mov_average(p2_v2_mcar)
#3
p3_v1_mcar_lwma <- weighted_mov_average(p3_v1_mcar)
p3_v2_mcar_lwma <- weighted_mov_average(p3_v2_mcar)
#4
p4_v1_mcar_lwma <- weighted_mov_average(p4_v1_mcar)
p4_v2_mcar_lwma <- weighted_mov_average(p4_v2_mcar)
#5
p5_v1_mcar_lwma <- weighted_mov_average(p5_v1_mcar)
p5_v2_mcar_lwma <- weighted_mov_average(p5_v2_mcar)
#6
p6_v1_mcar_lwma <- weighted_mov_average(p6_v1_mcar)
p6_v2_mcar_lwma <- weighted_mov_average(p6_v2_mcar)
#7
p7_v1_mcar_lwma <- weighted_mov_average(p7_v1_mcar)
p7_v2_mcar_lwma <- weighted_mov_average(p7_v2_mcar)
#8
p8_v1_mcar_lwma <- weighted_mov_average(p8_v1_mcar)
p8_v2_mcar_lwma <- weighted_mov_average(p8_v2_mcar)
#9
p9_v1_mcar_lwma <- weighted_mov_average(p9_v1_mcar)
p9_v2_mcar_lwma <- weighted_mov_average(p9_v2_mcar)
#10
p10_v1_mcar_lwma <- weighted_mov_average(p10_v1_mcar)
p10_v2_mcar_lwma <- weighted_mov_average(p10_v2_mcar)

# --------------------------------------
# TITLE: KINETICS PLOT BEFORE AND AFTER
# --------------------------------------
# -------------------------------
# Plot before and After 
# -------------------------------

#plot dots and line
plot_imputed_vs_original <- function(original, imputed, visit, type){
  #add source info 
  original_df <- original %>% mutate(Source = "Original")
  imputed_df <- imputed %>% mutate(Source = "Imputed")
  
  #combine data
  combined <- bind_rows(original_df, imputed_df)
  
  #ensure correct type
  combined <- combined %>%
    mutate(
      Visit = as.factor(Visit),
      Time_min = as.numeric(Time_min),
      Source = as.factor(Source)
    )
  
  #reshape into long
  long_df <- combined %>%
    pivot_longer(cols = 6:(ncol(combined) - 1), #last column source
                 names_to = "Metabolite",
                 values_to = "Concentration")
  
  #pateint list 
  patients <- unique(long_df$Patient)
  
  for (patient in patients){
    plot_data <- long_df %>%
      filter(Patient == patient)
    
    p <- ggplot(plot_data, aes(x = Time_min, y = Concentration, color = Source)) +
      geom_point(size = 2) +
      geom_line(aes(group = Source), linewidth =1) +
      facet_wrap(~ Metabolite, scales = "free_y", ncol = 6) +
      labs(
        title = paste(patient, ":", "Kintecs of all Metabolites (", visit, ", ", type, ")"),
        subtitle = "Original vs Imputed",
        x = "Time [min]",
        y = "Concentration [µM]"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14)
      ) +
      scale_color_manual(values = c("Original" = "darkblue", "Imputed" = "red"))
    
    print(p)
  }
}

# ----------------------------
# Part 1: Linear Interpolation
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Interpolation/MCAR_Interpolation_5thMV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p1_visit2, p1_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p2_visit2, p2_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p3_visit2, p3_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p4_visit2, p4_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p5_visit2, p5_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p6_visit2, p6_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p7_visit2, p7_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p8_visit2, p8_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p9_visit2, p9_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mcar_interpolation, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p10_visit2, p10_v2_mcar_interpolation, visit = "Visit 2", type = "MCAR")

dev.off()

# ----------------------------
# Part 2: Kalman Smoothing
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Kalman_Smoothing/MCAR_Kalman_5thMV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p1_visit2, p1_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p2_visit2, p2_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p3_visit2, p3_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p4_visit2, p4_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p5_visit2, p5_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p6_visit2, p6_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p7_visit2, p7_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p8_visit2, p8_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p9_visit2, p9_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mcar_kalman, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p10_visit2, p10_v2_mcar_kalman, visit = "Visit 2", type = "MCAR")

dev.off()


# ----------------------------
# Part 3: LWMA
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/LWMA/MCAR_LWMA_5thMV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p1_visit2, p1_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p2_visit2, p2_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p3_visit2, p3_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p4_visit2, p4_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p5_visit2, p5_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p6_visit2, p6_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p7_visit2, p7_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p8_visit2, p8_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p9_visit2, p9_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mcar_lwma, visit = "Visit 1", type = "MCAR")
plot_imputed_vs_original(p10_visit2, p10_v2_mcar_lwma, visit = "Visit 2", type = "MCAR")

dev.off()


# ------------------------
# TITLE: NRMSE Function
# ------------------------

#function to calculate nrmse 
calculate_nrsme <- function(original, imputed, method) {
  numeric_col_names <- names(original)[6:ncol(original)]
  
  nrmse_values <- sapply(numeric_col_names, function(col) {
    actual_val <- original[[col]]
    imputed_val <- imputed[[col]]
    
    valid_indices <- !is.na(actual_val) & !is.na(imputed_val)
    
    if (sum(valid_indices) > 2) {
      mse <- mean((actual_val[valid_indices] - imputed_val[valid_indices])^2)
      rmse <- sqrt(mse)
      norm_factor <- max(actual_val[valid_indices], na.rm = TRUE) - min(actual_val[valid_indices], na.rm = TRUE)
      
      if (norm_factor > 0) {
        return(rmse / norm_factor)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
  
  return(data.frame(
    Metabolite = numeric_col_names,
    Imputation_method = method,
    NRMSE = nrmse_values
  ))
}

# ----------------------------
# Part 1: Linear Interpolation
# ----------------------------

# --------------------------
# Part 1.1: NRMSE (MCAR)
# --------------------------

#call function to calcualte nrms
#interpolation
#p1
nrmse_interp_p1v1_mcar <- calculate_nrsme(p1_visit1, p1_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p1v2_mcar <- calculate_nrsme(p1_visit2, p1_v2_mcar_interpolation, method = "Linear Interpolation")
#p2
nrmse_interp_p2v1_mcar <- calculate_nrsme(p2_visit1, p2_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p2v2_mcar <- calculate_nrsme(p2_visit2, p2_v2_mcar_interpolation, method = "Linear Interpolation")
#p3
nrmse_interp_p3v1_mcar <- calculate_nrsme(p3_visit1, p3_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p3v2_mcar <- calculate_nrsme(p3_visit2, p3_v2_mcar_interpolation, method = "Linear Interpolation")
#p4
nrmse_interp_p4v1_mcar <- calculate_nrsme(p4_visit1, p4_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p4v2_mcar <- calculate_nrsme(p4_visit2, p4_v2_mcar_interpolation, method = "Linear Interpolation")
#p5
nrmse_interp_p5v1_mcar <- calculate_nrsme(p5_visit1, p5_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p5v2_mcar <- calculate_nrsme(p5_visit2, p5_v2_mcar_interpolation, method = "Linear Interpolation")
#p6
nrmse_interp_p6v1_mcar <- calculate_nrsme(p6_visit1, p6_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p6v2_mcar <- calculate_nrsme(p6_visit2, p6_v2_mcar_interpolation, method = "Linear Interpolation")
#p7
nrmse_interp_p7v1_mcar <- calculate_nrsme(p7_visit1, p7_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p7v2_mcar <- calculate_nrsme(p7_visit2, p7_v2_mcar_interpolation, method = "Linear Interpolation")
#p8
nrmse_interp_p8v1_mcar <- calculate_nrsme(p8_visit1, p8_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p8v2_mcar <- calculate_nrsme(p8_visit2, p8_v2_mcar_interpolation, method = "Linear Interpolation")
#p9
nrmse_interp_p9v1_mcar <- calculate_nrsme(p9_visit1, p9_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p9v2_mcar <- calculate_nrsme(p9_visit2, p9_v2_mcar_interpolation, method = "Linear Interpolation")
#p10
nrmse_interp_p10v1_mcar <- calculate_nrsme(p10_visit1, p10_v1_mcar_interpolation, method = "Linear Interpolation")
nrmse_interp_p10v2_mcar <- calculate_nrsme(p10_visit2, p10_v2_mcar_interpolation, method = "Linear Interpolation")

#combine visit 1 
nrmse_mcar_visit1 <- bind_rows(
  nrmse_interp_p1v1_mcar %>% mutate(Patient = "P1"),
  nrmse_interp_p2v1_mcar %>% mutate(Patient = "P2"),
  nrmse_interp_p3v1_mcar %>% mutate(Patient = "P3"),
  nrmse_interp_p4v1_mcar %>% mutate(Patient = "P4"),
  nrmse_interp_p5v1_mcar %>% mutate(Patient = "P5"),
  nrmse_interp_p6v1_mcar %>% mutate(Patient = "P6"),
  nrmse_interp_p7v1_mcar %>% mutate(Patient = "P7"),
  nrmse_interp_p8v1_mcar %>% mutate(Patient = "P8"),
  nrmse_interp_p9v1_mcar %>% mutate(Patient = "P9"),
  nrmse_interp_p10v1_mcar %>% mutate(Patient = "P10")
)

#combine visit 2
nrmse_mcar_visit2 <- bind_rows(
  nrmse_interp_p1v2_mcar %>% mutate(Patient = "P1"),
  nrmse_interp_p2v2_mcar %>% mutate(Patient = "P2"),
  nrmse_interp_p3v2_mcar %>% mutate(Patient = "P3"),
  nrmse_interp_p4v2_mcar %>% mutate(Patient = "P4"),
  nrmse_interp_p5v2_mcar %>% mutate(Patient = "P5"),
  nrmse_interp_p6v2_mcar %>% mutate(Patient = "P6"),
  nrmse_interp_p7v2_mcar %>% mutate(Patient = "P7"),
  nrmse_interp_p8v2_mcar %>% mutate(Patient = "P8"),
  nrmse_interp_p9v2_mcar %>% mutate(Patient = "P9"),
  nrmse_interp_p10v2_mcar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Interpolation/MCAR_Interpolation_5thMV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_mcar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1 (MCAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_mcar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2 (MCAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

dev.off()


# ----------------------------
# Part 2: Kalman
# ----------------------------

# --------------------------
# Part 2.1: NRMSE (MCAR)
# --------------------------

#call function to calcualte nrms
#kalman
#p1
nrmse_kalman_p1v1_mcar <- calculate_nrsme(p1_visit1, p1_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p1v2_mcar <- calculate_nrsme(p1_visit2, p1_v2_mcar_kalman, method = "Kalman Smoothing")
#p2
nrmse_kalman_p2v1_mcar <- calculate_nrsme(p2_visit1, p2_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p2v2_mcar <- calculate_nrsme(p2_visit2, p2_v2_mcar_kalman, method = "Kalman Smoothing")
#p3
nrmse_kalman_p3v1_mcar <- calculate_nrsme(p3_visit1, p3_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p3v2_mcar <- calculate_nrsme(p3_visit2, p3_v2_mcar_kalman, method = "Kalman Smoothing")
#p4
nrmse_kalman_p4v1_mcar <- calculate_nrsme(p4_visit1, p4_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p4v2_mcar <- calculate_nrsme(p4_visit2, p4_v2_mcar_kalman, method = "Kalman Smoothing")
#p5
nrmse_kalman_p5v1_mcar <- calculate_nrsme(p5_visit1, p5_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p5v2_mcar <- calculate_nrsme(p5_visit2, p5_v2_mcar_kalman, method = "Kalman Smoothing")
#p6
nrmse_kalman_p6v1_mcar <- calculate_nrsme(p6_visit1, p6_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p6v2_mcar <- calculate_nrsme(p6_visit2, p6_v2_mcar_kalman, method = "Kalman Smoothing")
#p7
nrmse_kalman_p7v1_mcar <- calculate_nrsme(p7_visit1, p7_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p7v2_mcar <- calculate_nrsme(p7_visit2, p7_v2_mcar_kalman, method = "Kalman Smoothing")
#p8
nrmse_kalman_p8v1_mcar <- calculate_nrsme(p8_visit1, p8_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p8v2_mcar <- calculate_nrsme(p8_visit2, p8_v2_mcar_kalman, method = "Kalman Smoothing")
#p9
nrmse_kalman_p9v1_mcar <- calculate_nrsme(p9_visit1, p9_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p9v2_mcar <- calculate_nrsme(p9_visit2, p9_v2_mcar_kalman, method = "Kalman Smoothing")
#p10
nrmse_kalman_p10v1_mcar <- calculate_nrsme(p10_visit1, p10_v1_mcar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p10v2_mcar <- calculate_nrsme(p10_visit2, p10_v2_mcar_kalman, method = "Kalman Smoothing")

#combine visit 1 
nrmse_kalman_mcar_visit1 <- bind_rows(
  nrmse_kalman_p1v1_mcar %>% mutate(Patient = "P1"),
  nrmse_kalman_p2v1_mcar %>% mutate(Patient = "P2"),
  nrmse_kalman_p3v1_mcar %>% mutate(Patient = "P3"),
  nrmse_kalman_p4v1_mcar %>% mutate(Patient = "P4"),
  nrmse_kalman_p5v1_mcar %>% mutate(Patient = "P5"),
  nrmse_kalman_p6v1_mcar %>% mutate(Patient = "P6"),
  nrmse_kalman_p7v1_mcar %>% mutate(Patient = "P7"),
  nrmse_kalman_p8v1_mcar %>% mutate(Patient = "P8"),
  nrmse_kalman_p9v1_mcar %>% mutate(Patient = "P9"),
  nrmse_kalman_p10v1_mcar %>% mutate(Patient = "P10")
)

#combine visit 2
nrmse_kalman_mcar_visit2 <- bind_rows(
  nrmse_kalman_p1v2_mcar %>% mutate(Patient = "P1"),
  nrmse_kalman_p2v2_mcar %>% mutate(Patient = "P2"),
  nrmse_kalman_p3v2_mcar %>% mutate(Patient = "P3"),
  nrmse_kalman_p4v2_mcar %>% mutate(Patient = "P4"),
  nrmse_kalman_p5v2_mcar %>% mutate(Patient = "P5"),
  nrmse_kalman_p6v2_mcar %>% mutate(Patient = "P6"),
  nrmse_kalman_p7v2_mcar %>% mutate(Patient = "P7"),
  nrmse_kalman_p8v2_mcar %>% mutate(Patient = "P8"),
  nrmse_kalman_p9v2_mcar %>% mutate(Patient = "P9"),
  nrmse_kalman_p10v2_mcar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Kalman_Smoothing/MCAR_Kalman_5thMV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_kalman_mcar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1 (MCAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_kalman_mcar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2 (MCAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

dev.off()


# ----------------------------
# Part 3: LWMA
# ----------------------------

# --------------------------
# Part 3.1: NRMSE (MCAR)
# --------------------------

#call function to calcualte nrms
#LWMA
#p1
nrmse_lwma_p1v1_mcar <- calculate_nrsme(p1_visit1, p1_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p1v2_mcar <- calculate_nrsme(p1_visit2, p1_v2_mcar_lwma, method = "LWMA")
#p2
nrmse_lwma_p2v1_mcar <- calculate_nrsme(p2_visit1, p2_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p2v2_mcar <- calculate_nrsme(p2_visit2, p2_v2_mcar_lwma, method = "LWMA")
#p3
nrmse_lwma_p3v1_mcar <- calculate_nrsme(p3_visit1, p3_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p3v2_mcar <- calculate_nrsme(p3_visit2, p3_v2_mcar_lwma, method = "LWMA")
#p4
nrmse_lwma_p4v1_mcar <- calculate_nrsme(p4_visit1, p4_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p4v2_mcar <- calculate_nrsme(p4_visit2, p4_v2_mcar_lwma, method = "LWMA")
#p5
nrmse_lwma_p5v1_mcar <- calculate_nrsme(p5_visit1, p5_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p5v2_mcar <- calculate_nrsme(p5_visit2, p5_v2_mcar_lwma, method = "LWMA")
#p6
nrmse_lwma_p6v1_mcar <- calculate_nrsme(p6_visit1, p6_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p6v2_mcar <- calculate_nrsme(p6_visit2, p6_v2_mcar_lwma, method = "LWMA")
#p7
nrmse_lwma_p7v1_mcar <- calculate_nrsme(p7_visit1, p7_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p7v2_mcar <- calculate_nrsme(p7_visit2, p7_v2_mcar_lwma, method = "LWMA")
#p8
nrmse_lwma_p8v1_mcar <- calculate_nrsme(p8_visit1, p8_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p8v2_mcar <- calculate_nrsme(p8_visit2, p8_v2_mcar_lwma, method = "LWMA")
#p9
nrmse_lwma_p9v1_mcar <- calculate_nrsme(p9_visit1, p9_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p9v2_mcar <- calculate_nrsme(p9_visit2, p9_v2_mcar_lwma, method = "LWMA")
#p10
nrmse_lwma_p10v1_mcar <- calculate_nrsme(p10_visit1, p10_v1_mcar_lwma, method = "LWMA")
nrmse_lwma_p10v2_mcar <- calculate_nrsme(p10_visit2, p10_v2_mcar_lwma, method = "LWMA")

#combine visit 1 
nrmse_lwma_mcar_visit1 <- bind_rows(
  nrmse_lwma_p1v1_mcar %>% mutate(Patient = "P1"),
  nrmse_lwma_p2v1_mcar %>% mutate(Patient = "P2"),
  nrmse_lwma_p3v1_mcar %>% mutate(Patient = "P3"),
  nrmse_lwma_p4v1_mcar %>% mutate(Patient = "P4"),
  nrmse_lwma_p5v1_mcar %>% mutate(Patient = "P5"),
  nrmse_lwma_p6v1_mcar %>% mutate(Patient = "P6"),
  nrmse_lwma_p7v1_mcar %>% mutate(Patient = "P7"),
  nrmse_lwma_p8v1_mcar %>% mutate(Patient = "P8"),
  nrmse_lwma_p9v1_mcar %>% mutate(Patient = "P9"),
  nrmse_lwma_p10v1_mcar %>% mutate(Patient = "P10")
)

#combine visit 2
nrmse_lwma_mcar_visit2 <- bind_rows(
  nrmse_lwma_p1v2_mcar %>% mutate(Patient = "P1"),
  nrmse_lwma_p2v2_mcar %>% mutate(Patient = "P2"),
  nrmse_lwma_p3v2_mcar %>% mutate(Patient = "P3"),
  nrmse_lwma_p4v2_mcar %>% mutate(Patient = "P4"),
  nrmse_lwma_p5v2_mcar %>% mutate(Patient = "P5"),
  nrmse_lwma_p6v2_mcar %>% mutate(Patient = "P6"),
  nrmse_lwma_p7v2_mcar %>% mutate(Patient = "P7"),
  nrmse_lwma_p8v2_mcar %>% mutate(Patient = "P8"),
  nrmse_lwma_p9v2_mcar %>% mutate(Patient = "P9"),
  nrmse_lwma_p10v2_mcar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/LWMA/MCAR_LWMA_5thMV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_lwma_mcar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1 (MCAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_lwma_mcar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2 (MCAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

dev.off()


# ----------------------------
# Part 4: All methods compared
# ----------------------------

#visit1
nrmse_visit1_tot <- bind_rows(
  nrmse_mcar_visit1,
  nrmse_kalman_mcar_visit1,
  nrmse_lwma_mcar_visit1
)

#visit2
nrmse_visit2_tot <- bind_rows(
  nrmse_mcar_visit2,
  nrmse_kalman_mcar_visit2,
  nrmse_lwma_mcar_visit2
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/NRMSE_MCAR_Imputation_methods.pdf", width = 16, height = 10)

#plot visit 1
ggplot(nrmse_visit1_tot, aes(x = Imputation_method, y = NRMSE, fill = Imputation_method)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, outlier.fill = "white") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "NRMSE of Imputed Values Only (Visit 1 - MCAR)",
    x = "Imputation Method",
    y = "Normalized RMSE"
  )

#plot visit 1
ggplot(nrmse_visit2_tot, aes(x = Imputation_method, y = NRMSE, fill = Imputation_method)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, outlier.fill = "white") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) +
  labs(
    title = "NRMSE of Imputed Values Only (Visit 2 - MCAR)",
    x = "Imputation Method",
    y = "Normalized RMSE"
  )

dev.off()

# ------------------------
# TITLE: AUC Function
# ------------------------

#function to calculate AUC
calculate_auc <- function(data){
  time <- data$Time_min
  aucs <- sapply(data[, 6:ncol(data)], function(metabolite){
    if (all(is.na(metabolite))){ #if all is NA skip completely
      return(NA)
    } else {
      return(trapz(time, metabolite)) #trapezodial AUC
    }
  })
  return(aucs)
}

# ----------------------------------
# Part 1: For 5th Value (T = 120 min) 
# -----------------------------------

#original AUC
#visit 1
auc_p1v1_original <- calculate_auc(p1_visit1)
auc_p2v1_original <- calculate_auc(p2_visit1)
auc_p3v1_original <- calculate_auc(p3_visit1)
auc_p4v1_original <- calculate_auc(p4_visit1)
auc_p5v1_original <- calculate_auc(p5_visit1)
auc_p6v1_original <- calculate_auc(p6_visit1)
auc_p7v1_original <- calculate_auc(p7_visit1)
auc_p8v1_original <- calculate_auc(p8_visit1)
auc_p9v1_original <- calculate_auc(p9_visit1)
auc_p10v1_original <- calculate_auc(p10_visit1)
#visit 2
auc_p1v2_original <- calculate_auc(p1_visit2)
auc_p2v2_original <- calculate_auc(p2_visit2)
auc_p3v2_original <- calculate_auc(p3_visit2)
auc_p4v2_original <- calculate_auc(p4_visit2)
auc_p5v2_original <- calculate_auc(p5_visit2)
auc_p6v2_original <- calculate_auc(p6_visit2)
auc_p7v2_original <- calculate_auc(p7_visit2)
auc_p8v2_original <- calculate_auc(p8_visit2)
auc_p9v2_original <- calculate_auc(p9_visit2)
auc_p10v2_original <- calculate_auc(p10_visit2)

#combine all 
#visit 1
original_visit1_auc <- bind_rows(
  auc_p1v1_original, 
  auc_p2v1_original,
  auc_p3v1_original,
  auc_p4v1_original,
  auc_p5v1_original,
  auc_p6v1_original,
  auc_p7v1_original,
  auc_p8v1_original,
  auc_p9v1_original,
  auc_p10v1_original
)
#visit 2
original_visit2_auc <- bind_rows(
  auc_p1v2_original, 
  auc_p2v2_original,
  auc_p3v2_original,
  auc_p4v2_original,
  auc_p5v2_original,
  auc_p6v2_original,
  auc_p7v2_original,
  auc_p8v2_original,
  auc_p9v2_original,
  auc_p10v2_original
)

#interpolation AUC
#visit 1
auc_p1v1_interpolation <- calculate_auc(p1_v1_mcar_interpolation)
auc_p2v1_interpolation <- calculate_auc(p2_v1_mcar_interpolation)
auc_p3v1_interpolation <- calculate_auc(p3_v1_mcar_interpolation)
auc_p4v1_interpolation <- calculate_auc(p4_v1_mcar_interpolation)
auc_p5v1_interpolation <- calculate_auc(p5_v1_mcar_interpolation)
auc_p6v1_interpolation <- calculate_auc(p6_v1_mcar_interpolation)
auc_p7v1_interpolation <- calculate_auc(p7_v1_mcar_interpolation)
auc_p8v1_interpolation <- calculate_auc(p8_v1_mcar_interpolation)
auc_p9v1_interpolation <- calculate_auc(p9_v1_mcar_interpolation)
auc_p10v1_interpolation <- calculate_auc(p10_v1_mcar_interpolation)
#visit 2
auc_p1v2_interpolation <- calculate_auc(p1_v2_mcar_interpolation)
auc_p2v2_interpolation <- calculate_auc(p2_v2_mcar_interpolation)
auc_p3v2_interpolation <- calculate_auc(p3_v2_mcar_interpolation)
auc_p4v2_interpolation <- calculate_auc(p4_v2_mcar_interpolation)
auc_p5v2_interpolation <- calculate_auc(p5_v2_mcar_interpolation)
auc_p6v2_interpolation <- calculate_auc(p6_v2_mcar_interpolation)
auc_p7v2_interpolation <- calculate_auc(p7_v2_mcar_interpolation)
auc_p8v2_interpolation <- calculate_auc(p8_v2_mcar_interpolation)
auc_p9v2_interpolation <- calculate_auc(p9_v2_mcar_interpolation)
auc_p10v2_interpolation <- calculate_auc(p10_v2_mcar_interpolation)

#combine
#visit 1
interpolation_visit1_auc <- bind_rows(
  auc_p1v1_interpolation, 
  auc_p2v1_interpolation,
  auc_p3v1_interpolation,
  auc_p4v1_interpolation,
  auc_p5v1_interpolation,
  auc_p6v1_interpolation,
  auc_p7v1_interpolation,
  auc_p8v1_interpolation,
  auc_p9v1_interpolation,
  auc_p10v1_interpolation
)
#visit 2
interpolation_visit2_auc <- bind_rows(
  auc_p1v2_interpolation, 
  auc_p2v2_interpolation,
  auc_p3v2_interpolation,
  auc_p4v2_interpolation,
  auc_p5v2_interpolation,
  auc_p6v2_interpolation,
  auc_p7v2_interpolation,
  auc_p8v2_interpolation,
  auc_p9v2_interpolation,
  auc_p10v2_interpolation
)

#kalman AUC
#visit 1
auc_p1v1_kalman <- calculate_auc(p1_v1_mcar_kalman)
auc_p2v1_kalman <- calculate_auc(p2_v1_mcar_kalman)
auc_p3v1_kalman <- calculate_auc(p3_v1_mcar_kalman)
auc_p4v1_kalman <- calculate_auc(p4_v1_mcar_kalman)
auc_p5v1_kalman <- calculate_auc(p5_v1_mcar_kalman)
auc_p6v1_kalman <- calculate_auc(p6_v1_mcar_kalman)
auc_p7v1_kalman <- calculate_auc(p7_v1_mcar_kalman)
auc_p8v1_kalman <- calculate_auc(p8_v1_mcar_kalman)
auc_p9v1_kalman <- calculate_auc(p9_v1_mcar_kalman)
auc_p10v1_kalman <- calculate_auc(p10_v1_mcar_kalman)
#visit 2
auc_p1v2_kalman <- calculate_auc(p1_v2_mcar_kalman)
auc_p2v2_kalman <- calculate_auc(p2_v2_mcar_kalman)
auc_p3v2_kalman <- calculate_auc(p3_v2_mcar_kalman)
auc_p4v2_kalman <- calculate_auc(p4_v2_mcar_kalman)
auc_p5v2_kalman <- calculate_auc(p5_v2_mcar_kalman)
auc_p6v2_kalman <- calculate_auc(p6_v2_mcar_kalman)
auc_p7v2_kalman <- calculate_auc(p7_v2_mcar_kalman)
auc_p8v2_kalman <- calculate_auc(p8_v2_mcar_kalman)
auc_p9v2_kalman <- calculate_auc(p9_v2_mcar_kalman)
auc_p10v2_kalman <- calculate_auc(p10_v2_mcar_kalman)

#combine
kalman_visit1_auc <- bind_rows(
  auc_p1v1_kalman, 
  auc_p2v1_kalman,
  auc_p3v1_kalman,
  auc_p4v1_kalman,
  auc_p5v1_kalman,
  auc_p6v1_kalman,
  auc_p7v1_kalman,
  auc_p8v1_kalman,
  auc_p9v1_kalman,
  auc_p10v1_kalman
)
#visit 2
kalman_visit2_auc <- bind_rows(
  auc_p1v2_kalman, 
  auc_p2v2_kalman,
  auc_p3v2_kalman,
  auc_p4v2_kalman,
  auc_p5v2_kalman,
  auc_p6v2_kalman,
  auc_p7v2_kalman,
  auc_p8v2_kalman,
  auc_p9v2_kalman,
  auc_p10v2_kalman
)

#LWMA AUC
#visit 1
auc_p1v1_lwma <- calculate_auc(p1_v1_mcar_lwma)
auc_p2v1_lwma <- calculate_auc(p2_v1_mcar_lwma)
auc_p3v1_lwma <- calculate_auc(p3_v1_mcar_lwma)
auc_p4v1_lwma <- calculate_auc(p4_v1_mcar_lwma)
auc_p5v1_lwma <- calculate_auc(p5_v1_mcar_lwma)
auc_p6v1_lwma <- calculate_auc(p6_v1_mcar_lwma)
auc_p7v1_lwma <- calculate_auc(p7_v1_mcar_lwma)
auc_p8v1_lwma <- calculate_auc(p8_v1_mcar_lwma)
auc_p9v1_lwma <- calculate_auc(p9_v1_mcar_lwma)
auc_p10v1_lwma <- calculate_auc(p10_v1_mcar_lwma)
#visit 2
auc_p1v2_lwma <- calculate_auc(p1_v2_mcar_lwma)
auc_p2v2_lwma <- calculate_auc(p2_v2_mcar_lwma)
auc_p3v2_lwma <- calculate_auc(p3_v2_mcar_lwma)
auc_p4v2_lwma <- calculate_auc(p4_v2_mcar_lwma)
auc_p5v2_lwma <- calculate_auc(p5_v2_mcar_lwma)
auc_p6v2_lwma <- calculate_auc(p6_v2_mcar_lwma)
auc_p7v2_lwma <- calculate_auc(p7_v2_mcar_lwma)
auc_p8v2_lwma <- calculate_auc(p8_v2_mcar_lwma)
auc_p9v2_lwma <- calculate_auc(p9_v2_mcar_lwma)
auc_p10v2_lwma <- calculate_auc(p10_v2_mcar_lwma)

#combine
#visit 1
lwma_visit1_auc <- bind_rows(
  auc_p1v1_lwma, 
  auc_p2v1_lwma,
  auc_p3v1_lwma,
  auc_p4v1_lwma,
  auc_p5v1_lwma,
  auc_p6v1_lwma,
  auc_p7v1_lwma,
  auc_p8v1_lwma,
  auc_p9v1_lwma,
  auc_p10v1_lwma
)
#visit 2
lwma_visit2_auc <- bind_rows(
  auc_p1v2_lwma, 
  auc_p2v2_lwma,
  auc_p3v2_lwma,
  auc_p4v2_lwma,
  auc_p5v2_lwma,
  auc_p6v2_lwma,
  auc_p7v2_lwma,
  auc_p8v2_lwma,
  auc_p8v2_lwma,
  auc_p10v2_lwma
)

#combine the datasets
#visit 1
visit1_auc_df <- bind_rows(
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p1v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p2v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p3v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p4v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p5v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p6v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p7v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p8v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p9v1_original)),
  data.frame(Method = "Original",      Visit = "Visit 1", stack(auc_p10v1_original)),
  
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p1v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p2v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p3v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p4v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p5v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p6v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p7v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p8v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p9v1_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 1", stack(auc_p10v1_interpolation)),
  
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p1v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p2v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p3v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p4v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p5v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p6v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p7v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p8v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p9v1_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 1", stack(auc_p10v1_kalman)),
  
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p1v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p2v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p3v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p4v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p5v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p6v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p7v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p8v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p9v1_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p10v1_lwma))
) %>% rename(AUC = values, Metabolite = ind)

#visit 2
visit2_auc_df <- bind_rows(
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p1v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p2v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p3v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p4v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p5v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p6v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p7v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p8v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p9v2_original)),
  data.frame(Method = "Original",      Visit = "Visit 2", stack(auc_p10v2_original)),
  
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p1v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p2v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p3v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p4v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p5v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p6v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p7v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p8v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p9v2_interpolation)),
  data.frame(Method = "Interpolation", Visit = "Visit 2", stack(auc_p10v2_interpolation)),
  
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p1v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p2v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p3v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p4v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p5v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p6v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p7v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p8v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p9v2_kalman)),
  data.frame(Method = "Kalman",        Visit = "Visit 2", stack(auc_p10v2_kalman)),
  
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p1v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p2v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p3v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p4v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p5v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p6v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p7v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p8v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p9v2_lwma)),
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p10v2_lwma))
) %>% rename(AUC = values, Metabolite = ind)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/AUC_Density_MCAR.pdf", width = 16, height = 10)

#now plot the density of AUC 
#Visit 1
ggplot(visit1_auc_df, aes(x = AUC, fill = Method, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~ Metabolite, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Visit 1: AUC Density per Metabolite",
    x = "AUC",
    y = "Density"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

#Visit 2
ggplot(visit2_auc_df, aes(x = AUC, fill = Method, color = Method)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~ Metabolite, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Visit 1: AUC Density per Metabolite",
    x = "AUC",
    y = "Density"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

dev.off()

