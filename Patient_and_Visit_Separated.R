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

#function to introduce 1 MCAR per dataframe randomly

#middle values only (one missing value)
MCAR_manipulation_middle <- function(data){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #filter eligible rows (30, 60 or 120)
  middle_rows <- which(data_copy$Time_min %in% c(30, 60, 120))
  
  for (col in colnames(data_copy[6:ncol(data_copy)])) {
    rand_row <- sample(middle_rows, 1)
    data_copy[rand_row, col] <- NA
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
p4_v1_mcar <- MCAR_manipulation_middle(p4_visit1)
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
# Part 3: MNAR Simulation
# --------------------------------------

#function to introduce 1 MNAR per dataframe (one missing value)
MNAR_manipulation_lowest <- function(data){
  #copy dataset to avoid modifying the original
  data_copy <- data

  #go through each column
  for (col in colnames(data_copy[6:ncol(data_copy)])) {
    min_row <- which.min(data_copy[[col]])
    data_copy[min_row, col] <- NA
  }

  return(data_copy)
}

# #function to introduce 2 MNAR per dataframe (one missing value)
# MNAR_manipulation_lowest <- function(data){
#   #copy dataset to avoid modifying the original
#   data_copy <- data
#   
#   #go through each column
#   for (col in colnames(data_copy[6:ncol(data_copy)])) {
#     #indices fot rhe two lowest value
#     min_rows <- order(data_copy[[col]], na.last = NA)[1:2]
#     data_copy[min_rows, col] <- NA
#   }
#   
#   return(data_copy)
# }


#call function
#p1
p1_v1_mnar <- MNAR_manipulation_lowest(p1_visit1)
p1_v2_mnar <- MNAR_manipulation_lowest(p1_visit2)
#p2
p2_v1_mnar <- MNAR_manipulation_lowest(p2_visit1)
p2_v2_mnar <- MNAR_manipulation_lowest(p2_visit2)
#p3
p3_v1_mnar <- MNAR_manipulation_lowest(p3_visit1)
p3_v2_mnar <- MNAR_manipulation_lowest(p3_visit2)
#p4
p4_v1_mnar <- MNAR_manipulation_lowest(p4_visit1)
p4_v2_mnar <- MNAR_manipulation_lowest(p4_visit2)
#p5
p5_v1_mnar <- MNAR_manipulation_lowest(p5_visit1)
p5_v2_mnar <- MNAR_manipulation_lowest(p5_visit2)
#p6
p6_v1_mnar <- MNAR_manipulation_lowest(p6_visit1)
p6_v2_mnar <- MNAR_manipulation_lowest(p6_visit2)
#p7
p7_v1_mnar <- MNAR_manipulation_lowest(p7_visit1)
p7_v2_mnar <- MNAR_manipulation_lowest(p7_visit2)
#p8
p8_v1_mnar <- MNAR_manipulation_lowest(p8_visit1)
p8_v2_mnar <- MNAR_manipulation_lowest(p8_visit2)
#p9
p9_v1_mnar <- MNAR_manipulation_lowest(p9_visit1)
p9_v2_mnar <- MNAR_manipulation_lowest(p9_visit2)
#p10
p10_v1_mnar <- MNAR_manipulation_lowest(p10_visit1)
p10_v2_mnar <- MNAR_manipulation_lowest(p10_visit2)


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
      return(na.approx(col, na.rm = FALSE)) #keep NA if interpolation not possible
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

#call interpolation function on mnar data
#p1
p1_v1_mnar_interpolation <- interpolate_missing(p1_v1_mnar)
p1_v2_mnar_interpolation <- interpolate_missing(p1_v2_mnar)
#p2
p2_v1_mnar_interpolation <- interpolate_missing(p2_v1_mnar)
p2_v2_mnar_interpolation <- interpolate_missing(p2_v2_mnar)
#p3
p3_v1_mnar_interpolation <- interpolate_missing(p3_v1_mnar)
p3_v2_mnar_interpolation <- interpolate_missing(p3_v2_mnar)
#p4
p4_v1_mnar_interpolation <- interpolate_missing(p4_v1_mnar)
p4_v2_mnar_interpolation <- interpolate_missing(p4_v2_mnar)
#p5
p5_v1_mnar_interpolation <- interpolate_missing(p5_v1_mnar)
p5_v2_mnar_interpolation <- interpolate_missing(p5_v2_mnar)
#p6
p6_v1_mnar_interpolation <- interpolate_missing(p6_v1_mnar)
p6_v2_mnar_interpolation <- interpolate_missing(p6_v2_mnar)
#p7
p7_v1_mnar_interpolation <- interpolate_missing(p7_v1_mnar)
p7_v2_mnar_interpolation <- interpolate_missing(p7_v2_mnar)
#p8
p8_v1_mnar_interpolation <- interpolate_missing(p8_v1_mnar)
p8_v2_mnar_interpolation <- interpolate_missing(p8_v2_mnar)
#p9
p9_v1_mnar_interpolation <- interpolate_missing(p9_v1_mnar)
p9_v2_mnar_interpolation <- interpolate_missing(p9_v2_mnar)
#p10
p10_v1_mnar_interpolation <- interpolate_missing(p10_v1_mnar)
p10_v2_mnar_interpolation <- interpolate_missing(p10_v2_mnar)

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

#MNAR
#p1
p1_v1_mnar_kalman <- kalman_imputation(p1_v1_mnar)
p1_v2_mnar_kalman <- kalman_imputation(p1_v2_mnar)
#2
p2_v1_mnar_kalman <- kalman_imputation(p2_v1_mnar)
p2_v2_mnar_kalman <- kalman_imputation(p2_v2_mnar)
#3
p3_v1_mnar_kalman <- kalman_imputation(p3_v1_mnar)
p3_v2_mnar_kalman <- kalman_imputation(p3_v2_mnar)
#4
p4_v1_mnar_kalman <- kalman_imputation(p4_v1_mnar)
p4_v2_mnar_kalman <- kalman_imputation(p4_v2_mnar)
#5
p5_v1_mnar_kalman <- kalman_imputation(p5_v1_mnar)
p5_v2_mnar_kalman <- kalman_imputation(p5_v2_mnar)
#6
p6_v1_mnar_kalman <- kalman_imputation(p6_v1_mnar)
p6_v2_mnar_kalman <- kalman_imputation(p6_v2_mnar)
#7
p7_v1_mnar_kalman <- kalman_imputation(p7_v1_mnar)
p7_v2_mnar_kalman <- kalman_imputation(p7_v2_mnar)
#8
p8_v1_mnar_kalman <- kalman_imputation(p8_v1_mnar)
p8_v2_mnar_kalman <- kalman_imputation(p8_v2_mnar)
#9
p9_v1_mnar_kalman <- kalman_imputation(p9_v1_mnar)
p9_v2_mnar_kalman <- kalman_imputation(p9_v2_mnar)
#10
p10_v1_mnar_kalman <- kalman_imputation(p10_v1_mnar)
p10_v2_mnar_kalman <- kalman_imputation(p10_v2_mnar)

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

#MNAR
#p1
p1_v1_mnar_lwma <- weighted_mov_average(p1_v1_mnar)
p1_v2_mnar_lwma <- weighted_mov_average(p1_v2_mnar)
#2
p2_v1_mnar_lwma <- weighted_mov_average(p2_v1_mnar)
p2_v2_mnar_lwma <- weighted_mov_average(p2_v2_mnar)
#3
p3_v1_mnar_lwma <- weighted_mov_average(p3_v1_mnar)
p3_v2_mnar_lwma <- weighted_mov_average(p3_v2_mnar)
#4
p4_v1_mnar_lwma <- weighted_mov_average(p4_v1_mnar)
p4_v2_mnar_lwma <- weighted_mov_average(p4_v2_mnar)
#5
p5_v1_mnar_lwma <- weighted_mov_average(p5_v1_mnar)
p5_v2_mnar_lwma <- weighted_mov_average(p5_v2_mnar)
#6
p6_v1_mnar_lwma <- weighted_mov_average(p6_v1_mnar)
p6_v2_mnar_lwma <- weighted_mov_average(p6_v2_mnar)
#7
p7_v1_mnar_lwma <- weighted_mov_average(p7_v1_mnar)
p7_v2_mnar_lwma <- weighted_mov_average(p7_v2_mnar)
#8
p8_v1_mnar_lwma <- weighted_mov_average(p8_v1_mnar)
p8_v2_mnar_lwma <- weighted_mov_average(p8_v2_mnar)
#9
p9_v1_mnar_lwma <- weighted_mov_average(p9_v1_mnar)
p9_v2_mnar_lwma <- weighted_mov_average(p9_v2_mnar)
#10
p10_v1_mnar_lwma <- weighted_mov_average(p10_v1_mnar)
p10_v2_mnar_lwma <- weighted_mov_average(p10_v2_mnar)


# --------------------------------------
# TITLE: KINETICS PLOT BEFORE AND AFTER
# --------------------------------------
# -------------------------------
# Plot before and After 
# -------------------------------

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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Interpolation/MCAR_Interpolation_1MV.pdf", width = 14, height = 10)

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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Interpolation/MNAR_Interpolation_1MV.pdf", width = 14, height = 10)

#MNAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p1_visit2, p1_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p2_visit2, p2_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p3_visit2, p3_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p4_visit2, p4_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p5_visit2, p5_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p6_visit2, p6_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p7_visit2, p7_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p8_visit2, p8_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p9_visit2, p9_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p10_visit2, p10_v2_mnar_interpolation, visit = "Visit 2", type = "MNAR")

dev.off()

# ----------------------------
# Part 2: Kalman Smoothing
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Kalman_Smoothing/MCAR_Kalman_1MV.pdf", width = 14, height = 10)

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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Kalman_Smoothing/MNAR_Kalman_1MV.pdf", width = 14, height = 10)

#MNAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p1_visit2, p1_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p2_visit2, p2_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p3_visit2, p3_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p4_visit2, p4_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p5_visit2, p5_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p6_visit2, p6_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p7_visit2, p7_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p8_visit2, p8_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p9_visit2, p9_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p10_visit2, p10_v2_mnar_kalman, visit = "Visit 2", type = "MNAR")

dev.off()

# ----------------------------
# Part 3: LWMA
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/LWMA/MCAR_LWMA_1MV.pdf", width = 14, height = 10)

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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/LWMA/MNAR_LWMA_1MV.pdf", width = 14, height = 10)

#MNAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p1_visit2, p1_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p2_visit2, p2_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p3_visit2, p3_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p4_visit2, p4_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p5_visit2, p5_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p6_visit2, p6_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p7_visit2, p7_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p8_visit2, p8_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p9_visit2, p9_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p10_visit2, p10_v2_mnar_lwma, visit = "Visit 2", type = "MNAR")

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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Interpolation/MCAR_Interpolation_1MV_NRMSE.pdf", width = 14, height = 10)

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

# -------------------------
# Part 1.2: NRMSE (MNAR)
# --------------------------

#call function to calcualte nrms
#interpolation
#p1
nrmse_interp_p1v1_mnar <- calculate_nrsme(p1_visit1, p1_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p1v2_mnar <- calculate_nrsme(p1_visit2, p1_v2_mnar_interpolation, method = "Linear Interpolation")
#p2
nrmse_interp_p2v1_mnar <- calculate_nrsme(p2_visit1, p2_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p2v2_mnar <- calculate_nrsme(p2_visit2, p2_v2_mnar_interpolation, method = "Linear Interpolation")
#p3
nrmse_interp_p3v1_mnar <- calculate_nrsme(p3_visit1, p3_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p3v2_mnar <- calculate_nrsme(p3_visit2, p3_v2_mnar_interpolation, method = "Linear Interpolation")
#p4
nrmse_interp_p4v1_mnar <- calculate_nrsme(p4_visit1, p4_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p4v2_mnar <- calculate_nrsme(p4_visit2, p4_v2_mnar_interpolation, method = "Linear Interpolation")
#p5
nrmse_interp_p5v1_mnar <- calculate_nrsme(p5_visit1, p5_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p5v2_mnar <- calculate_nrsme(p5_visit2, p5_v2_mnar_interpolation, method = "Linear Interpolation")
#p6
nrmse_interp_p6v1_mnar <- calculate_nrsme(p6_visit1, p6_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p6v2_mnar <- calculate_nrsme(p6_visit2, p6_v2_mnar_interpolation, method = "Linear Interpolation")
#p7
nrmse_interp_p7v1_mnar <- calculate_nrsme(p7_visit1, p7_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p7v2_mnar <- calculate_nrsme(p7_visit2, p7_v2_mnar_interpolation, method = "Linear Interpolation")
#p8
nrmse_interp_p8v1_mnar <- calculate_nrsme(p8_visit1, p8_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p8v2_mnar <- calculate_nrsme(p8_visit2, p8_v2_mnar_interpolation, method = "Linear Interpolation")
#p9
nrmse_interp_p9v1_mnar <- calculate_nrsme(p9_visit1, p9_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p9v2_mnar <- calculate_nrsme(p9_visit2, p9_v2_mnar_interpolation, method = "Linear Interpolation")
#p10
nrmse_interp_p10v1_mnar <- calculate_nrsme(p10_visit1, p10_v1_mnar_interpolation, method = "Linear Interpolation")
nrmse_interp_p10v2_mnar <- calculate_nrsme(p10_visit2, p10_v2_mnar_interpolation, method = "Linear Interpolation")

#combine visit 1 
nrmse_mnar_visit1 <- bind_rows(
  nrmse_interp_p1v1_mnar %>% mutate(Patient = "P1"),
  nrmse_interp_p2v1_mnar %>% mutate(Patient = "P2"),
  nrmse_interp_p3v1_mnar %>% mutate(Patient = "P3"),
  nrmse_interp_p4v1_mnar %>% mutate(Patient = "P4"),
  nrmse_interp_p5v1_mnar %>% mutate(Patient = "P5"),
  nrmse_interp_p6v1_mnar %>% mutate(Patient = "P6"),
  nrmse_interp_p7v1_mnar %>% mutate(Patient = "P7"),
  nrmse_interp_p8v1_mnar %>% mutate(Patient = "P8"),
  nrmse_interp_p9v1_mnar %>% mutate(Patient = "P9"),
  nrmse_interp_p10v1_mnar %>% mutate(Patient = "P10")
)

#combine visit 2
nrmse_mnar_visit2 <- bind_rows(
  nrmse_interp_p1v2_mnar %>% mutate(Patient = "P1"),
  nrmse_interp_p2v2_mnar %>% mutate(Patient = "P2"),
  nrmse_interp_p3v2_mnar %>% mutate(Patient = "P3"),
  nrmse_interp_p4v2_mnar %>% mutate(Patient = "P4"),
  nrmse_interp_p5v2_mnar %>% mutate(Patient = "P5"),
  nrmse_interp_p6v2_mnar %>% mutate(Patient = "P6"),
  nrmse_interp_p7v2_mnar %>% mutate(Patient = "P7"),
  nrmse_interp_p8v2_mnar %>% mutate(Patient = "P8"),
  nrmse_interp_p9v2_mnar %>% mutate(Patient = "P9"),
  nrmse_interp_p10v2_mnar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Interpolation/MNAR_Interpolation_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1 (MNAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2 (MNAR 1 MV in middle)",
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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Kalman_Smoothing/MCAR_Kalman_1MV_NRMSE.pdf", width = 14, height = 10)

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


# --------------------------
# Part 2.2: NRMSE (MNAR)
# --------------------------

#call function to calcualte nrms
#kalman
#p1
nrmse_kalman_p1v1_mnar <- calculate_nrsme(p1_visit1, p1_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p1v2_mnar <- calculate_nrsme(p1_visit2, p1_v2_mnar_kalman, method = "Kalman Smoothing")
#p2
nrmse_kalman_p2v1_mnar <- calculate_nrsme(p2_visit1, p2_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p2v2_mnar <- calculate_nrsme(p2_visit2, p2_v2_mnar_kalman, method = "Kalman Smoothing")
#p3
nrmse_kalman_p3v1_mnar <- calculate_nrsme(p3_visit1, p3_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p3v2_mnar <- calculate_nrsme(p3_visit2, p3_v2_mnar_kalman, method = "Kalman Smoothing")
#p4
nrmse_kalman_p4v1_mnar <- calculate_nrsme(p4_visit1, p4_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p4v2_mnar <- calculate_nrsme(p4_visit2, p4_v2_mnar_kalman, method = "Kalman Smoothing")
#p5
nrmse_kalman_p5v1_mnar <- calculate_nrsme(p5_visit1, p5_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p5v2_mnar <- calculate_nrsme(p5_visit2, p5_v2_mnar_kalman, method = "Kalman Smoothing")
#p6
nrmse_kalman_p6v1_mnar <- calculate_nrsme(p6_visit1, p6_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p6v2_mnar <- calculate_nrsme(p6_visit2, p6_v2_mnar_kalman, method = "Kalman Smoothing")
#p7
nrmse_kalman_p7v1_mnar <- calculate_nrsme(p7_visit1, p7_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p7v2_mnar <- calculate_nrsme(p7_visit2, p7_v2_mnar_kalman, method = "Kalman Smoothing")
#p8
nrmse_kalman_p8v1_mnar <- calculate_nrsme(p8_visit1, p8_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p8v2_mnar <- calculate_nrsme(p8_visit2, p8_v2_mnar_kalman, method = "Kalman Smoothing")
#p9
nrmse_kalman_p9v1_mnar <- calculate_nrsme(p9_visit1, p9_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p9v2_mnar <- calculate_nrsme(p9_visit2, p9_v2_mnar_kalman, method = "Kalman Smoothing")
#p10
nrmse_kalman_p10v1_mnar <- calculate_nrsme(p10_visit1, p10_v1_mnar_kalman, method = "Kalman Smoothing")
nrmse_kalman_p10v2_mnar <- calculate_nrsme(p10_visit2, p10_v2_mnar_kalman, method = "Kalman Smoothing")

#combine visit 1 
nrmse_kalman_mnar_visit1 <- bind_rows(
  nrmse_kalman_p1v1_mnar %>% mutate(Patient = "P1"),
  nrmse_kalman_p2v1_mnar %>% mutate(Patient = "P2"),
  nrmse_kalman_p3v1_mnar %>% mutate(Patient = "P3"),
  nrmse_kalman_p4v1_mnar %>% mutate(Patient = "P4"),
  nrmse_kalman_p5v1_mnar %>% mutate(Patient = "P5"),
  nrmse_kalman_p6v1_mnar %>% mutate(Patient = "P6"),
  nrmse_kalman_p7v1_mnar %>% mutate(Patient = "P7"),
  nrmse_kalman_p8v1_mnar %>% mutate(Patient = "P8"),
  nrmse_kalman_p9v1_mnar %>% mutate(Patient = "P9"),
  nrmse_kalman_p10v1_mnar %>% mutate(Patient = "P10")
)

#combine visit 2
nrmse_kalman_mnar_visit2 <- bind_rows(
  nrmse_kalman_p1v2_mnar %>% mutate(Patient = "P1"),
  nrmse_kalman_p2v2_mnar %>% mutate(Patient = "P2"),
  nrmse_kalman_p3v2_mnar %>% mutate(Patient = "P3"),
  nrmse_kalman_p4v2_mnar %>% mutate(Patient = "P4"),
  nrmse_kalman_p5v2_mnar %>% mutate(Patient = "P5"),
  nrmse_kalman_p6v2_mnar %>% mutate(Patient = "P6"),
  nrmse_kalman_p7v2_mnar %>% mutate(Patient = "P7"),
  nrmse_kalman_p8v2_mnar %>% mutate(Patient = "P8"),
  nrmse_kalman_p9v2_mnar %>% mutate(Patient = "P9"),
  nrmse_kalman_p10v2_mnar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/Kalman_Smoothing/MNAR_Kalman_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_kalman_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1 (MNAR 1 MV in middle)",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_kalman_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2 (MNAR 1 MV in middle)",
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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/LWMA/MCAR_LWMA_1MV_NRMSE.pdf", width = 14, height = 10)

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

# --------------------------
# Part 3.2: NRMSE (MNAR)
# --------------------------

#call function to calcualte nrms
#LWMA
#p1
nrmse_lwma_p1v1_mnar <- calculate_nrsme(p1_visit1, p1_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p1v2_mnar <- calculate_nrsme(p1_visit2, p1_v2_mnar_lwma, method = "LWMA")
#p2
nrmse_lwma_p2v1_mnar <- calculate_nrsme(p2_visit1, p2_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p2v2_mnar <- calculate_nrsme(p2_visit2, p2_v2_mnar_lwma, method = "LWMA")
#p3
nrmse_lwma_p3v1_mnar <- calculate_nrsme(p3_visit1, p3_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p3v2_mnar <- calculate_nrsme(p3_visit2, p3_v2_mnar_lwma, method = "LWMA")
#p4
nrmse_lwma_p4v1_mnar <- calculate_nrsme(p4_visit1, p4_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p4v2_mnar <- calculate_nrsme(p4_visit2, p4_v2_mnar_lwma, method = "LWMA")
#p5
nrmse_lwma_p5v1_mnar <- calculate_nrsme(p5_visit1, p5_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p5v2_mnar <- calculate_nrsme(p5_visit2, p5_v2_mnar_lwma, method = "LWMA")
#p6
nrmse_lwma_p6v1_mnar <- calculate_nrsme(p6_visit1, p6_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p6v2_mnar <- calculate_nrsme(p6_visit2, p6_v2_mnar_lwma, method = "LWMA")
#p7
nrmse_lwma_p7v1_mnar <- calculate_nrsme(p7_visit1, p7_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p7v2_mnar <- calculate_nrsme(p7_visit2, p7_v2_mnar_lwma, method = "LWMA")
#p8
nrmse_lwma_p8v1_mnar <- calculate_nrsme(p8_visit1, p8_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p8v2_mnar <- calculate_nrsme(p8_visit2, p8_v2_mnar_lwma, method = "LWMA")
#p9
nrmse_lwma_p9v1_mnar <- calculate_nrsme(p9_visit1, p9_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p9v2_mnar <- calculate_nrsme(p9_visit2, p9_v2_mnar_lwma, method = "LWMA")
#p10
nrmse_lwma_p10v1_mnar <- calculate_nrsme(p10_visit1, p10_v1_mnar_lwma, method = "LWMA")
nrmse_lwma_p10v2_mnar <- calculate_nrsme(p10_visit2, p10_v2_mnar_lwma, method = "LWMA")

#combine visit 1 
nrmse_lwma_mnar_visit1 <- bind_rows(
  nrmse_lwma_p1v1_mnar %>% mutate(Patient = "P1"),
  nrmse_lwma_p2v1_mnar %>% mutate(Patient = "P2"),
  nrmse_lwma_p3v1_mnar %>% mutate(Patient = "P3"),
  nrmse_lwma_p4v1_mnar %>% mutate(Patient = "P4"),
  nrmse_lwma_p5v1_mnar %>% mutate(Patient = "P5"),
  nrmse_lwma_p6v1_mnar %>% mutate(Patient = "P6"),
  nrmse_lwma_p7v1_mnar %>% mutate(Patient = "P7"),
  nrmse_lwma_p8v1_mnar %>% mutate(Patient = "P8"),
  nrmse_lwma_p9v1_mnar %>% mutate(Patient = "P9"),
  nrmse_lwma_p10v1_mnar %>% mutate(Patient = "P10")
)

#combine visit 2
nrmse_lwma_mnar_visit2 <- bind_rows(
  nrmse_lwma_p1v2_mnar %>% mutate(Patient = "P1"),
  nrmse_lwma_p2v2_mnar %>% mutate(Patient = "P2"),
  nrmse_lwma_p3v2_mnar %>% mutate(Patient = "P3"),
  nrmse_lwma_p4v2_mnar %>% mutate(Patient = "P4"),
  nrmse_lwma_p5v2_mnar %>% mutate(Patient = "P5"),
  nrmse_lwma_p6v2_mnar %>% mutate(Patient = "P6"),
  nrmse_lwma_p7v2_mnar %>% mutate(Patient = "P7"),
  nrmse_lwma_p8v2_mnar %>% mutate(Patient = "P8"),
  nrmse_lwma_p9v2_mnar %>% mutate(Patient = "P9"),
  nrmse_lwma_p10v2_mnar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/LWMA/MNAR_LWMA_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_lwma_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1 (MNAR 1 MV)",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_lwma_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2 (MNAR 1 MV)",
       y = "NRMSE", x = "Patient")

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

# ------------------------
# Part 1: Interpolation AUC
# ------------------------

#MCAR
#call function for interpolation 
#p1
auc_mcar_p1v1_original <- calculate_auc(p1_visit1)
auc_mcar_p1v1_interpolation <- calculate_auc(p1_v1_mcar_interpolation)
auc_mcar_p1v2_original <- calculate_auc(p1_visit2)
auc_mcar_p1v2_interpolation <- calculate_auc(p1_v2_mcar_interpolation)

#combine into dataframe for comparison
#p1
auc_mcar_combined_p1v1 <- data.frame(
  Metabolite = names(auc_mcar_p1v1_original),
  Original_AUC = auc_mcar_p1v1_original,
  Interpolation_AUC = auc_mcar_p1v1_interpolation,
  Difference = auc_mcar_p1v1_interpolation - auc_mcar_p1v1_original,
  Percent_Change_AUC = 100 *(auc_mcar_p1v1_interpolation - auc_mcar_p1v1_original) / auc_mcar_p1v1_original
)
auc_mcar_combined_p1v2 <- data.frame(
  Metabolite = names(auc_mcar_p1v2_original),
  Original_AUC = auc_mcar_p1v2_original,
  Interpolation_AUC = auc_mcar_p1v2_interpolation,
  Difference = auc_mcar_p1v2_interpolation - auc_mcar_p1v2_original,
  Percent_Change_AUC = 100 *(auc_mcar_p1v2_interpolation - auc_mcar_p1v2_original) / auc_mcar_p1v2_original
)

#plot
#p1
ggplot(auc_mcar_combined_p1v1, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(auc_mcar_combined_p1v2, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 

#MNAR
#call function for interpolation 
#p1
auc_mnar_p1v1_original <- calculate_auc(p1_visit1)
auc_mnar_p1v1_interpolation <- calculate_auc(p1_v1_mnar_interpolation)
auc_mnar_p1v2_original <- calculate_auc(p1_visit2)
auc_mnar_p1v2_interpolation <- calculate_auc(p1_v2_mnar_interpolation)

#combine into dataframe for comparison
#p1
auc_mnar_combined_p1v1 <- data.frame(
  Metabolite = names(auc_mnar_p1v1_original),
  Original_AUC = auc_mnar_p1v1_original,
  Interpolation_AUC = auc_mnar_p1v1_interpolation,
  Difference = auc_mnar_p1v1_interpolation - auc_mnar_p1v1_original,
  Percent_Change_AUC = 100 *(auc_mnar_p1v1_interpolation - auc_mnar_p1v1_original) / auc_mnar_p1v1_original
)

auc_mnar_combined_p1v2 <- data.frame(
  Metabolite = names(auc_mnar_p1v2_original),
  Original_AUC = auc_mnar_p1v2_original,
  Interpolation_AUC = auc_mnar_p1v2_interpolation,
  Difference = auc_mnar_p1v2_interpolation - auc_mnar_p1v2_original,
  Percent_Change_AUC = 100 *(auc_mnar_p1v2_interpolation - auc_mnar_p1v2_original) / auc_mnar_p1v2_original
)

#plot
#p1
ggplot(auc_mnar_combined_p1v1, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(auc_mnar_combined_p1v2, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ------------------------
# Part 2: Kalman AUC
# ------------------------

#MCAR
#call function for interpolation 
#p1
auc_mcar_p1v1_original <- calculate_auc(p1_visit1)
auc_mcar_p1v1_kalman <- calculate_auc(p1_v1_mcar_kalman)
auc_mcar_p1v2_original <- calculate_auc(p1_visit2)
auc_mcar_p1v2_kalman <- calculate_auc(p1_v2_mcar_kalman)

#combine into dataframe for comparison
#p1
auc_mcar_combined_p1v1_kalman <- data.frame(
  Metabolite = names(auc_mcar_p1v1_original),
  Original_AUC = auc_mcar_p1v1_original,
  Interpolation_AUC = auc_mcar_p1v1_kalman,
  Difference = auc_mcar_p1v1_kalman - auc_mcar_p1v1_original,
  Percent_Change_AUC = 100 *(auc_mcar_p1v1_kalman - auc_mcar_p1v1_original) / auc_mcar_p1v1_original
)
auc_mcar_combined_p1v2_kalman <- data.frame(
  Metabolite = names(auc_mcar_p1v2_original),
  Original_AUC = auc_mcar_p1v2_original,
  Interpolation_AUC = auc_mcar_p1v2_kalman,
  Difference = auc_mcar_p1v2_kalman - auc_mcar_p1v2_original,
  Percent_Change_AUC = 100 *(auc_mcar_p1v2_kalman - auc_mcar_p1v2_original) / auc_mcar_p1v2_original
)

#plot
#p1
ggplot(auc_mcar_combined_p1v1_kalman, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(auc_mcar_combined_p1v2_kalman, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#MNAR
#call function for interpolation 
#p1
auc_mnar_p1v1_original <- calculate_auc(p1_visit1)
auc_mnar_p1v1_kalman <- calculate_auc(p1_v1_mnar_kalman)
auc_mnar_p1v2_original <- calculate_auc(p1_visit2)
auc_mnar_p1v2_kalman <- calculate_auc(p1_v2_mcar_kalman)

#combine into dataframe for comparison
#p1
auc_mnar_combined_p1v1_kalman <- data.frame(
  Metabolite = names(auc_mnar_p1v1_original),
  Original_AUC = auc_mnar_p1v1_original,
  Interpolation_AUC = auc_mnar_p1v1_kalman,
  Difference = auc_mnar_p1v1_kalman - auc_mnar_p1v1_original,
  Percent_Change_AUC = 100 *(auc_mnar_p1v1_kalman - auc_mnar_p1v1_original) / auc_mnar_p1v1_original
)

auc_mnar_combined_p1v2_kalman <- data.frame(
  Metabolite = names(auc_mnar_p1v2_original),
  Original_AUC = auc_mnar_p1v2_original,
  Interpolation_AUC = auc_mnar_p1v2_kalman,
  Difference = auc_mnar_p1v2_kalman - auc_mnar_p1v2_original,
  Percent_Change_AUC = 100 *(auc_mnar_p1v2_kalman - auc_mnar_p1v2_original) / auc_mnar_p1v2_original
)

#plot
#p1
ggplot(auc_mnar_combined_p1v1_kalman, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(auc_mnar_combined_p1v2_kalman, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ------------------------
# Part 3: LWMA AUC
# ------------------------

#MCAR
#call function for interpolation 
#p1
auc_mcar_p1v1_original <- calculate_auc(p1_visit1)
auc_mcar_p1v1_lwma <- calculate_auc(p1_v1_mcar_lwma)
auc_mcar_p1v2_original <- calculate_auc(p1_visit2)
auc_mcar_p1v2_lwma <- calculate_auc(p1_v2_mcar_lwma)

#combine into dataframe for comparison
#p1
auc_mcar_combined_p1v1_lwma <- data.frame(
  Metabolite = names(auc_mcar_p1v1_original),
  Original_AUC = auc_mcar_p1v1_original,
  Interpolation_AUC = auc_mcar_p1v1_lwma,
  Difference = auc_mcar_p1v1_lwma - auc_mcar_p1v1_original,
  Percent_Change_AUC = 100 *(auc_mcar_p1v1_lwma - auc_mcar_p1v1_original) / auc_mcar_p1v1_original
)
auc_mcar_combined_p1v2_lwma <- data.frame(
  Metabolite = names(auc_mcar_p1v2_original),
  Original_AUC = auc_mcar_p1v2_original,
  Interpolation_AUC = auc_mcar_p1v2_lwma,
  Difference = auc_mcar_p1v2_lwma - auc_mcar_p1v2_original,
  Percent_Change_AUC = 100 *(auc_mcar_p1v2_lwma - auc_mcar_p1v2_original) / auc_mcar_p1v2_original
)

#plot
#p1
ggplot(auc_mcar_combined_p1v1_lwma, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(auc_mcar_combined_p1v2_lwma, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#MNAR
#call function for interpolation 
#p1
auc_mnar_p1v1_original <- calculate_auc(p1_visit1)
auc_mnar_p1v1_lwma <- calculate_auc(p1_v1_mnar_lwma)
auc_mnar_p1v2_original <- calculate_auc(p1_visit2)
auc_mnar_p1v2_lwma <- calculate_auc(p1_v2_mcar_lwma)

#combine into dataframe for comparison
#p1
auc_mnar_combined_p1v1_lwma <- data.frame(
  Metabolite = names(auc_mnar_p1v1_original),
  Original_AUC = auc_mnar_p1v1_original,
  Interpolation_AUC = auc_mnar_p1v1_lwma,
  Difference = auc_mnar_p1v1_lwma - auc_mnar_p1v1_original,
  Percent_Change_AUC = 100 *(auc_mnar_p1v1_lwma - auc_mnar_p1v1_original) / auc_mnar_p1v1_original
)

auc_mnar_combined_p1v2_lwma <- data.frame(
  Metabolite = names(auc_mnar_p1v2_original),
  Original_AUC = auc_mnar_p1v2_original,
  Interpolation_AUC = auc_mnar_p1v2_lwma,
  Difference = auc_mnar_p1v2_lwma - auc_mnar_p1v2_original,
  Percent_Change_AUC = 100 *(auc_mnar_p1v2_lwma - auc_mnar_p1v2_original) / auc_mnar_p1v2_original
)

#plot
#p1
ggplot(auc_mnar_combined_p1v1_lwma, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(auc_mnar_combined_p1v2_lwma, aes(x = Metabolite, y = Percent_Change_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  labs(
    title = "Percent Change in AUC after Interpolation (MCAR)",
    x = "Metabolite",
    y = "Percent Change in AUC (%)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))














