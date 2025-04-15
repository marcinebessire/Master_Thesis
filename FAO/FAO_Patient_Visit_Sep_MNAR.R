#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(zoo) #for interpolation
library(imputeTS) #for imputation methods
library(pracma) #for AUC calculation
library(missForest) #for RF
library(keras) #for LSTM
library(vegan) #for pca


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

# -------------------------------
# Part 2: Patient 4, missing t=120
# -------------------------------

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
p4_visit1_complete <- add_row_120(p4_visit1)

#assume linear interpolation for this time point 
p4_visit1_full <- interpolate_missing(p4_visit1_complete)

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
p4_v1_mnar <- MNAR_manipulation_lowest(p4_visit1_full)
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
interpolate_missing <- function(data) {
  data_copy <- data
  
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[, 6:ncol(data_copy)], function(col) {
    if (is.numeric(col)) {
      #interpolate and extrapolate both ends if needed
      return(na.approx(col, na.rm = FALSE, rule = 2)) #rule 2 allows extrapolation on both ends (otherwise stays NA)
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}


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
      return(na_ma(col, k = window, weighting = "exponential")) #give more value to the closer points and linearly decreasing
    } else {
      return(col)
    }
  })
  return(data_copy)
}

#call function for LWMA

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


# --------------------------------
# Part 4: LOESS + RF
# --------------------------------

# #LOESS (locally estimated scatterplot smoothing)
# #LOESS + Random Forest Imputation
# impute_loess_then_rf <- function(df, time_col = "Time_min", sd_threshold = 5) {
#   df_imputed <- df
#   metabolite_cols <- names(df)[6:ncol(df)]
#   
#   for (metabolite in metabolite_cols) {
#     message("LOESS for: ", metabolite)
#     
#     time <- df[[time_col]]
#     y <- df[[metabolite]]
#     na_indices <- which(is.na(y))
#     if (length(na_indices) == 0) next
#     
#     #require at least 3 non-missing values to fit LOESS
#     if (sum(!is.na(y)) < 3) {
#       message("-> Skipping LOESS: not enough non-missing values.")
#       next
#     }
#     
#     #calculate variability and choose span
#     variability <- sd(y, na.rm = TRUE)
#     span <- max(0.6, ifelse(variability < sd_threshold, 0.4, 0.75))
#     message("  -> Using span = ", span, " (SD = ", round(variability, 2), ")")
#     
#     #fit LOESS
#     df_non_na <- data.frame(time = time[!is.na(y)], y = y[!is.na(y)])
#     loess_fit <- tryCatch({ #tries to fit first-degress loess model (linear local fit)
#       loess(y ~ time, data = df_non_na, span = span, degree = 1, control = loess.control(surface = "direct"))
#     }, error = function(e) {
#       warning("  -> LOESS failed for ", metabolite, ": ", e$message)
#       return(NULL)
#     })
#     
#     #skip if model failed
#     if (is.null(loess_fit)) next
#     
#     #predict only for values within the fitting range
#     for (na_index in na_indices) {
#       t_missing <- time[na_index]
#       
#       #only predict if within time range (no extrapolation outsie)
#       if (t_missing >= min(df_non_na$time) && t_missing <= max(df_non_na$time)) {
#         predicted <- predict(loess_fit, newdata = data.frame(time = t_missing))
#         if (!is.na(predicted)) {
#           df_imputed[[metabolite]][na_index] <- predicted
#         } else {
#           message("  -> LOESS could not predict at time = ", t_missing)
#         }
#       } else {
#         message("  -> Time = ", t_missing, " is outside LOESS fitting range.")
#       }
#     }
#   }
#   
#   #random Forest to fill in any remaining NA
#   message("Running Random Forest refinement with missForest...")
#   rf_data <- df_imputed[, metabolite_cols]
#   rf_imputed <- missForest(rf_data)$ximp
#   df_imputed[, metabolite_cols] <- rf_imputed
#   
#   return(df_imputed)
# }

#LOESS + Random Forest Imputation with Log Transform
impute_loess_then_rf <- function(df, time_col = "Time_min", sd_threshold = 5) {
  df_imputed <- df
  metabolite_cols <- names(df)[6:ncol(df)]
  
  #apply log transform (avoid log(0))
  log_transformed <- df
  log_transformed[, metabolite_cols] <- log(df[, metabolite_cols] + 1)
  
  for (metabolite in metabolite_cols) {
    message("LOESS for: ", metabolite)
    
    time <- df[[time_col]]
    y <- log_transformed[[metabolite]]
    na_indices <- which(is.na(y))
    if (length(na_indices) == 0) next
    
    #require at least 3 non-missing values to fit LOESS
    if (sum(!is.na(y)) < 3) {
      message("  -> Skipping LOESS: not enough non-missing values.")
      next
    }
    
    #calculate variability and choose span
    variability <- sd(y, na.rm = TRUE)
    span <- max(0.6, ifelse(variability < sd_threshold, 0.4, 0.75))
    message("  -> Using span = ", span, " (SD = ", round(variability, 2), ")")
    
    #fit LOESS model
    df_non_na <- data.frame(time = time[!is.na(y)], y = y[!is.na(y)])
    loess_fit <- tryCatch({
      loess(y ~ time, data = df_non_na, span = span, degree = 1,
            control = loess.control(surface = "direct"))
    }, error = function(e) {
      warning("  -> LOESS failed for ", metabolite, ": ", e$message)
      return(NULL)
    })
    
    #skip if LOESS model failed
    if (is.null(loess_fit)) next
    
    #predict and fill NA values within LOESS range
    for (na_index in na_indices) {
      t_missing <- time[na_index]
      if (t_missing >= min(df_non_na$time) && t_missing <= max(df_non_na$time)) {
        predicted <- predict(loess_fit, newdata = data.frame(time = t_missing))
        if (!is.na(predicted)) {
          log_transformed[[metabolite]][na_index] <- predicted
        } else {
          message("  -> LOESS could not predict at time = ", t_missing)
        }
      } else {
        message("  -> Time = ", t_missing, " is outside LOESS fitting range.")
      }
    }
  }
  
  #RF on log-transformed data
  message("Running Random Forest refinement with missForest...")
  rf_data <- log_transformed[, metabolite_cols]
  rf_imputed <- missForest(rf_data)$ximp
  log_transformed[, metabolite_cols] <- rf_imputed
  
  #back-transform (inverse log)
  df_imputed[, metabolite_cols] <- exp(log_transformed[, metabolite_cols]) - 1
  
  return(df_imputed)
}


#call function for loess + rf imputation
#visit 1
p1_v1_loess <- impute_loess_then_rf(p1_v1_mnar)
p2_v1_loess <- impute_loess_then_rf(p2_v1_mnar)
p3_v1_loess <- impute_loess_then_rf(p3_v1_mnar)
p4_v1_loess <- impute_loess_then_rf(p4_v1_mnar)
p5_v1_loess <- impute_loess_then_rf(p5_v1_mnar)
p6_v1_loess <- impute_loess_then_rf(p6_v1_mnar)
p7_v1_loess <- impute_loess_then_rf(p7_v1_mnar)
p8_v1_loess <- impute_loess_then_rf(p8_v1_mnar)
p9_v1_loess <- impute_loess_then_rf(p9_v1_mnar)
p10_v1_loess <- impute_loess_then_rf(p10_v1_mnar)
#visit 2
p1_v2_loess <- impute_loess_then_rf(p1_v2_mnar)
p2_v2_loess <- impute_loess_then_rf(p2_v2_mnar)
p3_v2_loess <- impute_loess_then_rf(p3_v2_mnar)
p4_v2_loess <- impute_loess_then_rf(p4_v2_mnar)
p5_v2_loess <- impute_loess_then_rf(p5_v2_mnar)
p6_v2_loess <- impute_loess_then_rf(p6_v2_mnar)
p7_v2_loess <- impute_loess_then_rf(p7_v2_mnar)
p8_v2_loess <- impute_loess_then_rf(p8_v2_mnar)
p9_v2_loess <- impute_loess_then_rf(p9_v2_mnar)
p10_v2_loess <- impute_loess_then_rf(p10_v2_mnar)

# --------------------------------
# Part 5: LSTM
# --------------------------------

#step 1:create array/list of the dataframes 

#combine into one array
#visit1
all_mnar_list_v1 <- list(
  p1_v1_mnar, p2_v1_mnar, 
  p3_v1_mnar, p4_v1_mnar,
  p5_v1_mnar, p6_v1_mnar,
  p7_v1_mnar, p8_v1_mnar,
  p9_v1_mnar, p10_v1_mnar
)
#visit 2
all_mnar_list_v2 <- list(
  p1_v2_mnar, p2_v2_mnar, 
  p3_v2_mnar, p4_v2_mnar,
  p5_v2_mnar, p6_v2_mnar,
  p7_v2_mnar, p8_v2_mnar,
  p9_v2_mnar, p10_v2_mnar
)

#step 2:create function for LSTM imputation
impute_with_lstm_all_metabolites <- function(mnar_df_list, meta_start_col = 6, n_epochs = 200, batch_size = 4, verbose = 0) {
  #get metabolite columns
  metabolite_cols <- colnames(mnar_df_list[[1]])[meta_start_col:ncol(mnar_df_list[[1]])]
  
  #copy input list to update
  updated_list <- mnar_df_list
  
  for (metabolite in metabolite_cols) {
    message("Running LSTM imputation for: ", metabolite)
    
    #extract the sequences for the current metabolite
    sequences <- lapply(mnar_df_list, function(df) df[[metabolite]])
    sequences <- do.call(rbind, sequences)  #shape: patients x timepoints
    
    #replace NA with 0 for masking later
    X <- sequences
    X[is.na(X)] <- 0  #0 will be masked
    
    #reshape to 3D [samples, timesteps, features]
    X_array <- array(X, dim = c(nrow(X), ncol(X), 1))
    Y_array <- X_array  #autoencoder style: predict the full sequence
    
    #define LSTM model
    #bidericteion + deeper LSTM model
    model <- keras_model_sequential() %>% #linear stack of layers, each layer feeds into the next
      layer_masking(mask_value = 0, input_shape = c(ncol(X), 1)) %>% #tells the LTM to ignore (mask) any time step with value 0
      bidirectional(layer_lstm( #captures bothpast and future context
        units = 128,#larger capacity to model complex trends (128 memory cells per direction)
        return_sequences = TRUE, #output full sequence
        dropout = 0.3, #drop 30% of inputs (prevent overfitting)
        recurrent_dropout = 0.2 #drop 20% of recurrent state (adds regularization to memory connections)
      )) %>%
      layer_lstm( #second layer, stacked deeper on top of the first one
        units = 64, #smaller to progressively extract more compressed features
        return_sequences = TRUE,
        dropout = 0.3,
        recurrent_dropout = 0.2
      ) %>%
      layer_dense(units = 1) #single predicted value per timepoitn
    
    model %>% compile(
      loss = "mse",
      optimizer = "adam"
    )
    
    model %>% fit(
      x = X_array,
      y = Y_array,
      epochs = n_epochs,
      batch_size = batch_size,
      verbose = verbose
    )
    
    #predict missing values
    predicted <- model %>% predict(X_array)
    pred_matrix <- predicted[, , 1]  #shape: same as original
    
    #only replace the missing values
    for (i in seq_along(mnar_df_list)) {
      na_idx <- which(is.na(mnar_df_list[[i]][[metabolite]]))
      if (length(na_idx) > 0) {
        mnar_df_list[[i]][[metabolite]][na_idx] <- pred_matrix[i, na_idx]
      }
    }
    
    updated_list <- mnar_df_list  
  }
  
  return(updated_list)
}

#step 3:call function to get missing values via lstm
lstm_v1 <- impute_with_lstm_all_metabolites(all_mnar_list_v1)
lstm_v2 <- impute_with_lstm_all_metabolites(all_mnar_list_v2)

#step 4:fill those imputed values back into dataframe & create new dataframe
#lists of original MNAR dataframes for each visit
original_v1 <- list(p1_visit1, p2_visit1, p3_visit1, p4_visit1_full, p5_visit1,
                    p6_visit1, p7_visit1, p8_visit1, p9_visit1, p10_visit1)

original_v2 <- list(p1_visit2, p2_visit2, p3_visit2, p4_visit2, p5_visit2,
                    p6_visit2, p7_visit2, p8_visit2, p9_visit2, p10_visit2)

#create values names
output_names_v1 <- paste0("p", 1:10, "_v1_lstm")
output_names_v2 <- paste0("p", 1:10, "_v2_lstm")

#function to combine imputed values and original
combine_metadata_with_imputed <- function(original_df, imputed_df) {
  metadata <- original_df[, 1:5]
  metabolites <- imputed_df[, 6:ncol(imputed_df)]
  cbind(metadata, metabolites)
}

#loop to assign v1
for (i in seq_along(output_names_v1)) {
  assign(output_names_v1[i], combine_metadata_with_imputed(original_v1[[i]], lstm_v1[[i]]))
}

#loop to assign v2
for (i in seq_along(output_names_v2)) {
  assign(output_names_v2[i], combine_metadata_with_imputed(original_v2[[i]], lstm_v2[[i]]))
}

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


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/Interpolation/MNAR_Interpolation_1MV.pdf", width = 14, height = 10)

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
plot_imputed_vs_original(p4_visit1_full, p4_v1_mnar_interpolation, visit = "Visit 1", type = "MNAR")
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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/Kalman/MNAR_Kalman_1MV.pdf", width = 14, height = 10)

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
plot_imputed_vs_original(p4_visit1_full, p4_v1_mnar_kalman, visit = "Visit 1", type = "MNAR")
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


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/LWMA/MNAR_LWMA_1MV.pdf", width = 14, height = 10)

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
plot_imputed_vs_original(p4_visit1_full, p4_v1_mnar_lwma, visit = "Visit 1", type = "MNAR")
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


# ----------------------------
# Part 4: LEOSS + RF
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/LOESS_RF/MNAR_LOESS_RF_1MV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p1_visit2, p1_v2_loess, visit = "Visit 2", type = "MNAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p2_visit2, p2_v2_loess, visit = "Visit 2", type = "MNAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p3_visit2, p3_v2_loess, visit = "Visit 2", type = "MNAR")
#p4
plot_imputed_vs_original(p4_visit1_full, p4_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p4_visit2, p4_v2_loess, visit = "Visit 2", type = "MNAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p5_visit2, p5_v2_loess, visit = "Visit 2", type = "MNAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p6_visit2, p6_v2_loess, visit = "Visit 2", type = "MNAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p7_visit2, p7_v2_loess, visit = "Visit 2", type = "MNAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p8_visit2, p8_v2_loess, visit = "Visit 2", type = "MNAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p9_visit2, p9_v2_loess, visit = "Visit 2", type = "MNAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_loess, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p10_visit2, p10_v2_loess, visit = "Visit 2", type = "MNAR")

dev.off()


# ----------------------------
# Part 5: LSTM
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/LSTM/MNAR_LSTM_1MV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p1_visit2, p1_v2_lstm, visit = "Visit 2", type = "MNAR")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p2_visit2, p2_v2_lstm, visit = "Visit 2", type = "MNAR")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p3_visit2, p3_v2_lstm, visit = "Visit 2", type = "MNAR")
#p4
plot_imputed_vs_original(p4_visit1_full, p4_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p4_visit2, p4_v2_lstm, visit = "Visit 2", type = "MNAR")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p5_visit2, p5_v2_lstm, visit = "Visit 2", type = "MNAR")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p6_visit2, p6_v2_lstm, visit = "Visit 2", type = "MNAR")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p7_visit2, p7_v2_lstm, visit = "Visit 2", type = "MNAR")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p8_visit2, p8_v2_lstm, visit = "Visit 2", type = "MNAR")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p9_visit2, p9_v2_lstm, visit = "Visit 2", type = "MNAR")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_lstm, visit = "Visit 1", type = "MNAR")
plot_imputed_vs_original(p10_visit2, p10_v2_lstm, visit = "Visit 2", type = "MNAR")

dev.off()


# ------------------------
# TITLE: NRMSE Function
# ------------------------

#function to calculate nrmse (tot)
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
nrmse_interp_p4v1_mnar <- calculate_nrsme(p4_visit1_full, p4_v1_mnar_interpolation, method = "Linear Interpolation")
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


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/Interpolation/MNAR_Interpolation_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2",
       y = "NRMSE", x = "Patient")

dev.off()

# ----------------------------
# Part 2: Kalman
# ----------------------------

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
nrmse_kalman_p4v1_mnar <- calculate_nrsme(p4_visit1_full, p4_v1_mnar_kalman, method = "Kalman Smoothing")
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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/Kalman/MNAR_Kalman_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_kalman_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_kalman_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2",
       y = "NRMSE", x = "Patient")

dev.off()

# ----------------------------
# Part 3: LWMA
# ----------------------------

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
nrmse_lwma_p4v1_mnar <- calculate_nrsme(p4_visit1_full, p4_v1_mnar_lwma, method = "LWMA")
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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/LWMA/MNAR_LWMA_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_lwma_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_lwma_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2",
       y = "NRMSE", x = "Patient")

dev.off()


# ----------------------------
# Part 4: LOESS + RF
# ----------------------------

# --------------------------
# Part 4.1: NRMSE (MCAR)
# --------------------------

#call function to calcualte nrms
#LWMA
#p1
nrmse_loess_p1v1_mnar <- calculate_nrsme(p1_visit1, p1_v1_loess, method = "LOESS + RF")
nrmse_loess_p1v2_mnar <- calculate_nrsme(p1_visit2, p1_v2_loess, method = "LOESS + RF")
#p2
nrmse_loess_p2v1_mnar <- calculate_nrsme(p2_visit1, p2_v1_loess, method = "LOESS + RF")
nrmse_loess_p2v2_mnar <- calculate_nrsme(p2_visit2, p2_v2_loess, method = "LOESS + RF")
#p3
nrmse_loess_p3v1_mnar <- calculate_nrsme(p3_visit1, p3_v1_loess, method = "LOESS + RF")
nrmse_loess_p3v2_mnar <- calculate_nrsme(p3_visit2, p3_v2_loess, method = "LOESS + RF")
#p4
nrmse_loess_p4v1_mnar <- calculate_nrsme(p4_visit1_full, p4_v1_loess, method = "LOESS + RF")
nrmse_loess_p4v2_mnar <- calculate_nrsme(p4_visit2, p4_v2_loess, method = "LOESS + RF")
#p5
nrmse_loess_p5v1_mnar <- calculate_nrsme(p5_visit1, p5_v1_loess, method = "LOESS + RF")
nrmse_loess_p5v2_mnar <- calculate_nrsme(p5_visit2, p5_v2_loess, method = "LOESS + RF")
#p6
nrmse_loess_p6v1_mnar <- calculate_nrsme(p6_visit1, p6_v1_loess, method = "LOESS + RF")
nrmse_loess_p6v2_mnar <- calculate_nrsme(p6_visit2, p6_v2_loess, method = "LOESS + RF")
#p7
nrmse_loess_p7v1_mnar <- calculate_nrsme(p7_visit1, p7_v1_loess, method = "LOESS + RFS")
nrmse_loess_p7v2_mnar <- calculate_nrsme(p7_visit2, p7_v2_loess, method = "LOESS + RF")
#p8
nrmse_loess_p8v1_mnar <- calculate_nrsme(p8_visit1, p8_v1_loess, method = "LOESS + RF")
nrmse_loess_p8v2_mnar <- calculate_nrsme(p8_visit2, p8_v2_loess, method = "LOESS + RF")
#p9
nrmse_loess_p9v1_mnar <- calculate_nrsme(p9_visit1, p9_v1_loess, method = "LOESS + RF")
nrmse_loess_p9v2_mnar <- calculate_nrsme(p9_visit2, p9_v2_loess, method = "LOESS + RF")
#p10
nrmse_loess_p10v1_mnar <- calculate_nrsme(p10_visit1, p10_v1_loess, method = "LOESS + RF")
nrmse_loess_p10v2_mnar <- calculate_nrsme(p10_visit2, p10_v2_loess, method = "LOESS + RF")

#combine visit 1 
nrmse_loess_mnar_visit1 <- bind_rows(
  nrmse_loess_p1v1_mnar %>% mutate(Patient = "P1"),
  nrmse_loess_p2v1_mnar %>% mutate(Patient = "P2"),
  nrmse_loess_p3v1_mnar %>% mutate(Patient = "P3"),
  nrmse_loess_p4v1_mnar %>% mutate(Patient = "P4"),
  nrmse_loess_p5v1_mnar %>% mutate(Patient = "P5"),
  nrmse_loess_p6v1_mnar %>% mutate(Patient = "P6"),
  nrmse_loess_p7v1_mnar %>% mutate(Patient = "P7"),
  nrmse_loess_p8v1_mnar %>% mutate(Patient = "P8"),
  nrmse_loess_p9v1_mnar %>% mutate(Patient = "P9"),
  nrmse_loess_p10v1_mnar %>% mutate(Patient = "P10")
)


#combine visit 2
nrmse_loess_mnar_visit2 <- bind_rows(
  nrmse_loess_p1v2_mnar %>% mutate(Patient = "P1"),
  nrmse_loess_p2v2_mnar %>% mutate(Patient = "P2"),
  nrmse_loess_p3v2_mnar %>% mutate(Patient = "P3"),
  nrmse_loess_p4v2_mnar %>% mutate(Patient = "P4"),
  nrmse_loess_p5v2_mnar %>% mutate(Patient = "P5"),
  nrmse_loess_p6v2_mnar %>% mutate(Patient = "P6"),
  nrmse_loess_p7v2_mnar %>% mutate(Patient = "P7"),
  nrmse_loess_p8v2_mnar %>% mutate(Patient = "P8"),
  nrmse_loess_p9v2_mnar %>% mutate(Patient = "P9"),
  nrmse_loess_p10v2_mnar %>% mutate(Patient = "P10")
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/LOESS_RF/MNAR_LOESS_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_loess_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_loess_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2",
       y = "NRMSE", x = "Patient")

dev.off()



# ----------------------------
# Part 4: LSTM
# ----------------------------

# --------------------------
# Part 4.1: NRMSE
# --------------------------

#call function to calcualte nrms
#LWMA
#p1
nrmse_lstm_p1v1_mnar <- calculate_nrsme(p1_visit1, p1_v1_lstm, method = "LSTM")
nrmse_lstm_p1v2_mnar <- calculate_nrsme(p1_visit2, p1_v2_loess, method = "LSTM")
#p2
nrmse_lstm_p2v1_mnar <- calculate_nrsme(p2_visit1, p2_v1_lstm, method = "LSTM")
nrmse_lstm_p2v2_mnar <- calculate_nrsme(p2_visit2, p2_v2_loess, method = "LSTM")
#p3
nrmse_lstm_p3v1_mnar <- calculate_nrsme(p3_visit1, p3_v1_lstm, method = "LSTM")
nrmse_lstm_p3v2_mnar <- calculate_nrsme(p3_visit2, p3_v2_loess, method = "LSTM")
#p4
nrmse_lstm_p4v1_mnar <- calculate_nrsme(p4_visit1_full, p4_v1_lstm, method = "LSTM")
nrmse_lstm_p4v2_mnar <- calculate_nrsme(p4_visit2, p4_v2_loess, method = "LSTM")
#p5
nrmse_lstm_p5v1_mnar <- calculate_nrsme(p5_visit1, p5_v1_lstm, method = "LSTM")
nrmse_lstm_p5v2_mnar <- calculate_nrsme(p5_visit2, p5_v2_loess, method = "LSTM")
#p6
nrmse_lstm_p6v1_mnar <- calculate_nrsme(p6_visit1, p6_v1_lstm, method = "LSTM")
nrmse_lstm_p6v2_mnar <- calculate_nrsme(p6_visit2, p6_v2_loess, method = "LSTM")
#p7
nrmse_lstm_p7v1_mnar <- calculate_nrsme(p7_visit1, p7_v1_lstm, method = "LSTM")
nrmse_lstm_p7v2_mnar <- calculate_nrsme(p7_visit2, p7_v2_loess, method = "LSTM")
#p8
nrmse_lstm_p8v1_mnar <- calculate_nrsme(p8_visit1, p8_v1_lstm, method = "LSTM")
nrmse_lstm_p8v2_mnar <- calculate_nrsme(p8_visit2, p8_v2_loess, method = "LSTM")
#p9
nrmse_lstm_p9v1_mnar <- calculate_nrsme(p9_visit1, p9_v1_lstm, method = "LSTM")
nrmse_lstm_p9v2_mnar <- calculate_nrsme(p9_visit2, p9_v2_loess, method = "LSTM")
#p10
nrmse_lstm_p10v1_mnar <- calculate_nrsme(p10_visit1, p10_v1_lstm, method = "LSTM")
nrmse_lstm_p10v2_mnar <- calculate_nrsme(p10_visit2, p10_v2_loess, method = "LSTM")

#combine visit 1 
nrmse_lstm_mnar_visit1 <- bind_rows(
  nrmse_lstm_p1v1_mnar %>% mutate(Patient = "P1"),
  nrmse_lstm_p2v1_mnar %>% mutate(Patient = "P2"),
  nrmse_lstm_p3v1_mnar %>% mutate(Patient = "P3"),
  nrmse_lstm_p4v1_mnar %>% mutate(Patient = "P4"),
  nrmse_lstm_p5v1_mnar %>% mutate(Patient = "P5"),
  nrmse_lstm_p6v1_mnar %>% mutate(Patient = "P6"),
  nrmse_lstm_p7v1_mnar %>% mutate(Patient = "P7"),
  nrmse_lstm_p8v1_mnar %>% mutate(Patient = "P8"),
  nrmse_lstm_p9v1_mnar %>% mutate(Patient = "P9"),
  nrmse_lstm_p10v1_mnar %>% mutate(Patient = "P10")
)


#combine visit 2
nrmse_lstm_mnar_visit2 <- bind_rows(
  nrmse_lstm_p1v2_mnar %>% mutate(Patient = "P1"),
  nrmse_lstm_p2v2_mnar %>% mutate(Patient = "P2"),
  nrmse_lstm_p3v2_mnar %>% mutate(Patient = "P3"),
  nrmse_lstm_p4v2_mnar %>% mutate(Patient = "P4"),
  nrmse_lstm_p5v2_mnar %>% mutate(Patient = "P5"),
  nrmse_lstm_p6v2_mnar %>% mutate(Patient = "P6"),
  nrmse_lstm_p7v2_mnar %>% mutate(Patient = "P7"),
  nrmse_lstm_p8v2_mnar %>% mutate(Patient = "P8"),
  nrmse_lstm_p9v2_mnar %>% mutate(Patient = "P9"),
  nrmse_lstm_p10v2_mnar %>% mutate(Patient = "P10")
)


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/LSTM/MNAR_LSTM_1MV_NRMSE.pdf", width = 14, height = 10)

#plot visit 1
ggplot(nrmse_loess_mnar_visit1, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 1",
       y = "NRMSE", x = "Patient")

#plot visit 2
ggplot(nrmse_loess_mnar_visit2, aes(x = Patient, y = NRMSE)) +
  geom_boxplot(fill = "skyblue") +
  theme_minimal() +
  labs(title = "NRMSE per Patient: Visit 2",
       y = "NRMSE", x = "Patient")

dev.off()


# ----------------------------
# Part 4: All methods compared
# ----------------------------

#visit1
nrmse_visit1_tot <- bind_rows(
  nrmse_mnar_visit1,
  nrmse_kalman_mnar_visit1,
  nrmse_lwma_mnar_visit1,
  nrmse_loess_mnar_visit1,
  nrmse_lstm_mnar_visit1
)

#visit2
nrmse_visit2_tot <- bind_rows(
  nrmse_mnar_visit2,
  nrmse_kalman_mnar_visit2,
  nrmse_lwma_mnar_visit2,
  nrmse_loess_mnar_visit2,
  nrmse_lstm_mnar_visit2
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/NRMSE_MNAR_Imputation_methods.pdf", width = 16, height = 10)

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
    title = "NRMSE of Imputed Values Only (Visit 1)",
    x = "Imputation Method",
    y = "Normalized RMSE"
  ) +
  ylim(0,2)

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
    title = "NRMSE of Imputed Values Only (Visit 2)",
    x = "Imputation Method",
    y = "Normalized RMSE"
  ) +
  ylim(0,2)

ggplot(nrmse_visit1_tot, aes(x = reorder(Metabolite, NRMSE), y = NRMSE, fill = Imputation_method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "NRMSE per Metabolite (Visit 1 - All Methods)",
    x = "Metabolite",
    y = "NRMSE",
    fill = "Imputation Method"
  ) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_fill_brewer(palette = "Set2")

ggplot(nrmse_visit2_tot, aes(x = reorder(Metabolite, NRMSE), y = NRMSE, fill = Imputation_method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "NRMSE per Metabolite (Visit 1 - All Methods)",
    x = "Metabolite",
    y = "NRMSE",
    fill = "Imputation Method"
  ) +
  theme(axis.text.y = element_text(size = 10)) +
  scale_fill_brewer(palette = "Set2")

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
# Part 1: For MNAR (1MV)
# -----------------------------------

#original AUC
#visit 1
auc_p1v1_original <- calculate_auc(p1_visit1)
auc_p2v1_original <- calculate_auc(p2_visit1)
auc_p3v1_original <- calculate_auc(p3_visit1)
auc_p4v1_original <- calculate_auc(p4_visit1_full)
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
auc_p1v1_interpolation <- calculate_auc(p1_v1_mnar_interpolation)
auc_p2v1_interpolation <- calculate_auc(p2_v1_mnar_interpolation)
auc_p3v1_interpolation <- calculate_auc(p3_v1_mnar_interpolation)
auc_p4v1_interpolation <- calculate_auc(p4_v1_mnar_interpolation)
auc_p5v1_interpolation <- calculate_auc(p5_v1_mnar_interpolation)
auc_p6v1_interpolation <- calculate_auc(p6_v1_mnar_interpolation)
auc_p7v1_interpolation <- calculate_auc(p7_v1_mnar_interpolation)
auc_p8v1_interpolation <- calculate_auc(p8_v1_mnar_interpolation)
auc_p9v1_interpolation <- calculate_auc(p9_v1_mnar_interpolation)
auc_p10v1_interpolation <- calculate_auc(p10_v1_mnar_interpolation)
#visit 2
auc_p1v2_interpolation <- calculate_auc(p1_v2_mnar_interpolation)
auc_p2v2_interpolation <- calculate_auc(p2_v2_mnar_interpolation)
auc_p3v2_interpolation <- calculate_auc(p3_v2_mnar_interpolation)
auc_p4v2_interpolation <- calculate_auc(p4_v2_mnar_interpolation)
auc_p5v2_interpolation <- calculate_auc(p5_v2_mnar_interpolation)
auc_p6v2_interpolation <- calculate_auc(p6_v2_mnar_interpolation)
auc_p7v2_interpolation <- calculate_auc(p7_v2_mnar_interpolation)
auc_p8v2_interpolation <- calculate_auc(p8_v2_mnar_interpolation)
auc_p9v2_interpolation <- calculate_auc(p9_v2_mnar_interpolation)
auc_p10v2_interpolation <- calculate_auc(p10_v2_mnar_interpolation)

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
auc_p1v1_kalman <- calculate_auc(p1_v1_mnar_kalman)
auc_p2v1_kalman <- calculate_auc(p2_v1_mnar_kalman)
auc_p3v1_kalman <- calculate_auc(p3_v1_mnar_kalman)
auc_p4v1_kalman <- calculate_auc(p4_v1_mnar_kalman)
auc_p5v1_kalman <- calculate_auc(p5_v1_mnar_kalman)
auc_p6v1_kalman <- calculate_auc(p6_v1_mnar_kalman)
auc_p7v1_kalman <- calculate_auc(p7_v1_mnar_kalman)
auc_p8v1_kalman <- calculate_auc(p8_v1_mnar_kalman)
auc_p9v1_kalman <- calculate_auc(p9_v1_mnar_kalman)
auc_p10v1_kalman <- calculate_auc(p10_v1_mnar_kalman)
#visit 2
auc_p1v2_kalman <- calculate_auc(p1_v2_mnar_kalman)
auc_p2v2_kalman <- calculate_auc(p2_v2_mnar_kalman)
auc_p3v2_kalman <- calculate_auc(p3_v2_mnar_kalman)
auc_p4v2_kalman <- calculate_auc(p4_v2_mnar_kalman)
auc_p5v2_kalman <- calculate_auc(p5_v2_mnar_kalman)
auc_p6v2_kalman <- calculate_auc(p6_v2_mnar_kalman)
auc_p7v2_kalman <- calculate_auc(p7_v2_mnar_kalman)
auc_p8v2_kalman <- calculate_auc(p8_v2_mnar_kalman)
auc_p9v2_kalman <- calculate_auc(p9_v2_mnar_kalman)
auc_p10v2_kalman <- calculate_auc(p10_v2_mnar_kalman)

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
auc_p1v1_lwma <- calculate_auc(p1_v1_mnar_lwma)
auc_p2v1_lwma <- calculate_auc(p2_v1_mnar_lwma)
auc_p3v1_lwma <- calculate_auc(p3_v1_mnar_lwma)
auc_p4v1_lwma <- calculate_auc(p4_v1_mnar_lwma)
auc_p5v1_lwma <- calculate_auc(p5_v1_mnar_lwma)
auc_p6v1_lwma <- calculate_auc(p6_v1_mnar_lwma)
auc_p7v1_lwma <- calculate_auc(p7_v1_mnar_lwma)
auc_p8v1_lwma <- calculate_auc(p8_v1_mnar_lwma)
auc_p9v1_lwma <- calculate_auc(p9_v1_mnar_lwma)
auc_p10v1_lwma <- calculate_auc(p10_v1_mnar_lwma)
#visit 2
auc_p1v2_lwma <- calculate_auc(p1_v2_mnar_lwma)
auc_p2v2_lwma <- calculate_auc(p2_v2_mnar_lwma)
auc_p3v2_lwma <- calculate_auc(p3_v2_mnar_lwma)
auc_p4v2_lwma <- calculate_auc(p4_v2_mnar_lwma)
auc_p5v2_lwma <- calculate_auc(p5_v2_mnar_lwma)
auc_p6v2_lwma <- calculate_auc(p6_v2_mnar_lwma)
auc_p7v2_lwma <- calculate_auc(p7_v2_mnar_lwma)
auc_p8v2_lwma <- calculate_auc(p8_v2_mnar_lwma)
auc_p9v2_lwma <- calculate_auc(p9_v2_mnar_lwma)
auc_p10v2_lwma <- calculate_auc(p10_v2_mnar_lwma)

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


#LOESS AUC
#visit 1
auc_p1v1_loess <- calculate_auc(p1_v1_loess)
auc_p2v1_loess <- calculate_auc(p2_v1_loess)
auc_p3v1_loess <- calculate_auc(p3_v1_loess)
auc_p4v1_loess <- calculate_auc(p4_v1_loess)
auc_p5v1_loess <- calculate_auc(p5_v1_loess)
auc_p6v1_loess <- calculate_auc(p6_v1_loess)
auc_p7v1_loess <- calculate_auc(p7_v1_loess)
auc_p8v1_loess <- calculate_auc(p8_v1_loess)
auc_p9v1_loess <- calculate_auc(p9_v1_loess)
auc_p10v1_loess <- calculate_auc(p10_v1_loess)
#visit 2
auc_p1v2_loess <- calculate_auc(p1_v2_loess)
auc_p2v2_loess <- calculate_auc(p2_v2_loess)
auc_p3v2_loess <- calculate_auc(p3_v2_loess)
auc_p4v2_loess <- calculate_auc(p4_v2_loess)
auc_p5v2_loess <- calculate_auc(p5_v2_loess)
auc_p6v2_loess <- calculate_auc(p6_v2_loess)
auc_p7v2_loess <- calculate_auc(p7_v2_loess)
auc_p8v2_loess <- calculate_auc(p8_v2_loess)
auc_p9v2_loess <- calculate_auc(p9_v2_loess)
auc_p10v2_loess <- calculate_auc(p10_v2_loess)

#combine
#visit 1
loess_visit1_auc <- bind_rows(
  auc_p1v1_loess, 
  auc_p2v1_loess,
  auc_p3v1_loess,
  auc_p4v1_loess,
  auc_p5v1_loess,
  auc_p6v1_loess,
  auc_p7v1_loess,
  auc_p8v1_loess,
  auc_p9v1_loess,
  auc_p10v1_loess
)
#visit 2
loess_visit2_auc <- bind_rows(
  auc_p1v2_loess, 
  auc_p2v2_loess,
  auc_p3v2_loess,
  auc_p4v2_loess,
  auc_p5v2_loess,
  auc_p6v2_loess,
  auc_p7v2_loess,
  auc_p8v2_loess,
  auc_p9v2_loess,
  auc_p10v2_loess
)

#LSTM AUC
#visit 1
auc_p1v1_lstm <- calculate_auc(p1_v1_lstm)
auc_p2v1_lstm <- calculate_auc(p2_v1_lstm)
auc_p3v1_lstm <- calculate_auc(p3_v1_lstm)
auc_p4v1_lstm <- calculate_auc(p4_v1_lstm)
auc_p5v1_lstm <- calculate_auc(p5_v1_lstm)
auc_p6v1_lstm <- calculate_auc(p6_v1_lstm)
auc_p7v1_lstm <- calculate_auc(p7_v1_lstm)
auc_p8v1_lstm <- calculate_auc(p8_v1_lstm)
auc_p9v1_lstm <- calculate_auc(p9_v1_lstm)
auc_p10v1_lstm <- calculate_auc(p10_v1_lstm)
#visit 2
auc_p1v2_lstm <- calculate_auc(p1_v2_lstm)
auc_p2v2_lstm <- calculate_auc(p2_v2_lstm)
auc_p3v2_lstm <- calculate_auc(p3_v2_lstm)
auc_p4v2_lstm <- calculate_auc(p4_v2_lstm)
auc_p5v2_lstm <- calculate_auc(p5_v2_lstm)
auc_p6v2_lstm <- calculate_auc(p6_v2_lstm)
auc_p7v2_lstm <- calculate_auc(p7_v2_lstm)
auc_p8v2_lstm <- calculate_auc(p8_v2_lstm)
auc_p9v2_lstm <- calculate_auc(p9_v2_lstm)
auc_p10v2_lstm <- calculate_auc(p10_v2_lstm)

#combine
#visit 1
lstm_visit1_auc <- bind_rows(
  auc_p1v1_lstm, 
  auc_p2v1_lstm,
  auc_p3v1_lstm,
  auc_p4v1_lstm,
  auc_p5v1_lstm,
  auc_p6v1_lstm,
  auc_p7v1_lstm,
  auc_p8v1_lstm,
  auc_p9v1_lstm,
  auc_p10v1_lstm
)
#visit 2
lstm_visit2_auc <- bind_rows(
  auc_p1v2_lstm, 
  auc_p2v2_lstm,
  auc_p3v2_lstm,
  auc_p4v2_lstm,
  auc_p5v2_lstm,
  auc_p6v2_lstm,
  auc_p7v2_lstm,
  auc_p8v2_lstm,
  auc_p9v2_lstm,
  auc_p10v2_lstm
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
  data.frame(Method = "LWMA",          Visit = "Visit 1", stack(auc_p10v1_lwma)),
  
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p1v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p2v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p3v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p4v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p5v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p6v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p7v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p8v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p8v1_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 1", stack(auc_p10v1_loess)),
  
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p1v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p2v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p3v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p4v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p5v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p6v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p7v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p8v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p9v1_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 1", stack(auc_p10v1_lstm))
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
  data.frame(Method = "LWMA",          Visit = "Visit 2", stack(auc_p10v2_lwma)),
  
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p1v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p2v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p3v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p4v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p5v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p6v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p7v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p8v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p8v2_loess)),
  data.frame(Method = "LOESS + RF",          Visit = "Visit 2", stack(auc_p10v2_loess)),
  
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p1v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p2v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p3v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p4v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p5v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p6v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p7v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p8v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p9v2_lstm)),
  data.frame(Method = "LSTM",          Visit = "Visit 2", stack(auc_p10v2_lstm))
) %>% rename(AUC = values, Metabolite = ind)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/AUC_Density_MNAR.pdf", width = 16, height = 10)

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
    title = "Visit 2: AUC Density per Metabolite",
    x = "AUC",
    y = "Density"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

dev.off()

# ----------------------
# Part 2: Each separately
# -----------------------

#get imputation methods
methods_to_compare <- unique(visit1_auc_df$Method)
methods_to_compare <- methods_to_compare[methods_to_compare != "Original"]

#visit 1
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/AUC_Separated_V1.pdf", width = 16, height = 10)

for (method in methods_to_compare) {
  
  #subset original
  df_sub <- visit1_auc_df %>%
    filter(Method %in% c("Original", method))
  
  p <- ggplot(df_sub, aes(x = AUC, fill = Method, color = Method)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Visit 1: AUC Density - Original vs", method),
      x = "AUC",
      y = "Density"
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9)
    )
  
  print(p)
}

dev.off()

#visit 2
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/AUC_Separated_V2.pdf", width = 16, height = 10)

for (method in methods_to_compare) {
  
  #subset the original
  df_sub <- visit2_auc_df %>%
    filter(Method %in% c("Original", method))
  
  p <- ggplot(df_sub, aes(x = AUC, fill = Method, color = Method)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    facet_wrap(~ Metabolite, scales = "free") +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Visit 2: AUC Density - Original vs", method),
      x = "AUC",
      y = "Density"
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 9)
    )
  
  print(p)
}

dev.off()


# --------------------------------------------------------
# TITLE: Pearson Correlation between original and imputed 
# --------------------------------------------------------

#function to compute Pearson correlation for each metabolite
calculate_pearson_corr_visit <- function(original_df, imputed_df) {
  #check input dimensions
  stopifnot(nrow(original_df) == nrow(imputed_df))
  
  #metadata
  metabolite_cols <- colnames(original_df)[6:ncol(original_df)]
  
  results <- data.frame(
    Metabolite = character(),
    Pearson_Correlation = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (metabolite in metabolite_cols) {
    orig <- original_df[[metabolite]]
    imputed <- imputed_df[[metabolite]]
    
    #valid indice where both original and imputed have data
    valid_idx <- which(!is.na(orig) & !is.na(imputed))
    
    #nly compute correlation if enough valid points and non-zero variance
    if (length(valid_idx) > 2 &&
        sd(orig[valid_idx]) != 0 &&
        sd(imputed[valid_idx]) != 0) {
      
      corr <- cor(orig[valid_idx], imputed[valid_idx], method = "pearson")
    } else {
      corr <- NA
    }
    
    results <- rbind(results, data.frame(
      Metabolite = metabolite,
      Pearson_Correlation = corr
    ))
  }
  
  return(results)
}

#bar plot
plot_pearson_bar <- function(corr_df, method_name = "Interpolation", visit_label = "Visit 1") {
  p <- ggplot(corr_df, aes(x = reorder(Metabolite, Pearson_Correlation), y = Pearson_Correlation)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +  # flip for better readability
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Pearson Correlation per Metabolite (", visit_label, " - ", method_name, ")", sep = ""),
      x = "Metabolite",
      y = "Pearson Correlation"
    ) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen") + #strong correlation
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange") + #moderate correction
    ylim(-0.1, 1.05)
  
  print(p)
}

# --------------------
# Part 0: Original
# --------------------

#combine original data
#v1
original_v1 <- bind_rows(
  p1_visit1, p2_visit1, p3_visit1, p4_visit1_full, p5_visit1,
  p6_visit1, p7_visit1, p8_visit1, p9_visit1, p10_visit1
)

#v2
original_v2 <- bind_rows(
  p1_visit2, p2_visit2, p3_visit2, p4_visit2, p5_visit2,
  p6_visit2, p7_visit2, p8_visit2, p9_visit2, p10_visit2
)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/Pearson_Corr_MNAR.pdf", width = 16, height = 10)


# ----------------------
# Part 1: Interpolation
# ----------------------

#v1
interpolation_v1 <- bind_rows(
  p1_v1_mnar_interpolation, p2_v1_mnar_interpolation, p3_v1_mnar_interpolation,
  p4_v1_mnar_interpolation, p5_v1_mnar_interpolation, p6_v1_mnar_interpolation,
  p7_v1_mnar_interpolation, p8_v1_mnar_interpolation, p9_v1_mnar_interpolation,
  p10_v1_mnar_interpolation
)
#v2
interpolation_v2 <- bind_rows(
  p1_v2_mnar_interpolation, p2_v2_mnar_interpolation, p3_v2_mnar_interpolation,
  p4_v2_mnar_interpolation, p5_v2_mnar_interpolation, p6_v2_mnar_interpolation,
  p7_v2_mnar_interpolation, p8_v2_mnar_interpolation, p9_v2_mnar_interpolation,
  p10_v2_mnar_interpolation
)


#call function
pearson_results_v1_interp <- calculate_pearson_corr_visit(original_v1, interpolation_v1)
pearson_results_v2_interp <- calculate_pearson_corr_visit(original_v2, interpolation_v2)

#plot results
plot_pearson_bar(pearson_results_v1_interp, method_name = "Interpolation", visit_label = "Visit 1")
plot_pearson_bar(pearson_results_v2_interp, method_name = "Interpolation", visit_label = "Visit 2")

# ----------------------
# Part 2: Kalman
# ----------------------

#v1
kalman_v1 <- bind_rows(
  p1_v1_mnar_kalman, p2_v1_mnar_kalman, p3_v1_mnar_kalman,
  p4_v1_mnar_kalman, p5_v1_mnar_kalman, p6_v1_mnar_kalman,
  p7_v1_mnar_kalman, p8_v1_mnar_kalman, p9_v1_mnar_kalman,
  p10_v1_mnar_kalman
)
#v2
kalman_v2 <- bind_rows(
  p1_v2_mnar_kalman, p2_v2_mnar_kalman, p3_v2_mnar_kalman,
  p4_v2_mnar_kalman, p5_v2_mnar_kalman, p6_v2_mnar_kalman,
  p7_v2_mnar_kalman, p8_v2_mnar_kalman, p9_v2_mnar_kalman,
  p10_v2_mnar_kalman
)


#call function
pearson_results_v1_kalman <- calculate_pearson_corr_visit(original_v1, kalman_v1)
pearson_results_v2_kalman <- calculate_pearson_corr_visit(original_v2, kalman_v2)

#plot results
plot_pearson_bar(pearson_results_v1_kalman, method_name = "Kalman", visit_label = "Visit 1")
plot_pearson_bar(pearson_results_v2_kalman, method_name = "Kalman", visit_label = "Visit 2")


# ----------------------
# Part 3: LOESS + RF
# ----------------------

loess_v1 <- bind_rows(
  p1_v1_loess, p2_v1_loess, p3_v1_loess,
  p4_v1_loess, p5_v1_loess, p6_v1_loess,
  p7_v1_loess, p8_v1_loess, p9_v1_loess,
  p10_v1_loess
)

loess_v2 <- bind_rows(
  p1_v2_loess, p2_v2_loess, p3_v2_loess,
  p4_v2_loess, p5_v2_loess, p6_v2_loess,
  p7_v2_loess, p8_v2_loess, p9_v2_loess,
  p10_v2_loess
)

pearson_results_v1_loess <- calculate_pearson_corr_visit(original_v1, loess_v1)
pearson_results_v2_loess <- calculate_pearson_corr_visit(original_v1, loess_v2)

plot_pearson_bar(pearson_results_v1_loess, method_name = "LOESS + RF", visit_label = "Visit 1")
plot_pearson_bar(pearson_results_v2_loess, method_name = "LOESS + RF", visit_label = "Visit 2")

# ----------------------
# Part 4: LWMA
# ----------------------

lwma_v1 <- bind_rows(
  p1_v1_mnar_lwma, p2_v1_mnar_lwma, p3_v1_mnar_lwma,
  p4_v1_mnar_lwma, p5_v1_mnar_lwma, p6_v1_mnar_lwma,
  p7_v1_mnar_lwma, p8_v1_mnar_lwma, p9_v1_mnar_lwma,
  p10_v1_mnar_lwma
)

lwma_v2 <- bind_rows(
  p1_v2_mnar_lwma, p2_v2_mnar_lwma, p3_v2_mnar_lwma,
  p4_v2_mnar_lwma, p5_v2_mnar_lwma, p6_v2_mnar_lwma,
  p7_v2_mnar_lwma, p8_v2_mnar_lwma, p9_v2_mnar_lwma,
  p10_v2_mnar_lwma
)

pearson_results_v1_lwma <- calculate_pearson_corr_visit(original_v1, lwma_v1)
pearson_results_v2_lwma <- calculate_pearson_corr_visit(original_v2, lwma_v2)

plot_pearson_bar(pearson_results_v1_lwma, method_name = "LWMA", visit_label = "Visit 1")
plot_pearson_bar(pearson_results_v2_lwma, method_name = "LWMA", visit_label = "Visit 2")

# ----------------------
# Part 5: LSTM
# ----------------------

lstm_v1_combined <- bind_rows(
  p1_v1_lstm, p2_v1_lstm, p3_v1_lstm,
  p4_v1_lstm, p5_v1_lstm, p6_v1_lstm,
  p7_v1_lstm, p8_v1_lstm, p9_v1_lstm,
  p10_v1_lstm
)

lstm_v2_combined <- bind_rows(
  p1_v2_lstm, p2_v2_lstm, p3_v2_lstm,
  p4_v2_lstm, p5_v2_lstm, p6_v2_lstm,
  p7_v2_lstm, p8_v2_lstm, p9_v2_lstm,
  p10_v2_lstm
)

pearson_results_v1_lstm <- calculate_pearson_corr_visit(original_v1, lstm_v1_combined)
pearson_results_v2_lstm <- calculate_pearson_corr_visit(original_v2, lstm_v2_combined)

plot_pearson_bar(pearson_results_v1_lstm, method_name = "LSTM", visit_label = "Visit 1")
plot_pearson_bar(pearson_results_v2_lstm, method_name = "LSTM", visit_label = "Visit 2")

dev.off()


# --------------------
# Part 6: All together
# ---------------------


#put all results together
#visit 1
pearson_v1_all_methods <- bind_rows(
  mutate(pearson_results_v1_interp, Method = "Interpolation"),
  mutate(pearson_results_v1_kalman, Method = "Kalman"),
  mutate(pearson_results_v1_loess, Method = "LOESS + RF"),
  mutate(pearson_results_v1_lwma, Method = "LWMA"),
  mutate(pearson_results_v1_lstm, Method = "LSTM")
)

#visit 2
pearson_v2_all_methods <- bind_rows(
  mutate(pearson_results_v2_interp, Method = "Interpolation"),
  mutate(pearson_results_v2_kalman, Method = "Kalman"),
  mutate(pearson_results_v2_loess, Method = "LOESS + RF"),
  mutate(pearson_results_v2_lwma, Method = "LWMA"),
  mutate(pearson_results_v2_lstm, Method = "LSTM")
)


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/Pearson_Corr_separated.pdf", width = 16, height = 10)

ggplot(pearson_v1_all_methods, aes(x = reorder(Metabolite, Pearson_Correlation), y = Pearson_Correlation, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Visit 1: Pearson Correlation per Metabolite and Method",
    x = "Metabolite",
    y = "Pearson Correlation",
    fill = "Method"
  ) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange") +
  ylim(-0.1, 1.05)

ggplot(pearson_v2_all_methods, aes(x = reorder(Metabolite, Pearson_Correlation), y = Pearson_Correlation, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Visit 2: Pearson Correlation per Metabolite and Method",
    x = "Metabolite",
    y = "Pearson Correlation",
    fill = "Method"
  ) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "orange") +
  ylim(-0.1, 1.05)

dev.off()

# ------------------
# TITLE: PCA
# ------------------

#Function 1: reshape for PCA into wide format
reshape_for_pca <- function(df) {
  df <- df %>%
    mutate(Visit = gsub(" ", "", Visit)) %>% #Visit 1 to Visit1
    mutate(PatientVisit = paste0(Patient, "_", Visit)) #new column e.g. P1_Visit1
  
  metadata_cols <- c("ID", "Patient", "Date", "Time_min", "Visit", "PatientVisit") #metadatae
  metabolite_cols <- setdiff(names(df)[sapply(df, is.numeric)], metadata_cols) #metabolite data
  
  df %>%
    pivot_longer(cols = all_of(metabolite_cols), names_to = "Metabolite", values_to = "Value") %>% #each metabolite in one row (new columns are PatientVisit, Time_min, Metabolite, Value)
    mutate(Feature = paste0(Metabolite, "_T", Time_min)) %>% #new metabolite name e.g. Lysin_T0
    select(PatientVisit, Feature, Value) %>% 
    pivot_wider(names_from = Feature, values_from = Value) %>% #one row per PatientVisit, one column per Metabolite_Time
    arrange(PatientVisit)
}

#call function to reshape data
#v1
original_wide_v1 <- reshape_for_pca(original_v1)
interp_wide_v1 <- reshape_for_pca(interpolation_v1)
kalman_wide_v1 <- reshape_for_pca(kalman_v1)
lwma_wide_v1 <- reshape_for_pca(lwma_v1)
loess_wide_v1 <- reshape_for_pca(loess_v1)
lstm_wide_v1 <- reshape_for_pca(lstm_v1_combined)
#v2
original_wide_v2 <- reshape_for_pca(original_v2)
interp_wide_v2 <- reshape_for_pca(interpolation_v2)
kalman_wide_v2 <- reshape_for_pca(kalman_v2)
lwma_wide_v2 <- reshape_for_pca(lwma_v2)
loess_wide_v2 <- reshape_for_pca(loess_v2)
lstm_wide_v2 <- reshape_for_pca(lstm_v2_combined)

#Function 2: run PCA
#scaling function to handle columns with zero variance
scale_safe <- function(x) {
  apply(x, 2, function(col) {
    sd_col <- sd(col)
    if (sd_col == 0) return(rep(0, length(col)))  #constant columns
    (col - mean(col)) / sd_col
  })
}

#run PCA
run_pca <- function(df_wide) {
  #save row IDs
  row_ids <- df_wide$PatientVisit
  
  #separate features from ID
  data_mat <- df_wide %>%
    column_to_rownames("PatientVisit") %>%
    as.matrix()
  
  #identify zero-variance columns (because of LSTM)
  var_vec <- apply(data_mat, 2, var)
  zero_var_cols <- names(var_vec[var_vec == 0])
  
  #drop zero-variance columns
  data_for_pca <- data_mat[, var_vec > 0, drop = FALSE]
  
  #standardize (scale and center)
  scaled_data <- scale(data_for_pca)
  rownames(scaled_data) <- row_ids  #re-apply rownames 
  
  #run PCA
  pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
  
  #attach rownames to the PCA scores matrix
  rownames(pca_result$x) <- row_ids
  
  return(pca_result)
}


#call function to run PCA
#v1
original_pca_v1 <- run_pca(original_wide_v1)
interp_pca_v1 <- run_pca(interp_wide_v1)
kalman_pca_v1 <- run_pca(kalman_wide_v1)
lwma_pca_v1 <- run_pca(lwma_wide_v1)
loess_pca_v1 <- run_pca(loess_wide_v1)
lstm_pca_v1 <- run_pca(lstm_wide_v1)
#v2
original_pca_v2 <- run_pca(original_wide_v2)
interp_pca_v2 <- run_pca(interp_wide_v2)
kalman_pca_v2 <- run_pca(kalman_wide_v2)
lwma_pca_v2 <- run_pca(lwma_wide_v2)
loess_pca_v2 <- run_pca(loess_wide_v2)
lstm_pca_v2 <- run_pca(lstm_wide_v2)

#Function 3: compute mean distances 
#use PC1 and PC2 because they usually explain the most variance in the data, give a 2D comparison of structure, same components you use for PCA plots
compute_mean_distance <- function(pca1, pca2){
  scores1 <- as.data.frame(pca1$x[,1:2])
  scores2 <- as.data.frame(pca2$x[,1:2])
  
  common_rows <- intersect(rownames(scores1), rownames(scores2))
  
  scores1 <- scores1[common_rows, , drop = FALSE]
  scores2 <- scores2[common_rows, , drop = FALSE]
  
  distances <- sqrt(rowSums((scores1 - scores2)^2))
  mean(distances)
}


#call function to compare distances
#v1
distance_interp_v1 <- compute_mean_distance(original_pca_v1, interp_pca_v1)
distance_kalman_v1 <- compute_mean_distance(original_pca_v1, kalman_pca_v1)
distance_lwma_v1   <- compute_mean_distance(original_pca_v1, lwma_pca_v1)
distance_loess_v1  <- compute_mean_distance(original_pca_v1, loess_pca_v1)
distance_lstm_v1   <- compute_mean_distance(original_pca_v1, lstm_pca_v1)

#combine into one df
results_pca_v1 <- tibble(
  Method = c("Interpolation", "Kalman", "LWMA", "LOESS+RF", "LSTM"),
  Distance = c(distance_interp_v1, distance_kalman_v1, distance_lwma_v1, distance_loess_v1, distance_lstm_v1)
)

#v2
distance_interp_v2 <- compute_mean_distance(original_pca_v2, interp_pca_v2)
distance_kalman_v2 <- compute_mean_distance(original_pca_v2, kalman_pca_v2)
distance_lwma_v2   <- compute_mean_distance(original_pca_v2, lwma_pca_v2)
distance_loess_v2  <- compute_mean_distance(original_pca_v2, loess_pca_v2)
distance_lstm_v2   <- compute_mean_distance(original_pca_v2, lstm_pca_v2)

#combine into one df
results_pca_v2 <- tibble(
  Method = c("Interpolation", "Kalman", "LWMA", "LOESS+RF", "LSTM"),
  Distance = c(distance_interp_v2, distance_kalman_v2, distance_lwma_v2, distance_loess_v2, distance_lstm_v2)
)

#Function 4: plot PCA and compare the original with the imputed 
plot_pca_comparison <- function(pca_orig, pca_imp, method_label = "Imputed", visit_label = "Visit 1"){
  #get PCA scores
  scores_orig <- as.data.frame(pca_orig$x[,1:2]) %>%
    rownames_to_column("PatientVisit") %>%
    mutate(Type = "Original")
  
  scores_imp <- as.data.frame(pca_imp$x[,1:2]) %>%
    rownames_to_column("PatientVisit") %>%
    mutate(Type = method_label)
  
  #data for arrows
  arrows <- inner_join(scores_orig, scores_imp, by = "PatientVisit", suffix = c("_orig", "_imp"))
  
  
  ggplot() +
    #arrows between original and imptued
    geom_segment(data = arrows,
                 aes(x = PC1_orig, y = PC2_orig, xend = PC1_imp, yend = PC2_imp),
                 arrow = arrow(length = unit(0.2, "cm")), color = "gray60") +
    
    #original points + labels
    geom_point(data = scores_orig, aes(x = PC1, y = PC2, color = "Original"), size = 3) +
    geom_text(data = scores_orig, aes(x = PC1, y = PC2, label = PatientVisit), vjust = -0.6, size = 3) +
    
    #imputed points + labels
    geom_point(data = scores_imp, aes(x = PC1, y = PC2, color = method_label), size = 3) +
    geom_text(data = scores_imp, aes(x = PC1, y = PC2, label = PatientVisit), vjust = -0.6, size = 3) +
    
    labs(
      title = paste("PCA Comparison -", visit_label, ":", method_label),
      x = "PC1", y = "PC2", color = "Dataset"
    ) +
    theme_minimal()
}

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/FAO/PCA_comparison.pdf", width = 16, height = 10)

#call function to plot 
#v1
plot_pca_comparison(original_pca_v1, interp_pca_v1, method_label = "Interpolation", visit_label = "Visit 1")
plot_pca_comparison(original_pca_v1, kalman_pca_v1, method_label = "Kalman", visit_label = "Visit 1")
plot_pca_comparison(original_pca_v1, lwma_pca_v1, method_label = "LWMA", visit_label = "Visit 1")
plot_pca_comparison(original_pca_v1, loess_pca_v1, method_label = "LOESS + RF", visit_label = "Visit 1")
plot_pca_comparison(original_pca_v1, lstm_pca_v1, method_label = "LSTM", visit_label = "Visit 1")
#v2
plot_pca_comparison(original_pca_v2, interp_pca_v2, method_label = "Interpolation", visit_label = "Visit 2")
plot_pca_comparison(original_pca_v2, kalman_pca_v2, method_label = "Kalman", visit_label = "Visit 2")
plot_pca_comparison(original_pca_v2, lwma_pca_v2, method_label = "LWMA", visit_label = "Visit 2")
plot_pca_comparison(original_pca_v2, loess_pca_v2, method_label = "LOESS + RF", visit_label = "Visit 2")
plot_pca_comparison(original_pca_v2, lstm_pca_v2, method_label = "LSTM", visit_label = "Visit 2")

dev.off()

# ----------------------
# TITLE: Procrustes analysis
# ---------------------------

#Measures how similar two PCA configurations are overall

#visit 1
#interpolation
proc_interp_v1 <- procrustes(original_pca_v1$x, interp_pca_v1$x)
summary(proc_interp_v1)
#kalman
proc_kalman_v1 <- procrustes(original_pca_v1$x, kalman_pca_v1$x)
summary(proc_kalman_v1)
#LWMA
proc_lwma_v1 <- procrustes(original_pca_v1$x, lwma_pca_v1$x)
summary(proc_lwma_v1)
#LOESS + RF
proc_loess_v1 <- procrustes(original_pca_v1$x, loess_pca_v1$x)
summary(proc_loess_v1)
#LSTM
proc_lstm_v1 <- procrustes(original_pca_v1$x, lstm_pca_v1$x)
summary(proc_lstm_v1)

#visit 2
#interpolation
proc_interp_v2 <- procrustes(original_pca_v2$x, interp_pca_v2$x)
summary(proc_interp_v2)
#kalman
proc_kalman_v2 <- procrustes(original_pca_v2$x, kalman_pca_v2$x)
summary(proc_kalman_v2)
#LWMA
proc_lwma_v2 <- procrustes(original_pca_v2$x, lwma_pca_v2$x)
summary(proc_lwma_v2)
#LOESS + RF
proc_loess_v2 <- procrustes(original_pca_v2$x, loess_pca_v2$x)
summary(proc_loess_v2)
#LSTM
proc_lstm_v2 <- procrustes(original_pca_v2$x, lstm_pca_v2$x)
summary(proc_lstm_v2)









