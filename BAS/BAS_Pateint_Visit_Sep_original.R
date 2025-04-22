#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(zoo) #for Linear Interpolation
library(imputeTS) #for imputation methods
library(pracma) #for AUC calculation
library(missForest) #for RF
library(keras) #for LSTM


# --------------------------------------
# TITLE: PATIENT AND VISIT SEPARATED
# --------------------------------------

#load data
BAS_data <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/BAS_data.csv", check.names = FALSE) #34 metabolites

#count missing values and where they are
na_long <- BAS_data %>%
  pivot_longer(cols = 6:ncol(.), names_to = "Variable", values_to = "Value") %>%
  filter(is.na(Value)) %>%
  group_by(Patient, Visit, Variable) %>%
  summarise(Missing_Count = n(), .groups = "drop")

sum(na_long$Missing_Count) #total of 32 missing values 

# ------------------------
# TITLE: Data Preparation
# ------------------------
# ---------------------------------------------------
# Part 1: Split Dataframe according to Patient and Visit
# ----------------------------------------------------

p1_visit1 <- BAS_data[1:6,]
p1_visit2 <- BAS_data[7:12,]
p2_visit1 <- BAS_data[13:18,]
p2_visit2 <- BAS_data[19:24,]
p3_visit1 <- BAS_data[25:30,]
p3_visit2 <- BAS_data[31:36,]
p4_visit1 <- BAS_data[37:41,] #misses 120 min
p4_visit2 <- BAS_data[42:47,]
p5_visit1 <- BAS_data[48:53,]
p5_visit2 <- BAS_data[54:59,]
p6_visit1 <- BAS_data[60:65,]
p6_visit2 <- BAS_data[66:71,]
p7_visit1 <- BAS_data[72:77,]
p7_visit2 <- BAS_data[78:83,]
p8_visit1 <- BAS_data[84:89,]
p8_visit2 <- BAS_data[90:95,]
p9_visit1 <- BAS_data[96:101,]
p9_visit2 <- BAS_data[102:107,]
p10_visit1 <- BAS_data[108:113,]
p10_visit2 <- BAS_data[114:119,]

# -------------------------------
# Part 2: Patient 4 add row
# -------------------------------

#function to impute row T = 120 for patient 4
add_row_120 <- function(data){
  data_copy <- data
  
  #new row using metadata from existing row (e.g., row with Time_min == 60)
  template_row <- data_copy[which.min(abs(data_copy$Time_min - 120)), ]  #closest row (e.g., 60)
  
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

# ----------------------------------
# Part 3: Interpolate T=120 in p4_v1
# ----------------------------------

#function to interpolate only row T = 120
interpolate_only_row_120 <- function(data, time_col = "Time_min", target_time = 120) {
  data_copy <- data
  row_index <- which(data_copy[[time_col]] == target_time)
  if (length(row_index) == 0) stop("No row found with ", time_col, " == ", target_time)
  
  metabolite_cols <- names(data_copy)[6:ncol(data_copy)]
  
  for (col in metabolite_cols) {
    if (is.na(data_copy[[col]][row_index])) {
      valid_idx <- which(!is.na(data_copy[[col]]))
      if (length(valid_idx) >= 2) {
        interpolated <- approx(
          x = data_copy[[time_col]][valid_idx],
          y = data_copy[[col]][valid_idx],
          xout = target_time,
          rule = 2
        )$y
        data_copy[[col]][row_index] <- interpolated
      }
    }
  }
  
  return(data_copy)
}

#call function to interpolate row
p4_visit1 <-  interpolate_only_row_120(p4_visit1_complete)

# ---------------------------
# Part 4: MNAR Simulation
# ---------------------------

# #function to introduce 1 MNAR per dataframe (one missing value)
# MNAR_manipulation_lowest <- function(data){
#   #copy dataset to avoid modifying the original
#   data_copy <- data
#   
#   #go through each column
#   for (col in colnames(data_copy[6:ncol(data_copy)])) {
#     min_row <- which.min(data_copy[[col]])
#     data_copy[min_row, col] <- NA
#   }
#   
#   return(data_copy)
# }

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


# #call function
# #p1
# p1_v1_mnar <- MNAR_manipulation_lowest(p1_visit1)
# p1_v2_mnar <- MNAR_manipulation_lowest(p1_visit2)
# #p2
# p2_v1_mnar <- MNAR_manipulation_lowest(p2_visit1)
# p2_v2_mnar <- MNAR_manipulation_lowest(p2_visit2)
# #p3
# p3_v1_mnar <- MNAR_manipulation_lowest(p3_visit1)
# p3_v2_mnar <- MNAR_manipulation_lowest(p3_visit2)
# #p4
# p4_v1_mnar <- MNAR_manipulation_lowest(p4_visit1_full)
# p4_v2_mnar <- MNAR_manipulation_lowest(p4_visit2)
# #p5
# p5_v1_mnar <- MNAR_manipulation_lowest(p5_visit1)
# p5_v2_mnar <- MNAR_manipulation_lowest(p5_visit2)
# #p6
# p6_v1_mnar <- MNAR_manipulation_lowest(p6_visit1)
# p6_v2_mnar <- MNAR_manipulation_lowest(p6_visit2)
# #p7
# p7_v1_mnar <- MNAR_manipulation_lowest(p7_visit1)
# p7_v2_mnar <- MNAR_manipulation_lowest(p7_visit2)
# #p8
# p8_v1_mnar <- MNAR_manipulation_lowest(p8_visit1)
# p8_v2_mnar <- MNAR_manipulation_lowest(p8_visit2)
# #p9
# p9_v1_mnar <- MNAR_manipulation_lowest(p9_visit1)
# p9_v2_mnar <- MNAR_manipulation_lowest(p9_visit2)
# #p10
# p10_v1_mnar <- MNAR_manipulation_lowest(p10_visit1)
# p10_v2_mnar <- MNAR_manipulation_lowest(p10_visit2)

# --------------------------------------
# TITLE: IMPUTATION METHODS
# --------------------------------------
# --------------------------------------
# Part 1: Linear interpolation
# --------------------------------------

#interpolate missing data & fallback option
interpolate_missing <- function(data) {
  data_copy <- data
  
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[, 6:ncol(data_copy)], function(col) {
    if (is.numeric(col)) {
      #try interpolation
      interp <- na.approx(col, na.rm = FALSE, rule = 2)
      
      #alternative: if any NA still present, apply LOCF and then BOCF
      if (any(is.na(interp))) {
        interp <- na_locf(interp)  #forward fill
        interp <- na_locf(interp, option = "backward")  #backward fill
      }
      return(interp)
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}



#call interpolation function on mnar data
#p1
p1_v1_mnar_interpolation <- interpolate_missing(p1_visit1)
p1_v2_mnar_interpolation <- interpolate_missing(p1_visit2)
#p2
p2_v1_mnar_interpolation <- interpolate_missing(p2_visit1)
p2_v2_mnar_interpolation <- interpolate_missing(p2_visit2)
#p3
p3_v1_mnar_interpolation <- interpolate_missing(p3_visit1)
p3_v2_mnar_interpolation <- interpolate_missing(p3_visit2)
#p4
p4_v1_mnar_interpolation <- interpolate_missing(p4_visit1) 
p4_v2_mnar_interpolation <- interpolate_missing(p4_visit2)
#p5
p5_v1_mnar_interpolation <- interpolate_missing(p5_visit1)
p5_v2_mnar_interpolation <- interpolate_missing(p5_visit2)
#p6
p6_v1_mnar_interpolation <- interpolate_missing(p6_visit1)
p6_v2_mnar_interpolation <- interpolate_missing(p6_visit2)
#p7
p7_v1_mnar_interpolation <- interpolate_missing(p7_visit1)
p7_v2_mnar_interpolation <- interpolate_missing(p7_visit2)
#p8
p8_v1_mnar_interpolation <- interpolate_missing(p8_visit1)
p8_v2_mnar_interpolation <- interpolate_missing(p8_visit2)
#p9
p9_v1_mnar_interpolation <- interpolate_missing(p9_visit1)
p9_v2_mnar_interpolation <- interpolate_missing(p9_visit2)
#p10
p10_v1_mnar_interpolation <- interpolate_missing(p10_visit1)
p10_v2_mnar_interpolation <- interpolate_missing(p10_visit2)

# --------------------------------------
# Part 2: Kalman Smoothing
# --------------------------------------

#create function to impute MV with Kalman Smoothing
kalman_imputation_fallback <- function(data) {
  data_copy <- data
  
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[, 6:ncol(data_copy)], function(col) {
    if (is.numeric(col)) {
      non_na_count <- sum(!is.na(col))
      
      if (non_na_count >= 3) {
        #try Kalman
        tryCatch({
          return(na_kalman(col, model = "StructTS", smooth = TRUE))
        }, error = function(e) {
          message("Kalman failed – using interpolation: ", e$message)
          return(na.approx(col, na.rm = FALSE, rule = 2))
        })
      } else if (non_na_count >= 2) {
        #not enough for Kalman, but enough for interpolation
        return(na.approx(col, na.rm = FALSE, rule = 2))
      } else if (non_na_count == 1) {
        #only one value then use LOCF then BOCF
        message("Only one value – using LOCF + BOCF fallback")
        col <- na_locf(col, option = "locf", na_remaining = "rev")  #last observation carried forward
        col <- na_locf(col, option = "nocb", na_remaining = "rev")  #backward fill if still NA
        return(col)
      } else {
        #all values NA
        return(col)
      }
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}


#call function for kalman smoothing
#MNAR
#p1
p1_v1_mnar_kalman <- kalman_imputation_fallback(p1_visit1)
p1_v2_mnar_kalman <- kalman_imputation_fallback(p1_visit2)
#2
p2_v1_mnar_kalman <- kalman_imputation_fallback(p2_visit1)
p2_v2_mnar_kalman <- kalman_imputation_fallback(p2_visit2)
#3
p3_v1_mnar_kalman <- kalman_imputation_fallback(p3_visit1)
p3_v2_mnar_kalman <- kalman_imputation_fallback(p3_visit2)
#4
p4_v1_mnar_kalman <- kalman_imputation_fallback(p4_visit1)
p4_v2_mnar_kalman <- kalman_imputation_fallback(p4_visit2) #linear interpolation
#5
p5_v1_mnar_kalman <- kalman_imputation_fallback(p5_visit1) #linear interpolation
p5_v2_mnar_kalman <- kalman_imputation_fallback(p5_visit2) #linear interpolation
#6
p6_v1_mnar_kalman <- kalman_imputation_fallback(p6_visit1)
p6_v2_mnar_kalman <- kalman_imputation_fallback(p6_visit2)
#7
p7_v1_mnar_kalman <- kalman_imputation_fallback(p7_visit1)
p7_v2_mnar_kalman <- kalman_imputation_fallback(p7_visit2)
#8
p8_v1_mnar_kalman <- kalman_imputation_fallback(p8_visit1)
p8_v2_mnar_kalman <- kalman_imputation_fallback(p8_visit2)
#9
p9_v1_mnar_kalman <- kalman_imputation_fallback(p9_visit1)
p9_v2_mnar_kalman <- kalman_imputation_fallback(p9_visit1)
#10
p10_v1_mnar_kalman <- kalman_imputation_fallback(p10_visit1)
p10_v2_mnar_kalman <- kalman_imputation_fallback(p10_visit1)


# -----------------------------------------------------------
# Part 3: Weighted Moving Average (WMA)
# ----------------------------------------------------------

#function to make weighted moving average imputation
weighted_mov_average_fallback <- function(data, window = 3) {
  data_copy <- data
  
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[, 6:ncol(data_copy)], function(col) {
    if (is.numeric(col)) {
      non_na_count <- sum(!is.na(col))
      
      if (non_na_count >= 2) {
        #try WMA
        tryCatch({
          return(na_ma(col, k = window, weighting = "exponential")) #exponential
        }, error = function(e) {
          message("WMA failed – use linear interpolation: ", e$message)
          return(na.approx(col, na.rm = FALSE, rule = 2))
        })
      } else if (non_na_count == 1) {
        #only one value, try LOCF then BOCF
        message("Only one value – using LOCF + BOCF fallback")
        col <- na_locf(col, option = "locf", na_remaining = "rev")
        col <- na_locf(col, option = "nocb", na_remaining = "rev")
        return(col)
      } else {
        #all NAs
        return(col)
      }
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}


#call function for WMA
#MNAR
#p1
p1_v1_mnar_wma <- weighted_mov_average_fallback(p1_visit1)
p1_v2_mnar_wma <- weighted_mov_average_fallback(p1_visit2)
#2
p2_v1_mnar_wma <- weighted_mov_average_fallback(p2_visit1)
p2_v2_mnar_wma <- weighted_mov_average_fallback(p2_visit2)
#3
p3_v1_mnar_wma <- weighted_mov_average_fallback(p3_visit1)
p3_v2_mnar_wma <- weighted_mov_average_fallback(p3_visit2)
#4
p4_v1_mnar_wma <- weighted_mov_average_fallback(p4_visit1)
p4_v2_mnar_wma <- weighted_mov_average_fallback(p4_visit2) #linear interpolation
#5
p5_v1_mnar_wma <- weighted_mov_average_fallback(p5_visit1) #linear interpolation
p5_v2_mnar_wma <- weighted_mov_average_fallback(p5_visit2) #linear interpolation
#6
p6_v1_mnar_wma <- weighted_mov_average_fallback(p6_visit1)
p6_v2_mnar_wma <- weighted_mov_average_fallback(p6_visit2)
#7
p7_v1_mnar_wma <- weighted_mov_average_fallback(p7_visit1)
p7_v2_mnar_wma <- weighted_mov_average_fallback(p7_visit2)
#8
p8_v1_mnar_wma <- weighted_mov_average_fallback(p8_visit1)
p8_v2_mnar_wma <- weighted_mov_average_fallback(p8_visit2)
#9
p9_v1_mnar_wma <- weighted_mov_average_fallback(p9_visit1)
p9_v2_mnar_wma <- weighted_mov_average_fallback(p9_visit2)
#10
p10_v1_mnar_wma <- weighted_mov_average_fallback(p10_visit1)
p10_v2_mnar_wma <- weighted_mov_average_fallback(p10_visit2)


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
#       message("  -> Skipping LOESS: not enough non-missing values.")
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
#     loess_fit <- tryCatch({
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
#       #only predict if within time range
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
#   #random Forest 
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
  
  #step 1: Apply log transform (avoid log(0))
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
  
  #step 2: RF on log-transformed data
  message("Running Random Forest refinement with missForest...")
  rf_data <- log_transformed[, metabolite_cols]
  rf_imputed <- missForest(rf_data)$ximp
  log_transformed[, metabolite_cols] <- rf_imputed
  
  #step 3: back-transform (inverse log)
  df_imputed[, metabolite_cols] <- exp(log_transformed[, metabolite_cols]) - 1
  
  return(df_imputed)
}

#call function for loess + rf imputation
#visit 1
p1_v1_loess <- impute_loess_then_rf(p1_visit1)
p2_v1_loess <- impute_loess_then_rf(p2_visit1)
p3_v1_loess <- impute_loess_then_rf(p3_visit1)
p4_v1_loess <- impute_loess_then_rf(p4_visit1)
p5_v1_loess <- impute_loess_then_rf(p5_visit1)
p6_v1_loess <- impute_loess_then_rf(p6_visit1)
p7_v1_loess <- impute_loess_then_rf(p7_visit1)
p8_v1_loess <- impute_loess_then_rf(p8_visit1)
p9_v1_loess <- impute_loess_then_rf(p9_visit1)
p10_v1_loess <- impute_loess_then_rf(p10_visit1)
#visit 2
p1_v2_loess <- impute_loess_then_rf(p1_visit2)
p2_v2_loess <- impute_loess_then_rf(p2_visit2)
p3_v2_loess <- impute_loess_then_rf(p3_visit2)
p4_v2_loess <- impute_loess_then_rf(p4_visit2)
p5_v2_loess <- impute_loess_then_rf(p5_visit2)
p6_v2_loess <- impute_loess_then_rf(p6_visit2)
p7_v2_loess <- impute_loess_then_rf(p7_visit2)
p8_v2_loess <- impute_loess_then_rf(p8_visit2)
p9_v2_loess <- impute_loess_then_rf(p9_visit2)
p10_v2_loess <- impute_loess_then_rf(p10_visit2)

# --------------------------------
# Part 5: LSTM
# --------------------------------

#step 1:create array/list of the dataframes 
#combine into one array
#visit1
all_list_v1 <- list(
  p1_visit1, p2_visit1, 
  p3_visit1, p4_visit1,
  p5_visit1, p6_visit1,
  p7_visit1, p8_visit1,
  p9_visit1, p10_visit1
)
#visit 2
all_list_v2 <- list(
  p1_visit2, p2_visit2, 
  p3_visit2, p4_visit2,
  p5_visit2, p6_visit2,
  p7_visit2, p8_visit2,
  p9_visit2, p10_visit2
)

#step 2:function for LSTM imputation
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
    
    updated_list <- mnar_df_list  #save updated version
  }
  
  return(updated_list)
}

#step 3:call function to get missing values via lstm
lstm_v1 <- impute_with_lstm_all_metabolites(all_list_v1)
lstm_v2 <- impute_with_lstm_all_metabolites(all_list_v2)

#step 4:fill those imputed values from lsit back into dataframe & create new dataframe
#lists of original MNAR dataframes for each visit
original_v1 <- list(p1_visit1, p2_visit1, p3_visit1, p4_visit1, p5_visit1,
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


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/Interpolation/MNAR_Interpolation_1MV.pdf", width = 14, height = 10)

#MNAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p1_visit2, p1_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p2_visit2, p2_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p3_visit2, p3_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p4_visit2, p4_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p5_visit2, p5_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p6_visit2, p6_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p7_visit2, p7_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p8_visit2, p8_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p9_visit2, p9_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mnar_interpolation, visit = "Visit 1", type = "Linear Interpolation")
plot_imputed_vs_original(p10_visit2, p10_v2_mnar_interpolation, visit = "Visit 2", type = "Linear Interpolation")

dev.off()


# ----------------------------
# Part 2: Kalman Smoothing
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/Kalman/MNAR_Kalman_1MV.pdf", width = 14, height = 10)

#MNAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p1_visit2, p1_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p2_visit2, p2_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p3_visit2, p3_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p4_visit2, p4_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p5_visit2, p5_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p6_visit2, p6_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p7_visit2, p7_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p8_visit2, p8_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p9_visit2, p9_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mnar_kalman, visit = "Visit 1", type = "Kalman")
plot_imputed_vs_original(p10_visit2, p10_v2_mnar_kalman, visit = "Visit 2", type = "Kalman")

dev.off()

# ----------------------------
# Part 3: wma
# ----------------------------


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/WMA/MNAR_WMA_1MV.pdf", width = 14, height = 10)

#MNAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p1_visit2, p1_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p2_visit2, p2_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p3_visit2, p3_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p4_visit2, p4_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p5_visit2, p5_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p6_visit2, p6_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p7_visit2, p7_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p8_visit2, p8_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p9_visit2, p9_v2_mnar_wma, visit = "Visit 2", type = "WMA")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_mnar_wma, visit = "Visit 1", type = "WMA")
plot_imputed_vs_original(p10_visit2, p10_v2_mnar_wma, visit = "Visit 2", type = "WMA")

dev.off()


# ----------------------------
# Part 4: LEOSS + RF
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/LOESS_RF/MNAR_LOESS_RF_1MV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p1_visit2, p1_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p2_visit2, p2_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p3_visit2, p3_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p4_visit2, p4_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p5_visit2, p5_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p6_visit2, p6_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p7_visit2, p7_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p8_visit2, p8_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p9_visit2, p9_v2_loess, visit = "Visit 2", type = "LOESS + RF")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_loess, visit = "Visit 1", type = "LOESS + RF")
plot_imputed_vs_original(p10_visit2, p10_v2_loess, visit = "Visit 2", type = "LOESS + RF")

dev.off()


# ----------------------------
# Part 5: LSTM
# ----------------------------

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/LSTM/MNAR_LSTM_1MV.pdf", width = 14, height = 10)

#call function
#MCAR
#p1
plot_imputed_vs_original(p1_visit1, p1_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p1_visit2, p1_v2_lstm, visit = "Visit 2", type = "LSTM")
#p2
plot_imputed_vs_original(p2_visit1, p2_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p2_visit2, p2_v2_lstm, visit = "Visit 2", type = "LSTM")
#p3
plot_imputed_vs_original(p3_visit1, p3_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p3_visit2, p3_v2_lstm, visit = "Visit 2", type = "LSTM")
#p4
plot_imputed_vs_original(p4_visit1, p4_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p4_visit2, p4_v2_lstm, visit = "Visit 2", type = "LSTM")
#p5
plot_imputed_vs_original(p5_visit1, p5_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p5_visit2, p5_v2_lstm, visit = "Visit 2", type = "LSTM")
#p6
plot_imputed_vs_original(p6_visit1, p6_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p6_visit2, p6_v2_lstm, visit = "Visit 2", type = "LSTM")
#p7
plot_imputed_vs_original(p7_visit1, p7_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p7_visit2, p7_v2_lstm, visit = "Visit 2", type = "LSTM")
#p8
plot_imputed_vs_original(p8_visit1, p8_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p8_visit2, p8_v2_lstm, visit = "Visit 2", type = "LSTM")
#p9
plot_imputed_vs_original(p9_visit1, p9_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p9_visit2, p9_v2_lstm, visit = "Visit 2", type = "LSTM")
#p10
plot_imputed_vs_original(p10_visit1, p10_v1_lstm, visit = "Visit 1", type = "LSTM")
plot_imputed_vs_original(p10_visit2, p10_v2_lstm, visit = "Visit 2", type = "LSTMs")

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
# Part 1: All together
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

#WMA AUC
#visit 1
auc_p1v1_wma <- calculate_auc(p1_v1_mnar_wma)
auc_p2v1_wma <- calculate_auc(p2_v1_mnar_wma)
auc_p3v1_wma <- calculate_auc(p3_v1_mnar_wma)
auc_p4v1_wma <- calculate_auc(p4_v1_mnar_wma)
auc_p5v1_wma <- calculate_auc(p5_v1_mnar_wma)
auc_p6v1_wma <- calculate_auc(p6_v1_mnar_wma)
auc_p7v1_wma <- calculate_auc(p7_v1_mnar_wma)
auc_p8v1_wma <- calculate_auc(p8_v1_mnar_wma)
auc_p9v1_wma <- calculate_auc(p9_v1_mnar_wma)
auc_p10v1_wma <- calculate_auc(p10_v1_mnar_wma)
#visit 2
auc_p1v2_wma <- calculate_auc(p1_v2_mnar_wma)
auc_p2v2_wma <- calculate_auc(p2_v2_mnar_wma)
auc_p3v2_wma <- calculate_auc(p3_v2_mnar_wma)
auc_p4v2_wma <- calculate_auc(p4_v2_mnar_wma)
auc_p5v2_wma <- calculate_auc(p5_v2_mnar_wma)
auc_p6v2_wma <- calculate_auc(p6_v2_mnar_wma)
auc_p7v2_wma <- calculate_auc(p7_v2_mnar_wma)
auc_p8v2_wma <- calculate_auc(p8_v2_mnar_wma)
auc_p9v2_wma <- calculate_auc(p9_v2_mnar_wma)
auc_p10v2_wma <- calculate_auc(p10_v2_mnar_wma)

#combine
#visit 1
wma_visit1_auc <- bind_rows(
  auc_p1v1_wma, 
  auc_p2v1_wma,
  auc_p3v1_wma,
  auc_p4v1_wma,
  auc_p5v1_wma,
  auc_p6v1_wma,
  auc_p7v1_wma,
  auc_p8v1_wma,
  auc_p9v1_wma,
  auc_p10v1_wma
)
#visit 2
wma_visit2_auc <- bind_rows(
  auc_p1v2_wma, 
  auc_p2v2_wma,
  auc_p3v2_wma,
  auc_p4v2_wma,
  auc_p5v2_wma,
  auc_p6v2_wma,
  auc_p7v2_wma,
  auc_p8v2_wma,
  auc_p8v2_wma,
  auc_p10v2_wma
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
  
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p1v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p2v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p3v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p4v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p5v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p6v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p7v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p8v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p9v1_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 1", stack(auc_p10v1_interpolation)),
  
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
  
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p1v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p2v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p3v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p4v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p5v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p6v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p7v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p8v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p9v1_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 1", stack(auc_p10v1_wma)),
  
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
  
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p1v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p2v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p3v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p4v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p5v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p6v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p7v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p8v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p9v2_interpolation)),
  data.frame(Method = "Linear Interpolation", Visit = "Visit 2", stack(auc_p10v2_interpolation)),
  
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
  
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p1v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p2v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p3v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p4v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p5v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p6v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p7v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p8v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p9v2_wma)),
  data.frame(Method = "WMA",          Visit = "Visit 2", stack(auc_p10v2_wma)),
  
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

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/AUC_Density_MNAR.pdf", width = 16, height = 10)

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
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/AUC_Separated_V1.pdf", width = 16, height = 10)

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
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Patient_Visit_Separated/MNAR/BAS_original/AUC_Separated_V2.pdf", width = 16, height = 10)

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


