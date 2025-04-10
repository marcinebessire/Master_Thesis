library(keras)
library(dplyr)
library(tidyr)
library(tidyverse)

#step 1: create array/list of the dataframes 

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

#step 2: create function for LSTM imputation
impute_with_lstm_all_metabolites <- function(mnar_df_list, meta_start_col = 6, n_epochs = 200, batch_size = 4, verbose = 0) {
  #get metabolite columns
  metabolite_cols <- colnames(mnar_df_list[[1]])[meta_start_col:ncol(mnar_df_list[[1]])]
  
  #copy input list to update
  updated_list <- mnar_df_list
  
  for (metabolite in metabolite_cols) {
    message("Running LSTM imputation for: ", metabolite)
    
    #step 1: extract the sequences for the current metabolite
    sequences <- lapply(mnar_df_list, function(df) df[[metabolite]])
    sequences <- do.call(rbind, sequences)  # shape: patients x timepoints
    
    #step 2: replace NA with 0 for masking later
    X <- sequences
    X[is.na(X)] <- 0  # 0 will be masked
    
    #reshape to 3D: [samples, timesteps, features]
    X_array <- array(X, dim = c(nrow(X), ncol(X), 1))
    Y_array <- X_array  # Autoencoder style: predict the full sequence
    
    #step 3: define LSTM model
    model <- keras_model_sequential() %>%
      layer_masking(mask_value = 0, input_shape = c(ncol(X), 1)) %>%
      layer_lstm(units = 64, return_sequences = TRUE) %>%
      layer_dense(units = 1)
    
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
    
    #step 4: predict missing values
    predicted <- model %>% predict(X_array)
    pred_matrix <- predicted[, , 1]  #shape: same as original
    
    #step 5: only replace the missing values
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

#step 3: call function to get missing values via lstm
lstm_v1 <- impute_with_lstm_all_metabolites(all_mnar_list_v1)
lstm_v2 <- impute_with_lstm_all_metabolites(all_mnar_list_v2)

#step 4: fill those imputed values from lsit back into dataframe & create new dataframe
#lists of original MNAR dataframes for each visit
original_v1 <- list(p1_v1_mnar, p2_v1_mnar, p3_v1_mnar, p4_v1_mnar, p5_v1_mnar,
                    p6_v1_mnar, p7_v1_mnar, p8_v1_mnar, p9_v1_mnar, p10_v1_mnar)

original_v2 <- list(p1_v2_mnar, p2_v2_mnar, p3_v2_mnar, p4_v2_mnar, p5_v2_mnar,
                    p6_v2_mnar, p7_v2_mnar, p8_v2_mnar, p9_v2_mnar, p10_v2_mnar)

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


