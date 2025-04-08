library(keras)
library(dplyr)
library(tidyr)
library(tidyverse)

#step 1.) 

#combine into one array
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

#pick one metabolite to model fist
target_metabolite <- colnames(FAO_data[6]) #L-Isoleucine

#extract sequences for that metabolite
sequences <- lapply(all_mnar_list_v1, function(df) df[[target_metabolite]])
#make a matrix of shape 10x6
sequences <- do.call(rbind, sequences)

#step 2.) 

#prepare data for LSTM
X <- sequences
X[is.na(X)] <- 0 #set to 0 instead of NA

#reshape
X_array <- array(X, dim = c(nrow(X), ncol(X), 1))

Y_array <- X_array #Y and X are the same, we want to reconsturct full sequence

#step 3.) 

#build and train model 
model <- keras_model_sequential() %>%
  layer_masking(mask_value = 0, input_shape = c(6,1)) %>%
  layer_lstm(units = 64, return_sequences = TRUE) %>%
  layer_dense(units = 1)

model %>% compile(
  loss = "mse",
  optimizer  = "adam"
)

model %>% fit(X_array, Y_array, epochs = 200, batch_size = 4, verbose = 2)

#step 4.) 

#predict and impute the missing value
predicted <- model %>% predict(X_array)

#find the positions of missing values
na_positions <- which(is.na(sequences), arr.ind = TRUE)

#replace only those NAs
for (i in seq_len(nrow(na_positions))) {
  row <- na_positions[i, "row"]
  col <- na_positions[i, "col"]
  sequences[row, col] <- predicted[row, col, 1]
}



impute_with_lstm_all_metabolites <- function(mnar_df_list, n_epochs = 200, batch_size = 4, verbose = 0) {
  # Get metabolite columns (assumes they start from column 6)
  metabolite_cols <- colnames(mnar_df_list[[1]])[6:ncol(mnar_df_list[[1]])]
  
  # Initialize updated list
  updated_list <- mnar_df_list
  
  # Loop through each metabolite
  for (metabolite in metabolite_cols) {
    message("Running LSTM imputation for: ", metabolite)
    
    # Step 1: extract sequences into matrix
    sequences <- lapply(mnar_df_list, function(df) df[[metabolite]])
    sequences <- do.call(rbind, sequences)
    
    # Step 2: Prepare inputs (replace NA with 0)
    X <- sequences
    na_mask <- is.na(X)
    X[na_mask] <- 0
    
    # Reshape to 3D
    X_array <- array(X, dim = c(nrow(X), ncol(X), 1))
    Y_array <- X_array
    
    # Step 3: Build and train LSTM model
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
    
    # Step 4: Predict and insert imputed values
    predicted <- model %>% predict(X_array)
    pred_matrix <- predicted[, , 1]
    
    for (i in seq_along(mnar_df_list)) {
      na_idx <- which(is.na(mnar_df_list[[i]][[metabolite]]))
      if (length(na_idx) > 0) {
        mnar_df_list[[i]][[metabolite]][na_idx] <- pred_matrix[i, na_idx]
      }
    }
    
    # Save back
    updated_list <- mnar_df_list
  }
  
  return(updated_list)
}

lstm_v1 <- impute_with_lstm_all_metabolites(all_mnar_list_v1)

