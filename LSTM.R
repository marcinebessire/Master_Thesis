library(keras)

df_v1 <- FAO_v1_10pct_mcar
df_v2 <- FAO_v2_10pct_mcar

n_patients <- 10
n_timepoints <- 6
meta_start_col <- 6
metabolite_cols <- meta_start_col:ncol(df_v1)
metabolite_names <- names(df_v1)[metabolite_cols]

# Helper: Normalize and denormalize
normalize_matrix <- function(x) {
  mean_vals <- apply(x, 2, mean, na.rm = TRUE)
  sd_vals <- apply(x, 2, sd, na.rm = TRUE)
  norm <- scale(x, center = mean_vals, scale = sd_vals)
  list(norm = norm, mean = mean_vals, sd = sd_vals)
}

denormalize_matrix <- function(x, mean_vals, sd_vals) {
  sweep(x * sd_vals, 2, mean_vals, FUN = "+")
}

# Copy data for imputation
df_v1_imputed <- df_v1
df_v2_imputed <- df_v2

# Loop over each timepoint (1–6) separately per visit
for (tp in 1:n_timepoints) {
  message("Processing timepoint ", tp)
  
  # Extract data for this timepoint from both visits
  rows_v1 <- ((tp - 1) * n_patients + 1):(tp * n_patients)
  rows_v2 <- ((tp - 1) * n_patients + 1):(tp * n_patients)
  
  X_tp <- rbind(
    df_v1[rows_v1, metabolite_cols],
    df_v2[rows_v2, metabolite_cols]
  )
  
  # Normalize
  norm_obj <- normalize_matrix(X_tp)
  X_norm <- norm_obj$norm
  mean_vals <- norm_obj$mean
  sd_vals <- norm_obj$sd
  
  # Fill NAs with column mean (temporarily)
  X_filled <- apply(X_norm, 2, function(col) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })
  
  # Build and train autoencoder
  input_dim <- ncol(X_filled)
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = 'relu', input_shape = input_dim) %>%
    layer_dense(units = 32, activation = 'relu') %>%
    layer_dense(units = 64, activation = 'relu') %>%
    layer_dense(units = input_dim, activation = 'linear')
  
  model %>% compile(
    loss = "mean_squared_error",
    optimizer = optimizer_adam(learning_rate = 0.001)
  )
  
  model %>% fit(
    x = X_filled,
    y = X_filled,
    epochs = 100,
    batch_size = 8,
    verbose = 0
  )
  
  # Predict and denormalize
  X_pred <- model %>% predict(X_filled)
  X_imputed <- denormalize_matrix(X_pred, mean_vals, sd_vals)
  
  # Replace NAs only
  for (i in 1:nrow(X_tp)) {
    for (j in 1:ncol(X_tp)) {
      if (is.na(X_tp[i, j])) {
        X_tp[i, j] <- X_imputed[i, j]
      }
    }
  }
  
  # Put imputed values back into Visit 1 and Visit 2
  df_v1_imputed[rows_v1, metabolite_cols] <- X_tp[1:10, ]
  df_v2_imputed[rows_v2, metabolite_cols] <- X_tp[11:20, ]
}

# Now df_v1_imputed and df_v2_imputed have NAs imputed using other patients

