# LOESS + Random Forest Imputation
impute_loess_then_rf <- function(df, time_col = "Time_min", meta_start_col = 6, sd_threshold = 5) {
  library(missForest)
  
  df_imputed <- df
  metabolite_cols <- names(df)[meta_start_col:ncol(df)]
  
  # Step 1: LOESS Imputation
  for (metabolite in metabolite_cols) {
    message("LOESS for: ", metabolite)
    
    time <- df[[time_col]]
    y <- df[[metabolite]]
    na_indices <- which(is.na(y))
    if (length(na_indices) == 0) next
    
    # Adaptive span
    variability <- sd(y, na.rm = TRUE)
    span <- ifelse(variability < sd_threshold, 0.4, 0.75)
    message("  -> Using span = ", span, " (SD = ", round(variability, 2), ")")
    
    df_non_na <- data.frame(time = time[!is.na(y)], y = y[!is.na(y)])
    
    loess_fit <- tryCatch({
      loess(y ~ time, data = df_non_na, span = span, degree = 2)
    }, error = function(e) {
      warning("LOESS failed for ", metabolite, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(loess_fit)) {
      for (na_index in na_indices) {
        t_missing <- time[na_index]
        predicted <- predict(loess_fit, newdata = data.frame(time = t_missing))
        
        if (is.na(predicted)) {
          message("    -> LOESS failed; fallback to linear")
          predicted <- approx(x = df_non_na$time, y = df_non_na$y, xout = t_missing, rule = 2)$y
        }
        
        df_imputed[[metabolite]][na_index] <- predicted
      }
    }
  }
  
  # Step 2: Random Forest refinement
  message("Running Random Forest refinement with missForest...")
  rf_data <- df_imputed[, metabolite_cols]
  
  rf_imputed <- missForest(rf_data)$ximp  # Imputed matrix
  
  df_imputed[, metabolite_cols] <- rf_imputed  # Replace values
  
  return(df_imputed)
}

p1_v1_loess_rf <- impute_loess_then_rf(p1_v1_mcar)

plot_imputed_vs_original(p1_visit1, p1_v1_loess_rf, visit = "Visit 1", type = "MCAR")
