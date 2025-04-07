impute_loess_peak_aware <- function(df, time_col = "Time_min", meta_start_col = 6, sd_threshold = 5) {
  df_imputed <- df
  metabolite_cols <- names(df)[meta_start_col:ncol(df)]
  
  for (metabolite in metabolite_cols) {
    message("Fitting LOESS for: ", metabolite)
    
    time <- df[[time_col]]
    y <- df[[metabolite]]
    na_indices <- which(is.na(y))
    if (length(na_indices) == 0) next
    
    # Adjust LOESS span based on curve variability
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
    
    for (na_index in na_indices) {
      t_missing <- time[na_index]
      predicted <- NA
      
      # 1. Try LOESS
      if (!is.null(loess_fit)) {
        predicted <- predict(loess_fit, newdata = data.frame(time = t_missing))
      }
      
      # 2. Fallback: linear interpolation
      if (is.na(predicted)) {
        approx_fit <- approx(
          x = df_non_na$time,
          y = df_non_na$y,
          xout = t_missing,
          rule = 2  # allow extrapolation
        )
        predicted <- approx_fit$y
      }
      
      # 3. Peak-aware correction
      is_peak_candidate <- FALSE
      is_trough_candidate <- FALSE
      
      if (na_index > 1 && na_index < length(y)) {
        y_left <- y[na_index - 1]
        y_right <- y[na_index + 1]
        
        if (!is.na(y_left) && !is.na(y_right)) {
          is_peak_candidate <- y_left < predicted && y_right < predicted
          is_trough_candidate <- y_left > predicted && y_right > predicted
          
          if (is_peak_candidate || is_trough_candidate) {
            message("    -> Peak/trough detected at index ", na_index)
            
            t_window <- time[(na_index - 1):(na_index + 1)]
            y_window <- y[(na_index - 1):(na_index + 1)]
            df_window <- data.frame(t = t_window, y = y_window)
            
            parabola_fit <- tryCatch({
              lm(y ~ poly(t, 2, raw = TRUE), data = df_window)
            }, error = function(e) return(NULL))
            
            if (!is.null(parabola_fit)) {
              coefs <- coef(parabola_fit)
              a <- coefs[3]
              b <- coefs[2]
              t_vertex <- -b / (2 * a)
              t_vertex <- min(max(t_vertex, min(t_window)), max(t_window))  # clamp
              
              y_vertex <- predict(parabola_fit, newdata = data.frame(t = t_vertex))
              predicted <- y_vertex
              message("      ↪️  Parabola vertex used: ", round(predicted, 2))
            }
          }
        }
      }
      
      # 4. Boost underestimated peaks
      if (is_peak_candidate && na_index > 1 && na_index < length(y)) {
        y_left <- y[na_index - 1]
        y_right <- y[na_index + 1]
        peak_estimate <- max(c(y_left, y_right), na.rm = TRUE) + abs(y_left - y_right)
        
        if (!is.na(peak_estimate) && predicted < peak_estimate) {
          predicted <- peak_estimate
          message("      🔼 Boosted underestimated peak to: ", round(predicted, 2))
        }
      }
      
      # 5. Final fallback: use neighbor mean
      if (is.na(predicted)) {
        neighbors <- y[c(na_index - 1, na_index + 1)]
        predicted <- mean(neighbors, na.rm = TRUE)
        message("      ❗ Final fallback used: neighbor mean = ", round(predicted, 2))
      }
      
      df_imputed[[metabolite]][na_index] <- predicted
    }
  }
  
  return(df_imputed)
}


p1_v1_peak_aware <- impute_loess_peak_aware(p1_v1_mcar)

plot_imputed_vs_original(p1_visit1, p1_v1_peak_aware, visit = "Visit 1", type = "MCAR")
