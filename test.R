interpolate_missing2 <- function(data) {
  data_copy <- data
  time_vector <- data_copy$Time_min  # assuming this is your time column
  
  data_copy[, 6:ncol(data_copy)] <- lapply(data_copy[, 6:ncol(data_copy)], function(col) {
    if (is.numeric(col)) {
      return(interpolate_linear_with_slope(col, time = time_vector))
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}

p1_v1_interpol2 <- interpolate_missing2(p1_v1_mnar)
plot_imputed_vs_original(p1_visit1, p1_v1_interpol2, visit = "Visit 1", type = "MNAR")

