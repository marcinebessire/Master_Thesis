#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(missForest) #RF
library(zoo) #interpolation
library(imputeTS) #for imputation methods
library(pracma) #for AUC calculation
library(naniar) #to print missing values

# -------------------------------------------------
# Title: Separated based on Visit
# -------------------------------------------------

#load data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/FAO_data.csv", check.names = FALSE) #34 metabolites

FAO_data$Time_min <- as.numeric(FAO_data$Time_min)

# ------------------------------
# Part 1: Patient 4 add row 120
# ------------------------------

#copy data
FAO_data_full <- FAO_data

#ensure the row is added only once
if (!any(FAO_data_full$Patient == "P4" & FAO_data_full$Visit == "Visit 1" & FAO_data_full$Time_min == 120)) {
  
  #get all rows for Patient 4, Visit 1
  patient4_visit1 <- FAO_data_full %>% filter(Patient == "P4", Visit == "Visit 1")
  
  if (nrow(patient4_visit1) > 0) {
    #iuse the first row as a template for metadata
    new_row <- patient4_visit1[1, ]
    
    #set Time_min to 120
    new_row$Time_min <- 120
    
    #set all measurement columns (from 6 onward) to NA
    new_row[, 6:ncol(FAO_data_full)] <- NA
    
    #add the new row to the dataset
    FAO_data_full <- bind_rows(FAO_data_full, new_row)
  }
}

#get numeric columns (e.g., metabolites)
metabolite_cols <- names(FAO_data_full)[6:ncol(FAO_data_full)]

#extract Patient 4, Visit 1 (for interpoltation)
p4_v1 <- FAO_data_full %>%
  filter(Patient == "P4", Visit == "Visit 1") %>%
  arrange(Time_min)

#interpolate for each metabolite in this subset
p4_v1_interp <- p4_v1

for (col in metabolite_cols) {
  if (is.numeric(p4_v1[[col]])) {
    interpolated <- na.approx(p4_v1[[col]], x = p4_v1$Time_min, na.rm = FALSE, rule = 2)
    
    #only update the NA at Time = 120
    missing_idx <- which(p4_v1$Time_min == 120)
    if (length(missing_idx) == 1 && is.na(p4_v1[[col]][missing_idx])) {
      p4_v1_interp[[col]][missing_idx] <- interpolated[missing_idx]
    }
  }
}

#update the original full dataset
FAO_data_full <- FAO_data_full %>%
  anti_join(p4_v1 %>% filter(Time_min == 120), by = colnames(p4_v1)) %>%  
  bind_rows(p4_v1_interp %>% filter(Time_min == 120)) %>%                
  arrange(Patient, Visit, Time_min)

#sort by numeric patient number and then Date
FAO_data_full <- FAO_data_full %>%
  mutate(Patient_num = as.numeric(gsub("P", "", Patient))) %>%
  arrange(Patient_num, Date, Time_min) %>%
  select(-Patient_num)  


# ------------------------------
# Part 3: MNAR simulation
# ------------------------------

#function to simulate MNAR
MNAR_manipulation_visit <- function(data, missing_percentage){
  data_copy <- data
  
  #get Visit indices
  visit1_indices <- which(data_copy$Visit == "Visit 1")
  visit2_indices <- which(data_copy$Visit == "Visit 2")
  
  #for each measurement column (assuming starts at col 6)
  for (col in 6:ncol(data_copy)) {
    
    #visit 1
    if (length(visit1_indices) > 0) {
      values_v1 <- data_copy[visit1_indices, col]
      num_mv_v1 <- floor(length(values_v1) * missing_percentage)
      
      #get indices of the lowest values in Visit 1
      if (num_mv_v1 > 0) {
        lowest_v1 <- order(values_v1, na.last = NA)[1:num_mv_v1]
        rows_to_na_v1 <- visit1_indices[lowest_v1]
        data_copy[rows_to_na_v1, col] <- NA
      }
    }
    
    #visit 2
    if (length(visit2_indices) > 0) {
      values_v2 <- data_copy[visit2_indices, col]
      num_mv_v2 <- floor(length(values_v2) * missing_percentage)
      
      #get indices of the lowest values in Visit 2
      if (num_mv_v2 > 0) {
        lowest_v2 <- order(values_v2, na.last = NA)[1:num_mv_v2]
        rows_to_na_v2 <- visit2_indices[lowest_v2]
        data_copy[rows_to_na_v2, col] <- NA
      }
    }
  }
  
  return(data_copy)
}


#call function with different rates of mnar
FAO_10pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.10)
FAO_20pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.20)
FAO_25pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.25)
FAO_30pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.30)
FAO_35pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.35)
FAO_40pct_mnar <- MNAR_manipulation_visit(FAO_data_full, 0.40)

#split into two dataframe visit 1 and visit 2
#original
FAO_original_v1 <- FAO_data_full %>% filter(Visit == "Visit 1") #59
FAO_original_v2 <- FAO_data %>% filter(Visit == "Visit 2") #60
#10%
FAO_v1_10pct_mnar <- FAO_10pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_10pct_mnar <- FAO_10pct_mnar %>% filter(Visit == "Visit 2") #60
#20%
FAO_v1_20pct_mnar <- FAO_20pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_20pct_mnar <- FAO_20pct_mnar %>% filter(Visit == "Visit 2") #60
#25%
FAO_v1_25pct_mnar <- FAO_25pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_25pct_mnar <- FAO_25pct_mnar %>% filter(Visit == "Visit 2") #60
#30%
FAO_v1_30pct_mnar <- FAO_30pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_30pct_mnar <- FAO_30pct_mnar %>% filter(Visit == "Visit 2") #60
#35%
FAO_v1_35pct_mnar <- FAO_35pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_35pct_mnar <- FAO_35pct_mnar %>% filter(Visit == "Visit 2") #60
#40%
FAO_v1_40pct_mnar <- FAO_40pct_mnar %>% filter(Visit == "Visit 1") #60
FAO_v2_40pct_mnar <- FAO_40pct_mnar %>% filter(Visit == "Visit 2") #60

#plot all the missing value distirbution
#10%
vis_miss(FAO_v1_10pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
vis_miss(FAO_v2_10pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
#20%
vis_miss(FAO_v1_20pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
vis_miss(FAO_v2_20pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
#25%
vis_miss(FAO_v1_25pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
vis_miss(FAO_v2_25pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
#30%
vis_miss(FAO_v1_30pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
vis_miss(FAO_v2_30pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
#35%
vis_miss(FAO_v1_35pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
vis_miss(FAO_v2_35pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
#40%
vis_miss(FAO_v1_40pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])
vis_miss(FAO_v2_40pct_mnar[, 6:ncol(FAO_v1_10pct_mnar)])


# ------------------------------
# Part 3: Count Missing Values
# ------------------------------

#function to count NAs per Patient + Visit + Metabolite
get_missing_count_df <- function(data, missing_pct) {
  data %>%
    pivot_longer(cols = 6:ncol(.), names_to = "Metabolite", values_to = "Value") %>%
    filter(is.na(Value)) %>%
    group_by(Patient, Visit, Metabolite) %>%
    summarise(MissingCount = n(), .groups = "drop") %>%
    mutate(MissingPct = missing_pct)
}

#call function
#10%
missing_10 <- get_missing_count_df(FAO_10pct_mnar, 10)
missing_10pct_v1 <- get_missing_count_df(FAO_v1_10pct_mnar, 10)
missing_10pct_v2 <- get_missing_count_df(FAO_v2_10pct_mnar, 10)
#20%
missing_20 <- get_missing_count_df(FAO_20pct_mnar, 20)
missing_20pct_v1 <- get_missing_count_df(FAO_v1_20pct_mnar, 20)
missing_20pct_v2 <- get_missing_count_df(FAO_v2_20pct_mnar, 20)
#25%
missing_25 <- get_missing_count_df(FAO_25pct_mnar, 25)
missing_25pct_v1 <- get_missing_count_df(FAO_v1_25pct_mnar, 25)
missing_25pct_v2 <- get_missing_count_df(FAO_v2_25pct_mnar, 25)
#30%
missing_30 <- get_missing_count_df(FAO_30pct_mnar, 30)
missing_30pct_v1 <- get_missing_count_df(FAO_v1_30pct_mnar, 30)
missing_30pct_v2 <- get_missing_count_df(FAO_v2_30pct_mnar, 30)
#35%
missing_35 <- get_missing_count_df(FAO_35pct_mnar, 35)
missing_35pct_v1 <- get_missing_count_df(FAO_v1_35pct_mnar, 35)
missing_35pct_v2 <- get_missing_count_df(FAO_v2_35pct_mnar, 35)
#40%
missing_40 <- get_missing_count_df(FAO_40pct_mnar, 40)
missing_40pct_v1 <- get_missing_count_df(FAO_v1_40pct_mnar, 40)
missing_40pct_v2 <- get_missing_count_df(FAO_v2_40pct_mnar, 40)

#combine
all_missing_counts <- bind_rows(missing_10, missing_20, missing_25, missing_30, missing_35, missing_40)

#look into wide format
wide_missing <- all_missing_counts %>%
  pivot_wider(names_from = MissingPct, values_from = MissingCount, values_fill = 0)


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
      return(na.approx(col, na.rm = FALSE, rule = 2)) #uses the first non-NA value as a flat extrapolation (same for)
    } else {
      return(col)
    }
  })
  
  return(data_copy)
}

#call interpolation function on mcar data
#Visit1
v1_10pct_mnar_interpolation <- interpolate_missing(FAO_v1_10pct_mnar)
v1_20pct_mnar_interpolation <- interpolate_missing(FAO_v1_20pct_mnar)
v1_25pct_mnar_interpolation <- interpolate_missing(FAO_v1_25pct_mnar)
v1_30pct_mnar_interpolation <- interpolate_missing(FAO_v1_30pct_mnar)
v1_35pct_mnar_interpolation <- interpolate_missing(FAO_v1_35pct_mnar)
v1_40pct_mnar_interpolation <- interpolate_missing(FAO_v1_40pct_mnar)
#Visit2
v2_10pct_mnar_interpolation <- interpolate_missing(FAO_v2_10pct_mnar)
v2_20pct_mnar_interpolation <- interpolate_missing(FAO_v2_20pct_mnar)
v2_25pct_mnar_interpolation <- interpolate_missing(FAO_v2_25pct_mnar)
v2_30pct_mnar_interpolation <- interpolate_missing(FAO_v2_30pct_mnar)
v2_35pct_mnar_interpolation <- interpolate_missing(FAO_v2_35pct_mnar)
v2_40pct_mnar_interpolation <- interpolate_missing(FAO_v2_40pct_mnar)

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

#call kalman function
#Visit 1
v1_10pct_mnar_kalman <- kalman_imputation(FAO_v1_10pct_mnar)
v1_20pct_mnar_kalman <- kalman_imputation(FAO_v1_20pct_mnar)
v1_25pct_mnar_kalman <- kalman_imputation(FAO_v1_25pct_mnar)
v1_30pct_mnar_kalman <- kalman_imputation(FAO_v1_30pct_mnar)
v1_35pct_mnar_kalman <- kalman_imputation(FAO_v1_35pct_mnar)
v1_40pct_mnar_kalman <- kalman_imputation(FAO_v1_40pct_mnar)
#Visit 2
v2_10pct_mnar_kalman <- kalman_imputation(FAO_v2_10pct_mnar)
v2_20pct_mnar_kalman <- kalman_imputation(FAO_v2_20pct_mnar)
v2_25pct_mnar_kalman <- kalman_imputation(FAO_v2_25pct_mnar)
v2_30pct_mnar_kalman <- kalman_imputation(FAO_v2_30pct_mnar)
v2_35pct_mnar_kalman <- kalman_imputation(FAO_v2_35pct_mnar)
v2_40pct_mnar_kalman <- kalman_imputation(FAO_v2_40pct_mnar)

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

#call LWMA
#Visit 1
v1_10pct_mnar_lwma <- weighted_mov_average(FAO_v1_10pct_mnar)
v1_20pct_mnar_lwma <- weighted_mov_average(FAO_v1_20pct_mnar)
v1_25pct_mnar_lwma <- weighted_mov_average(FAO_v1_25pct_mnar)
v1_30pct_mnar_lwma <- weighted_mov_average(FAO_v1_30pct_mnar)
v1_35pct_mnar_lwma <- weighted_mov_average(FAO_v1_35pct_mnar)
v1_40pct_mnar_lwma <- weighted_mov_average(FAO_v1_40pct_mnar)
#Visit 2
v2_10pct_mnar_lwma <- weighted_mov_average(FAO_v2_10pct_mnar)
v2_20pct_mnar_lwma <- weighted_mov_average(FAO_v2_20pct_mnar)
v2_25pct_mnar_lwma <- weighted_mov_average(FAO_v2_25pct_mnar)
v2_30pct_mnar_lwma <- weighted_mov_average(FAO_v2_30pct_mnar)
v2_35pct_mnar_lwma <- weighted_mov_average(FAO_v2_35pct_mnar)
v2_40pct_mnar_lwma <- weighted_mov_average(FAO_v2_40pct_mnar)


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

#call kalman function
#Visit 1
v1_10pct_mcnar_kalman <- kalman_imputation(FAO_v1_10pct_mnar)
v1_20pct_mnar_kalman <- kalman_imputation(FAO_v1_20pct_mnar)
v1_25pct_mnar_kalman <- kalman_imputation(FAO_v1_25pct_mnar)
v1_30pct_mnar_kalman <- kalman_imputation(FAO_v1_30pct_mnar)
v1_35pct_mnar_kalman <- kalman_imputation(FAO_v1_35pct_mnar)
v1_40pct_mnar_kalman <- kalman_imputation(FAO_v1_40pct_mnar)
#Visit 2
v2_10pct_mnar_kalman <- kalman_imputation(FAO_v2_10pct_mnar)
v2_20pct_mnar_kalman <- kalman_imputation(FAO_v2_20pct_mnar)
v2_25pct_mnar_kalman <- kalman_imputation(FAO_v2_25pct_mnar)
v2_30pct_mnar_kalman <- kalman_imputation(FAO_v2_30pct_mnar)
v2_35pct_mnar_kalman <- kalman_imputation(FAO_v2_35pct_mnar)
v2_40pct_mnar_kalman <- kalman_imputation(FAO_v2_40pct_mnar)


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

#call LWMA
#Visit 1
v1_10pct_mnar_lwma <- weighted_mov_average(FAO_v1_10pct_mnar)
v1_20pct_mnar_lwma <- weighted_mov_average(FAO_v1_20pct_mnar)
v1_25pct_mnar_lwma <- weighted_mov_average(FAO_v1_25pct_mnar)
v1_30pct_mnar_lwma <- weighted_mov_average(FAO_v1_30pct_mnar)
v1_35pct_mnar_lwma <- weighted_mov_average(FAO_v1_35pct_mnar)
v1_40pct_mnar_lwma <- weighted_mov_average(FAO_v1_40pct_mnar)
#Visit 2
v2_10pct_mnar_lwma <- weighted_mov_average(FAO_v2_10pct_mnar)
v2_20pct_mnar_lwma <- weighted_mov_average(FAO_v2_20pct_mnar)
v2_25pct_mnar_lwma <- weighted_mov_average(FAO_v2_25pct_mnar)
v2_30pct_mnar_lwma <- weighted_mov_average(FAO_v2_30pct_mnar)
v2_35pct_mnar_lwma <- weighted_mov_average(FAO_v2_35pct_mnar)
v2_40pct_mnar_lwma <- weighted_mov_average(FAO_v2_40pct_mnar)

# --------------------------------
# Part 4: LOESS + RF
# --------------------------------

#LOESS (locally estimated scatterplot smoothing)
#LOESS + Random Forest Imputation
impute_loess_then_rf <- function(df, time_col = "Time_min", sd_threshold = 5) {
  df_imputed <- df
  metabolite_cols <- names(df)[6:ncol(df)]
  
  for (metabolite in metabolite_cols) {
    message("LOESS for: ", metabolite)
    
    time <- df[[time_col]]
    y <- df[[metabolite]]
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
    
    #fit LOESS
    df_non_na <- data.frame(time = time[!is.na(y)], y = y[!is.na(y)])
    loess_fit <- tryCatch({
      loess(y ~ time, data = df_non_na, span = span, degree = 1, control = loess.control(surface = "direct"))
    }, error = function(e) {
      warning("  -> LOESS failed for ", metabolite, ": ", e$message)
      return(NULL)
    })
    
    #skip if model failed
    if (is.null(loess_fit)) next
    
    #predict only for values within the fitting range
    for (na_index in na_indices) {
      t_missing <- time[na_index]
      
      #only predict if within time range
      if (t_missing >= min(df_non_na$time) && t_missing <= max(df_non_na$time)) {
        predicted <- predict(loess_fit, newdata = data.frame(time = t_missing))
        if (!is.na(predicted)) {
          df_imputed[[metabolite]][na_index] <- predicted
        } else {
          message("  -> LOESS could not predict at time = ", t_missing)
        }
      } else {
        message("  -> Time = ", t_missing, " is outside LOESS fitting range.")
      }
    }
  }
  
  #random Forest 
  message("Running Random Forest refinement with missForest...")
  rf_data <- df_imputed[, metabolite_cols]
  rf_imputed <- missForest(rf_data)$ximp
  df_imputed[, metabolite_cols] <- rf_imputed
  
  return(df_imputed)
}


#call function for loess+rf imputation
#visit 1
v1_10pct_mnar_loess <- impute_loess_then_rf(FAO_v1_10pct_mnar)
v1_20pct_mnar_loess <- impute_loess_then_rf(FAO_v1_20pct_mnar)
v1_25pct_mnar_loess <- impute_loess_then_rf(FAO_v1_25pct_mnar)
v1_30pct_mnar_loess <- impute_loess_then_rf(FAO_v1_30pct_mnar)
v1_35pct_mnar_loess <- impute_loess_then_rf(FAO_v1_35pct_mnar)
v1_40pct_mnar_loess <- impute_loess_then_rf(FAO_v1_40pct_mnar)
#visit 2
v2_10pct_mnar_loess <- impute_loess_then_rf(FAO_v2_10pct_mnar)
v2_20pct_mnar_loess <- impute_loess_then_rf(FAO_v2_20pct_mnar)
v2_25pct_mnar_loess <- impute_loess_then_rf(FAO_v2_25pct_mnar)
v2_30pct_mnar_loess <- impute_loess_then_rf(FAO_v2_30pct_mnar)
v2_35pct_mnar_loess <- impute_loess_then_rf(FAO_v2_35pct_mnar)
v2_40pct_mnar_loess <- impute_loess_then_rf(FAO_v2_40pct_mnar)


# --------------------------------------
# TITLE: KINETICS PLOT BEFORE AND AFTER
# --------------------------------------
# -------------------------------
# Plot before and After 
# -------------------------------

#plot dots and line
plot_imputed_vs_original <- function(data, visit, percent, file_path){
  #ensure correct type
  data <- data %>%
    mutate(
      Visit = as.factor(Visit),
      Time_min = as.numeric(Time_min),
      Method = as.factor(Method)
    )
  
  #long format
  long_df <- data %>%
    pivot_longer(cols = 6:(ncol(data) - 2), #last two cols: Method, MissingPct
                 names_to = "Metabolite",
                 values_to = "Concentration")
  
  #unique patients
  patients <- unique(long_df$Patient)
  
  #pdf
  pdf(file_path, width = 16, height = 10)
  
  for (patient in patients) {
    plot_data <- long_df %>% filter(Patient == patient)
    
    p <- ggplot(plot_data, aes(x = Time_min, y = Concentration, color = Method)) +
      geom_point(size = 2) +
      geom_line(aes(group = Method), linewidth = 1) +
      facet_wrap(~ Metabolite, scales = "free_y", ncol = 6) +
      labs(
        title = paste(patient, ":", "Kinetics of all Metabolites (", visit, ", ", percent, "%)"),
        subtitle = "Original vs Imputed Methods",
        x = "Time [min]",
        y = "Concentration [µM]"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14)
      )
    
    print(p)
  }
  
  dev.off()
}


# ----------------------------
# Part 1: Combine Visit 1
# ----------------------------

#10pct missingness
v1_10pct_all <- bind_rows(
  FAO_original_v1 %>% mutate(MissingPct = 10, Method = "Original"),
  v1_10pct_mnar_interpolation %>% mutate(MissingPct = 10, Method = "Interpolation"),
  v1_10pct_mnar_kalman %>% mutate(MissingPct = 10, Method = "Kalman"),
  v1_10pct_mnar_lwma %>% mutate(MissingPct = 10, Method = "LWMA"), 
  v1_10pct_mnar_loess %>% mutate(MissingPct = 10, Method = "LOESS-RF")
)

#20pct missingness
v1_20pct_all <- bind_rows(
  FAO_original_v1 %>% mutate(MissingPct = 20, Method = "Original"),
  v1_20pct_mnar_interpolation %>% mutate(MissingPct = 20, Method = "Interpolation"),
  v1_20pct_mnar_kalman %>% mutate(MissingPct = 20, Method = "Kalman"),
  v1_20pct_mnar_lwma %>% mutate(MissingPct = 20, Method = "LWMA"),
  v1_20pct_mnar_loess %>% mutate(MissingPct = 20, Method = "LOESS-RF")
)

#25pct missingness
v1_25pct_all <- bind_rows(
  FAO_original_v1 %>% mutate(MissingPct = 25, Method = "Original"),
  v1_25pct_mnar_interpolation %>% mutate(MissingPct = 25, Method = "Interpolation"),
  v1_25pct_mnar_kalman %>% mutate(MissingPct = 25, Method = "Kalman"),
  v1_25pct_mnar_lwma %>% mutate(MissingPct = 25, Method = "LWMA"),
  v1_25pct_mnar_loess %>% mutate(MissingPct = 25, Method = "LOESS-RF")
)

#30pct missingness
v1_30pct_all <- bind_rows(
  FAO_original_v1 %>% mutate(MissingPct = 30, Method = "Original"),
  v1_30pct_mnar_interpolation %>% mutate(MissingPct = 30, Method = "Interpolation"),
  v1_30pct_mnar_kalman %>% mutate(MissingPct = 30, Method = "Kalman"),
  v1_30pct_mnar_lwma %>% mutate(MissingPct = 30, Method = "LWMA"),
  v1_30pct_mnar_loess %>% mutate(MissingPct = 30, Method = "LOESS-RF")
)

#35pct missingness
v1_35pct_all <- bind_rows(
  FAO_original_v1 %>% mutate(MissingPct = 35, Method = "Original"),
  v1_35pct_mnar_interpolation %>% mutate(MissingPct = 35, Method = "Interpolation"),
  v1_35pct_mnar_kalman %>% mutate(MissingPct = 35, Method = "Kalman"),
  v1_35pct_mnar_lwma %>% mutate(MissingPct = 35, Method = "LWMA"),
  v1_35pct_mnar_loess %>% mutate(MissingPct = 35, Method = "LOESS-RF")
  
)

#40pct missingness
v1_40pct_all <- bind_rows(
  FAO_original_v1 %>% mutate(MissingPct = 40, Method = "Original"),
  v1_40pct_mnar_interpolation %>% mutate(MissingPct = 40, Method = "Interpolation"),
  v1_40pct_mnar_kalman %>% mutate(MissingPct = 40, Method = "Kalman"),
  v1_40pct_mnar_lwma %>% mutate(MissingPct = 40, Method = "LWMA"),
  v1_40pct_mnar_loess %>% mutate(MissingPct = 40, Method = "LOESS-RF")
)


#named list of all datasets
all_datasets_v1 <- list(
  "10" = v1_10pct_all,
  "20" = v1_20pct_all,
  "25" = v1_25pct_all,
  "30" = v1_30pct_all,
  "35" = v1_35pct_all,
  "40" = v1_40pct_all
)

#loop through and plot each
for (pct in names(all_datasets_v1)) {
  file_path <- paste0("/Users/marcinebessire/Desktop/Master_Thesis/Visit_Separated/MNAR/Visit1_", pct, "pct_imputation_plots.pdf")
  plot_imputed_vs_original(
    data = all_datasets_v1[[pct]],
    visit = "Visit1",
    percent = pct,
    file_path = file_path
  )
}

# ----------------------------
# Part 2: Combine Visit 2
# ----------------------------

#10pct missingness
v2_10pct_all <- bind_rows(
  FAO_original_v2 %>% mutate(MissingPct = 10, Method = "Original"),
  v2_10pct_mnar_interpolation %>% mutate(MissingPct = 10, Method = "Interpolation"),
  v2_10pct_mnar_kalman %>% mutate(MissingPct = 10, Method = "Kalman"),
  v2_10pct_mnar_lwma %>% mutate(MissingPct = 10, Method = "LWMA"),
  v2_10pct_mnar_loess %>% mutate(MissingPct = 10, Method = "LOESS-RF")
)

#20pct missingness
v2_20pct_all <- bind_rows(
  FAO_original_v2 %>% mutate(MissingPct = 20, Method = "Original"),
  v2_20pct_mnar_interpolation %>% mutate(MissingPct = 20, Method = "Interpolation"),
  v2_20pct_mnar_kalman %>% mutate(MissingPct = 20, Method = "Kalman"),
  v2_20pct_mnar_lwma %>% mutate(MissingPct = 20, Method = "LWMA"),
  v2_20pct_mnar_loess %>% mutate(MissingPct = 20, Method = "LOESS-RF")
)

#25pct missingness
v2_25pct_all <- bind_rows(
  FAO_original_v2 %>% mutate(MissingPct = 25, Method = "Original"),
  v2_25pct_mnar_interpolation %>% mutate(MissingPct = 25, Method = "Interpolation"),
  v2_25pct_mnar_kalman %>% mutate(MissingPct = 25, Method = "Kalman"),
  v2_25pct_mnar_lwma %>% mutate(MissingPct = 25, Method = "LWMA"),
  v2_25pct_mnar_loess %>% mutate(MissingPct = 25, Method = "LOESS-RF")
)

#30pct missingness
v2_30pct_all <- bind_rows(
  FAO_original_v2 %>% mutate(MissingPct = 30, Method = "Original"),
  v2_30pct_mnar_interpolation %>% mutate(MissingPct = 30, Method = "Interpolation"),
  v2_30pct_mnar_kalman %>% mutate(MissingPct = 30, Method = "Kalman"),
  v2_30pct_mnar_lwma %>% mutate(MissingPct = 30, Method = "LWMA"),
  v2_30pct_mnar_loess %>% mutate(MissingPct = 30, Method = "LOESS-RF")
)

#35pct missingness
v2_35pct_all <- bind_rows(
  FAO_original_v2 %>% mutate(MissingPct = 35, Method = "Original"),
  v2_35pct_mnar_interpolation %>% mutate(MissingPct = 35, Method = "Interpolation"),
  v2_35pct_mnar_kalman %>% mutate(MissingPct = 35, Method = "Kalman"),
  v2_35pct_mnar_lwma %>% mutate(MissingPct = 35, Method = "LWMA"),
  v2_35pct_mnar_loess %>% mutate(MissingPct = 35, Method = "LOESS-RF")
)

#40pct missingness
v2_40pct_all <- bind_rows(
  FAO_original_v2 %>% mutate(MissingPct = 40, Method = "Original"),
  v2_40pct_mnar_interpolation %>% mutate(MissingPct = 40, Method = "Interpolation"),
  v2_40pct_mnar_kalman %>% mutate(MissingPct = 40, Method = "Kalman"),
  v2_40pct_mnar_lwma %>% mutate(MissingPct = 40, Method = "LWMA"),
  v2_40pct_mnar_loess %>% mutate(MissingPct = 40, Method = "LOESS-RF")
)


#named list of all datasets
all_datasets_v2 <- list(
  "10" = v2_10pct_all,
  "20" = v2_20pct_all,
  "25" = v2_25pct_all,
  "30" = v2_30pct_all,
  "35" = v2_35pct_all,
  "40" = v2_40pct_all
)

#loop through and plot each
for (pct in names(all_datasets_v2)) {
  file_path <- paste0("/Users/marcinebessire/Desktop/Master_Thesis/Visit_Separated/MNAR/Visit2_", pct, "pct_imputation_plots.pdf")
  plot_imputed_vs_original(
    data = all_datasets_v2[[pct]],
    visit = "Visit1",
    percent = pct,
    file_path = file_path
  )
}

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

#Visit 1
nrmse_v1_10pct_interpolation <- calculate_nrsme(FAO_original_v1, v1_10pct_mnar_interpolation, method = "Interpolation")
nrmse_v1_20pct_interpolation <- calculate_nrsme(FAO_original_v1, v1_20pct_mnar_interpolation, method = "Interpolation")
nrmse_v1_25pct_interpolation <- calculate_nrsme(FAO_original_v1, v1_25pct_mnar_interpolation, method = "Interpolation")
nrmse_v1_30pct_interpolation <- calculate_nrsme(FAO_original_v1, v1_30pct_mnar_interpolation, method = "Interpolation")
nrmse_v1_35pct_interpolation <- calculate_nrsme(FAO_original_v1, v1_35pct_mnar_interpolation, method = "Interpolation")
nrmse_v1_40pct_interpolation <- calculate_nrsme(FAO_original_v1, v1_40pct_mnar_interpolation, method = "Interpolation")
#Visit 2
nrmse_v2_10pct_interpolation <- calculate_nrsme(FAO_original_v2, v2_10pct_mnar_interpolation, method = "Interpolation")
nrmse_v2_20pct_interpolation <- calculate_nrsme(FAO_original_v2, v2_20pct_mnar_interpolation, method = "Interpolation")
nrmse_v2_25pct_interpolation <- calculate_nrsme(FAO_original_v2, v2_25pct_mnar_interpolation, method = "Interpolation")
nrmse_v2_30pct_interpolation <- calculate_nrsme(FAO_original_v2, v2_30pct_mnar_interpolation, method = "Interpolation")
nrmse_v2_35pct_interpolation <- calculate_nrsme(FAO_original_v2, v2_35pct_mnar_interpolation, method = "Interpolation")
nrmse_v2_40pct_interpolation <- calculate_nrsme(FAO_original_v2, v2_40pct_mnar_interpolation, method = "Interpolation")

# ----------------------------
# Part 2: Kalmna
# ----------------------------

#Visit 1
nrmse_v1_10pct_kalman <- calculate_nrsme(FAO_original_v1, v1_10pct_mnar_kalman, method = "Kalman")
nrmse_v1_20pct_kalman <- calculate_nrsme(FAO_original_v1, v1_20pct_mnar_kalman, method = "Kalman")
nrmse_v1_25pct_kalman <- calculate_nrsme(FAO_original_v1, v1_25pct_mnar_kalman, method = "Kalman")
nrmse_v1_30pct_kalman <- calculate_nrsme(FAO_original_v1, v1_30pct_mnar_kalman, method = "Kalman")
nrmse_v1_35pct_kalman <- calculate_nrsme(FAO_original_v1, v1_35pct_mnar_kalman, method = "Kalman")
nrmse_v1_40pct_kalman <- calculate_nrsme(FAO_original_v1, v1_40pct_mnar_kalman, method = "Kalman")
#Visit 2
nrmse_v2_10pct_kalman <- calculate_nrsme(FAO_original_v2, v2_10pct_mnar_kalman, method = "Kalman")
nrmse_v2_20pct_kalman <- calculate_nrsme(FAO_original_v2, v2_20pct_mnar_kalman, method = "Kalman")
nrmse_v2_25pct_kalman <- calculate_nrsme(FAO_original_v2, v2_25pct_mnar_kalman, method = "Kalman")
nrmse_v2_30pct_kalman <- calculate_nrsme(FAO_original_v2, v2_30pct_mnar_kalman, method = "Kalman")
nrmse_v2_35pct_kalman <- calculate_nrsme(FAO_original_v2, v2_35pct_mnar_kalman, method = "Kalman")
nrmse_v2_40pct_kalman <- calculate_nrsme(FAO_original_v2, v2_40pct_mnar_kalman, method = "Kalman")

# ----------------------------
# Part 3:LWMA
# ----------------------------

#Visit 1
nrmse_v1_10pct_lwma <- calculate_nrsme(FAO_original_v1, v1_10pct_mnar_lwma, method = "LWMA")
nrmse_v1_20pct_lwma <- calculate_nrsme(FAO_original_v1, v1_20pct_mnar_lwma, method = "LWMA")
nrmse_v1_25pct_lwma <- calculate_nrsme(FAO_original_v1, v1_25pct_mnar_lwma, method = "LWMA")
nrmse_v1_30pct_lwma <- calculate_nrsme(FAO_original_v1, v1_30pct_mnar_lwma, method = "LWMA")
nrmse_v1_35pct_lwma <- calculate_nrsme(FAO_original_v1, v1_35pct_mnar_lwma, method = "LWMA")
nrmse_v1_40pct_lwma <- calculate_nrsme(FAO_original_v1, v1_40pct_mnar_lwma, method = "LWMA")
#Visit 2
nrmse_v2_10pct_lwma <- calculate_nrsme(FAO_original_v2, v2_10pct_mnar_lwma, method = "LWMA")
nrmse_v2_20pct_lwma <- calculate_nrsme(FAO_original_v2, v2_20pct_mnar_lwma, method = "LWMA")
nrmse_v2_25pct_lwma <- calculate_nrsme(FAO_original_v2, v2_25pct_mnar_lwma, method = "LWMA")
nrmse_v2_30pct_lwma <- calculate_nrsme(FAO_original_v2, v2_30pct_mnar_lwma, method = "LWMA")
nrmse_v2_35pct_lwma <- calculate_nrsme(FAO_original_v2, v2_35pct_mnar_lwma, method = "LWMA")
nrmse_v2_40pct_lwma <- calculate_nrsme(FAO_original_v2, v2_40pct_mnar_lwma, method = "LWMA")

# ----------------------------
# Part 4: LOESS-RF
# ----------------------------

#Visit 1
nrmse_v1_10pct_loess <- calculate_nrsme(FAO_original_v1, v1_10pct_mnar_loess, method = "LOESS-RF")
nrmse_v1_20pct_loess <- calculate_nrsme(FAO_original_v1, v1_20pct_mnar_loess, method = "LOESS-RF")
nrmse_v1_25pct_loess <- calculate_nrsme(FAO_original_v1, v1_25pct_mnar_loess, method = "LOESS-RF")
nrmse_v1_30pct_loess <- calculate_nrsme(FAO_original_v1, v1_30pct_mnar_loess, method = "LOESS-RF")
nrmse_v1_35pct_loess <- calculate_nrsme(FAO_original_v1, v1_35pct_mnar_loess, method = "LOESS-RF")
nrmse_v1_40pct_loess <- calculate_nrsme(FAO_original_v1, v1_40pct_mnar_loess, method = "LOESS-RF")
#Visit 2
nrmse_v2_10pct_loess <- calculate_nrsme(FAO_original_v2, v2_10pct_mnar_loess, method = "LOESS-RF")
nrmse_v2_20pct_loess <- calculate_nrsme(FAO_original_v2, v2_20pct_mnar_loess, method = "LOESS-RF")
nrmse_v2_25pct_loess <- calculate_nrsme(FAO_original_v2, v2_25pct_mnar_loess, method = "LOESS-RF")
nrmse_v2_30pct_loess <- calculate_nrsme(FAO_original_v2, v2_30pct_mnar_loess, method = "LOESS-RF")
nrmse_v2_35pct_loess <- calculate_nrsme(FAO_original_v2, v2_35pct_mnar_loess, method = "LOESS-RF")
nrmse_v2_40pct_loess <- calculate_nrsme(FAO_original_v2, v2_40pct_mnar_loess, method = "LOESS-RF")

# --------------------------------------
# Part 4: All methods compared (Visit 1)
# --------------------------------------

#visit1
nrmse_visit1_tot <- bind_rows(
  mutate(nrmse_v1_10pct_interpolation, percentage = 10),
  mutate(nrmse_v1_10pct_kalman, percentage = 10),
  mutate(nrmse_v1_10pct_lwma, percentage = 10),
  mutate(nrmse_v1_10pct_loess, percentage = 10),
  
  mutate(nrmse_v1_20pct_interpolation, percentage = 20),
  mutate(nrmse_v1_20pct_kalman, percentage = 20),
  mutate(nrmse_v1_20pct_lwma, percentage = 20),
  mutate(nrmse_v1_20pct_loess, percentage = 20),
  
  mutate(nrmse_v1_25pct_interpolation, percentage = 25),
  mutate(nrmse_v1_25pct_kalman, percentage = 25),
  mutate(nrmse_v1_25pct_lwma, percentage = 25),
  mutate(nrmse_v1_25pct_loess, percentage = 25),
  
  mutate(nrmse_v1_30pct_interpolation, percentage = 30),
  mutate(nrmse_v1_30pct_kalman, percentage = 30),
  mutate(nrmse_v1_30pct_lwma, percentage = 30),
  mutate(nrmse_v1_30pct_loess, percentage = 30),
  
  mutate(nrmse_v1_35pct_interpolation, percentage = 35),
  mutate(nrmse_v1_35pct_kalman, percentage = 35),
  mutate(nrmse_v1_35pct_lwma, percentage = 35),
  mutate(nrmse_v1_35pct_loess, percentage = 35),
  
  mutate(nrmse_v1_40pct_interpolation, percentage = 40),
  mutate(nrmse_v1_40pct_kalman, percentage = 40),
  mutate(nrmse_v1_40pct_lwma, percentage = 40),
  mutate(nrmse_v1_40pct_loess, percentage = 40)
)


#visit2
nrmse_visit2_tot <- bind_rows(
  mutate(nrmse_v2_10pct_interpolation, percentage = 10),
  mutate(nrmse_v2_10pct_kalman, percentage = 10),
  mutate(nrmse_v2_10pct_lwma, percentage = 10),
  mutate(nrmse_v2_10pct_loess, percentage = 10),
  
  mutate(nrmse_v2_20pct_interpolation, percentage = 20),
  mutate(nrmse_v2_20pct_kalman, percentage = 20),
  mutate(nrmse_v2_20pct_lwma, percentage = 20),
  mutate(nrmse_v2_20pct_loess, percentage = 20),
  
  mutate(nrmse_v2_25pct_interpolation, percentage = 25),
  mutate(nrmse_v2_25pct_kalman, percentage = 25),
  mutate(nrmse_v2_25pct_lwma, percentage = 25),
  mutate(nrmse_v2_25pct_loess, percentage = 25),
  
  mutate(nrmse_v2_30pct_interpolation, percentage = 30),
  mutate(nrmse_v2_30pct_kalman, percentage = 30),
  mutate(nrmse_v2_30pct_lwma, percentage = 30),
  mutate(nrmse_v2_30pct_loess, percentage = 30),
  
  mutate(nrmse_v2_35pct_interpolation, percentage = 35),
  mutate(nrmse_v2_35pct_kalman, percentage = 35),
  mutate(nrmse_v2_35pct_lwma, percentage = 35),
  mutate(nrmse_v2_35pct_loess, percentage = 35),
  
  mutate(nrmse_v2_40pct_interpolation, percentage = 40),
  mutate(nrmse_v2_40pct_kalman, percentage = 40),
  mutate(nrmse_v2_40pct_lwma, percentage = 40),
  mutate(nrmse_v2_40pct_loess, percentage = 40),
)


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Visit_Separated/MNAR/NRMSE_MCAR_Imputation_methods.pdf", width = 16, height = 10)

#plot visit 1
ggplot(nrmse_visit1_tot, aes(x = factor(percentage), y = NRMSE, fill = Imputation_method)) +
  geom_boxplot() +
  labs(
    title = "NRMSE by Imputation Method and MNAR Percentage (Visit 1)",
    x = "Missing Data Percentage (%)",
    y = "NRMSE"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(
    text = element_text(size = 12),
    legend.title = element_blank()
  )

#plot visit 2
ggplot(nrmse_visit2_tot, aes(x = factor(percentage), y = NRMSE, fill = Imputation_method)) +
  geom_boxplot() +
  labs(
    title = "NRMSE by Imputation Method and MNAR Percentage (Visit 2)",
    x = "Missing Data Percentage (%)",
    y = "NRMSE"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(
    text = element_text(size = 12),
    legend.title = element_blank()
  )

dev.off()




