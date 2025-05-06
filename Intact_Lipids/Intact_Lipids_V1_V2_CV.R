#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)
library(reshape2) #for melt
library(gtools) #for mixedsort
library(missForest)
library(imputeLCMD)
library(FSA) #for Dunns test

#open data 
data_full <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids_data.csv", check.names = FALSE)

# ----------------------
# TITLE: Data Processing
# ----------------------

# --------------------------
# Part 1: Remove Cols with NA
# --------------------------

#count how many columsn with NA
sum(colSums(is.na(data_full)) > 0) #92 

#metadata
metadata <- data_full[, 1:5]

#numeric data 
numeric_data <- data_full[, 6:ncol(data_full)] #total of 366 lipids

#remove columns with NAs
data_na_removed <- numeric_data[, colSums(is.na(numeric_data)) == 0] #remaining lipids are 274 
#92 were removed

#put back together 
data <- cbind(metadata, data_na_removed)

#separate based on visit (original)
data_original_v1 <- data %>% filter(Visit == "Visit 1")
data_original_v2 <- data %>% filter(Visit == "Visit 2")

# ---------------------------------
# Part 2: Plot Whole Distribution
# ---------------------------------

#function to plot the overall distribution of all numeric values
plot_overall_distribution <- function(original) {
  #numeric columns
  numeric_original <- original[, 6:ncol(original)]
  
  #convert the entire numeric data into a single vector
  all_values <- unlist(numeric_original, use.names = FALSE)
  
  #create a data frame for plotting
  df_all_values <- data.frame(Value = all_values)
  
  #plot the overall density distribution
  plot <- ggplot(df_all_values, aes(x = Value)) +
    geom_density(fill = "gray", alpha = 0.5, color = "black") +
    theme_minimal() +
    labs(title = "Overall Distribution of FAO Data",
         x = "Measurement",
         y = "Density") +
    xlim(-100,500)
  
  print(plot)
}

#call function
#original with NA
plot_overall_distribution(data_full)
#without NA
plot_overall_distribution(data)
#plot visit 1 data
plot_overall_distribution(data_original_v1)
#plot visit 2 data
plot_overall_distribution(data_original_v2)

# ---------------------
# TITLE: Calculate CV 
# ---------------------

# ---------------------
# Part 1: Calculate CV
# ---------------------

#function to calcualte CV
calculate_cv_dataframe <- function(data) {
  #numeric
  numeric_data <- data[, 6:ncol(data)]
  
  #calculate CV per column
  cv_values <- sapply(numeric_data, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
  })
  
  #create a CV results dataframe
  cv_df <- data.frame(
    Column = names(cv_values),
    CV = round(cv_values, 2)
  )
  
  return(cv_df)
}

#call function
cv_total <- calculate_cv_dataframe(data)

# ---------------------
# Part 2: Filter by CV
# ---------------------

#fitler by CV
filter_by_cv_threshold <- function(data, cv_threshold = 30) {
  #separate metadata and numeric parts
  metadata <- data[, 1:5]
  numeric_data <- data[, 6:ncol(data)]
  
  #calculate CV per column
  cv_values <- sapply(numeric_data, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
  })
  
  #get column names that meet CV threshold
  selected_columns <- names(cv_values)[cv_values <= cv_threshold & !is.na(cv_values)]
  
  #filter numeric data
  filtered_numeric_data <- numeric_data[, selected_columns, drop = FALSE]
  
  #return combined metadata and filtered numeric
  filtered_data <- cbind(metadata, filtered_numeric_data)
  return(filtered_data)
}

#process data of visit 1 and visit 2
#20% CV
data_cv_20 <- filter_by_cv_threshold(data, 20) #13 lipids
#25% CV
data_cv_25 <- filter_by_cv_threshold(data, 25) #29 lipids
#30% CV
data_cv_30 <- filter_by_cv_threshold(data, 30) #54 lipids
#35% CV
data_cv_35 <- filter_by_cv_threshold(data, 35) #85 lipids
#40% CV
data_cv_40 <- filter_by_cv_threshold(data, 40) #112 lipids
#50% CV
data_cv_50 <- filter_by_cv_threshold(data, 50) #166 lipids


# ---------------------------------
# TITLE: Introduce Missing Values
# ---------------------------------

# --------------------------
# Part 1: Simulate MNAR
# --------------------------

#MNAR because missing values depend on the observed data (lowest measurements are more likely missing)
#function to introduce MNAR missing values with balance across Visit 1 and Visit 2
MNAR_simulation <- function(data, missing_percentage){
  #copy dataset 
  data_copy <- data
  
  #get indices of visit 1 and visit 2
  v1_indices <- which(data_copy$Visit == "Visit 1")
  v2_indices <- which(data_copy$Visit == "Visit 2")
  
  #iterate through numeric cols
  for (col in 6:ncol(data_copy)){
    #skip non-numeric
    if (!is.numeric(data_copy[[col]])) next
    
    #number of NA to introduce based on missingness
    num_mv_tot <- round(nrow(data_copy) * missing_percentage)
    num_mv <- floor(num_mv_tot / 2) #even splitting between visit 1 and visit 2
    
    if (num_mv > 0 && length(v1_indices) > 0 && length(v2_indices) > 0){
      #visit 1
      v1_data <- data_copy[v1_indices, col]
      #find incies with lowest value
      v1_order <- order(v1_data, na.last = NA)
      v1_missing_indices <- v1_indices[v1_order[1:min(num_mv, length(v1_order))]]
      
      #visit 2
      v2_data <- data_copy[v2_indices, col]
      #find incies with lowest value
      v2_order <- order(v2_data, na.last = NA)
      v2_missing_indices <- v2_indices[v2_order[1:min(num_mv, length(v2_order))]]
      
      #set values to NA
      data_copy[v1_missing_indices, col] <- NA
      data_copy[v2_missing_indices, col] <- NA
    }
  }
  
  return(data_copy)
}
#call function to generate MNAR (evenly between visit 1 and visit 2)
#10% MV
data_cv20_10pct <- MNAR_simulation(data_cv_20, 0.1) #10% missing values
data_cv25_10pct <- MNAR_simulation(data_cv_25, 0.1) #10% missing values
data_cv30_10pct <- MNAR_simulation(data_cv_30, 0.1) #10% missing values
data_cv35_10pct <- MNAR_simulation(data_cv_35, 0.1) #10% missing values
data_cv40_10pct <- MNAR_simulation(data_cv_40, 0.1) #10% missing values
data_cv50_10pct <- MNAR_simulation(data_cv_50, 0.1) #10% missing values
#20% MV
data_cv20_20pct <- MNAR_simulation(data_cv_20, 0.2) #10% missing values
data_cv25_20pct <- MNAR_simulation(data_cv_25, 0.2) #10% missing values
data_cv30_20pct <- MNAR_simulation(data_cv_30, 0.2) #10% missing values
data_cv35_20pct <- MNAR_simulation(data_cv_35, 0.2) #10% missing values
data_cv40_20pct <- MNAR_simulation(data_cv_40, 0.2) #10% missing values
data_cv50_20pct <- MNAR_simulation(data_cv_50, 0.2) #10% missing values
#30% MV
data_cv20_30pct <- MNAR_simulation(data_cv_20, 0.3) #10% missing values
data_cv25_30pct <- MNAR_simulation(data_cv_25, 0.3) #10% missing values
data_cv30_30pct <- MNAR_simulation(data_cv_30, 0.3) #10% missing values
data_cv35_30pct <- MNAR_simulation(data_cv_35, 0.3) #10% missing values
data_cv40_30pct <- MNAR_simulation(data_cv_40, 0.3) #10% missing values
data_cv50_30pct <- MNAR_simulation(data_cv_50, 0.3) #10% missing values
#40% MV
data_cv20_40pct <- MNAR_simulation(data_cv_20, 0.4) #10% missing values
data_cv25_40pct <- MNAR_simulation(data_cv_25, 0.4) #10% missing values
data_cv30_40pct <- MNAR_simulation(data_cv_30, 0.4) #10% missing values
data_cv35_40pct <- MNAR_simulation(data_cv_35, 0.4) #10% missing values
data_cv40_40pct <- MNAR_simulation(data_cv_40, 0.4) #10% missing values
data_cv50_40pct <- MNAR_simulation(data_cv_50, 0.4) #10% missing values


# ---------------------------------
# TITLE: Prepare Data for Imputation
# ----------------------------------

# ------------------------------------
# Part 1: Separate data by Visit
# ------------------------------------

#original 
data_cv_20_v1 <- data_cv_20 %>% filter(Visit == "Visit 1")
data_cv_20_v2 <- data_cv_20 %>% filter(Visit == "Visit 2")
data_cv_25_v1 <- data_cv_25 %>% filter(Visit == "Visit 1")
data_cv_25_v2 <- data_cv_25 %>% filter(Visit == "Visit 2")
data_cv_30_v1 <- data_cv_30 %>% filter(Visit == "Visit 1")
data_cv_30_v2 <- data_cv_30 %>% filter(Visit == "Visit 2")
data_cv_35_v1 <- data_cv_35 %>% filter(Visit == "Visit 1")
data_cv_35_v2 <- data_cv_35 %>% filter(Visit == "Visit 2")
data_cv_40_v1 <- data_cv_40 %>% filter(Visit == "Visit 1")
data_cv_40_v2 <- data_cv_40 %>% filter(Visit == "Visit 2")

#10% MV
data_cv20_10pct_v1 <- data_cv20_10pct %>% filter(Visit == "Visit 1")
data_cv20_10pct_v2 <- data_cv20_10pct %>% filter(Visit == "Visit 2")
data_cv25_10pct_v1 <- data_cv25_10pct %>% filter(Visit == "Visit 1")
data_cv25_10pct_v2 <- data_cv25_10pct %>% filter(Visit == "Visit 2")
data_cv30_10pct_v1 <- data_cv30_10pct %>% filter(Visit == "Visit 1")
data_cv30_10pct_v2 <- data_cv30_10pct %>% filter(Visit == "Visit 2")
data_cv35_10pct_v1 <- data_cv35_10pct %>% filter(Visit == "Visit 1")
data_cv35_10pct_v2 <- data_cv35_10pct %>% filter(Visit == "Visit 2")
data_cv40_10pct_v1 <- data_cv40_10pct %>% filter(Visit == "Visit 1")
data_cv40_10pct_v2 <- data_cv40_10pct %>% filter(Visit == "Visit 2")
#20% MV
data_cv20_20pct_v1 <- data_cv20_20pct %>% filter(Visit == "Visit 1")
data_cv20_20pct_v2 <- data_cv20_20pct %>% filter(Visit == "Visit 2")
data_cv25_20pct_v1 <- data_cv25_20pct %>% filter(Visit == "Visit 1")
data_cv25_20pct_v2 <- data_cv25_20pct %>% filter(Visit == "Visit 2")
data_cv30_20pct_v1 <- data_cv30_20pct %>% filter(Visit == "Visit 1")
data_cv30_20pct_v2 <- data_cv30_20pct %>% filter(Visit == "Visit 2")
data_cv35_20pct_v1 <- data_cv35_20pct %>% filter(Visit == "Visit 1")
data_cv35_20pct_v2 <- data_cv35_20pct %>% filter(Visit == "Visit 2")
data_cv40_20pct_v1 <- data_cv40_20pct %>% filter(Visit == "Visit 1")
data_cv40_20pct_v2 <- data_cv40_20pct %>% filter(Visit == "Visit 2")
#30% MV
data_cv20_30pct_v1 <- data_cv20_30pct %>% filter(Visit == "Visit 1")
data_cv20_30pct_v2 <- data_cv20_30pct %>% filter(Visit == "Visit 2")
data_cv25_30pct_v1 <- data_cv25_30pct %>% filter(Visit == "Visit 1")
data_cv25_30pct_v2 <- data_cv25_30pct %>% filter(Visit == "Visit 2")
data_cv30_30pct_v1 <- data_cv30_30pct %>% filter(Visit == "Visit 1")
data_cv30_30pct_v2 <- data_cv30_30pct %>% filter(Visit == "Visit 2")
data_cv35_30pct_v1 <- data_cv35_30pct %>% filter(Visit == "Visit 1")
data_cv35_30pct_v2 <- data_cv35_30pct %>% filter(Visit == "Visit 2")
data_cv40_30pct_v1 <- data_cv40_30pct %>% filter(Visit == "Visit 1")
data_cv40_30pct_v2 <- data_cv40_30pct %>% filter(Visit == "Visit 2")
#40% MV
data_cv20_40pct_v1 <- data_cv20_40pct %>% filter(Visit == "Visit 1")
data_cv20_40pct_v2 <- data_cv20_40pct %>% filter(Visit == "Visit 2")
data_cv25_40pct_v1 <- data_cv25_40pct %>% filter(Visit == "Visit 1")
data_cv25_40pct_v2 <- data_cv25_40pct %>% filter(Visit == "Visit 2")
data_cv30_40pct_v1 <- data_cv30_40pct %>% filter(Visit == "Visit 1")
data_cv30_40pct_v2 <- data_cv30_40pct %>% filter(Visit == "Visit 2")
data_cv35_40pct_v1 <- data_cv35_40pct %>% filter(Visit == "Visit 1")
data_cv35_40pct_v2 <- data_cv35_40pct %>% filter(Visit == "Visit 2")
data_cv40_40pct_v1 <- data_cv40_40pct %>% filter(Visit == "Visit 1")
data_cv40_40pct_v2 <- data_cv40_40pct %>% filter(Visit == "Visit 2")


# --------------------------
# Part 2: Put all together
# --------------------------

#original
original_data_v1 <- list(
  cv20 = data_cv_20_v1,
  cv25 = data_cv_25_v1,
  cv30 = data_cv_30_v1,
  cv35 = data_cv_35_v1,
  cv40 = data_cv_40_v1
)

original_data_v2 <- list(
  cv20 = data_cv_20_v2,
  cv25 = data_cv_25_v2,
  cv30 = data_cv_30_v2,
  cv35 = data_cv_35_v2,
  cv40 = data_cv_40_v2
)


#mnar
mnar_data_v1 <- list(
  "cv20" = list(
    "10pct" = data_cv20_10pct_v1,
    "20pct" = data_cv20_20pct_v1,
    "30pct" = data_cv20_30pct_v1,
    "40pct" = data_cv20_40pct_v1
  ),
  "cv25" = list(
    "10pct" = data_cv25_10pct_v1,
    "20pct" = data_cv25_20pct_v1,
    "30pct" = data_cv25_30pct_v1,
    "40pct" = data_cv25_40pct_v1
  ),
  "cv30" = list(
    "10pct" = data_cv30_10pct_v1,
    "20pct" = data_cv30_20pct_v1,
    "30pct" = data_cv30_30pct_v1,
    "40pct" = data_cv30_40pct_v1
  ),
  "cv35" = list(
    "10pct" = data_cv35_10pct_v1,
    "20pct" = data_cv35_20pct_v1,
    "30pct" = data_cv35_30pct_v1,
    "40pct" = data_cv35_40pct_v1
  ),
  "cv40" = list(
    "10pct" = data_cv40_10pct_v1,
    "20pct" = data_cv40_20pct_v1,
    "30pct" = data_cv40_30pct_v1,
    "40pct" = data_cv40_40pct_v1
  )
)

mnar_data_v2 <- list(
  "cv20" = list(
    "10pct" = data_cv20_10pct_v2,
    "20pct" = data_cv20_20pct_v2,
    "30pct" = data_cv20_30pct_v2,
    "40pct" = data_cv20_40pct_v2
  ),
  "cv25" = list(
    "10pct" = data_cv25_10pct_v2,
    "20pct" = data_cv25_20pct_v2,
    "30pct" = data_cv25_30pct_v2,
    "40pct" = data_cv25_40pct_v2
  ),
  "cv30" = list(
    "10pct" = data_cv30_10pct_v2,
    "20pct" = data_cv30_20pct_v2,
    "30pct" = data_cv30_30pct_v2,
    "40pct" = data_cv30_40pct_v2
  ),
  "cv35" = list(
    "10pct" = data_cv35_10pct_v2,
    "20pct" = data_cv35_20pct_v2,
    "30pct" = data_cv35_30pct_v2,
    "40pct" = data_cv35_40pct_v2
  ),
  "cv40" = list(
    "10pct" = data_cv40_10pct_v2,
    "20pct" = data_cv40_20pct_v2,
    "30pct" = data_cv40_30pct_v2,
    "40pct" = data_cv40_40pct_v2
  )
)


# ------------------------------------
# Part 2: Visit 1 vs Visit 2 Plot 
# ------------------------------------

#reshape data into long format (use melt to convert dataframe from wide into long format)
data_long <- melt(data, id.vars = c("ID", "Patient", "Date", "Time_min", "Visit"),
                  variable.name = "Lipid", value.name = "Value")

#get unique lipids
lipids <- unique(data_long$Lipid)

#split into groups of 30 
lipid_groups <- split(lipids, ceiling(seq_along(lipids) / 26))

#open pdf device
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/Lipid_Comparison_Boxplots.pdf", width = 14, height = 10)

#generate boxplot for each lipid 
# Loop through each group and generate a plot
for (i in seq_along(lipid_groups)) {
  
  subset_data <- data_long %>% filter(Lipid %in% lipid_groups[[i]])
  
  p <- ggplot(subset_data, aes(x = Visit, y = Value, fill = Visit)) +
    geom_boxplot() +
    facet_wrap(~Lipid, scales = "free") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Comparison of Lipids Between Visits (Batch", i, ")"),
         x = "Visit",
         y = "Lipid Measurement") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold"))
  
  print(p)
}

dev.off()

# ------------------------------------
# TITLE: Imputation methods 
# ------------------------------------

# ------------------------------------
# Part 1: Half-min Imputation 
# ------------------------------------

#function for half-min imputation
half_min_imputation <- function(data){
  data_copy <- data 
  
  #metadata
  metadata <- data_copy[,1:5]
  
  #numeric
  num_data <- data_copy[, 6:ncol(data_copy)]
  
  #loop through each column/lipid
  for (col in names(num_data)){
    min_val <- min(num_data[[col]], na.rm = TRUE) #find min value
    num_data[[col]][is.na(num_data[[col]])] <- 0.5 * min_val
  }
  
  final_data <- cbind(metadata, num_data)
  
  return(final_data)
  
}

#half-min imputation to each dataset
#v1
halfmin_imputed_v1 <- lapply(mnar_data_v1, function(cv_group) {
  lapply(cv_group, half_min_imputation)
})
#v2
halfmin_imputed_v2 <- lapply(mnar_data_v2, function(cv_group) {
  lapply(cv_group, half_min_imputation)
})

#so see cv25 20 pct mnar: halfmin_imputed_v1[["cv25"]][["20pct"]] 

# --------------------------
# Part 2: KNN Imputation
# --------------------------

#function for KNN imputation
knn_imputation <- function(data){
  data_copy <- data
  
  #metadata
  metadata <- data_copy[, 1:5]
  
  #numeric 
  num_data <- data_copy[, 6:ncol(data_copy)]
  
  #first transfrom into matrix and perform knn
  imputed_data <- impute.knn(as.matrix(t(num_data)), rowmax = 0.5, colmax = 1)
  
  #transform back into dataframe
  imputed_df <- as.data.frame(t(imputed_data$data))
  
  final_data <- cbind(metadata, imputed_df)
  
  return(final_data)
  
}


#call function for KNN
#v1
knn_imputed_v1 <- lapply(mnar_data_v1, function(cv_group) {
  lapply(cv_group, knn_imputation)
})
#v2
knn_imputed_v2 <- lapply(mnar_data_v2, function(cv_group) {
  lapply(cv_group, knn_imputation)
})

# ------------------------------------
# Part 3: RF Imputation
# ------------------------------------

#function for random forest imputation
rf_imputation <- function(data){
  data_copy <- data
  
  #metadata
  metadata <- data_copy[,1:5]
  
  #numerci data
  num_data <- data_copy[,6:ncol(data_copy)]
  
  #apply missforest
  imputed_data <- missForest(num_data, maxiter = 10, ntree = 100)
  imputed_df <- as.data.frame(imputed_data$ximp) #ximp = imputed data matrix
  
  final_data <- cbind(metadata, imputed_df)
  
  return(final_data)
}


#call function for RF
#v1
rf_imputed_v1 <- lapply(mnar_data_v1, function(cv_group) {
  lapply(cv_group, rf_imputation)
})
#v2
rf_imputed_v2 <- lapply(mnar_data_v2, function(cv_group) {
  lapply(cv_group, rf_imputation)
})

# --------------------------
# Part 4: QRILC Imputation
# --------------------------

#function for QRILC imputation
qrilc_imputation <- function(data){
  #copy data
  data_copy <- data
  
  #select only numeric data
  numeric_data <- data_copy[, 6:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:5]
  
  #tranfser to log 
  log_data <- log(numeric_data + 1e-6)
  
  #impute with QRILC
  imputed_data <- impute.QRILC(log_data)[[1]] #returns list, extract imputed matrix
  
  #transform back 
  exp_imputed_data <- exp(imputed_data) - 1e-6
  
  #save as dataframe 
  imputed_df <- as.data.frame(exp_imputed_data)
  
  final_df <- cbind(meta_data, imputed_df)
  
  return(final_df)
}


#call function for QRILC
#v1
qrilc_imputed_v1 <- lapply(mnar_data_v1, function(cv_group) {
  lapply(cv_group, qrilc_imputation)
})
#v2
qrilc_imputed_v2 <- lapply(mnar_data_v2, function(cv_group) {
  lapply(cv_group, qrilc_imputation)
})


# ------------------------------------
# TITLE: Statistical Tests  
# ------------------------------------

# --------------------------
# Part 1: Shapiro-Wilk test
# --------------------------

#check normality with shapiro wilk test
shapiro_test <- function(data, cv_label = NA, pct_label = NA) {
  numeric_data <- data[, 6:ncol(data)]
  
  #shapiro-Wilk test
  shapiro_results <- apply(numeric_data, 2, function(x) shapiro.test(x)$p.value)
  adjusted_p_values <- p.adjust(shapiro_results, method = "BH")
  
  shapiro_df <- data.frame(
    Lipid = names(shapiro_results),
    p_value = shapiro_results,
    Adjusted_p_values = adjusted_p_values,
    CV = cv_label,
    Missingness = pct_label
  )
  
  significant_count <- sum(shapiro_df$Adjusted_p_values < 0.05, na.rm = TRUE)
  cat("\n---------------------------------------\n")
  cat("CV:", cv_label, "| Missingness:", pct_label, "\n")
  cat("Number of significant/non-normal Lipids: ", significant_count, "\n")
  cat("---------------------------------------\n")
  
  return(shapiro_df)
}


#call shapiro function
#original
shapiro_original <- shapiro_test(data_full) #113/366 
#removed NA
shapiro_original <- shapiro_test(data) #65/274 
#visit e
shapiro_original_v1 <- shapiro_test(data_original_v1) #7/274
shapiro_original_v2 <- shapiro_test(data_original_v2) #3/274

#Half-min
#v1
shapiro_halfmin_combined_v1 <- list()

for (cv_label in names(halfmin_imputed_v1)) {
  for (pct_label in names(halfmin_imputed_v1[[cv_label]])) {
    dataset <- halfmin_imputed_v1[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    shapiro_halfmin_combined_v1[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_halfmin_shapiro_v1 <- do.call(rbind, shapiro_halfmin_combined_v1)

#v2
shapiro_halfmin_combined_v2 <- list()

for (cv_label in names(halfmin_imputed_v2)) {
  for (pct_label in names(halfmin_imputed_v2[[cv_label]])) {
    dataset <- halfmin_imputed_v2[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    shapiro_halfmin_combined_v2[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_halfmin_shapiro_v2 <- do.call(rbind, shapiro_halfmin_combined_v2)

#KNN
#v1
shapiro_knn_combined_v1 <- list()

for (cv_label in names(knn_imputed_v1)) {
  for (pct_label in names(knn_imputed_v1[[cv_label]])) {
    dataset <- knn_imputed_v1[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    shapiro_knn_combined_v1[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_knn_shapiro_v1 <- do.call(rbind, shapiro_knn_combined_v1)

#v2
shapiro_knn_combined_v2 <- list()

for (cv_label in names(knn_imputed_v2)) {
  for (pct_label in names(knn_imputed_v2[[cv_label]])) {
    dataset <- knn_imputed_v2[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    shapiro_knn_combined_v2[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

final_knn_shapiro_v2 <- do.call(rbind, shapiro_knn_combined_v2)

#RF
#v1
shapiro_rf_combined_v1 <- list()

for (cv_label in names(rf_imputed_v1)) {
  for (pct_label in names(rf_imputed_v1[[cv_label]])) {
    dataset <- rf_imputed_v1[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    shapiro_rf_combined_v1[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_rf_shapiro_v1 <- do.call(rbind, shapiro_rf_combined_v1)

#v2
shapiro_rf_combined_v2 <- list()

for (cv_label in names(rf_imputed_v2)) {
  for (pct_label in names(rf_imputed_v2[[cv_label]])) {
    dataset <- rf_imputed_v2[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    shapiro_rf_combined_v2[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_rf_shapiro_v2 <- do.call(rbind, shapiro_rf_combined_v2)

#QRILC
#v1
qrilc_rf_combined_v1 <- list()

for (cv_label in names(qrilc_imputed_v1)) {
  for (pct_label in names(qrilc_imputed_v1[[cv_label]])) {
    dataset <- qrilc_imputed_v1[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    qrilc_rf_combined_v1[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_qrilc_shapiro_v1 <- do.call(rbind, qrilc_rf_combined_v1)

#v2
qrilc_rf_combined_v2 <- list()

for (cv_label in names(qrilc_imputed_v2)) {
  for (pct_label in names(qrilc_imputed_v2[[cv_label]])) {
    dataset <- qrilc_imputed_v2[[cv_label]][[pct_label]]
    
    result_df <- shapiro_test(dataset, cv_label, pct_label)
    qrilc_rf_combined_v2[[paste(cv_label, pct_label, sep = "_")]] <- result_df
  }
}

#combine all into one data.frame
final_qrilc_shapiro_v2 <- do.call(rbind, qrilc_rf_combined_v2)

#summarize the results
summarize_shapiro <- function(shapiro_df, method_name) {
  shapiro_df %>%
    group_by(CV, Missingness) %>%
    summarise(Non_Normal_Count = sum(Adjusted_p_values < 0.05, na.rm = TRUE)) %>%
    mutate(Method = method_name,
           Missingness = as.numeric(gsub("pct", "", Missingness))) %>%
    ungroup()
}

#summarize all the methods
#v1
summary_halfmin_v1 <- summarize_shapiro(final_halfmin_shapiro_v1, "Half-min")
summary_knn_v1     <- summarize_shapiro(final_knn_shapiro_v1, "KNN")
summary_rf_v1      <- summarize_shapiro(final_rf_shapiro_v1, "RF")
summary_qrilc_v1   <- summarize_shapiro(final_qrilc_shapiro_v1, "QRILC")

#original
summary_original_v1 <- data.frame(
  CV = NA,
  Missingness = 0,
  Non_Normal_Count = 7,
  Method = "Original"
)

#combine all
shapiro_summary1 <- bind_rows(
  summary_original_v1,
  summary_halfmin_v1,
  summary_knn_v1,
  summary_rf_v1,
  summary_qrilc_v1
)

#summarize all the methods
#v2
summary_halfmin_v2 <- summarize_shapiro(final_halfmin_shapiro_v2, "Half-min")
summary_knn_v2 <- summarize_shapiro(final_knn_shapiro_v2, "KNN")
summary_rf_v2 <- summarize_shapiro(final_rf_shapiro_v2, "RF")
summary_qrilc_v2 <- summarize_shapiro(final_qrilc_shapiro_v2, "QRILC")

#original
summary_original_v2 <- data.frame(
  CV = NA,
  Missingness = 0,
  Non_Normal_Count = 3,
  Method = "Original"
)

#combine all
shapiro_summary2 <- bind_rows(
  summary_original_v2,
  summary_halfmin_v2,
  summary_knn_v2,
  summary_rf_v2,
  summary_qrilc_v2
)

#prepare for plotting
shapiro_summary1_facet <- shapiro_summary1 %>%
  filter(!is.na(CV)) 

shapiro_summary2_facet <- shapiro_summary2 %>%
  filter(!is.na(CV)) 

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/Shapiro_Wilk.pdf", width = 14, height = 10)

#make a plot
ggplot(shapiro_summary1_facet, aes(x = Missingness, y = Non_Normal_Count, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = summary_original_v1$Non_Normal_Count, linetype = "dashed", color = "lightgreen", size = 1) +
  theme_minimal(base_size = 16) +
  facet_wrap(~CV) +
  labs(
    title = "Effect of Imputation on Normality of Lipids (Visit 1)",
    x = "Missingness Percentage (%)",
    y = "Count of Non-Normal Lipids",
    color = "Imputation Method"
  )

ggplot(shapiro_summary2_facet, aes(x = Missingness, y = Non_Normal_Count, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = summary_original_v1$Non_Normal_Count, linetype = "dashed", color = "lightgreen", size = 1) +
  theme_minimal(base_size = 16) +
  facet_wrap(~CV) +
  labs(
    title = "Effect of Imputation on Normality of Lipids (Visit 2)",
    x = "Missingness Percentage (%)",
    y = "Count of Non-Normal Lipids",
    color = "Imputation Method"
  )

dev.off()


# --------------------------
# Part 2: T-test (Visit 1 vs Visit 2)
# --------------------------

#function to make t-test between original and visit (should not be significantly different if good data)
visit_statistical_tests <- function(data, cv_label = NA, pct_label = NA){
  #numeric
  lipids <- colnames(data)[6:ncol(data)]
  
  #separate visit 1 and 2
  visit1 <- data %>% filter(Visit == "Visit 1") %>% select(Patient, all_of(lipids))
  visit2 <- data %>% filter(Visit == "Visit 2") %>% select(Patient, all_of(lipids))
  
  #merge visit 1 and visit 2 by patient
  paired_data <- merge(visit1, visit2, by = "Patient", suffixes = c("_V1", "_V2"))
  
  #convert all columns to numeric 
  data_v1 <- as.data.frame(lapply(paired_data[, paste0(lipids, "_V1")], as.numeric))
  data_v2 <- as.data.frame(lapply(paired_data[, paste0(lipids, "_V2")], as.numeric))
  
  #vector to store p-values 
  p_values_wilcox <- numeric(length(lipids))
  p_values_t_test <- numeric(length(lipids))
  
  
  #loop through each lipid and run t-test and wilcox
  for (i in seq_along(lipids)){
    lipid <- lipids[i]
    
    #wilcox
    wilcox_results <- wilcox.test(data_v1[, i], data_v2[,i], paired = TRUE, exact = FALSE)
    p_values_wilcox[i] <- wilcox_results$p.value
    #t-test
    t_test_result <- t.test(data_v1[,i], data_v2[,i], paired = TRUE)
    p_values_t_test[i] <- t_test_result$p.value
  }
  
  #adjust p-value with bonferroni-holms method
  adjusted_p_value_wilcox <- p.adjust(p_values_wilcox, method = "BH")
  adjusted_p_value_t_test <- p.adjust(p_values_t_test, method = "BH")
  
  #store results into dataframe
  results <- data.frame(
    Lipid = lipids,
    Wilcoxon_p_value = p_values_wilcox,
    Wilxocon_adjusted_p_value = adjusted_p_value_wilcox,
    T_test_p_value = p_values_t_test,
    T_Test_adjusted_p_value = adjusted_p_value_t_test
  )
  
  #significant threshold 
  alpha <- 0.05 
  
  #count significant lipids 
  wilcox_sign <- sum(results$Wilxocon_adjusted_p_value < alpha, na.rm = TRUE)
  t_test_sign <- sum(results$T_Test_adjusted_p_value < alpha, na.rm = TRUE)
  
  #print results
  cat("\n-------------------------------\n")
  cat("CV:", cv_label, "| Missingness:", pct_label, "\n")
  cat("Number of significant metabolites:\n")
  cat("- Wilcoxon test:", wilcox_sign, "\n")
  cat("- Paired t-test:", t_test_sign, "\n")
  cat("-------------------------------\n")
  
  
  
  #return results dataframe
  return(results)
}


#call function to compare Visit 1 vs. Visit 2
#original 
stats_original <- visit_statistical_tests(data)  #0 and 0

#Halfmin
halfmin_test_results <- list()

for (cv in names(halfmin_imputed_v1)) {
  halfmin_test_results[[cv]] <- list()
  
  for (pct in names(halfmin_imputed_v1[[cv]])) {
    combined_data <- bind_rows(
      halfmin_imputed_v1[[cv]][[pct]],
      halfmin_imputed_v2[[cv]][[pct]]
    )
    
    halfmin_test_results[[cv]][[pct]] <- visit_statistical_tests(
      combined_data,
      cv_label = cv,
      pct_label = pct
    )
  }
}

halfmin_test_summary <- bind_rows(
  lapply(names(halfmin_test_results), function(cv) {
    lapply(names(halfmin_test_results[[cv]]), function(pct) {
      res <- halfmin_test_results[[cv]][[pct]]
      data.frame(
        CV = cv,
        Missingness = as.numeric(gsub("pct", "", pct)),
        Method = "Half-min",
        Significant_t = sum(res$T_Test_adjusted_p_value < 0.05),
        Significant_wilcox = sum(res$Wilxocon_adjusted_p_value < 0.05)
      )
    }) %>% bind_rows()
  })
)


#KNN
knn_test_results <- list()

for (cv in names(knn_imputed_v1)) {
  knn_test_results[[cv]] <- list()
  
  for (pct in names(knn_imputed_v1[[cv]])) {
    combined_data <- bind_rows(
      knn_imputed_v1[[cv]][[pct]],
      knn_imputed_v2[[cv]][[pct]]
    )
    
    knn_test_results[[cv]][[pct]] <- visit_statistical_tests(
      combined_data,
      cv_label = cv,
      pct_label = pct
    )
  }
}

knn_test_summary <- bind_rows(
  lapply(names(knn_test_results), function(cv) {
    lapply(names(knn_test_results[[cv]]), function(pct) {
      res <- knn_test_results[[cv]][[pct]]
      data.frame(
        CV = cv,
        Missingness = as.numeric(gsub("pct", "", pct)),
        Method = "KNN",
        Significant_t = sum(res$T_Test_adjusted_p_value < 0.05),
        Significant_wilcox = sum(res$Wilxocon_adjusted_p_value < 0.05)
      )
    }) %>% bind_rows()
  })
)


#RF
rf_test_results <- list()

for (cv in names(rf_imputed_v1)) {
  rf_test_results[[cv]] <- list()
  
  for (pct in names(rf_imputed_v1[[cv]])) {
    combined_data <- bind_rows(
      rf_imputed_v1[[cv]][[pct]],
      rf_imputed_v2[[cv]][[pct]]
    )
    
    rf_test_results[[cv]][[pct]] <- visit_statistical_tests(
      combined_data,
      cv_label = cv,
      pct_label = pct
    )
  }
}

rf_test_summary <- bind_rows(
  lapply(names(rf_test_results), function(cv) {
    lapply(names(rf_test_results[[cv]]), function(pct) {
      res <- rf_test_results[[cv]][[pct]]
      data.frame(
        CV = cv,
        Missingness = as.numeric(gsub("pct", "", pct)),
        Method = "RF",
        Significant_t = sum(res$T_Test_adjusted_p_value < 0.05),
        Significant_wilcox = sum(res$Wilxocon_adjusted_p_value < 0.05)
      )
    }) %>% bind_rows()
  })
)


#QRILC
qrilc_test_results <- list()

for (cv in names(qrilc_imputed_v1)) {
  qrilc_test_results[[cv]] <- list()
  
  for (pct in names(qrilc_imputed_v1[[cv]])) {
    combined_data <- bind_rows(
      qrilc_imputed_v1[[cv]][[pct]],
      qrilc_imputed_v2[[cv]][[pct]]
    )
    
    qrilc_test_results[[cv]][[pct]] <- visit_statistical_tests(
      combined_data,
      cv_label = cv,
      pct_label = pct
    )
  }
}

qrilc_test_summary <- bind_rows(
  lapply(names(qrilc_test_results), function(cv) {
    lapply(names(qrilc_test_results[[cv]]), function(pct) {
      res <- qrilc_test_results[[cv]][[pct]]
      data.frame(
        CV = cv,
        Missingness = as.numeric(gsub("pct", "", pct)),
        Method = "QRILC",
        Significant_t = sum(res$T_Test_adjusted_p_value < 0.05),
        Significant_wilcox = sum(res$Wilxocon_adjusted_p_value < 0.05)
      )
    }) %>% bind_rows()
  })
)

#add original to each summary
halfmin_test_summary <- rbind(
  data.frame(
    CV = NA,
    Missingness = "Original",
    Method = "Half-min",
    Significant_t = 0,
    Significant_wilcox = 0
  ),
  halfmin_test_summary
)

knn_test_summary <- rbind(
  data.frame(
    CV = NA,
    Missingness = "Original",
    Method = "KNN",
    Significant_t = 0,
    Significant_wilcox = 0
  ),
  knn_test_summary
)

rf_test_summary <- rbind(
  data.frame(
    CV = NA,
    Missingness = "Original",
    Method = "RF",
    Significant_t = 0,
    Significant_wilcox = 0
  ),
  rf_test_summary
)

qrilc_test_summary <- rbind(
  data.frame(
    CV = NA,
    Missingness = "Original",
    Method = "QRILC",
    Significant_t = 0,
    Significant_wilcox = 0
  ),
  qrilc_test_summary
)

#combine overall stats 
combined_test_summary <- bind_rows(
  halfmin_test_summary,
  knn_test_summary,
  rf_test_summary,
  qrilc_test_summary
) %>%
  mutate(
    Missingness = as.character(Missingness),
    Missingness = ifelse(Missingness == "Original", "Original", paste0(Missingness, "%"))
  )

#pivot for plotting
stats_long <- combined_test_summary %>%
  pivot_longer(
    cols = c(Significant_t, Significant_wilcox),
    names_to = "Test",
    values_to = "Significant_Lipids"
  ) %>%
  mutate(
    Test = recode(Test,
                  Significant_t = "T-test",
                  Significant_wilcox = "Wilcoxon"),
    Missingness = factor(Missingness, levels = c("Original", "10%", "20%", "30%", "40%"))
  )

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/T_Test_Wilcox.pdf", width = 14, height = 10)

#ceate grouped bar plot
cv_levels <- stats_long %>%
  filter(!is.na(CV)) %>%
  pull(CV) %>%
  unique()

for (cv in cv_levels) {
  plot_data <- stats_long %>% filter(CV == cv | is.na(CV))
  
  p <- ggplot(plot_data, aes(x = Missingness, y = Significant_Lipids, fill = Test)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Significant_Lipids), 
              position = position_dodge(width = 0.9), 
              vjust = -0.3, size = 5, fontface = "bold") +  
    facet_wrap(~ Method, nrow = 1) +  
    labs(
      title = paste("Significant Lipids (T-test vs Wilcoxon) -", toupper(cv)),
      x = "Missingness Percentage",
      y = "Number of Significant Lipids",
      fill = "Statistical Test"
    ) +
    theme_minimal(base_size = 16) +  
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  
      axis.text.y = element_text(size = 12, face = "bold"),  
      axis.title.x = element_text(size = 14, face = "bold"),  
      axis.title.y = element_text(size = 14, face = "bold"),  
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      legend.title = element_text(size = 14, face = "bold"),  
      legend.text = element_text(size = 12), 
      strip.text = element_text(size = 14, face = "bold"), 
      panel.spacing = unit(2, "lines")
    )
  
  print(p)  # or save with ggsave if needed
}


dev.off()

# ---------------
# TITLE: NRMSE
# ---------------

#function to calcualte NRMSE
calculate_weighted_nrmse <- function(original, imputed, method, percentage){
  #numeric columns 
  numeric_col_names <- colnames(original)[6:ncol(original)]
  
  #calculate nrmse for each column
  nrmse_values <- sapply(numeric_col_names, function(col){
    actual_val <- original[[col]]
    imputed_val <- imputed[[col]]
    
    #ensure no missing value 
    valid_indices <- !is.na(actual_val) & !is.na(imputed_val)
    
    if (sum(valid_indices) > 2) { #if enough data
      mse <- mean((actual_val[valid_indices] - imputed_val[valid_indices])^2) #mean squared error
      rmse <- sqrt(mse) #root mean squared error
      norm_factor <- max(actual_val[valid_indices], na.rm = TRUE) - min(actual_val[valid_indices], na.rm = TRUE)
      
      if (norm_factor > 0){
        nrmse <- rmse / norm_factor
        weighted_nrmse <- nrmse * percentage
        return(weighted_nrmse)
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  })
  
  return(data.frame(
    Lipid = numeric_col_names,
    Imputation_Method = method,
    MNAR_proportion = (percentage * 100),
    Weighted_NRMSE = nrmse_values
  ))
}

#helper function to comppute weighted nrmse
compute_nrmse_all <- function(imputed_list, original_list, method_name, visit_label = "v1") {
  nrmse_results <- list()
  
  for (cv in names(imputed_list)) {
    nrmse_results[[cv]] <- list()
    
    for (pct in names(imputed_list[[cv]])) {
      imputed_data <- imputed_list[[cv]][[pct]]
      original_data <- original_list[[cv]]
      
      percentage_value <- as.numeric(gsub("pct", "", pct)) / 100
      
      nrmse_df <- calculate_weighted_nrmse(
        original = original_data,
        imputed = imputed_data,
        method = method_name,
        percentage = percentage_value
      )
      
      # Add metadata
      nrmse_df$CV <- cv
      nrmse_df$Missingness <- percentage_value * 100
      nrmse_df$Visit <- visit_label
      
      nrmse_results[[cv]][[pct]] <- nrmse_df
    }
  }
  
  # Combine all results
  bind_rows(
    lapply(nrmse_results, function(cv_group) {
      bind_rows(cv_group)
    })
  )
}

#call function to compute NRMSE
#Halfmin
#v1
nrmse_halfmin_v1 <- compute_nrmse_all(
  imputed_list = halfmin_imputed_v1,
  original_list = original_data_v1,
  method_name = "Half-min",
  visit_label = "v1"
)
#v1
nrmse_halfmin_v2 <- compute_nrmse_all(
  imputed_list = halfmin_imputed_v2,
  original_list = original_data_v2,
  method_name = "Half-min",
  visit_label = "v2"
)

#KNN
#v1
nrmse_knn_v1 <- compute_nrmse_all(
  imputed_list = knn_imputed_v1,
  original_list = original_data_v1,
  method_name = "KNN",
  visit_label = "v1"
)
#v1
nrmse_knn_v2 <- compute_nrmse_all(
  imputed_list = knn_imputed_v2,
  original_list = original_data_v2,
  method_name = "KNN",
  visit_label = "v2"
)

#RF
#v1
nrmse_rf_v1 <- compute_nrmse_all(
  imputed_list = rf_imputed_v1,
  original_list = original_data_v1,
  method_name = "RF",
  visit_label = "v1"
)
#v1
nrmse_rf_v2 <- compute_nrmse_all(
  imputed_list = rf_imputed_v2,
  original_list = original_data_v2,
  method_name = "RF",
  visit_label = "v2"
)

#QRILC
#v1
nrmse_qrilc_v1 <- compute_nrmse_all(
  imputed_list = qrilc_imputed_v1,
  original_list = original_data_v1,
  method_name = "QRILC",
  visit_label = "v1"
)
#v1
nrmse_qrilc_v2 <- compute_nrmse_all(
  imputed_list = qrilc_imputed_v2,
  original_list = original_data_v2,
  method_name = "QRILC",
  visit_label = "v2"
)

# --------------------
# Part 1: NRMSE Plot 
# --------------------

#v1
nrmse_data_v1 <- bind_rows(
  nrmse_halfmin_v1,
  nrmse_knn_v1,
  nrmse_rf_v1,
  nrmse_qrilc_v1 
)

#v2
nrmse_data_v2 <- bind_rows(
  nrmse_halfmin_v2,
  nrmse_knn_v2,
  nrmse_rf_v2,
  nrmse_qrilc_v2 
)

#convert MNAR to factor
nrmse_data_v1$MNAR_proportion <- factor(nrmse_data_v1$MNAR_proportion, levels = c(10, 20, 30, 40))
nrmse_data_v2$MNAR_proportion <- factor(nrmse_data_v2$MNAR_proportion, levels = c(10, 20, 30, 40))

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/NRMSE_v1.pdf", width = 14, height = 10)

cv_levels <- unique(nrmse_data_v1$CV)

#V1
for (cv in cv_levels) {
  
  plot_data <- nrmse_data_v1 %>% filter(CV == cv)
  
  p <- ggplot(plot_data, aes(x = MNAR_proportion, y = Weighted_NRMSE, fill = Imputation_Method)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    scale_fill_manual(values = c("lightblue", "orange", "blue", "magenta")) +
    labs(
      title = paste("Weighted NRMSE by Imputation Method -", toupper(cv), "(Visit 1)"),
      x = "MNAR Proportion (%)",
      y = "Weighted NRMSE",
      fill = "Imputation Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold")
    ) +
    ylim(0, 0.4)
  
  print(p)  
}

dev.off()

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/NRMSE_v2.pdf", width = 14, height = 10)

cv_levels <- unique(nrmse_data_v2$CV)

#V1
for (cv in cv_levels) {
  
  plot_data <- nrmse_data_v2 %>% filter(CV == cv)
  
  p <- ggplot(plot_data, aes(x = MNAR_proportion, y = Weighted_NRMSE, fill = Imputation_Method)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    scale_fill_manual(values = c("lightblue", "orange", "blue", "magenta")) +
    labs(
      title = paste("Weighted NRMSE by Imputation Method -", toupper(cv), "(Visit 1)"),
      x = "MNAR Proportion (%)",
      y = "Weighted NRMSE",
      fill = "Imputation Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold")
    ) +
    ylim(0, 0.4)
  
  print(p)  
}

dev.off()


# ---------------------------------------
# TITLE: Normalized Mean Difference (NMD)
# ---------------------------------------

#function to calculate normalized mean difference 
norm_mean_diff <- function(original, imputed, method, percentage, visit){
  #numeric columns
  numeric_original <- original[, 6:ncol(original)]
  numeric_imputed <- imputed[, 6:ncol(imputed)]
  
  #mean before imputation
  mean_before <- numeric_original %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Lipid", values_to = "Mean_Before")
  
  #mean after imputation
  mean_after <- numeric_imputed %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(cols = everything(), names_to = "Lipid", values_to = "Mean_After")
  
  #merge before and after
  mean_comparison <- merge(mean_before, mean_after, by = "Lipid")
  
  #compute NMD
  mean_comparison <- mean_comparison %>%
    mutate(Normalized_Difference = (Mean_After - Mean_Before) / Mean_Before)
  
  plot_title <- paste0("Normalized Difference with ", percentage, "% Missing Values using ", method, " Imputation ", "(Visit ", visit, ")")
  
  
  #plot the density of the normalized difference
  plot <- ggplot(mean_comparison, aes(x = Normalized_Difference)) +
    geom_density(fill = "blue", alpha = 0.4, color = "black") + 
    theme_minimal() +
    labs(title = plot_title,  
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red")  
  
  print(plot)  
  return(mean_comparison)  
}

#helper function
compute_nmd_all <- function(imputed_list, original_list, method_name, visit_label = "1") {
  all_nmd <- list()
  
  for (cv in names(imputed_list)) {
    all_nmd[[cv]] <- list()
    
    for (pct in names(imputed_list[[cv]])) {
      imputed_data <- imputed_list[[cv]][[pct]]
      original_data <- original_list[[cv]]
      
      percentage_value <- as.numeric(gsub("pct", "", pct))
      
      nmd_df <- norm_mean_diff(
        original = original_data,
        imputed = imputed_data,
        method = method_name,
        percentage = percentage_value,
        visit = visit_label
      )
      
      # Add metadata
      nmd_df$CV <- cv
      nmd_df$Missingness <- percentage_value
      nmd_df$Method <- method_name
      nmd_df$Visit <- paste0("v", visit_label)
      
      all_nmd[[cv]][[pct]] <- nmd_df
    }
  }
  
  # Combine into one dataframe
  bind_rows(lapply(all_nmd, bind_rows))
}


#call normalized difference function
#v1
nmd_halfmin_v1 <- compute_nmd_all(halfmin_imputed_v1, original_data_v1, method_name = "Half-min", visit_label = "1")
nmd_knn_v1     <- compute_nmd_all(knn_imputed_v1,     original_data_v1, method_name = "KNN",     visit_label = "1")
nmd_rf_v1      <- compute_nmd_all(rf_imputed_v1,      original_data_v1, method_name = "RF",      visit_label = "1")
nmd_qrilc_v1   <- compute_nmd_all(qrilc_imputed_v1,   original_data_v1, method_name = "QRILC",   visit_label = "1")

#v2
nmd_halfmin_v2 <- compute_nmd_all(halfmin_imputed_v2, original_data_v2, method_name = "Half-min", visit_label = "2")
nmd_knn_v2     <- compute_nmd_all(knn_imputed_v2,     original_data_v2, method_name = "KNN",     visit_label = "2")
nmd_rf_v2      <- compute_nmd_all(rf_imputed_v2,      original_data_v2, method_name = "RF",      visit_label = "2")
nmd_qrilc_v2   <- compute_nmd_all(qrilc_imputed_v2,   original_data_v2, method_name = "QRILC",   visit_label = "2")


# ---------------------------------------------
# Part 1: Plot of all NMD per Imputation Method
# --------------------------------------------

#combine all 
nmd_data_v1 <- bind_rows(nmd_halfmin_v1, nmd_knn_v1, nmd_rf_v1, nmd_qrilc_v1)
nmd_data_v2 <- bind_rows(nmd_halfmin_v2, nmd_knn_v2, nmd_rf_v2, nmd_qrilc_v2)

#CV and Percentage as factors
nmd_data_v1$CV <- factor(nmd_data_v1$CV, levels = c("cv20", "cv25", "cv30", "cv35", "cv40"))
nmd_data_v1$Percentage <- factor(paste0(nmd_data_v1$Missingness, "%"), levels = c("10%", "20%", "30%", "40%"))

nmd_data_v2$CV <- factor(nmd_data_v2$CV, levels = c("cv20", "cv25", "cv30", "cv35", "cv40"))
nmd_data_v2$Percentage <- factor(paste0(nmd_data_v2$Missingness, "%"), levels = c("10%", "20%", "30%", "40%"))

#plot for visit 1
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/NMD_v1.pdf", width = 14, height = 10)

cv_levels <- levels(nmd_data_v1$CV)

for (cv in cv_levels) {
  plot_data <- nmd_data_v1 %>% filter(CV == cv)
  
  p <- ggplot(plot_data, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(title = paste("Normalized Difference by Method -", toupper(cv), "(Visit 1)"),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~Method, scales = "free") +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  print(p)
}

dev.off()

#plot for visit 2
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/NMD_v2.pdf", width = 14, height = 10)

cv_levels <- levels(nmd_data_v2$CV)

for (cv in cv_levels) {
  plot_data <- nmd_data_v2 %>% filter(CV == cv)
  
  p <- ggplot(plot_data, aes(x = Normalized_Difference, fill = Percentage, color = Percentage)) +
    geom_density(alpha = 0.3) +
    theme_minimal() +
    labs(title = paste("Normalized Difference by Method -", toupper(cv), "(Visit 1)"),
         x = "Normalized Difference",
         y = "Density") +
    xlim(-0.2, 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~Method, scales = "free") +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  print(p)
}

dev.off()

# ------------------------------------
# TITLE: Distribution Plot (Whole Data)
# ------------------------------------

#function to plot distribution before and after imputation (entire dataset)
plot_whole_distribution <- function(original, imputed, method, percentage, visit, cv) {
  numeric_original <- original[, 6:ncol(original)]
  numeric_imputed <- imputed[, 6:ncol(imputed)]
  
  original_long <- numeric_original %>%
    pivot_longer(cols = everything(), names_to = "Lipid", values_to = "Value") %>%
    mutate(Data = "Original Data")
  
  imputed_long <- numeric_imputed %>%
    pivot_longer(cols = everything(), names_to = "Lipid", values_to = "Value") %>%
    mutate(Data = "Imputed Data")
  
  imputed_values <- numeric_original != numeric_imputed
  
  imputed_only_long <- numeric_imputed %>%
    as.data.frame() %>%
    replace(!imputed_values, NA) %>%
    pivot_longer(cols = everything(), names_to = "Lipid", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    mutate(Data = "Imputed Values")
  
  combined_data <- bind_rows(original_long, imputed_long, imputed_only_long)
  
  mean_data <- combined_data %>%
    group_by(Data) %>%
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(combined_data, aes(x = Value, fill = Data)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    labs(title = paste0("Overall Density: ", method, " Imputation (", percentage, "%, ", toupper(cv), ", Visit ", visit, ")"),
         x = "Value",
         y = "Density") +
    geom_vline(data = mean_data %>% filter(Data == "Original Data"),
               aes(xintercept = mean_value, color = "Original Mean"), linewidth = 0.5, linetype = "dashed") +
    geom_vline(data = mean_data %>% filter(Data == "Imputed Data"),
               aes(xintercept = mean_value, color = "Imputed Mean"), linewidth = 0.5, linetype = "dashed") +
    scale_color_manual(name = "Mean Values",
                       values = c("Original Mean" = "blue", "Imputed Mean" = "red")) +
    xlim(-100, 400) +
    theme(legend.position = "bottom")
  
  print(p)
}

#Halfmin
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/Distr_Halfmin.pdf", width = 14, height = 10)

#v1
for (cv in names(halfmin_imputed_v1)) {
  for (pct in names(halfmin_imputed_v1[[cv]])) {
    imputed_data <- halfmin_imputed_v1[[cv]][[pct]]
    original_data <- original_data_v1[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "Half-min",
      percentage = percentage_value,
      visit = 1,
      cv = cv
    )
  }
}

#v2
for (cv in names(halfmin_imputed_v2)) {
  for (pct in names(halfmin_imputed_v2[[cv]])) {
    imputed_data <- halfmin_imputed_v2[[cv]][[pct]]
    original_data <- original_data_v2[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "Half-min",
      percentage = percentage_value,
      visit = 2,
      cv = cv
    )
  }
}

dev.off()

#KNN
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/Distr_KNN.pdf", width = 14, height = 10)

#v1
for (cv in names(knn_imputed_v1)) {
  for (pct in names(knn_imputed_v1[[cv]])) {
    imputed_data <- knn_imputed_v1[[cv]][[pct]]
    original_data <- original_data_v1[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "KNN",
      percentage = percentage_value,
      visit = 1,
      cv = cv
    )
  }
}

#v2
for (cv in names(knn_imputed_v2)) {
  for (pct in names(knn_imputed_v2[[cv]])) {
    imputed_data <- knn_imputed_v2[[cv]][[pct]]
    original_data <- original_data_v2[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "KNN",
      percentage = percentage_value,
      visit = 2,
      cv = cv
    )
  }
}

dev.off()

#KNN
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/Distr_RF.pdf", width = 14, height = 10)

#v1
for (cv in names(rf_imputed_v1)) {
  for (pct in names(rf_imputed_v1[[cv]])) {
    imputed_data <- rf_imputed_v1[[cv]][[pct]]
    original_data <- original_data_v1[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "RF",
      percentage = percentage_value,
      visit = 1,
      cv = cv
    )
  }
}

#v2
for (cv in names(rf_imputed_v2)) {
  for (pct in names(rf_imputed_v2[[cv]])) {
    imputed_data <- rf_imputed_v2[[cv]][[pct]]
    original_data <- original_data_v2[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "RF",
      percentage = percentage_value,
      visit = 2,
      cv = cv
    )
  }
}

dev.off()

#QRILC
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/Distr_qrilc.pdf", width = 14, height = 10)

#v1
for (cv in names(qrilc_imputed_v1)) {
  for (pct in names(qrilc_imputed_v1[[cv]])) {
    imputed_data <- qrilc_imputed_v1[[cv]][[pct]]
    original_data <- original_data_v1[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "QRILC",
      percentage = percentage_value,
      visit = 1,
      cv = cv
    )
  }
}

#v2
for (cv in names(qrilc_imputed_v2)) {
  for (pct in names(qrilc_imputed_v2[[cv]])) {
    imputed_data <- qrilc_imputed_v2[[cv]][[pct]]
    original_data <- original_data_v2[[cv]]
    percentage_value <- as.numeric(gsub("pct", "", pct))
    
    plot_whole_distribution(
      original = original_data,
      imputed = imputed_data,
      method = "QRILC",
      percentage = percentage_value,
      visit = 2,
      cv = cv
    )
  }
}

dev.off()


# ------------------------------------
# TITLE: ANOVA
# ------------------------------------

#fit an ANOVA model with an interaction term between imputation method and missingness level 
#apply log to normalize data

#ANOVA 
#Visit 1
anova_all_v1 <- aov(log(Weighted_NRMSE) ~ Imputation_Method * MNAR_proportion * CV, data = nrmse_data_v1)
summary(anova_all_v1)

#Visit 2
anova_all_v2 <- aov(log(Weighted_NRMSE) ~ Imputation_Method * MNAR_proportion * CV, data = nrmse_data_v2)
summary(anova_all_v2)

# --------------------------------------------------------
# Part 1: Check Residuals and Normality for ANOVA result
# -------------------------------------------------------

#histogram of residuals (symmetric and bell suggests normality)
#v1
ggplot(data.frame(residuals = residuals(anova_all_v1)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue") +
  labs(title = "Histogram of Residuals (Visit 1)", x = "Residuals", y = "Frequency") +
  theme_minimal()
#v2
ggplot(data.frame(residuals = residuals(anova_all_v2)), aes(x = residuals)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue") +
  labs(title = "Histogram of Residuals (Visit 2)", x = "Residuals", y = "Frequency") +
  theme_minimal()


#QQ Plot
#plot residual against theoretical normal distirbutions
#red line shows perfect normal distirbution
#v1
qqnorm(residuals(anova_all_v1))
qqline(residuals(anova_all_v1), col = "red")
#v2
qqnorm(residuals(anova_all_v2))
qqline(residuals(anova_all_v2), col = "red")

#Tukey-Anscombe to check residuals vs fitted values
#x-axis = fitted values (predicted) and y-axis = residueals (error)
#residuals should be evenly spread

#v1
plot(fitted(anova_all_v1), resid(anova_all_v1), 
     main = "Tukey-Anscombe Plot (Visit 1)", 
     xlab = "Fitted", ylab = "Residuals", col = "blue")
#v2
plot(fitted(anova_all_v2), resid(anova_all_v2), 
     main = "Tukey-Anscombe Plot (Visit 2)",
     xlab = "Fitted", ylab = "Residuals", col = "blue")


# ------------------------------------
# Part 7: Kruskal-Wallis
# ------------------------------------

#perform kruskal-wallis test
#main effects only
#v1
kruskal.test(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data_v1)
kruskal.test(Weighted_NRMSE ~ MNAR_proportion, data = nrmse_data_v1)
kruskal.test(Weighted_NRMSE ~ CV, data = nrmse_data_v1)
#v2
kruskal.test(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data_v2)
kruskal.test(Weighted_NRMSE ~ MNAR_proportion, data = nrmse_data_v2)
kruskal.test(Weighted_NRMSE ~ CV, data = nrmse_data_v2)

# ------------------------------------
# Part 8: Dunn's Test for each imputation
# ------------------------------------

#perform Dunn's Test for pairwise comparison (BH correction for multiple testing)

#v1
dunn_method_v1 <- dunnTest(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data_v1, method = "bh")
dunn_mnar_v1   <- dunnTest(Weighted_NRMSE ~ MNAR_proportion, data = nrmse_data_v1, method = "bh")
dunn_cv_v1     <- dunnTest(Weighted_NRMSE ~ CV, data = nrmse_data_v1, method = "bh")

#v2
dunn_method_v2 <- dunnTest(Weighted_NRMSE ~ Imputation_Method, data = nrmse_data_v2, method = "bh")
dunn_mnar_v2   <- dunnTest(Weighted_NRMSE ~ MNAR_proportion, data = nrmse_data_v2, method = "bh")
dunn_cv_v2     <- dunnTest(Weighted_NRMSE ~ CV, data = nrmse_data_v2, method = "bh")

#view result tables
#v1
print(dunn_method_v1)
print(dunn_mnar_v1)
print(dunn_cv_v1)
#v2
print(dunn_method_v2)
print(dunn_mnar_v2)
print(dunn_cv_v2)

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/ANOVA_v1.pdf", width = 14, height = 10)

#plot results
#v1
ggplot(nrmse_data_v1, aes(x = Imputation_Method, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot() +
  labs(title = "Visit 1: NRMSE by Imputation Method", y = "Weighted NRMSE") +
  ylim(0,0.5) +
  theme_minimal()

ggplot(nrmse_data_v1, aes(x = MNAR_proportion, y = Weighted_NRMSE, fill = MNAR_proportion)) +
  geom_boxplot() +
  labs(title = "Visit 1: NRMSE by Missingness %", y = "Weighted NRMSE") +
  ylim(0,0.5) +
  theme_minimal()

ggplot(nrmse_data_v1, aes(x = CV, y = Weighted_NRMSE, fill = CV)) +
  geom_boxplot() +
  labs(title = "Visit 1: NRMSE by CV Threshold", y = "Weighted NRMSE") +
  ylim(0,0.5) +
  theme_minimal()

dev.off()

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/with_CV/ANOVA_v2.pdf", width = 14, height = 10)

#v2
ggplot(nrmse_data_v2, aes(x = Imputation_Method, y = Weighted_NRMSE, fill = Imputation_Method)) +
  geom_boxplot() +
  labs(title = "Visit 2: NRMSE by Imputation Method", y = "Weighted NRMSE") +
  ylim(0,0.5) +
  theme_minimal()

ggplot(nrmse_data_v2, aes(x = MNAR_proportion, y = Weighted_NRMSE, fill = MNAR_proportion)) +
  geom_boxplot() +
  labs(title = "Visit 2: NRMSE by Missingness %", y = "Weighted NRMSE") +
  ylim(0,0.5) +
  theme_minimal()

ggplot(nrmse_data_v2, aes(x = CV, y = Weighted_NRMSE, fill = CV)) +
  geom_boxplot() +
  labs(title = "Visit 2: NRMSE by CV Threshold", y = "Weighted NRMSE") +
  ylim(0,0.5) +
  theme_minimal()

dev.off()



