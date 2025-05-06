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
library(readxl)

#load data
data_whole <- read_excel("/Users/marcinebessire/Desktop/Master_Thesis/Papers_Data/FattyLiver_Pauling/mmc1.xlsx")

# --------------------------
# TITLE: Preprocessing
# --------------------------

# ------------------------
# Part 1: Data cleaning
# ------------------------

#remove unnecessary rows on top
data_whole <- data_whole[-c(1:9),]

#now transpose (patient as row and lipid as column) 
data_whole <- as.data.frame(t(data_whole))

#set first row as column names
colnames(data_whole) <- as.character(unlist(data_whole[1, ]))

#remove first row (now used as header)
data_whole <- data_whole[-1, ]

#reset row indices (1,2,3 etc)
rownames(data_whole) <- NULL 

#remove specie column and count NAs
data_whole <- data_whole[,-4]
data_whole <- data_whole[,-2]

#set data to numeric 
data_whole <- data_whole %>%
  mutate(across(3:ncol(.), as.numeric))

str(data_whole)

# ----------------------------------------
# Part 2: Remove Lipids that have NA
# ----------------------------------------

#remove cols with any NA
final_data <- data_whole %>%
  select(where(~ all(!is.na(.)))) #34 lipids left out of 316

# ----------------------------------------
# Part 3: Order by Condition
# ----------------------------------------

data_NC <- subset(final_data, Condition == "NC") #49 patients
rownames(data_NC) <- NULL 
data_HO <- subset(final_data, Condition == "HO") #51 patients
rownames(data_HO) <- NULL 
data_NAFL <- subset(final_data, Condition == "NAFL") #143 patients
rownames(data_NAFL) <- NULL 
data_NASH <- subset(final_data, Condition == "NASH") #94 patients
rownames(data_NASH) <- NULL 

# ---------------------
# TITLE: Calculate CV 
# ---------------------

# ---------------------
# Part 1: Calculate CV
# ---------------------

#function to calcualte CV
calculate_cv_dataframe <- function(data) {
  #numeric
  numeric_data <- data[, 3:ncol(data)]
  
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
cv_total <- calculate_cv_dataframe(final_data)
cv_NC <- calculate_cv_dataframe(data_NC)
cv_HO <- calculate_cv_dataframe(data_HO)
cv_NAFL <- calculate_cv_dataframe(data_NAFL)
cv_NASH <- calculate_cv_dataframe(data_NASH)

# ---------------------
# Part 2: Calculate SD
# ---------------------

#function to calcualte SD
calculate_sd_dataframe <- function(data) {
  #numeric
  numeric_data <- data[, 3:ncol(data)]
  
  #calculate sd per column
  sd_values <- sapply(numeric_data, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    sd(x, na.rm = TRUE)
  })
  
  #create a SD results dataframe
  sd_df <- data.frame(
    Column = names(sd_values),
    SD = round(sd_values, 2)
  )
  
  return(sd_df)
}

#call function
sd_total <- calculate_sd_dataframe(final_data)
sd_NC <- calculate_sd_dataframe(data_NC)
sd_HO <- calculate_sd_dataframe(data_HO)
sd_NAFL <- calculate_sd_dataframe(data_NAFL)
sd_NASH <- calculate_sd_dataframe(data_NASH)


# ---------------------------------
# TITLE: Plot Whole Distribution
# ---------------------------------

#function to plot the overall distribution of all numeric values
plot_overall_distribution <- function(original, name = NA) {
  #numeric columns
  numeric_original <- original[, 3:ncol(original)]
  
  #convert the entire numeric data into a single vector
  all_values <- unlist(numeric_original, use.names = FALSE)
  
  #create a data frame for plotting
  df_all_values <- data.frame(Value = all_values)
  
  #plot the overall density distribution
  plot <- ggplot(df_all_values, aes(x = Value)) +
    geom_density(fill = "gray", alpha = 0.5, color = "black") +
    theme_minimal() +
    labs(title = paste0("Overall Distribution of NAFLD Data (", name, ")"),
         x = "Measurement",
         y = "Density") +
    xlim(-10,50)
  
  print(plot)
}

#call function
#original with NA
plot_overall_distribution(data_whole, name = "Original")
#without data w/o NA
plot_overall_distribution(final_data, name = "After Filtering")
#NC data
plot_overall_distribution(data_NC, name = "NC")
#HO
plot_overall_distribution(data_HO, name = "HO")
#NAFL
plot_overall_distribution(data_NAFL, name = "NAFL")
#NASH
plot_overall_distribution(data_NASH, name = "NASH")

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
  
  #iterate through numeric cols
  for (col in 3:ncol(data_copy)){
    #skip non-numeric
    if (!is.numeric(data_copy[[col]])) next
    
    #number of NA to introduce based on missingness
    num_mv <- round(nrow(data_copy) * missing_percentage)
    
    if (num_mv > 0){
      col_data <- data_copy[[col]]
      #find indices with lowest value
      ordered_indices <- order(col_data, na.last = NA)
      #select indices with lwoest value
      missing_indices <- ordered_indices[1:min(num_mv, length(ordered_indices))]
      
      #set values to NA
      data_copy[missing_indices, col] <- NA
    }
  }
  
  return(data_copy)
}

#call function to generate MNAR (evenly between visit 1 and visit 2)
#NC
NC_10pct <- MNAR_simulation(data_NC, 0.10) #10% missing values (4 MV)
NC_15pct <- MNAR_simulation(data_NC, 0.15) #15% missing values (7 MV)
NC_20pct <- MNAR_simulation(data_NC, 0.20) #20% missing values (10 MV)
NC_25pct <- MNAR_simulation(data_NC, 0.25) #25% missing values (12 MV)
NC_30pct <- MNAR_simulation(data_NC, 0.30) #30% missing values (15 MV)
NC_40pct <- MNAR_simulation(data_NC, 0.40) #40% missing values (20 MV)
#HO
HO_10pct <- MNAR_simulation(data_HO, 0.10) #10% missing values (5 MV)
HO_15pct <- MNAR_simulation(data_HO, 0.15) #15% missing values (8 MV)
HO_20pct <- MNAR_simulation(data_HO, 0.20) #20% missing values (10 MV)
HO_25pct <- MNAR_simulation(data_HO, 0.25) #25% missing values (13 MV)
HO_30pct <- MNAR_simulation(data_HO, 0.30) #30% missing values (15 MV)
HO_40pct <- MNAR_simulation(data_HO, 0.40) #40% missing values (20 MV)
#NAFL
NAFL_10pct <- MNAR_simulation(data_NAFL, 0.10) #10% missing values (14 MV)
NAFL_15pct <- MNAR_simulation(data_NAFL, 0.15) #15% missing values (21 MV)
NAFL_20pct <- MNAR_simulation(data_NAFL, 0.20) #20% missing values (29 MV)
NAFL_25pct <- MNAR_simulation(data_NAFL, 0.25) #25% missing values (36 MV)
NAFL_30pct <- MNAR_simulation(data_NAFL, 0.30) #30% missing values (43 MV)
NAFL_40pct <- MNAR_simulation(data_NAFL, 0.40) #40% missing values (57 MV)
#NASH
NASH_10pct <- MNAR_simulation(data_NASH, 0.10) #10% missing values (9 MV)
NASH_15pct <- MNAR_simulation(data_NASH, 0.15) #15% missing values (14 MV)
NASH_20pct <- MNAR_simulation(data_NASH, 0.20) #20% missing values (19 MV)
NASH_25pct <- MNAR_simulation(data_NASH, 0.25) #25% missing values (24 MV)
NASH_30pct <- MNAR_simulation(data_NASH, 0.30) #30% missing values (28 MV)
NASH_40pct <- MNAR_simulation(data_NASH, 0.40) #40% missing values (38 MV)


#to count NAs: 
#sum(is.na(NASH_40pct[["Chol 27:0"]]))


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
  metadata <- data_copy[,1:2]
  
  #numeric
  num_data <- data_copy[, 3:ncol(data_copy)]
  
  #loop through each column/lipid
  for (col in names(num_data)){
    min_val <- min(num_data[[col]], na.rm = TRUE) #find min value
    num_data[[col]][is.na(num_data[[col]])] <- 0.5 * min_val
  }
  
  final_data <- cbind(metadata, num_data)
  
  return(final_data)
  
}

#call function for Half-min 
#NC 
halfmin_10pct_NC <- half_min_imputation(NC_10pct)
halfmin_15pct_NC <- half_min_imputation(NC_15pct)
halfmin_20pct_NC <- half_min_imputation(NC_20pct)
halfmin_25pct_NC <- half_min_imputation(NC_25pct)
halfmin_30pct_NC <- half_min_imputation(NC_30pct)
halfmin_40pct_NC <- half_min_imputation(NC_40pct)
#HO
halfmin_10pct_HO <- half_min_imputation(HO_10pct)
halfmin_15pct_HO <- half_min_imputation(HO_15pct)
halfmin_20pct_HO <- half_min_imputation(HO_20pct)
halfmin_25pct_HO <- half_min_imputation(HO_25pct)
halfmin_30pct_HO <- half_min_imputation(HO_30pct)
halfmin_40pct_HO <- half_min_imputation(HO_40pct)
#NAFL
halfmin_10pct_NAFL <- half_min_imputation(NAFL_10pct)
halfmin_15pct_NAFL <- half_min_imputation(NAFL_15pct)
halfmin_20pct_NAFL <- half_min_imputation(NAFL_20pct)
halfmin_25pct_NAFL <- half_min_imputation(NAFL_25pct)
halfmin_30pct_NAFL <- half_min_imputation(NAFL_30pct)
halfmin_40pct_NAFL <- half_min_imputation(NAFL_40pct)
#NASH
halfmin_10pct_NASH <- half_min_imputation(NASH_10pct)
halfmin_15pct_NASH <- half_min_imputation(NASH_15pct)
halfmin_20pct_NASH <- half_min_imputation(NASH_20pct)
halfmin_25pct_NASH <- half_min_imputation(NASH_25pct)
halfmin_30pct_NASH <- half_min_imputation(NASH_30pct)
halfmin_40pct_NASH <- half_min_imputation(NASH_40pct)

# --------------------------
# Part 2: KNN Imputation
# --------------------------

#function for KNN imputation
knn_imputation <- function(data){
  data_copy <- data
  
  #metadata
  metadata <- data_copy[, 1:2]
  
  #numeric 
  num_data <- data_copy[, 3:ncol(data_copy)]
  
  #first transfrom into matrix and perform knn
  imputed_data <- impute.knn(as.matrix(t(num_data)), rowmax = 0.5, colmax = 1)
  
  #transform back into dataframe
  imputed_df <- as.data.frame(t(imputed_data$data))
  
  final_data <- cbind(metadata, imputed_df)
  
  return(final_data)
  
}

#call function for KNN
#NC 
KNN_10pct_NC <- knn_imputation(NC_10pct)
KNN_15pct_NC <- knn_imputation(NC_15pct)
KNN_20pct_NC <- knn_imputation(NC_20pct)
KNN_25pct_NC <- knn_imputation(NC_25pct)
KNN_30pct_NC <- knn_imputation(NC_30pct)
KNN_40pct_NC <- knn_imputation(NC_40pct)
#HO
KNN_10pct_HO <- knn_imputation(HO_10pct)
KNN_15pct_HO <- knn_imputation(HO_15pct)
KNN_20pct_HO <- knn_imputation(HO_20pct)
KNN_25pct_HO <- knn_imputation(HO_25pct)
KNN_30pct_HO <- knn_imputation(HO_30pct)
KNN_40pct_HO <- knn_imputation(HO_40pct)
#NAFL
KNN_10pct_NAFL <- knn_imputation(NAFL_10pct)
KNN_15pct_NAFL <- knn_imputation(NAFL_15pct)
KNN_20pct_NAFL <- knn_imputation(NAFL_20pct)
KNN_25pct_NAFL <- knn_imputation(NAFL_25pct)
KNN_30pct_NAFL <- knn_imputation(NAFL_30pct)
KNN_40pct_NAFL <- knn_imputation(NAFL_40pct)
#NASH
KNN_10pct_NASH <- knn_imputation(NASH_10pct)
KNN_15pct_NASH <- knn_imputation(NASH_15pct)
KNN_20pct_NASH <- knn_imputation(NASH_20pct)
KNN_25pct_NASH <- knn_imputation(NASH_25pct)
KNN_30pct_NASH <- knn_imputation(NASH_30pct)
KNN_40pct_NASH <- knn_imputation(NASH_40pct)


# ------------------------------------
# Part 3: RF Imputation
# ------------------------------------

#function for random forest imputation
rf_imputation <- function(data){
  data_copy <- data
  
  #metadata
  metadata <- data_copy[,1:2]
  
  #numerci data
  num_data <- data_copy[,3:ncol(data_copy)]
  
  #apply missforest
  imputed_data <- missForest(num_data, maxiter = 10, ntree = 100)
  imputed_df <- as.data.frame(imputed_data$ximp) #ximp = imputed data matrix
  
  final_data <- cbind(metadata, imputed_df)
  
  return(final_data)
}

#call function for RF 
#NC 
RF_10pct_NC <- rf_imputation(NC_10pct)
RF_15pct_NC <- rf_imputation(NC_15pct)
RF_20pct_NC <- rf_imputation(NC_20pct)
RF_25pct_NC <- rf_imputation(NC_25pct)
RF_30pct_NC <- rf_imputation(NC_30pct)
RF_40pct_NC <- rf_imputation(NC_40pct)
#HO
RF_10pct_HO <- rf_imputation(HO_10pct)
RF_15pct_HO <- rf_imputation(HO_15pct)
RF_20pct_HO <- rf_imputation(HO_20pct)
RF_25pct_HO <- rf_imputation(HO_25pct)
RF_30pct_HO <- rf_imputation(HO_30pct)
RF_40pct_HO <- rf_imputation(HO_40pct)
#NAFL
RF_10pct_NAFL <- rf_imputation(NAFL_10pct)
RF_15pct_NAFL <- rf_imputation(NAFL_15pct)
RF_20pct_NAFL <- rf_imputation(NAFL_20pct)
RF_25pct_NAFL <- rf_imputation(NAFL_25pct)
RF_30pct_NAFL <- rf_imputation(NAFL_30pct)
RF_40pct_NAFL <- rf_imputation(NAFL_40pct)
#NASH
RF_10pct_NASH <- rf_imputation(NASH_10pct)
RF_15pct_NASH <- rf_imputation(NASH_15pct)
RF_20pct_NASH <- rf_imputation(NASH_20pct)
RF_25pct_NASH <- rf_imputation(NASH_25pct)
RF_30pct_NASH <- rf_imputation(NASH_30pct)
RF_40pct_NASH <- rf_imputation(NASH_40pct)

# --------------------------
# Part 4: QRILC Imputation
# --------------------------

#function for QRILC imputation
qrilc_imputation <- function(data){
  #copy data
  data_copy <- data
  
  #select only numeric data
  numeric_data <- data_copy[, 3:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:2]
  
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
#NC
QRILC_10pct_NC <- qrilc_imputation(NC_10pct)
QRILC_15pct_NC <- qrilc_imputation(NC_15pct)
QRILC_20pct_NC <- qrilc_imputation(NC_20pct)
QRILC_25pct_NC <- qrilc_imputation(NC_25pct)
QRILC_30pct_NC <- qrilc_imputation(NC_30pct)
QRILC_40pct_NC <- qrilc_imputation(NC_40pct)
#HO
QRILC_10pct_HO <- qrilc_imputation(HO_10pct)
QRILC_15pct_HO <- qrilc_imputation(HO_15pct)
QRILC_20pct_HO <- qrilc_imputation(HO_20pct)
QRILC_25pct_HO <- qrilc_imputation(HO_25pct)
QRILC_30pct_HO <- qrilc_imputation(HO_30pct)
QRILC_40pct_HO <- qrilc_imputation(HO_40pct)
#NAFL
QRILC_10pct_NAFL <- qrilc_imputation(NAFL_10pct)
QRILC_15pct_NAFL <- qrilc_imputation(NAFL_15pct)
QRILC_20pct_NAFL <- qrilc_imputation(NAFL_20pct)
QRILC_25pct_NAFL <- qrilc_imputation(NAFL_25pct)
QRILC_30pct_NAFL <- qrilc_imputation(NAFL_30pct)
QRILC_40pct_NAFL <- qrilc_imputation(NAFL_40pct)
#NASH
QRILC_10pct_NASH <- qrilc_imputation(NASH_10pct)
QRILC_15pct_NASH <- qrilc_imputation(NASH_15pct)
QRILC_20pct_NASH <- qrilc_imputation(NASH_20pct)
QRILC_25pct_NASH <- qrilc_imputation(NASH_25pct)
QRILC_30pct_NASH <- qrilc_imputation(NASH_30pct)
QRILC_40pct_NASH <- qrilc_imputation(NASH_40pct)







