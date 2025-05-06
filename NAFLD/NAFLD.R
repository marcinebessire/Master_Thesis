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

#combined 
halfmin_datasets <- list(
  NC_10 = halfmin_10pct_NC,
  NC_15 = halfmin_15pct_NC,
  NC_20 = halfmin_20pct_NC,
  NC_25 = halfmin_25pct_NC,
  NC_30 = halfmin_30pct_NC,
  NC_40 = halfmin_40pct_NC,
  
  HO_10 = halfmin_10pct_HO,
  HO_15 = halfmin_15pct_HO,
  HO_20 = halfmin_20pct_HO,
  HO_25 = halfmin_25pct_HO,
  HO_30 = halfmin_30pct_HO,
  HO_40 = halfmin_40pct_HO,
  
  NAFL_10 = halfmin_10pct_NAFL,
  NAFL_15 = halfmin_15pct_NAFL,
  NAFL_20 = halfmin_20pct_NAFL,
  NAFL_25 = halfmin_25pct_NAFL,
  NAFL_30 = halfmin_30pct_NAFL,
  NAFL_40 = halfmin_40pct_NAFL,
  
  NASH_10 = halfmin_10pct_NASH,
  NASH_15 = halfmin_15pct_NASH,
  NASH_20 = halfmin_20pct_NASH,
  NASH_25 = halfmin_25pct_NASH,
  NASH_30 = halfmin_30pct_NASH,
  NASH_40 = halfmin_40pct_NASH
)


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
#combined 
KNN_datasets <- list(
  NC_10 = KNN_10pct_NC,
  NC_15 = KNN_15pct_NC,
  NC_20 = KNN_20pct_NC,
  NC_25 = KNN_25pct_NC,
  NC_30 = KNN_30pct_NC,
  NC_40 = KNN_40pct_NC,
  
  HO_10 = KNN_10pct_HO,
  HO_15 = KNN_15pct_HO,
  HO_20 = KNN_20pct_HO,
  HO_25 = KNN_25pct_HO,
  HO_30 = KNN_30pct_HO,
  HO_40 = KNN_40pct_HO,
  
  NAFL_10 = KNN_10pct_NAFL,
  NAFL_15 = KNN_15pct_NAFL,
  NAFL_20 = KNN_20pct_NAFL,
  NAFL_25 = KNN_25pct_NAFL,
  NAFL_30 = KNN_30pct_NAFL,
  NAFL_40 = KNN_40pct_NAFL,
  
  NASH_10 = KNN_10pct_NASH,
  NASH_15 = KNN_15pct_NASH,
  NASH_20 = KNN_20pct_NASH,
  NASH_25 = KNN_25pct_NASH,
  NASH_30 = KNN_30pct_NASH,
  NASH_40 = KNN_40pct_NASH
)

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
#combined 
RF_datasets <- list(
  NC_10 = RF_10pct_NC,
  NC_15 = RF_15pct_NC,
  NC_20 = RF_20pct_NC,
  NC_25 = RF_25pct_NC,
  NC_30 = RF_30pct_NC,
  NC_40 = RF_40pct_NC,
  
  HO_10 = RF_10pct_HO,
  HO_15 = RF_15pct_HO,
  HO_20 = RF_20pct_HO,
  HO_25 = RF_25pct_HO,
  HO_30 = RF_30pct_HO,
  HO_40 = RF_40pct_HO,
  
  NAFL_10 = RF_10pct_NAFL,
  NAFL_15 = RF_15pct_NAFL,
  NAFL_20 = RF_20pct_NAFL,
  NAFL_25 = RF_25pct_NAFL,
  NAFL_30 = RF_30pct_NAFL,
  NAFL_40 = RF_40pct_NAFL,
  
  NASH_10 = RF_10pct_NASH,
  NASH_15 = RF_15pct_NASH,
  NASH_20 = RF_20pct_NASH,
  NASH_25 = RF_25pct_NASH,
  NASH_30 = RF_30pct_NASH,
  NASH_40 = RF_40pct_NASH
)

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
#combined 
QRILC_datasets <- list(
  NC_10 = QRILC_10pct_NC,
  NC_15 = QRILC_15pct_NC,
  NC_20 = QRILC_20pct_NC,
  NC_25 = QRILC_25pct_NC,
  NC_30 = QRILC_30pct_NC,
  NC_40 = QRILC_40pct_NC,
  
  HO_10 = QRILC_10pct_HO,
  HO_15 = QRILC_15pct_HO,
  HO_20 = QRILC_20pct_HO,
  HO_25 = QRILC_25pct_HO,
  HO_30 = QRILC_30pct_HO,
  HO_40 = QRILC_40pct_HO,
  
  NAFL_10 = QRILC_10pct_NAFL,
  NAFL_15 = QRILC_15pct_NAFL,
  NAFL_20 = QRILC_20pct_NAFL,
  NAFL_25 = QRILC_25pct_NAFL,
  NAFL_30 = QRILC_30pct_NAFL,
  NAFL_40 = QRILC_40pct_NAFL,
  
  NASH_10 = QRILC_10pct_NASH,
  NASH_15 = QRILC_15pct_NASH,
  NASH_20 = QRILC_20pct_NASH,
  NASH_25 = QRILC_25pct_NASH,
  NASH_30 = QRILC_30pct_NASH,
  NASH_40 = QRILC_40pct_NASH
)

# --------------------------
# Part 5: MICE Imputation
# --------------------------

#function for MICE imputation
mice_imputation <- function(data,  m = 5, seed = 123){
  data_copy <- data
  
  #select only numeric data
  numeric_data <- data_copy[, 3:ncol(data_copy)]
  
  #metadata
  meta_data <- data_copy[, 1:2]

  #impute with mice (m = 5 nr of multiple imputation, pmm =  Predictive mean matching, seed to make it reproducible)
  imputed <- mice(numeric_data, m = m, method = 'pmm', seed = seed)
  
  #extract first completed dataset
  completed_data <- complete(imputed, 1)
  
  # Combine metadata with imputed lipid data
  final_df <- cbind(meta_data, completed_data)
  
  return(final_df)
}

#call function for MICE
#NC
mice_10pct_NC <- mice_imputation(NC_10pct)
mice_15pct_NC <- mice_imputation(NC_15pct)
mice_20pct_NC <- mice_imputation(NC_20pct)
mice_25pct_NC <- mice_imputation(NC_25pct)
mice_30pct_NC <- mice_imputation(NC_30pct)
mice_40pct_NC <- mice_imputation(NC_40pct)
#HO
mice_10pct_HO <- mice_imputation(HO_10pct)
mice_15pct_HO <- mice_imputation(HO_15pct)
mice_20pct_HO <- mice_imputation(HO_20pct)
mice_25pct_HO <- mice_imputation(HO_25pct)
mice_30pct_HO <- mice_imputation(HO_30pct)
mice_40pct_HO <- mice_imputation(HO_40pct)
#NAFL
mice_10pct_NAFL <- mice_imputation(NAFL_10pct)
mice_15pct_NAFL <- mice_imputation(NAFL_15pct)
mice_20pct_NAFL <- mice_imputation(NAFL_20pct)
mice_25pct_NAFL <- mice_imputation(NAFL_25pct)
mice_30pct_NAFL <- mice_imputation(NAFL_30pct)
mice_40pct_NAFL <- mice_imputation(NAFL_40pct)
#NASH
mice_10pct_NASH <- mice_imputation(NASH_10pct)
mice_15pct_NASH <- mice_imputation(NASH_15pct)
mice_20pct_NASH <- mice_imputation(NASH_20pct)
mice_25pct_NASH <- mice_imputation(NASH_25pct)
mice_30pct_NASH <- mice_imputation(NASH_30pct)
mice_40pct_NASH <- mice_imputation(NASH_40pct)
#combined 
mice_datasets <- list(
  NC_10 = mice_10pct_NC,
  NC_15 = mice_15pct_NC,
  NC_20 = mice_20pct_NC,
  NC_25 = mice_25pct_NC,
  NC_30 = mice_30pct_NC,
  NC_40 = mice_40pct_NC,
  
  HO_10 = mice_10pct_HO,
  HO_15 = mice_15pct_HO,
  HO_20 = mice_20pct_HO,
  HO_25 = mice_25pct_HO,
  HO_30 = mice_30pct_HO,
  HO_40 = mice_40pct_HO,
  
  NAFL_10 = mice_10pct_NAFL,
  NAFL_15 = mice_15pct_NAFL,
  NAFL_20 = mice_20pct_NAFL,
  NAFL_25 = mice_25pct_NAFL,
  NAFL_30 = mice_30pct_NAFL,
  NAFL_40 = mice_40pct_NAFL,
  
  NASH_10 = mice_10pct_NASH,
  NASH_15 = mice_15pct_NASH,
  NASH_20 = mice_20pct_NASH,
  NASH_25 = mice_25pct_NASH,
  NASH_30 = mice_30pct_NASH,
  NASH_40 = mice_40pct_NASH
)


# ------------------------------------
# TITLE: Statistical Tests  
# ------------------------------------

# --------------------------
# Part 1: Shapiro-Wilk test
# --------------------------

#check normality with shapiro wilk test
shapiro_test <- function(data, label = "") {
  data_copy <- data 
  
  #numeric
  numeric_data <- data_copy[, 3:ncol(data_copy)]
  
  #shapiro test
  shapiro_results <- apply(numeric_data, 2, function(x) shapiro.test(x)$p.value)
  adjusted_p_values <- p.adjust(shapiro_results, method = "BH")
  
  shapiro_df <- data.frame(
    Lipid = names(shapiro_results),
    p_value = shapiro_results,
    Adjusted_p_values = adjusted_p_values
  )
  
  significant_count <- sum(shapiro_df$Adjusted_p_values < 0.05, na.rm = TRUE)
  
  cat("\n---------------------------------------\n")
  cat("Dataset: ", label, "\n")
  cat("Number of significant/non-normal Lipids: ", significant_count, "\n")
  cat("---------------------------------------\n")
  
  #return the count for plotting
  return(list(result = shapiro_df, non_normal_count = significant_count))
}


#call function for each data 
#original 
shapiro_orig_NC <- shapiro_test(data_NC, label = "Original") #34
#Halfmin
shapiro_res_halfmin <- lapply(names(halfmin_datasets), function(name) {
  res <- shapiro_test(halfmin_datasets[[name]], label = name)
  data.frame(Method = "Half-min", Dataset = name, Non_Normal_Count = res$non_normal_count)
})

#KNN
shapiro_res_KNN <- lapply(names(KNN_datasets), function(name) {
  res <- shapiro_test(KNN_datasets[[name]], label = name)
  data.frame(Method = "KNN", Dataset = name, Non_Normal_Count = res$non_normal_count)
})

#RF
shapiro_res_RF <- lapply(names(RF_datasets), function(name) {
  res <- shapiro_test(RF_datasets[[name]], label = name)
  data.frame(Method = "RF", Dataset = name, Non_Normal_Count = res$non_normal_count)
})

#QRILC
shapiro_res_QRILC <- lapply(names(QRILC_datasets), function(name) {
  res <- shapiro_test(QRILC_datasets[[name]], label = name)
  data.frame(Method = "QRILC", Dataset = name, Non_Normal_Count = res$non_normal_count)
})

#MICE
shapiro_res_mice <- lapply(names(mice_datasets), function(name) {
  res <- shapiro_test(mice_datasets[[name]], label = name)
  data.frame(Method = "MICE", Dataset = name, Non_Normal_Count = res$non_normal_count)
})

#combine all results of shapiro test
shapiro_summary <- bind_rows(
  shapiro_res_halfmin,
  shapiro_res_KNN,
  shapiro_res_RF,
  shapiro_res_QRILC,
  shapiro_res_mice
)

#extract number of missingness and condition
shapiro_summary <- shapiro_summary %>%
  mutate(
    Missingness = as.numeric(gsub(".*_(\\d+)$", "\\1", Dataset)),
    Condition = gsub("_(\\d+)$", "", Dataset)  # Extract 'NC', 'HO', etc.
  )

pdf("/Users/marcinebessire/Desktop/Master_Thesis/NAFLD/Shapiro.pdf", width = 14, height = 10)

#now plot 
ggplot(shapiro_summary, aes(x = factor(Missingness), y = Non_Normal_Count, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 34, linetype = "dashed", color = "darkgreen", size = 1) +
  facet_wrap(~ Method) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Non-Normal Lipids by Condition and Missingness (Faceted by Imputation Method)",
    x = "Missingness Percentage (%)",
    y = "Number of Non-Normal Lipids",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

dev.off()


# --------------------------- 
# Part 2: T-test and Wilcoxon
# ---------------------------

stats_test <- function(original_data, imputed_data, condition_label) {
  lipid_cols <- colnames(original_data)[3:ncol(original_data)]
  
  #run tests and store raw p-values
  raw_results <- lapply(lipid_cols, function(lipid) {
    orig_vals <- original_data[[lipid]]
    imp_vals  <- imputed_data[[lipid]]
    
    t_pval <- tryCatch(t.test(orig_vals, imp_vals)$p.value, error = function(e) NA)
    w_pval <- tryCatch(wilcox.test(orig_vals, imp_vals)$p.value, error = function(e) NA)
    
    data.frame(
      Lipid = lipid,
      TTest_p = t_pval,
      Wilcox_p = w_pval
    )
  })
  
  df <- do.call(rbind, raw_results)
  
  #apply BH adjustment
  df$TTest_p_adj   <- p.adjust(df$TTest_p, method = "BH")
  df$Wilcox_p_adj  <- p.adjust(df$Wilcox_p, method = "BH")
  df$Condition     <- condition_label
  
  return(df)
}

#combine original
original_datasets <- list(
  NC_10 = data_NC,
  NC_15 = data_NC,
  NC_20 = data_NC,
  NC_25 = data_NC,
  NC_30 = data_NC,
  NC_40 = data_NC,
  
  HO_10 = data_HO,
  HO_15 = data_HO,
  HO_20 = data_HO,
  HO_25 = data_HO,
  HO_30 = data_HO,
  HO_40 = data_HO,
  
  NAFL_10 = data_NAFL,
  NAFL_15 = data_NAFL,
  NAFL_20 = data_NAFL,
  NAFL_25 = data_NAFL,
  NAFL_30 = data_NAFL,
  NAFL_40 = data_NAFL,
  
  NASH_10 = data_NASH,
  NASH_15 = data_NASH,
  NASH_20 = data_NASH,
  NASH_25 = data_NASH,
  NASH_30 = data_NASH,
  NASH_40 = data_NASH
)

#call function for halfmin
#halfmin
results_halfmin <- lapply(names(halfmin_datasets), function(name) {
  original <- original_datasets[[name]]
  imputed <- halfmin_datasets[[name]]
  condition <- gsub("_\\d+$", "", name)
  missingness <- as.numeric(gsub(".*_(\\d+)$", "\\1", name))
  
  df <- stats_test(original, imputed, condition)
  df$Dataset <- name
  df$Missingness <- missingness
  df$Method <- "Half-min"
  return(df)
})

#save dataframe
results_halfmin_df <- do.call(rbind, results_halfmin)

#count singificant p-values
summary_tests_halfmin <- results_halfmin_df %>%
  mutate(
    TTest_signif = TTest_p_adj < 0.05,
    Wilcox_signif = Wilcox_p_adj < 0.05
  ) %>%
  group_by(Method, Condition, Missingness) %>%
  summarise(
    TTest_Significant = sum(TTest_signif, na.rm = TRUE),
    Wilcox_Significant = sum(Wilcox_signif, na.rm = TRUE),
    .groups = "drop"
  )

#knn
results_knn <- lapply(names(KNN_datasets), function(name) {
  original <- original_datasets[[name]]
  imputed <- KNN_datasets[[name]]
  condition <- gsub("_\\d+$", "", name)
  missingness <- as.numeric(gsub(".*_(\\d+)$", "\\1", name))
  
  df <- stats_test(original, imputed, condition)
  df$Dataset <- name
  df$Missingness <- missingness
  df$Method <- "KNN"
  return(df)
})

#save dataframe
results_KNN_df <- do.call(rbind, results_knn)

#knn
summary_tests_hknn <- results_KNN_df %>%
  mutate(
    TTest_signif = TTest_p_adj < 0.05,
    Wilcox_signif = Wilcox_p_adj < 0.05
  ) %>%
  group_by(Method, Condition, Missingness) %>%
  summarise(
    TTest_Significant = sum(TTest_signif, na.rm = TRUE),
    Wilcox_Significant = sum(Wilcox_signif, na.rm = TRUE),
    .groups = "drop"
  )


#RF
results_rf <- lapply(names(RF_datasets), function(name) {
  original <- original_datasets[[name]]
  imputed <- RF_datasets[[name]]
  condition <- gsub("_\\d+$", "", name)
  missingness <- as.numeric(gsub(".*_(\\d+)$", "\\1", name))
  
  df <- stats_test(original, imputed, condition)
  df$Dataset <- name
  df$Missingness <- missingness
  df$Method <- "RF"
  return(df)
})

#save dataframe
results_RF_df <- do.call(rbind, results_rf)

#count sign. p-values
summary_tests_RF <- results_RF_df %>%
  mutate(
    TTest_signif = TTest_p_adj < 0.05,
    Wilcox_signif = Wilcox_p_adj < 0.05
  ) %>%
  group_by(Method, Condition, Missingness) %>%
  summarise(
    TTest_Significant = sum(TTest_signif, na.rm = TRUE),
    Wilcox_Significant = sum(Wilcox_signif, na.rm = TRUE),
    .groups = "drop"
  )

#QRILC
results_qrilc <- lapply(names(QRILC_datasets), function(name) {
  original <- original_datasets[[name]]
  imputed <- QRILC_datasets[[name]]
  condition <- gsub("_\\d+$", "", name)
  missingness <- as.numeric(gsub(".*_(\\d+)$", "\\1", name))
  
  df <- stats_test(original, imputed, condition)
  df$Dataset <- name
  df$Missingness <- missingness
  df$Method <- "QRILC"
  return(df)
})

#save dataframe
results_QRILC_df <- do.call(rbind, results_qrilc)

#count sign. p-values
summary_tests_QRILC <- results_QRILC_df %>%
  mutate(
    TTest_signif = TTest_p_adj < 0.05,
    Wilcox_signif = Wilcox_p_adj < 0.05
  ) %>%
  group_by(Method, Condition, Missingness) %>%
  summarise(
    TTest_Significant = sum(TTest_signif, na.rm = TRUE),
    Wilcox_Significant = sum(Wilcox_signif, na.rm = TRUE),
    .groups = "drop"
  )

#MICE
results_mice <- lapply(names(mice_datasets), function(name) {
  original <- original_datasets[[name]]
  imputed <- mice_datasets[[name]]
  condition <- gsub("_\\d+$", "", name)
  missingness <- as.numeric(gsub(".*_(\\d+)$", "\\1", name))
  
  df <- stats_test(original, imputed, condition)
  df$Dataset <- name
  df$Missingness <- missingness
  df$Method <- "MICE"
  return(df)
})

#save dataframe
results_mice_df <- do.call(rbind, results_mice)

#count sign. p-values
summary_tests_mice <- results_mice_df %>%
  mutate(
    TTest_signif = TTest_p_adj < 0.05,
    Wilcox_signif = Wilcox_p_adj < 0.05
  ) %>%
  group_by(Method, Condition, Missingness) %>%
  summarise(
    TTest_Significant = sum(TTest_signif, na.rm = TRUE),
    Wilcox_Significant = sum(Wilcox_signif, na.rm = TRUE),
    .groups = "drop"
  )


