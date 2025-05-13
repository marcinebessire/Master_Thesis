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
library(VIM) #for knn
library(caTools)  # for trapz function

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

#check if numeric conversion worked
str(data_whole)

# ------------------------------------
# TITLE: Imputation methods 
# ------------------------------------

# ------------------------------------
# Part 1: Half-min Imputation 
# ------------------------------------

#function for half-min imputation
halfmin_imputation <- function(data){
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

#call function to impute missing values
data_halfmin <- halfmin_imputation(data_whole)

# --------------------------
# Part 2: KNN Imputation
# --------------------------

#function for KNN imputation
knn_imputation <- function(data, k = 10) {
  data_copy <- data
  
  #metadata
  metadata <- data_copy[, 1:2]
  num_data <- data_copy[, 3:ncol(data_copy)]
  
  #combine numeric data and dummy ID (required by VIM kNN to preserve structure)
  num_data$ID_temp__ <- 1:nrow(num_data)
  
  #apply KNN imputation
  imputed_data <- kNN(num_data, variable = colnames(num_data)[1:(ncol(num_data)-1)], 
                      k = k, imp_var = FALSE)
  
  #rm temporary ID
  imputed_data$ID_temp__ <- NULL
  
  #combine back with metadata
  final_data <- cbind(metadata, imputed_data)
  
  return(final_data)
}

#call function for KNN
data_knn <- knn_imputation(data_whole)

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
data_rf <- rf_imputation(data_whole)

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

#call function to impute missing values
data_qrilc <- qrilc_imputation(data_whole)

# ------------------
# TITLE: PCA
# ------------------

#function for PCA
pca_func <- function(data, method){
  #numeric 
  numeric <- data[, 3:ncol(data)]
  
  #metadata
  metadata <- data[, 1:2]
  
  #perform PCA
  pca_res <- prcomp(numeric, center = TRUE, scale. = TRUE)

  #combine scores with metadata
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Condition <- metadata$Condition  #for coloring
  
  #scatter plot: Patients that are far apart have different lipid patterns
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA of Lipidomics Data", x = "PC1", y = "PC2")
  
  print(p1)
  
  #scree plot: x-axis = PC1, PC2 etc., y-axis = proportion of variance
  plot(pca_res, type = "lines", main = "Scree Plot")
  
}

#call pca function
pca_halfmin <- pca_func(data_halfmin)
pca_knn <- pca_func(data_knn)
pca_rf <- pca_func(data_rf)
pca_qrilc <- pca_func(data_qrilc)








