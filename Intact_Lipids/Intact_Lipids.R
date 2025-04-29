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
data_10pct <- MNAR_simulation(data, 0.1) #10% missing values
data_20pct <- MNAR_simulation(data, 0.2) #20% missing values
data_30pct <- MNAR_simulation(data, 0.3) #30% missing values
data_40pct <- MNAR_simulation(data, 0.4) #40% missing values

# --------------------------
# Part 2: Plot MNAR
# --------------------------

#plot missing values
plot_missing_values <- function(data, title) {
  #convert data to long format
  data_long <- melt(data, id.vars = c("ID", "Patient", "Date", "Time_min", "Visit"))
  
  #create a column indicating missing values
  data_long$Missing <- ifelse(is.na(data_long$value), "Missing", "Present")
  
  #correctly order Patinet IDs 
  sorted_patient <- mixedsort(unique(data_long$Patient))  
  data_long$Patient <- factor(data_long$Patient, levels = sorted_patient)  #apply order
  
  #plot missing values using a heatmap
  ggplot(data_long, aes(x = variable, y = Patient, fill = Missing)) +
    geom_tile(color = "grey") +  #add borders for clarity
    facet_wrap(~Visit, ncol = 2) +  #separate by Visit 1 and Visit 2
    scale_fill_manual(values = c("Present" = "white", "Missing" = "red")) +
    theme_minimal(base_size = 16) +  #set global font size
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),  
      axis.text.y = element_text(size = 12, face = "bold"),  
      axis.title.x = element_text(size = 16, face = "bold"),  
      axis.title.y = element_text(size = 16, face = "bold"),  
      strip.text = element_text(size = 18, face = "bold"),  
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
      legend.title = element_text(size = 16, face = "bold"),  
      legend.text = element_text(size = 14)  
    ) +
    labs(title = title, x = "Lipid", y = "Patient ID", fill = "Data Status")
}

#call function for plotting
plot_missing_values(data_10pct, "Missing Data Pattern (10% MNAR)")
plot_missing_values(data_20pct, "Missing Data Pattern (20% MNAR)")
plot_missing_values(data_30pct, "Missing Data Pattern (30% MNAR)")
plot_missing_values(data_40pct, "Missing Data Pattern (40% MNAR)")

# ---------------------------------
# TITLE: Prepare Data for Imputation
# ----------------------------------

# ------------------------------------
# Part 1: Separate data by Visit
# ------------------------------------

#original 
data_original_v1 <- data %>% filter(Visit == "Visit 1")
data_original_v2 <- data %>% filter(Visit == "Visit 2")
#simulated
data_10pct_v1 <- data_10pct %>% filter(Visit == "Visit 1")
data_10pct_v2 <- data_10pct %>% filter(Visit == "Visit 2")
data_20pct_v1 <- data_20pct %>% filter(Visit == "Visit 1")
data_20pct_v2 <- data_20pct %>% filter(Visit == "Visit 2")
data_30pct_v1 <- data_30pct %>% filter(Visit == "Visit 1")
data_30pct_v2 <- data_30pct %>% filter(Visit == "Visit 2")
data_40pct_v1 <- data_40pct %>% filter(Visit == "Visit 1")
data_40pct_v2 <- data_40pct %>% filter(Visit == "Visit 2")


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
pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/Lipid_Comparison_Boxplots.pdf", width = 14, height = 10)

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

#call function for Half-min imputation
#v1
halfmin_10pct_v1 <- half_min_imputation(data_10pct_v1)
halfmin_20pct_v1 <- half_min_imputation(data_20pct_v1)
halfmin_30pct_v1 <- half_min_imputation(data_30pct_v1)
halfmin_40pct_v1 <- half_min_imputation(data_40pct_v1)
#v2
halfmin_10pct_v2 <- half_min_imputation(data_10pct_v2)
halfmin_20pct_v2 <- half_min_imputation(data_20pct_v2)
halfmin_30pct_v2 <- half_min_imputation(data_30pct_v2)
halfmin_40pct_v2 <- half_min_imputation(data_40pct_v2)

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
#Visit 1
knn_10pct_v1 <- knn_imputation(data_10pct_v1)
knn_20pct_v1 <- knn_imputation(data_20pct_v1)
knn_30pct_v1 <- knn_imputation(data_30pct_v1)
knn_40pct_v1 <- knn_imputation(data_40pct_v1)
#Visit 2
knn_10pct_v2 <- knn_imputation(data_10pct_v2)
knn_20pct_v2 <- knn_imputation(data_20pct_v2)
knn_30pct_v2 <- knn_imputation(data_30pct_v2)
knn_40pct_v2 <- knn_imputation(data_40pct_v2)

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
#Visit 1
rf_10pct_v1 <- rf_imputation(data_10pct_v1)
rf_20pct_v1 <- rf_imputation(data_20pct_v1)
rf_30pct_v1 <- rf_imputation(data_30pct_v1)
rf_40pct_v1 <- rf_imputation(data_40pct_v1)
#Visit 2
rf_10pct_v2 <- rf_imputation(data_10pct_v2)
rf_20pct_v2 <- rf_imputation(data_20pct_v2)
rf_30pct_v2 <- rf_imputation(data_30pct_v2)
rf_40pct_v2 <- rf_imputation(data_40pct_v2)


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
#Visit 1
qrilc_10pct_v1 <- qrilc_imputation(data_10pct_v1)
qrilc_20pct_v1 <- qrilc_imputation(data_20pct_v1)
qrilc_30pct_v1 <- qrilc_imputation(data_30pct_v1)
qrilc_40pct_v1 <- qrilc_imputation(data_40pct_v1)
#Visit 2
qrilc_10pct_v2 <- qrilc_imputation(data_10pct_v2)
qrilc_20pct_v2 <- qrilc_imputation(data_20pct_v2)
qrilc_30pct_v2 <- qrilc_imputation(data_30pct_v2)
qrilc_40pct_v2 <- qrilc_imputation(data_40pct_v2)

# ---------------------------------
# TITLE: Combine datasets into one
# ----------------------------------

#make funciton to merge and order dataframes for each Imputation
merge_and_order <- function(data1, data2) {
  #combine datasets
  merged_data <- rbind(data1, data2)
  
  #convert Participant to a properly sorted factor
  merged_data$Patient <- factor(merged_data$Patient,
                                    levels = mixedsort(unique(merged_data$Patient)))
 
  #order by patient and date
  merged_data <- merged_data[order(merged_data$Patient, merged_data$Date), ]
  
  
  #reset row numbering to ensure sequential indices
  row.names(merged_data) <- NULL  
  
  return(merged_data)
}

#call function
#Halfmin
halfmin_10pct_tot <- merge_and_order(halfmin_10pct_v1, halfmin_10pct_v2)
halfmin_20pct_tot <- merge_and_order(halfmin_20pct_v1, halfmin_20pct_v2)
halfmin_30pct_tot <- merge_and_order(halfmin_30pct_v1, halfmin_30pct_v2)
halfmin_40pct_tot <- merge_and_order(halfmin_40pct_v1, halfmin_40pct_v2)
#KNN
knn_10pct_tot <- merge_and_order(knn_10pct_v1, knn_10pct_v2)
knn_20pct_tot <- merge_and_order(knn_20pct_v1, knn_20pct_v2)
knn_30pct_tot <- merge_and_order(knn_30pct_v1, knn_30pct_v2)
knn_40pct_tot <- merge_and_order(knn_40pct_v1, knn_40pct_v2)
#RF
rf_10pct_tot <- merge_and_order(rf_10pct_v1, rf_10pct_v2)
rf_20pct_tot <- merge_and_order(rf_20pct_v1, rf_20pct_v2)
rf_30pct_tot <- merge_and_order(rf_30pct_v1, rf_30pct_v2)
rf_40pct_tot <- merge_and_order(rf_40pct_v1, rf_40pct_v2)
#QRILC
qrilc_10pct_tot <- merge_and_order(qrilc_10pct_v1, qrilc_10pct_v2)
qrilc_20pct_tot <- merge_and_order(qrilc_20pct_v1, qrilc_20pct_v2)
qrilc_30pct_tot <- merge_and_order(qrilc_30pct_v1, qrilc_30pct_v2)
qrilc_40pct_tot <- merge_and_order(qrilc_40pct_v1, qrilc_40pct_v2)

# ------------------------------------
# TITLE: Statistical Tests  
# ------------------------------------

# --------------------------
# Part 1: Shapiro-Wilk test
# --------------------------

#check normality with shapiro wilk test
shapiro_test <- function(data){
  data_copy <- data 
  
  #numeric
  numeric_data <- data_copy[, 6:ncol(data)]
  
  #get shapiro restuls 
  shapiro_results <- apply(numeric_data, 2, function(x) shapiro.test(x)$p.value)
  
  #adjust p value
  adjusted_p_values <- p.adjust(shapiro_results, method = "BH")
  
  #convert to dataframe
  shapiro_df <- data.frame(
    Lipid = names(shapiro_results),
    p_value = shapiro_results,
    Adjusted_p_values = adjusted_p_values
  )
  
  #print how many are sign (if p-value < 0.05 then not normal distribution)
  significant_count <- sum(shapiro_df$Adjusted_p_values < 0.05, na.rm = TRUE)
  cat("\n---------------------------------------\n")
  cat("Number of singificant/non-normal Lipids: ", significant_count, "\n")
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
#visit 1
shapiro1_halfmin10 <- shapiro_test(halfmin_10pct_v1) #5/274
shapiro1_halfmin20 <- shapiro_test(halfmin_20pct_v1) #2/274
shapiro1_halfmin30 <- shapiro_test(halfmin_30pct_v1) #2/274
shapiro1_halfmin40 <- shapiro_test(halfmin_40pct_v1) #5/274
#visit 2
shapiro2_halfmin10 <- shapiro_test(halfmin_10pct_v2) #3/274
shapiro2_halfmin20 <- shapiro_test(halfmin_20pct_v2) #3/274
shapiro2_halfmin30 <- shapiro_test(halfmin_30pct_v2) #6/274
shapiro2_halfmin40 <- shapiro_test(halfmin_40pct_v2) #6/274

#knn
#visit 1
shapiro1_knn10 <- shapiro_test(knn_10pct_v1) #96/274
shapiro1_knn20 <- shapiro_test(knn_20pct_v1) #109/274
shapiro1_knn30 <- shapiro_test(knn_30pct_v1) #162/274
shapiro1_knn40 <- shapiro_test(knn_40pct_v1) #179/274
#visit 2
shapiro2_knn10 <- shapiro_test(knn_10pct_v2) #83/274
shapiro2_knn20 <- shapiro_test(knn_20pct_v2) #103/274
shapiro2_knn30 <- shapiro_test(knn_30pct_v2) #108/274
shapiro2_knn40 <- shapiro_test(knn_40pct_v2) #122/274

#RF
#visit 1
shaprio1_rf10 <- shapiro_test(rf_10pct_v1) #61/274
shaprio1_rf20 <- shapiro_test(rf_20pct_v1) #72/274
shaprio1_rf30 <- shapiro_test(rf_30pct_v1) #106/274
shaprio1_rf40 <- shapiro_test(rf_40pct_v1) #119/274
#visit 2
shaprio2_rf10 <- shapiro_test(rf_10pct_v2) #43/274
shaprio2_rf20 <- shapiro_test(rf_20pct_v2) #78/274
shaprio2_rf30 <- shapiro_test(rf_30pct_v2) #71/274
shaprio2_rf40 <- shapiro_test(rf_40pct_v2) #67/274

#QRILC
#visit1
shapiro1_qrilc10 <- shapiro_test(qrilc_10pct_v1) #7/274
shapiro1_qrilc20 <- shapiro_test(qrilc_20pct_v1) #38/174
shapiro1_qrilc30 <- shapiro_test(qrilc_30pct_v1) #19/274
shapiro1_qrilc40 <- shapiro_test(qrilc_40pct_v1) #63/274
#visit2
shapiro2_qrilc10 <- shapiro_test(qrilc_10pct_v2) #3/274
shapiro2_qrilc20 <- shapiro_test(qrilc_20pct_v2) #3/174
shapiro2_qrilc30 <- shapiro_test(qrilc_30pct_v2) #3/274
shapiro2_qrilc40 <- shapiro_test(qrilc_40pct_v2) #0/274

pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/Shapiro_Wilk.pdf", width = 14, height = 10)

#combine all data visit 1
shapiro_summary1 <- data.frame(
  Method = c(
    rep("Original", 1),
    rep("Half-min", 4),
    rep("KNN", 4),
    rep("RF", 4),
    rep("QRILC", 4)
  ),
  Missingness = c(
    0, #original
    
    10, 20, 30, 40,  #halfmin
    10, 20, 30, 40,  #knn
    10, 20, 30, 40,  #RF
    10, 20, 30, 40   #QRILC
  ),
  Non_Normal_Count = c(
    7, #original
    
    5, 2, 2, 5,  #Halfmin
    96, 109, 162, 179,     #knn
    61, 72, 106, 119,     #RF
    7, 38, 19, 63     #QRILC
  )
)

#make ggplot (bar plot)
ggplot(shapiro_summary1, aes(x = factor(Missingness), y = Non_Normal_Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Non-Normally Distributed Lipids (Visit 1)",
       x = "Missingness Percentage (%)",
       y = "Count of Non-Normal Lipids",
       fill = "Imputation Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#extract the Non_Normal_Count for the "Original" method
original_non_normal1 <- shapiro_summary1 %>%
  filter(Method == "Original") %>%
  pull(Non_Normal_Count)

#plot (dot plot with line)
ggplot(shapiro_summary1, aes(x = Missingness, y = Non_Normal_Count, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = original_non_normal1, linetype = "dashed", color = "lightgreen", size = 1) + 
  theme_minimal(base_size = 16) +  
  labs(
    title = "Effect of Imputation on Normality of Lipids (Visit 1)",
    x = "Missingness Percentage (%)",
    y = "Count of Non-Normal Lipids",
    color = "Imputation Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  
    axis.text.y = element_text(size = 14, face = "bold"),  
    axis.title.x = element_text(size = 16, face = "bold"), 
    axis.title.y = element_text(size = 16, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
    legend.title = element_text(size = 16, face = "bold"),  
    legend.text = element_text(size = 14) 
  )

#combine all data visit 2
shapiro_summary2 <- data.frame(
  Method = c(
    rep("Original", 1),
    rep("Half-min", 4),
    rep("KNN", 4),
    rep("RF", 4),
    rep("QRILC", 4)
  ),
  Missingness = c(
    0, #original
    
    10, 20, 30, 40,  #halfmin
    10, 20, 30, 40,  #knn
    10, 20, 30, 40,  #RF
    10, 20, 30, 40   #QRILC
  ),
  Non_Normal_Count = c(
    3, #original
    
    3, 3, 6, 6,  #Halfmin
    83, 103, 108, 122,     #knn
    43, 78, 71, 67,     #RF
    3, 3, 3, 0     #QRILC
  )
)


#make ggplot (bar plot)
ggplot(shapiro_summary2, aes(x = factor(Missingness), y = Non_Normal_Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Number of Non-Normally Distributed Lipids (Visit 2)",
       x = "Missingness Percentage (%)",
       y = "Count of Non-Normal Lipids",
       fill = "Imputation Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#extract the Non_Normal_Count for the "Original" method
original_non_normal2 <- shapiro_summary2 %>%
  filter(Method == "Original") %>%
  pull(Non_Normal_Count)

#plot (dot plot with line)
ggplot(shapiro_summary2, aes(x = Missingness, y = Non_Normal_Count, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = original_non_normal2, linetype = "dashed", color = "lightgreen", size = 1) + 
  theme_minimal(base_size = 16) +  
  labs(
    title = "Effect of Imputation on Normality of Lipids (Visit 2)",
    x = "Missingness Percentage (%)",
    y = "Count of Non-Normal Lipids",
    color = "Imputation Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  
    axis.text.y = element_text(size = 14, face = "bold"),  
    axis.title.x = element_text(size = 16, face = "bold"), 
    axis.title.y = element_text(size = 16, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
    legend.title = element_text(size = 16, face = "bold"),  
    legend.text = element_text(size = 14) 
  )

dev.off()

# --------------------------
# Part 2: T-test (Visit 1 vs Visit 2)
# --------------------------

#function to make t-test between original and visit (should not be significantly different if good data)
visit_statistical_tests <- function(data){
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
  cat("Number of significant metabolites:\n")
  cat("- Wilcoxon test:", wilcox_sign, "\n")
  cat("- Paired t-test:", t_test_sign, "\n")
  cat("-------------------------------\n")
  
  
  #return results dataframe
  return(results)
}


#call function to compare Visit 1 vs. Visit 2
#original 
stats_origina <- visit_statistical_tests(data)  #0 and 0
#Halfmin
visit_Halfmin10pct_res <- visit_statistical_tests(halfmin_10pct_tot) #0 and 0
visit_Halfmin20pct_res <- visit_statistical_tests(halfmin_20pct_tot) #0 and 0
visit_Halfmin30pct_res <- visit_statistical_tests(halfmin_30pct_tot) #0 and 0
visit_Halfmin40pct_res <- visit_statistical_tests(halfmin_40pct_tot) #0 and 0
#KNN
visit_KNN10pct_res <- visit_statistical_tests(knn_10pct_tot) #0 W and 5 T
visit_KNN20pct_res <- visit_statistical_tests(knn_20pct_tot) #0 W and 5 T
visit_KNN30pct_res <- visit_statistical_tests(knn_30pct_tot) #0 W and 14 T
visit_KNN40pct_res <- visit_statistical_tests(knn_40pct_tot) #0 W and 5 T
#RF
visit_RF10pct_res <- visit_statistical_tests(rf_10pct_tot) #0 W and 6 T
visit_RF20pct_res <- visit_statistical_tests(rf_20pct_tot) #0 W and 27 T
visit_RF30pct_res <- visit_statistical_tests(rf_30pct_tot) #0 W and 50 T
visit_RF40pct_res <- visit_statistical_tests(rf_40pct_tot) #50 W and 67 T
#QRILC
visit_QRILC10pct_res <- visit_statistical_tests(qrilc_10pct_tot) #0 and 0
visit_QRILC20pct_res <- visit_statistical_tests(qrilc_20pct_tot) #0 and 0
visit_QRILC30pct_res <- visit_statistical_tests(qrilc_30pct_tot) #0 and 14
visit_QRILC40pct_res <- visit_statistical_tests(qrilc_40pct_tot) #0 and 19

#create data frame for number of significant metabolites
stats_overall <- data.frame(
  Method = rep(c("Halfmin", "Halfmin", "Halfmin", "Halfmin", "Halfmin",
                 "KNN", "KNN", "KNN", "KNN", "KNN",
                 "RF", "RF", "RF", "RF", "RF",
                 "QRILC", "QRILC", "QRILC", "QRILC", "QRILC"), each = 1),
  Percentage = c("Original", "10%", "20%", "30%", "40%",
                 "Original", "10%", "20%", "30%", "40%",
                 "Original", "10%", "20%", "30%", "40%",
                 "Original", "10%", "20%", "30%", "40%"),
  Wilcoxon = c(0, 0, 0, 0, 0,   
               0, 0, 0, 0, 0,
               0, 0, 0, 0, 50,  
               0, 0, 0, 0, 0),  
  TTest = c(0, 0, 0, 0, 0,  
            0, 5, 5, 14, 5,  
            0, 6, 27, 50, 67,  
            0, 0, 0, 14, 19)  
)

#convert percentage to factor 
stats_overall$Percentage <- factor(stats_overall$Percentage, levels = c("Original", "10%", "20%", "30%", "40%"))


#reshape into long format
stats_long <- stats_overall %>%
  pivot_longer(cols = c("Wilcoxon", "TTest"), names_to = "Test", values_to = "Significant_Lipids")


pdf("/Users/marcinebessire/Desktop/Master_Thesis/Intact_Lipids/T_Test_Wilcox.pdf", width = 14, height = 10)

#ceate grouped bar plot
ggplot(stats_long, aes(x = Percentage, y = Significant_Lipids, fill = Test)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Significant_Lipids), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 6, fontface = "bold") +  
  facet_wrap(~ Method, scales = "free_x", nrow = 2) +  
  labs(title = "Number of Significant Lipids (Visit 1 vs Visit 2)",
       x = "Missingness Percentage",
       y = "Number of Significant Lipids",
       fill = "Statistical Test") +
  theme_minimal(base_size = 16) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  
    axis.text.y = element_text(size = 14, face = "bold"),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.title.y = element_text(size = 16, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
    legend.title = element_text(size = 16, face = "bold"),  
    legend.text = element_text(size = 14), 
    strip.text = element_text(size = 16, face = "bold"), 
    panel.spacing = unit(2, "lines")  
  )

dev.off()














