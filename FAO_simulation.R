#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

#load data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/FAO_data.csv", check.names = FALSE) #34 metabolites

# --------------------------------------
# Part 1: Plot kinetic distribution
# --------------------------------------

#reshape into long format
FAO_long <- FAO_data %>%
  pivot_longer(cols=6:ncol(.),
               names_to = "Metabolite",
               values_to = "Concentration")

#loop through each patient and metabolite
unique_patients <- unique(FAO_long$Patient)

for (patient in unique_patients) {
    
    plot_data <- FAO_long %>%
      filter(Patient == patient)
    
    p <- ggplot(plot_data, aes(x = Time_min, y = Concentration, color = Visit)) +
      geom_point(size = 2) +
      geom_line(aes(group = Visit), linewidth = 1) +
      facet_wrap(~ Metabolite, scales = "free_y", ncol = 6) +
      labs(
        title = paste(patient, ":", "Kinetics of All Metabolites"),
        subtitle = "Visit 1 vs visit 2",
        x = "Time [min]",
        y = "Concentration [µM]"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom", 
        plot.title = element_text(face = "bold", size = 14)
      ) +
      scale_color_manual(values = c("Visit 1" = "lightblue", "Visit 2" = "orange"))
    
    #show plot
    plot(p)
    
    #press enter to go to the next plot 
    readline(prompt = "Press [Enter] to see the next plot")
}


# -------------------------------------------------
# Part 2: Split Dataframe according to Visit only
# -------------------------------------------------

# -------------------------------------------------
# Part 2.1: MCAR
# -------------------------------------------------

#function to generate missing values per visit dataframe
MCAR_manipulation_visit <- function(data, missing_percentage){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #get the Visit 1 and Visit 2 indices
  visit1_indices <- which(data_copy$Visit == "Visit 1")
  visit2_indices <- which(data_copy$Visit == "Visit 2")
  
  #iterate through each numeric column
  for (col in 6:ncol(data_copy)) {
    #get number of missing values to introduce (rounded down to ensure pairing)
    num_mv_total <- round(nrow(data_copy) * missing_percentage)
    num_mv <- floor(num_mv_total / 2)  #even split between Visit 1 and Visit 2
    
    if (num_mv > 0 && length(visit1_indices) > 0 && length(visit2_indices) > 0) {
      #randomly choose indices for missing values in each visit group
      missing_indices_v1 <- sample(visit1_indices, num_mv, replace = FALSE)
      missing_indices_v2 <- sample(visit2_indices, num_mv, replace = FALSE)
      
      #set selected values to NA
      data_copy[missing_indices_v1, col] <- NA
      data_copy[missing_indices_v2, col] <- NA
    }
  }
  
  return(data_copy)
}

#call function
FAO_visit_mcar <- MCAR_manipulation_visit(FAO_data, 0.2)

#split into two dataframe visit 1 and visit 2
FAO_visit1_mcar <- FAO_visit_mcar %>% filter(Visit == "Visit 1") #60
FAO_visit2_mcar <- FAO_visit_mcar %>% filter(Visit == "Visit 2") #59

# -------------------------------------------------
# Part 2.1: MNAR
# -------------------------------------------------

#function to generate missing values per visit dataframe
MNAR_manipulation_visit <- function(data, missing_percentage){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #get the Visit 1 and Visit 2 indices
  visit1_indices <- which(data_copy$Visit == "Visit 1")
  visit2_indices <- which(data_copy$Visit == "Visit 2")
  
  #iterate through each numeric column
  for (col in 6:ncol(data_copy)) {
    #get number of missing values to introduce (rounded down to ensure pairing)
    num_mv_total <- round(nrow(data_copy) * missing_percentage)
    num_mv <- floor(num_mv_total / 2)  #even split between Visit 1 and Visit 2
    
    if (num_mv > 0 && length(visit1_indices) >= num_mv && length(visit2_indices) >= num_mv) {
      #column values for each visit
      col_v1 <- data_copy[visit1_indices, col]
      col_v2 <- data_copy[visit2_indices, col]
      
      #indices with lowest value between each group
      lowest_v1_indices <- visit1_indices[order(col_v1, na.last = NA)][1:num_mv]
      lowest_v2_indices <- visit2_indices[order(col_v2, na.last = NA)][1:num_mv]
      
      #set to NA
      data_copy[lowest_v1_indices, col] <- NA
      data_copy[lowest_v2_indices, col] <- NA
     }
  }
  
  return(data_copy)
}

#call function
FAO_visit_mnar <- MNAR_manipulation_visit(FAO_data, 0.2)

#split into two dataframe visit 1 and visit 2
FAO_visit1_mnar <- FAO_visit_mnar %>% filter(Visit == "Visit 1") #60
FAO_visit2_mnar <- FAO_visit_mnar %>% filter(Visit == "Visit 2") #59

# --------------------------------------
# Part 3: Patient and Visit separated
# --------------------------------------

# ---------------------------------------------------
# Part 3.1: Split Dataframe according to Patient and Visit
# ----------------------------------------------------

p1_visit1 <- FAO_data[1:6,]
p1_visit2 <- FAO_data[7:12,]
p2_visit1 <- FAO_data[13:18,]
p2_visit2 <- FAO_data[19:24,]
p3_visit1 <- FAO_data[25:30,]
p3_visit2 <- FAO_data[31:36,]
p4_visit1 <- FAO_data[37:41,] #misses 120 min
p4_visit2 <- FAO_data[42:47,]
p5_visit1 <- FAO_data[48:53,]
p5_visit2 <- FAO_data[54:59,]
p6_visit1 <- FAO_data[60:65,]
p6_visit2 <- FAO_data[66:71,]
p7_visit1 <- FAO_data[72:77,]
p7_visit2 <- FAO_data[78:83,]
p8_visit1 <- FAO_data[84:89,]
p8_visit2 <- FAO_data[90:95,]
p9_visit1 <- FAO_data[96:101,]
p9_visit2 <- FAO_data[102:107,]
p10_visit1 <- FAO_data[108:113,]
p10_visit2 <- FAO_data[114:119,]

# --------------------------------------
# Part 3.2: MCAR 
# --------------------------------------

#function to introduce 1 MCAR per dataframe randomly

#middle values only (one missing value)
MCAR_manipulation_middle <- function(data){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #filter eligible rows (30, 60 or 120)
  middle_rows <- which(data_copy$Time_min %in% c(30, 60, 120))
  
  for (col in colnames(data_copy[6:ncol(data_copy)])) {
    rand_row <- sample(middle_rows, 1)
    data_copy[rand_row, col] <- NA
  }
  
  return(data_copy)
}

#call function
#p1
p1_v1_mcar <- MCAR_manipulation_middle(p1_visit1)
p1_v2_mcar <- MCAR_manipulation_middle(p1_visit2)
#p2
p2_v1_mcar <- MCAR_manipulation_middle(p2_visit1)
p2_v2_mcar <- MCAR_manipulation_middle(p2_visit2)
#p3
p3_v1_mcar <- MCAR_manipulation_middle(p3_visit1)
p3_v2_mcar <- MCAR_manipulation_middle(p3_visit2)
#p4
p4_v1_mcar <- MCAR_manipulation_middle(p4_visit1)
p4_v2_mcar <- MCAR_manipulation_middle(p4_visit2)
#p5
p5_v1_mcar <- MCAR_manipulation_middle(p5_visit1)
p5_v2_mcar <- MCAR_manipulation_middle(p5_visit2)
#p6
p6_v1_mcar <- MCAR_manipulation_middle(p6_visit1)
p6_v2_mcar <- MCAR_manipulation_middle(p6_visit2)
#p7
p7_v1_mcar <- MCAR_manipulation_middle(p7_visit1)
p7_v2_mcar <- MCAR_manipulation_middle(p7_visit2)
#p8
p8_v1_mcar <- MCAR_manipulation_middle(p8_visit1)
p8_v2_mcar <- MCAR_manipulation_middle(p8_visit2)
#p9
p9_v1_mcar <- MCAR_manipulation_middle(p9_visit1)
p9_v2_mcar <- MCAR_manipulation_middle(p9_visit2)
#p10
p10_v1_mcar <- MCAR_manipulation_middle(p10_visit1)
p10_v2_mcar <- MCAR_manipulation_middle(p10_visit2)

# --------------------------------------
# Part 3.3: MNAR
# --------------------------------------

#function to introduce 1 MNAR per dataframe (one missing value)
MNAR_manipulation_lowest <- function(data){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  #go through each column
  for (col in colnames(data_copy[6:ncol(data_copy)])) {
    min_row <- which.min(data_copy[[col]])
    data_copy[min_row, col] <- NA
  }
  
  return(data_copy)
}

#call function
#p1
p1_v1_mnar <- MNAR_manipulation_lowest(p1_visit1)
p1_v2_mnar <- MNAR_manipulation_lowest(p1_visit2)
#p2
p2_v1_mnar <- MNAR_manipulation_lowest(p2_visit1)
p2_v2_mnar <- MNAR_manipulation_lowest(p2_visit2)
#p3
p3_v1_mnar <- MNAR_manipulation_lowest(p3_visit1)
p3_v2_mnar <- MNAR_manipulation_lowest(p3_visit2)
#p4
p4_v1_mnar <- MNAR_manipulation_lowest(p4_visit1)
p4_v2_mnar <- MNAR_manipulation_lowest(p4_visit2)
#p5
p5_v1_mnar <- MNAR_manipulation_lowest(p5_visit1)
p5_v2_mnar <- MNAR_manipulation_lowest(p5_visit2)
#p6
p6_v1_mnar <- MNAR_manipulation_lowest(p6_visit1)
p6_v2_mnar <- MNAR_manipulation_lowest(p6_visit2)
#p7
p7_v1_mnar <- MNAR_manipulation_lowest(p7_visit1)
p7_v2_mnar <- MNAR_manipulation_lowest(p7_visit2)
#p8
p8_v1_mnar <- MNAR_manipulation_lowest(p8_visit1)
p8_v2_mnar <- MNAR_manipulation_lowest(p8_visit2)
#p9
p9_v1_mnar <- MNAR_manipulation_lowest(p9_visit1)
p9_v2_mnar <- MNAR_manipulation_lowest(p9_visit2)
#p10
p10_v1_mcar <- MNAR_manipulation_lowest(p10_visit1)
p10_v2_mcar <- MNAR_manipulation_lowest(p10_visit2)


# --------------------------------------
# Part 4: Whole dataset
# --------------------------------------

# --------------------------------------
# Part 4.1: MCAR 
# --------------------------------------

#function to introduce 1 MCAR per dataframe randomly

#middle values only 
MCAR_manipulation_middle_whole <- function(data, missing_percentage){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  for (col in colnames(data_copy[6:ncol(data_copy)])) {
    #calculate how many mv depending on percentage
    num_mv <- round(nrow(data_copy) * missing_percentage)

    #randomly sample rows in middle
    selected_rows <- sample(1:nrow(data_copy), num_mv, replace = FALSE)
      
    #set to NA
    data_copy[selected_rows, col] <- NA
      
  
  }
  
  return(data_copy)
}

#call function
FAO_10pct_mcar <- MCAR_manipulation_middle_whole(FAO_data, 0.1)
FAO_20pct_mcar <- MCAR_manipulation_middle_whole(FAO_data, 0.2)
FAO_25pct_mcar <- MCAR_manipulation_middle_whole(FAO_data, 0.25)
FAO_30pct_mcar <- MCAR_manipulation_middle_whole(FAO_data, 0.3)
FAO_35pct_mcar <- MCAR_manipulation_middle_whole(FAO_data, 0.35)
FAO_40pct_mcar <- MCAR_manipulation_middle_whole(FAO_data, 0.4)


