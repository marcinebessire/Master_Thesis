#load necessary libraries 
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(openxlsx)

# --------------------------------------
# Title: Whole dataset 
# --------------------------------------

#load data
FAO_data <- read.csv("/Users/marcinebessire/Desktop/Master_Thesis/FAO_data.csv", check.names = FALSE) #34 metabolites

# --------------------------------------
# Part 1: Plot kinetic distribution
# --------------------------------------

#convert visit and time to factor and numeric
FAO_data <- FAO_data %>%
  mutate(
    Visit = as.factor(Visit),
    Time_min = as.numeric(Time_min)
  )

#reshape into long format
FAO_long <- FAO_data %>%
  pivot_longer(cols = 6:ncol(FAO_data), #metabolite columns only
               names_to = "Metabolite",
               values_to = "Concentration")

#loop through each patient and metabolite
patients <- unique(FAO_long$Patient)

#open a PDF device
pdf("/Users/marcinebessire/Desktop/Master_Thesis/No_Separation/FAO_original_Kinetics.pdf", width = 14, height = 10)


for (patient in patients) {

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
}

#close pdf device
dev.off()


# --------------------------------------
# Part 2: MCAR Simulation
# --------------------------------------

#function to introduce MCAR per dataframe randomly
MCAR_manipulation <- function(data, missing_percentage){
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

# --------------------------------------
# Part 3: MNAR Simulation
# --------------------------------------

#function to introduce MNAR by choosing the lowest values
MNAR_manipulation <- function(data, missing_percentage){
  #copy dataset to avoid modifying the original
  data_copy <- data
  
  for (col in colnames(data_copy[6:ncol(data_copy)])) {
    #sort rows by metbaolite value (ascending)
    ordered_indices <- order(data_copy[[col]], na.last = NA)
    
    #calculate how many mv depending on percentage
    num_mv <- round(nrow(data_copy) * missing_percentage)
    
    #randomly sample rows in middle
    selected_rows <- ordered_indices[1:num_mv]
    
    #set to NA
    data_copy[selected_rows, col] <- NA
    
    
  }
  
  return(data_copy)
}

#call function
FAO_10pct_mnar <- MNAR_manipulation(FAO_data, 0.1)
FAO_20pct_mnar <- MNAR_manipulation(FAO_data, 0.2)
FAO_25pct_mnar <- MNAR_manipulation(FAO_data, 0.25)
FAO_30pct_mnar <- MNAR_manipulation(FAO_data, 0.3)
FAO_35pct_mnar <- MNAR_manipulation(FAO_data, 0.35)
FAO_40pct_mnar <- MNAR_manipulation(FAO_data, 0.4)


