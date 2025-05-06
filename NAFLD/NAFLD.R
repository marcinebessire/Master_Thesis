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

#remove specie column
data_whole <- data_whole[,-4]

#set data to numeric 
data_whole <- data_whole %>%
  mutate(across(4:ncol(.), as.numeric))

#also set column 2 (number of NAs) as numeric
data_whole[[2]] <- as.numeric(as.character(data_whole[[2]]))

str(data_whole)

# ----------------------------------------
# Part 2: Remove Rows/Patient that have NA
# ----------------------------------------

data_whole <- data_whole %>%
  drop_na()



