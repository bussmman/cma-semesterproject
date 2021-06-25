############################################################ 
## Detecting potential mating through trajectory analysis ##
## of fertile male and female wild boars                  ##
############################################################

####################################################
## Step 1: Loading the required data and packages ##
####################################################

## Load packages
library(ComputationalMovementAnalysisData)
                      # to load the wild boar data
library(ggplot2)      # to visualize data
library(readr)        # to import tabular data (e.g. csv)
library(dplyr)        # to manipulate (tabular) data
library(sf)           # to handle spatial vector data
library(terra)        # To handle raster data
library(lubridate)    # To handle dates and times
#library(zoo)          # To smoothen using moving window functions
#library(SimilarityMeasures)   # Similarity measures

## Load the wild boar data and metadata
boar = wildschwein_BE
boar_meta = wildschwein_metadata

## Load the background map
map_BE = terra::rast("pk100_BE.tif")

################################
## Step 2: Data preprocessing ##
################################

## Select the data for the mating season 2015/16
# Add a date only column to the data frame "boar"
boar = boar %>%
  mutate(Date = as.Date(boar$DatetimeUTC))

# select the data based on the dates of the mating season
season = boar[boar$Date >= "2015-11-01"
              & boar$Date <= "2016-04-30",]

## Identify the fertile males and females that were
## tracked during the mating season
# Get the names, sex and weight of all wild boars
boar_weight = boar_meta[boar_meta$TierName %in% season$TierName,
              c("TierName","Sex","Gewicht")]

# Order them by weight
boar_weight = boar_weight[order(boar_weight$Gewicht),]

# Divide them up according to sex
males = boar_weight[boar_weight$Sex=="m",]
females = boar_weight[boar_weight$Sex=="f",]
unique(males$TierName)                          # 2 males
unique(females$TierName)                        # 6 females
