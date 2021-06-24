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

################################
## Step 2: Data preprocessing ##
################################

## Identify the fertile males and females
# Get the names, sex and weight of all wild boars in BE
boar_weight = boar_meta[boar_meta$TierName %in% boar$TierName,
              c("TierName","Sex","Gewicht")]

# Order them by weight
boar_weight = boar_weight[order(boar_weight$Gewicht),]

# Divide them up according to sex
males = boar_weight[boar_weight$Sex=="m",]      # 6 males
females = boar_weight[boar_weight$Sex=="f",]    # 13 females
