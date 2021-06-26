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
library(terra)        # to handle raster data
library(lubridate)    # to handle dates and times
library(tmap)         # to create thematic maps
#library(zoo)          # to smoothen using moving window functions
#library(SimilarityMeasures)   # Similarity measures

## Load the wild boar data and metadata
boar = wildschwein_BE
boar_meta = wildschwein_metadata

## Load the background map
map_BE = terra::rast("pk100_BE.tif")

## Create functions that are needed later on
# Relate the tracked location to the mating season
loc = function(dat,ref){
  l = length(dat$DatetimeUTC)
  name = deparse(substitute(dat))
  animal = ref[ref$TierName==name,]
  
  for (i in 1:l){
    if (dat$DatetimeUTC[i] %in% animal$DatetimeUTC){
      dat$E[i] = animal[animal$DatetimeUTC==dat$DatetimeUTC[i],"E"]
      dat$N[i] = animal[animal$DatetimeUTC==dat$DatetimeUTC[i],"N"]
    } else {
      dat$E[i] = NA
      dat$N[i] = NA
    }
  }
  return(dat)
}

################################
## Step 2: Data preprocessing ##
################################

## Define the mating season
starttime = as_datetime("2015-11-01 00:00:00 UTC")
endtime = as_datetime("2016-04-30 23:45:00 UTC")

## Select the data for the mating season
# Round the times to 15 minutes to be able to match them
boar$DatetimeUTC = lubridate::round_date(
  boar$DatetimeUTC,"15 minutes")

# # Add a date only column to the data frame "boar"
# boar = boar %>%
#   mutate(Date = as.Date(boar$DatetimeUTC))

# select the data based on the dates of the mating season
season = boar[boar$DatetimeUTC >= starttime
              & boar$DatetimeUTC <= endtime,]

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

## Identify the tracked locations during mating season
## for each male and female wild boar
# Males
Amos = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Amos$E = NA
Amos$N = NA
Amos = loc(Amos,season)

Ueli = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Ueli$E = NA
Ueli$N = NA
Ueli = loc(Ueli,season)

# Females
Venus = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Venus$E = NA
Venus$N = NA
Venus = loc(Venus,season)

Gaby = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Gaby$E = NA
Gaby$N = NA
Gaby = loc(Gaby,season)

Miriam = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Miriam$E = NA
Miriam$N = NA
Miriam = loc(Miriam,season)

Caroline = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Caroline$E = NA
Caroline$N = NA
Caroline = loc(Caroline,season)

Evelin = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Evelin$E = NA
Evelin$N = NA
Evelin = loc(Evelin,season)

Frida = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Frida$E = NA
Frida$N = NA
Frida = loc(Frida,season)

#################################
## Step 3: Modelling procedure ##
#################################

