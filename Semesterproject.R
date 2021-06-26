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
      dat$E[i] = animal[animal$DatetimeUTC==dat$DatetimeUTC[i],5]
      dat$N[i] = animal[animal$DatetimeUTC==dat$DatetimeUTC[i],6]
    } else {
      dat$E[i] = NA
      dat$N[i] = NA
    }
  }
  return(dat)
}

# Calculate the Euclidian distance between a male and a
# female wild boar at all times
euc_dis = function(male,female){
  male = male %>%
    mutate(n=sqrt((male$E-female$E)^2+(male$N-female$N)^2))
  
  names(male)[names(male)=="n"] = deparse(substitute(female))
  return(male)
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
Amos$E[lengths(Amos$E)!=1] = NA
Amos$E = as.double(unlist(Amos$E))
Amos$N[lengths(Amos$N)!=1] = NA
Amos$N = as.double(unlist(Amos$N))

Ueli = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Ueli$E = NA
Ueli$N = NA
Ueli = loc(Ueli,season)
Ueli$E[lengths(Ueli$E)!=1] = NA
Ueli$E = as.double(unlist(Ueli$E))
Ueli$N[lengths(Ueli$N)!=1] = NA
Ueli$N = as.double(unlist(Ueli$N))

# Females
Venus = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Venus$E = NA
Venus$N = NA
Venus = loc(Venus,season)
Venus$E[lengths(Venus$E)!=1] = NA
Venus$E = as.double(unlist(Venus$E))
Venus$N[lengths(Venus$N)!=1] = NA
Venus$N = as.double(unlist(Venus$N))

Gaby = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Gaby$E = NA
Gaby$N = NA
Gaby = loc(Gaby,season)
Gaby$E[lengths(Gaby$E)!=1] = NA
Gaby$E = as.double(unlist(Gaby$E))
Gaby$N[lengths(Gaby$N)!=1] = NA
Gaby$N = as.double(unlist(Gaby$N))

Miriam = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Miriam$E = NA
Miriam$N = NA
Miriam = loc(Miriam,season)
Miriam$E[lengths(Miriam$E)!=1] = NA
Miriam$E = as.double(unlist(Miriam$E))
Miriam$N[lengths(Miriam$N)!=1] = NA
Miriam$N = as.double(unlist(Miriam$N))

Caroline = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Caroline$E = NA
Caroline$N = NA
Caroline = loc(Caroline,season)
Caroline$E[lengths(Caroline$E)!=1] = NA
Caroline$E = as.double(unlist(Caroline$E))
Caroline$N[lengths(Caroline$N)!=1] = NA
Caroline$N = as.double(unlist(Caroline$N))

Evelin = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Evelin$E = NA
Evelin$N = NA
Evelin = loc(Evelin,season)
Evelin$E[lengths(Evelin$E)!=1] = NA
Evelin$E = as.double(unlist(Evelin$E))
Evelin$N[lengths(Evelin$N)!=1] = NA
Evelin$N = as.double(unlist(Evelin$N))

Frida = data.frame(DatetimeUTC=seq(starttime,endtime,900))
Frida$E = NA
Frida$N = NA
Frida = loc(Frida,season)
Frida$E[lengths(Frida$E)!=1] = NA
Frida$E = as.double(unlist(Frida$E))
Frida$N[lengths(Frida$N)!=1] = NA
Frida$N = as.double(unlist(Frida$N))

#################################
## Step 3: Modelling procedure ##
#################################

## Calculate the Euclidian distance between all the females
## and the male at all times for each male
# Amos
Amos = euc_dis(Amos,Venus)
Amos = euc_dis(Amos,Gaby)
Amos = euc_dis(Amos,Miriam)
Amos = euc_dis(Amos,Caroline)
Amos = euc_dis(Amos,Evelin)
Amos = euc_dis(Amos,Frida)

# Ueli
Ueli = euc_dis(Ueli,Venus)
Ueli = euc_dis(Ueli,Gaby)
Ueli = euc_dis(Ueli,Miriam)
Ueli = euc_dis(Ueli,Caroline)
Ueli = euc_dis(Ueli,Evelin)
Ueli = euc_dis(Ueli,Frida)

## Check for the minimal distance that occured between all
## the females and the male at all times for each male
# Amos
min(Amos$Venus, na.rm=TRUE)
min(Amos$Gaby, na.rm=TRUE)
min(Amos$Miriam, na.rm=TRUE)      # very likely meeting
min(Amos$Caroline, na.rm=TRUE)    # very likely meeting
min(Amos$Evelin, na.rm=TRUE) 
min(Amos$Frida, na.rm=TRUE)       # likely meeting

# Ueli
min(Ueli$Venus, na.rm=TRUE)
min(Ueli$Gaby, na.rm=TRUE)        # no non-NA value!
min(Ueli$Miriam, na.rm=TRUE)      # very likely meeting
min(Ueli$Caroline, na.rm=TRUE)    # very likely meeting
min(Ueli$Evelin, na.rm=TRUE) 
min(Ueli$Frida, na.rm=TRUE)       # very likely meeting
