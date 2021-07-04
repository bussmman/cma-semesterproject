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
library(gridExtra)
#library(raster)
#library(zoo)          # to smoothen using moving window functions
#library(SimilarityMeasures)   # Similarity measures

## Load the wild boar data and metadata
boar = wildschwein_BE
boar_meta = wildschwein_metadata

## Load the background map
map_BE = terra::rast("pk100_BE.tif")
# map_BE = raster("pk100_BE.tif")
# map_BE <- as(map_BE, "SpatialPixelsDataFrame")
# map_BE <- as.data.frame(map_BE)

## Create functions that are needed later on
# Relate the tracked location to the mating season
loc = function(dat,ref){
  l = length(dat$DatetimeUTC)
  name = deparse(substitute(dat))
  animal = ref[ref$TierName==name,]
  
  for (i in 1:l){
    if (dat$DatetimeUTC[i] %in% animal$DatetimeUTC_rounded){
      dat$E[i] = animal[animal$DatetimeUTC_rounded==dat$DatetimeUTC[i],5]
      dat$N[i] = animal[animal$DatetimeUTC_rounded==dat$DatetimeUTC[i],6]
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

# Subset a dataset for the time of a meeting
traj = function(male,female,dat,threshold){
  name = deparse(substitute(female))
  
  t = male$DatetimeUTC[which(male[,name] <= threshold)]
  dat[dat$DatetimeUTC_rounded %in% t,]
}

# Define Trajectory ID for each segment
traj_ID = function(dat){
  l = length(dat$DatetimeUTC_rounded)
  ID = 1
  dat$timelag = lead(dat$DatetimeUTC_rounded)
  dat$timediff = dat$timelag - dat$DatetimeUTC_rounded
  
  dat = dat %>%
    mutate(Traj_ID = NA)
  
  for (i in 1:l){
    td = as.numeric(dat$timediff[i], units="secs")
    td[is.na(td)] = 0
    
    if (td != 900){
      dat$Traj_ID[i] = ID
      ID = ID + 1
    } else {
      dat$Traj_ID[i] = ID
    }
  }
  return(dat)
}

# Omit trajectories with a length smaller than a set limit
omit = function(dat,lim){
  dat %>%
    group_by(Traj_ID) %>%
    mutate(n_pos = n()) %>%
    filter(n_pos > lim)
}

# Plot the trajectories of the boars during a potential
# mating event
traj_plot = function(dat,male){
  d0 = as.character(min(dat$DatetimeUTC_rounded))
  d1 = as.character(max(dat$DatetimeUTC_rounded))
  
  p1 = ggplot(dat,aes(x=E, y=N, group=interaction(TierName,Traj_ID))) +
    ggtitle(paste(d0,"-",d1,sep=" ")) +
    geom_path(aes(color=TierName),show.legend = FALSE) +
    geom_point(aes(color=TierName),show.legend = FALSE) +
    xlim(min(dat$E[dat$TierName==male])-500,
         max(dat$E[dat$TierName==male])+500) +
    ylim(min(dat$N[dat$TierName==male])-500,
         max(dat$N[dat$TierName==male])+500)
  
  p2 = ggplot(dat,aes(x=E, y=N, group=interaction(TierName,Traj_ID))) +
    ggtitle(paste(d0,"-",d1,sep=" ")) +
    geom_path(aes(color=TierName)) +
    geom_point(aes(color=TierName)) +
    xlim(min(dat$E[dat$TierName==male])-50,
         max(dat$E[dat$TierName==male])+50) +
    ylim(min(dat$N[dat$TierName==male])-50,
         max(dat$N[dat$TierName==male])+50)
  
  grid.arrange(p1, p2, ncol=2)
}

################################
## Step 2: Data preprocessing ##
################################

## Define the mating season
starttime = as_datetime("2015-11-01 00:00:00 UTC")
endtime = as_datetime("2016-04-30 23:45:00 UTC")

## Select the data for the mating season
# Round the times to 15 minutes to be able to match them
boar = boar %>%
  mutate(DatetimeUTC_rounded
         = lubridate::round_date(boar$DatetimeUTC,"15 minutes"))

# # Add a date only column to the data frame "boar"
# boar = boar %>%
#   mutate(Date = as.Date(boar$DatetimeUTC))

# select the data based on the dates of the mating season
season = boar[boar$DatetimeUTC_rounded >= starttime
              & boar$DatetimeUTC_rounded <= endtime,]

# Calculate the following position and timestamp for each
# data point
season = season %>%
  group_by(CollarID) %>%
  mutate(timelag=lead(DatetimeUTC_rounded),
         E_lag=lead(E),N_lag=lead(N),
         timediff=timelag-DatetimeUTC_rounded)
season = ungroup(season)

# Interpolate the E und N positions to the rounded Datetime
season$E = season$E + (
  as.numeric(season$DatetimeUTC-season$DatetimeUTC_rounded)
  /as.numeric(season$timediff)*(season$E_lag-season$E))

season$N = season$N + (
  as.numeric(season$DatetimeUTC-season$DatetimeUTC_rounded)
  /as.numeric(season$timediff)*(season$N_lag-season$N))

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

# Plot the duration of tracking for each wild boar during
# the mating season
season %>%
  ggplot(aes(DatetimeUTC, TierName, colour = TierName)) +
  geom_point(show.legend = FALSE)

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
## the females and the male for each male
# Amos
min(Amos$Venus, na.rm=TRUE)
min(Amos$Gaby, na.rm=TRUE)
min(Amos$Miriam, na.rm=TRUE)      # very likely meeting
min(Amos$Caroline, na.rm=TRUE)    # very likely meeting
min(Amos$Evelin, na.rm=TRUE) 
min(Amos$Frida, na.rm=TRUE)

# Ueli
min(Ueli$Venus, na.rm=TRUE)
min(Ueli$Gaby, na.rm=TRUE)        # no non-NA value!
min(Ueli$Miriam, na.rm=TRUE)      # very likely meeting
min(Ueli$Caroline, na.rm=TRUE)    # very likely meeting
min(Ueli$Evelin, na.rm=TRUE) 
min(Ueli$Frida, na.rm=TRUE)       # very likely meeting

## Examine all the "very likely meetings": Define a euclidian
## distance threshold (50m) and examine the trajectories
## of all wild boars at the times, when the "meeting couple"
## was closer to each other than the set threshold

# Amos & Miriam (Threshold = 50m)
AM50 = traj(Amos,Miriam,season,50)
AM50 = traj_ID(AM50)
AM50 = omit(AM50,1)
traj_plot(AM50,"Amos")

# Meeting 1: 2015-11-15 18:15:00 - 2015-11-15 18:30:00

# # Threshold = 100m
# AM100 = traj(Amos,Miriam,season,100)
# AM100 = traj_ID(AM100)
# AM100 = omit(AM100,3)
# traj_plot(AM100,"Amos")

# Amos & Caroline (Threshold = 50m)
AC50 = traj(Amos,Caroline,season,50)
AC50 = traj_ID(AC50)
AC50 = omit(AC50,1)
traj_plot(AC50,"Amos")

# Meeting 2: 2015-11-14 08:45:00 - 2015-11-14 15:00:00
# Meeting 3: 2015-11-14 15:30:00 - 2015-11-14 17:15:00
# There is only one point missing between meeting 2 and 3
# --> counts as one meeting (Meeting 2)

# # Threshold = 100m
# AC100 = traj(Amos,Caroline,season,100)
# AC100 = traj_ID(AC100)
# AC100 = omit(AC100,1)
# traj_plot(AC100,"Amos")

# Ueli & Miriam (Threshold = 50m)
UM50 = traj(Ueli,Miriam,season,50)
UM50 = traj_ID(UM50)
UM50 = omit(UM50,1)
traj_plot(UM50,"Ueli")

# Meeting 3: 2016-02-03 00:00:00 - 2016-02-03 00:30:00

# # Threshold = 100m
# UM100 = traj(Ueli,Miriam,season,100)
# UM100 = traj_ID(UM100)
# UM100 = omit(UM100,1)
# traj_plot(UM100,"Ueli")

# Ueli & Caroline (Threshold = 50m)
UC50 = traj(Ueli,Caroline,season,50)
UC50 = traj_ID(UC50)
UC50 = omit(UC50,1)
traj_plot(UC50,"Ueli")

# Meeting 4: 2016-03-24 20:45:00 - 2016-03-24 21:00:00
# Meeting 5: 2016-03-24 21:30:00 - 2016-03-24 21:45:00
# There is only one point missing between meeting 4 and 5
# --> counts as one meeting (Meeting 4)
# Meeting 5: 2016-04-08 20:45:00 - 2016-04-08 21:30:00

# # Threshold = 100m
# UC100 = traj(Ueli,Caroline,season,100)
# UC100 = traj_ID(UC100)
# UC100 = omit(UC100,1)
# traj_plot(UC100,"Ueli")

# Ueli & Frida (Threshold = 50m)
UF50 = traj(Ueli,Frida,season,50)
UF50 = traj_ID(UF50)
UF50 = omit(UF50,1)
traj_plot(UF50,"Ueli")

# Meeting 6: 2016-02-01 06:45:00 - 2016-02-01 15:15:00
# Meeting 7: 2016-02-02 23:45:00 - 2016-02-03 00:15:00
# --> Corresponds to Meeting 3 (according to the timestamp)
# --> Refered to as Meeting 3!
# Meeting 7: 2016-04-11 00:45:00 - 2016-04-11 01:00:00
# Meeting 8: 2016-04-11 23:30:00 - 2016-04-12 00:45:00
# Meeting 9: 2016-04-12 20:30:00 - 2016-04-12 23:15:00
#Meeting 10: 2016-04-13 07:45:00 - 2016-04-13 09:00:00
#Meeting 11: 2016-04-28 22:30:00 - 2016-04-28 23:00:00

# # Threshold = 100m
# UF100 = traj(Ueli,Frida,season,100)
# UF100 = traj_ID(UF100)
# UF100 = omit(UF100,1)
# traj_plot(UF100,"Ueli")

## Plot selected Meeting events of interest
# Meeting 2: Amos and Caroline
M2 = AC50
traj_plot(M2,"Amos")

# Meeting 3: Ueli, Miriam and Frida
M3 = season[season$DatetimeUTC_rounded
            %in% seq(as_datetime("2016-02-02 23:45:00 UTC"),
                     as_datetime("2016-02-03 00:30:00 UTC"),
                     by="15 min"),]
M3 = traj_ID(M3)
traj_plot(M3,"Ueli")

# Meeting 4: Ueli and Caroline
M4 = season[season$DatetimeUTC_rounded
             %in% seq(as_datetime("2016-03-24 20:45:00 UTC"),
                      as_datetime("2016-03-24 21:45:00 UTC"),
                      by="15 min"),]
M4 = traj_ID(M4)
traj_plot(M4,"Ueli")

# Meeting 5: Ueli and Caroline
M5 = omit(UC50,max(UC50$n_pos)-1)
traj_plot(M5,"Ueli")

# Meeting 6: Ueli and Frida
M6 = omit(UF50,max(UF50$n_pos)-1)
traj_plot(M6,"Ueli")

# Meeting 10: Ueli and Frida
M10 = season[season$DatetimeUTC_rounded
            %in% seq(as_datetime("2016-04-13 07:45:00 UTC"),
                     as_datetime("2016-04-13 09:00:00 UTC"),
                     by="15 min"),]
M10 = traj_ID(M10)
traj_plot(M10,"Ueli")

## Add Environmental Context to the Meeting Events
# Turn Meeting Events into sf-objects
# Meeting 2
M2_sf = M2[M2$TierName %in% c("Amos","Caroline"),]
M2_sf = st_as_sf(M2_sf, coords = c("E", "N"), crs = 2056)
M2_sf$Event = "M2"
M2_sf$n_pos = NULL

# Meeting 3
M3_sf = M3[M3$TierName %in% c("Ueli","Miriam","Frida"),]
M3_sf = st_as_sf(M3_sf, coords = c("E", "N"), crs = 2056)
M3_sf$Event = "M3"

# Meeting 4
M4_sf = M4[M4$TierName %in% c("Ueli","Caroline"),]
M4_sf = st_as_sf(M4_sf, coords = c("E", "N"), crs = 2056)
M4_sf$Event = "M4"

# Meeting 5
M5_sf = M5[M5$TierName %in% c("Ueli","Caroline"),]
M5_sf = st_as_sf(M5_sf, coords = c("E", "N"), crs = 2056)
M5_sf$Event = "M5"
M5_sf$n_pos = NULL

# Meeting 6
M6_sf = M6[M6$TierName %in% c("Ueli","Frida"),]
M6_sf = st_as_sf(M6_sf, coords = c("E", "N"), crs = 2056)
M6_sf$Event = "M6"
M6_sf$n_pos = NULL

# Meeting 10
M10_sf = M10[M10$TierName %in% c("Ueli","Frida"),]
M10_sf = st_as_sf(M10_sf, coords = c("E", "N"), crs = 2056)
M10_sf$Event = "M10"

# Merge Meeting Events together in one data frame
M_tot = rbind(M2_sf,M3_sf,M4_sf,M5_sf,M6_sf,M10_sf)
M_tot = summarise(group_by(M_tot,Event))
M_cvx_h = st_convex_hull(M_tot)

# Plot the convex hull of the meeting events (only including
# the locations of the wild boars that actually met)
#tmap_mode("view")
tm_shape(map_BE) + 
  tm_rgb() +
  tm_shape(M_cvx_h) +
  tm_polygons(col = "Event",alpha = 0.4,border.col = "red") +
  tm_legend(bg.color = "white")
