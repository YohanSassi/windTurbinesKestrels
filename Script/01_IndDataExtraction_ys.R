#___________________________________________________________________________________________________________________________________#
###
###       
###       Project: Relative probability of selecton by Lesser Kestrels in Aumelas wind farm
###
###        Data come from wild Lesser Kestrels equipped by LPO, deposited on MoveBank plateform.
### 
###       !!! Script used to extract GPS locations of Lesser kestrels, filter and save them !!!
###
###       Author: Yohan Sassi 
###               yohan.sassi@cefe.cnrs.fr
###
#___________________________________________________________________________________________________________________________________#



rm(list = ls())

#________________________________________________________________
#
#### Package needed      ####
#
#________________________________________________________________

library(dplyr)
library(lubridate) 
library(amt)
library(stringr)
library("move")
library(sf)
library(sp)
library(rgeos)
library(raster)



#Homemade functions & needed files
source("Functions/df_to_spatial_df_utm.R") #To convert long/lat to UTM



#________________________________________________________________
#
#### Data preparation, cleaning / filtering      ####
#
#________________________________________________________________


# Parameters
loginStored <- move::movebankLogin(username="XXX", password="XXX") # Store login requirement for MOVEBANK
studyName <- "Falco naumanni Aumelas wind farm interactions [ID_PROG 311]" # Study name
#year <- c("2017", "2018", "2019","2020") #years of study

indNames <- unique(move::getMovebankAnimals(study = studyName, login = loginStored)["local_identifier"] %>% 
  filter(!grepl("2016|2022|2023",.$local_identifier)) %>% 
  pull()) # Names of individual in the study for the year considered


# -- LOOP ON EACH INDIVIDUAL
for (i in 1:length(indNames)) {
  
  # Download data for an individual for a specific period - from movebank
  indData <- move::getMovebankLocationData(study = studyName,
                                           sensorID = "GPS",
                                           login = loginStored,
                                           animalName = indNames[i],
                                           timestamp_start = "20170101000000000",
                                           timestamp_end = "20211231235959000")
  
  # Remove outlier & extract only interesting columns (reduce the weight of the dataset)
  indData <- indData[which(indData$visible == "true"),c(1,3:6,11,16)]
  
  # Add column of date & take the list
  indData$timestamp <- strptime(indData$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  indData$date <- strptime(indData$timestamp, format = "%Y-%m-%d")
  
  # Keep data only between May and July, and between 07:00 and 18:00 UTC so between 09:00 and 20:00 local time
  indData <- indData %>% 
    filter(lubridate::month(date) >= 5 & lubridate::month(date) <= 7) %>% 
    filter(lubridate::hour(timestamp) >= 7 & lubridate::hour(timestamp) < 18) %>% 
    filter(is.na(location.lat) == FALSE & is.na(location.long) == FALSE) # Be sure that no NA in lat/long
  
  if (nrow(indData) == 0){ # If no data corresponding to the filtering move to the next individuals
    next
  } else {
    
    #Transform coordinates to UTM
    df_to_spatial_df_utm(dataframe=indData, 
                         column_longitude="location.long", 
                         column_latitude="location.lat", 
                         utm_zone=31, hemisphere="N")
    
    # Add the UTM coordinates to the initial df
    indData <- cbind(indData, dataframe.sp_utm@coords)
    
    # Keep data in a 10km buffer around the center of the MCP
    CRS.31N <- CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs") # CRS of the project (proj4string argument)
    baryCoordinates <- SpatialPoints(coords = cbind(550421,4820652), proj4string = CRS.31N) # UTM 31N Latitude & longitude of the MCP barycenter in an sp object
    baryMCPBuffer <- gBuffer(baryCoordinates,  width=10000, quadsegs = 100) # 10km buffer
    lon.lat.matrix <- matrix(data = c(indData$longitude,indData$latitude), nrow = length(indData$longitude), ncol = 2) # create lon lat matrix to create an sp.dataframe
    sp.df <- SpatialPointsDataFrame(lon.lat.matrix, data = indData, proj4string = CRS.31N) # create an sp dataframe for the sp::over function
    indData <- indData[!is.na(sp::over(sp.df, baryMCPBuffer)),]
     
    # Visual check
    # plot(SpatialPoints(cbind(indData[,c("longitude","latitude")]), proj4string = CRS.31N))
    # plot(baryMCPBuffer, col = "red", add = TRUE)
    # plot(SpatialPoints(cbind(indData2[,c("longitude","latitude")])), col = "green", add = TRUE)
    
    # Free memory
    rm(dataframe.sp,dataframe.sp_utm)
    gc()
    
    if (nrow(indData) == 0){ # Again, if no data into the buffer then move to next individual
      next
    } else {
      
      # Extract all dates this individual was present
      allDates <- unique(indData$date)
      
      # -- LOOP ON EACH DAY FOR THIS INDIVIDUAL
      for (d in 1:length(allDates)){ 
        
        # Data for a specific day
        wantindData <- indData[which(indData$date == allDates[d]),]
        
        # Convert into track_xyt for amt
        wantindData_Track <- make_track(wantindData,
                                        longitude,
                                        latitude,
                                        timestamp)
        
        ## Subsample to have a point every 10min or the closest
        # if no points found in 10mn then it take the following point and start again from this point
        wantindData_TrackResampled <- track_resample(wantindData_Track, 
                                                     rate = minutes(10), 
                                                     tolerance = minutes(0), 
                                                     start = 1)
        # Concatenate the re-sampled GPS fixes
        if (d == 1){
          allSubsampled <- wantindData_TrackResampled 
        } else {
          allSubsampled <- rbind(allSubsampled, wantindData_TrackResampled)
        }
      }
      
      # Keep only the points corresponding to the re-sampled paths
      indData <- indData[which(indData$longitude %in% allSubsampled$x_ &
                                 indData$latitude %in% allSubsampled$y_ &
                                 indData$timestamp %in% allSubsampled$t_),]
      
      # Save the individual file
      hyfenindices <- stringr::str_locate_all(pattern ='-', indNames[i]) # locate the - in the name
      cutIndName <- substr(indNames[i], 
                           as.integer(hyfenindices[[1]][4,1])-3, 
                           as.integer(hyfenindices[[1]][4,1])+7) # cut around the - to have inside the []
      
      # Save as R.Data
      # save(indData,
      #      file = paste("./Data/FilteredIndData/indData_allyears_",cutIndName,".Rdata", sep=""))

      # Free memory
      rm(allSubsampled)
      gc()
      
    }
  }
}
