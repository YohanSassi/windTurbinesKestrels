#___________________________________________________________________________________________________________________________________#
###
###       
###       Project: Relative probability of selecton by Lesser Kestrels in Aumelas wind farm
###
###        Data come from wild Lesser Kestrels equipped by LPO, deposited on MoveBank plateform.
### 
###       !!! Script associate environmental variables to individual tracking data
###           + sample pseudo-false GPS locations in MCP to fit RSF !!!
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

library("oce")
library(lubridate)
library(dplyr)
library(doParallel)
library(plyr)
library(tictoc)
library(amt)
library(stringr)
library(sf)
library(sp)
library(rgeos)
library(raster)


#Homemade functions & needed files
source("Functions/df_to_spatial_df_utm.R") #To convert long/lat to UTM
source("Functions/WindRose_fct.R") #To plot windrose



#________________________________________________________________
#
#### Parameters & required files      ####
#
#________________________________________________________________

# List of individual data
indData_fls <- list.files(path = "./Data/FilteredIndData", full.names = TRUE, recursive = FALSE)

# Meteorological data
meteoData <- read.csv("./Data/MeteoFrance_StAndre_Sangonis_2017-2021 mai-juil_6-19hUTC.csv", header = TRUE) # Meteo file

# Identify days with problems of wind measurement
weatherStationStopDays <- meteoData %>% 
  dplyr::select(., -c(outdoor_temp_avg)) %>% # Remove useless columns - full NA
  dplyr::mutate(timestamp = strptime(as.POSIXct(.$Date) + seconds(0), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")) %>% # Convert string to datetime in GMT
  filter(lubridate::month(timestamp) >= 5 & lubridate::month(timestamp) <= 7) %>% # filter out years, months and hours not considered in the study
  filter(lubridate::hour(timestamp) >= 7 & lubridate::hour(timestamp) <= 19) %>% 
  filter(lubridate::year(timestamp) >= 2017 & lubridate::year(timestamp) <= 2021) %>% 
  filter_all(any_vars(is.na(.))) # Keep only days where NA in one of the columns

weatherStationStopDays <- unique(lubridate::date(weatherStationStopDays$Date))
  
# Filter the meteoData
meteoData <- meteoData %>%
  dplyr::select(., -c(outdoor_temp_avg)) %>% # Remove useless columns - full NA
  dplyr::mutate(timestamp = strptime(as.POSIXct(.$Date) + seconds(0), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")) %>% # Convert string to datetime in GMT
  filter(lubridate::month(timestamp) >= 5 & lubridate::month(timestamp) <= 7) %>% # filter out years, months and hours not considered in the study
  filter(lubridate::hour(timestamp) >= 7 & lubridate::hour(timestamp) <= 19) %>% 
  filter(lubridate::year(timestamp) >= 2017 & lubridate::year(timestamp) <= 2021) %>% 
  filter_all(~ !is.na(.))

## Exploration WindRose
plot.windrose(data = meteoData,
              spd = meteoData$wind_speed_avg,
              dir = meteoData$wind_direction_avg,
              spdmax = 16,
              spdmin = 0,
              spdseq = seq(0,16,2),
              dirres = 45) + 
  theme_bw() +
  xlab("") + 
  ylab("Frequency") +
  theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 14),
              panel.grid.minor = element_line(colour = "grey93"),
              panel.grid.major = element_line(colour = "grey93"),
              strip.background = element_rect(colour = "white",
                                              fill = "white"),
              strip.text = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 14))





# Raster of slope, exposition and landcover (projection: UTM31N)
expoRaster <- raster("./Data/expo_31N.tif") # Exposition
slopeRaster <- raster("./Data/pente_31N.tif") # Slope
landcoverRaster <- raster("./Data/landcover.tif") # Landcover

allRaster <- stack(expoRaster, slopeRaster, landcoverRaster) # Stack raster together 

# Table of colony / bird ID association
colonyTable <-  read.csv("./Data/colony.csv", sep = ";", header = TRUE) 
colonyTable$colony <- str_replace_all(colonyTable$colony, " ", "_") # Remove spaces in names of colony
colonyTable$logger <- gsub(pattern= "LAR15 (PIC26)", replacement= "LAR15", x = colonyTable$logger, fixed = T) # change LR15 logger id

# Projection UTM31N stored in object
CRS.31N <- CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs") # CRS of the project (proj4string argument)

# Load shapefile to have the MCP limits
MCPshape <- read_sf('./Data/MCP_300buffer_31N.shp')
MCPpoly <- as(st_geometry(MCPshape), "Spatial")
MCPpoly@proj4string <- CRS.31N

# Load 200m buffered sleeping sites
SleepingSites <- read_sf("./Data/Dortoirs_et_nids_2016-2020_buff200_31N.shp")
SleepingSites <- as(st_geometry(SleepingSites), "Spatial")
SleepingSites@proj4string <- CRS.31N


# Load table with turbines locations + add latlong coordinates
windTurbineTable <-  read.csv("./Data/locWindTurbinesAumelas.csv", header = TRUE)




#________________________________________________________________
#
#### Association of environmental variables to individual data ####
#
#________________________________________________________________


# If you have multiple cores you can run this in parallel
detectCores()
doParallel::registerDoParallel(detectCores()-2)


# -- LOOP ON EACH INDIVIDUAL
tic()
for (i in 1:length(indData_fls)){
  
  
  # Load individual data
  load(indData_fls[i])
  
  # Filter data to remove those present in days where weather station was not working
  if(any(unique(lubridate::date(indData$date)) %in% weatherStationStopDays)){
    
    indData <- indData %>% 
      filter(!(lubridate::date(.$date) %in% weatherStationStopDays)) # keep data corresponding to days not in weatherStationStopDays
  }

  # Convert into a dataframe to then make a track - needed for sampling environmental variables
  forTrack <- indData %>% 
    dplyr::rename(x = longitude,
                  y = latitude,
                  t = timestamp,
                  ID = individual.id)
  
  # Convert in track 
  forTrack.xyt <- mk_track(forTrack,
                           x,
                           y,
                           all_cols = TRUE)
  forTrack.xyt$x_ <- as.numeric(forTrack.xyt$x_)
  forTrack.xyt$y_ <- as.numeric(forTrack.xyt$y_)
  
  
  # Add additional information:
  # Landscape features -> exposition, slope and landcover
  # meteorological data -> wind direction and speed
  # Orographic updraft velocity
  # Colony of origin
  # Loc in MCP or not
  IndData_allVars <- forTrack.xyt %>% 
    extract_covariates(allRaster) %>% # Extract simultaneously expo, slope and landcover
    dplyr::mutate(expo_rad = (.$expo_31N*pi) / 180, # convert expo in radian
                  slope_rad = (.$pente_31N*pi) / 180, # convert slope in radian
                  colony = colonyTable[str_which(forTrack.xyt$individual.local.identifier[1], colonyTable$logger),"colony"][[1]]) %>% # Add colony
    bind_cols(do.call("rbind",llply(1:nrow(.), function(a){ # Need to do each point as round_date can end up with same timestamp 
      
      # Is the loc in the MCP? 
      inMCP = ifelse(is.na(sp::over(SpatialPoints(forTrack.xyt[a,c("x_","y_")], proj4string = CRS.31N), MCPpoly)[[1]]), 0, 1)
      
      # Is the loc in a SleepingSite?
      inSleepingSite = ifelse(is.na(sp::over(SpatialPoints(forTrack.xyt[a,c("x_","y_")], proj4string = CRS.31N), SleepingSites)[[1]]), 0, 1)
      
      # Wind speed [m/s]
      windspeed <- meteoData[which(meteoData$timestamp == lubridate::round_date(.$t[a], "hour")),"wind_speed_avg"][[1]]
      
      # Wind direction [radian]
      winddir <- (meteoData[which(meteoData$timestamp == lubridate::round_date(.$t[a], "hour")),"wind_direction_avg"][[1]]*pi)/180
      
      # Estimation of orographic updraft velocity W0
      coef <- ifelse(.$slope_rad[a] != 0, sin(.$slope_rad[a])*cos(winddir - .$expo_rad[a]), 0)
      W0 <- windspeed * coef
      
      return(c(inMCP ,inSleepingSite, windspeed, winddir, W0))
    },.parallel = TRUE))) %>% 
    dplyr::rename(inMCP = ...17, # Rename columns
                  inSleepingSite = ...18,
                  windspeed = ...19, 
                  winddir_rad = ...20,
                  W0 = ...21) 
  
  
  # Extract name and save both df in separate files
  hyfenindices <- stringr::str_locate_all(pattern ='-', IndData_allVars$individual.local.identifier[1]) # locate the - in the name
  cutIndName <- substr(IndData_allVars$individual.local.identifier[1], 
                       as.integer(hyfenindices[[1]][4,1])-3, 
                       as.integer(hyfenindices[[1]][4,1])+7) # cut around the - to have inside the []
  
  
  if(i == 6){ IndData_allVars$colony <- "Saint-P_de_M"} # Correct colony for PIC27
  
  # Save individual dataset
  save(IndData_allVars, # This will be used for GLMM
       file = paste("./Data/IndData_VarAssociated/IndData_allVars_",cutIndName,".Rdata", sep=""))
  
}
toc()











#________________________________________________________________
#
#### Prepare data for RSFs     ####
#
#________________________________________________________________


# Quantile of Tramontane wind speed
meteoData %>% 
  filter(wind_direction_avg >= 270 & wind_direction_avg <= 360) %>% 
  reframe(Tram = c(quantile(wind_speed_avg,0.10), quantile(wind_speed_avg,0.5), quantile(wind_speed_avg,0.90)))


# Quantile of sea breeze wind speed
meteoData %>% 
  filter(wind_direction_avg >= 90 & wind_direction_avg <= 179) %>% 
  reframe(SeaBreeze = c(quantile(wind_speed_avg,0.10), quantile(wind_speed_avg,0.5), quantile(wind_speed_avg,0.90)))


# Object combining wind information for filtering
WindNames <- c("Tramontane", "SeaBreeze")

WindInfo <- vector("list")
WindInfo[[1]] <- list(name = WindNames[1],
                      angles = c(270,360),
                      speed = c(4.4, 7.5, 10.3))
WindInfo[[2]] <- list(name = WindNames[2],
                      angles = c(80,160),
                      speed = c(1.8, 4.2, 6.3))




#Data set will all individuals combined
allFiles <- list.files("./Data/IndData_VarAssociated" , all.files = F, full.names = T)

for (i in 1:length(allFiles)){
  
  load(allFiles[i])
  
  if (i == 1){
    allLoc <- IndData_allVars
  } else {
    allLoc <- rbind(allLoc, IndData_allVars)
  }
} 


# Filtering to have only data inMCP, from the considered colonies
# Convert in track
allLoc <- allLoc %>% 
  filter(.$colony %in% c("Saint-P_de_M", "Saint-Pargoire", "Villeveyrac")) %>% # Filter out individuals of other colony that have less than 5 loc in MCP
  filter(.$inMCP == 1) 


for(v in 1:2){ #Loop on main wind direction
  
  forRSF_subwind <- allLoc %>% 
    filter(winddir_rad >= (WindInfo[[v]]$angles[1]*pi)/180 & winddir_rad <= (WindInfo[[v]]$angles[2]*pi)/180)
  
  
  for(q in 1:3){ #Loop on quantile for wind speed, q1 = 10% (Low), q2 = median, q3 = 90% (High)
    
    #Filter to keep data with wind speed +/- 1m/s around the considered quantile
    forRSF <- forRSF_subwind %>% 
      rowwise() %>%
      filter(any(between(windspeed, WindInfo[[v]]$speed[q] - 1, WindInfo[[v]]$speed[q] + 1))) %>% 
      mk_track(., x_, y_, all_cols = TRUE)
    
    
    # Sample random points in MCP associated to specific time
    # Attribute the same value as the used loc
    # Except that expo, slope, W0 and distance to turbine
    allLoc_forRSF <- as.data.frame(do.call(rbind,
                                           llply(1:nrow(forRSF), function(b){
                                             
                                             rdmPoints <- random_points(MCPshape, n=10) %>% # Take 10 more random point in MCP
                                               bind_cols(forRSF[b,c("event.id","t")]) %>% 
                                               bind_cols(data.frame(oce::utm2lonlat(.$x_, .$y_, zone = 31, hemisphere = "N", km = FALSE))) %>% # Add the long/lat to match de strcture
                                               dplyr::rename(location.long = longitude, location.lat = latitude) %>% # Rename to match structure
                                               bind_cols(forRSF[b,c("ID","individual.local.identifier","ground.speed","date")]) %>% 
                                               extract_covariates(allRaster) %>% # Extract simultaneously expo, slope and landcover
                                               dplyr::mutate(expo_rad = (.$expo_31N*pi) / 180, # convert expo in radian
                                                             slope_rad = (.$pente_31N*pi) / 180, # convert slope in radian
                                                             colony = forRSF[b,"colony"][[1]], # Add colony
                                                             inMCP = forRSF[b,"inMCP"][[1]], 
                                                             windspeed = forRSF[b,"windspeed"][[1]], # Add windspeed
                                                             winddir_rad = forRSF[b,"winddir_rad"][[1]]) %>% # Add winddir
                                               bind_cols(do.call("rbind",llply(1:nrow(.), function(c){ # Need to do each point as round_date can end up with same timestamp
                                                 
                                                 # Estimation of orographic updraft velocity W0
                                                 coef <- ifelse(.$slope_rad[c] != 0, sin(.$slope_rad[c])*cos(.$winddir_rad[c] - .$expo_rad[c]), 0)
                                                 W0 <- .$windspeed[c] * coef
                                                 
                                                 return(W0)
                                               },.parallel = TRUE))) %>% 
                                               dplyr::rename(W0 = ...21) %>% 
                                               bind_rows(cbind(case_ = TRUE, forRSF[b,])) %>% # Add the corresponding used loc
                                               bind_cols(do.call("rbind",llply(1:nrow(.), function(d){ # Need to do each point separately otherwise the min distance is among all loc
                                                 
                                                 minDistToWT = min(pointDistance(cbind(.$location.long[d],.$location.lat[d]),
                                                                                 cbind(windTurbineTable$x, windTurbineTable$y), lonlat = TRUE))
                                                 return(minDistToWT)
                                               },.parallel = TRUE))) %>% 
                                               dplyr::rename(minDistToTurbine = ...23)
                                           },.parallel = TRUE))) %>% 
      dplyr::mutate(w = ifelse(.$case_ == TRUE, 1, 5000))
    
    # Save file
    save(allLoc_forRSF, # This will be used for RSF
         file = paste0("./Data/IndData_forRSF/allLoc_forRSF_",WindNames[v],"_q",q,".Rdata"))
    
  }
}
