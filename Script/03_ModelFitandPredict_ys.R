#___________________________________________________________________________________________________________________________________#
###
###       
###       Project: Relative probability of selecton by Lesser Kestrels in Aumelas wind farm
###
###        Data come from wild Lesser Kestrels equipped by LPO, deposited on MoveBank plateform.
### 
###       !!! Script used to model presence and relative selection of habitat of birds within MCP !!!
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


library(glmmTMB)
library(dplyr)
library(lubridate)
library(ggeffects)
library(mgcv)
library(jtools)
library(RColorBrewer)
library(ggplot2)
library(raster)
library(sf)
library(cowplot)


#Homemade functions & needed files
source("Functions/df_to_spatial_df_utm.R") #To convert long/lat to UTM
rad <- function(angle) {angle*pi/180} # convert degrees angle in radians

#________________________________________________________________
#
#### > GLMM model    ####
#
#________________________________________________________________

# Load the data sets prepared for each individual and combine them into a full one
allFiles <- list.files("./Data/IndData_VarAssociated" , all.files = F, full.names = T)

for (i in 1:length(allFiles)){
  
  load(allFiles[i])
  
  if (i == 1){
    allLoc <- IndData_allVars
  } else {
    allLoc <- rbind(allLoc, IndData_allVars)
  }
} #Dataset will all individual combined


allLoc <- allLoc %>% 
  mutate(month = factor(month(.$date)), # Add variable month
         landcover = factor(landcover), # Convert variable to factor
         colony = factor(colony),
         inMCP = factor(inMCP),
         W0 = scale(W0, center = T, scale = T)) %>% 
  filter(.$colony %in% c("Saint-P_de_M", "Saint-Pargoire", "Villeveyrac")) %>% # Filter out individuals of other colony that have less than 5 loc in MCP
  filter(.$landcover != "1") %>% # Urban area
  filter(.$inSleepingSite != "1") # Remove location from sleeping sites


# check if all individuals have data in the MCP 
isMCP <- allLoc %>% filter(.$inMCP == "1")
indiv.out <- setdiff(unique(allLoc$individual.local.identifier), unique(isMCP$individual.local.identifier))
# individuals from allLoc are not in isMCP -> do not have locations in MCP -> filter out the two individuals
allLoc <- allLoc %>% 
  filter(!.$individual.local.identifier %in% indiv.out) # delete individuals which don't have data in MCP



# Fit full model
fit_full <- glmmTMB(inMCP ~ W0 + month + colony
               + (1|individual.local.identifier),
               family=binomial(), 
               data = allLoc)


summary(fit_full)

### Fit check ---
modDrop <- drop1(fit_full, test = "Chisq")


### final model
fit <- glmmTMB(inMCP ~ W0 + month
                    + (1|individual.local.identifier),
                    family=binomial(), 
                    data = allLoc)

fit_simresid <- DHARMa::simulateResiduals(fit, plot=T) # simulation of residuals to check for oversdispersion (graph)

plot(fit_simresid)
testOutliers(fit_simresid, type = "default", plot = T)

### Summary tables ---
sjPlot::tab_model(fit,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

summary(fit)


#________________________________________________________________
#
#### > Predictions from GLMM    ####
#
#________________________________________________________________

# Load the data sets prepared for each individual and combine them into a full one
allFiles <- list.files("./Data/IndData_VarAssociated" , all.files = F, full.names = T)

for (i in 1:length(allFiles)){
  
  load(allFiles[i])
  
  if (i == 1){
    allLoc <- IndData_allVars
  } else {
    allLoc <- rbind(allLoc, IndData_allVars)
  }
} #Dataset will all individual combined


allLoc <- allLoc %>% 
  mutate(month = factor(month(.$date)),
         inMCP = factor(inMCP)) %>% 
  filter(.$colony %in% c("Saint-P_de_M", "Saint-Pargoire", "Villeveyrac")) %>% # Filter out individuals of other colony that have less than 5 loc in MCP
  filter(.$landcover != "1") %>% # Urban area
  filter(.$inSleepingSite != "1") 


### Fit best model (i.e. without colony)
final_mod <- glmmTMB(inMCP ~ W0 + month
               + (1|individual.local.identifier),
               family=binomial(), 
               data = allLoc)

summary(final_mod)
performance::check_collinearity(final_mod)


### Create and fill df for predictions plot
forProba <- data.frame(bind = c(1,2,3,4,5))
forProba$proba1 <- rep(NA,nrow(forProba))
forProba$ntot1 <-rep(NA,nrow(forProba))
forProba$month <- rep(NA,nrow(forProba))

forProba <- rbind(forProba,forProba,forProba) # For the three months
c = 0 # Counter

for (m in 5:7){
  
  for(i in 1:5){
    
    c = c+1
    forProba$month[c] <- m
    
    if (i == 1){
      dat <- allLoc[which(allLoc$W0 <= 1 & allLoc$month == m),]
      forProba$proba1[c] <- nrow(dat[which(dat$inMCP == 1),]) / nrow(dat)
      forProba$ntot1[c] <- nrow(dat)
      
    } else if (i == 5){
      dat <- allLoc[which(allLoc$W0 >= 4 & allLoc$month == m),]
      forProba$proba1[c] <- nrow(dat[which(dat$inMCP == 1),]) / nrow(dat)
      forProba$ntot1[c] <- nrow(dat)
      
    } else {

      dat <- allLoc[which(allLoc$W0 > forProba$bind[i-1] & allLoc$W0 < forProba$bind[i] & allLoc$month == m),]
      forProba$proba1[c] <- nrow(dat[which(dat$inMCP == 1),]) / nrow(dat)
      forProba$ntot1[c] <- nrow(dat)
      
    }
  }
}


### Plot predicted probability of presence in MCP --- 
mydf <- ggpredict(final_mod, type = "fe", back.transform = TRUE, ci.lvl = 0.95, terms = c("W0[0:6]","month"))

plotInMCP <- ggplot() +
  geom_line(data = mydf, aes(x, predicted, colour = group)) + 
  geom_ribbon(data = mydf, aes(x, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, show.legend = FALSE) +
  geom_point(aes(x = bind, y = proba1, size = ntot1, colour = as.factor(month)), data = forProba, show.legend = FALSE) +
  xlab("W0 [m/s]") + 
  ylab("Probability of presence in WEF-MCP") +  
  xlim(0,6) + 
  ylim(0,0.7) + 
  theme_light() +
  expand_limits(x = 15, y = 0) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) +
  scale_colour_discrete(name="Month",
                          breaks=c("5", "6", "7"),
                          labels=c("May", "June", "July"))





#________________________________________________________________
#
#### > RSF models     ####
#
#________________________________________________________________


# Load the data sets prepared for each individual and combine them into a full one
allFiles <- list.files("./Data/IndData_forRSF" , all.files = F, full.names = T)



# RSF on Sea Breeze - q1
load(allFiles[1])

# Prepare data type (add month variable and change to factors)
allLoc_forRSF <- allLoc_forRSF %>%
  mutate(month = factor(month(.$date)),
         minDistToTurbine = scale(minDistToTurbine, center = T, scale = T), # Scaled
         W0 = scale(W0, center = T, scale = T), # Scaled
         W0sq = scale(W0*W0, center = T, scale = T)) # Scaled the W0 squared

# Fit RSF model
rsf_SBq1 <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                  + (1|individual.local.identifier),
                  family=binomial(),
                  data = allLoc_forRSF,
                  weights = w)


sjPlot::tab_model(rsf_SBq1,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")


# RSF on Sea Breeze - q2
load(allFiles[2])

allLoc_forRSF <- allLoc_forRSF %>%
  mutate(month = factor(month(.$date)),
         minDistToTurbine = scale(minDistToTurbine, center = T, scale = T), # Scaled
         W0 = scale(W0, center = T, scale = T), # Scaled
         W0sq = scale(W0*W0, center = T, scale = T)) # Scaled the W0 squared

rsf_SBq2 <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                    + (1|individual.local.identifier),
                    family=binomial(),
                    data = allLoc_forRSF,
                    weights = w)


sjPlot::tab_model(rsf_SBq2,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

# RSF on Sea Breeze - q3
load(allFiles[3])

allLoc_forRSF <- allLoc_forRSF %>%
  mutate(month = factor(month(.$date)),
         minDistToTurbine = scale(minDistToTurbine, center = T, scale = T), # Scaled
         W0 = scale(W0, center = T, scale = T), # Scaled
         W0sq = scale(W0*W0, center = T, scale = T)) # Scaled the W0 squared

rsf_SBq3 <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                    + (1|individual.local.identifier),
                    family=binomial(),
                    data = allLoc_forRSF,
                    weights = w)



sjPlot::tab_model(rsf_SBq3,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

# RSF on Tramontane - q1
load(allFiles[4])

allLoc_forRSF <- allLoc_forRSF %>%
  mutate(month = factor(month(.$date)),
         minDistToTurbine = scale(minDistToTurbine, center = T, scale = T), # Scaled
         W0 = scale(W0, center = T, scale = T), # Scaled
         W0sq = scale(W0*W0, center = T, scale = T)) # Scaled the W0 squared

rsf_Tramq1 <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                    + (1|individual.local.identifier),
                    family=binomial(),
                    data = allLoc_forRSF,
                    weights = w)


sjPlot::tab_model(rsf_Tramq1,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

# RSF on Tramontane - q2
load(allFiles[5])

allLoc_forRSF <- allLoc_forRSF %>%
  mutate(month = factor(month(.$date)),
         minDistToTurbine = scale(minDistToTurbine, center = T, scale = T), # Scaled
         W0 = scale(W0, center = T, scale = T), # Scaled
         W0sq = scale(W0*W0, center = T, scale = T)) # Scaled the W0 squared

rsf_Tramq2 <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                      + (1|individual.local.identifier),
                      family=binomial(),
                      data = allLoc_forRSF,
                      weights = w)


sjPlot::tab_model(rsf_Tramq2,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

# RSF on Tramontane - q3
load(allFiles[6])

allLoc_forRSF <- allLoc_forRSF %>%
  mutate(month = factor(month(.$date)),
         minDistToTurbine = scale(minDistToTurbine, center = T, scale = T), # Scaled
         W0 = scale(W0, center = T, scale = T), # Scaled
         W0sq = scale(W0*W0, center = T, scale = T)) # Scaled the W0 squared

rsf_Tramq3 <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                      + (1|individual.local.identifier),
                      family=binomial(),
                      data = allLoc_forRSF,
                      weights = w)


sjPlot::tab_model(rsf_Tramq3,
                  transform = NULL,
                  show.est = TRUE,
                  string.est = "Estimate")

# Forest plot of each models
plot_summs(rsf_SBq1,
           rsf_SBq2,
           rsf_SBq3,
           rsf_Tramq1,
           rsf_Tramq2,
           rsf_Tramq3,
           ci_level = 0.95,
           point.shape = 19,
           colors = c("#3690C0", "#0570B0", "#034E7B", "#41AE76", "#238B45","#005824"), #brewer.pal(6, "Set2"),
           omit.coefs = c("sd__(Intercept)", "(Intercept)"),
           model.names = c("Sea breeze - 10th","Sea breeze - med","Sea breeze - 90th",
                           "Tramontane - 10th","Tramontane - med","Tramontane - 90th"), 
           coefs = c("W0" = "W0", 
                     "Squared W0" = "W0sq",
                     "Min distance to turbine" = "minDistToTurbine")) + 
  labs(x = "\n Estimates \n ") +
  scale_x_continuous(limits=c(-2,+3.2),
                     breaks = seq(-2,3,0.5)) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_line(colour = "grey93"),
        panel.grid.major = element_line(colour = "grey93"),
        strip.background = element_rect(colour = "white",
                                        fill = "white"),
        strip.text = element_text(face = "bold", size = 14)) +
  ylab("")









#________________________________________________________________
#
#### > Predict from RSF     ####
#
#________________________________________________________________


# see: https://terpconnect.umd.edu/~egurarie/research/NWT/Step09_RSF_PartIV.html#3_predict_over_a_raster


# Load shapefile to have the MCP limits
CRS.31N <- CRS("+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs +type=crs") # CRS of the project (proj4string argument)
MCPshape <- read_sf('./Data/MCP_300buffer_31N.shp')
MCPpoly <- as(st_geometry(MCPshape), "Spatial")
MCPpoly@proj4string <- CRS.31N



## First load exposition and slope rasters to create expo and slope rasters in radians at MCP dimensions

#### 1- Exposition raster ####

# mutate expo raster into radians
expoRaster <- raster("./Data/expo_31N.tif") # Exposition
MCP_expo <- mask(expoRaster, MCPpoly)
plot(MCP_expo, xlim= c(547000,553000), ylim=c(4819000,4822500)) #visual check

#save raster
#writeRaster(MCP_expo, filename = "./Output/MCP_expo.tif", format = "GTiff", progress = "text")

# Convert to radian
MCP_expo_rad <- calc(MCP_expo, fun = rad) # MCP
names(MCP_expo_rad) <- "exposition"

#save raster
#writeRaster(MCP_expo_rad, filename = "./Outputs/MCP_expo_rad.tif", format = "GTiff", progress = "text")




#### 2- Slope raster ####
#Repeate for the slope raster

# mutate slope raster into radians
slopeRaster <- raster("./Data/pente_31N.tif") # Slope
MCP_slope <- mask(slopeRaster, MCPpoly)
plot(MCP_slope, xlim= c(547000,553000), ylim=c(4819000,4822500)) #visual check

#save raster
#writeRaster(MCP_slope, filename = "./Data/MCP_slope.tif", format = "GTiff", progress = "text")

# Convert to radian
MCP_slope_rad <- calc(MCP_slope, fun = rad) # MCP
names(MCP_slope_rad) <- "slope"

#save raster
#writeRaster(MCP_slope_rad, filename = "./Outputs/MCP_slope_rad.tif", format = "GTiff", progress = "text")




#### 3- Tramontane rasters ####

# Object combining wind information
WindNames <- c("Tramontane", "SeaBreeze")

WindInfo <- vector("list")
WindInfo[[1]] <- list(name = WindNames[1],
                      angles = 303.2, # mean angle Tram
                      speed = c(4.4, 7.5, 10.3)) #10%, med, 90% of wind speed distribution
WindInfo[[2]] <- list(name = WindNames[2],
                      angles = 132.2, # mean angle Sea breeze
                      speed = c(1.8, 4.2, 6.3))


## Create a rasters with mean tramontane angle (i.e northeast wind) in radian at MCP dimensions.
## and wind speed at 10th, median and 90th percentile

# Raster Tram mean direction at MCP dimensions
tram_mean_dir_rad <- raster(extent(MCP_expo_rad), nrow = nrow(MCP_expo_rad), ncol = ncol(MCP_expo_rad))
values(tram_mean_dir_rad) <- WindInfo[[1]]$angles %>% rad(.)
tram_mean_dir_rad <- mask(tram_mean_dir_rad, MCP_expo_rad) 
names(tram_mean_dir_rad) <- "wind_direction"
plot(tram_mean_dir_rad, xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check

#save raster
#writeRaster(tram_mean_dir_rad, filename = "./Outputs/MCP_Tram_meanDir.tif", format = "GTiff", progress = "text")



## Create rasters with tramontane wind speed at at 10th, median and 90th percentile (m/s) at MCP dimensions.
tram_speed <- raster(extent(MCP_expo_rad), nrow = nrow(MCP_expo_rad), ncol = ncol(MCP_expo_rad))

for (s in 1:3){
  
  values(tram_speed) <- WindInfo[[1]]$speed[s]
  tram_speed <- mask(tram_speed, MCP_expo_rad)
  #plot(tram_speed, xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check
  
  #save raster
  writeRaster(tram_speed, filename = paste0("./Outputs/MCP_Tram_speed_q",s,".tif"), format = "GTiff", progress = "text")
}



#### 4- Sea breeze (SB) rasters ####

# Raster SB mean direction at MCP dimensions
SB_mean_dir_rad <- raster(extent(MCP_expo_rad), nrow = nrow(MCP_expo_rad), ncol = ncol(MCP_expo_rad))
values(SB_mean_dir_rad) <- WindInfo[[2]]$angles %>% rad(.)
SB_mean_dir_rad <- mask(SB_mean_dir_rad, MCP_expo_rad) 
names(SB_mean_dir_rad) <- "wind_direction"
#plot(SB_mean_dir_rad, xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check

#save raster
#writeRaster(SB_mean_dir_rad, filename = "./Outputs/MCP_SB_meanDir.tif", format = "GTiff", progress = "text")



## Create rasters with SB wind speed at 10th, median and 90th percentile (m/s) at MCP dimensions.
SB_speed <- raster(extent(MCP_expo_rad), nrow = nrow(MCP_expo_rad), ncol = ncol(MCP_expo_rad))

for (s in 1:3){
  
  values(SB_speed) <- WindInfo[[2]]$speed[s]
  SB_speed <- mask(SB_speed, MCP_expo_rad)
  #plot(tram_speed, xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check
  
  #save raster
  writeRaster(SB_speed, filename = paste0("./Outputs/MCP_SB_speed_q",s,".tif"), format = "GTiff", progress = "text")
}






#### 5- W0 & W0sq rasters  ####
### compute the orographic updraft velocity W0 on raster for each conditions

# For the tramontane
coef <- function(teta,alpha,beta) {sin(teta)*cos(alpha-beta)}
windNames <- c("Tram", "SB")
rastName <- c("./Outputs/MCP_Tram_meanDir.tif","./Outputs/MCP_SB_meanDir.tif")
rastWindNames <- c("./Outputs/MCP_Tram_speed_q1.tif","./Outputs/MCP_Tram_speed_q2.tif","./Outputs/MCP_Tram_speed_q3.tif",
                   "./Outputs/MCP_SB_speed_q1.tif","./Outputs/MCP_SB_speed_q2.tif","./Outputs/MCP_SB_speed_q3.tif")
c=0

for (v in 1:2){
  
  # Load the mean direction of wind
  MCP_windDir <- raster(rastName[v])
  
  for (s in 1:3){ # Loop on at 10th, median and 90th percentile wind speed for each wind
    c <- c+1 #counter
    
    #Load the wind speed raster
    MCP_windspeed <- raster(rastWindNames[c])
    
    # W0
    MCP_W0 <- overlay(MCP_slope_rad, MCP_windDir, MCP_expo_rad, fun = coef) %>%
              overlay(.,MCP_windspeed, fun = function(x,y) x*y)
    names(MCP_W0) <- "W0"
    
    #plot(MCP_W0, xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check
    
    # W0sq
    MCP_W00 <- MCP_W0 %>% reclassify(c(-Inf,0,0))
    MCP_W0sq <- calc(MCP_W00, fun=function(x){x^2}) # create W0 quadratic raster
    names(MCP_W0sq) <- "W0sq" # rename to match name var expected by the model
    
    #plot(MCP_W0sq , xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check

    # Scaling for future predictions
    MCP_W0 <- calc(MCP_W0, fun = scale)
    MCP_W0sq <- calc(MCP_W0sq, fun = scale) 
    
    #save raster
    writeRaster(MCP_W0, filename = paste0("./Outputs/MCP_W0_",windNames[v],"_q",s,".tif"), format = "GTiff", progress = "text")
    writeRaster(MCP_W0sq, filename = paste0("./Outputs/MCP_W0sq_",windNames[v],"_q",s,".tif"), format = "GTiff", progress = "text")
    
  }
}





#### 6- minDistToTurbine raster ####

## Create a raster with the pixel distance to the nearest wind turbine
# Load table with turbines locations + add latlong coordinates
windTurbineTable <-  read.csv("./Data/locWindTurbinesAumelas.csv", header = TRUE)

#Transform coordinates to UTM
df_to_spatial_df_utm(dataframe=windTurbineTable, 
                     column_longitude="x", 
                     column_latitude="y", 
                     utm_zone=31, hemisphere="N")


windturbines.31N <- cbind(windTurbineTable, dataframe.sp_utm@coords) %>%
  dplyr::select(-gid, -id_zde, -id_poste, -x_icpe, -y_icpe, -sys_coord,-alt_base, 
                -etat_mat, -date_real, -date_demol,-en_service,-srce_geom,-d_service,
                -e_adm_mat, -dept,-obs,-xcoord, -ycoord, -i_pc_anc, -c_mat_diff, -n_mat_diff) %>%
  dplyr::rename(lon.31N = longitude,
                lat.31N = latitude) # dataframe used later for windturbine plots


# Create a new raster with same dimensions 
new_raster <- raster(extent(MCP_expo_rad), resolution = c(24.98211, 24.98211), crs = CRS.31N)

# compute distances 
MCP_minDistToTurbine <- distanceFromPoints(new_raster, dataframe.sp_utm)
#plot(MCP_minDistToTurbine, xlim= c(547000,553000), ylim=c(4819000,4822500)) # Visual check

# Scale the raster
MCP_minDistToTurbine <- calc(MCP_minDistToTurbine, fun = scale) 
names(MCP_minDistToTurbine) <- "minDistToTurbine"

# Save raster
# writeRaster(MCP_minDistToTurbine, filename = paste0("./Outputs/MCP_minDistToTurbine.tif"), format = "GTiff", progress = "text")




#### 7- Predictions from RSF ####

# Prepare names list to call them in the loop
MCP_W0_Names <- c("./Outputs/MCP_W0_Tram_q1.tif","./Outputs/MCP_W0_Tram_q2.tif","./Outputs/MCP_W0_Tram_q3.tif",
                  "./Outputs/MCP_W0_SB_q1.tif","./Outputs/MCP_W0_SB_q2.tif","./Outputs/MCP_W0_SB_q3.tif")

MCP_W0sq_Names <- c("./Outputs/MCP_W0sq_Tram_q1.tif","./Outputs/MCP_W0sq_Tram_q2.tif","./Outputs/MCP_W0sq_Tram_q3.tif",
                    "./Outputs/MCP_W0sq_SB_q1.tif","./Outputs/MCP_W0sq_SB_q2.tif","./Outputs/MCP_W0sq_SB_q3.tif")

indData_names <- c("./Data/IndData_forRSF/allLoc_forRSF_Tramontane_q1.Rdata",
                   "./Data/IndData_forRSF/allLoc_forRSF_Tramontane_q2.Rdata",
                   "./Data/IndData_forRSF/allLoc_forRSF_Tramontane_q3.Rdata",
                   "./Data/IndData_forRSF/allLoc_forRSF_SeaBreeze_q1.Rdata",
                   "./Data/IndData_forRSF/allLoc_forRSF_SeaBreeze_q2.Rdata",
                   "./Data/IndData_forRSF/allLoc_forRSF_SeaBreeze_q3.Rdata")

# Raster distance
MCP_minDistToTurbine <- raster("./Outputs/MCP_minDistToTurbine.tif")
names(MCP_minDistToTurbine) <- "minDistToTurbine"

# False raster individual -> necessary to fit model for predictions
MCP_individual.local.identifier <- raster(extent(MCP_minDistToTurbine),resolution = c(24.98211, 24.98211), crs = CRS.31N)
values(MCP_individual.local.identifier) <- as.factor("a")
names(MCP_individual.local.identifier) <-  "individual.local.identifier"

windName <- "Tram"
c = 0

for (i in 1:6){ # Loop to predict in each condition of wind
  
  # Load raster values W0 / W0sq and stack them with minDist and ind
  MCP_W0 <- raster(MCP_W0_Names[i])
  names(MCP_W0) <- "W0"
  
  MCP_W0sq <- raster(MCP_W0sq_Names[i])
  names(MCP_W0sq) <- "W0sq"
  
  all_rasters <- stack(MCP_W0, MCP_W0sq, MCP_minDistToTurbine, MCP_individual.local.identifier)
  
  
  # Load data for specific scenario
  # run the model with only 80% of the data (train_set) for visual cross validation (20% remaining)
  load(indData_names[i])
  indData_train <- allLoc_forRSF %>% 
    filter(event.id %in% sample(unique(.$event.id),round(length(unique(.$event.id))*0.8))) %>% #80% of the dataset
    mutate(month = factor(month(.$date)),
           minDistToTurbine = scale(minDistToTurbine, center = T, scale = T),
           W0 = scale(W0, center = T, scale = T),
           W0sq = scale(W0*W0, center = T, scale = T))
  
  indData_test <- allLoc_forRSF %>% 
    filter(!event.id %in% unique(indData_train$event.id)) %>%  #20% remaining for visual check
    mutate(month = factor(month(.$date)),
           minDistToTurbine = scale(minDistToTurbine, center = T, scale = T),
           W0 = scale(W0, center = T, scale = T),
           W0sq = scale(W0*W0, center = T, scale = T))
  
  # model
  rsf_mod <- glmmTMB(case_ ~ W0 + W0sq + minDistToTurbine
                     + (1|individual.local.identifier),
                     family=binomial(),
                     data = indData_train)
  
  # prediction under this scenario
  RSF.prediction <- predict(all_rasters, rsf_mod, type = "response", re.form = NA, allow.new.levels=T) # re.form = NA to set random effect to 0 (indiv local identifier = null)
  #plot(RSF.prediction, xlim= c(547000,553000), ylim=c(4819000,4822500)) Visual check
  
  RSF.prediction_norm <- RSF.prediction/max(values(RSF.prediction), na.rm = TRUE) # Rescale proba between 0-1
  # plot(RSF.prediction_norm, xlim= c(547000,553000), ylim=c(4819000,4822500)) #Visual check
  # plot(dataframe.sp_utm, add = TRUE)
  
  if(i == 4){c = 0}
  if(i >= 4){windName <- "SB"}
  c = c+1
  
  # Save raster
  writeRaster(RSF.prediction_norm, filename = paste0("./Outputs/RSF.prediction_",windName,"_q",c,".tif"), format = "GTiff", progress = "text")
  save(indData_test, file = paste0("./Outputs/IndData_test_",windName,"_q",c,".RData"))
}




#### 8- Prediction maps and risky turbines ####

# change to spatial point data frame (SPDF) - Needed after
windturbines.31N_spdf <- SpatialPointsDataFrame(cbind(windturbines.31N$lon.31N, windturbines.31N$lat.31N), 
                                                data = windturbines.31N)



##### Map TRAMONTANE 10th ####
RSF.prediction_norm <- raster("./Outputs/RSF.prediction_Tram_q1.tif")
load("./Outputs/IndData_test_Tram_q1.RData")

### Turbines at risks p > 50
# change raster to data for ggplot
data.pred <- rasterToPoints(RSF.prediction_norm) %>%
              data.frame(.) %>%
              dplyr::rename(predictions = layer,
                   longitude = x,
                   latitude = y)


#Look for the most risky wind turbines #
# Extract maximum selection probabilities in a 100m buffer around each turbine (200m?)
prob_max <- raster::extract(RSF.prediction_norm,  # raster layer
                            windturbines.31N_spdf,# wt SPDF for buffers
                            buffer = 50, # buffer size, units depend on CRS
                            fun=max,      # what to value to extract
                            df=TRUE)      # return a dataframe 

# grab the names of the plots from the centroid_spdf
prob_max <- prob_max[,-1] %>% # delete first column
  bind_cols(windturbines.31N_spdf$numero) %>% # Binding wind turbine numbers
  dplyr::rename(max_prob = ...1,
         WT_ID = ...2)  #fix the column names

windturbines.31N_spdf <- as.data.frame(windturbines.31N_spdf)
windturbines.31N_spdf <- merge(windturbines.31N_spdf, prob_max, by.x = 'numero', by.y = 'WT_ID')



#### Visualize risky turbines - 50%
tram_q1_50pct <- ggplot() +
                geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
                geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.50,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
                geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.50,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
                geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
                theme_bw() +
                scale_fill_viridis_c() +
                xlab("") +
                ylab("Latitude") + 
                annotate("text", x=548700, y=4822500, label= "Wind speed = 4.4 m/s") +
                theme(legend.position='none',
                      axis.text.x = element_text(size = 10),
                      axis.text.y = element_text(size = 10),
                      axis.title = element_text(size = 12),
                      strip.text = element_text(face = "bold", size = 14))


#### Visualize risky turbines - 75%
tram_q1_75pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.75,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.75,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("Longitude") +
  ylab("Latitude") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 4.4 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))




##### Map TRAMONTANE med ####
RSF.prediction_norm <- raster("./Outputs/RSF.prediction_Tram_q2.tif")
load("./Outputs/IndData_test_Tram_q2.Rdata")

### Turbines at risks p > 50
# change raster to data for ggplot
data.pred <- rasterToPoints(RSF.prediction_norm) %>%
  data.frame(.) %>%
  dplyr::rename(predictions = layer,
                longitude = x,
                latitude = y)


#Look for the most risky wind turbines #
# Extract maximum selection probabilities in a 100m buffer around each turbine (200m?)
windturbines.31N_spdf <- SpatialPointsDataFrame(cbind(windturbines.31N$lon.31N, windturbines.31N$lat.31N), 
                                                data = windturbines.31N)

prob_max <- raster::extract(RSF.prediction_norm,  # raster layer
                            windturbines.31N_spdf,# wt SPDF for buffers
                            buffer = 50, # buffer size, units depend on CRS
                            fun=max,      # what to value to extract
                            df=TRUE)      # return a dataframe 

# grab the names of the plots from the centroid_spdf
prob_max <- prob_max[,-1] %>% # delete first column
  bind_cols(windturbines.31N_spdf$numero) %>% # Binding wind turbine numbers
  dplyr::rename(max_prob = ...1,
                WT_ID = ...2)  #fix the column names

windturbines.31N_spdf <- as.data.frame(windturbines.31N_spdf)
windturbines.31N_spdf <- merge(windturbines.31N_spdf, prob_max, by.x = 'numero', by.y = 'WT_ID')



#### Visualize risky turbines - 50%
tram_q2_50pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.50,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.50,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 7.5 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


#### Visualize risky turbines - 75%
tram_q2_75pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.75,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.75,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("Longitude") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 7.5 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))



##### Map TRAMONTANE 90th ####
RSF.prediction_norm <- raster("./Outputs/RSF.prediction_Tram_q3.tif")
load("./Outputs/IndData_test_Tram_q3.Rdata")

### Turbines at risks p > 50
# change raster to data for ggplot
data.pred <- rasterToPoints(RSF.prediction_norm) %>%
  data.frame(.) %>%
  dplyr::rename(predictions = layer,
                longitude = x,
                latitude = y)


#Look for the most risky wind turbines #
# Extract maximum selection probabilities in a 100m buffer around each turbine (200m?)
windturbines.31N_spdf <- SpatialPointsDataFrame(cbind(windturbines.31N$lon.31N, windturbines.31N$lat.31N), 
                                                data = windturbines.31N)

prob_max <- raster::extract(RSF.prediction_norm,  # raster layer
                            windturbines.31N_spdf,# wt SPDF for buffers
                            buffer = 50, # buffer size, units depend on CRS
                            fun=max,      # what to value to extract
                            df=TRUE)      # return a dataframe 

# grab the names of the plots from the centroid_spdf
prob_max <- prob_max[,-1] %>% # delete first column
  bind_cols(windturbines.31N_spdf$numero) %>% # Binding wind turbine numbers
  dplyr::rename(max_prob = ...1,
                WT_ID = ...2)  #fix the column names

windturbines.31N_spdf <- as.data.frame(windturbines.31N_spdf)
windturbines.31N_spdf <- merge(windturbines.31N_spdf, prob_max, by.x = 'numero', by.y = 'WT_ID')



#### Visualize risky turbines - 50%
tram_q3_50pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.50,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.50,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 10.3 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


#### Visualize risky turbines - 75%
tram_q3_75pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.75,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.75,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("Longitude") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 10.3 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


### Combined plot of prediction map - Tram---
plot_grid(tram_q1_50pct,tram_q2_50pct,tram_q3_50pct, 
          tram_q1_75pct,tram_q2_75pct,tram_q3_75pct,
          labels=c("A","B","C","D","E","F"), 
          ncol = 3, nrow = 2)
















##### Map SB 10th ####
RSF.prediction_norm <- raster("./Outputs/RSF.prediction_SB_q1.tif")
load("./Outputs/IndData_test_SB_q1.Rdata")

### Turbines at risks p > 50
# change raster to data for ggplot
data.pred <- rasterToPoints(RSF.prediction_norm) %>%
  data.frame(.) %>%
  dplyr::rename(predictions = layer,
                longitude = x,
                latitude = y)


#Look for the most risky wind turbines #
# Extract maximum selection probabilities in a 100m buffer around each turbine (200m?)
windturbines.31N_spdf <- SpatialPointsDataFrame(cbind(windturbines.31N$lon.31N, windturbines.31N$lat.31N), 
                                                data = windturbines.31N)

prob_max <- raster::extract(RSF.prediction_norm,  # raster layer
                            windturbines.31N_spdf,# wt SPDF for buffers
                            buffer = 50, # buffer size, units depend on CRS
                            fun=max,      # what to value to extract
                            df=TRUE)      # return a dataframe 

# grab the names of the plots from the centroid_spdf
prob_max <- prob_max[,-1] %>% # delete first column
  bind_cols(windturbines.31N_spdf$numero) %>% # Binding wind turbine numbers
  dplyr::rename(max_prob = ...1,
                WT_ID = ...2)  #fix the column names

windturbines.31N_spdf <- as.data.frame(windturbines.31N_spdf)
windturbines.31N_spdf <- merge(windturbines.31N_spdf, prob_max, by.x = 'numero', by.y = 'WT_ID')



#### Visualize risky turbines - 50%
SB_q1_50pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.50,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.50,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("") +
  ylab("Latitude") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 1.8 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


#### Visualize risky turbines - 75%
SB_q1_75pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.75,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.75,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("Longitude") +
  ylab("Latitude") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 1.8 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))




##### Map SB med ####
RSF.prediction_norm <- raster("./Outputs/RSF.prediction_SB_q2.tif")
load("./Outputs/IndData_test_SB_q2.Rdata")

### Turbines at risks p > 50
# change raster to data for ggplot
data.pred <- rasterToPoints(RSF.prediction_norm) %>%
  data.frame(.) %>%
  dplyr::rename(predictions = layer,
                longitude = x,
                latitude = y)


#Look for the most risky wind turbines #
# Extract maximum selection probabilities in a 100m buffer around each turbine (200m?)
windturbines.31N_spdf <- SpatialPointsDataFrame(cbind(windturbines.31N$lon.31N, windturbines.31N$lat.31N), 
                                                data = windturbines.31N)

prob_max <- raster::extract(RSF.prediction_norm,  # raster layer
                            windturbines.31N_spdf,# wt SPDF for buffers
                            buffer = 50, # buffer size, units depend on CRS
                            fun=max,      # what to value to extract
                            df=TRUE)      # return a dataframe 

# grab the names of the plots from the centroid_spdf
prob_max <- prob_max[,-1] %>% # delete first column
  bind_cols(windturbines.31N_spdf$numero) %>% # Binding wind turbine numbers
  dplyr::rename(max_prob = ...1,
                WT_ID = ...2)  #fix the column names

windturbines.31N_spdf <- as.data.frame(windturbines.31N_spdf)
windturbines.31N_spdf <- merge(windturbines.31N_spdf, prob_max, by.x = 'numero', by.y = 'WT_ID')



#### Visualize risky turbines - 50%
SB_q2_50pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.50,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.50,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 4.2 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


#### Visualize risky turbines - 75%
SB_q2_75pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.75,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.75,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("Longitude") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 4.2 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))



##### Map SB 90th ####
RSF.prediction_norm <- raster("./Outputs/RSF.prediction_SB_q3.tif")
load("./Outputs/IndData_test_SB_q3.Rdata")

### Turbines at risks p > 50
# change raster to data for ggplot
data.pred <- rasterToPoints(RSF.prediction_norm) %>%
  data.frame(.) %>%
  dplyr::rename(predictions = layer,
                longitude = x,
                latitude = y)


#Look for the most risky wind turbines #
# Extract maximum selection probabilities in a 100m buffer around each turbine (200m?)
windturbines.31N_spdf <- SpatialPointsDataFrame(cbind(windturbines.31N$lon.31N, windturbines.31N$lat.31N), 
                                                data = windturbines.31N)

prob_max <- raster::extract(RSF.prediction_norm,  # raster layer
                            windturbines.31N_spdf,# wt SPDF for buffers
                            buffer = 50, # buffer size, units depend on CRS
                            fun=max,      # what to value to extract
                            df=TRUE)      # return a dataframe 

# grab the names of the plots from the centroid_spdf
prob_max <- prob_max[,-1] %>% # delete first column
  bind_cols(windturbines.31N_spdf$numero) %>% # Binding wind turbine numbers
  dplyr::rename(max_prob = ...1,
                WT_ID = ...2)  #fix the column names

windturbines.31N_spdf <- as.data.frame(windturbines.31N_spdf)
windturbines.31N_spdf <- merge(windturbines.31N_spdf, prob_max, by.x = 'numero', by.y = 'WT_ID')



#### Visualize risky turbines - 50%
SB_q3_50pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.50,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.50,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 6.3 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


#### Visualize risky turbines - 75%
SB_q3_75pct <- ggplot() +
  geom_raster(data = data.pred, aes(x = longitude, y = latitude, fill = predictions)) +
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob < 0.75,], aes(x = lon.31N, y = lat.31N), color = "white", shape = 9 ,size = 2)+
  geom_point(data = windturbines.31N_spdf[windturbines.31N_spdf$max_prob >= 0.75,], aes(x = lon.31N, y = lat.31N), fill = "#FF0000", color = "#FF0000", shape = 23 ,size = 3)+
  geom_point(data = indData_test[which(indData_test$case_ == TRUE),], aes(x = x_, y = y_), color = "#990000", size = 1)+
  theme_bw() +
  scale_fill_viridis_c() +
  xlab("Longitude") +
  ylab("") + 
  annotate("text", x=548700, y=4822500, label= "Wind speed = 6.3 m/s") +
  theme(legend.position='none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14))


### Combined plot of prediction map - SB ---
plot_grid(SB_q1_50pct,SB_q2_50pct,SB_q3_50pct, 
          SB_q1_75pct,SB_q2_75pct,SB_q3_75pct,
          labels=c("A","B","C","D","E","F"), 
          ncol = 3, nrow = 2)

