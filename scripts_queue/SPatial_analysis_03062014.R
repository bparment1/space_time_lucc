####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc          #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbour to predict.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 02/07/2014 
#DATE MODIFIED: 03/06/2014
#Version: 2
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones             
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast)
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)

###### Functions used in this script

#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"

#set up the working directory
setwd(in_dir)

moore_window <- file.path(in_dir,"moore_window.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84

CRS_interp <- proj_modis_str

out_suffix <-"03062014" #output suffix for the files that are masked for quality and for 

########### START SCRIPT #################

## Read in data
moore_w <- raster(moore_window) 
projection(moore_w) <- CRS_interp

#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
#reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013.rst
#reg_var_list <- list.files(pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")

r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)

projection(r_stack) <- CRS_interp #this is modis sinusoidal
r_var_w <- crop(r_stack,moore_w) #clip (or "window") raster image to Moore subset area

r_var_dates <- subset(r_var_w,152:156) #hurricane took place between image 153 and 154
r_var_dates <- mask(r_var_dates,moore_w) #mask areas
plot(r_var_dates)

r1 <- subset(r_var_dates,1) #first date, i.e. DOY 209 (image 152)
r2 <- subset(r_var_dates,2) #second date, i.e. DOY 225 (image 153)

r_NA <-r1+r2 #this contains NA to mask values...

r1 <- mask(r1,r_NA) #use r_NA to mask both alyer
r2 <- mask(r2,r_NA)

r_s <- stack(r1,r2) #will contain the values for teh spatial model
names(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
layerNames(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster

r_s <- projectRaster(r_s,crs=CRS_WGS84) #project to latlong

#convert to polygon...for input in model
r_poly <- rasterToPolygons(r_s, fun=NULL, n=4, na.rm=TRUE, 
                        digits=12, dissolve=FALSE)

r_poly3 <- r_poly[-c(8392,11455,24685),] #these three poly that have no neighbours, remove them
                                         #note that this will change depending on input images
r_poly3$UNIQID <- 1:nrow(r_poly3) #assign unique ID
#get neighbours

writeOGR(r_poly3,layer="r_poly3",dsn=in_dir,driver="ESRI Shapefile")

r_nb <-poly2nb(r_poly3,row.names=r_poly3$UNIQID,queen=TRUE)
#sam.listw<-nb2listw(sample1.neigh, style="W",zero.policy=TRUE)
#r_nb2 <- r_nb
#r_nb2 <- r_nb[-c(8392,11455,24685)]

r_listw<-nb2listw(r_nb, style="W") #get the list of weights

# SAR model: currently not working with R 2.14 and R 3.02
sam.sar <- errorsarlm(v2 ~ v1, listw=r_listw, 
                      data=r_poly3,tol.solve=1e-36)
summary(sam.sar)

# AR model
sam.ar <- lagsarlm(v2 ~ v1, listw=r_listw, data=r_poly3,tol.solve=1e-36)
summary(sam.ar)


