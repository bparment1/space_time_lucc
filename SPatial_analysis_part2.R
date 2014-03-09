####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc          #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbour to predict.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 03/09/2014
#Version: 1
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
#in_dir <-"home/parmentier/Data/Space_Time"
in_dir <- "C:/Users/mmmillones/Dropbox/Space_Time/"

#set up the working directory
setwd(in_dir)

moore_window <- file.path(in_dir,"moore_window.rst")
test_shp_fname <-"C:/Users/mmmillones/Dropbox/Space_Time/GEODA_Analysis/TEST_SPATIAL_ONE.shp"
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

out_suffix <-"03092014" #output suffix for the files that are masked for quality and for 

########### START SCRIPT #################

## Read in data
moore_w <- raster(moore_window) 
projection(moore_w) <- CRS_interp
test_shp <- readOGR(dsn=dirname(test_shp_fname),
                    layer=gsub(extension(basename(test_shp_fname)),"",basename(test_shp_fname)))

#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
#reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013.rst
#reg_var_list <- list.files(pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")

r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
#levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)

projection(r_stack) <- CRS_interp #this is modis sinusoidal
r_var <- subset(r_stack,153:154) #date before hurricane is 153
r_clip <- moore_w

test_df <- as.data.frame(test_shp)
lm_mod<-lm(v2~W_V1+v1,data=test_df) #tested model
test_shp$pred_R <- lm_mod$fitted.values
test_shp$pred_R_residuals <- lm_mod$residuals
proj4string(test_shp) <-CRS_WGS84

r_var_w <- crop(r_var,r_clip) #clip (or "window") raster image to Moore subset area
r_var_w <- mask(r_var_w,r_clip) #mask areas
#plot(r_var_dates)
  
r1 <- subset(r_var_w,1) #first date, i.e. DOY 209 (image 152)
r2 <- subset(r_var_w,2) #second date, i.e. DOY 225 (image 153)
  
r_NA <-r1+r2 #this contains NA to mask values...
r_NA <- r_NA > -9999

r1 <- mask(r1,r_NA) #use r_NA to mask both alyer
r2 <- mask(r2,r_NA)
  
r_s <- stack(r1,r2) #will contain the values for teh spatial model
#names(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
names(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
  
if(!is.null(proj_str)){
    r_s <- projectRaster(r_s,crs=CRS_WGS84) #project to latlong
}
  
r_spat_pred <- rasterize(test_shp,r_s,field="pred_R") #this is the prediction from lm model

plot(r_spat_pred,subset(r_s,2)) #quick visualization...

#Now write out gal and shp files
out_dir <-in_dir
out_fname <-file.path(out_dir,paste("test_shp_in_R2","_",out_suffix,".shp",sep=""))  
writeOGR(test_shp,dsn= dirname(out_fname),layer= sub(".shp","",basename(out_fname)), 
         driver="ESRI Shapefile",overwrite_layer=TRUE)

names(test_shp)
writeRaster(r_spat_pred,filename="test_pred_R.tif") #writeRaster(r_spat_pred,filename="test_pred_R.rst" it is not working in this versionof R yet worked for previous ones
# to use GDAL
writeRaster(r_spat_pred,filename="test_pred_R.rst")
system("gdal_tranlate_test_R.tif test_pred_R.rst")
            # to translante a stack of raster files to rst from tiff
            # system("gdal_translate *.tif *.rst")