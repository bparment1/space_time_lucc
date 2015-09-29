####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc          #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbour to predict.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 02/07/2014 
#DATE MODIFIED: 03/07/2014
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

create_sp_poly_spatial_reg <- function(r_var,r_clip=NULL,proj_str=NULL,out_suffix,out_dir="."){
  #Function to remove NA accross two raster layers and create no empty neighbour list
  #inputs: 
  #1)r_var: stack of two raster layer
  #2)r_clip: layer used to clip, if null there is no clipping
  #3)proj_str: to reproject if needed, if null no reprojection
  #4)out_dir: output directory path
  #5)out_suffix: suffix added to output file name
  #outputs: list of four elements 
  #r_listw: list of weights (Queen)
  #r_poly_name: shapefile name screeened for NA and no neighbours features
  #r_nb_name: neighbour object
  #zero_nb: feature polygon with no neighbours that were removed

  ## START function
  if(!is.null(r_clip)){
      r_var_w <- crop(r_var,r_clip) #clip (or "window") raster image to Moore subset area
      r_var_w <- mask(r_var_w,r_clip) #mask areas
  }else{
    r_var_w <- r_var
  }
  #plot(r_var_dates)
  
  r1 <- subset(r_var_w,1) #first date, i.e. DOY 209 (image 152)
  r2 <- subset(r_var_w,2) #second date, i.e. DOY 225 (image 153)
  
  r_NA <-r1+r2 #this contains NA to mask values...
  
  r1 <- mask(r1,r_NA) #use r_NA to mask both alyer
  r2 <- mask(r2,r_NA)
  
  r_s <- stack(r1,r2) #will contain the values for teh spatial model
  names(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
  #layerNames(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
  
  if(!is.null(proj_str)){
    r_s <- projectRaster(r_s,crs=CRS_WGS84) #project to latlong
  }
  
  #convert to polygon...for input in model
  r_poly <- rasterToPolygons(r_s, fun=NULL, n=4, na.rm=TRUE, 
                             digits=12, dissolve=FALSE)
  
  #Now find polygons wiht no neighbours..
  r_poly$UNIQID <- 1:nrow(r_poly) #assign unique ID
  r_nb <-poly2nb(r_poly,row.names=r_poly$UNIQID,queen=TRUE)
  
  ID_selected <- which(card(r_nb)==0)  #Find entities with zero neighbour i.e. cardinality==0
  
  if(length(ID_selected>0)){
      r_nb_sub <- subset(r_nb, subset=card(r_nb) > 0)
      r_poly_sub <- r_poly[-ID_selected,] #these three poly that have no neighbours, remove them
      #note that this will change depending on input images
      r_poly_sub$UNIQID <- 1:nrow(r_poly_sub) #assign unique ID
      r_poly <- r_poly_sub
      r_nb <- r_nb_sub
  }
    
  r_listw<-nb2listw(r_nb, style="W") #get the list of weights
  
  #Now write out gal and shp files
  out_fname <-file.path(out_dir,paste("r_poly","_",out_suffix,".shp",sep=""))  
  writeOGR(r_poly,dsn= dirname(out_fname),layer= sub(".shp","",basename(out_fname)), 
           driver="ESRI Shapefile",overwrite_layer=TRUE)
  #writeOGR(r_poly_sub,layer="r_poly",dsn=out_dir,driver="ESRI Shapefile")

  #save gal write.gal
  out_gal_fname <-file.path(out_dir,paste("r_poly","_",out_suffix,".gal",sep=""))  
  write.nb.gal(nb=r_nb, file=out_gal_fname, oldstyle=FALSE, shpfile=out_fname)

  nb_obj<- list(r_listw,out_fname,out_gal_fname,ID_selected)
  names(nb_obj) <- c("r_listw","r_poly_name","r_nb_name","zero_nb")
  
  return(nb_obj)
  
}

create_sp_poly_one_date_spatial_reg <- function(r_var,r_clip=NULL,proj_str=NULL,out_suffix,out_dir="."){
  #Function to remove NA accross two raster layers and create no empty neighbour list
  #inputs: 
  #1)r_var: stack of two raster layer
  #2)r_clip: layer used to clip, if null there is no clipping
  #3)proj_str: to reproject if needed, if null no reprojection
  #4)out_dir: output directory path
  #5)out_suffix: suffix added to output file name
  #outputs: list of four elements 
  #r_listw: list of weights (Queen)
  #r_poly_name: shapefile name screeened for NA and no neighbours features
  #r_nb_name: neighbour object
  #zero_nb: feature polygon with no neighbours that were removed

  ## START function
  if(!is.null(r_clip)){
      r_var_w <- crop(r_var,r_clip) #clip (or "window") raster image to Moore subset area
      r_var_w <- mask(r_var_w,r_clip) #mask areas
  }else{
    r_var_w <- r_var
  }
  #plot(r_var_dates)

  #convert to polygon...for input in model
  r_poly <- rasterToPolygons(r_var_w, fun=NULL, n=4, na.rm=TRUE, 
                             digits=12, dissolve=FALSE)
  
  #Now find polygons wiht no neighbours..
  r_poly$UNIQID <- 1:nrow(r_poly) #assign unique ID
  r_nb <-poly2nb(r_poly,row.names=r_poly$UNIQID,queen=TRUE)
  
  ID_selected <- which(card(r_nb)==0)  #Find entities with zero neighbour i.e. cardinality==0
  
  if(length(ID_selected>0)){
      r_nb_sub <- subset(r_nb, subset=card(r_nb) > 0)
      r_poly_sub <- r_poly[-ID_selected,] #these three poly that have no neighbours, remove them
      #note that this will change depending on input images
      r_poly_sub$UNIQID <- 1:nrow(r_poly_sub) #assign unique ID
      r_poly <- r_poly_sub
      r_nb <- r_nb_sub
  }
    
  r_listw<-nb2listw(r_nb, style="W") #get the list of weights
  
  #Now write out gal and shp files
  out_fname <-file.path(out_dir,paste("r_poly","_",out_suffix,".shp",sep=""))  
  writeOGR(r_poly,dsn= dirname(out_fname),layer= sub(".shp","",basename(out_fname)), 
           driver="ESRI Shapefile",overwrite_layer=TRUE)
  #writeOGR(r_poly_sub,layer="r_poly",dsn=out_dir,driver="ESRI Shapefile")

  #save gal write.gal
  out_gal_fname <-file.path(out_dir,paste("r_poly","_",out_suffix,".gal",sep=""))  
  write.nb.gal(nb=r_nb, file=out_gal_fname, oldstyle=FALSE, shpfile=out_fname)

  nb_obj<- list(r_listw,out_fname,out_gal_fname,ID_selected)
  names(nb_obj) <- c("r_listw","r_poly_name","r_nb_name","zero_nb")
  
  return(nb_obj)
  
}


#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
#in_dir <-"home/parmentier/Data/Space_Time"
in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"
out_dir <- in_dir
#set up the working directory
setwd(in_dir)

moore_window <- file.path(in_dir,"moore_window.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84

CRS_interp <- proj_modis_str

out_suffix <-"03072014" #output suffix for the files that are masked for quality and for 

########### START SCRIPT #################

## Read in data
raster_w <- raster(moore_window) 
projection(raster_w) <- CRS_interp

#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
#reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013.rst
#reg_var_list <- list.files(pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")

r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
projection(r_stack) <- CRS_interp #this is modis sinusoidal

#levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)

###### PART 1: screened out for NA and no neighbour entitty

#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
r_var <- subset(r_stack,152:153)
r_clip <- raster_w
proj_str<- CRS_WGS84 
out_suffix_s <- paste("mt2_mt1",out_suffix,sep="_")
#out_dir and out_suffix set earlier

nb_obj_for_pred_t_minus1 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

### Prepare dataset 2 for function: date t-1 (mt1) and t+1 (pt1) before and after hurricane (t0)
#this is for prediction t+1
r_var <- subset(r_stack,153:154)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("mt1_pt1",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_t1 <- create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

### Prepare dataset 3 for function: date t+1 (pt1) and t+2 (pt2) before and after hurricane (t0)
#this is for prediction t+2
r_var <- subset(r_stack,154:155)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("pt1_pt2",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_t2 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

### Prepare dataset 4 for function: date t+2 (pt2) and t+3 (pt3) before and after hurricane (t0)
#this is for prediction t+3
r_var <- subset(r_stack,155:156)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("pt2_pt3",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_t3 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

###### PART 2: run spatial regression models

### Now run spatial regression models

mod_spat_param_list <- list(nb_obj_for_pred_t_minus1,nb_obj_for_pred_t1,nb_obj_for_pred_t2,nb_obj_for_pred_t3)

#does not run right now

#for(i in 1:length(mod_spat_param_list)){
#  r_poly_fname <- mod_spat_param_list[[i]]$r_poly_name
#  r_listw <- mod_spat_param_list[[i]]$r_listw
#    
#  #read...
#  r_poly <- readOGR(dsn=)
#  r_poly <- readOGR(dsn=dirname(r_poly_fname),
#                    layer=gsub(extension(basename(r_poly_fname)),"",basename(r_poly_fname))) 
#  #readd
#  # SAR model: currently not working with R 2.14 and R 3.02
#  sam.sar <- errorsarlm(v2 ~ v1, listw=r_listw, 
#                      data=r_poly,tol.solve=1e-36)
#  summary(sam.sar)

#  # AR model
#  sam.ar <- lagsarlm(v2 ~ v1, listw=r_listw, data=r_poly,tol.solve=1e-36)
#  summary(sam.ar)

#}

#### PART 3: CREATE DATA FOR SAME DATE SPATIAL AUTO REG



### Prepare dataset 4 for function: date t+2 (pt2) and t+3 (pt3) before and after hurricane (t0)
#this is 152
r_var <- subset(r_stack,152)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("one_date_152",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_152 <-create_sp_poly_one_date_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

#this is 153
r_var <- subset(r_stack,153)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("one_date_153",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_153 <-create_sp_poly_one_date_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

#this is 154
r_var <- subset(r_stack,154)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("one_date_154",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_154 <-create_sp_poly_one_date_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

#this is 155
r_var <- subset(r_stack,155)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("one_date_155",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_155 <-create_sp_poly_one_date_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

#this is 156
r_var <- subset(r_stack,156)
r_clip <- raster_w
proj_str<- CRS_WGS84 
#out_dir and out_suffix set earlier
out_suffix_s <- paste("one_date_156",out_suffix,sep="_")
#undebug(create_sp_poly_spatial_reg)

nb_obj_for_pred_156 <-create_sp_poly_one_date_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)
