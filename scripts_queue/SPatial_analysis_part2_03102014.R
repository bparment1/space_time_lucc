####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbour to predict NDVI values in the MOORE EDGY region.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 03/10/2014
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

### Make a function to uniformize NA accros given dates!!!

create_uniform_mask_fun <- function(r_var,proj_str=NULL,r_clip=NULL){
  #Function to create uniform mask accross raster layer stack
  
  if(!is.null(r_clip)){
      r_var_w <- crop(r_var,r_clip) #clip (or "window") raster image to Moore subset area
      r_var_w <- mask(r_var_w,r_clip) #mask areas
  }else{
    r_var_w <- r_var
  }
  
  r_NA_mask <-calc(r_var_w,sum,na.rm=FALSE)
  if(!is.null(proj_str)){
    r_NA_mask <- projectRaster(r_NA_mask,crs=proj_str) #project to latlong
  }
  return(r_NA_mask)
}

predict_var_spatial_reg <-function(i,list_param){

  #Extract parameters/arguments
  test_shp_fname <- list_param$list_shp[i]
  out_dir  <- list_param$out_dir
  r_ref_s    <- list_param$r_ref #if NULL, no image is created
  proj_str <- list_param$proj_str
  list_models <- list_param$list_models
  out_suffix <- list_param$out_suffix[i]
  file_format <- list_param$file_format
  
  #### START SCRIPT
  
  list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!  
  formula <-list_formulas[[i]]
  
  r_ref <- subset(r_ref_s,i)
  test_shp <- readOGR(dsn=dirname(test_shp_fname),
                    layer=gsub(extension(basename(test_shp_fname)),"",basename(test_shp_fname)))

  test_df <- as.data.frame(test_shp)
  #lm_mod<-lm(v2~W_V1+v1,data=test_df) #tested model
  
  lm_mod <- try(lm(formula,data=test_df)) #tested model

  test_shp$pred_R <- lm_mod$fitted.values
  test_shp$pred_R_res <- lm_mod$residuals
  proj4string(test_shp) <- proj_str

  r_spat_pred <- rasterize(test_shp,r_ref,field="pred_R") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name <- paste("r_spat_pred","_",out_suffix,file_format,sep="")
  writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)
  #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
  r_spat_res <- rasterize(test_shp,r_ref,field="pred_R_res") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name2 <- paste("r_spat_res","_",out_suffix,file_format,sep="")
  writeRaster(r_spat_res,filename=file.path(out_dir,raster_name2),overwrite=TRUE)

  #Now write out gal and shp files

  out_fname <-file.path(out_dir,paste("test_shp_in_R","_",out_suffix,".shp",sep=""))  
  writeOGR(test_shp,dsn= dirname(out_fname),layer= gsub(".shp","",basename(out_fname)), 
         driver="ESRI Shapefile",overwrite_layer=TRUE)
  
  #could add mod objet
  spat_reg_obj <- list(lm_mod,out_fname,raster_name,raster_name2)
  names(spat_reg_obj) <- c("lm_mod","shp","raster_pred","raster_res")
  save(spat_reg_obj,file= file.path(out_dir,paste("spat_reg_obj","_",out_suffix,".RData",sep="")))

  return(spat_reg_obj)
}

#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"/home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"

#set up the working directory
setwd(in_dir)

moore_window <- file.path(in_dir,"moore_window.rst")
test_shp_path <- file.path(in_dir,"GEODA_Analysis") #path to get to the files for spatial pred
#  "/Users/benoitparmentier/Dropbox/Data/Space_Time/GEODA_Analysis/TEST_SPATIAL_ONE.shp"
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

out_suffix <-"03102014" #output suffix for the files that are masked for quality and for 

########### START SCRIPT #################

## Read in data
moore_w <- raster(moore_window) 
projection(moore_w) <- CRS_interp

#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
#reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013.rst
#reg_var_list <- list.files(pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")

r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
#levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)

projection(r_stack) <- CRS_interp #this is modis sinusoidal
r_var <- subset(r_stack,153:154) #date before hurricane is 153
r_clip <- moore_w
NAvalue(r_clip) <- 0

r_NA_mask_1 <- create_uniform_mask_fun(r_var=subset(r_stack,152:153),proj_str=CRS_WGS84,r_clip=r_clip)
r_NA_mask_2 <- create_uniform_mask_fun(r_var=subset(r_stack,153:154),proj_str=CRS_WGS84,r_clip=r_clip)
r_NA_mask_3 <- create_uniform_mask_fun(r_var=subset(r_stack,154:155),proj_str=CRS_WGS84,r_clip=r_clip)
r_NA_mask_4 <- create_uniform_mask_fun(r_var=subset(r_stack,155:156),proj_str=CRS_WGS84,r_clip=r_clip)

#r4 <- create_uniform_mask_fun(r_var,proj_str=NULL,r_clip=NULL)

list_shp <- list.files(path=test_shp_path,pattern="TEST_SPATIAL.*.shp",full.names=TRUE)
list_shp <- list_shp[c(2,4,3,1)]
r_ref <- stack(r_NA_mask_1,r_NA_mask_2,r_NA_mask_3,r_NA_mask_4)

#Add option to run other stuff than lm!!
list_models<-c("v2~W_V1+v1",
               "v2~W_V1+v1",  
               "v2~W_V1+v1",
               "v2~W_V1+v1")

## function to predict...
out_suffix_s <- paste(c("_1","_2","_3","_4"),out_suffix,sep="_")
list_param_spat_reg <-list(list_shp, out_dir,r_ref,proj_str,list_models,out_suffix_s,file_format)
names(list_param_spat_reg) <- c("list_shp", "out_dir","r_ref","proj_str","list_models","out_suffix","file_format")

#debug(predict_var_spatial_reg)
spat_reg_pred_obj1 <- predict_var_spatial_reg(1,list_param=list_param_spat_reg) 
spat_reg_pred_obj2 <- predict_var_spatial_reg(2,list_param=list_param_spat_reg) 
spat_reg_pred_obj3 <- predict_var_spatial_reg(3,list_param=list_param_spat_reg) 
spat_reg_pred_obj4 <- predict_var_spatial_reg(4,list_param=list_param_spat_reg) 

#list_spat_reg_obj <- lapply(1:length(list_models),FUN=predict_var_spatial,list_param=list_param_spat_reg)

r_pred_spat <- stack(list.files(path=out_dir,pattern=paste(out_suffix,file_format,sep="")))
plot(r_pred_s,colNA="black")

r_pred_spat <- subset(r_pred_spat,1:4)
#r_pred_s <- subset(r_pred_s,1:4)

##### END OF SCRIPT ########