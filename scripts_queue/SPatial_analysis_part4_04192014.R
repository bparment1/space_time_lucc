####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbour to predict NDVI values in the MOORE EDGY region.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 04/19/2014
#Version: 3
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
library(sphet)

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

create_dir_fun <- function(out_dir,out_suffix){
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    out_dir <- file.path(out_dir,out_name)
  }
  #create if does not exists
  if(!file.exists(out_dir)){
    dir.create(out_dir)
  }
  return(out_dir)
}

#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"/home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"
out_dir <- in_dir

input_data <- file.path(in_dir,"EDGY_dat_spdf_04072014.txt")

moore_window <- file.path(in_dir,"moore_window.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

file_format <- ".tif"
NA_value <- "-9999"
out_suffix <-"_predictions_04192014" #output suffix for the files that are masked for quality and for 
create_out_dir_param=TRUE

########### START SCRIPT #################

#set up the working directory
#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

## Read in data
moore_w <- raster(moore_window)
projection(moore_w) <- proj_modis_str
rast_ref <- projectRaster(from=moore_w,crs=CRS_WGS84,method="ngb") #Check that it is using ngb

#reads input NDVI time series images
reg_var_list <- list.files(path=file.path(in_dir,"moore_NDVI_wgs84"),
                           pattern="moore_reg.*.MOD13A2_A.*04072014.*.tif$",full.name=TRUE)
#moore_reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif

r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
#levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)

r_var <- subset(r_stack,153:154) #date before hurricane is 153
r_clip <- rast_ref

s_dat_var <- read.table(input_data,header=T,row.names=1,sep=",") #first  columns = row.names
#s_dat_var_poly <- readOGR(dsn=dirname(input_data),
#                          layer=exentension(basename(input_data)))
coordinates(s_dat_var) <- s_dat_var[,c("r_x","r_y")]
proj4string(s_dat_var) <- CRS_WGS84
#s_dat_var_poly <- readOGR(dsn=dirname(input_data),
#                    layer=gsub(extension(basename(input_data)),"",basename(input_data)))

pix_id_r <- rasterize(s_dat_var,rast_ref,"pix_id_r")

#pix_id_r <- mask(pix_id_r,rast_ref)
pix_id_poly <- rasterToPolygons(pix_id_r)# {raster}
reg_list_nb <-poly2nb(pix_id_poly,pix_id_poly$pix_id_r,queen=TRUE)
reg_listw_b <-nb2listw(reg_list_nb, style="B",zero.policy=TRUE)
reg_listw_w <- nb2listw(reg_list_nb, style="W",zero.policy=TRUE)

#EDGY_dat_spdf_04062014
# SAR model
data_reg <- as.data.frame(s_dat_var) #errorsarlm needs a dataframe!!!
NDVI_selected="NDVI_153"
data_reg$NDVI <- data_reg[[NDVI_selected]]

sam.esar <- errorsarlm(NDVI~1, listw=reg_listw_w, 
                      data=data_reg,na.action=na.omit,zero.policy=TRUE,
                      tol.solve=1e-36)
sam.esar <- errorsarlm(NDVI~ 1, listw=reg_listw_w, 
                       data=data_reg,na.action=na.omit,zero.policy=TRUE,
                       tol.solve=1e-36)

#summary(sam.sar) had an error message
sam.lsar <- lagsarlm(NDVI~ NDVI, listw=reg_listw_w, 
                      data=data_reg,na.action=na.omit,zero.policy=TRUE,
                      tol.solve=1e-36)

#summary(esar1f)
system.time(esar1f <- spautolm(NDVI~ NDVI,listw=reg_listw_w, 
                               data=data_reg, family="SAR", method="eigen", na.action=na.omit,zero.policy=TRUE,
                               tol.solve=1e-36,verbose=TRUE))
summary(esar1f)
#r_NA_mask_1 <- create_uniform_mask_fun(r_var=subset(r_stack,152),proj_str=CRS_WGS84,r_clip=r_clip)
#r_NA_mask_2 <- create_uniform_mask_fun(r_var=subset(r_stack,153:154),proj_str=CRS_WGS84,r_clip=r_clip)
#r4 <- create_uniform_mask_fun(r_var,proj_str=NULL,r_clip=NULL)


### CLEAN OUT AND SCREEN NA and list of neighbours

#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
r_var <- subset(r_stack,153)
r_clip <- rast_ref
proj_str<- CRS_WGS84 
out_suffix_s <- paste("d153",out_suffix,sep="_")
#out_dir and out_suffix set earlier

nb_obj_for_pred_t_153 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)

r_poly_name <- nb_obj_for_pred_t_153$r_poly_name
reg_listw_w <- nb_obj_for_pred_t_153$r_listw

data_reg <- readOGR(dsn=dirname(r_poly_name),
                    layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))
data_reg <- as.data.frame(data_reg)

res<- spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
            het = TRUE, verbose=TRUE)

res<- spreg(v1 ~ v2,data=data_reg, listw= reg_listw_w, model="error",   
            het = TRUE, verbose=TRUE)

res<- spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
            het = TRUE, verbose=TRUE)

res<- spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
            het = FALSE, verbose=TRUE)




# #Run regression with spatial lag variable...
# list_models<-c("v2~ W_V1",
#                "v2~ W_V1",  
#                "v2~ W_V1",
#                "v2~ W_V1")
# 
# ## function to predict...
# out_suffix_s <- paste(c("_spat_lag_1","_spat_lag_2","_spat_lag_3","_spat_lag_4"),out_suffix,sep="_")
# list_param_spat_reg <-list(list_shp, out_dir,r_ref,proj_str,list_models,out_suffix_s,file_format)
# names(list_param_spat_reg) <- c("list_shp", "out_dir","r_ref","proj_str","list_models","out_suffix","file_format")
# 
# #debug(predict_var_spatial_reg)
# #spat_reg_pred_obj1 <- predict_var_spatial_reg(1,list_param=list_param_spat_reg) 
# #spat_reg_pred_obj2 <- predict_var_spatial_reg(2,list_param=list_param_spat_reg) 
# #spat_reg_pred_obj3 <- predict_var_spatial_reg(3,list_param=list_param_spat_reg) 
# #spat_reg_pred_obj4 <- predict_var_spatial_reg(4,list_param=list_param_spat_reg) 
# 
# #list_spat_reg_obj <- lapply(1:length(list_models),FUN=predict_var_spatial_reg,list_param=list_param_spat_reg)
# list_spat_reg_obj <- mclapply(1:length(list_models),FUN=predict_var_spatial_reg,list_param=list_param_spat_reg,mc.preschedule=FALSE,mc.cores = 4)
# 
# r_pred_spat <- stack(list.files(path=out_dir,pattern=paste(out_suffix_s,".*",file_format,sep="")))
# plot(r_pred_spat,colNA="black")
# 
# r_pred_spat <- subset(r_pred_spat,1:4)
# #r_pred_s <- subset(r_pred_s,1:4)

##### END OF SCRIPT ########