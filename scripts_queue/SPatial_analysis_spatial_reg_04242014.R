####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 04/23/2014
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
#library(BMS) #contains hex2bin and bin2hex
#library(bitops)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
#library(forecast)
#library(xts)
#library(zoo)
#library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
#library(sphet)

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
  #names(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
  layerNames(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
  
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

predict_temp_reg_fun <-function(i,list_param){
  #Extract parameters/arguments
  #test_shp_fname <- list_param$list_shp[i]
  out_dir  <- list_param$out_dir
  r_ref_s    <- list_param$r_var #if NULL, no image is created
  r_clip     <- list_param$r_clip
  proj_str <- list_param$proj_str
  list_models <- list_param$list_models
  out_suffix <- list_param$out_suffix[i]
  file_format <- list_param$file_format
  
  #### START SCRIPT
  
  if(!is.null(list_models)){
    list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!  
    formula <-list_formulas[[i]]
  }
  #r_ref_s <- subset(r_ref_s,i)
  #out_dir and out_suffix set earlier
  n_pred<- i+1
  r_var2 <- subset(r_ref_s,i:n_pred)
  r_ref_s <- crop(r_var2,r_clip)
  
  data_reg2_spdf <- as(r_ref_s,"SpatialPointsDataFrame")
  names(data_reg2_spdf) <- c("t1","t2")
  data_reg2 <- as.data.frame(data_reg2_spdf)
  data_reg2 <- na.omit(data_reg2) #remove NA...this reduces the number of observations
  lm_mod <- lm(t2 ~ t1, data=data_reg2)
  #summary(lm_mod)
  
  #Predicted values and
  data_reg2$lm_temp_pred <- lm_mod$fitted.values
  data_reg2$lm_temp_res <- lm_mod$residuals
  
  coordinates(data_reg2) <- c("x","y")
  proj4string(data_reg2) <- proj_str
  
  r_temp_pred <- rasterize(data_reg2,rast_ref,field="lm_temp_pred") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name <- paste("r_temp_pred","_",out_suffix,file_format,sep="")
  writeRaster(r_temp_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)
  
  #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
  r_temp_res <- rasterize(data_reg2,r_ref,field="lm_temp_res") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name2 <- paste("r_temp_res","_",out_suffix,file_format,sep="")
  writeRaster(r_temp_res,filename=file.path(out_dir,raster_name2),overwrite=TRUE)
  
  #could add mod objet
  temp_reg_obj <- list(lm_mod,raster_name,raster_name2)
  names(temp_reg_obj) <- c("lm_mod","raster_pred","raster_res")
  save(temp_reg_obj,file= file.path(out_dir,paste("temp_reg_obj","_",out_suffix,".RData",sep="")))
  
  return(temp_reg_obj)
  
}

predict_spat_reg_fun <- function(i,list_param){
  
  #Extract parameters/arguments
  #test_shp_fname <- list_param$list_shp[i]
  out_dir  <- list_param$out_dir
  r_ref_s    <- list_param$r_var #if NULL, no image is created
  r_clip     <- list_param$r_clip
  proj_str <- list_param$proj_str
  list_models <- list_param$list_models
  out_suffix <- list_param$out_suffix[i]
  file_format <- list_param$file_format
  
  #### START SCRIPT
  
  if(!is.null(list_models)){
    list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!  
    formula <-list_formulas[[i]]
  }
  r_ref_s <- subset(r_ref_s,i)
  #out_dir and out_suffix set earlier
  
  nb_obj_for_pred_t <-create_sp_poly_spatial_reg(r_ref_s,r_clip,proj_str,out_suffix=out_suffix,out_dir)
  
  r_poly_name <- nb_obj_for_pred_t$r_poly_name
  reg_listw_w <- nb_obj_for_pred_t$r_listw
  
  data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                           layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))
  data_reg <- as.data.frame(data_reg_spdf)
  
  ## Add options later to choose the model type: lagsar,esar,spreg,lm etc.
  if(mle==TRUE){
    sam.esar <- try(errorsarlm(v1~ 1, listw=reg_listw_w, 
                               data=data_reg,na.action=na.omit,zero.policy=TRUE,
                               tol.solve=1e-36))
  }
  if(gmm==TRUE){
    #lm_mod <- try(lm(formula,data=test_df)) #tested model
    sam.esar<- spreg(v1 ~ v2,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE)
    
  }
  
  #summary(sam.esar)
  #Predicted values and
  data_reg_spdf$spat_reg_pred <- sam.esar$fitted.values
  data_reg_spdf$spat_reg_res <- sam.esar$residuals
  
  r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name <- paste("r_spat_pred","_",out_suffix,file_format,sep="")
  writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)
  
  #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
  r_spat_res <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_res") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name2 <- paste("r_spat_res","_",out_suffix,file_format,sep="")
  writeRaster(r_spat_res,filename=file.path(out_dir,raster_name2),overwrite=TRUE)
  
  #could add mod objet
  spat_reg_obj <- list(sam.esar,r_poly_name,reg_listw_w,raster_name,raster_name2)
  names(spat_reg_obj) <- c("esar_mod","r_poly_name","reg_listw_w","raster_pred","raster_res")
  save(spat_reg_obj,file= file.path(out_dir,paste("spat_reg_obj","_",out_suffix,".RData",sep="")))
  
  return(spat_reg_obj)
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

load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}


calc_ac_stat_fun <- function(r_pred_s,r_var_s,r_zones,file_format=".tif",out_suffix){
  #Parameters:
  #Input:
  #r_pred_s: raster stack of layers predictions (for example predicted NDVI)
  #r_var_s: raster stack of layers actual values (for example observed NDVI)
  #r_zones: raster defining zones of relevance to accuracy (e.g.hurricane winds zones)
  #
  
  ##Functions used
  rmse_fun <-function(x){sqrt(mean(x^2,na.rm=TRUE))}
  mae_fun <-function(x){mean(abs(x),na.rm=TRUE)}
  
  ##Start script
  
  #Accuracy/errors by zones
  r_res_s <- r_pred_s - r_var_s #residuals, stack of raster layers
  mse_zones_tb <- zonal(r_res_s^2,r_zones,fun="mean") #mean square error
  mae_zones_tb <- zonal(abs(r_res_s),r_zones,fun="mean") #absolute error
  rmse_zones_tb <- sqrt(mse_zones_tb) #root mean square error
  
  #Overall Accuracy/errors 
  
  mae_tb <- cellStats(abs(r_res_s),mean) #calculate MAE for layer stack
  rmse_tb <- sqrt(cellStats((r_res_s)^2,mean)) #calculate rmse overall
  #write out residuals rasters
  r_res_s <- writeRaster(r_res_s,filename=file.path(out_dir,"r_res_s.tif"),bylayer=TRUE,
                         suffix=paste(1:nlayers(r_res_s),out_suffix,sep="_"),overwrite=TRUE)
  ac_obj <- list(mae_tb,rmse_tb,mae_zones_tb,rmse_zones_tb)
  names(ac_obj) <- c("mae_tb","rmse_tb","mae_zones_tb","rmse_zones_tb")
  
  return(ac_obj)
}

#####  Parameters and argument set up ###########

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
#in_dir <-"/home/parmentier/Data/Space_Time"
in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"
out_dir <-  in_dir

moore_window <- file.path(in_dir,"R_data_geoprocessing_workshop_04232014","window_test4.rst")
winds_zones_fname <- file.path(in_dir,"R_data_geoprocessing_workshop_04232014","00_windzones_moore_sin.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

file_format <- ".rst"
NA_value <- "-9999"
out_suffix <-"_predictions_04232014" #output suffix for the files that are masked for quality and for 
create_out_dir_param=TRUE

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

## Read in data defining the area of study
moore_w <- raster(moore_window) #read raster file
moore_w <- moore_w > -999 #reclass image
rast_ref <- moore_w #create ref image
freq(rast_ref) #frequency of values in rast_ref
rast_ref[rast_ref==0] <- NA #assign NA for zero val

winds_r <- raster(winds_zones_fname) #read raster file

projection(winds_r) <- proj_modis_str #assign projection

#reproject data to latlong WGS83 (EPSG4326)
winds_wgs84 <- projectRaster(from=winds_r,crs=CRS_WGS84,method="ngb") #Check that it is using ngb

#reads input NDVI time series images
reg_var_list <- list.files(path=file.path(in_dir,"moore_NDVI_wgs84"),
                           pattern="moore_reg.*.MOD13A2_A.*04062014.*.rst$",full.name=TRUE)
#moore_reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif

#Create a stack of layers
r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)
plot(r_stack,y=1:2)

################ PART II : RUN SPATIAL REGRESSION ############

#Now mask and crop layer
r_var <- subset(r_stack,153:154) #date before hurricane is 153
r_clip <- rast_ref
r_var <- crop(r_var,rast_ref)

pix_id_r <- rast_ref
values(pix_id_r) <- 1:ncell(rast_ref)
pix_id_r <- mask(pix_id_r,rast_ref)

#pix_id_r <- mask(pix_id_r,rast_ref)
pix_id_poly <- rasterToPolygons(pix_id_r)# {raster}
reg_list_nb <-poly2nb(pix_id_poly,pix_id_poly$pix_id_r,queen=TRUE)
reg_listw_b <-nb2listw(reg_list_nb, style="B",zero.policy=TRUE)
reg_listw_w <- nb2listw(reg_list_nb, style="W",zero.policy=TRUE)

#EDGY_dat_spdf_04062014
# SAR model

data_reg <- as(r_var,"SpatialPointsDataFrame")
data_reg <- as.data.frame(data_reg) #errorsarlm needs a dataframe!!!
names(data_reg)[1:2] <- c("NDVI_153","NDVI_154")
  
NDVI_selected="NDVI_153"
data_reg$NDVI <- data_reg[[NDVI_selected]]

## NOW RUN THE SPATIAL MODEL...

sam.esar <- errorsarlm(NDVI~1, listw=reg_listw_w, 
                      data=data_reg,na.action=na.omit,zero.policy=TRUE,
                      tol.solve=1e-36)
#This does not work because of NA!!

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

data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                    layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))
data_reg <- as.data.frame(data_reg_spdf)

sam.esar <- errorsarlm(v1~ 1, listw=reg_listw_w, 
                       data=data_reg,na.action=na.omit,zero.policy=TRUE,
                       tol.solve=1e-36)
summary(sam.esar)
#Predicted values and
data_reg_spdf$spat_reg_pred <- sam.esar$fitted.values
data_reg_spdf$spat_reg_res <- sam.esar$residuals

################### PART III RUN TEMPORAL MODEL USING LM ########

#### NOW PREDICTION USING PREVIOUS TIME STEP 

r_var2 <- subset(r_stack,152:153)
r_var2_w <- crop(r_var2,rast_ref)

data_reg2_spdf <- as(r_var2_w,"SpatialPointsDataFrame")
names(data_reg2_spdf) <- c("t1","t2")
data_reg2 <- as.data.frame(data_reg2_spdf)
data_reg2 <- na.omit(data_reg2) #remove NA...this reduces the number of observations
lm_mod <- lm(t2 ~ t1, data=data_reg2)
summary(lm_mod)

#Predicted values and
data_reg2$lm_temp_pred <- lm_mod$fitted.values
data_reg2$lm_temp_res <- lm_mod$residuals

coordinates(data_reg2) <- c("x","y")
proj4string(data_reg2) <- CRS_WGS84

### NOW CREATE THE IMAGES BACK ...

#proj4string(test_shp) <- proj_str

r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
#file_format <- ".rst"
raster_name <- paste("r_spat_pred","_",out_suffix,file_format,sep="")
writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)

r_temp_pred <- rasterize(data_reg2,rast_ref,field="lm_temp_pred") #this is the prediction from lm model
#file_format <- ".rst"
raster_name <- paste("r_temp_pred","_",out_suffix,file_format,sep="")
writeRaster(r_temp_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)

## This can be repeated in a loop to produce predictions and compare the actual to predicted

##############################################################################################
############## PART IV PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER 4 time steps ####

r_spat_var <- subset(r_stack,152:155)
list_models <-NULL
out_suffix_s <- paste("t_",152:155,out_suffix,sep="")
list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format")
n_pred <- nlayers(r_spat_var)
#debug(predict_spat_reg_fun)
test_spat <- predict_spat_reg_fun(1,list_param_spat_reg)

pred_spat_l <-lapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg)

#pred_temp_l <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)

spat_pred_rast <- stack(lapply(pred_spat_l,FUN=function(x){x$raster_pred}))

## Predict using temporal info: time steps 153,154,155,156

r_temp_var <- subset(r_stack,151:155)
list_models <-NULL
out_suffix_s <- paste("t_",152:155,out_suffix,sep="")
list_param_temp_reg <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,file_format)
names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format")
n_pred <- nlayers(r_temp_var) -1
#undebug(predict_temp_reg_fun)
test_temp <- predict_temp_reg_fun(1,list_param_temp_reg)

pred_temp_l <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)

temp_pred_rast <- stack(lapply(pred_temp_l,FUN=function(x){x$raster_pred}))


############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

r_huric_w <- subset(r_stack,152:155)
r_huric_w <- crop(r_huric_w,rast_ref)

res_temp_s <- temp_pred_rast - r_huric_w
res_spat_s <- spat_pred_rast _ r_huric_w

out_suffix_s <- paste("temp_",out_suffix,sep="_")
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=r_pred_temp,r_var_s=r_huric_w,r_zones=r_winds_m,
                                file_format=file_format,out_suffix=out_suffix_s)
out_suffix_s <- paste("spat_",out_suffix,sep="_")  
ac_spat_obj <- calc_ac_stat_fun(r_pred_s=r_pred_spat,r_var_s=r_huric_w,r_zones=r_winds_m,
                                file_format=file_format,out_suffix=out_suffix_s)

##### END OF SCRIPT ########