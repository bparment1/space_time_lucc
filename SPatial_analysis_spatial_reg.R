####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 06/19/2014
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
library(sphet) #contains spreg

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

##This function needs to be modified...it is still using two dates format...!!
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
  names(r_s)<- c("v1","v2") #use layerNames(...) for earlier version of Raster
  
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
#Predict using the previous date and OLS
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
  proj4string(data_reg2) <- projection(r_ref_s)
  
  r_temp_pred <- rasterize(data_reg2,r_ref_s,field="lm_temp_pred") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name <- paste("r_temp_pred","_",out_suffix,file_format,sep="")
  writeRaster(r_temp_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)
  
  #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
  r_temp_res <- rasterize(data_reg2,r_ref_s,field="lm_temp_res") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name2 <- paste("r_temp_res","_",out_suffix,file_format,sep="")
  writeRaster(r_temp_res,filename=file.path(out_dir,raster_name2),overwrite=TRUE)
  
  #Adding mod object, the ARIMA model will be different...function...most likely
  temp_reg_obj <- list(lm_mod,raster_name,raster_name2)
  names(temp_reg_obj) <- c("lm_mod","raster_pred","raster_res")
  save(temp_reg_obj,file= file.path(out_dir,paste("temp_reg_obj","_",out_suffix,".RData",sep="")))
  
  return(temp_reg_obj)
  
}
#Predict using lag or error model...
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
  estimator <- list_param$estimator
    
  #### START SCRIPT
  
  if(!is.null(list_models)){
    list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!  
    formula <-list_formulas[[i]]
  }
  r_ref_s <- subset(r_ref_s,i) #subset the relevant layer
  #out_dir and out_suffix set earlier
  
  nb_obj_for_pred_t <- create_sp_poly_spatial_reg(r_ref_s,r_clip,proj_str,out_suffix=out_suffix,out_dir)
  
  r_poly_name <- nb_obj_for_pred_t$r_poly_name
  reg_listw_w <- nb_obj_for_pred_t$r_listw
  
  data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                           layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))
  data_reg <- as.data.frame(data_reg_spdf)
  
  ## Add options later to choose the model type: lagsar,esar,spreg,lm etc.
  ##Can also add a loop to use all of them and output table???
  if(estimator=="mle"){ #maximum likelihood estimator is used for the 
    spat_mod <- try(errorsarlm(v1~ 1, listw=reg_listw_w, 
                               data=data_reg,na.action=na.omit,zero.policy=TRUE,
                               tol.solve=1e-36))
  }
  if(estimator=="gmm"){ #generalized method of moments: this is not available old packages...
    #lm_mod <- try(lm(formula,data=test_df)) #tested model
    spat_mod <- try(spreg(v1 ~ v2,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))
  }
  #res<- spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
  #            het = TRUE, verbose=TRUE)
  #ordinary least square, this estimator should not be used for spatial reg but is  here for didactic purposes (course):
  #ie. values of standard errors can be compared with other...
  if(estimator=="ols"){  
    #lm_mod <- try(lm(formula,data=test_df)) #tested model
    v2 <- lag(reg_listw_w,var=data_reg$v2) #creating a lag variable...#maybe this should be in the cleaning out function?
    data_reg$v3 <- data_reg$v2 #this is the lag variable...
    data_reg$v3 <- v2
    spat_mod <- try(lm(v1 ~ v2,data=data_reg))
  }
  
  #summary(sam.esar)
  #Predicted values and
  data_reg_spdf$spat_reg_pred <- spat_mod$fitted.values #should work for any of the three model
  data_reg_spdf$spat_reg_res <- spat_mod$residuals
  
  r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name <- paste("r_spat_pred_",estimator,"_",out_suffix,file_format,sep="")
  writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)
  
  #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
  r_spat_res <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_res") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name2 <- paste("r_spat_res_",estimator,"_",out_suffix,file_format,sep="")
  writeRaster(r_spat_res,filename=file.path(out_dir,raster_name2),overwrite=TRUE)
  
  #Return object contains model fitted, input and output rasters
  spat_reg_obj <- list(spat_mod,r_poly_name,reg_listw_w,raster_name,raster_name2)
  names(spat_reg_obj) <- c("spat_mod","r_poly_name","reg_listw_w","raster_pred","raster_res")
  save(spat_reg_obj,file= file.path(out_dir,paste("spat_reg_obj_",estimator,"_",out_suffix,".RData",sep="")))
  
  return(spat_reg_obj)
}

create_dir_fun <- function(out_dir,out_suffix){
  #if out_suffix is not null then append out_suffix string
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

rasterize_df_fun <- function(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format=".rst",NA_flag_val=-9999){
  data_spdf <- data_tb
  coordinates(data_spdf) <- cbind(data_spdf[[coord_names[1]]],data_spdf[[coord_names[2]]])
  proj4string(data_spdf) <- proj_str

  data_pix <- as(data_spdf,"SpatialPixelsDataFrame")
  data_grid <- as(d,"SpatialGridDataFrame") #making it a regural grid
  r_ref <- raster(data_grid) #this is the ref image
  rast_list <- vector("list",length=ncol(data_tb))
  
  for(i in 1:(ncol(data_tb))){
    field_name <- names(data_tb)[i]
    r <-rasterize(data_spdf,r_ref,field_name)
    data_name<-paste("r_",field_name,sep="") #can add more later...
    #raster_name<-paste(data_name,out_names[j],".tif", sep="")
    raster_name<-paste(data_name,out_suffix,file_format, sep="")
  
    writeRaster(r, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name),overwrite=TRUE)
    #Writing the data in a raster file format...
    rast_list[i] <-file.path(out_dir,raster_name)
  }
  return(unlist(rast_list))
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
  mse_zones_tb <- zonal(r_res_s^2,r_zones,stat="mean") #mean square error
  mae_zones_tb <- zonal(abs(r_res_s),r_zones,stat="mean") #absolute error
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

in_dir<-"~/Data/Space_beats_time/R_Workshop_April2014"
#in_dir <-"/home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Google Drive/R_Workshop_April2014"
#out_dir <-  "/Users/benoitparmentier/Google Drive/R_Workshop_April2014"

moore_window <- file.path(in_dir,"window_test4.rst")
winds_zones_fname <- file.path(in_dir,"00_windzones_moore_sin.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

file_format <- ".rst"
NA_value <- -9999
NA_flag_val <- NA_value
out_suffix <-"_predictions_06192014" #output suffix for the files that are masked for quality and for 
create_out_dir_param=TRUE

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_tb <-read.table(data_fname,sep=",")

#### Make this a function...

coord_names <- c("XCoord","YCoord")
l_rast <- rasterize_df_fun(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format,NA_flag_val)
#debug(rasterize_df_fun)

reg_var_list <- l_rast[6:18]
r_stack <- stack(reg_var_list)
plot(1:13,data_tb[600,6:18])

r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
#reads input NDVI time series images
reg_var_list <- list.files(path=file.path(in_dir,"moore_NDVI_wgs84"),
                           pattern="moore_reg.*.MOD13A2_A.*04062014.*.rst$",full.name=TRUE)
#moore_reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif

#Create a stack of layers
#r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
levelplot(r_stack,layers=5:6) #show first 2 images (half a year more or less)
plot(r_stack,y=5:6)

################ PART II : RUN SPATIAL REGRESSION ############

#Now mask and crop layer
r_var <- subset(r_stack,5:6) #date before hurricane is 153
#r_clip <- rast_ref
#r_var <- crop(r_var,rast_ref) #crop/clip image using another reference image
rast_ref <- subset(r_stack,1)
rast_ref <- rast_ref != NA_flag_val
pix_id_r <- rast_ref
values(pix_id_r) <- 1:ncell(rast_ref) #create an image with pixel id for every observation
pix_id_r <- mask(pix_id_r,rast_ref) #3854 pixels from which 779 pixels are NA

#pix_id_r <- mask(pix_id_r,rast_ref)
pix_id_poly <- rasterToPolygons(pix_id_r) #convert raster to polygon (raster package option)
reg_list_nb <-poly2nb(pix_id_poly,pix_id_poly$pix_id_r,queen=TRUE) #list of neighbours generated from poly
reg_listw_b <-nb2listw(reg_list_nb, style="B",zero.policy=TRUE) #weight list using binary style
reg_listw_w <- nb2listw(reg_list_nb, style="W",zero.policy=TRUE) #weight list using row centering
#Note that zero.policy allows for tracking features with zero neighbours
#EDGY_dat_spdf_04062014
# SAR model

data_reg <- as(r_var,"SpatialPointsDataFrame")
data_reg <- as.data.frame(data_reg) #errorsarlm needs a dataframe!!!
names(data_reg)[1:2] <- c("pop_2004","pop_2005")
  
var_selected="pop_2005" #this is the variable modeled...
data_reg$var <- data_reg[[var_selected]]

## NOW RUN THE SPATIAL MODEL...

sam.esar <- errorsarlm(var~1, listw=reg_listw_w, 
                      data=data_reg,na.action=na.omit,zero.policy=TRUE,
                      tol.solve=1e-36)
#This does not work because some pixels have neighbours with NA!!

### CLEAN OUT AND SCREEN NA and list of neighbours
#Let's use our function we created to clean out neighbours
#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
r_var <- subset(r_stack,6) #this is the date we want to use to run the spatial regression
r_clip <- rast_ref #this is the image defining the study area
proj_str<- NULL #SRS/CRS projection system
out_suffix_s <- paste("d2005",out_suffix,sep="_")
#out_dir and out_suffix set earlier

nb_obj_for_pred_t_2005 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)
names(nb_obj_for_pred_t_2005) #show the structure of object (which is made up of a list)
r_poly_name <- nb_obj_for_pred_t_2005$r_poly_name #name of the shapefile that has been cleaned out
reg_listw_w <- nb_obj_for_pred_t_2005$r_listw #list of weights for cleaned out shapefile
#Use OGR to load the screened out data: note that we have now 2858 features with 3 fields
data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                    layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))
data_reg <- as.data.frame(data_reg_spdf) #convert spdf to df

#Now run the spatial regression, since there are no covariates, the error and lag models are equivalent
#This takes a bit less than 3minutes for the dataset containing 2858 polygons
sam.esar <- errorsarlm(v1~ 1, listw=reg_listw_w, 
                       data=data_reg,na.action=na.omit,zero.policy=TRUE,
                       tol.solve=1e-36) #tol.solve use in matrix operations
summary(sam.esar)
#Predicted values and
data_reg_spdf$spat_reg_pred <- sam.esar$fitted.values
data_reg_spdf$spat_reg_res <- sam.esar$residuals

################### PART III RUN TEMPORAL MODEL USING LM ########

#### NOW PREDICTION USING PREVIOUS TIME STEP 

r_var2 <- subset(r_stack,5:6)
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

###########################################
############## PART IV: Produce images from the individual predictions using time and space ####

### NOW CREATE THE IMAGES BACK ...

r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
#file_format <- ".rst"
raster_name <- paste("r_spat_pred","_",out_suffix,file_format,sep="")
writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)

r_temp_pred <- rasterize(data_reg2,rast_ref,field="lm_temp_pred") #this is the prediction from lm model
#file_format <- ".rst"
raster_name <- paste("r_temp_pred","_",out_suffix,file_format,sep="")
writeRaster(r_temp_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)

## This can be repeated in a loop to produce predictions and compare the actual to predicted
r_pred <- stack(r_spat_pred,r_temp_pred) #stack of predictions
names(r_pred) <- c("spatial pred","temporal pred") #change layerNames to names when using R 3.0 or above
plot(stack(r_spat_pred,r_temp_pred))
levelplot(r_pred,regions.col=rev(terrain.colors(255)),main="Var predictions after hurricane")

#### Examining difference between predictions
r_dif <- r_spat_pred - r_temp_pred
plot(r_dif)
hist(r_dif)

##############################################################################################
############## PART V PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER 4 time steps ####

##This times will we use an automated function to generate predictions over 4 dates

### Predict using spatial regression
r_spat_var <- subset(r_stack,2:13) #predict before (step 152) and three dates after (step 153)
list_models <- NULL
proj_str <- NULL #if null the raster images are not reprojected
out_suffix_s <- paste("t_",2001:2012,out_suffix,sep="") #use mle estimation, name "mle" is included in the output automatically?
estimator <- "mle"

list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator")
n_pred <- nlayers(r_spat_var)
## First predict for one date for testing and comparison
#debug(predict_spat_reg_fun)
#use mle estimator
testd1_spat_mle <- predict_spat_reg_fun(1,list_param_spat_reg)
#use generalized methods of moment as estimator: won't work with version of 2.14 and old sphet package from Feb 2012
list_param_spat_reg$estimator <- "gmm"
#use ols as estimator...this is added as a dictatic form for teaching, not to be used for research
testd1_spat_gmm <- predict_spat_reg_fun(1,list_param_spat_reg)
list_param_spat_reg$estimator <- "ols"
testd1_spat_ols <- predict_spat_reg_fun(1,list_param_spat_reg)

## Now predict for four dates using "mle"
list_param_spat_reg$estimator <- "mle"
pred_spat_mle <-lapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg)
#Use parallel processing on MAC and Linux/Unix systems
#pred_spat_mle <-mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,,mc.preschedule=FALSE,mc.cores = 2)
#list_param_spat_reg$estimator <- "gmm"
#pred_spat_gmm <-lapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg)

spat_pred_rast <- stack(lapply(pred_spat_mle,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast <- stack(lapply(pred_spat_mle,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast) #view the four predictions using mle spatial reg.

## Predict using temporal info: time steps 2005..2012

r_temp_var <- subset(r_stack,1:13) #need date 151 because it relies on the previous date in contrast to spat reg
list_models <-NULL
out_suffix_s <- paste("t_",2001:2012,out_suffix,sep="")
list_param_temp_reg <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,file_format)
names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format")
n_pred <- nlayers(r_temp_var) -1
#debug(predict_temp_reg_fun)
#test_temp <- predict_temp_reg_fun(1,list_param_temp_reg)

pred_temp_lm <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)

temp_pred_rast <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
temp_res_rast <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
levelplot(temp_pred_rast) #view the four predictions using mle spatial reg.
projection(temp_pred_rast) <- CRS_WGS84
projection(spat_pred_rast) <- CRS_WGS84

############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

r_huric_w <- subset(r_stack,2:13)
r_huric_w <- crop(r_huric_w,rast_ref)

#reproject data to latlong WGS84 (EPSG4326)
#r_winds_m <- projectRaster(from=r_winds,res_temp_s,method="ngb") #Check that it is using ngb

#r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
res_temp_s <- temp_pred_rast - r_huric_w
res_spat_s <- spat_pred_rast - r_huric_w

out_suffix_s <- paste("temp_",out_suffix,sep="_")
#debug(calc_ac_stat_fun)
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast,r_var_s=r_huric_w,r_zones=rast_ref,
                                file_format=file_format,out_suffix=out_suffix_s)
out_suffix_s <- paste("spat_",out_suffix,sep="_")  
ac_spat_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,r_var_s=r_huric_w,r_zones=rast_ref,
                                file_format=file_format,out_suffix=out_suffix_s)

#mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg","temp")
mae_tot_tb$time <- 1:12

plot(spat_reg ~ time, type="b",col="magenta",data=mae_tot_tb,ylim=c(0,1000))
lines(temp ~ time, type="b",col="cyan",data=mae_tot_tb)
write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
legend("topleft",legend=c("spat","temp"),col=c("magenta","cyan"),lty=1)
title("Overall MAE for spatial and temporal models")

#### BY ZONES ASSESSMENT

#mae_zones_tb <- rbind(ac_spat_obj$mae_zones_tb[2:3,],
#                      ac_temp_obj$mae_zones_tb[2:3,])
#mae_zones_tb <- as.data.frame(mae_zones_tb)
#mae_zones_tb$method <- c("spat_reg","spat_reg","temp","temp")
#names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")

#write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

#Very quick and dirty plot
#time <-1:4
#x <- as.numeric(mae_zones_tb[1,2:5])
#plot(x~time,, type="b",col="magenta",lty=1,ylim=c(400,2000),ylab="MAE for NDVI")
#x <- as.numeric(mae_zones_tb[2,2:5])
#lines(x~time, type="b",lty=2,col="magenta")
#add temporal
x <- as.numeric(mae_zones_tb[3,2:5]) #zone 4
lines(x~time,, type="b",col="cyan",lty=1,ylim=c(400,2000))
x <- as.numeric(mae_zones_tb[4,2:5]) #zone 5
lines(x~time, type="b",lty=2,col="cyan")
legend("topleft",legend=c("spat zone 4","spat zone 5","temp zone 4","temp zone 5"),
        col=c("magenta","magenta","cyan","cyan"),lty=c(1,2,1,2))
title("MAE per wind zones for spatial and temporal models")

### more advanced plot to fix later....
#mae_val <- (as.vector(as.matrix(mae_zones_tb[,2:5])))
#avg_ac_tb <- as.data.frame(mae_val)

# avg_ac_tb$metric <- rep("mae",length(mae_val))
# avg_ac_tb$zones <- rep(c(3,4,5),4)
# avg_ac_tb$time <- c(rep(1,6),rep(2,6),rep(3,6),rep(4,6))
# avg_ac_tb$method <- rep(c("spat_reg","spat_reg","spat_reg","arima","arima","arima"),4)
# names(avg_ac_tb)[1]<- "ac_val"
# names_panel_plot <- c("time -1","time +1","time +2","time +3")
# p <- xyplot(ac_val~zones|time, # |set up pannels using method_interp
#             group=method,data=avg_ac_tb, #group by model (covariates)
#             main="Average MAE by winds zones and time step ",
#             type="b",as.table=TRUE,
#             #strip = strip.custom(factor.levels=c("time -1","time +1","time +2","time +3")),
#             strip=strip.custom(factor.levels=names_panel_plot),
#             index.cond=list(c(1,2,3,4)), #this provides the order of the panels)
#             pch=1:length(avg_ac_tb$method),
#             par.settings=list(superpose.symbol = list(
#               pch=1:length(avg_ac_tb$method))),
#             auto.key=list(columns=1,space="right",title="Model",cex=1),
#             #auto.key=list(columns=5),
#             xlab="Winds zones",
#             ylab="MAE for NDVI")
# print(p)

################### END OF SCRIPT ##################