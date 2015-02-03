####################################    Space Time Analyses Project   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script functions to produce predictions for the dates following the Hurricane Dean event.       
#The script uses spatial regression with weight matrix to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 09/25/2014
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Space beats time paper 
#################################################################################################

#This script currently contains 7 functions:

#[1] "calc_ac_stat_fun" : compute MAE for predictions wiht or without regions/zones
#[2] "create_dir_fun"   : create an output directory                    
#[3] "create_sp_poly_spatial_reg" : create a list of weights and polygon file for spatial regression
#[4] "create_uniform_mask_fun" : harmonize NA values given input layers with different valid values            
#[5] "load_obj" : load R object                          
#[6] "predict_spat_reg_fun" : function to perform spatial regresssion prediction
#[7] "predict_temp_reg_fun : function to preforma temperoral  prediction       

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
#library(forecast) #ARIMA forecasting
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
  #This might need to be changed!!!
  if(estimator=="gmm"){ #generalized method of moments: this is not available old packages...
    #lm_mod <- try(lm(formula,data=test_df)) #tested model
    v2 <- rnorm(nrow(data_reg))  #see if this work...
    data_reg$v2 <- v2 - mean(v2)
    spat_mod <- try(spreg(v1 ~ v2,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))
    #spat_mod <- try(spreg(v1 ~ v2,data=data_reg, listw= reg_listw_w, model="error",   
    #                 het = TRUE, verbose=TRUE))
  }
  #res<- spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
  #            het = TRUE, verbose=TRUE)
  #ordinary least square, this estimator should not be used for spatial reg but is  here for didactic purposes (course):
  #ie. values of standard errors can be compared with other...
  if(estimator=="ols"){  
    #lm_mod <- try(lm(formula,data=test_df)) #tested model
    v2 <- lag.listw(reg_listw_w,var=data_reg$v2) #creating a lag variable...#maybe this should be in the cleaning out function?
    data_reg$v3 <- data_reg$v2 #this is the lag variable...
    #This is the error!!! changed on 08/06/2014
    data_reg$v2 <- v2
    spat_mod <- try(lm(v1 ~ v2,data=data_reg))
  }
  
  #summary(sam.esar)
  #Predicted values and
  if(estimator!="gmm"){
    data_reg_spdf$spat_reg_pred <- spat_mod$fitted.values #should work for any of the three model
  }else{
    data_reg_spdf$spat_reg_pred <- spat_mod$residuals + data_reg_spdf$v1
  }
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

#Fuction to rasterize a table with coordinates and variables...,maybe add option for ref image??
#Make this more efficient!!!
rasterize_df_fun <- function(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format=".rst",NA_flag_val=-9999,tolerance_val= 0.000120005){
  data_spdf <- data_tb
  coordinates(data_spdf) <- cbind(data_spdf[[coord_names[1]]],data_spdf[[coord_names[2]]])
  proj4string(data_spdf) <- proj_str

  data_pix <- try(as(data_spdf,"SpatialPixelsDataFrame"))
  #tolerance_val <- 0.000120005 
  #tolerance_val <- 0.000856898
  if(inherits(data_pix,"try-error")){
      data_pix <- SpatialPixelsDataFrame(data_spdf, data=data_spdf@data, tolerance=tolerance_val) 
  }
  
  #test <- as(data_spdf,"SpatialPixelsDataFrame")

  # set up an 'empty' raster, here via an extent object derived from your data
  #e <- extent(s100[,1:2])
  #e <- e + 1000 # add this as all y's are the same

  #r <- raster(e, ncol=10, nrow=2)
  # or r <- raster(xmn=, xmx=,  ...

  data_grid <- as(data_pix,"SpatialGridDataFrame") #making it a regural grid
  r_ref <- raster(data_grid) #this is the ref image
  rast_list <- vector("list",length=ncol(data_tb))
  
  for(i in 1:(ncol(data_tb))){
    field_name <- names(data_tb)[i]
    var <- as.numeric(data_spdf[[field_name]])
    data_spdf$var  <- var
    #r <-rasterize(data_spdf,r_ref,field_name)
    r <-rasterize(data_spdf,r_ref,"var",NAflag=NA_flag_val,fun=mean) #prolem with NA in NDVI!!

    data_name<-paste("r_",field_name,sep="") #can add more later...
    #raster_name<-paste(data_name,out_names[j],".tif", sep="")
    raster_name<-paste(data_name,"_",out_suffix,file_format, sep="")
  
    writeRaster(r, NAflag=NA_flag_val,filename=file.path(out_dir,raster_name),overwrite=TRUE)
    #Writing the data in a raster file format...
    rast_list[i] <-file.path(out_dir,raster_name)
  }
  return(unlist(rast_list))
}

################### END OF SCRIPT ##################