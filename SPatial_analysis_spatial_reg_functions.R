####################################    Space Time Analyses Project   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script functions to produce predictions for the dates following the Hurricane Dean event.       
#The script uses spatial regression with weight matrix to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 05/208/2015
#Version: 2
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Space beats time paper
#PROJECT: Geocomputation and AAG 2015
#TO DO:
# modify the rasterize_df_fun function to allow ref image
# add the ARIMA method...in the temporal pred?
#
#################################################################################################

#This script currently contains 7 functions:

#[1] "calc_ac_stat_fun" : compute MAE for predictions wiht or without regions/zones
#[2] "create_dir_fun"   : create an output directory                    
#[3] "create_sp_poly_spatial_reg" : create a list of weights and polygon file for spatial regression
#[4] "create_uniform_mask_fun" : harmonize NA values given input layers with different valid values            
#[5] "load_obj" : load R object                          
#[6] "predict_spat_reg_fun" : function to perform spatial regresssion prediction
#[7] "predict_temp_reg_fun : function to preform temporal  predictions       
#[8] "rasterize_df_fun" : create raster from textfile with location information
#
#Add arima functions
#[9] "convert_arima_pred_to_raster"         
#[10] "extract_arima_mod_info"              
#[11] "pixel_ts_arima_predict"              
#[12] "raster_NA_image"                     
#[13] "raster_ts_arima"                      
#[14] "raster_ts_arima_predict"                         

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast) #ARIMA forecasting
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
library(sphet) #contains spreg
library(BMS) #contains hex2bin and bin2hex
library(bitops)

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

###ARIMA RELATED FUNCTIONS...need to be improved a lot!!

raster_ts_arima<-function(pixel,na.rm=T,arima_order){
  arima_obj<-arima(pixel,order=arima_order)
  a<-as.numeric(coef(arima_obj)[1]) 
  return(a)
}

#This takes a pixel time series ... extracted from a stack

raster_ts_arima_predict <- function(pixel,na.rm=T,arima_order=NULL,n_ahead=2){
  if(is.null(arima_order)){
    arima_mod <- auto.arima(pixel)
    p_arima<- try(predict(arima_mod,n.ahead=n_ahead))
  }else{
    arima_mod<-arima(pixel,order=arima_order)
    p_arima<- try(predict(arima_mod,n.ahead=n_ahead))
  }
  if (!inherits(p_arima,"try-error")){
    y<- t(as.data.frame(p_arima$pred)) #this makes a matrix...should probably be a data.frame
    y_error <- t(as.data.frame(p_arima$se)) #this contains the standard errrors related to the value predicted by arima
    #y_error <- rep(0,n_ahead)
  }
  if (inherits(p_arima,"try-error")){
    y<- rep(NA,n_ahead)
    y_error <- rep(1,n_ahead)
  }
  pred_obj <- list(y,y_error)
  names(pred_obj)<-c("pred","error")
                                  
  return(pred_obj)
}

#This is the main function!!!

pixel_ts_arima_predict <- function(i,list_param){
  #store coordinates???
  #change function to check for coordinates...it is better to use a data.frame...
  #V1 should be var and x, y can be coordinates...
  #extract parameters                                  
  pixel <-list_param$pix_val[i]
  arima_order <-list_param$arima_order
  n_ahead <- list_param$n_ahead
  out_dir <- list_param$out_dir
  out_suffix <- list_param$out_suffix
  pixel_xy <- try(list_param$df_xy[i,])
  na.rm=T
  
  #Start

  if(is.null(arima_order)){
    arima_mod <- auto.arima(pixel)
    p_arima<- try(predict(arima_mod,n.ahead=n_ahead))
  }else{
    arima_mod<-arima(pixel,order=arima_order)
    p_arima<- try(predict(arima_mod,n.ahead=n_ahead))
  }
  if (!inherits(p_arima,"try-error")){
    y<- t(as.data.frame(p_arima$pred)) #this makes a matrix...should probably be a data.frame
    y_error <- t(as.data.frame(p_arima$se)) #this contains the standard errrors related to the value predicted by arima
    #y_error <- rep(0,n_ahead)
  }
  if (inherits(p_arima,"try-error")){
    y<- rep(NA,n_ahead)
    y_error <- rep(1,n_ahead)
  }
  
  arima_mod_filename <- file.path(out_dir,paste("arima_mod_","pixel_",i,"_",out_suffix,".RData",sep=""))
  
  save(arima_mod,file=arima_mod_filename)             
  
  ## Prepare object to return...
  #note that the arima mod object could be included but is not at the time being
  #if pixel time series is geographically referenced
  #if(!is.null(pixel_xy)){
  if(!inherits(p_arima,"try-error")){
    pred_obj <- list(y,y_error,arima_mod_filename,pixel_xy)
    names(pred_obj)<-c("pred","error","arima_mod_filename","pixel_xy")
  }else{
    pred_obj <- list(y,y_error,arima_mod_filename)
    names(pred_obj)<-c("pred","error","arima_mod_filename")
  }
  
  return(pred_obj)
}

raster_NA_image <- function(r_stack){
  list_r_NA <- vector("list",length=nlayers(r_stack))
  for (i in 1:nlayers(r_stack)){
    r <- subset(r_stack,i)
    r_NA <- is.na(r)
    list_r_NA[[i]] <- r_NA
  }
  return(list_r_NA)
}

convert_arima_pred_to_raster <- function(i,list_param){
  #This function produces a raster image from ariam pred obj
  #Read in the parameters...
  r_ref <-list_param$r_ref
  ttx <- list_param$ttx
  file_format <- list_param$file_format
  out_dir <-list_param$out_dir
  out_suffix <- list_param$out_suffix
  out_rastname <-list_param$out_rastnames[i]
  file_format <- list_param$file_format
  NA_flag_val <- list_param$NA_flag_val
  
  #start script
  #pred_t <- lapply(ttx,FUN=function(x){x$pred[i]})
  #error_t <- lapply(ttx,FUN=function(x){x$error[i]})
  
  l_r <-vector("list", length=2)
  l_r[[1]]<-lapply(ttx,FUN=function(x){x$pred[i]})
  l_r[[2]] <- lapply(ttx,FUN=function(x){x$error[i]})
  
  for (j in 1:2){
    tt_dat <- do.call(rbind,l_r[[j]])
    tt_dat <- as.data.frame(tt_dat)
    pred_t <-as(r_ref,"SpatialPointsDataFrame")
    pred_t <- as.data.frame(pred_t)
    pred_t <- cbind(pred_t,tt_dat)
    coordinates(pred_t) <- cbind(pred_t$x,pred_t$y)
    raster_pred <- rasterize(pred_t,r_ref,"V1",fun=mean)
    l_r[[j]] <- raster_pred
  }
  
  #tmp_name <- extension(out_rastname)
  #modify output name to for error image
  tmp_name <- unlist(strsplit(out_rastname,"_"))
  nb<- length(tmp_name)
  tmp_name <-paste(paste(tmp_name[1:(nb-1)],collapse="_"),
             "error",tmp_name[nb],sep="_")
  writeRaster( l_r[[2]],NAflag=NA_flag_val,
              filename=file.path(out_dir,tmp_name),
              overwrite=TRUE)  
  writeRaster( l_r[[1]],NAflag=NA_flag_val,
              filename=file.path(out_dir,out_rastname),
              overwrite=TRUE)  
  return(list(out_rastname,tmp_name))
}

extract_arima_mod_info <- function(i,list_param){
  fname <- list_param$arima_mod_name[i]
  arima_mod <- load_obj(fname)
  #summary(arima_mod)
  #coef(arima_mod)
  arima_specification <- arima_mod$arma
  arima_coef <-  coef(arima_mod)
  #http://stackoverflow.com/questions/19483952/how-to-extract-integration-order-d-from-auto-arima
  #a$arma[length(a$arma)-1] is the order d
  #[1] 2 0 0 0 1 0 0
  #A compact form of the specification, as a vector giving the number of AR (1), MA (2), 
  #seasonal AR (3) and seasonal MA coefficients (4), 
  #plus the period (5) and the number of non-seasonal (6) and seasonal differences (7).
  
  return(list(arima_specification,arima_coef))
} 

#Still in process, apply function pixel by pixel in a block read from disk for stack or brick raster
#not used here because we loose the pixel id for the object for ARIMA, this needs to be work out at some point!!!

#    readBloackRaster <- function(r_var,r_mask,bs_read=NULL,out_rast_fname=NULL,pixelFun=NULL,list_param=NULL){
#      
#      colNo <- ncol(r_var) #number of column
#      rowNo <- nrow(r_var)
#      #out <- raster(x,1)
#      #bs <- blockSize(out)#
#      if(is.null(bs_read)){
#        bs <- blockSize(r_var)
#        bsno <- bs$n #number of lbocks to read
#      }
#      #bs <- blockSize(r_var)
#      if(!is.null(out_rast_fname)){
#        out_rast <- writeStart(out_rast_fname, out_rast, overwrite=TRUE)
#      }
#      
#      for (i in 1:bsno){
#        v <- getValuesBlock(r_var, 
#                           row=bs$row[i], #starting row
#                           nrows=bs$nrows, 
#                           col=1, 
#                           ncols=colNo)#, 
#                           #format='')
#        v_out <- v
#        if(!is.null(pixFun)){
#          v_out <- lapply(1:nrow(v),FUN=pixelFun,list_param=list_param) #assumes that this is one value per pixel
#          out_rast <- writeValues(out_rast, v_out, bs$row[i])
#
#          }
#        }  
#      return(v_out)
#    }


#Predict using the previous date and OLS
predict_temp_reg_fun <-function(i,list_param){
  #Extract parameters/arguments
  out_dir  <- list_param$out_dir
  r_ref_s    <- list_param$r_var #if NULL, no image is created, this is the reference image
  #list_param$ <- rast_ref
  r_clip     <- list_param$r_clip
  proj_str <- list_param$proj_str
  list_models <- list_param$list_models
  file_format <- list_param$file_format
  estimator <- list_param$estimator
  estimation_method <- list_param$estimation_method #currently used only for mle from errorsarlm
  NA_flag_val <- list_param$NA_flag_val
  
  if(estimator== "arima"){
    out_suffix <- list_param$out_suffix[1]
  }else{
    out_suffix <- list_param$out_suffix[1] #changed this for now...
  }
  #ARIMA specific
  num_cores <- list_param$num_cores #paraallelization in space.this should be done by row or til enot by pixel!!!!
  time_step <- list_param$time_step #this is the time step for which to start the arima model with
  n_pred_ahead <- list_param$n_pred_ahead
  r_stack <- list_param$r_stack
  arima_order <- list_param$arima_order
  
  #### START SCRIPT
  
  if(!is.null(list_models)){
    list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!  
    formula <-list_formulas[[i]]
  }
  #r_ref_s <- subset(r_ref_s,i)
  #out_dir and out_suffix set earlier
  
  if(estimation_method=="ols"){
    
    n_pred<- i+1
    r_var2 <- subset(r_ref_s,i:n_pred)
    r_ref_s <- crop(r_var2,r_clip)
    
    data_reg2_spdf <- as(r_ref_s,"SpatialPointsDataFrame")
    names(data_reg2_spdf) <- c("t1","t2")
    data_reg2 <- as.data.frame(data_reg2_spdf)
    data_reg2 <- na.omit(data_reg2) #remove NA...this reduces the number of observations
    temp_mod <-try(lm(t2 ~ t1, data=data_reg2))
    #summary(lm_mod)
  
    #Predicted values and
    data_reg2$temp_pred <- temp_mod$fitted.values
    data_reg2$temp_res <- temp_mod$residuals
  
    coordinates(data_reg2) <- c("x","y")
    proj4string(data_reg2) <- projection(r_ref_s)
  
    r_temp_pred <- rasterize(data_reg2,r_ref_s,field="temp_pred") #this is the prediction from lm model
    r_temp_res <- rasterize(data_reg2,r_ref_s,field="temp_res") #this is the prediction from lm model
    
    #can change later to have t_step
    raster_name_pred <- paste(paste("r_temp_pred","_",estimator,"_",estimation_method,"_",n_pred,sep=""),"_",out_suffix,file_format,sep="")
    raster_name_res <- paste(paste("r_temp_res","_",estimator,"_",estimation_method,"_",n_pred,sep=""),"_",out_suffix,file_format,sep="")

    
    #raster_name_pred <- paste("r_temp_pred","_",estimator,"_",estimation_method,"_",out_suffix,file_format,sep="")
    writeRaster(r_temp_pred,filename=file.path(out_dir,raster_name_pred),overwrite=TRUE)
  
    #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
    #file_format <- ".rst"
    #raster_name_res <- paste("r_temp_res","_",estimator,"_",estimation_method,"_",out_suffix,file_format,sep="")
    writeRaster(r_temp_res,filename=file.path(out_dir,raster_name_res),overwrite=TRUE)

  }
 
  if(estimation_method=="arima"){
    
    #This will be made a function callled from here!!!
    #also add line by line reading to avoid reading all file in memory
    
    #PARAM
    #r_stack
    #n_pred_ahead <- 4 #number of temporal ARIMA predictions ahead..
    #time_step
    # file_format <- ".rst"
    #NA_flag_val <- -9999
    #r_ref_s
    #r_clip
    #out_dir
    #out_suffix
    
    ### start of function
    
    if(!is.null(r_clip)){
      #r_stack2<-r_stack #this is should not be in memoor!!!
      r_stack <- crop(r_stack,r_clip)
      r_ref_s <- crop(r_ref_s,r_clip)
    }
    
    n_start <- c(time_step) +1
    n_end   <- c(time_step)+n_pred_ahead
    r_obs_s <- subset(r_stack,n_start:n_end) #stack of observed layers
    
    #r1 <- subset(r_obs_s,1)
    #xy <-coordinates(r_obs_s)  #get x and y projected coordinates...
    #CRS_interp<-proj4string(r1)
    r_x <-init(r_clip,v="x")
    r_y <-init(r_clip,v="y")
    names(r_x) <- "x"
    names(r_y) <- "y"
    r_stack <- stack(r_stack,r_x,r_y)
    r_stack <- mask(r_stack,r_clip)
    #rm(r1)

<<<<<<< HEAD
    #Very inefficient, will be changed to avoid reading in memory: problem to be sloved
    #readBlockRaster() see earlier
=======
    #Very inefficient, will be changed to avoid reading in memory
    
>>>>>>> 976db84ba0842bafd991171971cd76173e60cd19
    
    pix_val <- as(r_stack,"SpatialPointsDataFrame") #this will be changed later...to read line by line!!!!
    pix_val2 <- as.data.frame(pix_val)
    df_xy <- pix_val2[,c("x","y")]
    pix_val2 <-  pix_val2[,1:time_step] #152
    pix_val2 <- as.data.frame(t(as.matrix(pix_val2 )))#dim 152x26,616

    ### Should add a window option to subset the pixels time series
    #

    ## Now prepare to call pixel based predictions for arima:
    
    out_suffix_s <- paste("arima","_",out_suffix,sep="") #can modify name of output suffix
    out_dir_arima <- create_dir_fun(out_dir,out_suffix_s) #arima models will be stored here

    #pixel <-list_param$pix_val[i]
    #arima_order <-list_param$arima_order
    #n_ahead <- list_param$n_ahead
    #out_dir <- list_param$out_dir
    #out_suffix <- list_param$out_suffix
    #na.rm=T
    #list_param_predict_arima_2 <- list(pix_val=pix_val2,na.rm=T,arima_order=NULL,n_ahead=n_pred_ahead)
    #adde coordinates: df_xy
    list_param_predict_arima_2 <- list(pix_val=pix_val2,arima_order=arima_order,n_ahead=n_pred_ahead,out_dir=out_dir_arima,out_suffix=out_suffix,na.rm=T,df_xy=df_xy)

    #undebug(pixel_ts_arima_predict)
    #test_pix_obj <- pixel_ts_arima_predict(1,list_param=list_param_predict_arima_2)
    arima_pixel_pred_obj <- mclapply(1:length(pix_val2), FUN=pixel_ts_arima_predict,list_param=list_param_predict_arima_2,mc.preschedule=FALSE,mc.cores = num_cores) 
   
    #save this in a separate folder!!!
   
    save(arima_pixel_pred_obj,file=paste("arima_pixel_pred_obj","_",out_suffix,".RData",sep=""))

    #r_ref_s
   
    #### Now convert predictions to raster images
   
    #out_rastnames <- paste(paste("NDVI_pred_mooore_auto",1:n_pred_ahead,sep="_"),"_",out_suffix,file_format,sep="")
    #out_rastnames <- paste(paste("r_temp_pred_arima_",time_step:c(time_step+n_pred_ahead),sep="_"),"_",out_suffix,file_format,sep="")
    #out_rastnames <- paste(paste("r_temp_pred_arima_",1:n_pred_ahead,sep="_"),"_",out_suffix,file_format,sep="")
    #out_rastnames_res <- paste(paste("r_temp_res_arima_",1:n_pred_ahead,sep="_"),"_",out_suffix,file_format,sep="")
    
    #can change later to have t_step
    raster_name_pred <- paste(paste("r_temp_pred","_",estimator,"_",estimation_method,"_",1:n_pred_ahead,sep=""),"_",out_suffix,file_format,sep="")
    raster_name_res <- paste(paste("r_temp_res","_",estimator,"_",estimation_method,"_",1:n_pred_ahead,sep=""),"_",out_suffix,file_format,sep="")

    #list_param_arima_convert <- list(r_ref_s,arima_pixel_pred_obj,file_format,out_dir,raster_name_pred,file_format,NA_flag_val)
    list_param_arima_convert <- list(r_clip,arima_pixel_pred_obj,file_format,out_dir,raster_name_pred,file_format,NA_flag_val)
   
    names(list_param_arima_convert) <- c("r_ref","ttx","file_format","out_dir","out_rastnames","file_format","NA_flag_val")

    #debug(convert_arima_pred_to_raster)
    ## Convert predicted values to raster...
    #pred_t_l<-lapply(1:1,FUN=convert_arima_pred_to_raster,list_param=list_param_arima_convert) #,mc.preschedule=FALSE,mc.cores = num_cores)

    pred_t_l<-lapply(1:n_pred_ahead,FUN=convert_arima_pred_to_raster,list_param=list_param_arima_convert) #,mc.preschedule=FALSE,mc.cores = num_cores)
    #pred_t_l<-mclapply(1:n_pred_ahead,FUN=convert_arima_pred_to_raster,list_param=list_param_arima_convert,mc.preschedule=FALSE,mc.cores = num_cores)

    #arima_pixel_pred_obj <- mclapply(1:length(pix_val2), FUN=pixel_ts_arima_predict,list_param=list_param_predict_arima_2,mc.preschedule=FALSE,mc.cores = num_cores) 

    pred_t_l <-unlist(pred_t_l)
    r_temp_pred  <- stack(pred_t_l[-c(grep(pattern="error",pred_t_l))])
    r_temp_error <- stack(pred_t_l[c(grep(pattern="error",pred_t_l))])
    r_temp_res <- r_temp_pred - r_obs_s
   
    names(r_temp_res) <- raster_name_res
   
    #raster_name_pred <- out_rastnames
    #raster_name_res <-  out_rastnames_res
    #writeRaster(r_temp_res,bylayer=T)
    writeRaster(r_temp_res,filename=file.path(out_dir,raster_name_res),NAflag=NA_flag_val,bylayer=T,overwrite=TRUE)
   
    temp_mod <- NA #too many to store?
     
  }
  
   #### PREPARE OBJECT TO RETURN

  #Adding mod object, the ARIMA model will be different...function...most likely
  temp_reg_obj <- list(temp_mod,file.path(out_dir,raster_name_pred),file.path(out_dir,raster_name_res))
  names(temp_reg_obj) <- c("temp_mod","raster_pred","raster_res")
  save(temp_reg_obj,file= file.path(out_dir,paste("temp_reg_obj","_",estimator,"_",estimation_method,"_",out_suffix,".RData",sep="")))
  
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
  estimation_method <- list_param$estimation_method #currently used only for mle from errorsarlm
  
  #### START SCRIPT
  
  #Formula option not in full use yet...
  if(!is.null(list_models)){
    list_formulas<-lapply(list_models,as.formula,env=.GlobalEnv) #mulitple arguments passed to lapply!!  
    formula <-list_formulas[[i]]
  }
  #else set formula to default??
  
  r_ref_s <- subset(r_ref_s,i) #subset the relevant layer, i.e. date from i is the relevant timestep
  #out_dir and out_suffix set earlier
  
  nb_obj_for_pred_t <- create_sp_poly_spatial_reg(r_ref_s,r_clip,proj_str,out_suffix=out_suffix,out_dir)
  
  r_poly_name <- nb_obj_for_pred_t$r_poly_name
  reg_listw_w <- nb_obj_for_pred_t$r_listw
  
  data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                           layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))
  data_reg <- as.data.frame(data_reg_spdf)
  
  ## Add options later to choose the model type: lagsar,esar,spreg,lm etc.
  ##Can also add a loop to use all of them and output table???
  if(estimator=="mle"){ #maximum likelihood estimator is used with eigen method as default, if not other method is used if not null
    if(is.null(estimation_method)){
      estimation_method <- "eigen"
    }
    spat_mod <- try(errorsarlm(v1~ 1, listw=reg_listw_w, 
                               data=data_reg,method=estimation_method,na.action=na.omit,zero.policy=TRUE,
                               tol.solve=1e-36))
  }
  #This might need to be changed!!!
  if(estimator=="gmm"){ #generalized method of moments: this is not available old packages...
    if(is.null(estimation_method)){
      estimation_method <- "gs2slshac"
    }

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
    if(is.null(estimation_method)){
      estimation_method <- "ols"
    }

    #lm_mod <- try(lm(formula,data=test_df)) #tested model
    v2 <- lag.listw(reg_listw_w,var=data_reg$v2) #creating a lag variable...#maybe this should be in the cleaning out function?
    data_reg$v3 <- data_reg$v2 #this is the lag variable...
    #This is the error!!! changed on 08/06/2014
    data_reg$v2 <- v2
    spat_mod <- try(lm(v1 ~ v2,data=data_reg))
  }
  
  #### PREPARE OBJECT TO RETURN
  
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
  raster_name <- paste("r_spat_pred_",estimator,"_",estimation_method,"_",out_suffix,file_format,sep="")
  writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)
  
  #plot(r_spat_pred,subset(r_s,2)) #quick visualization...
  r_spat_res <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_res") #this is the prediction from lm model
  #file_format <- ".rst"
  raster_name2 <- paste("r_spat_res_",estimator,"_",estimation_method,"_",out_suffix,file_format,sep="")
  writeRaster(r_spat_res,filename=file.path(out_dir,raster_name2),overwrite=TRUE)
  
  #Return object contains model fitted, input and output rasters
  spat_reg_obj <- list(spat_mod,r_poly_name,reg_listw_w,raster_name,raster_name2)
  names(spat_reg_obj) <- c("spat_mod","r_poly_name","reg_listw_w","raster_pred","raster_res")
  save(spat_reg_obj,file= file.path(out_dir,paste("spat_reg_obj_",estimator,"_",estimation_method,"_",out_suffix,".RData",sep="")))
  
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
  #Purpose: Calculate accuracy statistics for given regions/zones of the study area
  #Statistics are MAE (Mean Absolute Error) and RMSE(Root Mean Square Error)
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
  sd_mae_zones_tb <- zonal(abs(r_res_s),r_zones,stat="sd") #absolute error
  #sd_mse_zones_tb <- zonal(r_res_s^2,r_zones,stat="sd") #mean square error
  
  rmse_zones_tb <- cbind(mse_zones_tb[,1],sqrt(mse_zones_tb[,2:dim(mse_zones_tb)[2]])) #root mean square error
  colnames(rmse_zones_tb)[1] <- c("zone")
  #Overall Accuracy/errors 
  
  mae_tb <- cellStats(abs(r_res_s),mean) #calculate MAE for layer stack
  rmse_tb <- sqrt(cellStats((r_res_s)^2,mean)) #calculate rmse overall
  #sd_rmse_zones_tb <- zonal(abs(r_res_s),r_zones,stat="sd") #absolute error

  
  #write out residuals rasters
  r_res_s <- writeRaster(r_res_s,filename=file.path(out_dir,"r_res_s.tif"),bylayer=TRUE,
                         suffix=paste(1:nlayers(r_res_s),out_suffix,sep="_"),overwrite=TRUE)
  ac_obj <- list(mae_tb,rmse_tb,mae_zones_tb,rmse_zones_tb,sd_mae_zones_tb)
  names(ac_obj) <- c("mae_tb","rmse_tb","mae_zones_tb","rmse_zones_tb","sd_mae_zones_tb")
  
  return(ac_obj)
}

#Fuction to rasterize a table with coordinates and variables...,maybe add option for ref image??
#Make this more efficient!!!
rasterize_df_fun <- function(data_tb,coord_names,proj_str,out_suffix,num_cores=0,rast_ref=NULL,out_dir=".",file_format=".rst",NA_flag_val=-9999,tolerance_val= 0.000120005){
  #This function creates a stack of raster images from a data frame table.
  #If not reference raster image is provided, the function generate its own raster grid.
  #The generated grid is a function of the distance between points and the tolrance index
  #
  #Inputs: 
  #1)data)tb: data.frame with data
  #2)coord_names: names of columns containing x,y coordinates
  #3)proj_str: to reproject if needed, if null no reprojection
  #4)out_suffix: suffix added to output file name
  #5)num_cores: number of cores to use, if zero the no parralization (not implemented yet )
  #6)rast_ref: 
  #4)out_dir: output directory path
  #8)file_format: raster format used in writing out files
  #9)NA_flag_val: values of NA
  #10)tolerance_val: spatial resolution tolerance used in the grid generation
  #
  #Outputs: list of raster files that are written out 
  #r_listw: list of weights (Queen)
  #r_poly_name: shapefile name screeened for NA and no neighbours features
  #r_nb_name: neighbour object
  #zero_nb: feature polygon with no neighbours that were removed

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
  
  #this should be parallelized!!!
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
  
  ##
  return(unlist(rast_list))
}

################### END OF SCRIPT ##################