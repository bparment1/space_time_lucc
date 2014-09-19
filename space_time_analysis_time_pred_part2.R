####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: ARIMA          #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#Predictions are done with ARIMA model leveraging temporal correlation.                        
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/27/2013 
#DATE MODIFIED: 04/21/2014
#Version: 4
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones             
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
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

###### Functions used in this script

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
    y<- t(as.data.frame(p_arima$pred))
    y_error <- rep(0,n_ahead)
  }
  if (inherits(p_arima,"try-error")){
    y<- rep(NA,n_ahead)
    y_error <- rep(1,n_ahead)
  }
  pred_obj <- list(y,y_error)
  names(pred_obj)<-c("pred","error")
                                  
  return(pred_obj)
}

raster_ts_arima_predict <- function(i,list_param){
  
  #extract parameters                                  
  pixel <-list_param$pix_val[i]
  arima_order <-list_param$arima_order
  n_ahead <- list_param$n_ahead
  out_suffix <- list_param$out_suffix
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
    y<- t(as.data.frame(p_arima$pred))
    y_error <- rep(0,n_ahead)
  }
  if (inherits(p_arima,"try-error")){
    y<- rep(NA,n_ahead)
    y_error <- rep(1,n_ahead)
  }
  pred_obj <- list(y,y_error)
  names(pred_obj)<-c("pred","error")
  save(arima_mod,file=paste("arima_mod_",i,out_suffix,".RData",sep=""))             
  
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

################### Parameters and arguments #################

## Location on Benoit's laptop (Dropbox)
#in_dir<- "/Users/benoitparmentier/Dropbox/Data/Space_Time"
in_dir<- "/home/parmentier/Data/Space_Time/"
#out_dir<- "/Users/benoitparmentier/Dropbox/Data/Space_Time"
out_dir<- "/home/parmentier/Data/Space_Time"
out_suffix <-"04212014" #output suffix for the files that are masked for quality and for 

moore_window <- file.path(in_dir,"moore_window.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

function_analyses_paper <- "MODIS_and_raster_processing_functions_04172014.R"
#script_path <- "~/Dropbox/Data/NCEAS/git_space_time_lucc/scripts_queue" #path to script functions
script_path <- file.path(in_dir,"R") #path to script functions

source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.
#This is the shape file of outline of the study area    
ref_rast_name<-file.path(in_dir,"/reg_input_yucatan/gyrs_sin_mask_1km_windowed.rst")  #local raster name defining resolution, exent: oregon

infile_modis_grid <- file.path(in_dir,"/reg_input_yucatan/modis_sinusoidal_grid_world.shp")

create_out_dir_param <- TRUE #create output directory using previously set out_dir and/or out_suffix
#It is an input/output of the covariate script
infile_reg_outline <- "~/Data/Space_Time/GYRS_MX_trisate_sin_windowed.shp"  #input region outline defined by polygon: Oregon
## Other specific parameters
NA_flag_val<- -9999
file_format <- ".rst" #problem wiht writing IDRISI for hte time being

############################  START SCRIPT ###################

## Read in data
moore_w <- raster(moore_window)
projection(moore_w) <- proj_modis_str
rast_ref <- projectRaster(from=moore_w,crs=CRS_WGS84,method="ngb") #Check that it is using ngb

#reads input NDVI time series images
reg_var_list <- list.files(path=file.path(in_dir,"moore_NDVI_wgs84"),
                           pattern="moore_reg.*.MOD13A2_A.*04072014.*.tif$",full.name=TRUE)
#moore_reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif

r_stack <- stack(reg_var_list)
levelplot(r_stack,layers=1:12) #show first 12 images (half a year more or less)

#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#out_dir <- getwd() #to allow raster package to work...

################# PART 1: TEMPORAL PREDICTIONS ###################

#Hurricane August 17, 2007
day_event<-strftime(as.Date("2007.08.17",format="%Y.%m.%d"),"%j")
#153,154
names(r_stack)
grep(paste("2007",day_event,sep=""),names(r_stack))

#r_huric_w <- mask(r_huric_w,moore_w)

r_stack_w <- mask(r_stack,rast_ref)

#### Now do Moore area!! about 31,091 pixels

pix_val <- as(r_stack_w,"SpatialPointsDataFrame")
#pix_val <- as.data.frame(r_stack_w)
#moore_dat <- as(moore_w,"SpatialPointsDataFrame")
r_huric_w <- subset(r_stack_w,152:156) #date before hurricane and  after

#pix_val <- extract(r_stack,moore_dat,df=TRUE)
#pix_val <- t(pix_val[,1:153])
pix_val2 <- as.data.frame(pix_val)
pix_val2 <-  pix_val2[,1:152] 
pix_val2 <- as.data.frame(t(as.matrix(pix_val2 )))#dim 152x26,616
#
#ttx2 <- lapply(pix_val,FUN=raster_ts_arima_predict,na.rm=T,arima_order=NULL,n_ahead=2)
# <- mclapply(pix_val,,)

n_pred_ahead <- 4 #number of temporal ARIMA predictions ahead..

## Now prepare predictions: should this a be a function?

list_param_predict_arima_2 <- list(pix_val=pix_val2,na.rm=T,arima_order=NULL,n_ahead=n_pred_ahead)
#debug(raster_ts_arima_predict)
#tmp_val <- raster_ts_arima_predict(1,list_param_predict_arima_2)
ttx2 <- mclapply(1:length(pix_val2), FUN=raster_ts_arima_predict,list_param=list_param_predict_arima_2,mc.preschedule=FALSE,mc.cores = 11) 
save(ttx2,file=paste("raster_ts_arima_predict_obj",out_suffix,".RData",sep=""))

ttx3 <- load_obj("raster_ts_arima_predict_obj04212014.RData")
#ttx2 <- mclapply(1:length(pix_val), FUN=raster_ts_arima_predict,list_param=list_param_predict_arima_2,mc.preschedule=FALSE,mc.cores = 12) #This is the end bracket from mclapply(...) statement
#ttx2 <- lapply(1:length(pix_val2), FUN=raster_ts_arima_predict,list_param=list_param_predict_arima_2) #This is the end bracket from mclapply(...) statement

r_ref <- rast_ref
file_format <- ".rst"
NA_flag_val <- -9999

out_rastnames <- paste(paste("NDVI_pred_mooore_auto",1:n_pred_ahead,sep="_"),"_",out_suffix,file_format,sep="")
list_param_arima_convert <- list(rast_ref,ttx3,file_format,out_dir,out_rastnames,file_format,NA_flag_val)
names(list_param_arima_convert) <- c("r_ref","ttx","file_format","out_dir","out_rastnames","file_format","NA_flag_val")

#debug(convert_arima_pred_to_raster)
## Convert predicted values to raster...

pred_t_l<-lapply(1:n_pred_ahead,FUN=convert_arima_pred_to_raster,list_param=list_param_arima_convert)

pred_t_l <-unlist(pred_t_l)
r_pred  <- stack(pred_t_l[-c(grep(pattern="error",pred_t_l))])
r_error <- stack(pred_t_l[c(grep(pattern="error",pred_t_l))])

r_huric_w <- subset(r_stack,153:156)
#r_huric_w <- crop(r_huric_w,moore_w)

r_t0_pred <- stack(subset(r_huric_w,1),r_pred_t1,r_pred_t2)
names(r_t0_pred) <- c("NDVI_t_0","NDVI_pred_t_1","NDVI_pred_t_2")

#### Analyses of Model fitted...

l_arima_mod <- list.files(path=out_dir,pattern="arima_mod.*.RData",full.names=T)

list_param_extract_arima <- list(arima_mod_name=l_arima_mod)

#debug(extract_arima_mod_info)
test <- extract_arima_mod_info(1,list_param_extract_arima)

#tmp_val <- raster_ts_arima_predict(1,list_param_predict_arima_2)
l_arima_info <- mclapply(1:length(l_arima_mod), FUN=extract_arima_mod_info,list_param=list_param_extract_arima,
                         mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement


save(l_arima_info,file=paste("l_arima_info_obj",out_suffix,".RData",sep=""))
#l_arima_info <-load_obj("l_arima_mod_info_obj04212014.RData")
#mod_specification <- mclapply(l_arima_info[1:11],FUN=function(i){i[[1]]},mc.preschedule=FALSE,mc.cores = 11)
#test_name <- list.files(".","l_arima_mod_info_obj.*.RData")
mod_specification <- mclapply(l_arima_info,FUN=function(i){i[[1]]},mc.preschedule=FALSE,mc.cores = 11)
mod_specification2 <- mod_specification[1:26216]
#[1] 2 0 0 0 1 0 0
#A compact form of the specification, as a vector giving the number of AR (1), MA (2), 
#seasonal AR (3) and seasonal MA coefficients (4), 
#plus the period (5) and the number of non-seasonal (6) and seasonal differences (7).

fitted_mod_spec <- as.data.frame(do.call(rbind,mod_specification2))
#should add  id or coordinates!!!!so you can map it!!
names(fitted_mod_spec) <- c("AR","MA","S_AR","S_MA","P","D","S_D")
arima_mod_str <- lapply(1:nrow(fitted_mod_spec),FUN=function(i){paste(c(fitted_mod_spec[i,]),collapse="_")})
fitted_mod_spec$mod <- as.character(arima_mod_str)
fitted_mod_spec$mod_fac <- as.factor(fitted_mod_spec$mod)

#read.table(file=paste("fitted_mod_spec_tb","_",out_suffix,".txt",)

write.table(fitted_mod_spec,file=paste("fitted_mod_spec_tb","_",out_suffix,".txt",col.names=T,row.names=F,sep=","))

## Now plot
            
p1 <- histogram(fitted_mod_spec$AR) #how to present everything at once?? should this be transposed?
p2 <- histogram(fitted_mod_spec$MA)
p3 <- histogram(fitted_mod_spec$P)
p4 <- histogram(fitted_mod_spec$D)
p3 <- histogram(fitted_mod_spec$mod_fac)

#p3 <- histogram(fitted_mod_spec$mod)

p3 <- histogram(fitted_mod_spec$mod_fac)
dat_arima <- fitted_mod_spec
coordinates(dat_arima) <- coordinates(pix_val)

#rasterize(dat_arima)
r_ar <- rasterize(dat_arima,rast_ref,field="AR") #this is the prediction from lm model
r_ma <- rasterize(dat_arima,rast_ref,field="MA") #this is the prediction from lm model
r_p <- rasterize(dat_arima,rast_ref,field="P") #this is the prediction from lm model
r_d <- rasterize(dat_arima,rast_ref,field="D") #this is the prediction from lm model

r_arima_s <- stack(r_ar,r_ma,r_p,r_d)
names(r_arima_s) <- c("AR","MA","P","D")
plot(r_arima_s)

plot(r_ar,col=c("red","blue","green","grey","black"),main="AR")
plot(r_ma,col=c("red","blue","green","grey","black"),main="MA")

################### END OF SCRIPT ################