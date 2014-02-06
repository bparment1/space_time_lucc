####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: ARIMA          #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#                        
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/27/2013 
#DATE MODIFIED: 01/25/2014
#Version: 2
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
  out_rastname <-list_param$out_rastnames[i]
  file_format <- list_param$file_format
  NA_flag_val <- list_param$NA_flag_val
  
  #start script
  pred_t <- lapply(ttx,FUN=function(x){x$pred[i]})
  tt_dat <- do.call(rbind,pred_t)
  tt_dat <- as.data.frame(tt_dat)
  pred_t <-as(r_ref,"SpatialPointsDataFrame")
  pred_t <- as.data.frame(pred_t)
  pred_t <- cbind(pred_t,tt_dat)
  coordinates(pred_t) <- cbind(pred_t$x,pred_t$y)
  raster_pred <- rasterize(pred_t,r_ref,"V1",fun=mean)
  writeRaster(raster_pred,NAflag=NA_flag_val,filename=file.path(out_dir,out_rastname),
              overwrite=TRUE)  
  return(out_rastname)
}

#freq_r_stack(r_stack){
#  
#}

####### Parameters and arguments

## Location on Benoit's laptop (Dropbox)
in_dir<- "/Users/benoitparmentier/Dropbox/Data/Space_Time"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/Space_Time"
setwd(out_dir)

function_analyses_paper <- "MODIS_and_raster_processing_functions_01252014.R"
script_path <- "/Users/Parmentier/Dropbox/Data/NCEAS/git_space_time_lucc/scripts_queue" #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
infile_reg_outline <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/GYRS_MX_trisate_sin_windowed.shp"  #input region outline defined by polygon: Oregon

#ref_rast_name<-""  #local raster name defining resolution, exent, local projection--. set on the fly?? 
ref_rast_name<-"/Users/benoitparmentier/Dropbox/Data/Space_Time/gyrs_sin_mask_1km_windowed.rst"  #local raster name defining resolution, exent: oregon
ref_samp4_name <-"/Users/benoitparmentier/Dropbox/Data/Space_Time/reg_Sample4.rst"
ref_EDGY_name <-"/Users/benoitparmentier/Dropbox/Data/Space_Time/reg_EDGY_mask_sin_1km.rst"
ref_egg_rings_gyr_name <-"/Users/benoitparmentier/Dropbox/Data/Space_Time/reg_egg_rings_gyr.rst"
infile_modis_grid<-"/Users/benoitparmentier/Dropbox/Data/Space_Time/modis_sinusoidal_grid_world.shp" #modis grid tiling system, global
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_interp <- proj_modis_str

out_suffix <-"01252014" #output suffix for the files that are masked for quality and for 
    
## Other specific parameters
NA_flag_val<- -9999
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
r_stack <- stack(reg_var_list)
levelplot(r_stack,layers=1:12) #show first 12 images (half a year more or less)

ref_rast_r  <- raster(ref_rast_name)
ref_samp4_r <- raster(ref_samp4_name)
mask_EDGY_r <- raster(ref_EDGY_name) #create raster image
egg_rings_gyr_r <- raster(ref_egg_rings_gyr_name)

projection(mask_EDGY_r) <- proj_modis_str #assign projection coord defined earlier
projection(ref_samp4_r) <- proj_modis_str
projection(ref_rast_r) <- proj_modis_str
projection(egg_rings_gyr_r) <- proj_modis_str

r_x <- ref_rast_r #create raster image that will contain x coordinates
r_y <- ref_rast_r #create raster image that will contain y coordiates
values(r_x) <- coordinates(ref_rast_r)[,1] #set values in raster image to x coord
values(r_y) <- coordinates(ref_rast_r)[,2] #set values in raster image to y coord
pix_id_r <- ref_rast_r
values(pix_id_r) <- 1:ncell(ref_rast_r)
s_dat_var <-stack(pix_id_r,r_x,r_y,egg_rings_gyr_r,mask_EDGY_r,ref_samp4_r)
layerNames(s_dat_var) <- c("pix_id_r","r_x","r_y","egg_rings_gyr_r",
                           "mask_EDGY_r","ref_samp4_r")
projection(s_dat_var) <- proj_modis_str
plot(s_dat_var)


################# PART 2: TEMPORAL PREDICTIONS ###################

#Hurricane August 17, 2007
day_event<-strftime(as.Date("2007.08.17",format="%Y.%m.%d"),"%j")
#153,154
layerNames(r_stack)
grep(paste("2007",day_event,sep=""),layerNames(r_stack))

### NOW ANALYSES WITH TIME AND SPACE...

## get small subset 1

# now extract these pixels and fit and arim for before period of hurricane
r_w <-raster("cropped_area.rst")
r_w_spdf <-as(r_w,"SpatialPointsDataFrame")
pix_val <- extract(r_stack,r_w_spdf,df=TRUE)

plot(pix_val[1,],type="b")
abline(v=153,col="red") # Find out fire event in the area...there are some very low
#values that are not related to the hurricane...
plot(pix_val[1,139:161],type="b")
abline(v=(154-139),col="red")

levelplot(r_stack,layer=152:155)
plot(subset(r_stack,154),legend=F)
plot(egg_rings_gyr_r,add=T,legend=F)
plot(ref_samp4_r,add=T,col=c("blue","red","black","yellow"))
plot(r_w,add=T,col="cyan",legend=F)
r_w_spdf<-as(r_w,"SpatialPointsDataFrame")

### Testing ARIMA model fitting with and without seasonality
#arima(USAccDeaths, order = c(0,1,1), seasonal = list(order=c(0,1,1)),
#plot(USAccDeaths)     
ts_x <- ts(data=pix_val[1,1:153],frequency=23) #frequency determines period...
frequency(ts_x) #frequency parameters is used in arima model fitting and spectrum
acf(ts_x,na.action=na.pass)
pacf(ts_x,na.action=na.pass)
acf(pix_val[1,1:153],na.action=na.pass) #wihtout a ts object
pacf(pix_val[1,1:153],na.action=na.pass)

arima_mod <- auto.arima(ts_x)
#arima_mod <- auto.arima(pix_val[1,1:153])
arima_mod1 <- auto.arima(pix_val[1,1:153])

p_arima<-predict(arima_mod,n.ahead=2)
p_arima1<-predict(arima_mod1,n.ahead=2)

#arima_mod <- arima(pix_val[1,1:153], order = c(1,1,1), 
#                   seasonal = list(order=c(0,1,1)),
            
plot(pix_val[1,152:155],type="b")
lines(c(pix_val[1,152:153],p_arima$pred),type="b",col="red") #prediction with seasonality
lines(c(pix_val[1,152:153],p_arima1$pred),type="b",col="blue")

raster_ts_arima(pix_val[1,1:153],na.rm=T,c(1,0,0))                        
raster_ts_arima(pix_val[1,1:153],na.rm=T,arima_order=c(0,0,2),n_ahead=2)                        

#Transform data to only include 1 to 153!!
pix <- as.data.frame(t(pix_val[,1:153]))

#acf(pix[1:153],na.action=na.pass)

tt<-raster_ts_arima_predict(pix[,1],na.rm=T,arima_order=NULL,n_ahead=2)
ttx<-lapply(pix,FUN=raster_ts_arima_predict,na.rm=T,arima_order=NULL,n_ahead=2)

tt_dat<-do.call(rbind,ttx)
class(tt_dat)

## Convert predicted values to raster...
r_pred_t1 <- r_w
r_pred_t2 <- r_w
values(r_pred_t1) <- tt_dat[,1]
values(r_pred_t2) <- tt_dat[,2]

r_huric_w <- subset(r_stack,153:155)
r_huric_w <- crop(r_huric_w,r_w)

r_t0_pred <- stack(subset(r_huric_w,1),r_pred_t1,r_pred_t2)
layerNames(r_t0_pred) <- c("NDVI_t_0","NDVI_pred_t_1","NDVI_pred_t_2")
dif_pred  <- r_t0_pred - r_huric_w 

#mfrow(c(2,3))
temp.colors <- colorRampPalette(c('blue', 'white', 'red'))
plot(stack(r_huric_w,r_t0_pred)) #compare visually predicted and actual NDVI
plot(dif_pred,col=temp.colors(25),colNA="black") #compare visually predicted and actual NDVI

sd_vals <- cellStats(dif_pred,"sd")
mean_vals <- cellStats(dif_pred,"mean")
#mode_vals <- cellStats(dif_pred,"mode")

hist(dif_pred)
dif_pred1 <-subset(dif_pred,2)
high_res<- dif_pred1>2*sd_vals[2] # select high residuals...
freq(high_res)

#There is large number of NA values after the hurricane
r_NA<-stack(raster_NA_image(r_huric_w))
freq(subset(r_NA,1)) # 10 NA cells
freq(subset(r_NA,2)) # 85 NA cells, after hurricane event...
freq(subset(r_NA,3)) # 26 NA cells

### apply to Ejido 2....

ref_samp4_spdf<- as(ref_samp4_r,"SpatialPointsDataFrame")
pix_val <- extract(r_stack,ref_samp4_spdf,df=TRUE)
#pix_val <- t(pix_val[,1:153])
pix_samp4 <- as.data.frame(t(pix_val[,1:153]))


#acf(pix[1:153],na.action=na.pass)
tt <- raster_ts_arima_predict(pix_samp4[,4],na.rm=T,arima_order=NULL,n_ahead=2)
ttx<-lapply(pix_samp4,FUN=raster_ts_arima_predict,na.rm=T,arima_order=NULL,n_ahead=2)

## Now of ejido

#ejido2 <- as.data.frame(ref_samp4_spdf)
#ejido2 <- ejido2[ejido2$reg_Sample4==2,]
ejido <- ref_samp4_spdf[ref_samp4_spdf$reg_Sample4==4,]

pix_val <- extract(r_stack,ejido,df=TRUE)
#pix_val <- t(pix_val[,1:153])
pix_ejido <- as.data.frame(t(pix_val[,1:153]))

ttx<-lapply(pix_ejido,FUN=raster_ts_arima_predict,na.rm=T,arima_order=NULL,n_ahead=2)
#ttx<-lapply(pix,   FUN=raster_ts_arima_predict,na.rm=T,arima_order=NULL,n_ahead=2)

pred_error <-lapply(ttx,FUN=function(x){x$error[1]})
tt_dat_error<-do.call(rbind,pred_error)
sum(tt_dat_error)

pred_t1 <-lapply(ttx,FUN=function(x){x$pred[1]})
tt_dat1<-do.call(rbind,pred_t1)

pred_t2 <-lapply(ttx,FUN=function(x){x$pred[2]})
tt_dat2<-do.call(rbind,pred_t2)

i<-1
sref_samp4_spdf
r_ref <- ref_samp4_r
raster_name <- "test.raster"
file_format <- ".rst"
NA_flag_val <- -9999

out_rastnames <- c("test1.rst","test2.rst")
list_param_arima_convert <- list(r_ref,ttx,file_format,out_dir,out_rastnames,file_format,NA_flag_val)
names(list_param_arima_convert) <- c("r_ref","ttx","file_format","out_dir","out_rastnames","file_format","NA_flag_val")

debug(convert_arima_pred_to_raster)
convert_arima_pred_to_raster(1,list_param=list_param_arima_convert)



## Convert predicted values to raster...
r_pred_t1 <- ref_samp4_r #timestep 1
r_pred_t2 <- ref_samp4_r #timestep 2
values(r_pred_t1) <- tt_dat1
values(r_pred_t2) <- tt_dat2

r_huric_w <- subset(r_stack,153:155)
r_huric_w <- crop(r_huric_w,r_w)

r_t0_pred <- stack(subset(r_huric_w,1),r_pred_t1,r_pred_t2)
layerNames(r_t0_pred) <- c("NDVI_t_0","NDVI_pred_t_1","NDVI_pred_t_2")
dif_pred  <- r_t0_pred - r_huric_w 
