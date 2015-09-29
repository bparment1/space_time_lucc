###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(gtools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast)
library(xts)
library(zoo)
library(lubridate)

### Parameters and arguments

in_dir<- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan"
out_dir<- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan"
setwd(out_dir)

function_analyses_paper <-"MODIS_and_raster_processing_functions_09242013.R"
script_path<-in_dir #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

#infile_reg_outline=""  #input region outline defined by polygon: none for Venezuela
#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
#infile_reg_outline <- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/GYRS_MX_tri-state_latlong.shp"  #input region outline defined by polygon: Oregon
infile_reg_outline <- "/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/GYRS_MX_trisate_sin_windowed.shp"  #input region outline defined by polygon: Oregon

#ref_rast_name<-""  #local raster name defining resolution, exent, local projection--. set on the fly?? 
ref_rast_name<-"/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/gyrs_sin_mask_1km_windowed.rst"  #local raster name defining resolution, exent: oregon
#/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/GYRS_MX_trisate_sin_windowed.rst
#ref_rast_name<-"/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/GYRS_Municipios/GYRS_sin_mask_1km.rst"  #local raster name defining resolution, exent: oregon
ref_samp4_name <-"/Users/Parmentier/Google Drive/GA_group/Sample4.rst"
ref_EDGY_name <-"/Users/Parmentier/Google Drive/GA_group/EDGY_mask_sin_1km.rst"
infile_modis_grid<-"/Users/Parmentier/Documents/Benoit/Space_Time_Yucatan/modis_sinusoidal_grid_world.shp" #modis grid tiling system, global
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +units=m +no_defs";
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_interp <- proj_modis_str

out_suffix <-"09242013" #output suffix for the files that are masked for quality and for 
    
## Other specific parameters
NA_flag_val<- -9999
reg_var_list <- list.files(path=in_dir,pattern=paste("reg_mosaiced_MOD13A2_.*.",out_suffix,".rst$",sep=""))

r_stack <-stack(reg_var_list)

levelplot(r_stack,layers=1:12) #show first 12 images (half a year more or less)

reg_mask<-raster(ref_rast_name)
projection(reg_mask)<-CRS_interp


start_date<-"2001.01.01"
end_date <- "2012.12.31"
end_date <- "2002.12.31"

st <- as.Date(start_date,format="%Y.%m.%d")
en <- as.Date(end_date,format="%Y.%m.%d")
ll <- seq.Date(st, en, by="1 day")
ll <- seq.Date(st, en, by="16 day")
#dates_queried <- format(ll,"%Y.%m.%d")

generate_dates_modis <-function(start_date,end_date,step_date){
  library(xts)
  library(zoo)
  library(lubridate)
  
  st <- as.Date(start_date,format="%Y-%m-%d")
  en <- as.Date(end_date,format="%Y-%m-%d")
  year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
  }
  
  dates_modis <-as.Date(unlist((ll_list))) 
  return(dates_modis)
}
#debug(generate_dates_modis)
test_dates <-generate_dates_modis(start_date="2001-01-01",end_date="2002-12-31",step_date="16 day")
#Finish this part...write out
#names_hdf<-as.character(unlist(strsplit(x=basename(hdf), split="[.]")))
#product<-names_hdf[1]
#char_nb<-length(names_hdf)-2
#names_hdf <- names_hdf[1:char_nb]
#raster_name <- paste(paste(names_hdf,collapse="_"),"_",out_suffix,file_format,sep="")
leap_year(2012)
t44<-stack(reg_var_list[1:46])
ndvi_ts <- setZ(t44, idx)

ndvi_ts

x <- zApply(ndvi_ts, by=as.yearqtr, fun=mean, name="quarters") #aggregate times series by quarter
#names(SISmm) <- month.abb
x <- zApply(ndvi_ts, by=as.yearmon, fun=mean, name="month") #aggregate time series by month
#x <- zApply(ndvi_ts, by="month",fun=mean,name="month") #overall montlhy mean mean

x <- zApply(ndvi_ts, by="day",fun=mean,name="overall mean") #overall mean
x <- zApply(ndvi_ts, by=c(1,24),fun=mean,name="overall mean") #overall mean
r_date<-getZ(ndvi_ts)
#x <- apply.daily(ndvi_ts,FUN=mean) does not work
plot(x)

r_date <-generate_dates_modis(start_date="2001-01-01",end_date="2012-12-31",step_date="16 day")
r_date <-zoo(data.frame(r_index=1:length(r_date)),r_date)
month(r_date) #lubridate package 
mday(r_date) #find day of the mont!!
yday(r_date)
#format(r_date,"%m")
r_stack <-stack(reg_var_list)
setZ(r_stack,r_date)

input_selected <- r_date$r_index[day(r_date)==1] #January (first month)

input_selected <- r_date$r_index[month(r_date)==1] #January (first month)

clim_calc <- function(input_selected,r_stack,r_date=NULL){
  library(xts)
  library(zoo)
  library(lubridate)
  
  #r_date <- getZ(r_stack)
  #r_date <-zoo(data.frame(r_index=1:length(r_date)),r_date)
  
  #input_selected<-c(1,24)
  test <-calc(subset(r_stack,input_selected),mean,na.rm=TRUE)
  is_not_na <- function(x){sum(!is.na(x))}
  nobs <-calc(subset(r_stack,input_selected),is_not_na)
  r_obj<-stack(test,nobs)
  layerNames(r_obj)<-c("mean","nobs")
  return(r_obj)
}
#pix<-as.vector(pixel)
pixel<-as.vector(r_stack[300,500,]) #extract pixel at row 300 and column 500)
plot(pixel,type="b")
arima_order<-c(1,0,0) #order p,d,q for AR, difference, MA
arima_obj<-arima(pixel,order=arima_order)
auto.arima(pixel)

acf(as.ts(pixel))

#Take a subset...
raster_ts_arima<-function(pixel,na.rm=T,arima_order){
  arima_obj<-arima(pixel,order=arima_order)
  a<-as.numeric(coef(arima_obj)[1]) 
  return(a)
}
raster_ts_arima_fun<-function(pixel,list_param){
  arima_order <-list_param$arima_order
  arima_obj<-arima(pixel,order=arima_order)
  a<-as.numeric(coef(arima_obj)[1]) 
  return(a)
}


raster_test_fun<-function(pixel,na.rm){
  x<-sum(pixel,na.rm=T)
  #x<-sum(pixel)
  return(x)
}

tx<-raster_test_fun(pixel)
tx<-raster_ts_arima(pixel,na.rm=T,arima_order=arima_order)

### GET A SMALLER IMAGE FOR TESTING...
r1<-subset(r_stack,1)
plot(r1)
click(r1,id=TRUE,xy=TRUE,cell=TRUE)
cell_id1 <-127559
plot(r1)
e<-drawExtent()

r_crop<-crop(r1,e)

r_e<-crop(r_stack,e)
writeRaster(r_e,NAflag=NA_flag_val,filename="test_subset_10172013.tif",overwrite=TRUE)  
arima_order<-c(1,0,0) #order p,d,q for AR, difference, MA
test44 <-calc(r_e,fun=raster_ts_arima)
test44 <-calc(r_e,fun=raster_test_fun)
test45 <-calc(r_e,fun=sum,na.rm=T)
test45<-stackApply(r_e,indices=rep(1,length=nlayers(r_e)),fun=raster_test_fun,na.rm=T)
list_param<-list(arima_order=arima_order)
test45<-stackApply(r_e,indices=rep(1,length=nlayers(r_e)),raster_ts_arima_fun,list_param=list_param,na.rm=T)
test45<-stackApply(r_e,indices=rep(1,length=nlayers(r_e)),raster_ts_arima_fun,na.rm=T,list_param=list_param)
test45<-stackApply(r_e,indices=rep(1,length=nlayers(r_e)),raster_ts_arima,na.rm=T,arima_order=arima_order)

plot(test44)

test<-clim_calc(input_selected,r_stack)

x <- zApply(r_stack, by=input_selected,fun=mean,name="overall mean") #overall mean

test2<-subset(ndvi_ts,1) + subset(ndvi_ts,24)
#apply.daily(x, FUN, ...)
#apply.weekly(x, FUN, ...)
#apply.monthly(x, FUN, ...)
#apply.quarterly(x, FUN, ...)
#apply.yearly(x, FUN, ...)
#use auto.arima, arima,

url <- "~/Datos/Cressie/"
sst.dat = read.table(paste(url, "SST011970_032003.dat", sep=''), header = FALSE) 
sst.ll = read.table(paste(url, "SSTlonlat.dat", sep=''), header = FALSE)

spSST <- SpatialPointsDataFrame(sst.ll, sst.dat)
gridded(spSST) <- TRUE
proj4string(spSST) = "+proj=longlat +datum=WGS84"
SST <- brick(spSST)

idx33 <- seq(as.Date('1970-01-01'), as.Date('2003-03-01'), by='month')
idx33 <- as.yearmon(idx33)
SST <- setZ(SST, idx)
names(SST) <- as.character(idx)
hovmoller(SST, contour=FALSE, panel=panel.levelplot.raster,
          yscale.components=yscale.raster.subticks,
          interpolate=TRUE, par.settings=RdBuTheme)

### test function
f4_tst <- function(x, a, filename) {
  out <- raster(x)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE)
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
    v <- v + a
    out <- writeValues(out, v, bs$row[i])
  }
  out <- writeStop(out)
  return(out)
}
t44<- f4_tst(r1, 5, filename='test44.grd')

#blockSize(s)

t44<- f4_tst(r_e, na.rm=T,arima_order, filename='test44.grd',reg_mask)
debug(f4_tst)
t44<- f4_tst(r_stack, na.rm=T,arima_order, filename='test44.grd',reg_mask)


f4_tst <- function(x, na.rm=T,arima_order, filename,mask_rast) {
  raster_ts_arima<-function(pixel,na.rm=T,arima_order){
    arima_obj<-try(arima(pixel,order=arima_order))
    if(inherits(arima_obj,"try-error")){
      a<-NA
    }else{
      a<-as.numeric(coef(arima_obj)[1])
    }
    return(a)
  }
  raster_test_fun<-function(pixel){
    x<-sum(pixel,na.rm=T)
    #x<-sum(pixel)
    return(x)
  }
  
  col_no<-ncol(x)
  out <- raster(x,1)
  bs <- blockSize(out)
  out <- writeStart(out, filename, overwrite=TRUE)
  #for (i in 1:bs$n) {
  for (i in 1:bs$n) {
    v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
    v_mask <- getValues(mask_rast, row=bs$row[i], nrows=bs$nrows[i] )
    v_out<- vector("numeric",length=nrow(v))
    for (j in 1:nrow(v)){
      
      if(is.na(v_mask[j])){
        v_out[j]<-NA
      }else{
        v_out[j]<-raster_ts_arima(as.vector(v[j,]),na.rm=T,arima_order)
      }
      #v[j,]<-raster_test_fun(as.vector(v[j,]))
    }
    out <- writeValues(out, v_out, bs$row[i])
  }
  
  out <- writeStop(out)
  return(out)
}


#Select area that was heavily affected...(10pix) --> matrix 276*10
#selecct area that was not heavily affected...(10pix)
#or pixel that underwent different land cover change transitions.


