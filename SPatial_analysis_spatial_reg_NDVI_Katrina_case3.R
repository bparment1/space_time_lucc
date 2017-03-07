####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict Light data values in the Hurricane Katrina New Orleans region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 03/07/2017
#Version: 3
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#PROJECT: AAG 2015 in Chicago, with Marco Millones
#PROJECT: Geocomputation conference in Dallas with Marco Millones
#
#COMMENTS: - Testing alternative methods to eigen for spatial predictions: "Chebyshev" on new light data
#         - clean up and organize code to be more general for any dataset
#TO DO:
# - add confidence interval around reg coef: this is important!!
# - add variance around the MAE values in the accuracy assessment
# - modify parallelization so that it works both on windows and linux/macos
# - automation to call from the terminal/shell
#
#
#COMMIT: more changes to ARIMA models
#
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
library(zoo)  #original times series function functionalies and objects/classes
library(xts)  #extension of time series objects and functionalities
library(forecast) #arima and other time series methods
library(lubridate) #date and time handling tools
library(colorRamps) #contains matlab.like color palette
library(rgeos) #spatial analysis, topological and geometric operations e.g. interesect, union, contain etc.
library(sphet) #spatial analyis, regression eg.contains spreg for gmm estimation

###### Functions used in this script

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_03072017_functions.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_01092016.R" #PARAM 1
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir <- "~/Data/Space_beats_time/case3data/" #lights/table" #PARAM3
#in_dir <- "/home/parmentier/Data/Space_beats_time/Case2_data_NDVI/"
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/Case2_data_NDVI/"

#in_dir <- "~/Data/Space_beats_time/case3data/lights/table"
#in_dir <- "~/Data/Space_beats_time/Case1a_data"
#out_dir <- "/home/parmentier/Data/Space_beats_time/outputs"
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs"

#moore_window <- file.path(in_dir,"window_test4.rst")
#winds_zones_fname <- file.path(in_dir,"00_windzones_moore_sin.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"NDVI_Katrina_03032017" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
#data_fname <- file.path(in_dir,"lights/table","Kat_lights.txt") #PARAM 10
#data_fname <- file.path(in_dir,"output_Katrina_04082015","dat_reg_var_list_NDVI_Katrina_04082015.txt")
data_fname <- file.path(in_dir,"dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt")


coord_names <- c("x","y") #PARAM 11
#coord_names <- c("Long","Lat") #PARAM 11
#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12

var_names <- 1:230 #PARAM 13 #Data is stored in the columns 3 to 22
#num_cores <- 11 #PARAM 14
num_cores <- 4 #PARAM 14

n_time_event <- 108 #PARAM 15 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
time_window_selected <- var_names #PARAM 16: use alll dates for now
time_window_selected <- 100:116 #PARAM 16: use alll dates for now

re_initialize_arima <- T #PARAM 17, use re-initialization ie apply arima model with one step forward at each time step
previous_step <- T #PARAM 18

#date_range1 <- c("2001.01.01","2012.12.31") #EDGY DEAN
#date_range2 <- c("1992.01.01","2013.12.31") #Light Katrina: annual
date_range3 <- c("2001.01.01","2010.12.31") #NDVI Katrina

#dates1 <- generate_dates_by_step(date_range1[1],date_range1[2],16)$dates
#dates2 <- unique(year(generate_dates_by_step(date_range2[1],date_range2[2],1)$dates)) #extract year
dates3 <- generate_dates_by_step(date_range3[1],date_range3[2],16)$dates #NDVI Katrina

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

##specific processing done for srm
#r_dem <- raster("/data/project/layers/commons/data_workflow/inputs/dem-cgiar-srtm-1km-tif/srtm_1km.tif")
#r_dem_Katrina<-crop(r_dem,r_stack)
#writeRaster(r_dem_Katrina,file.path(in_dir,"r_srtm_Katrina.rst"))
#library(ggmap)
### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

data_tb <-read.table(data_fname,sep=",",header=T)
#data_tb$ezone_c[data_tb$Elev_zone 

#add elevation?
#attach(data_tb)
#data_tb$ezone_c[Elev_Zone < 3] <- 1
#data_tb$ezone_c[Elev_Zone == 3] <- NA  
#data_tb$ezone_c[Elev_Zone == 4] <- 2  
#detach(data_tb)

#### Make this a function...that will run automatically the predictions
## from lines 152-265

#Transform table text file into a raster image

#This function is very slow and inefficienct, needs improvement (add parallelization)
l_rast <- rasterize_df_fun(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format,NA_flag_val,tolerance_val=0.000120005)

#debug(rasterize_df_fun)
s_raster <- stack(l_rast) #stack with all the variables
projection(s_raster) <- CRS_reg
names(s_raster) <- names(data_tb)                
r_FID <- subset(s_raster,1) #Assumes ID or reference image is the first image of the stack

##Figure 1: reference layer

plot(r_FID,main=l_rast[1])
head(freq(r_FID))

##Figure 2: zonal layer

plot(subset(s_raster,zonal_colnames),main=zonal_colnames)

reg_var_list <- l_rast[var_names] #only select population raster
r_stack <- stack(reg_var_list)
projection(r_stack) <- CRS_reg
names(r_stack) <- names(data_tb)[var_names]

#Later rerport this basic information in  a text file 
dim(r_stack) #34x49x230 
ncell(s_raster) #1666
freq(r_FID,value=NA) #122
ncell(s_raster) - freq(r_FID,value=NA) #1544
res(r_stack) #about 1km

#Automate this step?

## Figure 3: visualization of time series

#projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
levelplot(r_stack,layers=1:4,col.regions=matlab.like(125)) #show first four images (half a year more or less)
plot(r_stack,y=1:4)

## Figure 4: visualization of event in time series

n_time_event #108
dates3[108]
#[1] "2005-08-29"

## plottting after and before event
n_before <- n_time_event - 5
n_after <- n_time_event + 5

levelplot(r_stack,layers=n_before:n_after,col.regions=matlab.like(125))

## Figure 5: visualization of histogram of event in time series
histogram(subset(r_stack,n_before:n_after))

## Figure 6: Profile for the time series

mean_vals <- colMeans(data_tb[,var_names],na.rm=T)
pixval <- data_tb[800,var_names]
#pix300 <- data_tb[300,var_names]

plot(1:length(var_names),mean_vals,type="b",ylab="var",xlab="time step",
     main="Average variable and pixel 800 profile")
lines(1:length(var_names),pixval,type="b",ylab="var",xlab="time step",col=c("red"),
     main="Average variable and pixel 800 profile")
legend("bottomleft",legend=c("Overall average VAR ","PIX 800 "),
        col=c("black","red"),lty=c(1,2))
abline(v=n_time_event,col="blue")
## Figure 5: visualization of histogram of event in time series

#title("Average pop per elevation zones (observed data)")
## By zone/strata

r_zonal <- subset(s_raster,zonal_colnames)
zones_tb_avg<- zonal(r_stack,r_zonal,fun='mean')

zones_avg_df <- as.data.frame(zones_tb_avg)
n_zones <- length(unique(zones_avg_df$zone))

n_time <- ncol(zones_avg_df) -1
#pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")

mydata<- zones_avg_df
#dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
dd <- do.call(make.groups, mydata) 
#dd$lag <- mydata$lag 
dd$zones <- mydata$zones
#dd$method <- mydata$method
#drop first few rows that contain no data but zones...
n_start <-n_zones +1
dd <- dd[n_start:nrow(dd),]
tmp_time <-unlist(lapply(1:n_time,FUN=function(i){rep(i,n_zones)}))
dd$time  <- tmp_time
dd$zones <- mydata$zone #use recycle rule

xyplot(data~time |zones,#group=method,
       data=dd,type="b",xlab="time",ylab="VAR",
       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2)
)

xyplot(data~time,group=zones,
       data=dd,type="b",xlab="time",ylab="VAR",
       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2),
      main="Average by zones for VAR"
)

write.table(zones_avg_df,file=paste("zones_avg_df","_",out_suffix,".txt",sep=""))
write.table(dd,file=paste("zones_avg_df_long_table","_",out_suffix,".txt",sep=""))

##############################################################################################
############## PART V PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER MULTIPLE time steps ####

##This times will we use an automated function to generate predictions over multiple dates

##########################
#### RUN FOR SELECTED DATES and the three methods..

###########SPATIAL METHODS
## Now predict for four dates using "mle": this does not work for such large area such as EDGY!!

### Predict using spatial regression: this should be a master function...
#r_spat_var <- subset(r_stack,139:161) #predict before (step 152) and three dates after (step 153)
### Important:
#As input!!
#r_spat_var should contain one image before the on being predicted, i.e. to predict step 100, you need at least 
#a stack of raster 99 and raster 100.

#num_cores_tmp <- 11
time_step_start <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
time_step_end <- n_time_event + 8
time_step_subset <- time_step_start -1 #use 99
time_window_selected <- time_step_subset:time_step_end

time_window_predicted <- time_step_start:time_step_end #100 to 116
r_spat_var <- subset(s_raster,time_window_selected) #predict before and after event

rast_ref <- subset(r_stack,1)
#rast_ref <- rast_ref != NA_flag_val
#pix_id_r <- rast_ref
#values(pix_id_r) <- 1:ncell(rast_ref) #create an image with pixel id for every observation
#pix_id_r <- mask(pix_id_r,rast_ref) #3854 pixels from which 779 pixels are NA

########
### TEST SPATIAL Prediction on one date....
                
### CLEAN OUT AND SCREEN NA and list of neighbours
#Let's use our function we created to clean out neighbours
#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
#r_var <- subset(r_stack,n_time_event) #this is the date we want to use to run the spatial regression
r_clip <- rast_ref #this is the image defining the study area

### MLE EIGEN Using previous step first

list_models <- NULL
proj_str <- NULL #if null the raster images are not reprojected
#the ouput suffix was wrong, needs to be 153!!!
#Use 100 to 116
#out_suffix_s <- paste("t_",100:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!

#estimator <- "mle"
estimator <- "mle"
estimation_method <- "eigen"
#estimation_method <- "LU"
#estimation_method <- "Chebyshev"
#estimator <- gmm

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_03022017_functions.R" #PARAM 1
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

previous_step <- T #Use t-1 neighbour and t-1 value to predict t0

out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!

#list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,previous_step)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator","estimation_method","previous_step")
n_pred <- nlayers(r_spat_var) - 1 # minus one because the first date is not predicted
#debug(predict_spat_reg_fun)
#predict_spat_reg_fun(1,list_param=list_param_spat_reg)

pred_spat_mle_eigen_with_previous  <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,
                                               list_param=list_param_spat_reg,
                                               mc.preschedule=FALSE,
                                               mc.cores = num_cores)

save(pred_spat_mle_eigen_with_previous,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,
                                                      "_","with_previous_",out_suffix,".RData",sep="")))

#pred_spat_mle_eigen: extract raster images from object
spat_pred_rast_mle_eigen_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast_mle_eigen_with_previous <- stack(lapply(pred_spat_mle_eigen_with_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast_mle_eigen_with_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_eigen_with_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

### MLE EIGEN not using previous step first

previous_step <- F

#list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
list_param_spat_reg<- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,previous_step)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator","estimation_method","previous_step")
n_pred <- nlayers(r_spat_var) -1


#debug(predict_spat_reg_fun)
#predict_spat_reg_fun(1,list_param=list_param_spat_reg)

#debug(predict_spat_reg_fun)
#predict_spat_reg_fun(1,list_param=list_param_spat_reg)

pred_spat_mle_eigen_no_previous <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,
                                            list_param=list_param_spat_reg,
                                            mc.preschedule=FALSE,
                                            mc.cores = num_cores)

save(pred_spat_mle_eigen_no_previous,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,"_",
                                                                  "_","no_previous_",out_suffix,".RData",sep="")))

#pred_spat_mle_eigen: extract raster images from object
spat_pred_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast_mle_eigen_no_previous <- stack(lapply(pred_spat_mle_eigen_no_previous,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast_mle_eigen_no_previous,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_eigen_no_previous,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

###########TEMPORAL METHODS
## Predict using temporal info: time steps 

### Use ARIMA!!! with time before

#out_dir
#r_clip
#proj_str
#file_format

r_temp_var <- subset(s_raster,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
list_models <-NULL

#the ouput suffix was wrong, needs to be 153!!!
#Use 100 to 116
out_suffix_s <- paste("t_",100:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
estimator <- "lm"
estimation_method <-"ols"

#ARIMA specific

num_cores_tmp <- num_cores
time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
time_step_subset <- time_step -1 #use 99

rast_ref <- subset(s_raster,1) #first image ID
r_stack_arima <- mask(r_stack,rast_ref)

#r_stack <- r_stack_arima
arima_order <- NULL
n_pred <- nlayers(r_temp_var) -1
n_pred_ahead <- n_pred

list_param_temp_reg_lm <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
names(list_param_temp_reg_lm) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
#debug(predict_temp_reg_fun)
#test_temp <- predict_temp_reg_fun(14,list_param_temp_reg_lm)
#plot(raster(test_temp$raster_pred),main=basename(test_temp$raster_pred))


#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

pred_temp_lm <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg_lm) 
#note that ojbect saved for each time step: i refers to the time step but prediciton is one step ahead!!
#-rw-r--r-- 1 parmentier parmentier      601 May 15 10:26 r_temp_res_lm_ols_17_t_100_NDVI_Katrina_05152015.rdc
#-rw-r--r-- 1 parmentier parmentier  9734856 May 15 10:26 temp_reg_obj_lm_ols_t_16_t_100_NDVI_Katrina_05152015.RData

save(pred_temp_lm,file= file.path(out_dir, paste("pred_temp_lm_",out_suffix,".RData",sep="")))

#########
##Run ARIMA models

#There are 954 pixels predicted
#not run here...

r_temp_var <- subset(s_raster,time_window_selected) # relies on the previous date in contrast to spat reg, this is r_var param
list_models <-NULL
#the ouput suffix was wrong, needs to be 153!!!
#Use 100 to 116
#out_suffix_s <- paste("t_",100:length(time_window_selected),"_",out_suffix,sep="")#this should really be automated!!!
out_suffix_s <- paste("t_",time_window_predicted,"_",out_suffix,sep="")#this should really be automated!!!
num_cores_tmp <- num_cores

#ARIMA specific parameters

time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
n_pred <- nlayers(r_temp_var) -1
n_pred_ahead <- n_pred
rast_ref <- subset(s_raster,1) #first image ID
r_stack_arima <- mask(r_stack,rast_ref)
#r_stack <- r_stack_arima
arima_order <- NULL

estimator <- "arima"
estimation_method <-"arima"
r_clip_tmp <- NULL
r_clip_tmp <- rast_ref

list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")

#debug(predict_temp_reg_fun)

if(re_initialize_arima==T){
   l_pred_temp_arima <- vector("list",length=n_pred)

  for(i in 1:n_pred){
    n_pred_tmp <- 1
    n_pred_ahead <- n_pred_tmp
    #time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99

    #n_time_pred_start <- 99 + i
    n_time_pred_start <- time_window_selected[1]  + i
    time_step <-n_time_pred_start
    #n_time_pred_end <- n_time_pred_start + length(time_window_selected)
    n_time_pred_end <- time_window_selected[length(time_window_selected)]
    
    out_suffix_s <- paste("t_",n_time_pred_start:n_time_pred_end,"_",out_suffix,sep="")#this should really be automated!!!
    list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
    names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
    #debug(predict_temp_reg_fun)
    pred_temp_arima <- predict_temp_reg_fun(i,list_param_temp_reg) #only one date predicted...four step ahead
    l_pred_temp_arima[[i]] <- pred_temp_arima  #only one date predicted...one step ahead
  }
  save(l_pred_temp_arima,file= file.path(out_dir, paste("l_pred_temp_arima_",out_suffix,".RData",sep="")))

}else{
  time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
  #n_pred <- 16
  n_pred <- length(time_window_selected) -1
  n_pred_ahead <- n_pred

  #n_time_pred_start <- 100 
  n_time_pred_start <- time_window_selected[1] 

  n_time_pred_end <- n_time_pred_start + length(time_window_selected)
    
  out_suffix_s <- paste("t_",n_time_pred_start:n_time_pred_end,"_",out_suffix,sep="")#this should really be automated!!!
  list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
   names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")

  #debug(predict_temp_reg_fun)
  pred_temp_arima <- predict_temp_reg_fun(1,list_param_temp_reg) #only one date predicted...four step ahead

}
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#using 11 cores
#pred_temp_arima <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)

#pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_153_EDGY_predictions_03182015.RData")
#pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_100_NDVI_Katrina_04102015.RData")

#extract_files <- function(i,x){obj<-x[[i]];obj$raster_pred}
#debug(extract_files)
#extract_files(1,l_pred_temp_arima)
r_temp_pred_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=l_pred_temp_arima)))
r_temp_res_rast_arima <- stack(unlist(lapply(1:length(l_pred_temp_arima),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=l_pred_temp_arima)))

r_temp_pred_rast_lm <- stack(unlist(lapply(1:length(pred_temp_lm),FUN=function(i,x){obj<-x[[i]];obj$raster_pred},x=pred_temp_lm)))
r_temp_res_rast_lm <- stack(unlist(lapply(1:length(pred_temp_lm),FUN=function(i,x){obj<-x[[i]];obj$raster_res},x=pred_temp_lm)))

levelplot(r_temp_pred_rast_lm,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(r_temp_res_rast_lm,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(r_temp_pred_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(r_temp_pred_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.

#levelplot(r_temp_res_rast_arima,col.regions=matlab.like(255),main="Var residuals after hurricane")

projection(r_temp_pred_rast_arima) <- CRS_reg

#r_temp_pred_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
#r_temp_res_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
#levelplot(r_temp_pred_rast_lm,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
#levelplot(r_temp_res_rast_lm,col.regions=rev(terrain.colors(255)),main="Var residuals after hurricane")

levelplot(spat_pred_rast_mle_eigen,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_eigen,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_eigen,col.regions=matlab.like(255)) #view the four predictions using mle spatial reg.

projection(r_temp_pred_rast_arima) <- CRS_reg
projection(r_temp_res_rast_arima) <- CRS_reg
projection(r_temp_res_rast_lm) <- CRS_reg
projection(r_temp_pred_rast_lm) <- CRS_reg

             
############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

projection(spat_pred_rast_mle_eigen) <- CRS_reg
projection(r_temp_pred_rast_arima) <- CRS_reg
projection(r_temp_pred_rast_lm) <- CRS_reg

#select specific model:
spat_pred_rast <- subset(spat_pred_rast_mle_eigen,2:length(time_window_selected)) #spatial model prediction: 101 to 116
#temp_pred_rast <- r_temp_pred_rast_lm #temporal model prediction
time_selected <- length(time_window_selected)-1
temp_pred_rast_arima <- subset(r_temp_pred_rast_arima,1:time_selected) #temporal model prediction
temp_pred_rast_lm <- r_temp_pred_rast_lm #temporal model prediction, lm, from time step 2 to 17

rast_zonal <- subset(s_raster,zonal_colnames)

projection(rast_ref) <- CRS_reg
projection(rast_zonal) <- CRS_reg
projection(s_raster) <- CRS_reg
projection(spat_pred_rast) <- CRS_reg #this should be CRS_interp
projection(temp_pred_rast_arima) <- CRS_reg
projection(temp_pred_rast_lm) <- CRS_reg

r_huric_obs <- subset(s_raster,time_window_selected[-1])

#r_huric_w <- crop(r_huric_w,rast_ref)
levelplot(r_huric_obs,col.regions=matlab.like(25))
levelplot(temp_pred_rast_arima,col.regions=matlab.like(25))

#r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
res_temp_s_lm <- temp_pred_rast_lm - r_huric_obs
res_temp_s_arima <- temp_pred_rast_arima - r_huric_obs
res_spat_s <- spat_pred_rast - r_huric_obs

names(res_temp_s_lm) <- sub("pred","res",names(res_temp_s_lm))
names(res_temp_s_arima) <- sub("pred","res",names(res_temp_s_arima))

#names(res_spat_mle_s) <- sub("pred","res",names(res_spat_mle_s))
names(res_spat_s) <- sub("pred","res",names(res_spat_s))
#names(res_temp_s) <- paste("r_res_s_",1:nlayers(res_temp_s),"_",out_suffix,sep="")

r_results <- stack(s_raster,rast_zonal,r_temp_pred_rast_arima,r_temp_pred_rast_lm,
                   spat_pred_rast_mle_eigen,spat_pred_rast_mle_Chebyshev,res_temp_s_arima,res_temp_s_lm,res_spat_s)

dat_out <- as.data.frame(r_results)
dat_out <- na.omit(dat_out)
write.table(dat_out,file=paste("dat_out_",out_suffix,".txt",sep=""),
            row.names=F,sep=",",col.names=T)

### Now accuracy assessment using MAE

out_suffix_s <- paste("temp_arima_",out_suffix,sep="_")

#undebug(calc_ac_stat_fun)
ac_temp_arima_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast_arima,r_var_s=r_huric_obs,r_zones=rast_zonal,
                                file_format=file_format,out_suffix=out_suffix_s)

out_suffix_s <- paste("temp_lm_",out_suffix,sep="_")

ac_temp_lm_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast_lm,r_var_s=r_huric_obs,r_zones=rast_zonal,
                                file_format=file_format,out_suffix=out_suffix_s)

out_suffix_s <- paste("spat_mle",out_suffix,sep="_")  
ac_spat_mle_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,r_var_s=r_huric_obs,r_zones=rast_zonal,
                                file_format=file_format,out_suffix=out_suffix_s)

#mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
#mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_mle_obj$mae_tb,ac_temp_arima_obj$mae_tb,ac_temp_lm_obj$mae_tb))

mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg","temp_arima","temp_lm")
#mae_tot_tb$time <- 2:nrow(mae_tot_tb)
#mae_tot_tb$time <- 2:n_pred
mae_tot_tb$time <- 2:n_pred

#mae_tot_tb$time <- 2:17

y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp_arima,mae_tot_tb$temp_lm))
plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range)
lines(temp_arima ~ time, type="b",col="magenta",data=mae_tot_tb)
lines(temp_lm ~ time, type="b",col="red",data=mae_tot_tb)

legend("topleft",legend=c("spat","temp_arima","temp_lm"),col=c("cyan","magenta","red"),lty=1)
title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!
write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))

#### BY ZONES ASSESSMENT: Ok now it is general so it should be part of the function...

#mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb[1:3,],
#                      ac_temp_obj$mae_zones_tb[1:3,])
mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb,
                      ac_temp_arima_obj$mae_zones_tb,ac_temp_lm_obj$mae_zones)

mae_zones_tb <- as.data.frame(mae_zones_tb)

n_zones <- length(unique(mae_zones_tb$zone))

mae_zones_tb$method <- c(rep("spat_reg",n_zones),rep("temp_arima_reg",n_zones),rep("temp_lm_reg",n_zones))

n_time <- ncol(mae_zones_tb) -1
pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")
#names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")
names(mae_zones_tb) <- pred_names

write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
#drop first few rows that contain no data but zones...
n_start <-n_zones*3 +1 #3 because we have 3 methods...
#n_start <-n_zones*2 +1
dd <- dd[n_start:nrow(dd),]

dd$zones <- mydata$zones
dd$method <- mydata$method

dd$zones <- mydata$zone #use recycle rule

#Note that to get the title correct, one needs to add
xyplot(data~which | as.factor(zones) ,group=method,data=dd,type="b",xlab="time",ylab="VAR",
      strip = strip.custom(factor.levels=unique(as.character(dd$zones))), #fix this!!!
      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2)
)

histogram(r_zonal)

################### END OF SCRIPT ##################