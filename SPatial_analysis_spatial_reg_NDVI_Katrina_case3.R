####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict Light data values in the Hurricane Katrina New Orleans region. 
#Temporal predictions use OLS with the image of the previous time or the ARIMA method.
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 04/14/2015
#Version: 3
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#PROJECT: AAG 2015 with Marco Millones
#PROJECT: Geocomputation with Marco Millones
#
#COMMENTS: - Testing alternative methods to eigen for spatial predictions: "Chebyshev" on new light data
#         - clean up and organize code to be more general for any dataset
#TO DO:
# - add ARIMA function with parallelization: in process, must be modified for increased efficiency
# - add confidence interval around reg coef: this is important!!
# - add variance around the MAE values in the accuracy assessment
# - modify parallelization so that it works both on windows and linux/macos
# - automation to call from the terminal/shell
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

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_04102015_functions.R" #PARAM 1
script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir <- "~/Data/Space_beats_time/case3data/" #lights/table" #PARAM3
in_dir <- "/home/parmentier/Data/Space_beats_time/Case2_data_NDVI/"
#in_dir <- "~/Data/Space_beats_time/case3data/lights/table"
#in_dir <- "~/Data/Space_beats_time/Case1a_data"
#in_dir_NDVI <- file.path(in_dir,"moore_NDVI_wgs84") #contains NDVI 

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
out_suffix <-"NDVI_Katrina_04102015" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
#data_fname <- file.path(in_dir,"lights/table","Kat_lights.txt") #PARAM 10
#data_fname <- file.path(in_dir,"output_Katrina_04082015","dat_reg_var_list_NDVI_Katrina_04082015.txt")
data_fname <- file.path(in_dir,"dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt")

#coord_names <- c("Long","Lat") #PARAM 11
coord_names <- c("x","y") #PARAM 11

#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12

var_names <- 1:230 #PARAM 13 #Data is stored in the columns 3 to 22
num_cores <- 11 #PARAM 14
n_time_event <- 108 #PARAM 15 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
time_window_selected <- var_names #PARAM 16: use alll dates for now
time_window_selected <- 100:116 #PARAM 16: use alll dates for now
  
################# START SCRIPT ###############################

##specific processing done for srm
#r_dem <- raster("/data/project/layers/commons/data_workflow/inputs/dem-cgiar-srtm-1km-tif/srtm_1km.tif")
#r_dem_Katrina<-crop(r_dem,r_stack)
#writeRaster(r_dem_Katrina,file.path(in_dir,"r_srtm_Katrina.rst"))
library(ggmap)
### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#set up the working directory
#Create output directory

out_dir <- dirname(in_dir) #output will be created in the input dir
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

#Transform table text file into a raster image
#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

#This function is very slow and inefficienct, needs improvement (add parallelization)
l_rast <- rasterize_df_fun(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format,NA_flag_val,tolerance_val=0.000120005)

#debug(rasterize_df_fun)
s_raster <- stack(l_rast) #stack with all the variables
projection(s_raster) <- CRS_reg
names(s_raster) <- names(data_tb)                

r_FID <- subset(s_raster,1) #Assumes ID or reference image is the first image of the stack

##Figure 1: reference layer

plot(r_FID,main=l_rast[1])

freq(r_FID)

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
#freq(subset(s_raster,"NDVI_1"),value=NA) #
res(r_stack)

#Automate this step?

## Figure 3: visualization of time series

#projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
levelplot(r_stack,layers=1:4,col.regions=matlab.like(125)) #show first four images (half a year more or less)
plot(r_stack,y=1:4)

## Figure 4: visualization of event in time series

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
zones_tb_avg<- zonal(r_stack,r_zonal,stat='mean')

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

################ PART II : RUN SPATIAL REGRESSION ############

#Now mask and crop layer
r_var <- subset(r_stack,c(n_time_event,n_time_event+1)) #date before and after event

#r_clip <- rast_ref
#r_var <- crop(r_var,rast_ref) #crop/clip image using another reference image
rast_ref <- subset(r_stack,1)
rast_ref <- rast_ref != NA_flag_val
pix_id_r <- rast_ref
values(pix_id_r) <- 1:ncell(rast_ref) #create an image with pixel id for every observation
pix_id_r <- mask(pix_id_r,rast_ref) #3854 pixels from which 779 pixels are NA

########
### TEST SPATIAL Prediction on one date....
                
### CLEAN OUT AND SCREEN NA and list of neighbours
#Let's use our function we created to clean out neighbours
#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
r_var <- subset(r_stack,n_time_event) #this is the date we want to use to run the spatial regression
r_clip <- rast_ref #this is the image defining the study area

proj_str<- NULL #SRS/CRS projection system
out_suffix_s <- paste("d",n_time_event,out_suffix,sep="_")
#out_dir and out_suffix set earlier

nb_obj_for_pred_t_event <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)
names(nb_obj_for_pred_t_event) #show the structure of object (which is made up of a list)
r_poly_name <- nb_obj_for_pred_t_event$r_poly_name #name of the shapefile that has been cleaned out
reg_listw_w <- nb_obj_for_pred_t_event$r_listw #list of weights for cleaned out shapefile
#Use OGR to load the screened out data: note that we have now 2858 features with 3 fields
data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                    layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))

data_reg <- as.data.frame(data_reg_spdf) #convert spdf to df
#test a fake covariate
#Now run the spatial regression, since there are no covariates, the error and lag models are equivalent
#This takes a bit less than 3minutes for the dataset containing 2858 polygons
sam_esar_eigen <- errorsarlm(v1~ 1, listw=reg_listw_w,method="eigen",
                       data=data_reg,na.action=na.omit,zero.policy=TRUE,
                       tol.solve=1e-36) #tol.solve use in matrix operations
summary(sam_esar_eigen)
v5 <- rnorm(nrow(data_reg))
data_reg$v5 <- v5 - mean(v5) 
sam.esar2 <- errorsarlm(v1~ v5, listw=reg_listw_w, 
                       data=data_reg,na.action=na.omit,zero.policy=TRUE,
                       tol.solve=1e-36) #tol.solve use in matrix operations
summary(sam.esar2) #using randomfake variable as  covariate

print(coef(sam_esar_eigen))
print(coef(sam.esar2)) #almost equal could set the seed to make it reproducib

#Predicted values and
data_reg_spdf$spat_reg_pred <- sam_esar_eigen$fitted.values
data_reg_spdf$spat_reg_res <- sam_esar_eigen$residuals

spat_mod <- try(spautolm(v1 ~ 1,data=data_reg, listw= reg_listw_w,na.action=na.omit,zero.policy=TRUE,
                               tol.solve=1e-36))
spat_mod_spreg <- try(spreg(v1 ~ v2,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))

#res<-gstslshet(v1 ~ v2 , data=data_reg, listw=reg_listw_w)
spat_mod_spreg2 <- try(spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))

lm_mod <- try(lm(v1 ~ v5,data=data_reg))
#lm_mod2 <- try(lm(v1 ~ v4,data=data_reg))
mean(data_reg$v2)
mean(data_reg$v5)

all.equal(mean(data_reg$v2),as.numeric(coef(lm_mod)[1])) #ok this work...
#A work around may be to include a standar normal random variable!!


################### PART III RUN TEMPORAL MODEL USING LM ########

#### NOW PREDICTION USING PREVIOUS TIME STEP 

n_before <- n_time_event-1
r_var2 <- subset(r_stack,n_before:n_time_event)
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
############## PART IV:Produce images from the individual predictions using time and space ####

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
plot(r_pred)

levelplot(r_pred,regions.col=rev(terrain.colors(255)),main="Var predictions after disturbance event")
levelplot(r_pred,col.regions=matlab.like(25),main="Var predictions after disturbance event")

#### Examining difference between predictions
r_dif <- r_spat_pred - r_temp_pred
plot(r_dif,main="Difference between spatial and temporal models")
hist(r_dif,main="Difference between spatial and temporal models")

##############################################################################################
############## PART V PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER MULTIPLE time steps ####

##This times will we use an automated function to generate predictions over multiple dates

##########################
#### RUN FOR SELECTED DATES and the three methods..

###########SPATIAL METHODS
## Now predict for four dates using "mle": this does not work for such large area such as EDGY!!

### Predict using spatial regression: this should be a master function...
#r_spat_var <- subset(r_stack,139:161) #predict before (step 152) and three dates after (step 153)

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

### MLE Chebyshev

list_models <- NULL
proj_str <- NULL #if null the raster images are not reprojected
out_suffix_s <- paste("t_",1:length(time_window_selected),"_",out_suffix,sep="") #note that we set the time step in the output suffix!!!
#estimator <- "mle"
estimator <- "mle"
estimation_method <- "Chebyshev"
#estimation_method <- "LU"

#list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator","estimation_method")
n_pred <- nlayers(r_spat_var)

pred_spat_mle_chebyshev <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,mc.preschedule=FALSE,mc.cores = num_cores)
save(pred_spat_mle_chebyshev,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,"_",out_suffix,".RData",sep="")))

#pred_spat_mle_chebyshev: extract raster images from object
spat_pred_rast_mle_Chebyshev <- stack(lapply(pred_spat_mle_chebyshev,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast_mle_Chebyshev <- stack(lapply(pred_spat_mle_chebyshev,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast_mle_Chebyshev,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_Chebyshev,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

### MLE eigen
estimator <- "mle"
estimation_method <- "eigen"

#list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator","estimation_method")
n_pred <- nlayers(r_spat_var)

pred_spat_mle_eigen <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,mc.preschedule=FALSE,mc.cores = num_cores)
save(pred_spat_mle_eigen,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,"_",out_suffix,".RData",sep="")))

#pred_spat_mle_chebyshev: extract raster images from object
spat_pred_rast_mle_eigen <- stack(lapply(pred_spat_mle_eigen,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast_mle_eigen <- stack(lapply(pred_spat_mle_eigen,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast_mle_eigen,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_eigen,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

### OLS for didactive purpose

list_param_spat_reg$estimator <- "ols"
list_param_spat_reg$estimation_method <- "ols"

pred_spat_ols <- mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,mc.preschedule=FALSE,mc.cores = num_cores)
save(pred_spat_ols,file=file.path(out_dir,paste("pred_spat_",estimator,"_",estimation_method,"_",out_suffix,".RData",sep="")))

#pred_spat_mle_chebyshev: extract raster images from object
spat_pred_rast_ols <- stack(lapply(pred_spat_ols,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast_ols <- stack(lapply(pred_spat_ols,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast_ols,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_ols,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

## Alternative gmm etc.


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

num_cores_tmp <- 11
time_step <- n_time_event - 8 #this is the time step for which to start the arima model with, start at 99
n_pred_ahead <- 16
rast_ref <- subset(s_raster,1) #first image ID
r_stack_arima <- mask(r_stack,rast_ref)

#r_stack <- r_stack_arima
arima_order <- NULL

#list_param_temp_reg <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
#                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
#names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
#                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
#n_pred <- nlayers(r_temp_var) -1
#undebug(predict_temp_reg_fun)
#test_temp <- predict_temp_reg_fun(14,list_param_temp_reg)
#plot(raster(test_temp$raster_pred),main=basename(test_temp$raster_pred))

#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#pred_temp_lm <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg) 


##ARIMA

#not run here...

estimator <- "arima"
estimation_method <-"arima"
r_clip_tmp <- NULL
r_clip_tmp <- rast_ref

num_cores_tmp <- 11

list_param_temp_reg <- list(out_dir,r_temp_var,r_clip_tmp,proj_str,list_models,out_suffix_s,file_format,estimator,estimation_method,
                            num_cores_tmp,time_step,n_pred_ahead,r_stack_arima,arima_order,NA_flag_val)
names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format","estimator","estimation_method",
                            "num_cores","time_step","n_pred_ahead","r_stack","arima_order","NA_flag_val")
#n_pred <- nlayers(r_temp_var) -1
n_pred <- 16
#debug(predict_temp_reg_fun)
pred_temp_arima <- predict_temp_reg_fun(1,list_param_temp_reg) #only one date predicted...four step ahead
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#using 11 cores
#pred_temp_arima <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)

#pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_153_EDGY_predictions_03182015.RData")
pred_temp_arima <- load_obj("temp_reg_obj_arima_arima_t_100_NDVI_Katrina_04102015.RData")

r_temp_pred_rast_arima <- stack(pred_temp_arima$raster_pred)
r_temp_res_rast_arima <- stack(pred_temp_arima$raster_res)

#r_temp_pred_rast_arima <- stack(lapply(pred_temp_arima,FUN=function(x){x$raster_pred}))
#r_temp_res_rast_arima <- stack(lapply(pred_temp_arima,FUN=function(x){x$raster_res}))
levelplot(r_temp_pred_rast_arima,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
#levelplot(r_temp_res_rast_arima,col.regions=matlab.like(255),main="Var residuals after hurricane")

projection(r_temp_pred_rast_arima) <- CRS_WGS84

#r_temp_pred_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
#r_temp_res_rast_lm <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
#levelplot(r_temp_pred_rast_lm,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
#levelplot(r_temp_res_rast_lm,col.regions=rev(terrain.colors(255)),main="Var residuals after hurricane")

levelplot(spat_pred_rast_mle_eigen,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_mle_eigen,col.regions=rev(terrain.colors(255))) #view the four predictions using mle spatial reg.

projection(r_temp_pred_rast_lm) <- CRS_reg
projection(r_temp_res_rast_lm) <- CRS_reg

####### Comparing coefficients

## Extract spatial coefficients

#pred_spat_mle
#pred_spat_mle_Chebyshev_test <-load_obj("pred_spat_mle_chebyshev_EDGY_predictions_03092015.RData")
l_coef_mle_chebyshev <- lapply(pred_spat_mle_chebyshev,FUN=function(x){coef(x$spat_mod)})
tb_coef_mle_chebyshev <- as.data.frame(do.call(rbind,l_coef_mle_chebyshev))
tb_coef_mle_chebyshev$estimation_method <- "Chebyshev"

l_coef_mle_eigen <- lapply(pred_spat_mle_eigen,FUN=function(x){coef(x$spat_mod)})
tb_coef_mle_eigen <- as.data.frame(do.call(rbind,l_coef_mle_eigen))
tb_coef_mle_eigen$estimation_method <- "eigen"

#pred_spat_mle_Matrix_test <-load_obj("pred_spat_mle_Matrix_EDGY_predictions_03092015.RData")
#l_coef_mle_Matrix <- lapply(pred_spat_mle_Matrix_test,FUN=function(x){coef(x$spat_mod)})
#tb_coef_mle_Matrix <- as.data.frame(do.call(rbind,l_coef_mle_Matrix))
#tb_coef_mle_Matrix$estimation_method <- "Matrix"

#tb_coef_mle <- as.data.frame(do.call(rbind,list(tb_coef_mle_Chebyshev,tb_coef_mle_LU,tb_coef_mle_Matrix,tb_coef_mle_MC)))
tb_coef_mle <- as.data.frame(do.call(rbind,list(tb_coef_mle_chebyshev,tb_coef_mle_eigen))) #,tb_coef_mle_Matrix,tb_coef_mle_MC)))

tb_coef_mle$v2 <- NA                
tb_coef_mle <- tb_coef_mle[,c(2,1,4,3)]
names(tb_coef_mle)<- c("(Intercept)","rho","v2","estimation_method")
tb_coef_mle$time <- 1:n_pred             
tb_coef_mle$estimator <- "mle"
tb_coef_mle$method <- paste(tb_coef_mle$estimator,tb_coef_mle$estimation_method,sep="_")                

#View(tb_coef_mle)
head(tb_coef_mle)

#pred_spat_gmm
#pred_spat_gmm_test <-load_obj("pred_spat_gmm_EDGY_predictions_03092015.RData")
#l_coef_gmm <- lapply(pred_spat_gmm,FUN=function(x){coef(x$spat_mod)})
#tb_coef_gmm <- as.data.frame(t(do.call(cbind,(l_coef_gmm))))
#tb_coef_gmm <- tb_coef_gmm[,c(1,3,2)]
#names(tb_coef_gmm)<- c("(Intercept)","rho","v2")
#tb_coef_gmm$estimation_method <- "gmm"
#tb_coef_gmm$time <- 1:4                
#tb_coef_gmm$estimator <- "gmm"
#tb_coef_gmm$method <- "gmm"                


#tb_coef_mle  <- tb_coef_mle[,c(1,3,2)]
#names(tb_coef_mle)<- c("(Intercept)","rho")
#tb_coef_mle$v2<- NA
#tb_coef_mle$time <- 2001:2012
#tb_coef_mle$method <- "mle"
#pred_spat_ols <-load_obj("pred_spat_ols_EDGY_predictions_03092015.RData")
#pred_spat_ols
l_coef_ols <- lapply(pred_spat_ols,FUN=function(x){coef(x$spat_mod)})

tb_coef_ols <- as.data.frame(do.call(rbind,l_coef_ols))
tb_coef_ols <- tb_coef_ols[,c(1,2)]
names(tb_coef_ols)<- c("(Intercept)","rho")
tb_coef_ols$v2 <- NA   
tb_coef_ols$estimation_method <- "ols"
tb_coef_ols$time <- 1:n_pred #20 for this dataset                
tb_coef_ols$estimator <- "ols"
tb_coef_ols$method <- "ols"                


#tb_coef_sas_file <- file.path(in_dir,"EDGY_coeficients_SAS.csv")
#tb_coef_sas <- read.table(tb_coef_sas_file,sep=";",header=T)
#names(tb_coef_sas) <- names(tb_coef_gmm)

#names(tb_coef_sas)

names(tb_coef_mle)
names(tb_coef_ols)
#names(tb_coef_gmm)

#head(tb_coef_sas)


#tb_coef_method <- rbind(tb_coef_mle,tb_coef_gmm,tb_coef_ols)
#tb_coef_method <- rbind(tb_coef_gmm,tb_coef_sas,tb_coef_ols,tb_coef_mle)
tb_coef_method <- rbind(tb_coef_ols,tb_coef_mle)

xyplot(rho~time,groups=method,data=tb_coef_method,type="b",
                 auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2),
                main="Comparison of rho coefficients with different methods "
)

write.table(tb_coef_method,paste("tb_coef_method",out_suffix,".txt",sep=""),row.names=F,col.names=T)                

############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

#select specific model:
spat_pred_rast <- subset(spat_pred_rast_mle_eigen,2:length(time_window_selected)) #spatial model prediction
temp_pred_rast <- r_temp_pred_rast_lm #temporal model prediction
temp_pred_rast <- r_temp_pred_rast_arima #temporal model prediction

rast_zonal <- subset(s_raster,zonal_colnames)

projection(rast_ref) <- CRS_reg
projection(rast_zonal) <- CRS_reg
projection(s_raster) <- CRS_reg
projection(spat_pred_rast) <- CRS_reg #this should be CRS_interp
projection(temp_pred_rast) <- CRS_reg
projection(spat_pred_rast_mle_Chebyshev) <- CRS_reg
projection(spat_pred_rast_mle_eigen) <- CRS_reg
projection(r_temp_pred_rast_lm) <- CRS_reg

#r_temp_s <- r_temp_pred_rast_arima #Now temporal predicitons based on ARIMA!!!
#temp_pred_rast <- r_temp_s
#projection(temp_pred_rast) <- CRS_WGS84

#projection(spat_pred_rast) <- CRS_WGS84
#projection(spat_pred_rast_gmm) <- CRS_WGS84
#levelplot(spat_pred_rast_gmm,col.regions=matlab.like(25))

r_huric_w <- subset(s_raster,time_window_selected[-1])
#r_huric_w <- crop(r_huric_w,rast_ref)

#r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
res_temp_s <- temp_pred_rast - r_huric_w
res_spat_s <- spat_pred_rast - r_huric_w
#pred_spat_mle_Chebyshev_test <-load_obj("pred_spat_mle_chebyshev_EDGY_predictions_03092015.RData")

#r_res_s_1_temp__EDGY_predictions_03182015.tif
names(res_temp_s) <- sub("pred","res",names(res_temp_s))
#names(res_spat_mle_s) <- sub("pred","res",names(res_spat_mle_s))
names(res_spat_s) <- sub("pred","res",names(res_spat_s))
#names(res_temp_s) <- paste("r_res_s_",1:nlayers(res_temp_s),"_",out_suffix,sep="")
#debug(calc_ac_stat_fun)
#r_zones <- raster(l_rast[22])
#reproject data to latlong WGS84 (EPSG4326)
#r_winds <- raster(winds_zones_fname)
#projection(r_winds) <- proj_modis_str 
#r_winds_m <- projectRaster(from=r_winds,res_temp_s,method="ngb") #Check that it is using ngb
             
#r_in <-stack(l_rast)
#r_results <- stack(s_raster,temp_pred_rast,spat_pred_rast_mle,spat_pred_rast_gmm,res_temp_s,res_spat_mle_s,res_spat_gmm_s)
#r_results <- stack(s_raster,r_winds_m,temp_pred_rast,spat_pred_rast_gmm,res_temp_s,res_spat_gmm_s)
r_results <- stack(s_raster,rast_zonal,r_temp_pred_rast_lm,spat_pred_rast_mle_eigen,spat_pred_rast_mle_Chebyshev,res_temp_s,res_spat_s)

dat_out <- as.data.frame(r_results)
dat_out <- na.omit(dat_out)
write.table(dat_out,file=paste("dat_out_",out_suffix,".txt",sep=""),
            row.names=F,sep=",",col.names=T)

### Now accuracy assessment using MAE

out_suffix_s <- paste("temp_",out_suffix,sep="_")
#debug(calc_ac_stat_fun)
#projection(rast_ref) <- CRS_WGS84
#projection(spat_pred_rast_mle) <- CRS_WGS84
#projection(z_winds) <- CRS_WGS84 #making sure proj4 representation of projections are the same

undebug(calc_ac_stat_fun)
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast,r_var_s=r_huric_w,r_zones=rast_zonal,
                                file_format=file_format,out_suffix=out_suffix_s)
#out_suffix_s <- paste("spat_",out_suffix,sep="_")  
#ac_spat_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,r_var_s=r_huric_w,r_zones=rast_ref,
#                                file_format=file_format,out_suffix=out_suffix_s)

out_suffix_s <- paste("spat_mle",out_suffix,sep="_")  
ac_spat_mle_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,r_var_s=r_huric_w,r_zones=rast_zonal,
                                file_format=file_format,out_suffix=out_suffix_s)

#mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
#mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_mle_obj$mae_tb,ac_temp_obj$mae_tb))

mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg","temp")
#mae_tot_tb$time <- 2:nrow(mae_tot_tb)
mae_tot_tb$time <- 2:n_time

y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range)
lines(temp ~ time, type="b",col="magenta",data=mae_tot_tb)
legend("topleft",legend=c("spat","temp"),col=c("cyan","magenta"),lty=1)
title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!
write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))

#### BY ZONES ASSESSMENT: Ok now it is general so it should be part of the function...

#mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb[1:3,],
#                      ac_temp_obj$mae_zones_tb[1:3,])
mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb,
                      ac_temp_obj$mae_zones_tb)

mae_zones_tb <- as.data.frame(mae_zones_tb)

n_zones <- length(unique(mae_zones_tb$zone))

mae_zones_tb$method <- c(rep("spat_reg",n_zones),rep("temp_reg",n_zones))

n_time <- ncol(mae_zones_tb) -1
pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")
#names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")
names(mae_zones_tb) <- pred_names

write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
#drop first few rows that contain no data but zones...
n_start <-n_zones*2 +1
dd <- dd[n_start:nrow(dd),]
dd$zones <- mydata$zone #use recycle rule

xyplot(data~which |zones,group=method,data=dd,type="b",xlab="time",ylab="VAR",
       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2)
)


################### END OF SCRIPT ##################