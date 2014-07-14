####################################    Space Time Analyses Project   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces predictions for the dates following the Hurricane Dean event.       
#The script uses spatial regression with weight matrix to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 07/17/2014
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to spatial regression with R 
#PROJECT: Space beats time paper 
#PROJECT: William and Mary Geospatial Pattern Analyses, William and Mary GIS Summer class
#################################################################################################

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

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_07142014_functions.R"

script_path <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/R_workshop_WM_04232014" #path to script
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
#Using source, reads the script containing the functions, loads functions in the workspace/enviroment
#making them available for the user.

#####  Parameters and argument set up ###########

in_dir<-"~/Data/Space_beats_time/R_Workshop_April2014"
#in_dir <-"/home/parmentier/Data/Space_Time"
#in_dir <- "/Users/benoitparmentier/Google Drive/R_Workshop_April2014"
#out_dir <-  "/Users/benoitparmentier/Google Drive/R_Workshop_April2014"

moore_window <- file.path(in_dir,"window_test4.rst") #this is the raster used to crop the input data
winds_zones_fname <- file.path(in_dir,"00_windzones_moore_sin.rst") #wind zone for the hurricane

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

file_format <- ".rst" #output raster format, all raster will be IDRISI format
NA_value <- "-9999" # No data values used in the raster images
out_suffix <-"_predictions_07142014" #output suffix for the files and directory created
create_out_dir_param=TRUE #if true, a new output directory is created...

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

## Read in data defining the area of study
moore_w <- raster(moore_window) #read raster file and create a RasterLayer object
moore_w <- moore_w > -9999 #reclass image (this creates a boolean image with 1 and 0)
rast_ref <- moore_w #create ref image, this is a smaller area than the NDVI images
projection(rast_ref) <- CRS_WGS84 #assign spatial reference informaiton
freq(rast_ref) #frequency of values in rast_ref
rast_ref[rast_ref==0] <- NA #assign NA for zero val

r_winds <- raster(winds_zones_fname) #read raster file about winds zone this is sinusoidal projection

projection(r_winds) <- proj_modis_str #assign projection

#list input NDVI time series images: This is a stack of 276 raster images...
reg_var_list <- list.files(path=file.path(in_dir,"moore_NDVI_wgs84"),
                           pattern="moore_reg.*.MOD13A2_A.*04062014.*.rst$",full.name=TRUE)
#moore_reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif

#Create a stack of layers
r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
levelplot(r_stack,layers=1:2) #show first 2 images (half a year more or less)
plot(r_stack,y=1:2)

################ PART II : RUN SPATIAL REGRESSION for one date ############

#Now mask and crop layer
r_var <- subset(r_stack,153:154) #date before hurricane is 153
r_clip <- rast_ref #the reference image is used for cropping the large area into a smaller area
r_var <- crop(r_var,rast_ref) #crop/clip image using another reference image

pix_id_r <- rast_ref #make a copy of the raster image
values(pix_id_r) <- 1:ncell(rast_ref) #create an image with pixel id for every observation
pix_id_r <- mask(pix_id_r,rast_ref) #3854 pixels from which 779 pixels are NA

### Convert raster to vector and calculate neighbours weights
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
names(data_reg)[1:2] <- c("NDVI_153","NDVI_154")
  
NDVI_selected="NDVI_153" #this is the variable modeled...
data_reg$NDVI <- data_reg[[NDVI_selected]] #selecte a specific column in the data.frame an assign it to a new column

## NOW RUN THE SPATIAL MODEL...

sam.esar <- errorsarlm(NDVI~1, listw=reg_listw_w, 
                      data=data_reg,na.action=na.omit,zero.policy=TRUE,
                      tol.solve=1e-36)
#This does not work because some pixels have no neighbours or neighbours with NA!!

### CLEAN OUT AND SCREEN NA and list of neighbours
#Let's use our function we created to clean out neighbours
#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
r_var <- subset(r_stack,153) #this is the date we want to use to run the spatial regression
r_clip <- rast_ref #this is the image defining the study area
proj_str<- NULL #SRS/CRS projection system
out_suffix_s <- paste("d153",out_suffix,sep="_") #suffix used
#out_dir and out_suffix set earlier

#Used the functions from the script sourced file to screen NA and create a neighbour object
nb_obj_for_pred_t_153 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)
names(nb_obj_for_pred_t_153) #show the structure of object (which is made up of a list)
r_poly_name <- nb_obj_for_pred_t_153$r_poly_name #name of the shapefile that has been cleaned out
reg_listw_w <- nb_obj_for_pred_t_153$r_listw #list of weights for cleaned out shapefile
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
#Predicted values and residuals
data_reg_spdf$spat_reg_pred <- sam.esar$fitted.values #store fitted/predicted values
data_reg_spdf$spat_reg_res <- sam.esar$residuals #store residuals values

################### PART III: RUN TEMPORAL MODEL USING LM for one date ########

#### NOW PREDICTION USING PREVIOUS TIME STEP 

r_var2 <- subset(r_stack,152:153)
r_var2_w <- crop(r_var2,rast_ref)

data_reg2_spdf <- as(r_var2_w,"SpatialPointsDataFrame") #convert raster image to spatial point data.frame
names(data_reg2_spdf) <- c("t1","t2") #assign name
data_reg2 <- as.data.frame(data_reg2_spdf)
data_reg2 <- na.omit(data_reg2) #remove NA...this reduces the number of observations
lm_mod <- lm(t2 ~ t1, data=data_reg2)
summary(lm_mod)

#Extract predicted values and residuals
data_reg2$lm_temp_pred <- lm_mod$fitted.values
data_reg2$lm_temp_res <- lm_mod$residuals

#Convert data.frame back to a spatial point data.frame
coordinates(data_reg2) <- c("x","y")
proj4string(data_reg2) <- CRS_WGS84

###########################################
############## PART IV: Produce images from the individual predictions using time and space ####

### NOW CREATE RASTER IMAGES FROM THE PREDICTIONS ...

#Convert the SpatialPolygonsDataFrame field containing the spatial prediction into a raster image
r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
#file_format <- ".rst"
raster_name <- paste("r_spat_pred","_",out_suffix,file_format,sep="")
writeRaster(r_spat_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)

#Convert the SpatialPolygonsDataFrame field containing the temporal prediction into a raster image
r_temp_pred <- rasterize(data_reg2,rast_ref,field="lm_temp_pred") #this is the prediction from lm model
#file_format <- ".rst"
raster_name <- paste("r_temp_pred","_",out_suffix,file_format,sep="") #this is the output name for the raster created
writeRaster(r_temp_pred,filename=file.path(out_dir,raster_name),overwrite=TRUE)

## This can be repeated in a loop to produce predictions and compare the actual to predicted
r_pred <- stack(r_spat_pred,r_temp_pred) #stack of predictions
names(r_pred) <- c("spatial pred","temporal pred") #change layerNames to names when using R 3.0 or above
plot(stack(r_spat_pred,r_temp_pred))
levelplot(r_pred,regions.col=rev(terrain.colors(255)),main="NDVI predictions after hurricane")

#### Examining difference between temporal and spatial predictions
r_dif <- r_spat_pred - r_temp_pred
plot(r_dif)
hist(r_dif)

##############################################################################################
############## PART V PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER 4 time steps ####

##This times will we use an automated function to generate predictions over 4 dates.

### Predict using spatial regression function, first prepare the input arguments in a list
r_spat_var <- subset(r_stack,152:155) #predict before (step 152) and three dates after (step 153)
list_models <- NULL #this parameter is not used
proj_str <- NULL #if null the raster images are not reprojected
out_suffix_s <- paste("t_",152:155,out_suffix,sep="") #use mle estimation, name "mle" is included in the output automatically?
estimator <- "mle" # we wish to use the Maximum likelihood estimator

list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator")
n_pred <- nlayers(r_spat_var)
## First predict for one date for testing and comparison
#debug(predict_spat_reg_fun)
#use mle estimator
testd1_spat_mle <- predict_spat_reg_fun(1,list_param_spat_reg) 

# We will not use GMM in the William and Mary lab
#use generalized methods of moment as estimator: won't work with version of 2.14 and old sphet package from Feb 2012
#list_param_spat_reg$estimator <- "gmm"
#use ols as estimator...this is added as a dictatic form for teaching, not to be used for research
#testd1_spat_gmm <- predict_spat_reg_fun(1,list_param_spat_reg)

list_param_spat_reg$estimator <- "ols" #use the OLS estimator now
testd1_spat_ols <- predict_spat_reg_fun(1,list_param_spat_reg) #check if it runs for one date

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

## Predict using temporal info: time steps 153,154,155,156

r_temp_var <- subset(r_stack,151:155) #need date 151 because it relies on the previous date in contrast to spat reg
list_models <-NULL      #this parameter is set to null
out_suffix_s <- paste("t_",152:155,out_suffix,sep="") #predicted time steps
list_param_temp_reg <- list(out_dir,r_temp_var,r_clip,proj_str,list_models,out_suffix_s,file_format)
names(list_param_temp_reg) <- c("out_dir","r_var","r_clip","proj_str","list_models","out_suffix_s","file_format")
n_pred <- nlayers(r_temp_var) -1
#debug(predict_temp_reg_fun)
#test_temp <- predict_temp_reg_fun(1,list_param_temp_reg)

pred_temp_lm <- lapply(1:n_pred,FUN=predict_temp_reg_fun,list_param=list_param_temp_reg)

temp_pred_rast <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_pred}))
temp_res_rast <- stack(lapply(pred_temp_lm,FUN=function(x){x$raster_res}))
levelplot(temp_pred_rast) #view the four predictions using mle spatial reg.

############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

r_huric_w <- subset(r_stack,152:155) #
r_huric_w <- crop(r_huric_w,rast_ref)

#r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
res_temp_s <- temp_pred_rast - r_huric_w #temporal residuals
res_spat_s <- spat_pred_rast - r_huric_w #spatial residuals

#Winds data is in MODIS sinusoidal projection so reproject data to latlong WGS84 (EPSG4326)
r_winds_m <- projectRaster(from=r_winds,res_temp_s,method="ngb") #Check that it is using ngb

##prepare output name using a suffix
out_suffix_s <- paste("temp_",out_suffix,sep="_")
#debug(calc_ac_stat_fun)
#Use function provided to compute the MAE for each time step
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast,r_var_s=r_huric_w,r_zones=r_winds_m,
                                file_format=file_format,out_suffix=out_suffix_s)
out_suffix_s <- paste("spat_",out_suffix,sep="_")  
ac_spat_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,r_var_s=r_huric_w,r_zones=r_winds_m,
                                file_format=file_format,out_suffix=out_suffix_s)

#mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg","temp")
mae_tot_tb$time <- 1:4

plot(spat_reg ~ time, type="b",col="magenta",data=mae_tot_tb,ylim=c(800,2000))
lines(temp ~ time, type="b",col="cyan",data=mae_tot_tb)
write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
legend("topleft",legend=c("spat","temp"),col=c("magenta","cyan"),lty=1)
title("Overall MAE for spatial and temporal models")

#### BY ZONES ASSESSMENT

mae_zones_tb <- rbind(ac_spat_obj$mae_zones_tb[2:3,],
                      ac_temp_obj$mae_zones_tb[2:3,])
mae_zones_tb <- as.data.frame(mae_zones_tb)
mae_zones_tb$method <- c("spat_reg","spat_reg","temp","temp")
names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")

write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

#Very quick and dirty plot
time <-1:4
x <- as.numeric(mae_zones_tb[1,2:5])
plot(x~time,, type="b",col="magenta",lty=1,ylim=c(400,2000),ylab="MAE for NDVI")
x <- as.numeric(mae_zones_tb[2,2:5])
lines(x~time, type="b",lty=2,col="magenta")
#add temporal
x <- as.numeric(mae_zones_tb[3,2:5]) #zone 4
lines(x~time,, type="b",col="cyan",lty=1,ylim=c(400,2000))
x <- as.numeric(mae_zones_tb[4,2:5]) #zone 5
lines(x~time, type="b",lty=2,col="cyan")
legend("topleft",legend=c("spat zone 4","spat zone 5","temp zone 4","temp zone 5"),
        col=c("magenta","magenta","cyan","cyan"),lty=c(1,2,1,2))
title("MAE per wind zones for spatial and temporal models")

################### END OF SCRIPT ##################