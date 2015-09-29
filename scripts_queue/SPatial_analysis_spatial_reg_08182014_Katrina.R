####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 08/18/2014
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

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_08172014_functions.R"

script_path <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/R_workshop_WM_04232014" #path to script
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.


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
out_suffix <-"_predictions_08182014" #output suffix for the files that are masked for quality and for 
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
data_tb <-read.table(data_fname,sep=",",header=T)

#### Make this a function...

coord_names <- c("XCoord","YCoord")
l_rast <- rasterize_df_fun(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format,NA_flag_val)
#debug(rasterize_df_fun)

reg_var_list <- l_rast[6:18]
r_stack <- stack(reg_var_list)
plot(1:13,data_tb[600,6:18],type="l")

#r_spat_pred <- rasterize(data_reg_spdf,rast_ref,field="spat_reg_pred") #this is the prediction from lm model
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
plot(r_pred)

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