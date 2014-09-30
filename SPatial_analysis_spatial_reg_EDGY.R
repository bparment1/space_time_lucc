####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 09/30/2014
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
library(forecast)
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
library(sphet) #contains spreg

###### Functions used in this script

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_09252014_functions.R"
script_path <- "/home/parmentier/Data/Space_Time/R" #path to script
#script_path <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/R_workshop_WM_04232014" #path to script
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
out_suffix <-"EDGY_predictions_09302014" #output suffix for the files that are masked for quality and for 
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
#EDGY_dat_spdf_04072014.txt
#data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname <- file.path("/home/parmentier/Data/Space_Time","EDGY_dat_spdf_04072014.txt")

data_tb <-read.table(data_fname,sep=",",header=T)

#### Make this a function...

#Transform table text file into a raster image
#coord_names <- c("XCoord","YCoord") #for Katrina
coord_names <- c("r_x","r_y") #for EDGY

#l_rast <- rasterize_df_fun(data_tb[,1:5],coord_names,proj_str,out_suffix,out_dir=".",file_format,NA_flag_val)
l_rast <- rasterize_df_fun(data_tb,coord_names,proj_str,out_suffix,out_dir=".",file_format,NA_flag_val)

#debug(rasterize_df_fun)
s_raster <- stack(l_rast) #stack with all the variables
names(s_raster) <- names(data_tb)                

r_FID <- raster(l_rast[1])
plot(r_FID,main="Pixel ID")
freq(r_FID)

reg_var_list <- l_rast[6:281] #only select population raster
r_stack <- stack(reg_var_list)
names(r_stack) <- paste("NDVI",1:nlayers(r_stack),sep="_") #2000:2012

#debug(rasterize_df_fun)
#r_FID <- raster(l_rast[3])
#plot(r_FID,main="Pixel ID")
#freq(r_FID)
#reg_var_list <- l_rast[6:18] #only select population raster
#r_stack <- stack(reg_var_list)
#names(r_stack) <- 2000:2012

### Quick plots

#projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
levelplot(r_stack,layers=153:154,col.regions=matlab.like(125)) #show first 2 images (half a year more or less)
plot(r_stack,y=153:154)
histogram(subset(r_stack,153:154))

plot(1:276,data_tb[300,6:281],type="b", main="Population temporal  profile at pixel 300 and pixel 200 (red)",
     ylim=c(0,10000),ylab="Population",xlab="Year")#This is for line 600 in the table
abline(v=154,lty=3)
lines(1:276,data_tb[200,6:281],col="red",type="b", main="Population temporal  profile at pixel 600",
     ylab="Population",xlab="Year")#This is for line 600 in the table
#Plot only year 2007
plot(1:23,data_tb[8000,139:161],type="b", main="Population temporal  profile at pixel 300 and pixel 200 (red)",
     ylim=c(2000,10000),ylab="Population",xlab="Year")#This is for line 600 in the table
abline(v=16,lty=3)
lines(1:23,data_tb[20000,139:161],col="red",type="b", main="Population temporal  profile at pixel 600",
     ylab="Population",xlab="Year")#This is for line 600 in the table

data_tb$FID[40000] #FID 599
data_tb$FID[200] #FID 599

data_tb

## Read in data defining the area of study
#moore_w <- raster(moore_window) #read raster file and create a RasterLayer object
#don't use window...?
#moore_w <- moore_w > -9999 #reclass image (this creates a boolean image with 1 and 0)
#rast_ref <- moore_w #create ref image, this is a smaller area than the NDVI images
rast_ref <- subset(s_raster)
projection(rast_ref) <- CRS_WGS84 #assign spatial reference informaiton
freq(rast_ref) #frequency of values in rast_ref
rast_ref[rast_ref==0] <- NA #assign NA for zero val

r_winds <- raster(winds_zones_fname) #read raster file about winds zone this is sinusoidal projection

projection(r_winds) <- proj_modis_str #assign projection

#list input NDVI time series images: This is a stack of 276 raster images...
#ndvi_list <- list.files(path=file.path(in_dir,"moore_NDVI_wgs84"),
#                           pattern="moore_reg.*.MOD13A2_A.*04062014.*.rst$",full.name=TRUE)
#moore_reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif

#Create a stack of layers
r_stack <- stack(reg_var_list) #276 images from 2001 to 2012 included
projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
levelplot(r_stack,layers=153:154,col.regions=rev(terrain.colors(255))) #show first 2 images (half a year more or less)
plot(r_stack,y=1:2)

#projection(r_stack) <- CRS_WGS84 #making sure the same  CRS format is used
#levelplot(r_stack,layers=5:6,col.regions=matlab.like(125)) #show first 2 images (half a year more or less)
#plot(r_stack,y=5:6)

#plot(2000:2012,data_tb[600,6:18],type="b", main="Population temporal  profile at pixel 600 and pixel 300 (red)",
#     ylab="Population",xlab="Year")#This is for line 600 in the table
#abline(v=2005,lty=3)
#lines(2000:2012,data_tb[200,6:18],col="red",type="b", main="Population temporal  profile at pixel 600",
#     ylab="Population",xlab="Year")#This is for line 600 in the table

#data_tb$FID[600] #FID 599
#data_tb$FID[300] #FID 599

#data_tb
################ PART II : RUN SPATIAL REGRESSION ############

#Now mask and crop layer
r_var <- subset(r_stack,153:154) #date before hurricane is 153
#r_clip <- rast_ref
#r_var <- crop(r_var,rast_ref) #crop/clip image using another reference image
rast_ref <- subset(r_stack,1)
rast_ref <- rast_ref != NA_flag_val
pix_id_r <- rast_ref
values(pix_id_r) <- 1:ncell(rast_ref) #create an image with pixel id for every observation
pix_id_r <- mask(pix_id_r,rast_ref) #3854 pixels from which 779 pixels are NA

### CLEAN OUT AND SCREEN NA and list of neighbours
#Let's use our function we created to clean out neighbours
#Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
#this is for prediction t -1 
#list_timesteps <- 
#for(i in 1:length(list_timesteps)){
#r_var <- subset(r_stack,6) #this is the date we want to use to run the spatial regression
#r_clip <- rast_ref #this is the image defining the study area
#proj_str<- NULL #SRS/CRS projection system
#out_suffix_s <- paste("d2005",out_suffix,sep="_")

r_var <- subset(r_stack,154) #this is the date we want to use to run the spatial regression
r_clip <- rast_ref #this is the image defining the study area
proj_str<- NULL #SRS/CRS projection system
out_suffix_s <- paste("t_154",out_suffix,sep="_")
#out_dir and out_suffix set earlier

nb_obj_for_pred_t_154 <-create_sp_poly_spatial_reg(r_var,r_clip,proj_str,out_suffix=out_suffix_s,out_dir)
names(nb_obj_for_pred_t_154) #show the structure of object (which is made up of a list)
r_poly_name <- nb_obj_for_pred_t_154$r_poly_name #name of the shapefile that has been cleaned out
reg_listw_w <- nb_obj_for_pred_t_154$r_listw #list of weights for cleaned out shapefile
#Use OGR to load the screened out data: note that we have now 2858 features with 3 fields
data_reg_spdf <- readOGR(dsn=dirname(r_poly_name),
                    layer=gsub(extension(basename(r_poly_name)),"",basename(r_poly_name)))

data_reg <- as.data.frame(data_reg_spdf) #convert spdf to df

#Now run the spatial regression, since there are no covariates, the error and lag models are equivalent
v5 <- rnorm(nrow(data_reg))
data_reg$v5 <- v5 - mean(v5) 

#sam.esar3 <- errorsarlm(v1~ v5, listw=reg_listw_w, 
#                       data=data_reg,na.action=na.omit,zero.policy=TRUE,
#                       tol.solve=1e-36) #tol.solve use in matrix operations

#data_reg_spdf$spat_reg_pred <- sam.esar$fitted.values
#data_reg_spdf$spat_reg_res <- sam.esar$residuals

spat_mod_spreg3 <- try(spreg(v1 ~ v5,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))
data_reg_spdf$spat_reg_pred <- data_reg$v1 + spat_mod_spreg3$residuals
data_reg_spdf$spat_reg_res <- spat_mod_spreg3$residuals

#lm_mod <- try(lm(v1 ~ v5,data=data_reg))
#mean(data_reg$v5)
#all.equal(mean(data_reg$v2),as.numeric(coef(lm_mod)[1])) #ok this work...

################### PART III RUN TEMPORAL MODEL USING LM ########

#### NOW PREDICTION USING PREVIOUS TIME STEP 

r_var2 <- subset(r_stack,153:154)
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
levelplot(r_pred,col.regions=matlab.like(25),main="Var predictions after hurricane")

#### Examining difference between predictions
r_dif <- r_spat_pred - r_temp_pred
plot(r_dif,main="Difference between spatial and temporal models")
hist(r_dif,main="Difference between spatial and temporal models")

##############################################################################################
############## PART V PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER 4 time steps ####

##This times will we use an automated function to generate predictions over 4 dates

### Predict using spatial regression
#r_spat_var <- subset(r_stack,139:161) #predict before (step 152) and three dates after (step 153)
r_spat_var <- subset(r_stack,153:156) #predict before (step 152) and three dates after (step 153)

list_models <- NULL
proj_str <- NULL #if null the raster images are not reprojected
out_suffix_s <- paste("t_",153:156,out_suffix,sep="") #use mle estimation, name "mle" is included in the output automatically?
#estimator <- "mle"
estimator <- "mle"

list_param_spat_reg <- list(out_dir,r_spat_var,r_clip,proj_str,list_models,out_suffix_s,file_format,estimator)
names(list_param_spat_reg) <- c("out_dir","r_var_spat","r_clip","proj_str","list_models","out_suffix","file_format","estimator")
n_pred <- nlayers(r_spat_var)
## First predict for one date for testing and comparison
#debug(predict_spat_reg_fun)
#use mle estimator
#testd1_spat_mle <- predict_spat_reg_fun(1,list_param_spat_reg)
#use generalized methods of moment as estimator: won't work with version of 2.14 and old sphet package from Feb 2012
list_param_spat_reg$estimator <- "gmm"
#use ols as estimator...this is added as a dictatic form for teaching, not to be used for research
#debug(predict_spat_reg_fun)
testd1_spat_gmm <- predict_spat_reg_fun(1,list_param_spat_reg)

list_param_spat_reg$estimator <- "ols"
testd1_spat_ols <- predict_spat_reg_fun(1,list_param_spat_reg)

## Now predict for four dates using "mle"
#list_param_spat_reg$estimator <- "mle"
#pred_spat_mle <-lapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg)
#Use parallel processing on MAC and Linux/Unix systems
#pred_spat_mle <-mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,,mc.preschedule=FALSE,mc.cores = 2)
list_param_spat_reg$estimator <- "gmm"
#pred_spat_gmm <-lapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg)
pred_spat_gmm <-mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,mc.preschedule=FALSE,mc.cores = 4)

list_param_spat_reg$estimator <- "ols"
#pred_spat_ols <-lapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg)
pred_spat_ols <-mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,mc.preschedule=FALSE,mc.cores = 4)

#spat_pred_rast <- stack(lapply(pred_spat_mle,FUN=function(x){x$raster_pred})) #get stack of predicted images
#spat_res_rast <- stack(lapply(pred_spat_mle,FUN=function(x){x$raster_res})) #get stack of predicted images
#levelplot(spat_pred_rast) #view the four predictions using mle spatial reg.
#levelplot(spat_res_rast,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

spat_pred_rast_gmm <- stack(lapply(pred_spat_gmm,FUN=function(x){x$raster_pred})) #get stack of predicted images
spat_res_rast_gmm <- stack(lapply(pred_spat_gmm,FUN=function(x){x$raster_res})) #get stack of predicted images
levelplot(spat_pred_rast_gmm) #view the four predictions using mle spatial reg.
levelplot(spat_res_rast_gmm,col.regions=matlab.like(25)) #view the four predictions using mle spatial reg.

## Predict using temporal info: time steps 152..156

r_temp_var <- subset(r_stack,152:156) #need date 151 because it relies on the previous date in contrast to spat reg
list_models <-NULL
out_suffix_s <- paste("t_",152:156,out_suffix,sep="")
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
projection(spat_pred_rast_gmm) <- CRS_WGS84

## Extract spatial coefficients

#pred_spat_mle
#l_coef_mle <- lapply(pred_spat_mle,FUN=function(x){coef(x$spat_mod)})
#pred_spat_gmm
l_coef_gmm <- lapply(pred_spat_gmm,FUN=function(x){coef(x$spat_mod)})
#pred_spat_ols
l_coef_ols <- lapply(pred_spat_ols,FUN=function(x){coef(x$spat_mod)})

#tb_coef_mle <- as.data.frame(do.call(rbind,l_coef_mle))
#tb_coef_mle  <- tb_coef_mle[,c(2,1)]
#names(tb_coef_mle)<- c("(Intercept)","rho")
#tb_coef_mle$v2<- NA
#tb_coef_mle$time <- 2001:2012
#tb_coef_mle$method <- "mle"
tb_coef_gmm <- as.data.frame(t(do.call(cbind,(l_coef_gmm))))
tb_coef_gmm <- tb_coef_gmm[,c(1,3,2)]
names(tb_coef_gmm)<- c("(Intercept)","rho","v2")
tb_coef_gmm$time <- 1:4
tb_coef_gmm$method <- "gmm"             
tb_coef_ols <- as.data.frame(do.call(rbind,l_coef_ols))
tb_coef_ols <- tb_coef_ols[,c(1,2)]
names(tb_coef_ols)<- c("(Intercept)","rho")
tb_coef_ols$v2 <- NA                
tb_coef_ols$time <- 1:4                
tb_coef_ols$method <- "ols"                

#tb_coef_method <- rbind(tb_coef_mle,tb_coef_gmm,tb_coef_ols)
tb_coef_method <- rbind(tb_coef_gmm,tb_coef_ols)

xyplot(rho~time,groups=method,data=tb_coef_method,type="b",
                 auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2),
                main="Comparison of rho coefficients with different methods for 2001-2012"
)

write.table(tb_coef_method,paste("tb_coef_method",out_suffix,".txt",sep=""),row.names=F,col.names=T)                

############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

r_huric_w <- subset(r_stack,153:156)
#r_huric_w <- crop(r_huric_w,rast_ref)

#r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
res_temp_s <- temp_pred_rast - r_huric_w
res_spat_gmm_s <- spat_pred_rast_gmm - r_huric_w

names(res_temp_s) <- sub("predictions","residuals",names(res_temp_s))
#names(res_spat_mle_s) <- sub("predictions","residuals",names(res_spat_mle_s))
names(res_spat_gmm_s) <- sub("predictions","residuals",names(res_spat_gmm_s))

#debug(calc_ac_stat_fun)
projection(rast_ref) <- CRS_WGS84
#r_zones <- raster(l_rast[22])
z_winds <- subset(s_raster,"z_winds")
projection(r_winds) <- CRS_WGS84
#reproject data to latlong WGS84 (EPSG4326)
r_winds <- raster(winds_zones_fname)
projection(r_winds) <- proj_modis_str 
r_winds_m <- projectRaster(from=r_winds,res_temp_s,method="ngb") #Check that it is using ngb
             
#r_in <-stack(l_rast)
projection(s_raster) <- CRS_WGS84
#r_results <- stack(s_raster,temp_pred_rast,spat_pred_rast_mle,spat_pred_rast_gmm,res_temp_s,res_spat_mle_s,res_spat_gmm_s)
r_results <- stack(s_raster,r_winds_m,temp_pred_rast,spat_pred_rast_gmm,res_temp_s,res_spat_gmm_s)

dat_out <- as.data.frame(r_results)
dat_out <- na.omit(dat_out)
write.table(dat_out,file=paste("dat_out_",out_suffix,".txt",sep=""),
            row.names=F,sep=",",col.names=T)

### Now accuracy assessment using MAE

out_suffix_s <- paste("temp_",out_suffix,sep="_")
#debug(calc_ac_stat_fun)
projection(rast_ref) <- CRS_WGS84
ac_temp_obj <- calc_ac_stat_fun(r_pred_s=temp_pred_rast,r_var_s=r_huric_w,r_zones=z_winds,
                                file_format=file_format,out_suffix=out_suffix_s)
#out_suffix_s <- paste("spat_",out_suffix,sep="_")  
#ac_spat_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast,r_var_s=r_huric_w,r_zones=rast_ref,
#                                file_format=file_format,out_suffix=out_suffix_s)

out_suffix_s <- paste("spat_gmm",out_suffix,sep="_")  
ac_spat_gmm_obj <- calc_ac_stat_fun(r_pred_s=spat_pred_rast_gmm,r_var_s=r_huric_w,r_zones=z_winds,
                                file_format=file_format,out_suffix=out_suffix_s)

#mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
#mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_gmm_obj$mae_tb,ac_temp_obj$mae_tb))

mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg","temp")
mae_tot_tb$time <- 1:4

plot(spat_reg ~ time, type="b",col="magenta",data=mae_tot_tb,ylim=c(0,1800))
lines(temp ~ time, type="b",col="cyan",data=mae_tot_tb)
write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
legend("topleft",legend=c("spat","temp"),col=c("magenta","cyan"),lty=1)
title("Overall MAE for spatial and temporal models")


#### BY ZONES ASSESSMENT

mae_zones_tb <- rbind(ac_spat_gmm_obj$mae_zones_tb[1:3,],
                      ac_temp_obj$mae_zones_tb[1:3,])
mae_zones_tb <- as.data.frame(mae_zones_tb)
mae_zones_tb$method <- c("spat_reg","spat_reg","spat_reg","temp","temp","temp")
names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")

write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
#drop first four rows
dd <- dd[7:nrow(dd),]

xyplot(data~which |zones,group=method,data=dd,type="b",xlab="year",ylab="NDVI",
       strip = strip.custom(factor.levels=c("z3","z4","z5")),
      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                     border = FALSE, lines = TRUE,cex=1.2)
)

#Very quick and dirty plot
time <-1:4
x <- as.numeric(mae_zones_tb[1,3:5])
plot(x~time,, type="b",col="magenta",lty=1,ylim=c(400,2000),ylab="MAE for NDVI")
#x <- as.numeric(mae_zones_tb[2,2:5])
#lines(x~time, type="b",lty=2,col="magenta")
#add temporal
#x <- as.numeric(mae_zones_tb[3,2:5]) #zone 4
#lines(x~time,, type="b",col="cyan",lty=1,ylim=c(400,2000))
#x <- as.numeric(mae_zones_tb[4,2:5]) #zone 5
#lines(x~time, type="b",lty=2,col="cyan")
#legend("topleft",legend=c("spat zone 4","spat zone 5","temp zone 4","temp zone 5"),
#        col=c("magenta","magenta","cyan","cyan"),lty=c(1,2,1,2))
#title("MAE per wind zones for spatial and temporal models")

### Use ARIMA predions instead of lm temporal pred...

#...add here



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