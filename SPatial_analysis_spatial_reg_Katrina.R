####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc -PART 2  #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbours to predict NDVI values in the MOORE EDGY region. 
#Temporal predictions use OLS with the image of the previous time step rather than ARIMA.
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 03/09/2014 
#DATE MODIFIED: 09/25/2014
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
# library(BMS) #contains hex2bin and bin2hex
# library(bitops)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
# library(forecast)
# library(xts)
# library(zoo)
# library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
library(sphet) #contains spreg
library(reshape) #Data format and type transformation

###### Functions used in this script

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_09252014_functions.R"
script_path <- "/home/parmentier/Data/Space_Time/R" #path to script
# script_path <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/R_workshop_WM_04232014" #path to script

# source all functions used in this script 
source(file.path(script_path, function_spatial_regression_analyses))  

#####  Parameters and argument set up ###########

in_dir <- "~/Data/Space_beats_time/R_Workshop_April2014"
# in_dir <-"/home/parmentier/Data/Space_Time"
# in_dir <- "/Users/benoitparmentier/Google Drive/R_Workshop_April2014"
# out_dir <-  "/Users/benoitparmentier/Google Drive/R_Workshop_April2014"

moore_window <- file.path(in_dir ,"window_test4.rst")
winds_zones_fname <- file.path(in_dir, "00_windzones_moore_sin.rst")

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
# CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_st r<- CRS_WGS84
CRS_interp <- proj_modis_str

file_format <- ".rst"
NA_value <- -9999
NA_flag_val <- NA_value
out_suffix <- "_predictions_09252014" #output suffix for the files that are masked for quality and for 
create_out_dir_param = TRUE

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
# set up the working directory
# Create output directory

out_dir <- in_dir  # output will be created in the input dir
out_suffix_s <- out_suffix  # can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir, out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

# data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014",
#                         "Katrina_Output_CSV - Katrina_pop.csv")

data_fname <- file.path("~/Data/Space_beats_time/stu/Katrina/run2/csv", "Katrina2.csv")
# /Users/benoitparmentier/Google Drive/Space_beats_time/stu/Katrina/run2/csv

data_tb <- read.table(data_fname, sep = ",", header = T)
# data_tb$ezone_c[data_tb$Elev_zone 

attach(data_tb)
data_tb$ezone_c[Elev_Zone < 3] <- 1
data_tb$ezone_c[Elev_Zone == 3] <- NA  
data_tb$ezone_c[Elev_Zone == 4] <- 2  
detach(data_tb)

#### Make this a function...

# Transform table text file into a raster image
coord_names <- c("XCoord", "YCoord")
# coord_names <- c("POINT_X1","POINT_Y1")
l_rast <- rasterize_df_fun(data_tb, coord_names, proj_str, out_suffix, 
                           out_dir = ".", file_format, NA_flag_val, 
                           tolerance_val = 0.000120005)

# debug(rasterize_df_fun)
s_raster <- stack(l_rast)  # stack with all the variables
names(s_raster) <- names(data_tb)                

r_FID <- raster(l_rast[4])
plot(r_FID, main = "Pixel ID")
freq(r_FID)

reg_var_list <- l_rast[7:19]  # only select population raster
r_stack <- stack(reg_var_list)
names(r_stack) <- 2000:2012

# projection(r_stack) <- CRS_WGS84 
# making sure the same CRS format is used

# show first 2 images (half a year more or less)
levelplot(r_stack, layers = 5:6, col.regions = matlab.like(125))  
plot(r_stack, y = 5:6)
histogram(subset(r_stack, 5:6))

plot(2000:2012, data_tb[300 ,7:19], type = "b", 
     main = "Population temporal profile at pixel 300 and pixel 200 (red)",
     ylim = c(0, 1600), ylab = "Population", xlab = "Year")  # This is for line 600 in the table
abline(v=2005, lty=3)
lines(2000:2012, data_tb[200, 7:19], col = "red", type = "b", 
      main = "Population temporal profile at pixel 600",
      ylab = "Population", xlab = "Year")  # This is for line 600 in the table

data_tb$FID[300]  # FID 599
data_tb$FID[200]  # FID 599

data_tb

mean_vals <- colMeans(data_tb[,7:19], na.rm = T)
plot(2000:2012, mean_vals, type = "b", ylab = "pop", xlab = "year", 
     ylim = c(0,1600), main = "Average Population and pixel 200 profile")
lines(2000:2012, data_tb[200, 7:19], col = "red", type = "b", 
      main = "Population temporal  profile at pixel 600",
      ylab = "Population",xlab="Year")  # This is for line 600 in the table
legend("topright", legend = c("Overall average pop", "Pixel 200 pop"),
        col = c("black", "red"), lty = c(1,2))
# title("Average pop per elevation zones (observed data)")

data_tb

mean_by_zones <- aggregate(.~Elev_Zone, data = data_tb[,c(7:19,22)], mean)
plot(2000:2012, as.numeric(mean_by_zones[1,2:14]), type = "b",
     col = "black", lty = 1, ylim = c(0, 2000),
     ylab = "pop", xlab = "year")
lines(2000:2012,as.numeric(mean_by_zones[2, 2:14]), type = "b", 
      col = "red", lty = 2)
lines(2000:2012,as.numeric(mean_by_zones[3, 2:14]), type = "b", 
      col = "blue", lty = 3)
lines(2000:2012,as.numeric(mean_by_zones[4, 2:14]), type = "b", 
      col = "green", lty = 4)

legend("topright", legend =c ("zone 1", "zone 2", "zone 3", "zone 4"),
        col = c("black", "red", "blue", "green"), lty = c(1, 2, 3, 4))
title("Average pop per elevation zones (observed data)")

# get mean for all column of data.frame!!
mean_by_ezone_c <- aggregate(.~ezone_c, data = data_tb[,c(7:19, 23)], mean) 
plot(2000:2012, as.numeric(mean_by_ezone_c[1, 2:14]), type = "b",
     col = "black", lty = 1, ylim = c(0, 2000),
     ylab = "pop", xlab = "year")
lines(2000:2012, as.numeric(mean_by_ezone_c[2, 2:14]), type = "b", 
      col = "red", lty = 2)

legend("topright", legend = c("zone 1: low", "zone 2: high"),
        col = c("black", "red", "blue", "green"), lty = c(1,2,3,4))
title("Average pop per elevation ezones c (observed data)")

elev_s <- subset(s_raster, c("Elev_Zone", "ezone_c"))
levelplot(elev_s)             

################ PART II : RUN SPATIAL REGRESSION ############

# Now mask and crop layer
r_var <- subset(r_stack, 5:6)  # date before hurricane is 153
# r_clip <- rast_ref
# r_var <- crop(r_var,rast_ref)  # crop/clip image using another reference image
rast_ref <- subset(r_stack, 1)
rast_ref <- rast_ref != NA_flag_val
pix_id_r <- rast_ref
values(pix_id_r) <- 1:ncell(rast_ref)  # create an image with pixel id for every observation
pix_id_r <- mask(pix_id_r, rast_ref)  # 3854 pixels from which 779 pixels are NA

########
### TEST SPATIAL Prediction on one date....
                
### CLEAN OUT AND SCREEN NA and list of neighbours
# Let's use our function we created to clean out neighbours
# Prepare dataset 1 for function: date t-2 (mt2) and t-1 (mt1) before hurricane
# this is for prediction t -1 
r_var <- subset(r_stack,6)  # this is the date we want to use to run the spatial regression
r_clip <- rast_ref  # this is the image defining the study area
proj_str<- NULL  # SRS/CRS projection system
out_suffix_s <- paste("d2005", out_suffix, sep = "_")
# out_dir and out_suffix set earlier

nb_obj_for_pred_t_2005 <- create_sp_poly_spatial_reg(r_var, r_clip, proj_str,
                                                     out_suffix = out_suffix_s, out_dir)
names(nb_obj_for_pred_t_2005)  # show the structure of object (which is made up of a list)
r_poly_name <- nb_obj_for_pred_t_2005$r_poly_name  # name of the shapefile that has been cleaned out
reg_listw_w <- nb_obj_for_pred_t_2005$r_listw  # list of weights for cleaned out shapefile
# Use OGR to load the screened out data: note that we have now 2858 features with 3 fields
data_reg_spdf <- readOGR(dsn = dirname(r_poly_name),
                    layer = gsub(extension(basename(r_poly_name)), "", 
                                 basename(r_poly_name)))

data_reg <- as.data.frame(data_reg_spdf)  # convert spdf to df
# test a fake covariate
data_reg$v4 <- 1

# Now run the spatial regression, since there are no covariates, the error and lag models are equivalent
# This takes a bit less than 3minutes for the dataset containing 2858 polygons
# tol.solve use in matrix operations
sam.esar <- errorsarlm(v1~1, listw = reg_listw_w, data = data_reg, 
                       na.action = na.omit, zero.policy = TRUE, tol.solve=1e-36)  
summary(sam.esar)
v5 <- rnorm(nrow(data_reg))
data_reg$v5 <- v5 - mean(v5) 
sam.esar2 <- errorsarlm(v1~ v4, listw = reg_listw_w, 
                       data = data_reg, na.action = na.omit, zero.policy = TRUE,
                       tol.solve = 1e-36)  
summary(sam.esar)

# data_reg$v5 <- rnorm(nrow(data_reg))
sam.esar3 <- errorsarlm(v1~ v5, listw = reg_listw_w, 
                       data = data_reg, na.action = na.omit, zero.policy = TRUE,
                       tol.solve = 1e-36) 
print(coef(sam.esar))
# almost equal could set the seed to make it reproduction
print(coef(sam.esar3)) 
# Predicted values 
data_reg_spdf$spat_reg_pred <- sam.esar$fitted.values
data_reg_spdf$spat_reg_res <- sam.esar$residuals

spat_mod <- try(spautolm(v1 ~ 1,data = data_reg, listw = reg_listw_w, 
                         na.action = na.omit, zero.policy = TRUE, tol.solve = 1e-36))
spat_mod_spreg <- try(spreg(v1 ~ v2, data = data_reg, listw = reg_listw_w, 
                            model = "error", het = TRUE, verbose = TRUE))

# res<-gstslshet(v1 ~ v2 , data=data_reg, listw=reg_listw_w)
spat_mod_spreg2 <- try(spreg(v1 ~ 1,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))

# Should we the var with mean removed?                
data_reg$v3 <- data_reg$v1 - mean(data_reg$v1)
spat_mod_spreg2 <- try(spreg(v1 ~ v4,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))

spat_mod_spreg3 <- try(spreg(v1 ~ v5,data=data_reg, listw= reg_listw_w, model="error",   
                     het = TRUE, verbose=TRUE))

lm_mod <- try(lm(v1 ~ v5,data=data_reg))
lm_mod2 <- try(lm(v1 ~ v4,data=data_reg))
mean(data_reg$v2)
mean(data_reg$v5)

all.equal(mean(data_reg$v2),as.numeric(coef(lm_mod)[1]))  # ok this work...
# A work around may be to include a standar normal random variable!!

################### PART III RUN TEMPORAL MODEL USING LM ########

#### NOW PREDICTION USING PREVIOUS TIME STEP 

r_var2 <- subset(r_stack, 5:6)
r_var2_w <- crop(r_var2, rast_ref)

data_reg2_spdf <- as(r_var2_w, "SpatialPointsDataFrame")
names(data_reg2_spdf) <- c("t1", "t2")
data_reg2 <- as.data.frame(data_reg2_spdf)
data_reg2 <- na.omit(data_reg2)  # remove NA...this reduces the number of observations
lm_mod <- lm(t2 ~ t1, data = data_reg2)
summary(lm_mod)

# Predicted values and
data_reg2$lm_temp_pred <- lm_mod$fitted.values
data_reg2$lm_temp_res <- lm_mod$residuals

coordinates(data_reg2) <- c("x", "y")
proj4string(data_reg2) <- CRS_WGS84

########################################### 
############## PART IV:Produce images from the individual predictions using time and space ####

### NOW CREATE THE IMAGES BACK ...

r_spat_pred <- rasterize(data_reg_spdf, rast_ref, field = "spat_reg_pred")  # this is the prediction from lm model
# file_format <- ".rst"
raster_name <- paste("r_spat_pred", "_", out_suffix, file_format, sep = "")
writeRaster(r_spat_pred, filename = file.path(out_dir, raster_name), overwrite = TRUE)

r_temp_pred <- rasterize(data_reg2, rast_ref, field = "lm_temp_pred")  # this is the prediction from lm model
# file_format <- ".rst"
raster_name <- paste("r_temp_pred", "_", out_suffix, file_format, sep = "")
writeRaster(r_temp_pred, filename = file.path(out_dir, raster_name), overwrite = TRUE)

## This can be repeated in a loop to produce predictions and compare the actual to predicted
r_pred <- stack(r_spat_pred, r_temp_pred)  # stack of predictions
# change layerNames to names when using R 3.0 or above
names(r_pred) <- c("spatial pred", "temporal pred")  
plot(r_pred)

levelplot(r_pred, regions.col = rev(terrain.colors(255)), main = "Var predictions after hurricane")
levelplot(r_pred, col.regions = matlab.like(25), main = "Var predictions after hurricane")

#### Examining difference between predictions
r_dif <- r_spat_pred - r_temp_pred
plot(r_dif, main = "Difference between spatial and temporal models")
hist(r_dif, main = "Difference between spatial and temporal models")

##############################################################################################
############## PART V PREDICT MODELS FOR USING TEMP AND SPAT REGRESSION OVER 4 time steps ####

##This times will we use an automated function to generate predictions over 4 dates

### Predict using spatial regression
# predict before (step 152) and three dates after (step 153)
r_spat_var <- subset(r_stack, 2:13)  
list_models <- NULL
proj_str <- NULL  # if null the raster images are not reprojected
out_suffix_s <- paste("t_", 2001:2012, out_suffix, sep = "") 
# use mle estimation, name "mle" is included in the output automatically?
estimator <- "mle"

list_param_spat_reg <- list(out_dir ,r_spat_var, r_clip, proj_str, list_models,
                            out_suffix_s, file_format, estimator)
names(list_param_spat_reg) <- c("out_dir", "r_var_spat", "r_clip", "proj_str", 
                                "list_models", "out_suffix", "file_format", "estimator")
n_pred <- nlayers(r_spat_var)
## First predict for one date for testing and comparison
# debug(predict_spat_reg_fun)
# use mle estimator
testd1_spat_mle <- predict_spat_reg_fun(1, list_param_spat_reg)
# use generalized methods of moment as estimator
# won't work with version of 2.14 and old sphet package from Feb 2012
list_param_spat_reg$estimator <- "gmm"
# use ols as estimator...this is added as a dictatic form for teaching, not to be used for research
# debug(predict_spat_reg_fun)
testd1_spat_gmm <- predict_spat_reg_fun(1, list_param_spat_reg)
testd1_spat_gmm <- predict_spat_reg_fun(5, list_param_spat_reg) #2005

# debug(predict_spat_reg_fun)                
list_param_spat_reg$estimator <- "ols"
testd1_spat_ols <- predict_spat_reg_fun(1, list_param_spat_reg)

## Now predict for four dates using "mle"
list_param_spat_reg$estimator <- "mle"
pred_spat_mle <- lapply(1:n_pred, FUN = predict_spat_reg_fun, 
                       list_param = list_param_spat_reg)
# Use parallel processing on MAC and Linux/Unix systems
# pred_spat_mle <-mclapply(1:n_pred,FUN=predict_spat_reg_fun,
#                          list_param=list_param_spat_reg,,
#                           mc.preschedule=FALSE,mc.cores = 2)

list_param_spat_reg$estimator <- "gmm"
pred_spat_gmm <- lapply(1:n_pred, FUN = predict_spat_reg_fun, 
                       list_param = list_param_spat_reg)

list_param_spat_reg$estimator <- "ols"
pred_spat_ols <- lapply(1:n_pred, FUN = predict_spat_reg_fun, 
                       list_param = list_param_spat_reg)

# get stack of predicted images
spat_pred_rast_mle <- stack(lapply(pred_spat_mle, FUN = function(x){x$raster_pred}))
spat_res_rast_mle <- stack(lapply(pred_spat_mle, FUN = function(x){x$raster_res})) 
# view the four predictions using mle spatial reg.
levelplot(spat_pred_rast_mle)  
levelplot(spat_res_rast_mle, col.regions = matlab.like(25))  

#get stack of predicted images
spat_pred_rast_gmm <- stack(lapply(pred_spat_gmm, FUN = function(x){x$raster_pred})) 
spat_res_rast_gmm <- stack(lapply(pred_spat_gmm, FUN = function(x){x$raster_res}))
#view the four predictions using mle spatial reg.
levelplot(spat_pred_rast_gmm) 
levelplot(spat_res_rast_gmm, col.regions = matlab.like(25)) 

## Predict using temporal info: time steps 2005..2012

# need date 151 because it relies on the previous date in contrast to spat reg
r_temp_var <- subset(r_stack, 1:13)  
list_models <- NULL
out_suffix_s <- paste("t_", 2001:2012, out_suffix, sep = "")
list_param_temp_reg <- list(out_dir, r_temp_var, r_clip, proj_str, list_models, 
                            out_suffix_s, file_format)
names(list_param_temp_reg) <- c("out_dir", "r_var", "r_clip", "proj_str", 
                                "list_models", "out_suffix_s", "file_format")
n_pred <- nlayers(r_temp_var) - 1
# debug(predict_temp_reg_fun)
# test_temp <- predict_temp_reg_fun(1,list_param_temp_reg)

pred_temp_lm <- lapply(1:n_pred, FUN = predict_temp_reg_fun, 
                       list_param = list_param_temp_reg)

temp_pred_rast <- stack(lapply(pred_temp_lm, FUN = function(x){x$raster_pred}))
temp_res_rast <- stack(lapply(pred_temp_lm, FUN = function(x){x$raster_res}))
levelplot(temp_pred_rast)  # view the four predictions using mle spatial reg
projection(temp_pred_rast) <- CRS_WGS84
projection(spat_pred_rast_mle) <- CRS_WGS84
projection(spat_pred_rast_gmm) <- CRS_WGS84

## Extract spatial coefficients

# pred_spat_mle
l_coef_mle <- lapply(pred_spat_mle, FUN=function(x){coef(x$spat_mod)})
# pred_spat_gmm
l_coef_gmm <- lapply(pred_spat_gmm, FUN=function(x){coef(x$spat_mod)})
# pred_spat_ols
l_coef_ols <- lapply(pred_spat_ols, FUN=function(x){coef(x$spat_mod)})

tb_coef_mle <- as.data.frame(do.call(rbind, l_coef_mle))
tb_coef_mle  <- tb_coef_mle[,c(2,1)]
names(tb_coef_mle)<- c("(Intercept)", "rho")
tb_coef_mle$v2<- NA
tb_coef_mle$time <- 2001:2012
tb_coef_mle$method <- "mle"
tb_coef_gmm <- as.data.frame(t(do.call(cbind, (l_coef_gmm))))
tb_coef_gmm <- tb_coef_gmm[,c(1,3,2)]
names(tb_coef_gmm)<- c("(Intercept)", "rho", "v2")
tb_coef_gmm$time <- 2001:2012
tb_coef_gmm$method <- "gmm"             
tb_coef_ols <- as.data.frame(do.call(rbind, l_coef_ols))
tb_coef_ols <- tb_coef_ols[,c(1,2)]
names(tb_coef_ols)<- c("(Intercept)", "rho")
tb_coef_ols$v2 <- NA                
tb_coef_ols$time <- 2001:2012                
tb_coef_ols$method <- "ols"                

tb_coef_method <- rbind(tb_coef_mle, tb_coef_gmm, tb_coef_ols)

xyplot(rho~time,groups=method,data=tb_coef_method,type="b",
       auto.key = list("topright", corner = c(0,1),  # col=c("black","red"),
       border = FALSE, lines = TRUE, cex = 1.2),
       main="Comparison of rho coefficients with different methods for 2001-2012")

write.table(tb_coef_method,paste("tb_coef_method", out_suffix, ".txt", sep = ""),
            row.names = F, col.names = T)                

############ PART V COMPARE MODELS IN PREDICTION ACCURACY #################

r_huric_w <- subset(r_stack,2:13)
r_huric_w <- crop(r_huric_w, rast_ref)

# reproject data to latlong WGS84 (EPSG4326)
# r_winds_m <- projectRaster(from=r_winds,res_temp_s,method="ngb") 
# Check that it is using ngb

# r_winds_m <- crop(winds_wgs84,res_temp_s) #small test window
res_temp_s <- temp_pred_rast - r_huric_w
res_spat_mle_s <- spat_pred_rast_mle - r_huric_w
res_spat_gmm_s <- spat_pred_rast_gmm  - r_huric_w
names(temp_pred_rast)
names(spat_pred_rast_mle)
names(res_temp_s) <- sub("predictions", "residuals", names(res_temp_s))
names(res_spat_mle_s) <- sub("predictions", "residuals", names(res_spat_mle_s))
names(res_spat_gmm_s) <- sub("predictions", "residuals", names(res_spat_gmm_s))

# debug(calc_ac_stat_fun)
projection(rast_ref) <- CRS_WGS84
# r_zones <- raster(l_rast[22])
r_elev_zones <- subset(s_raster, "Elev_Zone")
projection(r_elev_zones) <- CRS_WGS84

r_ezones_c <- subset(s_raster, "ezone_c")
projection(r_ezones_c) <- CRS_WGS84
                
# r_in <-stack(l_rast)
projection(s_raster) <- CRS_WGS84
r_results <- stack(s_raster, temp_pred_rast, spat_pred_rast_mle, spat_pred_rast_gmm, 
                   res_temp_s, res_spat_mle_s, res_spat_gmm_s)

dat_out <- as.data.frame(r_results)
dat_out <- na.omit(dat_out)
write.table(dat_out,file=paste("dat_out_", out_suffix, ".txt", sep = ""),
            row.names = F, sep = ",", col.names = T)

# Use elevation categories  as zones
out_suffix_s <- paste("temp_", out_suffix, sep = "_")
ac_temp_obj <- calc_ac_stat_fun(r_pred_s = temp_pred_rast, r_var_s = r_huric_w,
                                r_zones = r_elev_zones, file_format = file_format, 
                                out_suffix = out_suffix_s)
out_suffix_s <- paste("spat_", out_suffix, sep = "_")  
ac_spat_obj <- calc_ac_stat_fun(r_pred_s = spat_pred_rast, r_var_s = r_huric_w, 
                                r_zones = r_elev_zones, file_format = file_format, 
                                out_suffix = out_suffix_s)
                
# Use modified elevation categories  as zones: ezones c                
out_suffix_s <- paste("temp2_", out_suffix, sep = "_")
ac_temp_obj2 <- calc_ac_stat_fun(r_pred_s = temp_pred_rast, r_var_s = r_huric_w, 
                                 r_zones = r_ezones_c, file_format = file_format,
                                 out_suffix = out_suffix_s)
out_suffix_s <- paste("spat2_", out_suffix, sep = "_")  
ac_spat_obj2 <- calc_ac_stat_fun(r_pred_s = spat_pred_rast, r_var_s = r_huric_w, 
                                 r_zones = r_ezones_c, file_format = file_format,
                                 out_suffix = out_suffix_s)

out_suffix_s <- paste("spat_gmm2", out_suffix,sep="_")  
ac_spat_gmm_obj2 <- calc_ac_stat_fun(r_pred_s = spat_pred_rast_gmm, r_var_s = r_huric_w, 
                                     r_zones = r_ezones_c, file_format = file_format,
                                     out_suffix = out_suffix_s)

# mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c("spat_reg", "temp")
mae_tot_tb$time <- 1:12

plot(spat_reg ~ time, type = "b",col = "magenta", data = mae_tot_tb,
     ylim = c(0,1000), ylab = "MAE", xlab = "time step")
lines(temp ~ time, type = "b",col = "cyan", data = mae_tot_tb)
legend("topright" ,legend = c("spat","temp"), col = c("magenta","cyan"), lty = 1)
title("Overall MAE for spatial and temporal models")

write.table(mae_tot_tb, file = paste("mae_tot_tb", "_", out_suffix, ".txt", 
                                     sep = ""))

# mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
mae_tot_tb2 <- (cbind(ac_spat_obj2$mae_tb, ac_temp_obj2$mae_tb))
mae_tot_tb2 <- as.data.frame(mae_tot_tb2)
row.names(mae_tot_tb2) <- NULL
names(mae_tot_tb2)<- c("spat_reg", "temp")
mae_tot_tb2$time <- 1:12

plot(spat_reg ~ time, type = "b",col = "magenta", data = mae_tot_tb2,
     ylim = c(0,1000), ylab = "MAE", xlab = "time step")
lines(temp ~ time, type = "b",col = "cyan", data = mae_tot_tb2)
legend("topright", legend = c("spat","temp"), col = c("magenta","cyan"),lty = 1)
title("Overall MAE for spatial and temporal models")

write.table(mae_tot_tb2, file = paste("mae_tot_tb2", "_", out_suffix, ".txt",
                                      sep = ""))

##Now add GMM

# mae_tot_tb <- t(rbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb))
# gmm_mae_tot_tb <- (cbind(ac_spat_obj$mae_tb,ac_temp_obj$mae_tb,ac_spat_gmm_obj$mae_tb))
# mae_tot_tb <- as.data.frame(gmm_mae_tot_tb)
# row.names(mae_tot_tb) <- NULL
# names(mae_tot_tb)<- c("spat_reg","temp","gmm")
# mae_tot_tb$time <- 1:12

# plot(spat_reg ~ time, type="b",col="magenta",data=mae_tot_tb,ylim=c(0,1000))
# lines(temp ~ time, type="b",col="cyan",data=mae_tot_tb)
# lines(gmm ~ time, type="b",col="cyan",data=mae_tot_tb)

# write.table(mae_tot_tb,file=paste("mae_tot_tb","_",out_suffix,".txt",sep=""))
# legend("topleft",legend=c("spat","temp"),col=c("magenta","cyan"),lty=1)
# title("Overall MAE for spatial and temporal models")

#### BY ZONES ASSESSMENT

mae_zones_tb <- rbind(ac_spat_obj$mae_zones_tb[1:4,],
                      ac_temp_obj$mae_zones_tb[1:4,])
mae_zones_tb <- as.data.frame(mae_zones_tb)
mae_zones_tb$method <- c("spat_reg", "spat_reg", "spat_reg", "spat_reg",
                         "temp", "temp", "temp", "temp")

# names(mae_zones_tb) <- c("zones",paste("d",2001:2012,sep="_"),"method")
names(mae_zones_tb) <- c("zones", 2001:2012, "method")

write.table(mae_zones_tb, file=paste("mae_zones_tb", "_", out_suffix, ".txt", sep = ""))

mae_zones_tb2 <- rbind(ac_spat_obj2$mae_zones_tb[1:2,],
                      ac_temp_obj2$mae_zones_tb[1:2,])
mae_zones_tb2 <- as.data.frame(mae_zones_tb2)
mae_zones_tb2$method <- c("spat_reg", "spat_reg", "temp", "temp")

# names(mae_zones_tb) <- c("zones",paste("d",2001:2012,sep="_"),"method")
names(mae_zones_tb2) <- c("zones", 2001:2012, "method")

#prepare to automate the plotting of   all columns

### Use Elev Zone
                
# Should make this a function...
mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
# dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
# drop first four rows
dd <- dd[9:nrow(dd),]

xyplot(data~which | zones, group=method, data = dd, type = "b", xlab = "year", 
       ylab = "pop", strip = strip.custom(factor.levels = c("z1","z2","z3","z4")))

# Should make this a function...USING EZONE_C!! RECLASSIFIED DATA
mydata2<- mae_zones_tb2
dd2 <- do.call(make.groups, mydata2[,-ncol(mydata2)]) 
# dd$lag <- mydata$lag 
dd2$zones <- mydata2$zones
dd2$method <- mydata2$method
# drop first four rows
dd2 <- dd2[5:nrow(dd2),]

xyplot(data~which | zones, group = method, data = dd2, type = "b", xlab = "year",
       ylab = "pop", strip = strip.custom(factor.levels=c("z1","z2")))
#             strip=strip.custom(factor.levels=names_panel_plot),

# xyplot(data~which |zones+method,dd,type="b",xlab="year",ylab="lab")
# xyplot(data~which,groups=zones+method,dd,type="b")
# xyplot(data~which, groups=method+as.character(zones),dd,type="b")

plot(1:12, as.numeric(mae_zones_tb[1,2:13]), type = "b", col = "black", 
     ylim = c(0,1000))
lines(1:12, as.numeric(mae_zones_tb[2,2:13]), type = "b", col = "red")
lines(1:12, as.numeric(mae_zones_tb[3,2:13]), type = "b", col = "blue")
lines(1:12 ,as.numeric(mae_zones_tb[4,2:13]), type = "b", col = "green")

# Zone 1: spatial and temp
plot(1:12, as.numeric(mae_zones_tb[1,2:13]), type = "b", col = "magenta",
     ylim=c(0,1000), lty=1, ylab = "pop", xlab = "year")
lines(1:12,as.numeric(mae_zones_tb[5,2:13]), type = "b", col ="cyan",
      lty=2)
legend("topleft", legend = c("spat zone 1", "temp zone 1"),
        col = c("magenta", "cyan"), lty=c(1,2))

# Zone 2: spatial and temp
plot(1:12,as.numeric(mae_zones_tb[2,2:13]), type = "b", col = "magenta",
     ylim=c(0,1000), lty=1, xlab="year", ylab = "year")
lines(1:12,as.numeric(mae_zones_tb[6,2:13]), type = "b",col = "cyan",
      lty = 2)
legend("topleft", legend=c("spat zone 2", "temp zone 2"),
        col = c("magenta","cyan"), lty = c(1,2))

# Zone 3: spatial and temp
plot(1:12,as.numeric(mae_zones_tb[3,2:13]), type = "b", col = "magenta",
     ylim = c(0,1000), lty = 1)
lines(1:12,as.numeric(mae_zones_tb[7,2:13]), type = "b", col = "cyan",
      lty = 2)
legend("topleft", legend = c("spat zone 3","temp zone 3"),
        col = c("magenta","cyan"), lty = c(1,2))

# Zone 4: spatial and temp
plot(1:12, as.numeric(mae_zones_tb[4,2:13]), type = "b", col = "magenta",
     ylim = c(0,1000), lty = 1)
lines(1:12,as.numeric(mae_zones_tb[8,2:13]), type = "b", col = "cyan",
      lty = 2)
legend("topleft", legend = c("spat zone 4","temp zone 4"),
        col = c("magenta","cyan"), lty = c(1,2))

# Very quick and dirty plot
# time <-1:12
# x <- as.numeric(mae_zones_tb[1,2:5])

#mod_pat<-glob2rx("d_*")
#var_pat<-grep(mod_pat,names(mae_zones_tb),value=TRUE) # using grep with "value" extracts the matching names
#var_pat <- names(mae_zones_tb)
#tb_melt<-melt(mae_zones_tb,
#              measure=var_pat,
#              id=c("zones","method"),
#              na.rm=F)

#tb_melt<-melt(mae_zones_tb,
#              id=var_pat,
#              measure=var_pat,
#              na.rm=F)


#plot(x~time,, type="b",col="magenta",lty=1,ylim=c(400,2000),ylab="MAE for NDVI")
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

#dat_out__predictions_09092014.txt

#pdf("SampleGraph.pdf",width=7,height=5)
#plot(r_in,1)

#dev.off()

################### END OF SCRIPT ##################