####################################    Space Time Analyses PAPER   #######################################
###########################################  Figures production  #######################################
#This script produces figures based on outputs from earlier analyses for the space beats time framework.

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/20/2015 
#DATE MODIFIED: 02/23/2018
#Version: 1
#PROJECT: SBT framework - Book chapter with Rita results
#COMMENTS: 
#         
#COMMIT: changes to output figures of temporal profiles
#TO DO:
#produce jpegs with the tfw files...so that they can be opened in a transparency mode in a GIS directly overlaid.
#
#COMMITS: figures for RITA results
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

###### Functions used in this script sourced from other files

#function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_11242015_functions.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_02102018.R" #PARAM 1
function_space_and_time_assessment <- "space_and_time_assessment_functions_02232018b.R" #PARAM 1

script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts" #path on bpy50 #PARAM 2
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path on Atlas
#source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_space_and_time_assessment)) #source all functions used in this script 1.

#zonal stratum for NLU
##### Functions used in this script 

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

#####  Parameters and argument set up ###########

#ARG1: input dir
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time" 
#ARG2: output dir
out_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs" #bpy50 laptop
#ARG3: proj4 coordinates system info for modis projection
proj_str <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#ARG
file_format <- ".rst" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"sbt_book_figures_02202018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
datafname <- NULL #Need to update this
coord_names <- "x;y" #PARAM 9
zonal_colnames <- "r_zonal_rev" #PARAM 12
var_names <- "1;230" #PARAM 10 #Data is stored in the columns 3 to 22
num_cores <- "4" #PARAM 11
n_time_event <- "110;7;7" #PARAM 12 #timestep corresponding to the event for obs,spat pred, temp pred 
time_window_selected <- "105;114" #PARAM 13: use alll dates for now
previous_step <- TRUE #PARAM 14
date_range <- "2001.01.01;2010.12.31" #date
#Closest date to the even for each example:
date_event <- "2005-09-24" #Hurricane RITA


cat_name <- NULL
r_ref <- NULL
temp_fname <- NULL #Need to update this
spat_fname <- NULL #Need to update this
method_space <- "mle;eigen" #method for space and time used to predict space and time respectively
method_time <- "arima;arima;TRUE"
pixel_index <- NULL
mosaic_dir <- NULL
mosaic_out_suffix <- NULL

###########
date_range_str <- unlist(strsplit(date_range,";"))
dates <- generate_dates_by_step(date_range_str[1],date_range_str[2],16)$dates

#method_space <- df_args[19,index_val] 
#method_time <- df_args[20,index_val] 
#take a look at:
#https://github.com/bparment1/space_time_lucc/blob/master/space_and_time_assessment.R

#### Clean this up: Need to this general for any region!!!
#Latest relevant folders, bpy50 laptop
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_tiles_combined_NDVI_Rita_11062017" #EDGY Dean
in_dir1a <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_tile_1_NDVI_Rita_11062017" #EDGY Dean
in_dir1b <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017" #EDGY Dean
in_dir_ref <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI/rev_project_output"

### This will need to be changed: this should be a file with list of files to read from tile dir
data_fname1a <- file.path(in_dir1a,"dat_out_tile_1_NDVI_Rita_11062017.txt")
data_fname1b <- file.path(in_dir1b,"dat_out_tile_2_NDVI_Rita_11062017.txt")
data_fname_mae_zone_tb1a <- file.path(in_dir1a,"mae_zones_tb_tile_1_NDVI_Rita_11062017.txt")
data_fname_mae_zone_tb1b <- file.path(in_dir1b,"mae_zones_tb_tile_2_NDVI_Rita_11062017.txt")

data_fname_mae_tot_tb1a <- file.path(in_dir1a,"mae_zones_tb_tile_1_NDVI_Rita_11062017.txt")
data_fname_mae_tot_tb1b <- file.path(in_dir1b,"mae_zones_tb_tile_2_NDVI_Rita_11062017.txt")

###### constant
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1

################# START SCRIPT ###############################

#set up the working directory
#Create output directory

#out_dir <- in_dir #output will be created in the input dir
out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

############### START OF SCRIPT ##########

##Figure 1:    Concept for SBT (outside R) (Marco)
##Figure 2:    Study areas (outside R) (Stu)
##Figure 3:    Strata: Zonal areas maps: Generate here
##Figure 4:    Average Temporal profiles in a year by zones and overall to show impact of events on the variable
##Figure 5:    Spatial patterns: Maps of Observed, predicted, residuals
##Figure 8:    Temporal MAE patterns for event over x dates

# e.g. Figure 4a: temporal profiles for Dean NDVI

#################################################
## PART 1: Read the datasets ####

###  Load data if needed:
# if(class(r_temp_pred)=="character"){
#   #/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_2
#   r_temp_pred_df <- read.table(r_temp_pred,stringsAsFactors = F)
#   lf <-r_temp_pred_df[,1] 
#   r_temp_pred <- stack(lf)
# }
# if(class(r_spat_pred)=="character"){
#   r_spat_pred_df <- read.table(r_spat_pred,stringsAsFactors = F)
#   lf <-r_spat_pred_df[,1] 
#   r_spat_pred <- stack(lf)
# }
# if(class(s_raster)=="character"){
#   df_lf <- read.table(s_raster,sep=",",stringsAsFactors = F,header=T) 
#   s_raster <- stack(df_lf[,1])
# }

### Observed, predicted and residulas data for DEAN case study (1)

#r_temp_pred1 <- list.files(path=in_dir1,
#                           pattern="r_temp_pred_arima_arima_.*._tile_1_NDVI_Rita_11062017.tif",
#                           full.names=T)
#r_temp_pred2 <- list.files(path=in_dir2,
#                           pattern="r_temp_pred_arima_arima_.*._tile_2_NDVI_Rita_11062017.tif",
#                           full.names=T)

r_temp_pred <- stack(list.files(path=in_dir,
                           pattern="r_temp_pred_arima_arima_.*.tif",
                           full.names=T))
r_spat_pred <- stack(list.files(path=in_dir,
                          pattern="r_spat_pred_mle_eigen_no_previous_step_.*.tif",
                          full.names=T))

#r_spat_pred <- list.files(path=in_dir,
#                          pattern="r_spat_pred_mle_eigen_no_previous_step_.*.tif",
#                          full.names=T)

r_var <- stack(mixedsort(list.files(path=in_dir_ref,paste0(file_format,"$"),full.names=T))) #input raster images for the study area (276 images)

#### Prepare data: time window subset:

n_time_event <- as.numeric(unlist(strsplit(n_time_event,";")))

n_time_event_obs <- n_time_event[1]
n_time_event_spat <- n_time_event[2]
n_time_event_temp <- n_time_event[3]
n_time_subset <- c(n_time_event_obs -1 ,n_time_event_obs + 2)
steps_subset <- n_time_subset[1]:n_time_subset[2] 
r_obs <- subset(r_var,steps_subset)
#r_spat <- subset(r_spat_pred,5:8)
#r_temp <- subset(r_temp_pred,5:8)

n_time_subset <- c(n_time_event_spat - 1 ,n_time_event_spat + 2)
steps_subset <- n_time_subset[1]:n_time_subset[2] 

#### Will change here using input info
r_spat <- subset(r_spat_pred,5:8)
r_temp <- subset(r_temp_pred,5:8)

r_res_spat <- r_spat - r_obs
r_res_temp <- r_temp - r_obs
r_zonal <- subset(r_var,zonal_colnames)

############## PANEL MAPS
names_layers_obs <- c("Observed NDVI T-1","Observed NDVI T+1","Observed NDVI T+2","Observed NDVI T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted  T+3")
names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

##### First spatial pattern
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))
pix_res_height <- 500
pix_res_width <- 480

p <- levelplot(r_obs, margin=FALSE,
                ylab=NULL,xlab=NULL,
                scales=list(draw=FALSE), #don't draw geographic coordinates
                par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                    par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                names.attr= names_layers_obs,
                #main=paste(names_layers[i],"NDVI",sep=" "),
                col.regions=palette_colors,at=seq(-3000,10000,by=1000))

png_filename <- paste("Figure","_2_observed_",out_suffix,".png", sep="")
png(png_filename,
    height=pix_res_height*layout_m[2],width=pix_res_width*layout_m[1])
print(p) #to plot in a loop!!  
dev.off()

##spatial plot

layout_m <- c(4,1) #
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))
pix_res_height <- 500
pix_res_width <- 480

p_spat <- levelplot(r_spat, margin=FALSE,
                    ylab=NULL,xlab=NULL,
                    scales=list(draw=FALSE), #don't draw geographic coordinates
                    par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                        par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                    par.strip.text=list(font=2,cex=2),
                    layout= layout_m,
                    zlim=c(-3000,8000),
                    names.attr= names_layers_pred_spat,
                    #main=paste(names_layers[i],"NDVI",sep=" "),
                    col.regions=palette_colors,at=seq(-3000,10000,by=1000))

png_filename <- paste("Figure","_spatial_prediction_",out_suffix,".png", sep="")

png(png_filename,
    height=pix_res_height*layout_m[2],width=pix_res_width*layout_m[1])
print(p_spat) #to plot in a loop!!  

dev.off()

##temp plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_temp <- levelplot(r_temp, margin=FALSE,
                    ylab=NULL,xlab=NULL,
                    scales=list(draw=FALSE), #don't draw geographic coordinates
                    par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                        par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                    par.strip.text=list(font=2,cex=2),
                    layout= layout_m,
                    names.attr= names_layers_pred_temp,
                    #main=paste(names_layers[i],"NDVI",sep=" "),
                    col.regions=palette_colors,
                    at=seq(-3000,10000,by=1000))

png_filename <- paste("Figure","_time_prediction_",out_suffix,".png", sep="")

png(png_filename,
    height=480*layout_m[2],width=480*layout_m[1])

print(p_temp) #to plot in a loop!!  

dev.off()

##combined plots obs,spat,temp
layout_m <- c(4,3)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

r_all_var <- stack(r_obs,r_spat,r_temp)


p_all_var <- levelplot(r_all_var, margin=FALSE,
                       ylab=NULL,xlab=NULL,
                       scales=list(draw=FALSE), #don't draw geographic coordinates
                       par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                           par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                       par.strip.text=list(font=2,cex=2),
                       layout= layout_m,
                       names.attr= names_layers_all,
                       col.regions=palette_colors,
                       at=seq(-3000,10000,by=1000))

png_filename <- paste("Figure","_2_","combined_all_obs_time_space_pred_",out_suffix,".png", sep="")
png(png_filename,
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_var) #to plot in a loop!!  
dev.off()

#### Now plot residuals for NDVI

layout_m <- c(4,2)
no_brks <- 255
#palette_colors <- (matlab.like(no_brks))
palette_colors <- matlab.like(no_brks)

r_all_res <- stack(r_res_spat,r_res_temp)
names_layers_res_spat <- c("Spatial residuals T-1","Spatial residuals T+1","Spatial residuals T+2","Spatial residuals T+3")
names_layers_res_temp <- c("Temporal residuals T-1","Temporal residuals T+1","Temporal residuals T+2","Temporal residuals T+3")
names_layers_all_res <- c(names_layers_res_spat,names_layers_res_temp)

p_all_res <- levelplot(r_all_res, margin=FALSE,
                       ylab=NULL,xlab=NULL,
                       scales=list(draw=FALSE), #don't draw geographic coordinates
                       par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                           par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                       par.strip.text=list(font=2,cex=2),
                       layout= layout_m,
                       names.attr=names_layers_all_res,
                       col.regions=matlab.like(255),
                       at=seq(-3000,5000,by=200))
#col.regions=palette_colors)

png(paste("Figure","_2_","combined_all_residuals_models_time_and_space_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_res) #to plot in a loop!!  

dev.off()


###################

#######################################################
########## Figure 3:  Strata: Zonal areas maps

############
##Figure 3a: Zonal figure

freq_zonal_df <- as.data.frame(freq(r_zonal,merge=T))

n_zones <- sum(as.numeric(!is.na(freq_zonal_df$value)))
col_pal_all <- c("red","blue","green","brown","violet") #used in all the areas, use palette from earlier code

if(is.null(cat_name)){
  cat_name <- rev(paste("Zone",1:n_zones))
}

title_str <- "FEMA flood zones"

### Find number of zones from freq tb
layout_m <- c(1,1)
png(paste("Figure","_3a_","zonal_variable_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

#cat_name <- rev(c("Zone 3","Zone 4","Zone 5"))
#cat_name <- rev(c("Zone 2","Zone 1"))
#par(xpd = FALSE)
col_pal <- col_pal_all[1:n_zones]
plot(r_zonal,col=col_pal,legend=F)
#plot(r_zonal2,col=c("red","green"))EDG

#par(xpd = TRUE)
#"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".
legend("topleft", legend = cat_name, fill = rev(col_pal), 
       title="Zones",
       cex = 0.9, 
       #inset = 0.9,
       bty="n")
title(title_str)
dev.off()

############################################################
###Figure 4:  Average Temporal profiles overall for the time series under study
## This illustrate the change (dip) directly after the Hurricane event

#debug(accuracy_space_time_calc)
data_fname <- r_var
#r_ref <- NULL

#debug(accuracy_space_time_calc)
accuracy_space_and_time_obj <- accuracy_space_time_calc(r_temp_pred=r_temp_pred,
                                 r_spat_pred=r_spat_pred,
                                 s_raster=data_fname,
                                 proj_str=proj_str,
                                 time_window_selected=time_window_selected,
                                 n_time_event=n_time_event_obs,
                                 r_zonal=zonal_colnames,
                                 method_space=method_space,
                                 method_time=method_time,
                                 r_ref=r_ref,
                                 out_suffix=out_suffix,
                                 var_names=var_names,
                                 NA_flag_val=NA_flag_val,
                                 file_format=file_format, 
                                 date_range=date_range,
                                 out_dir=out_dir,
                                 create_out_dir_param=create_out_dir_param)

mae_tot_tb <- accuracy_space_and_time_obj$mae_tot_tb
mae_zones_tb <- accuracy_space_and_time_obj$mae_zones_tb

dim(mae_zones_tb)
dim(mae_tot_tb)

method_time <- unlist(strsplit(method_time,";"))
method_space <- unlist(strsplit(method_space,";"))

name_method_time <- paste0("temp_",method_time[1],"_",method_time[2])
name_method_space <- paste0("spat_",method_space[1],"_",method_space[2])
mae_tot_tb <- as.data.frame(mae_tot_tb)
row.names(mae_tot_tb) <- NULL
names(mae_tot_tb)<- c(name_method_space,name_method_time)
mae_tot_tb$time <- 1:nrow(mae_tot_tb)

## Input params for the plot:
#https://stackoverflow.com/questions/5506046/how-do-i-put-more-space-between-the-axis-labels-and-axis-title-in-an-r-boxplot

y_range <- range(cbind(mae_tot_tb[[name_method_space]],mae_tot_tb[[name_method_time]]))
legend_val <- c("tempor model","spatial model")
  
res_pix <- 960
col_mfrow<- 1
row_mfrow<- 0.7
png(filename=paste("Figure_temporal_profiles_MAE_tot_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

temp_formula_str <- as.formula(paste0(paste0(name_method_time," ~ ","time")))
spat_formula_str <- as.formula(paste0(paste0(name_method_space," ~ ","time")))

#par(mgp=c(0,1,0))
## margin for side 2 is 7 lines in size
#op <- par(mar=c(5,7,4,2) +0.1) ## default is c(5,4,4,2) + 0.1
op <- par(mar=c(5,7,6,2) +0.1) ## default is c(5,4,4,2) + 0.1, 4 is 4+32 for top

plot(temp_formula_str,type="b",
     col="cyan",
     pch=16, #filled in circles...
     lwd=3,
     data=mae_tot_tb,
     ylim=y_range,
     ylab="Mean Absolute Error (MAE)",
     xaxt="n",
     xlab="Time",
     cex.lab=2,
     cex.axis=1.5,
     cex=1.5)
lines(spat_formula_str, 
      type="b",
      lwd=3,
      pch=16, #filled in circles
      col="magenta",
      data=mae_tot_tb)

x_labels <- c("T-5","T-4","T-3","T-2","T-1","T+1","T+2","T+3","T+4","T+5")
axis(1, at = 1:10, labels = x_labels , cex.axis = 1.5)
#axis(2, cex.axis = 2)

legend("topleft",
       legend=legend_val,
       col=c("cyan","magenta"),
       lty=1,
       lwd=3,
       cex=1.5,
       bty="n")
title("Overall MAE for spatial and temporal models",cex.main=2.3,font=2) #Note that the results are different than for ARIMA!!!
par(op)
dev.off()

############################################################
###Figure 5:  Average Temporal profiles by zones for the time series under study

#### BY ZONES ASSESSMENT: Ok now it is general so it should be part of the function...

#mae_zones_tb <- rbind(ac_spat_mle_obj$mae_zones_tb[1:3,],
#                      ac_temp_obj$mae_zones_tb[1:3,])
#mae_zones_tb <- rbind(ac_spat_mle_eigen_no_previous_obj$mae_zones_tb,
#                      ac_spat_mle_eigen_with_previous_obj$mae_zones_tb,
#                      ac_temp_obj$mae_zones_tb)

mae_zones_tb <- accuracy_space_and_time_obj$mae_zones_tb

n_zones <- length(unique(mae_zones_tb$zone))


#mae_zones_tb$method <- c(rep("spat_reg_no",n_zones),rep("temp_with_reg",n_zones),rep("temp_arima_reg",n_zones))

#n_time <- ncol(mae_zones_tb) -1
#pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")
#names(mae_zones_tb) <- c("zones","pred1","pred2","pred3","pred4","method")
#names(mae_zones_tb) <- pred_names

#write.table(mae_zones_tb,file=paste("mae_zones_tb","_",out_suffix,".txt",sep=""))

mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
#drop first few rows that contain no data but zones...
n_methods <- 2 #2 because we have 2 methods... (one for space and one for time)
n_start <-n_zones*2 +1 
#n_start <-n_zones*2 +1
dd <- dd[n_start:nrow(dd),]

dd$zones <- mydata$zones
dd$method <- mydata$method

dd$zones <- mydata$zone #use recycle rule

#Note that to get the title correct, one needs to add
res_pix<-960
col_mfrow<-1
row_mfrow<-1
png(filename=paste("Figure1_ref_layer_time_1_",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

p <- xyplot(data~which | as.factor(zones) ,group=method,data=dd,type="b",xlab="time",ylab="VAR",
            strip = strip.custom(factor.levels=unique(as.character(dd$zones))), #fix this!!!
            auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
                            border = FALSE, lines = TRUE,cex=1.2))
try(print(p))

dev.off()
#histogram(r_zonal)

#####################end of added on 02/10/2018



#View(dat_out1a)
#zones_tb_avg<- zonal(r_var1,r_zonal1,fun='mean')

#zones_avg_df <- as.data.frame(zones_tb_avg)
#n_zones <- length(unique(zones_avg_df$zone))

#n_time <- ncol(zones_avg_df) - 1
#pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")

############
##Figure 4a: Dean case for a year?

#mean_vals <- colMeans(data_tb1[,6:281],na.rm=T)
#col_pal <- c("black",col_pal_all)

#layout_m <- c(1.5,1.1)

#png(paste("Figure","_4a_","average_temporal_profiles_by_zones_subset","NDVI_Dean_",out_suffix,".png", sep=""),
#    height=520*layout_m[2],width=520*layout_m[1])

#plot(1:23,mean_vals[139:161],type="b",pch=1,col=col_pal[1],ylim=c(3000,10000),ylab="NDVI",xlab="Time step (16 days)") #overall
# lines(1:23,zones_avg_df[1,139:161],type="b",pch=2,col=col_pal[2]) #zone 3
# lines(1:23,zones_avg_df[2,139:161],type="b",pch=3,col=col_pal[3]) #zone 4
# lines(1:23,zones_avg_df[3,139:161],type="b",pch=4,col=col_pal[4]) #zone 5
# abline(v=15.5,lty="dashed")
# legend("topleft",legend=c("Overall","zone 3","zone 4", "zone 5"),
#         cex=1.2, col=col_pal,bty="n",
#         lty=1,pch=1:4)
# legend("topright",legend=c("hurricane event"),cex=1.2,lty="dashed",bty="n")
# title("Overall and zonal averages NDVI in the Dean study area",cex=1.8, font=2)
# 
# dev.off()

############
##Figure 4b: Light use NLU for Katrina the whole time series

zones_tb_avg2<- zonal(r_var2,r_zonal2,fun='mean')

zones_avg_df2 <- as.data.frame(zones_tb_avg2)
n_zones <- length(unique(zones_avg_df2$zone))

n_time <- ncol(zones_avg_df2) -1
#pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")

##Subset to a year...for clarity

#zones_avg_df

mean_vals <- colMeans(data_tb2[,3:24],na.rm=T)
#pixval <- data_tb[800,var_names]
#pix300 <- data_tb[300,var_names]
col_pal <- c(col_pal_all[1:2])

layout_m <- c(1.5,1)

png(paste("Figure","_4b_","average_temporal_profiles_by_zones_subset","_light_data_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

plot(1:22,mean_vals[1:22],type="b",pch=1,col=col_pal[1],ylim=c(50,65),
     ylab="Light Intensity",xlab="Time step (annual)") #overall
lines(1:22,zones_avg_df2[1,2:23],type="b",pch=2,col=col_pal[2]) #zone 4
lines(1:22,zones_avg_df2[2,2:23],type="b",pch=3,col=col_pal[3]) #zone 5
#lines(1:22,zones_avg_df[3,2:23],type="b",pch=4,col="blue") #zone 6
abline(v=13.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 1","zone 2"),
       cex=0.8, col=col_pal,bty="n",
       lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=0.8,lty="dashed",bty="n")
title("Average light in the Katrina study area and by zones",cex=1.6, font=2)

dev.off()

############
##Figure 4c: Katrina case  temporal profiles

zones_tb_avg3<- zonal(r_var3,r_zonal3,fun='mean')
zones_avg_df3 <- as.data.frame(zones_tb_avg3)
n_zones <- length(unique(zones_avg_df3$zone))

mean_vals <- colMeans(data_tb3[,6:281],na.rm=T)
col_pal <- c("black",col_pal_all)

layout_m <- c(1.5,1.1)

png(paste("Figure","_4c_","average_temporal_profiles_by_zones_subset_","NDVI_Katrina_",out_suffix,".png", sep=""),
    height=520*layout_m[2],width=520*layout_m[1])

plot(1:23,mean_vals[139:161],type="b",pch=1,col=col_pal[1],ylim=c(3000,10000),ylab="NDVI",xlab="Time step (16 days)") #overall
lines(1:23,zones_avg_df3[1,139:161],type="b",pch=2,col=col_pal[2]) #zone 3
lines(1:23,zones_avg_df3[2,139:161],type="b",pch=3,col=col_pal[3]) #zone 4
lines(1:23,zones_avg_df3[3,139:161],type="b",pch=4,col=col_pal[4]) #zone 5
abline(v=15.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 3","zone 4", "zone 5"),
       cex=1.2, col=col_pal,bty="n",
       lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=1.2,lty="dashed",bty="n")
title("Overall and zonal averages NDVI in the Katrina Study area")
dev.off()


#########################################################################
##### Figure 8:  Temporal MAE patterns Dean four dates
##### Accuracy assessment by MAE for Dean data

#mae_zones_tb1 <- read.table(file.path(in_dir1,"mae_zones_tb_EDGY_predictions_03182015.txt"))
#mae_tot_tb1 <- read.table(file.path(in_dir1,"mae_tot_tb_EDGY_predictions_03182015.txt"))

##Make this a function for any variable
#layout_m <- c(1,1.2)

#Arguments
plot_filename <- paste("Figure_8_accuracy_","mae","_by_zone_and_four_timesteps","_","NDVI_Dean_",out_suffix,".png", sep="")
var_name <- "NDVI"
event_timestep <- 1.5
pix_res <- 480
input_data_df <- mae_zones_tb1
layout_m <- c(1,3)

#source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.

#debug(plot_by_zones_and_timestep_fun)
plot_zones_obj1 <- plot_by_zones_and_timestep_fun(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m)

## Now do the same for total...
#Arguments
plot_filename <- paste("Figure_8_accuracy_","mae","_by_total_and_four_timesteps","_","NDVI_Dean_",out_suffix,".png", sep="")
var_name <- "NDVI"
event_timestep <- 1.5
pix_res <- 520
input_data_df <- mae_tot_tb1
layout_m <- c(1,1)

#debug(plot_by_tot_and_timestep_fun)
plot_tot_obj1 <- plot_by_tot_and_timestep_fun(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m)

#################################################################
#### Figure 10:  Temporal MAE patterns Katrina NDVI four dates

#Arguments
plot_filename <- paste("Figure_10_accuracy_","mae","_by_zone_and_four_timesteps","_","NDVI_Katrina_",out_suffix,".png", sep="")
var_name <- "NDVI"
event_timestep <- 1.5
pix_res <- 480
#input_data_df <- mae_zones_tb3
layout_m <- c(1,3)

##Format the data first to zoom in the four relevant time steps +zones+method columns
input_data_df <- mae_zones_tb3[,c(1,8:11,18)]
names(input_data_df)[1] <- "zones"
#replace temp_reg by temp!!
input_data_df$method <- as.character(input_data_df$method)
input_data_df$method[input_data_df$method=="temp_reg"]<- "temp"
#source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.

#debug(plot_by_zones_and_timestep_fun)
plot_zones_obj3 <- plot_by_zones_and_timestep_fun(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m)


## Now do the same for total...
#Arguments
plot_filename <- paste("Figure_10_accuracy_","mae","_by_total_and_four_timesteps","_","NDVI_Katrina",out_suffix,".png", sep="")
var_name <- "NDVI"
event_timestep <- 1.5
pix_res <- 520
input_data_df <- mae_tot_tb3[c(7:10),]
layout_m <- c(1,1)

#debug(plot_by_tot_and_timestep_fun)
plot_tot_obj3 <- plot_by_tot_and_timestep_fun(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m)

###################################################
#FIGURE 1: temporal profiles for a specific time period (year here for NDVI)

### Figure for Katrina


#########################
######### NDVI Katrina

start_date <- "2005-01-01"
end_date<- "2005-12-31"
dates <- dates3
date_event <- date_event3
n_time_event <- n_time_event3
data_tb <- data_tb3
r_var <- r_var3
r_zonal <- r_zonal3
var_name <- "NDVI"
y_range <- c(2200,10000)
x_label<- "Dates 16-day time step"
title_str <- "Average NDVI for year 2005 in the Katrina study area and by zones"
out_suffix_str <- paste("NDVI_Katrina_",out_suffix,sep="")
out_dir

#start_date,end_date,dates,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str
#debug(plot_temporal_time_series_profile_by_zones)

plot_temporal_time_series_profile_by_zones(start_date,end_date,dates,date_event,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str)


################################### END OF SCRIPT #######################################


