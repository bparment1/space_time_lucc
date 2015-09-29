####################################    Space Time Analyses PAPER   #######################################
###########################################  Figures production  #######################################
#This script produces figures based on outputs from earlier analyses for the space beats time project.
#The script uses spatial regression and temporal model (ARIMA and temporal OLS) values predicted for various. 
#Current case studies include:

#-raster NDVI MODIS data in the Yucatan after hurricane Dean
#-raster population data in New Orleans after hurricane Katrina
#-raster DSMP ligth data in New Orleans after hurricane Katrina
#-raster NDVI MODIS data in New Orleans after hurricane Katrina

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/20/2015 
#DATE MODIFIED: 04/20/2015
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#PROJECT: Chicago AAG 2015 with Marco Millones
#PROJECT: Dallas Geocomputation 2015 with Marco Millones
#
#COMMENTS: - 
#         - 
#TO DO:
#
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

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_04182015_functions.R" #PARAM 1
script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir <- "~/Data/Space_beats_time/case3data/" #lights/table" #PARAM3
in_dir <- "/home/parmentier/Data/Space_beats_time/"
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
out_suffix <-"geocomputation_figures_04202015" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#data_fname <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
#data_fname <- file.path(in_dir,"lights/table","Kat_lights.txt") #PARAM 10
#data_fname <- file.path(in_dir,"output_Katrina_04082015","dat_reg_var_list_NDVI_Katrina_04082015.txt")
#data_fname <- file.path(in_dir,"dat_reg2_var_list_NDVI_NDVI_Katrina_04102015.txt")

#coord_names <- c("Long","Lat") #PARAM 11
coord_names <- c("x","y") #PARAM 11

#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

#zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12

#var_names <- 1:230 #PARAM 13 #Data is stored in the columns 3 to 22
#num_cores <- 11 #PARAM 14
#n_time_event <- 108 #PARAM 15 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
#time_window_selected <- var_names #PARAM 16: use alll dates for now
#time_window_selected <- 100:116 #PARAM 16: use alll dates for now

#Latest relevant folders...
in_dir1 <- "/home/parmentier/Data/Space_beats_time/output_EDGY_predictions_03182015" #EDGY Dean
in_dir2 <- "/home/parmentier/Data/Space_beats_time/output__predictions_09252014" #pop Katrina
in_dir3 <- "/home/parmentier/Data/Space_beats_time/output_light_Katrina_03222015" #light Katrina
in_dir4 <- "/home/parmentier/Data/Space_beats_time/output_geocomputation_04202015" #NDVI Katrina

data_fname1 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname2 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname3 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname4 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")

data_tb1 <- read.table(data_fname,sep=",",header=T)
data_tb2 <- read.table(data_fname,sep=",",header=T)
data_tb3 <- read.table(data_fname,sep=",",header=T)
data_tb4 <- read.table(data_fname,sep=",",header=T)


################# START SCRIPT ###############################

##specific processing done for srm
#r_dem <- raster("/data/project/layers/commons/data_workflow/inputs/dem-cgiar-srtm-1km-tif/srtm_1km.tif")
#r_dem_Katrina<-crop(r_dem,r_stack)
#writeRaster(r_dem_Katrina,file.path(in_dir,"r_srtm_Katrina.rst"))
#library(ggmap)
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

############### START OF SCRIPT ##########

###############################
##### Figure 1

r_var1 <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst",full.names=T)))
n_time_event <- 154
n_before <- n_time_event - 1
n_after <- n_time_event + 2

r_var <- subset(r_var1,n_before:n_after) 
#names(r_var1) <- c("T\-1","T\+1","T\+2","T\+3")
#levelplot(r_var1,layers=n_before:n_after,col.regions=rev(terrain.colors(255))

no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

#names_layers <- c("Monthly Climatology","Daily deviation","Daily prediction")
names_layers <- c("T-1","T+1","T+2","T+3")

fig_nb <- c("1_t153","1_t154","1_t155","1_t156")
list_p <- vector("list",length=length(names_layers))
i<-1 # for testing

for(i in 1:nlayers(r_var)){
  
  p <- levelplot(r_var,layers=i, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
  png(paste("Figure",fig_nb[i],"_paper_climatology_daily_deviation_daily_prediction_model7_levelplot_",out_suffix,".png", sep=""),
    height=480*1.4,width=480*1.4)
  print(p) #to plot in a loop!!  
  dev.off()

}

p1 <- levelplot(r_var, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
png(paste("Figure","_1_","combined_152_155_paper_climatology_daily_deviation_daily_prediction_model7_levelplot_",out_suffix,".png", sep=""),
    height=480*2.8,width=480*2.8)
print(p1) #to plot in a loop!!  
dev.off()

###############################
##### Figure 2

r_spat_res_mle_Chebyshev_t_154EDGY_predictions_03182015.rst

r_temp_pred_lm_ols_t_153_EDGY_predictions_03182015.rst

r_var1 <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst",full.names=T)))


################### END OF SCRIPT ##################