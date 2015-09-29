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
#DATE MODIFIED: 04/23/2015
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
out_suffix <-"AAG_figures_04232015" #output suffix for the files and ouptu folder #PARAM 8
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
in_dir4 <- "/home/parmentier/Data/Space_beats_time/output_NDVI_Katrina_04182015" #NDVI Katrina

data_fname1 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname2 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname3 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")
data_fname4 <- file.path("/home/parmentier/Data/Space_beats_time/R_Workshop_April2014","Katrina_Output_CSV - Katrina_pop.csv")

data_tb1 <- read.table(data_fname,sep=",",header=T)
data_tb2 <- read.table(data_fname,sep=",",header=T)
data_tb3 <- read.table(data_fname,sep=",",header=T)
data_tb4 <- read.table(data_fname,sep=",",header=T)


################# START SCRIPT ###############################

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
##### Figure 1: zonal contrast based on average of observed values

r_var1 <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst",full.names=T)))
r_zonal <- raster(list.files(path=in_dir1,pattern="r_z_winds_EDGY_.*.rst$",full.names=T))
#zonal stat proffile...

zones_tb_avg<- zonal(r_var1,r_zonal,stat='mean')

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

p_zones_avg <- xyplot(data~time ,group=zones,
               data=dd,type="b",xlab="time",ylab="VAR",
               auto.key=list(columns=1,space="right",title="Zones",cex=1),
               border = FALSE, lines = TRUE,cex=1.2
)

png(paste("Figure","_1a_","average_temporal_profiles_by_zones_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_zones_avg)

update(p_zones_avg,panel=function(...){ 
        panel.xyplot(...) 
        panel.abline(v=154) 
} ) 

dev.off()

##Subset to a year...for clarity

png(paste("Figure","_1b_","average_temporal_profiles_by_zones_subset",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_zones_avg)

update(p_zones_avg,panel=function(...){ 
        panel.xyplot(...) 
        panel.abline(v=154) 
}) 

dev.off()

###############################
##### Figure 2: Comparisons of observed, predicted temp, predicted spat

r_var1 <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst",full.names=T)))

spat_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_pred_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))
#temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_lm_ols_t_.*._EDGY_predictions_03182015.rst",full.names=T)))

spat_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
temp_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))

n_time_event <- 154
n_before <- n_time_event - 1
n_after <- n_time_event + 2

r_obs <- subset(r_var1,n_before:n_after) 
#names(r_var1) <- c("T\-1","T\+1","T\+2","T\+3")
#levelplot(r_var1,layers=n_before:n_after,col.regions=rev(terrain.colors(255))

no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

#names_layers <- c("Monthly Climatology","Daily deviation","Daily prediction")
names_layers <- c("T-1","T+1","T+2","T+3")

fig_nb <- c("2_t153","2_t154","2_t155","2_t156")
list_p <- vector("list",length=length(names_layers))
i<-1 # for testing

for(i in 1:nlayers(r_var)){
  
  p <- levelplot(r_var,layers=i, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
  png(paste("Figure",fig_nb[i],"_observed_NDVI_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*1.4,width=480*1.4)
  print(p) #to plot in a loop!!  
  dev.off()

}

##### Combined figures EDGY

layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p1 <- levelplot(r_obs, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
png(paste("Figure","_2_","combined_observed_NDVI_152_155_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p1) #to plot in a loop!!  
dev.off()

##spatial plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_spat <- levelplot(spat_pred_rast, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 zlim=c(-3000,8000),
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
png(paste("Figure","_2_","combined_spat_NDVI_152_155_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_spat) #to plot in a loop!!  

dev.off()

##temp plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_temp <- levelplot(temp_pred_rast, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                zlim=c(-3000,8000),
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
png(paste("Figure","_2_","combined_temp_NDVI_152_155_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_temp) #to plot in a loop!!  

dev.off()

##combined plots obs,spat,temp
layout_m <- c(4,3)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

r_all_var <- stack(r_obs,spat_pred_rast,temp_pred_rast)

p_all_var <- levelplot(r_all_var, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                zlim=c(-3000,8000),
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)

png(paste("Figure","_2_","combined_all_NDVI_153_156_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_var) #to plot in a loop!!  

dev.off()


####

#spat_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
#temp_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))

layout_m <- c(4,2)
no_brks <- 255
palette_colors <- (matlab.like(no_brks)

r_all_res <- stack(spat_res_rast,temp_res_rast)

p_all_res <- levelplot(r_all_res, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                zlim=c(-3000,8000),
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors)

png(paste("Figure","_2_","combined_all_res_152_155_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_res) #to plot in a loop!!  

dev.off()

###############################
##### Figure 3: accuracy assessment by MAE

mae_zones_tb <- read.table(file.path(in_dir1,"mae_zones_tb_EDGY_predictions_03182015.txt"))
mae_tot_tb <- read.table(file.path(in_dir1,"mae_tot_tb_EDGY_predictions_03182015.txt"))

### mae for total region
layout_m <- c(1.4,1.4)
png(paste("Figure_3a_accuracy_","mae","_EDGY_","by_tot_and_timestep","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range)
lines(temp ~ time, type="b",col="magenta",data=mae_tot_tb)
legend("topleft",legend=c("spat","temp"),col=c("cyan","magenta"),lty=1)
title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!

print(p)
dev.off()

### mae by zones
mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
#drop first four rows
dd <- dd[7:nrow(dd),]

layout_m <- c(1.4,1.4)

p<- xyplot(data~which |zones,group=method,data=dd,type="b",xlab="year",ylab="NDVI",
       strip = strip.custom(factor.levels=aPs.character(unique(dd$zones))),
       #strip = strip.custom(factor.levels=as.character(unique(dd$which))),
       auto.key=list(columns=1,space="right",title="Model",cex=1),
                     border = FALSE, lines = TRUE,cex=1.2)

png(paste("Figure_3b_accuracy_","mae","_EDGY_","by_timestep_and_zones","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
print(p)
dev.off()

layout_m <- c(1.4,1.4)
png(paste("Figure_3b_accuracy_","mae","_EDGY_","by_zones_and_timestep","_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p <- xyplot(data~zones |which,group=method,data=dd,type="b",xlab="zones",ylab="NDVI",
      strip = strip.custom(factor.levels=as.character(unique(dd$which))),
      #auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
      auto.key=list(columns=1,space="right",title="Model",cex=1),
      border = FALSE, lines = TRUE,cex=1.2
)

print(p)
dev.off()


################### END OF SCRIPT ##################