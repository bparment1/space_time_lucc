####################################    Space Time Analyses PAPER   #######################################
###########################################  Figures production  #######################################
#This script produces figures based on outputs from earlier analyses for the space beats time project.
#The script uses spatial regression and temporal model (ARIMA and temporal OLS) values predicted for various. 
#Current case studies include:

#-raster NDVI MODIS data in the Yucatan after hurricane Rita

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/20/2015 
#DATE MODIFIED: 12/17/2017
#Version: 1
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones            
#PROJECT: Workshop for William and Mary: an intro to geoprocessing with R 
#PROJECT: Chicago AAG 2015 with Marco Millones
#PROJECT: Dallas Geocomputation 2015 with Marco Millones, Springer book
#
#COMMENTS: - 
#         - 
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

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_11242015_functions.R" #PARAM 1
function_paper_figures_analyses <- "space_beats_time_sbt_paper_figures_functions_01292016.R" #PARAM 1

script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts" #path on bpy50 #PARAM 2
#script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path on Atlas
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.
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


in_dir <- "/home/bparmentier/Google Drive/Space_beats_time" #bpy50 laptop
#in_dir <- "/home/parmentier/Data/Space_beats_time"

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"sbt_book_figures_12172017" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#### Clean this up: Need to this general for any region!!!
#Latest relevant folders, bpy50 laptop
in_dir <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_tiles_combined_NDVI_Rita_11062017" #EDGY Dean
in_dir1a <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_tile_1_NDVI_Rita_11062017" #EDGY Dean
in_dir1b <- "/home/bparmentier/Google Drive/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017" #EDGY Dean
in_dir_ref <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Rita_NDVI/rev_project_output"

data_fname1a <- file.path(in_dir1a,"dat_out_tile_1_NDVI_Rita_11062017.txt")
data_fname1b <- file.path(in_dir1b,"dat_out_tile_2_NDVI_Rita_11062017.txt")

data_tb1a <- read.table(data_fname1a,sep=",",header=T) #EDGY DEAN
data_tb1b <- read.table(data_fname1b,sep=",",header=T) #EDGY DEAN

mae_zones_tb1a <- read.table(file.path(in_dir1a,"mae_zones_tb_tile_1_NDVI_Rita_11062017.txt"))
mae_zones_tb1b <- read.table(file.path(in_dir1b,"mae_zones_tb_tile_2_NDVI_Rita_11062017.txt"))

mae_tot_tb1a <- read.table(file.path(in_dir1a,"mae_zones_tb_tile_1_NDVI_Rita_11062017.txt"))
mae_tot_tb1b <- read.table(file.path(in_dir1b,"mae_zones_tb_tile_2_NDVI_Rita_11062017.txt"))

## Parameters that vary from case studies to case studies...

#coord_names <- c("Long","Lat") #PARAM 11
#coord_names <- c("x","y") #PARAM 11
#coord_names <- c("XCoord","YCoord")
#coord_names <- c("POINT_X1","POINT_Y1")

#zonal_colnames <- "r_srtm_Katrina_rec2" #PARAM 12

#var_names <- 1:230 #PARAM 13 #Data is stored in the columns 3 to 22
#num_cores <- 11 #PARAM 14
#n_time_event <- 108 #PARAM 15 #this is the timestep corresponding to the event ie Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orelans
#time_window_selected <- var_names #PARAM 16: use alll dates for now
#time_window_selected <- 100:116 #PARAM 16: use alll dates for now

date_range <- c("2001.01.01","2010.12.31") #date

dates <- generate_dates_by_step(date_range[1],date_range[2],16)$dates
#Closest date to the even for each example:
n_time_event1 <- 154 #PARAM 15 # #timestep for Hurricane Katrina (Aug 23- Aub 31 2005): 235-243 DOY, storm surge Aug 29 in New Orleans

date_event1 <- "2005-09-24" #Dean

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

##Figure 1:    Concept for SBT (outside R) (Marco)
##Figure 2:    Study areas (outside R) (Stu)
##Figure 3:    Strata: Zonal areas maps (Stu)
##Figure 4:    Average Temporal profiles in a year by zones and overall to show impact of events on the variable
##Figure 5:    Spatial patterns Dean: Maps of Observed, predicted, residuals
##Figure 6:    Spatial patterns Katrina NLU: Maps of Observed, predicted, residuals
##Figure 7:    Spatial patterns Katrina NDVI: Maps of Observed, predicted, residuals
##Figure 8:    Temporal MAE patterns Dean four dates
##Figure 9:    Temporal MAE patterns Katrina NLU for four dates
##Figure 10:   Temporal MAE patterns Katrina NDVI four dates

#Figure index
#a: Dean NDVI
#b: Katrina NLU
#c: Katrina NDVI
#
# e.g. Figure 4a: temporal profiles for Dean NDVI

#################################################
## PART 1: Read the datasets ####

r_var1 <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst$",full.names=T)))
#r_var2 <- stack(mixedsort(list.files(path=in_dir2,"r_F.*.light_Katrina_03222015.rst$",full.names=T)))
#r_var3 <- stack(mixedsort(list.files(path=in_dir3,"r_reg2_NDVI.*.rst$",full.names=T)))

r_zonal1 <- raster(list.files(path=in_dir1,pattern="r_z_winds_EDGY_.*.rst$",full.names=T))
#r_zonal2 <- raster(file.path(in_dir2,"r_high1_lowminus1_light_Katrina_03222015.rst"))
#r_zonal3 <- raster(file.path(in_dir3,"r_r_srtm_Katrina_rec2_NDVI_Katrina_04182015.rst"))
#r_zonal3 <- raster(list.files(path=in_dir3,pattern="r_r_srtm_.*._rec2_NDVI_Katrina_04182015.rst$",full.names=T))

### Observed, predicted and residulas data for DEAN case study (1)

r_var <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst$",full.names=T))) #input raster images for the study area (276 images)
spat_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_pred_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst$",full.names=T)))
temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_arima_arima_.*._EDGY_predictions_03182015.rst$",full.names=T)))
spat_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst$",full.names=T)))
temp_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst$",full.names=T)))

#######################################################
########## Figure 3:  Strata: Zonal areas maps

col_pal_all <- c("red","blue","green") #used in all the areas

############
##Figure 3a: Dean case zonal stat

layout_m <- c(1,1)
png(paste("Figure","_3a_","zonal_variable_","Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

cat_name <- rev(c("Zone 3","Zone 4","Zone 5"))
par(xpd = FALSE)
col_pal <- col_pal_all
plot(r_zonal1,col=col_pal,legend=F)
#plot(r_zonal2,col=c("red","green"))EDG

par(xpd = TRUE)
legend(x = -87.4, y = 19.2, legend = cat_name, fill = rev(col_pal), 
       title="Zones",
       cex = 0.9, inset = 0.9,bty="n")

dev.off()

############
##Figure 3b: Katrina case zonal map , NLU

layout_m <- c(1,1)
png(paste("Figure","_3c_","zonal_variable_","light_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

col_pal <- col_pal_all[2]
cat_name <- rev(c("high","low"))
par(xpd = FALSE)
plot(r_zonal2,col=c("red","brown"),legend=F)
#plot(r_zonal2,col=c("red","green"))

par(xpd = TRUE)
legend(x = -89.88, y = 30, legend = cat_name, fill = rev(col_pal), 
       title="Zones",
       cex = 0.9, inset = 0.9,bty="n")

dev.off()


############
##Figure 3c: Katrina case zonal map 

r_zonal3 <- raster(file.path(in_dir3,"r_r_srtm_Katrina_rec2_NDVI_Katrina_04182015.rst"))
#r_r_srtm_Katrina_NDVI_Katrina_04182015.rst

col_pal <- col_pal_all
layout_m <- c(1,1)
png(paste("Figure","_3c_","zonal_variable_","NDVI_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

cat_name <- rev(c("high","medium","low"))
par(xpd = FALSE)
plot(r_zonal3,col=col_pal,legend=F)
#plot(r_zonal2,col=c("red","green"))

par(xpd = TRUE)
legend(x = -89.72, y = 30.10, legend = cat_name, fill = rev(col_pal), 
       title="Zones",
       cex = 0.9, inset = 0.9,bty="n")

dev.off()

############################# This is repeated?? 
############################################################
###Figure 4:  Average Temporal profiles overall for the time series under study
## This illustrate the change (dip) directly after the Hurricane event

zones_tb_avg<- zonal(r_var1,r_zonal1,fun='mean')

zones_avg_df <- as.data.frame(zones_tb_avg)
n_zones <- length(unique(zones_avg_df$zone))

n_time <- ncol(zones_avg_df) -1
#pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")


############
##Figure 4a: Dean case for a year?

mean_vals <- colMeans(data_tb1[,6:281],na.rm=T)
col_pal <- c("black",col_pal_all)

layout_m <- c(1.5,1.1)

png(paste("Figure","_4a_","average_temporal_profiles_by_zones_subset","NDVI_Dean_",out_suffix,".png", sep=""),
    height=520*layout_m[2],width=520*layout_m[1])

plot(1:23,mean_vals[139:161],type="b",pch=1,col=col_pal[1],ylim=c(3000,10000),ylab="NDVI",xlab="Time step (16 days)") #overall
lines(1:23,zones_avg_df[1,139:161],type="b",pch=2,col=col_pal[2]) #zone 3
lines(1:23,zones_avg_df[2,139:161],type="b",pch=3,col=col_pal[3]) #zone 4
lines(1:23,zones_avg_df[3,139:161],type="b",pch=4,col=col_pal[4]) #zone 5
abline(v=15.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 3","zone 4", "zone 5"),
        cex=1.2, col=col_pal,bty="n",
        lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=1.2,lty="dashed",bty="n")
title("Overall and zonal averages NDVI in the Dean study area",cex=1.8, font=2)

dev.off()

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


##########################################################################
##### Figure 5: Comparisons of observed, predicted temp, predicted spat: EDGY case study



#####################################################################
##### Figure 7: Comparisons of observed, predicted temp, predicted spat: Katrina NDVI case study

#r_var3 <- stack(mixedsort(list.files(path=in_dir3,"r_reg2_NDVI.*.rst",full.names=T)))

#spat_pred_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_spat_pred_mle_eigen_t_.*.NDVI_Katrina_04182015.rst",full.names=T)))
#temp_pred_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_temp_pred_arima_arima_.*.NDVI_Katrina_04182015.rst",full.names=T)))
#temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_lm_ols_t_.*._EDGY_predictions_03182015.rst",full.names=T)))

#r_spat_res_mle_eigen_t_15_NDVI_Katrina_04182015

#spat_res_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_spat_res_mle_eigen_t_.*.NDVI_Katrina_04182015.rst",full.names=T)))
#temp_res_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_temp_res_arima_arima_.*.NDVI_Katrina_04182015.rst",full.names=T)))

n_time_event <- 108
n_before <- n_time_event - 1
n_after <- n_time_event + 2

r_obs3 <- subset(r_var3,n_before:n_after) 
spat_pred_rast3 <- subset(spat_pred_rast3,8:11) 
temp_pred_rast3 <- subset(temp_pred_rast3,7:10) 
spat_res_rast3 <- subset(spat_res_rast3,8:11) 
temp_res_rast3 <- subset(temp_res_rast3,7:10) 

#names(r_var1) <- c("T\-1","T\+1","T\+2","T\+3")
#levelplot(r_var1,layers=n_before:n_after,col.regions=rev(terrain.colors(255))

no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

#names_layers <- c("Monthly Climatology","Daily deviation","Daily prediction")
names_layers <- c("T-1","T+1","T+2","T+3")
names_layers_obs <- c("Observed NDVI T-1","Observed NDVI T+1","Observed NDVI T+2","Observed NDVI T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted  T+3")
#names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

fig_nb <- c("2_t107","2_t108","2_t109","2_t110")
list_p <- vector("list",length=length(names_layers))
i<-1 # for testing

for(i in 1:nlayers(r_obs3)){
  
  p <- levelplot(r_obs3,layers=i, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                     par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 main=paste(names_layers_obs[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
  png(paste("Figure",fig_nb[i],"_observed_NDVI_Katrina_paper_space_beats_time_",out_suffix,".png", sep=""),
      height=480*1.4,width=480*1.4)
  print(p) #to plot in a loop!!  
  dev.off()
  
}

##### Combined figures NDVI Katrina

#levelplot(meot_rast_m,main=title_plot, ylab=NULL,xlab=NULL,,par.settings = list(axis.text = list(font = 2, cex = 1.5),
#                      par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),
#                      par.strip.text=list(font=2,cex=1.5),
#                      col.regions=temp.colors,at=seq(-1,1,by=0.02))

layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p1 <- levelplot(r_obs3, margin=FALSE,
                ylab=NULL,xlab=NULL,
                par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                    par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                names.attr= names_layers_obs,
                #main=paste(names_layers[i],"NDVI",sep=" "),
                col.regions=palette_colors,at=seq(-3000,10000,by=0.02))

png(paste("Figure","_2_","combined_observed_NDVI_Katina_1107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p1) #to plot in a loop!!  
dev.off()

##spatial plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_spat <- levelplot(spat_pred_rast3, margin=FALSE,
                    ylab=NULL,xlab=NULL,
                    par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                        par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                    par.strip.text=list(font=2,cex=2),
                    layout= layout_m,
                    zlim=c(-3000,8000),
                    names.attr= names_layers_pred_spat,
                    #main=paste(names_layers[i],"NDVI",sep=" "),
                    col.regions=palette_colors,at=seq(-3000,10000,by=0.02))

png(paste("Figure","_2_","combined_spat_NDVI_Katrina_107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_spat) #to plot in a loop!!  

dev.off()

##temp plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_temp <- levelplot(temp_pred_rast3, margin=FALSE,
                    ylab=NULL,xlab=NULL,
                    par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                        par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                    par.strip.text=list(font=2,cex=2),
                    layout= layout_m,
                    names.attr= names_layers_pred_temp,
                    #main=paste(names_layers[i],"NDVI",sep=" "),
                    col.regions=palette_colors,
                    at=seq(-3000,10000,by=0.02))

png(paste("Figure","_2_","combined_temp_NDVI_Katrina_107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_temp) #to plot in a loop!!  

dev.off()

##combined plots obs,spat,temp
layout_m <- c(4,3)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

r_all_var3 <- stack(r_obs3,spat_pred_rast3,temp_pred_rast3)

names_layers_obs <- c("Observed NDVI T-1","Observed NDVI T+1","Observed NDVI T+2","Observed NDVI T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted  T+3")
names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

p_all_var <- levelplot(r_all_var3, margin=FALSE,
                       ylab=NULL,xlab=NULL,
                       par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                           par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                       par.strip.text=list(font=2,cex=2),
                       layout= layout_m,
                       names.attr= names_layers_all,
                       col.regions=palette_colors,
                       at=seq(-3000,10000,by=0.02))

png(paste("Figure","_2_","combined_all_NDVI_Katrina_107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_var) #to plot in a loop!!  

dev.off()


#### Now plot residuals for NDVI

#spat_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
#temp_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))

layout_m <- c(4,2)
no_brks <- 255
#palette_colors <- (matlab.like(no_brks))
palette_colors <- matlab.like(no_brks)

r_all_res3 <- stack(spat_res_rast3,temp_res_rast3)
names_layers_res_spat <- c("Spatial residuals T-1","Spatial residuals T+1","Spatial residuals T+2","Spatial residuals T+3")
names_layers_res_temp <- c("Temporal residuals T-1","Temporal residuals T+1","Temporal residuals T+2","Temporal residuals T+3")
names_layers_all_res <- c(names_layers_res_spat,names_layers_res_temp)

p_all_res <- levelplot(r_all_res3, margin=FALSE,
                       ylab=NULL,xlab=NULL,
                       par.settings = list(axis.text = list(font = 2, cex = 1.5),
                                           par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                       par.strip.text=list(font=2,cex=2),
                       layout= layout_m,
                       names.attr=names_layers_all_res,
                       col.regions=matlab.like(255))
#col.regions=palette_colors)

png(paste("Figure","_2_","combined_all_res_107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_res) #to plot in a loop!!  

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

##########################################################################
##### Figure 9:  Temporal MAE patterns Katrina light four dates
##### accuracy assessment by MAE for Katrina light data

#Arguments
plot_filename <- paste("Figure_9_accuracy_","mae","_by_zone_and_four_timesteps","_","NLU_Katrina_",out_suffix,".png", sep="")
var_name <- "NLU"
event_timestep <- 1.5
pix_res <- 480
#input_data_df <- mae_zones_tb3
layout_m <- c(1,2)

##Format the data first to zoom in the four relevant time steps +zones+method columns
input_data_df <- mae_zones_tb2[,c(1,13:16,21)]
names(input_data_df)[1] <- "zones"
#replace temp_reg by temp!!
input_data_df$method <- as.character(input_data_df$method)
input_data_df$method[input_data_df$method=="temp_reg"]<- "temp"
#source(file.path(script_path,function_paper_figures_analyses)) #source all functions used in this script 1.

#debug(plot_by_zones_and_timestep_fun)

plot_zones_obj2 <- plot_by_zones_and_timestep_fun(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m)

## Now do the same for total...
#Arguments
plot_filename <- paste("Figure_9_accuracy_","mae","_by_total_and_four_timesteps","_","NLU_Katrina",out_suffix,".png", sep="")
var_name <- "NLU"
event_timestep <- 1.5
pix_res <- 520
input_data_df <- mae_tot_tb2[12:15,]
layout_m <- c(1,1)

#debug(plot_by_tot_and_timestep_fun)
plot_tot_obj2 <- plot_by_tot_and_timestep_fun(plot_filename,var_name,event_timestep,pix_res,input_data_df,layout_m)


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

####################################
### Figure for EDGY

#> dates[153]
#[1] "2007-08-13"

start_date <- "2007-01-01"
end_date<- "2007-12-31"
dates <- dates1
n_time_event <- n_time_event1 #the effect of the event on the data are not observed on the 
                              #actual date of the event!!! So let's add another variable for
                              #the exact event date itself.
date_event1 <- "2007-08-21" #August 21, Hurrifcane Dean landfall in Yucatan
date_event <- date_event1
data_tb <- data_tb1
r_var <- r_var1
r_zonal <- r_zonal1
var_name <- "NDVI"
y_range <- c(2200,10000)
x_label<- "Dates 16-day time step"
title_str <- "Average NDVI for year 2007 in the Dean study area and by zones"
out_suffix_str <- paste("NDVI_Dean_",out_suffix,sep="")
out_dir

#start_date,end_date,dates,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str
#debug(plot_temporal_time_series_profile_by_zones)

plot_temporal_time_series_profile_by_zones(start_date,end_date,dates,date_event,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str)

######### NDVI NLU
#date_range2
start_date <- "1992-01-01"
end_date<- "2013-12-31"
dates <- dates2
date_event <- date_event2
n_time_event <- n_time_event2
data_tb <- data_tb2
r_var <- r_var2
r_zonal <- r_zonal2
var_name <- "NLU"
y_range <- c(50,65)
x_label<- "Dates annual time step"
title_str <- "Average NLU for year 1992-2013 in the Katrina study area and by zones"
out_suffix_str <- paste("NLU_Katrina_",out_suffix,sep="")
out_dir

#start_date,end_date,dates,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str
#debug(plot_temporal_time_series_profile_by_zones)

plot_temporal_time_series_profile_by_zones(start_date,end_date,dates,date_event,n_time_event,data_tb,r_var,r_zonal,var_name,y_range,x_label,title_str,out_dir,out_suffix_str)

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


################### END OF SCRIPT ##################


