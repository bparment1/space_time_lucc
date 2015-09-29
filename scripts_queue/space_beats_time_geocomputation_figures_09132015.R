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
#DATE MODIFIED: 09/19/2015
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

###### Functions used in this script sourced from other files

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_05172015_functions.R" #PARAM 1
script_path <- "/home/parmentier/Data/Space_beats_time/sbt_scripts" #path to script #PARAM 2
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

##### Functions used in this script 

generate_dates_modis <-function(start_date,end_date,step_date){
  library(xts)
  library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(strftime(st,"%Y"),strftime(en,"%Y")) #extract year
  
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
  }
  
  dates_modis <-as.Date(unlist((ll_list))) 
  return(dates_modis)
}

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

#debug(generate_dates_modis)
#test_dates <-generate_dates_modis(start_date="2001-01-01",end_date="2010-12-31",step_date="16 day")
#start_date<-"2001.01.01"
#end_date <- "2012.12.31"
#end_date <- "2002.12.31"

#st <- as.Date(start_date,format="%Y.%m.%d")
#en <- as.Date(end_date,format="%Y.%m.%d")
#ll <- seq.Date(st, en, by="1 day")
#ll <- seq.Date(st, en, by="16 day")
#dates_queried <- format(ll,"%Y.%m.%d")

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
out_suffix <-"geocomputation_figures_09132015" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

#Latest relevant folders...
in_dir1 <- "/home/parmentier/Data/Space_beats_time/output_EDGY_predictions_03182015" #EDGY Dean
#in_dir2 <- "/home/parmentier/Data/Space_beats_time/output__predictions_09252014" #pop Katrina
in_dir2 <- "/home/parmentier/Data/Space_beats_time/output_light_Katrina_03222015" #light Katrina
in_dir3 <- "/home/parmentier/Data/Space_beats_time/output_NDVI_Katrina_04182015" #NDVI Katrina

data_fname1 <- file.path(in_dir1,"dat_out_EDGY_predictions_03182015.txt")
#data_fname2 <- file.path(in_dir2,"dat_out__predictions_09252014")
data_fname2 <- file.path(in_dir2,"dat_out_light_Katrina_03222015.txt")
data_fname3 <- file.path(in_dir3,"dat_out_NDVI_Katrina_04182015.txt")

data_tb1 <- read.table(data_fname1,sep=",",header=T) #EDGY DEAN
#data_tb2 <- read.table(data_fname,sep=",",header=T) #pop Katrina landscan
data_tb2 <- read.table(data_fname2,sep=",",header=T) #light Katrina
data_tb3 <- read.table(data_fname3,sep=",",header=T) #NDVI Katrina

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

#########################################
########### Case 1: NDVI EDGY data

###############################
##### Figure 1: zonal contrast based on average of observed values: EDGY area

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

layout_m <- c(1.8,1)

p_zones_avg <- xyplot(data~time ,group=zones,
               data=dd,type="b",xlab="Time step",ylab="NDVI",
               auto.key=list(columns=1,space="right",title="Zones",cex=1,pch=16),
               main="Average NDVI in the Dean study area and by zones",
               border = FALSE, lines = TRUE,cex=0.6,lwd=0.6,pch=16
)

png(paste("Figure","_1a_","average_temporal_profiles_by_zones_","NDVI_Dean_",out_suffix,".png", sep=""),
    height=550*layout_m[2],width=550*layout_m[1])
print(p_zones_avg)

update(p_zones_avg,panel=function(...){ 
        panel.xyplot(...) 
        panel.abline(v=154) 
} ) 

dev.off()

##Subset to a year...for clarity

#zones_avg_df

mean_vals <- colMeans(data_tb1[,6:281],na.rm=T)
#pixval <- data_tb[800,var_names]
#pix300 <- data_tb[300,var_names]
layout_m <- c(1.5,1.1)

png(paste("Figure","_1b_","average_temporal_profiles_by_zones_subset","NDVI_Dean_",out_suffix,".png", sep=""),
    height=520*layout_m[2],width=520*layout_m[1])

plot(1:23,mean_vals[139:161],type="b",pch=1,col="red",ylim=c(3000,10000),ylab="NDVI",xlab="Time step (16 days)") #overall
lines(1:23,zones_avg_df[1,139:161],type="b",pch=2,col="black") #zone 4
lines(1:23,zones_avg_df[2,139:161],type="b",pch=3,col="green") #zone 5
lines(1:23,zones_avg_df[3,139:161],type="b",pch=4,col="blue") #zone 6
abline(v=15.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 3","zone 4", "zone 5"),
        cex=1.2, col=c("red","black","green","blue"),bty="n",
        lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=1.2,lty="dashed",bty="n")
title("Overall and zonal averages NDVI in the Dean study area",cex=1.8, font=2)

dev.off()

###############################
##### Figure 2: Comparisons of observed, predicted temp, predicted spat: EDGY case study

r_var1 <- stack(mixedsort(list.files(path=in_dir1,"r_NDVI.*.rst",full.names=T)))

spat_pred_rast1 <- stack(mixedsort(list.files(path=in_dir1,"r_spat_pred_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
temp_pred_rast1 <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))
#temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_lm_ols_t_.*._EDGY_predictions_03182015.rst",full.names=T)))

spat_res_rast1 <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
temp_res_rast1 <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))

n_time_event <- 154
n_before <- n_time_event - 1
n_after <- n_time_event + 2

r_obs1 <- subset(r_var1,n_before:n_after) 
#names(r_var1) <- c("T\-1","T\+1","T\+2","T\+3")
#levelplot(r_var1,layers=n_before:n_after,col.regions=rev(terrain.colors(255))

no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

#names_layers <- c("Monthly Climatology","Daily deviation","Daily prediction")
names_layers <- c("T-1","T+1","T+2","T+3")
names_layers_obs <- c("Observed NDVI T-1","Observed NDVI T+1","Observed NDVI T+2","Observed NDVI T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted T+3")
names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

fig_nb <- c("2_t153","2_t154","2_t155","2_t156")
list_p <- vector("list",length=length(names_layers))
i<-1 # for testing

for(i in 1:nlayers(r_obs1)){
  
  p <- levelplot(r_obs1,layers=i, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 main=paste(names_layers_obs[i],sep=" "),
                 col.regions=palette_colors)
  
  png(paste("Figure",fig_nb[i],"_observed_NDVI_paper_space_beats_time_","NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*1.4,width=480*1.4)
  print(p) #to plot in a loop!!  
  dev.off()

}

##### Combined figures EDGY

#levelplot(meot_rast_m,main=title_plot, ylab=NULL,xlab=NULL,,par.settings = list(axis.text = list(font = 2, cex = 1.5),
#                      par.main.text=list(font=2,cex=2.2),strip.background=list(col="white")),
#                      par.strip.text=list(font=2,cex=1.5),
#                      col.regions=temp.colors,at=seq(-1,1,by=0.02))

layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p1 <- levelplot(r_obs1, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr = names_layers_obs,
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors,at=seq(-3000,10000,by=0.02))
  
png(paste("Figure","_2_","combined_observed_paper_space_beats_time_","NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p1) #to plot in a loop!!  
dev.off()

##spatial plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_spat <- levelplot(spat_res_rast1, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr = names_layers_pred_spat,
                 col.regions=palette_colors,at=seq(-3000,10000,by=0.02))
  
png(paste("Figure","_2_","combined_spat_paper_space_beats_time_NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_spat) #to plot in a loop!!  

dev.off()

##temp plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_temp <- levelplot(temp_pred_rast1, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr = names_layers_pred_temp,
                 col.regions=palette_colors,
                 at=seq(-3000,10000,by=0.02))
  
png(paste("Figure","_2_","combined_temp_paper_space_beats_time_NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_temp) #to plot in a loop!!  

dev.off()

##Full figure combined plots of obs data,spat predictions,temp predictions
layout_m <- c(4,3)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

r_all_var1 <- stack(r_obs1,spat_pred_rast1,temp_pred_rast1)

p_all_var <- levelplot(r_all_var1, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr = names_layers_all,
                 col.regions=palette_colors,
                 at=seq(-3000,10000,by=0.02))

png(paste("Figure","_2_","combined_all_paper_space_beats_time_NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_var) #to plot in a loop!!  

dev.off()

#### NOW PLOT RESIDUALS FROM DEAN NDVI PREDICTIONS ####

#spat_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
#temp_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))
names_layers_res_spat <- c("Spatial residuals T-1","Spatial residuals T+1","Spatial residuals T+2","Spatial residuals T+3")
names_layers_res_temp <- c("Temporal residuals T-1","Temporal residuals T+1","Temporal residuals T+2","Temporal residuals T+3")
names_layers_all_res <- c(names_layers_res_spat,names_layers_res_temp)


layout_m <- c(4,2)
no_brks <- 255
#palette_colors <- (matlab.like(no_brks))
palette_colors <- matlab.like(no_brks)

r_all_res1 <- stack(spat_res_rast1,temp_res_rast1)

p_all_res <- levelplot(r_all_res1, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                names.attr=names_layers_all_res,
                col.regions=matlab.like(255))
                #col.regions=palette_colors)

png(paste("Figure","_2_","combined_all_res_paper_space_beats_time_NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_res) #to plot in a loop!!  

dev.off()

############################### DEAN case study
##### Figure 3: accuracy assessment by MAE

mae_zones_tb <- read.table(file.path(in_dir1,"mae_zones_tb_EDGY_predictions_03182015.txt"))
mae_tot_tb <- read.table(file.path(in_dir1,"mae_tot_tb_EDGY_predictions_03182015.txt"))

### mae for total region
layout_m <- c(1,1.2)
png(paste("Figure_3a_accuracy_","mae","_by_tot_and_timestep","_","NDVI_Dean_",out_suffix,".png", sep=""),
    height=560*layout_m[1],width=560*layout_m[2])

y_range <- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
y_range <- c(y_range[1],y_range[2]+100)
#xlab_tick <- mae_tot_tb$time
xlab_tick <- c("T-1","T+1","T+2","T+3")
x_tick_position <- mae_tot_tb$time

plot(spat_reg ~ time, type="l",col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "NDVI",xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6)
points(spat_reg ~ time, col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "NDVI",xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6,cex.pch=1.6)
lines(temp ~ time, col="magenta",lwd=3,pch=16,data=mae_tot_tb)
points(temp ~ time, col="magenta",lwd=3,pch=16,data=mae_tot_tb,cex.pch=1.6)
axis(1,at= x_tick_position,labels=xlab_tick)

abline(v=1.5,lty="dashed")
#legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n")
legend("top",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.4)
#legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.4)

legend("topright",legend=c("spatial","temporal"),col=c("cyan","magenta"),lty=1,lwd=3,bty="n",cex=1.6)
title("Overall MAE for Spatial and Temporal models", cex.main=2) #Note that the results are different than for ARIMA!!!

dev.off()

### mae by zones
mydata <- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
#drop first four rows
dd <- dd[7:nrow(dd),]


#xyplot(data~which |zones,group=method,data=dd,type="b",xlab="time",ylab="VAR",
#       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
#      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
#                     border = FALSE, lines = TRUE,cex=1.2)
#)

#xyplot(data~which |zones,group=method,data=dd,type="b",xlab="time",ylab="VAR",
#+       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
#+      auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
#+                     border = FALSE, lines = TRUE,cex=1.2)

layout_m <- c(1,3)

p<- xyplot(data~which |zones,group=method,data=dd,type="b",xlab="Time Step",ylab="NDVI",
       strip = strip.custom(factor.levels=as.character(unique(dd$zones))),
       #strip = strip.custom(factor.levels=as.character(unique(dd$which))),
       auto.key=list(columns=1,space="right",title="Model",cex=1),
                     border = FALSE, lines = TRUE,cex=1.2)

png(paste("Figure_3b_accuracy_","mae","by_timestep_and_zones","_NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
print(p)
dev.off()

layout_m <- c(1,3)
png(paste("Figure_3b_accuracy_","mae","_by_zones_and_timestep","_NDVI_Dean_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p <- xyplot(data~zones |which,group=method,data=dd,type="b",xlab="zones",ylab="NDVI",
      strip = strip.custom(factor.levels=as.character(unique(dd$which))),
      #auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
      auto.key=list(columns=1,space="right",title="Model",cex=1),
      border = FALSE, lines = TRUE,cex=1.2
)

print(p)
dev.off()


###
layout_m <- c(1,3)

#layout_m <- c(1,1.2)
png(paste("Figure_3b_accuracy_","mae","_by_zone_and_timestep","_","NDVI_Dean_",out_suffix,".png", sep=""),
    height=560*layout_m[1],width=560*layout_m[2])

input_data <- mae_zones_tb
for ( i in 3:5){
  col_names <- unique(input_data$method)
  input_data <- subset(input_data,zones==3)
  input_data <- subset(input_data,select=c(names(input_data)!="zones"))
  input_data <- subset(input_data,select=c(names(input_data)!="method"))
  
  names(input_data) <- col_names
  
  y_range <- range(cbind(input_data$spat_reg,input_data$temp))
  y_range <- c(y_range[1],y_range[2]+100)
  #xlab_tick <- mae_tot_tb$time
  xlab_tick <- c("T-1","T+1","T+2","T+3")
  x_tick_position <- mae_tot_tb$time

  plot(spat_reg ~ time, type="l",col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "NDVI",xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6)
  points(spat_reg ~ time, col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "NDVI",xlab="Time step",xaxt="n",lwd=3,pch=16,cex.lab=1.6,cex.pch=1.6)
  lines(temp ~ time, col="magenta",lwd=3,pch=16,data=mae_tot_tb)
  points(temp ~ time, col="magenta",lwd=3,pch=16,data=mae_tot_tb,cex.pch=1.6)
  axis(1,at= x_tick_position,labels=xlab_tick)

  abline(v=1.5,lty="dashed")
  #legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n")
  legend("top",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.4)
  #legend("topleft",legend=c("hurricane event"),lty="dashed",bty="n",cex=1.4)

  legend("topright",legend=c("spatial","temporal"),col=c("cyan","magenta"),lty=1,lwd=3,bty="n",cex=1.6)
  title("Overall MAE for Spatial and Temporal models", cex.main=2) #Note that the results are different than for ARIMA!!!

}

dev.off()



######################################
########## Case 2: light data Katrina

## Elevation figures:

r_zonal2 <- raster(file.path(in_dir2,"r_high1_lowminus1_light_Katrina_03222015.rst"))

layout_m <- c(1,1)
png(paste("Figure","_0_","zonal_variable_","light_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

cat_name <- rev(c("high","low"))
par(xpd = FALSE)
plot(r_zonal2,col=c("black","red"),legend=F)
#plot(r_zonal2,col=c("red","green"))

par(xpd = TRUE)
legend(x = -89.88, y = 30, legend = cat_name, fill = rev(c("red","black")), 
       title="Zones",
    cex = 0.9, inset = 0.9,bty="n")

dev.off()

###############################
##### Figure 4: zonal contrast based on average of observed values: EDGY area

r_var2 <- stack(mixedsort(list.files(path=in_dir2,"r_F.*.light_Katrina_03222015.rst",full.names=T)))
#r_F101992_light_Katrina_03222015.rst


#r_high1_lowminus1_light_Katrina_03222015
#zonal stat proffile...

zones_tb_avg2<- zonal(r_var2,r_zonal2,stat='mean')

zones_avg_df <- as.data.frame(zones_tb_avg2)
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

layout_m <- c(1.5,1)

p_zones_avg <- xyplot(data~time ,group=zones,
               data=dd,type="b",xlab="Time step",ylab="Light Intensity",
               auto.key=list(columns=1,space="right",title="Zones",cex=1),
               main="Average light in the Katrina study area and by zones",
               border = FALSE, lines = TRUE,cex=1.2
)

png(paste("Figure","_1a_","average_temporal_profiles_by_zones_","light_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_zones_avg)

update(p_zones_avg,panel=function(...){ 
        panel.xyplot(...) 
        panel.abline(v=154) 
} ) 

dev.off()

##Subset to a year...for clarity

#zones_avg_df

mean_vals <- colMeans(data_tb2[,3:24],na.rm=T)
#pixval <- data_tb[800,var_names]
#pix300 <- data_tb[300,var_names]
layout_m <- c(1.5,1)

png(paste("Figure","_1b_","average_temporal_profiles_by_zones_subset","_light_data_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

plot(1:22,mean_vals[1:22],type="b",pch=1,col="red",ylim=c(50,65),
     ylab="Light Intensity",xlab="Time step (annual)") #overall
lines(1:22,zones_avg_df[1,2:23],type="b",pch=2,col="black") #zone 4
lines(1:22,zones_avg_df[2,2:23],type="b",pch=3,col="green") #zone 5
#lines(1:22,zones_avg_df[3,2:23],type="b",pch=4,col="blue") #zone 6
abline(v=13.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 1","zone 2"),
        cex=0.8, col=c("red","black","green"),bty="n",
        lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=0.8,lty="dashed",bty="n")
title("Average light in the Katrina study area and by zones",cex=1.6, font=2)

dev.off()


###############################
##### Figure 6: Comparisons of observed, predicted temp, predicted spat: EDGY case study

#r_var3 <- stack(mixedsort(list.files(path=in_dir3,"r_reg2_NDVI.*.rst",full.names=T)))
r_var2 <- stack(mixedsort(list.files(path=in_dir2,"r_F.*.light_Katrina_03222015.rst",full.names=T)))

spat_pred_rast2 <- stack(mixedsort(list.files(path=in_dir2,"r_spat_pred_mle_eigen_t_.*.light_Katrina_03222015.rst",full.names=T)))
temp_pred_rast2 <- stack(mixedsort(list.files(path=in_dir2,"r_temp_pred_lm_ols_t.*.light_Katrina_03222015.rst",full.names=T)))
#r_temp_res_lm_ols_t_20_light_Katrina_03222015
#r_spat_pred_ols_ols_t_1_light_Katrina_03222015.rst
#temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_lm_ols_t_.*._EDGY_predictions_03182015.rst",full.names=T)))

#r_spat_res_mle_eigen_t_15_NDVI_Katrina_04182015

spat_res_rast2 <- stack(mixedsort(list.files(path=in_dir2,"r_spat_res_mle_eigen_t_.*.light_Katrina_03222015.rst",full.names=T)))
temp_res_rast2 <- stack(mixedsort(list.files(path=in_dir2,"r_temp_res_lm_ols_t.*.light_Katrina_03222015.rst",full.names=T)))

n_time_event <- 14
n_before <- n_time_event - 1
n_after <- n_time_event + 2

r_obs2 <- subset(r_var2,n_before:n_after) 
spat_pred_rast2 <- subset(spat_pred_rast2,13:16) 
temp_pred_rast2 <- subset(temp_pred_rast2,12:15) 
spat_res_rast2 <- subset(spat_res_rast2,13:16) 
temp_res_rast2 <- subset(temp_res_rast2,12:15) 

no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

#names_layers <- c("Monthly Climatology","Daily deviation","Daily prediction")
names_layers <- c("T-1","T+1","T+2","T+3")
names_layers_obs <- c("Observed light T-1","Observed light T+1","Observed light T+2","Observed light T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted  T+3")
#names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

names_layers_obs <- c("Observed light T-1","Observed light T+1","Observed light T+2","Observed light T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted  T+3")
names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

fig_nb <- c("2_t13","2_t14","2_t15","2_t16")
list_p <- vector("list",length=length(names_layers))
i<-1 # for testing

for(i in 1:nlayers(r_obs2)){
  
  p <- levelplot(r_obs2,layers=i, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),par.strip.text=list(font=2,cex=2),
                 main=paste(names_layers_obs[i],"NDVI",sep=" "),
                 col.regions=palette_colors)
  
  png(paste("Figure",fig_nb[i],"_observed_light_Katrina_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*1.4,width=480*1.4)
  print(p) #to plot in a loop!!  
  dev.off()

}

##### Combined figures light Katrina

layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p1 <- levelplot(r_obs2, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr= names_layers_obs,
                 col.regions=palette_colors,at=seq(-3000,10000,by=0.02))
  
png(paste("Figure","_2_","combined_observed_light_Katina_1107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p1) #to plot in a loop!!  
dev.off()

##spatial plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_spat <- levelplot(spat_pred_rast2, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 zlim=c(-3000,8000),
                 names.attr= names_layers_pred_spat,
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors,at=seq(-3000,10000,by=0.02))
  
png(paste("Figure","_2_","combined_spat_light_Katrina_107_110_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_spat) #to plot in a loop!!  

dev.off()

##temp plot
layout_m <- c(4,1)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

p_temp <- levelplot(temp_pred_rast2, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr= names_layers_pred_temp,
                 #main=paste(names_layers[i],"NDVI",sep=" "),
                 col.regions=palette_colors,
                 ,at=seq(-3000,10000,by=0.02))
  
png(paste("Figure","_2_","combined_temp_light_Katrina_13_16_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_temp) #to plot in a loop!!  

dev.off()

##combined plots obs,spat,temp
layout_m <- c(4,3)
no_brks <- 255
palette_colors <- rev(terrain.colors(no_brks))

r_all_var2 <- stack(r_obs2,spat_pred_rast2,temp_pred_rast2)

names_layers_obs <- c("Observed light T-1","Observed light T+1","Observed light T+2","Observed light T+3")
names_layers_pred_spat <- c("Spatial predicted T-1","Spatial predicted T+1","Spatial predicted T+2","Spatial predicted T+3")
names_layers_pred_temp <- c("Temporal predicted T-1","Temporal predicted T+1","Temporal predicted T+2","Temporal predicted  T+3")
names_layers_all <- c(names_layers_obs,names_layers_pred_spat,names_layers_pred_temp)

p_all_var <- levelplot(r_all_var2, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                 par.strip.text=list(font=2,cex=2),
                 layout= layout_m,
                 names.attr= names_layers_all,
                 col.regions=palette_colors,
                 at=seq(-3000,10000,by=0.02))

png(paste("Figure","_2_","combined_all_light_Katrina_13_16_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_var) #to plot in a loop!!  

dev.off()

####

#spat_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_spat_res_mle_Chebyshev_t_.*.EDGY_predictions_03182015.rst",full.names=T)))
#temp_res_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_res_arima_arima_.*._EDGY_predictions_03182015.rst",full.names=T)))

layout_m <- c(4,2)
no_brks <- 255
#palette_colors <- (matlab.like(no_brks))
palette_colors <- matlab.like(no_brks)

r_all_res2 <- stack(spat_res_rast2,temp_res_rast2)
names_layers_res_spat <- c("Spatial residuals T-1","Spatial residuals T+1","Spatial residuals T+2","Spatial residuals T+3")
names_layers_res_temp <- c("Temporal residuals T-1","Temporal residuals T+1","Temporal residuals T+2","Temporal residuals T+3")
names_layers_all_res <- c(names_layers_res_spat,names_layers_res_temp)

p_all_res <- levelplot(r_all_res2, margin=FALSE,
                 ylab=NULL,xlab=NULL,
                 par.settings = list(axis.text = list(font = 2, cex = 1.5),
                              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
                par.strip.text=list(font=2,cex=2),
                layout= layout_m,
                names.attr=names_layers_all_res,
                col.regions=matlab.like(255))
                #col.regions=palette_colors)

png(paste("Figure","_2_","combined_all_res_light_Katrina_13_16_paper_space_beats_time_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])
print(p_all_res) #to plot in a loop!!  

dev.off()

###############################
##### Figure 6: accuracy assessment by MAE for Katrina light data

mae_zones_tb <- read.table(file.path(in_dir2,"mae_zones_tb_light_Katrina_03222015.txt"))
mae_tot_tb <- read.table(file.path(in_dir2,"mae_tot_tb_light_Katrina_03222015.txt"))

### mae for total region
layout_m <- c(,1.2)
png(paste("Figure_3a_accuracy_","mae","_by_tot_and_timestep","_","light_Katrina_",out_suffix,".png", sep=""),
    height=520*layout_m[1],width=520*layout_m[2])

#y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
#plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range,ylab="MAE",xlab="Time step (annual)")
#lines(temp ~ time, type="b",col="magenta",data=mae_tot_tb)

y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
xlab_tick <- mae_tot_tb$time
#xlab_tick <- c("T-1","T+1","T+2","T+3")
x_tick_position <- mae_tot_tb$time

plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "Light Intensity",xlab="Time step",xaxt="n",pch=16,cex.lab=1.2)
lines(temp ~ time, type="b",col="magenta",pch=16,data=mae_tot_tb)
axis(1,at= x_tick_position,labels=xlab_tick)
legend("topleft",legend=c("spat","temp"),col=c("cyan","magenta"),lty=1,bty="n")
title("Overall MAE for spatial and temporal models") #Note that the results are different than for ARIMA!!!

#print(p)
dev.off()

### mae by zones
#mydata<- mae_zones_tb
#dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
#dd$lag <- mydata$lag 
#dd$zones <- mydata$zones
#dd$method <- mydata$method
#drop first four rows
#dd <- dd[7:nrow(dd),]

n_zones <- length(unique(mae_zones_tb$zone))
n_time <- ncol(mae_zones_tb) -1
pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")

mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
#drop first few rows that contain no data but zones...
n_start <-n_zones*2 +1
dd <- dd[n_start:nrow(dd),]
dd$zones <- mydata$zone #use recycle rule

xyplot(data~which |zones,group=method,data=dd,type="b",xlab="Time step (annual)",ylab="MAE for Light intensity",
       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
       auto.key=list(columns=1,space="right",title="Model",cex=1),
       border = FALSE, lines = TRUE,cex=1.2)

#layout_m <- c(1.4,1.4)
layout_m <- c(1,2)


p<- xyplot(data~which |zones,group=method,data=dd,type="b",xlab="Time step (annual)",ylab="MAE for Light intensity",
       strip = strip.custom(factor.levels=as.character(unique(dd$zones))),
       #strip = strip.custom(factor.levels=as.character(unique(dd$which))),
       auto.key=list(columns=1,space="right",title="Model",cex=1),
                     border = FALSE, lines = TRUE,cex=1.2)

png(paste("Figure_3b_accuracy_","mae","by_timestep_and_zones","_light_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
print(p)
dev.off()

layout_m <- c(1.4,1.4)
png(paste("Figure_3b_accuracy_","mae","_by_zones_and_timestep","_light_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p <- xyplot(data~zones |which,group=method,data=dd,type="b",xlab="Zones",ylab="MAE for Light intensity",
      strip = strip.custom(factor.levels=as.character(unique(dd$which))),
      #auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
      auto.key=list(columns=1,space="right",title="Model",cex=1),
      border = FALSE, lines = TRUE,cex=1.2
)

print(p)
dev.off()

######################################
########## CASE 3: NDVI Katrina

r_zonal3 <- raster(file.path(in_dir3,"r_r_srtm_Katrina_rec2_NDVI_Katrina_04182015.rst"))
#r_r_srtm_Katrina_NDVI_Katrina_04182015.rst

layout_m <- c(1,1)
png(paste("Figure","_0_","zonal_variable_","NDVI_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

cat_name <- rev(c("high","medium","low"))
par(xpd = FALSE)
plot(r_zonal3,col=c("black","yellow","red"),legend=F)
#plot(r_zonal2,col=c("red","green"))

par(xpd = TRUE)
legend(x = -89.72, y = 30.10, legend = cat_name, fill = rev(c("red","yellow","black")), 
       title="Zones",
    cex = 0.9, inset = 0.9,bty="n")

dev.off()

###############################
##### Figure 7: zonal contrast based on average of observed values: EDGY area

r_var3 <- stack(mixedsort(list.files(path=in_dir3,"r_reg2_NDVI.*.rst",full.names=T)))

r_zonal3 <- raster(list.files(path=in_dir3,pattern="r_r_srtm_.*._rec2_NDVI_Katrina_04182015.rst$",full.names=T))
#r_high1_lowminus1_light_Katrina_03222015
#zonal stat proffile...

zones_tb_avg3<- zonal(r_var3,r_zonal3,stat='mean')

zones_avg_df <- as.data.frame(zones_tb_avg3)
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

layout_m <- c(1.8,1)

p_zones_avg <- xyplot(data~time ,group=zones,
               data=dd,type="b",xlab="Time step",ylab="NDVI",
               auto.key=list(columns=1,space="right",title="Zones",cex=1,pch=16),
               main="Average NDVI in the Katrina study area and by zones",
               border = FALSE, lines = TRUE,cex=0.6,lwd=0.6,pch=16
)


png(paste("Figure","_1a_","average_temporal_profiles_by_zones_","NDVI_Katrina_",out_suffix,".png", sep=""),
    height=550*layout_m[2],width=550*layout_m[1])
print(p_zones_avg)

update(p_zones_avg,panel=function(...){ 
        panel.xyplot(...) 
        panel.abline(v=108) 
} ) 

dev.off()

##Subset to a year...for clarity

#zones_avg_df

mean_vals <- colMeans(data_tb3[,1:230],na.rm=T)
#pixval <- data_tb[800,var_names]
#pix300 <- data_tb[300,var_names]
layout_m <- c(1.5,1)

png(paste("Figure","_1b_","average_temporal_profiles_by_zones_subset","_NDVI_data_",out_suffix,".png", sep=""),
    height=480*layout_m[2],width=480*layout_m[1])

plot(1:230,mean_vals[1:230],type="b",pch=1,col="red",ylim=c(2000,8000),
     ylab="Light Intensity",xlab="Time step (annual)") #overall
lines(1:230,zones_avg_df[1,2:231],type="b",pch=2,col="black") #zone 4
lines(1:230,zones_avg_df[2,2:231],type="b",pch=3,col="green") #zone 5
lines(1:230,zones_avg_df[3,2:231],type="b",pch=4,col="blue") #zone 6
abline(v=107.5,lty="dashed")
legend("topleft",legend=c("Overall","zone 1","zone 2","zone 3"),
        cex=0.8, col=c("red","black","green","blue"),bty="n",
        lty=1,pch=1:4)
legend("topright",legend=c("hurricane event"),cex=0.8,lty="dashed",bty="n")
title("Average NDVI in the Katrina study area and by zones",cex=1.6, font=2)

dev.off()

###############################
##### Figure 8: Comparisons of observed, predicted temp, predicted spat: EDGY case study

r_var3 <- stack(mixedsort(list.files(path=in_dir3,"r_reg2_NDVI.*.rst",full.names=T)))

spat_pred_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_spat_pred_mle_eigen_t_.*.NDVI_Katrina_04182015.rst",full.names=T)))
temp_pred_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_temp_pred_arima_arima_.*.NDVI_Katrina_04182015.rst",full.names=T)))
#temp_pred_rast <- stack(mixedsort(list.files(path=in_dir1,"r_temp_pred_lm_ols_t_.*._EDGY_predictions_03182015.rst",full.names=T)))

#r_spat_res_mle_eigen_t_15_NDVI_Katrina_04182015

spat_res_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_spat_res_mle_eigen_t_.*.NDVI_Katrina_04182015.rst",full.names=T)))
temp_res_rast3 <- stack(mixedsort(list.files(path=in_dir3,"r_temp_res_arima_arima_.*.NDVI_Katrina_04182015.rst",full.names=T)))

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
                 ,at=seq(-3000,10000,by=0.02))
  
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

############################
######## 

mae_zones_tb <- read.table(file.path(in_dir3,"mae_zones_tb_NDVI_Katrina_04182015.txt"))
#mae_tot_tb_NDVI_Katrina_04182015
mae_tot_tb <- read.table(file.path(in_dir3,"mae_tot_tb_NDVI_Katrina_04182015.txt"))

### mae for total region
layout_m <- c(1,1.2)
png(paste("Figure_3a_accuracy_","mae","_by_tot_and_timestep","_","NDVI_Katrina_",out_suffix,".png", sep=""),
    height=560*layout_m[1],width=560*layout_m[2])

#y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
#plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range,ylab="MAE",xlab="Time step (annual)")
#lines(temp ~ time, type="b",col="magenta",data=mae_tot_tb)

y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
xlab_tick <- mae_tot_tb$time
#xlab_tick <- c("T-1","T+1","T+2","T+3")
x_tick_position <- mae_tot_tb$time

plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "NDVI",xlab="Time step",xaxt="n",pch=16,cex.lab=1.2)
lines(temp ~ time, type="b",col="magenta",pch=16,data=mae_tot_tb)
axis(1,at= x_tick_position,labels=xlab_tick)

legend("topleft",legend=c("spat","temp"),col=c("cyan","magenta"),pch=16,lty=1,bty="n")
title("Overall MAE for Spatial and Temporal models") #Note that the results are different than for ARIMA!!!

#print(p)
dev.off()

### mae for total region
layout_m <- c(1,1.2)
png(paste("Figure_3a_accuracy_","mae","_by_tot_and_timestep","_","NDVI_Dean_",out_suffix,".png", sep=""),
    height=560*layout_m[1],width=560*layout_m[2])

y_range<- range(cbind(mae_tot_tb$spat_reg,mae_tot_tb$temp))
#xlab_tick <- mae_tot_tb$time
xlab_tick <- c("T-1","T+1","T+2","T+3")
x_tick_position <- mae_tot_tb$time

plot(spat_reg ~ time, type="b",col="cyan",data=mae_tot_tb,ylim=y_range,
     ylab= "NDVI",xlab="Time step",xaxt="n",pch=16,cex.lab=1.2)
lines(temp ~ time, type="b",col="magenta",pch=16,data=mae_tot_tb)
axis(1,at= x_tick_position,labels=xlab_tick)

legend("topleft",legend=c("spat","temp"),col=c("cyan","magenta"),pch=16,lty=1,bty="n")
title("Overall MAE for Spatial and Temporal models") #Note that the results are different than for ARIMA!!!
dev.off()

### mae by zones

n_zones <- length(unique(mae_zones_tb$zone))
n_time <- ncol(mae_zones_tb) -1
pred_names <- c("zone",paste("t",2:n_time,sep="_"),"method")

mydata<- mae_zones_tb
dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
dd$lag <- mydata$lag 
dd$zones <- mydata$zones
dd$method <- mydata$method
#drop first few rows that contain no data but zones...
n_start <-n_zones*2 +1
dd <- dd[n_start:nrow(dd),]
dd$zones <- mydata$zone #use recycle rule

xyplot(data~which |zones,group=method,data=dd,type="b",xlab="Time step (16 day)",ylab="MAE for NDVI intensity",
       #strip = strip.custom(factor.levels=c("z3","z4","z5")), #fix this!!!
       auto.key=list(columns=1,space="right",title="Model",cex=1),
       border = FALSE, lines = TRUE,cex=1.2)

layout_m <- c(1.4,1.4)
layout_m <- c(1,3)

p<- xyplot(data~which |zones,group=method,data=dd,type="b",xlab="Time step (16 day)",ylab="MAE for NDVI intensity",
       strip = strip.custom(factor.levels=as.character(unique(dd$zones))),
       #strip = strip.custom(factor.levels=as.character(unique(dd$which))),
       auto.key=list(columns=1,space="right",title="Model",cex=1),
                     border = FALSE, lines = TRUE,cex=1.2)

png(paste("Figure_3b_accuracy_","mae","by_timestep_and_zones","_NDVI_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
print(p)
dev.off()

layout_m <- c(1.4,1.4)
#layout_m <- c(1,3)

png(paste("Figure_3b_accuracy_","mae","_by_zones_and_timestep","_NDVI_Katrina_",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])

p <- xyplot(data~zones |which,group=method,data=dd,type="b",xlab="Zones",ylab="MAE for NDVI",
      strip = strip.custom(factor.levels=as.character(unique(dd$which))),
      #auto.key = list("topright", corner = c(0,1),# col=c("black","red"),
      auto.key=list(columns=1,space="right",title="Model",cex=1),
      border = FALSE, lines = TRUE,cex=1.2
)

print(p)
dev.off()

#######

################### END OF SCRIPT ##################