############################    SBT comparison of results with pysal   #######################################
################################  Generation of predictions with python  #######################################
#This script produces spatial prediction using spatial regression.
#AUTHORS: Benoit Parmentier                                           
#DATE CREATED:  01/09/2015
#DATE MODIFIED: 01/12/2015
#Version: 1
#PROJECT: SBT
#################################################################################################

###Loading R library and packages                                                      

library(sp)  #Spatial objects definition
library(spdep) #Spatial objects functions for analyses
library(rasterVis) #Raster visualization
library(raster) #Raster objects definition and function for analyses
library(rgdal) #GDAL binding for R

library(rgeos) #GEOS binding for R
library(gtools) #general additional tools
library(maptools) #mapping tools
library(colorRamps) #Palette/coloramp for display,contains matlab.like color palette

###### Functions used in this script


#function_TOC_creation <- "TOC.R"
#function_plot_TOC <- "TOCplot.R"

#script_path <- "C:/Users/parmentier/Dropbox/Data/TOC/My Run" #path to script
#source(file.path(script_path,function_TOC_creation)) #source all functions used in this script 1.
#Using source, reads the script containing the functions, loads functions in the workspace/enviroment
#making them available for the user.

#####  Parameters and argument set up ###########

in_dir <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/pysal_test/"
in_dir1 <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/output__predictions_09252014"

data_test <-read.table(file.path(in_dir,"pysal_res_2004_09252014.csv"),sep=",")
names(data_test) <- c("ID","res_2004","obs_2004","rnd_val","pred_2004")
dim(data_test)

mean(data_test$rnd_val)
sd(data_test$rnd_val)

lf_csv <-mixedsort(list.files(path=in_dir,pattern="pysal_res_.*._09252014.csv")) #,full.names=T)

lf_shp <-mixedsort(list.files(path=in_dir1,pattern="^r_poly_t.*.shp")) #,full.names=T)
#
data_reg <- readOGR(in_dir1,gsub(".shp","",lf_shp[i]))
head(data_reg)
head(data_test)

data_tmp<-cbind(data_reg,data_test)
head(data_tmp)

coordinates(data_tmp) <- coordinates(data_reg)
mean(data_tmp$obs_2004)
mae_val <- mean(abs(data_tmp$res_2004))
lf_rst <- list.files(in_dir1,pattern=".rst",full.names=T)
r <- raster(lf_rst[1])

r_res <- rasterize(data_tmp,r,"res_2004")
r_pred <- rasterize(data_tmp,r,"pred_2004")
r_obs <- rasterize(data_tmp,r,"obs_2004")


#link  to shapefile...
#calculate MAE per date
#plot MAE and compare to other methods
#do difference compared to mle by dates
#check random number

#compare maps...