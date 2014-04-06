####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: Preparation of Data  #######################################
#This script prepares data for the Space-Time Symposium at Williams and Mary.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 04/06/2014 
#DATE MODIFIED: 04/06/2014
#Version: 1
#PROJECT: GLP Conference Berlin and Space time project, YUCATAN CASE STUDY with Marco Millones             
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
library(forecast)
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)

###### Functions used in this script

#Insert here
create_output_dir <- function(i,l_out_dir,out_suffix=NULL){
  out_dir <-l_out_dir[[i]]
  
  if (!file.exists(out_dir)){
  dir.create(out_dir)
  #} else{
  #  out_path <-paste(out_path..)
  }
}

##### Parameters and arguments

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"/home/parmentier/Data/Space_Time" #On Atlas server...
#in_dir <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/"
out_dir <- in_dir
#set up the working directory
setwd(in_dir)

moore_window <- file.path(in_dir,"moore_window.rst")
#test_shp_path <- file.path(in_dir,"GEODA_Analysis") #path to get to the files for spatial pred
#  "/Users/benoitparmentier/Dropbox/Data/Space_Time/GEODA_Analysis/TEST_SPATIAL_ONE.shp"
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
proj_str<- CRS_WGS84
CRS_interp <- proj_modis_str

function_analyses_paper <- "MODIS_and_raster_processing_functions_04062014.R"
script_path <- file.path(in_dir,"R") #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
infile_reg_outline <- file.path(in_dir,"GYRS_MX_trisate_sin_windowed.shp")  #input region outline defined by polygon: Oregon

infile_modis_grid<-file.path(in_dir,"modis_sinusoidal_grid_world.shp") #modis grid tiling system, global
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_interp <- proj_modis_str

file_format <- ".tif"
ref_sejidos_name <- file.path(in_dir,"00_sejidos_group_sel5_ids.rst")
ref_winds_name <- file.path(in_dir,"00_windzones_moore_sin.rst")

out_suffix <-"04072014" #output suffix for the files that are masked for quality and for 
    
## Other specific parameters
NA_flag_val<- -9999

############ PART I: READ IN MODIS PROCESSED DATA ###############

#Read raster stack
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
r_stack <- stack(reg_var_list)
projection(r_stack) <- proj_modis_str

#levelplot(r_stack,layers=1:12) #show first 12 images (half a year more or less)

### We need to reproject the raster images...

#first create folder to store reprojected data:

l_out_dir <- list(file.path(in_dir,"moore_NDVI_wgs84"))

create_output_dir(1,l_out_dir)
moore_window <- file.path(in_dir,"moore_window.rst")
projection(moore_w) <- proj_modis_str
ref_rast <- projectRaster(from=moore_w,crs=CRS_WGS84,method="ngb") #Check that it is using ngb

##Create output names for region
out_suffix_var <-paste(out_suffix,file_format,sep="") 
list_var_names <-extension(reg_var_list,"")
var_list_outnames <- change_names_file_list(list_var_names,out_suffix_var,"moore_",file_format,out_path=l_out_dir[[1]])     

j<-1
#list_param_create_region<-list(j,raster_name=list_var_mosaiced,reg_ref_rast=ref_rast,out_rast_name=var_list_outnames)
list_param_create_region<-list(j,reg_var_list,ref_rast,var_list_outnames,NA_flag_val,proj_modis_str)
names(list_param_create_region) <-c("j","raster_name","reg_ref_rast","out_rast_name","NA_flag_val","input_proj_str")

#debug(create__m_raster_region)
#test<-create__m_raster_region(1,list_param_create_region)

#reg_var_moore_list <-mclapply(1:11, list_param=list_param_create_region, create__m_raster_region,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
reg_var_list_moore_list <-mclapply(1:length(var_list_outnames), list_param=list_param_create_region, create__m_raster_region,mc.preschedule=FALSE,mc.cores = 11) #This is the end bracket from mclapply(...) statement
#reg_var_list <-lapply(1:length(var_list_outnames), list_param=list_param_create_region, create__m_raster_region) 

#test<-stack(reg_var_moore_list[1:4])
#plot(test)

############ PART 2: READ IN ADDITIONAL VARIABLES FOR SPACE-TIME ANALYSES ###############

#ref_sejidos_name <- "~/Data/Space_Time/00_sejidos_group_sel5_ids.rst"
#ref_winds_name <- "~/Data/Space_Time/00_windzones_moore_sin.rst"


r_clip <- ref_rast
NAvalue(r_clip) <- 0

r_winds <- raster(ref_winds_name)
r_sejidos <- raster(ref_sejidos_name)
projection(r_winds) <- proj_modis_str
projection(r_sejidos) <- proj_modis_str

r_winds_m <- projectRaster(r_winds,ref_rast,method="ngb") #Check that it is using ngb
r_sejidos_m <- projectRaster(r_sejidos,ref_rast,method="ngb") #Check that it is using ngb
#add state

r_x <- ref_rast #create raster image that will contain x coordinates
r_y <- ref_rast #create raster image that will contain y coordiates
values(r_x) <- coordinates(ref_rast)[,1] #set values in raster image to x coord
values(r_y) <- coordinates(ref_rast)[,2] #set values in raster image to y coord
pix_id_r <- ref_rast
values(pix_id_r) <- 1:ncell(ref_rast)

s_dat_var <-stack(pix_id_r,r_x,r_y,r_winds_m,r_sejidos_m)
names(s_dat_var) <- c("pix_id_r","r_x","r_y","z_winds","sejidos")

projection(s_dat_var) <- CRS_WGS84
plot(s_dat_var)


## Now add the NDVI data for only the EDGY/moore region of Yucatan

reg_var_m <- stack(reg_var_moore_list)

s_dat_var <- stack(s_dat_var,reg_var_m)

r_clip <- ref_rast
NAvalue(r_clip) <- 0

s_dat_var <- mask(s_dat_var,r_clip)

EDGY_dat_spdf <- as(s_dat_var,"SpatialPointsDataFrame") #create a SpatialPointsDataFrame
NDVI_names <- paste("NDVI_",1:276,sep="")
names(EDGY_dat_spdf) <- c(c("pix_id_r","r_x","r_y","z_winds","sejidos"),NDVI_names)
#change names!!
#Save spdf as .RData object
save(EDGY_dat_spdf,file= file.path(out_dir,paste("EDGY_dat_spdf_",out_suffix,".RData",sep="")))

#Save spdf as shapefile, note that column names can be truncated
outfile1<-file.path(out_dir,paste("EDGY_dat_spdf","_",out_suffix,".shp",sep=""))
writeOGR(EDGY_dat_spdf,dsn= dirname(outfile1),layer= sub(".shp","",basename(outfile1)), driver="ESRI Shapefile",overwrite_layer=TRUE)

#Save spdf as delimited csv text file.
outfile1<-file.path(out_dir,paste("EDGY_dat_spdf","_",out_suffix,".txt",sep=""))
write.table(as.data.frame(EDGY_dat_spdf),file=outfile1,sep=",")

############### END OF SCRIPT ##################