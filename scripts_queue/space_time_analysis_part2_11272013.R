###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(BMS) #contains hex2bin and bin2hex
library(bitops)
library(gtools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast)
library(xts)
library(zoo)
library(lubridate)

### Parameters and arguments

in_dir<- "/Users/benoitparmentier/Dropbox/Data/Space_Time"
out_dir<- "/Users/benoitparmentier/Dropbox/Data/Space_Time"
setwd(out_dir)

function_analyses_paper <- "MODIS_and_raster_processing_functions_10182013.R"
script_path <- in_dir #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

#This is the shape file of outline of the study area                                                      #It is an input/output of the covariate script
infile_reg_outline <- "/Users/benoitparmentier/Dropbox/Data/Space_Time/GYRS_MX_trisate_sin_windowed.shp"  #input region outline defined by polygon: Oregon

#ref_rast_name<-""  #local raster name defining resolution, exent, local projection--. set on the fly?? 
ref_rast_name<-"/Users/benoitparmentier/Dropbox/Data/Space_Time/gyrs_sin_mask_1km_windowed.rst"  #local raster name defining resolution, exent: oregon
ref_samp4_name <-"/Users/benoitparmentier/Dropbox/Data/Space_Time/reg_Sample4.rst"
ref_EDGY_name <-"/Users/benoitparmentier/Dropbox/Data/Space_Time/reg_EDGY_mask_sin_1km.rst"
ref_egg_rings_gyr_name <-"/Users/benoitparmentier/Dropbox/Data/Space_Time/reg_egg_rings_gyr.rst"
infile_modis_grid<-"/Users/benoitparmentier/Dropbox/Data/Space_Time/modis_sinusoidal_grid_world.shp" #modis grid tiling system, global
proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_interp <- proj_modis_str

out_suffix <-"11272013" #output suffix for the files that are masked for quality and for 
    
## Other specific parameters
NA_flag_val<- -9999
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
r_stack <- stack(reg_var_list)
levelplot(r_stack,layers=1:12) #show first 12 images (half a year more or less)

ref_rast_r  <- raster(ref_rast_name)
ref_samp4_r <- raster(ref_samp4_name)
mask_EDGY_r <- raster(ref_EDGY_name) #create raster image
egg_rings_gyr_r <- raster(ref_egg_rings_gyr_name)

projection(mask_EDGY_r) <- proj_modis_str #assign projection coord defined earlier
projection(ref_samp4_r) <- proj_modis_str
projection(ref_rast_r) <- proj_modis_str
projection(egg_rings_gyr_r) <- proj_modis_str

r_x <- ref_rast_r #create raster image that will contain x coordinates
r_y <- ref_rast_r #create raster image that will contain y coordiates
values(r_x) <- coordinates(ref_rast_r)[,1] #set values in raster image to x coord
values(r_y) <- coordinates(ref_rast_r)[,2] #set values in raster image to y coord
pix_id_r <- ref_rast_r
values(pix_id_r) <- 1:ncell(ref_rast_r)
s_dat_var <-stack(pix_id_r,r_x,r_y,egg_rings_gyr_r,mask_EDGY_r,ref_samp4_r)
layerNames(s_dat_var) <- c("pix_id_r","r_x","r_y","egg_rings_gyr_r",
                           "mask_EDGY_r","ref_samp4_r")
projection(s_dat_var) <- proj_modis_str
plot(s_dat_var)

## Now extract data for only the EDGY region of Yucatan
EDGY_spdf <- as(mask_EDGY_r,"SpatialPointsDataFrame") #create a SpatialPointsDataFrame
data_EDGY<- extract(r_stack,EDGY_spdf) #extract pixels with NDVI in EDGY area in a matrix
s_dat_var_EDGY <- extract(s_dat_var,EDGY_spdf) #extract pixels with attributes in EDGY area in a matrix
EDGY_dat_spdf<- cbind(data_EDGY,s_dat_var_EDGY) #Add columns from matrix to EDGY 

#EDGY_dat_spdf<- EDGY_spdf

#Save spdf as .RData object
save(EDGY_dat_spdf,file= file.path(out_dir,paste("EDGY_dat_spdf_",out_suffix,".RData",sep="")))

#Save spdf as shapefile, note that column names can be truncated
outfile1<-file.path(out_dir,paste("EDGY_dat_spdf","_",out_suffix,".shp",sep=""))
writeOGR(EDGY_dat_spdf,dsn= dirname(outfile1),layer= sub(".shp","",basename(outfile1)), driver="ESRI Shapefile",overwrite_layer=TRUE)

#Save spdf as delimited csv text file.
outfile1<-file.path(out_dir,paste("EDGY_dat_spdf","_",out_suffix,".txt",sep=""))
write.table(as.data.frame(EDGY_dat_spdf),file=outfile1,sep=",")

#Hurricane August 17, 2007

### NOW ANALYSES WITH TIME AND SPACE...

#filename<-sub(".shp","",infile_reg_outline)             #Removing the extension from file.
#interp_area <- readOGR(dsn=dirname(filename),basename(filename))
#CRS_interp<-proj4string(interp_area)         #Storing the coordinate information: geographic coordinates longlat WGS84

