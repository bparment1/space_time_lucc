####################################    Space Time Analyses PAPER   #######################################
############################  Yucatan case study: SAR, SARMA etc          #######################################
#This script produces a prediction for the dates following the Hurricane event.       
#The script uses spatial neighbour to predict.                        
#AUTHORS: Marco Millones and Benoit Parmentier                                             
#DATE CREATED: 02/07/2014 
#DATE MODIFIED: 03/02/2014
#Version: 2
#PROJECT: GLP Conference Berlin,YUCATAN CASE STUDY with Marco Millones             
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
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

#####  Parameters

#in_dir<-"C:/Users/mmmillones/Dropbox/Space_Time"
in_dir <-"home/parmentier/Data/Space_Time"

#set up the working directory
setwd(in_dir)

Moore_extent_name <- "~/Data/Space_Time/00_moore_clipped_sin_reduced.rst"    
moore_window <- "~/Data/Space_Time/moore_window.rst"

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84

CRS_interp <- proj_modis_str


########### START SCRIPT #################

## Read in data
moore_r <- raster(Moore_extent_name)
moore_w <- raster(moore_window)
projection(moore_w) <- CRS_interp
moore_wgs84_w <- projectRaster(moore_w,crs=CRS(CRS_WGS84),
                               method="ngb")
#loads shapefile data (from shapefile)

moore_w_poly <- rasterToPolygons(moore_wgs84_w, fun=NULL, n=4, na.rm=TRUE, 
                        digits=12, dissolve=FALSE)

moore_w_poly$UNIQID <- 1:nrow(moore_w_poly)

#moore_wgs84_w <- as(moore_wgs84_w,"SpatialGridDataFrame")
#list_nb_moore <- cell2nb(moore_wgs84_w,"queen")
sample1 <- readOGR(dsn=in_dir,layer="sample1") #didn't work
#sample1 <- moore_w_poly
#creates a aditional random variable
sample1$var<-rnorm(n=nrow(sample1),mean=100000,sd=25000) #creates new variables

#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
r_stack <- stack(reg_var_list)
levelplot(r_stack,layers=1:2) #show first 12 images (half a year more or less)

ref_EDGY_name <- file.path(in_dir,"reg_EDGY_mask_sin_1km.rst")
mask_EDGY_r <- raster(ref_EDGY_name) #create raster image
projection(mask_EDGY_r) <- proj_modis_str #assign projection coord defined earlier

## Now extract data for only the EDGY region of Yucatan
EDGY_spdf <- as(mask_EDGY_r,"SpatialPointsDataFrame") #create a SpatialPointsDataFrame
data_EDGY<- extract(r_stack,EDGY_spdf) #extract pixels with NDVI in EDGY area in a matrix

#centroids_sp <- as.data.frame(centroids)
#coordinates(centroids_sp) <- cbind(centroids[,1],centroids[,2])
#proj4string(centroids_sp) <- proj4string(sample1)
#mask_sample1 <- spTransform(centroids_sp,CRS(CRS_interp))

#extract for the subset of EDGY
#ED_spdf <- as(mask_EDGY_r,"SpatialPointsDataFrame") #create a SpatialPointsDataFrame
#

data_sample1 <- extract(r_stack,mask_sample1,df=T,sp=T) #extract pixels with NDVI in EDGY area in a matrix
#Join data to spatial polygon data polygon
sample1$NDVId225 <- data_sample1[,153]
sample1$NDVId241 <- data_sample1[,154]

windows(record=TRUE)

#Univartiate diagnostics of dependent variable
# create original raw histogram and qqnormal plot
hist(sample1$REC_NUM)
qqnorm(sample1$REC_NUM)

#transorm raw data(log transformation) DOESNT WORK
Tsample1$REC_NUM <- log(sample1$REC_NUM+0)

#create tranformed histogram
hist(sample1$trspb)
qqnorm(sample1$trspb)

#calculate shapiro wilk ryan joiner tests of normality
shapiro.test(sample1$REC_NUM)
shapiro.test(sample1$REC_NUM)

#calculate levene test for constant variance (runs but doesnt work yet)
levene.test(sample1$REC_NUM, sample1$REC_NUM,"mean")

#creates the neighborhood file (will take time)

sample1.neigh<-poly2nb(sample1,row.names=sample1$UNIQID,queen=TRUE)
centroid<-cbind(sample1$MU,sample1$MV)
centroids<-coordinates(sample1)
sample1.neigh[[1]]
sam.listb<-nb2listw(sample1.neigh, style="B")
sam.listw<-nb2listw(sample1.neigh, style="W")
sam.listb$weights[[1]]
n<-length(sample1.neigh)
connect.b<-matrix(rep(0,n*n),n,n)
for(i in 1:n){
  connect.b[i,sample1.neigh[[i]]]<-sam.listb$weights[[i]]
}

#calculate moran coefficent,geary ratio, moran scatter plot
moran.test(sample1$REC_NUM, sam.listw)
moran.plot(sample1$REC_NUM, sam.listw)
geary.test(sample1$REC_NUM, sam.listw)


#calculate moran scatter plot manually
wx <- lag.listw(sam.listw, sample1$REC_NUM)
plot(sample1$REC_NUM, wx)
cor(sample1$REC_NUM, wx)

#plots the map with the IDs 
centroid<-cbind(sample1$MU,sample1$MV)
plot(sample1)
text(centroid[,1],centroid[,2],label=sample1$ID)

#plots the map with the neighborhood structure
plot(sample1)
plot(sample1.neigh,centroid,add=T,col="blue")


# OLS model

sam.lm <- lm(REC_NUM ~ MV, data=sample1)
sam.lm <- lm(NDVId225 ~ NDVId225, data=sample1)
summary(sam.lm)

# SAR model
#sam.sar <- errorsarlm(REC_NUM ~ MV, listw=sam.listw, data=sample1)
#summary(sam.sar) had an error message

sam.sar <- errorsarlm(NDVId225 ~ NDVId225, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.sar)

sam.sar <- errorsarlm(NDVId225 ~ var, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.sar)

sam.sar <- errorsarlm(REC_NUM ~ var, listw=sam.listw, data=sample1)
summary(sam.sar)

# AR model
sam.ar <- lagsarlm(NDVId225 ~ NDVId225, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.ar)

#### create back raster imagge

data_sample1$pred_d225 <- sam.ar$fitted.value

predict(sam.ar,sample1$NDVId225,listw=sam.listw)

# Krigging

library(gstat)
syr.v <- variogram(trspb ~ 1, syracuse)
plot(syr.v)
syr.vf.exp <- fit.variogram(syr.v, vgm(1,"Exp", 3000, 0.2))
plot(syr.v, syr.vf.exp)

syr.reg <- spsample(syracuse, 10000, type="regular")
syr.grid <- SpatialPixels(syr.reg)
ok.exp <- krige(trspb ~ 1, syracuse, syr.grid, syr.vf.exp)
color.pal <- colorRampPalette(c("dark red","orange","light Yellow"))
color.palr <- colorRampPalette(c("light yellow","orange","dark red"))
spplot(ok.exp["var1.pred"], col.regions=color.pal)
spplot(ok.exp["var1.var"], col.regions=color.palr)
