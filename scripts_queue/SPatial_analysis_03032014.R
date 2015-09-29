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
#centroids_sp <- gCentroid(moore_w_poly,byid=TRUE)

#moore_wgs84_w <- as(moore_wgs84_w,"SpatialGridDataFrame")
#list_nb_moore <- cell2nb(moore_wgs84_w,"queen")
sample1 <- readOGR(dsn=in_dir,layer="sample1") #didn't work
#sample1 <- moore_w_poly
centroids_sp <- gCentroid(sample1,byid=TRUE) #centroids of poly

mask_sample1 <- spTransform(centroids_sp,CRS(CRS_interp))
#reads input NDVI time series images
reg_var_list <- list.files(path=in_dir,pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")
#reg_mosaiced_MOD13A2_A2012353__005_1_km_16_days_NDVI_09242013_09242013.rst
#reg_var_list <- list.files(pattern="reg_mosaiced_MOD13A2_A.*.__005_1_km_16_days_NDVI_09242013_09242013.rst$")

r_stack <- stack(reg_var_list)
levelplot(r_stack,layers=1:2) #show first 12 images (half a year more or less)

data_sample1 <- extract(r_stack,mask_sample1,df=T,sp=T) #extract pixels with NDVI in EDGY area in a matrix
#Join data to spatial polygon data polygon

projection(r_stack) <- CRS_interp
r_var_w <- crop(r_stack,moore_w)
#sample1 <- moore_w_poly 

r_var_dates <- subset(r_var_w,152:156)
r_var_dates <- mask(r_var_dates,moore_w)
plot(r_var_dates)

r1 <- subset(r_var_dates,1)
r2 <- subset(r_var_dates,2)

r_NA <-r1+r2 #this contains NA to mask values...

r1 <- mask(r1,r_NA)
r2 <- mask(r2,r_NA)

r_s <- stack(r1,r2)
names(r_s)<- c("v1","v2")
r_s <- projectRaster(r_s,crs=CRS_WGS84)

r_poly <- rasterToPolygons(r_s, fun=NULL, n=4, na.rm=TRUE, 
                        digits=12, dissolve=FALSE)

r_poly$UNIQID <- 1:nrow(r_poly)
r_nb <-poly2nb(r_poly,row.names=r_poly$UNIQID,queen=TRUE)
#sam.listw<-nb2listw(sample1.neigh, style="W",zero.policy=TRUE)

r_listw<-nb2listw(r_nb, style="W",zero.policy=TRUE)

# SAR model
sam.sar <- errorsarlm(v2 ~ v1, listw=r_listw, 
                      data=r_poly,tol.solve=1e-36,zero.policy=TRUE)
summary(sam.sar)

sam.sar <- errorsarlm(NDVId225 ~ var, listw=sam.listw, 
                      data=sample1)

sam.sar <- errorsarlm(NDVId225 ~ var, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.sar)

sam.sar <- errorsarlm(REC_NUM ~ var, listw=sam.listw, data=sample1)
summary(sam.sar)

# AR model
sam.ar <- lagsarlm(NDVId225 ~ NDVId225, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.ar)



######

sample1$var <- rnorm(n=nrow(sample1),mean=100000,sd=25000) #creates new variables

data_df <- as.data.frame(data_sample1)
#Should use merge with ID!!
sample1$NDVId209 <- data_df[,152] #Before hurricane
sample1$NDVId225 <- data_df[,153] #Before hurricane
sample1$NDVId241 <- data_df[,154] #After hurricane

#windows(record=TRUE)

#Univartiate diagnostics of dependent variable
# create original raw histogram and qqnormal plot
hist(sample1$REC_NUM)
qqnorm(sample1$REC_NUM)

#transorm raw data(log transformation) DOESNT WORK
sample1$REC_NUM <- log(sample1$REC_NUM+0)

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
)
#sam.listw<-nb2listw(sample1.neigh, style="W",zero.policy=TRUE)
sam.listw<-nb2listw(sample1.neigh, style="W")

sam.listb$weights[[1]]

sam.listw$weights[[1]]

n<-length(samplcentroids_sp <- gCentroid(sample1,byid=TRUE) #centroids of poly
#e1.neigh)
connect.b<-matrix(rep(0,n*n),n,n)

for(i in 1:n){
  connect.b[i,sample1.neigh[[i]]]<-sam.listb$weights[[i]]
}

#calculate moran coefficent,geary ratio, moran scatter plot
moran.test(sample1$NDVId225, sam.listw)
moran.plot(sample1$NDVId225, sam.listw)
geary.test(sample1$NDVId225, sam.listw)

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
sam.sar <- errorsarlm(REC_NUM ~ MV, listw=sam.listw, 
                      data=sample1,tol.solve=1e-36)
summary(sam.sar) #had an error message

sam.sar <- errorsarlm(NDVId225 ~ NDVId225, listw=sam.listw, 
                      data=sample1,tol.solve=1e-36,zero.policy=TRUE)
summary(sam.sar)

sam.sar <- errorsarlm(NDVId225 ~ var, listw=sam.listw, 
                      data=sample1)

sam.sar <- errorsarlm(NDVId225 ~ var, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.sar)

sam.sar <- errorsarlm(REC_NUM ~ var, listw=sam.listw, data=sample1)
summary(sam.sar)

# AR model
sam.ar <- lagsarlm(NDVId225 ~ NDVId225, listw=sam.listw, data=sample1,tol.solve=1e-36)
summary(sam.ar)

#### create back raster imagge

data_sample1$pred_d225 <- sam.ar$fitted.value
data_sample1$pred_d225 <- sam.sar$fited.value

r_pred <-rasterize(x=data_sample1,y=)
#predict(sam.ar,sample1$NDVId225,listw=sam.listw)


