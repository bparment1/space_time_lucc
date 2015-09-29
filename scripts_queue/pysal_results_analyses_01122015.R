############################    SBT comparison of results with pysal   #######################################
################################  Generation of predictions with python  #######################################
#This script produces spatial prediction using spatial regression.
#AUTHORS: Benoit Parmentier                                           
#DATE CREATED:  01/09/2015
#DATE MODIFIED: 01/16/2015
#Version: 1
#PROJECT: Space beats time project

### TO DO LIST:
#compare coefficients
#check random number

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
library(gridExtra)

###### Functions used in this script

pysal_pred_conversion_fun <- function(i,list_param){

  #####
  #This function uses output from pysal script 
  #Inputs:
  #lf_csv: list of files with predicted values and observed from pysal script
  #lf_shp: list of files with shapefiles of location for observed values
  #date: year or date time step predicted
  #out_dir
  #out_suffix : output suffix
  
  ### BEGIN ####
  
  fname <- list_param$lf_csv[i]
  shp_name <- list_param$lf_shp[i]
  date<- list_param$l_date[i]
  out_dir <- list_param$out_dir
  out_suffix_str <- list_param$out_suffix
  r <- list_param$r #raster reference
  
  #data_pred <-read.table(file.path(in_dir,fname,sep=","))
  data_pred <-read.table(fname,sep=",")


  names(data_pred) <- c("ID",
                      paste("res_",date,sep=""),
                      paste("obs_",date,sep=""),
                      paste("rnd_val_",date,sep=""),
                      paste("pred_",date,sep=""))
  #names(data_test) <- c("ID","res_2004","obs_2004","rnd_val","pred_2004")

  dim(data_pred)

  #data_reg <- readOGR(in_dir1,gsub(".shp","",shp_name))
  data_reg <- readOGR(dirname(shp_name),gsub(".shp","",basename(shp_name)))

  head(data_reg)
  head(data_pred)

  data_val<-cbind(data_reg,data_pred)
  head(data_val)

  coordinates(data_val) <- coordinates(data_reg)
  proj4string(data_val) <- proj4string(data_reg)
  #mean(data_tmp$obs_2004)
  mae_val <- mean(abs(data_val[[paste("res_",date,sep="")]]))

  r_res <- rasterize(data_val,r,paste("res_",date,sep=""))
  r_pred <- rasterize(data_val,r,paste("pred_",date,sep=""))
  r_rnd <- rasterize(data_val,r,paste("rnd_val_",date,sep=""))
  r_obs <- rasterize(data_val,r,paste("obs_",date,sep=""))
  
  r_stack <- stack(r_obs,r_pred,r_rnd,r_obs)

  names(r_stack)<-c(paste("res_",date,sep=""),paste("pred_",date,sep=""),
                    paste("rnd_val_",date,sep=""),paste("obs_",date,sep=""))

  raster_name <- paste("r_stack_",date,"_",out_suffix_str,".tif",sep="")
  
  writeRaster(r_stack,filename=file.path(out_dir,raster_name),bylayer=F,overwrite=T)
  #add part to save projected file??
  #return reg_outline!!!
  conv_obj <-list(r_stack,mae_val,file.path(out_dir,raster_name),names(r_stack))
  names(conv_obj) <-c("r_stack","mae_val","raster_name","r_stack_name")
  save(conv_obj,file= file.path(out_dir,paste("conv_obj","_",date,"_",out_suffix_str,".RData",sep="")))

  return(conv_obj)
}

function_spatial_regression_analyses <- "SPatial_analysis_spatial_reg_09252014_functions.R"
script_path <- "/home/parmentier/Data/Space_Time/R" #path to script
#script_path <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/R_workshop_WM_04232014" #path to script
source(file.path(script_path,function_spatial_regression_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/pysal_test/"
in_dir1 <- "/home/parmentier/Data/Space_beats_time/R_Workshop_April2014/output__predictions_09252014"

setwd(in_dir)

lf_csv <-mixedsort(list.files(path=in_dir,pattern="pysal_res_.*._09252014.csv",full.names=T)) #,full.names=T)
lf_shp <-mixedsort(list.files(path=in_dir1,pattern="^r_poly_t.*.shp",full.names=T)) #,full.names=T)
out_dir <- in_dir
l_date <- 2001:2012
out_suffix <- "01212014"

lf_rst <- mixedsort(list.files(path=in_dir1,pattern=".rst",full.names=T))
r<- raster(lf_rst[1])

list_param_pred_conversion <- list(lf_csv,lf_shp,l_date,r,out_dir,out_suffix)
names(list_param_pred_conversion) <- c("lf_csv","lf_shp","l_date","r","out_dir","out_suffix")
#debug(pysal_pred_conversion_fun)
pysal_pred_conversion_fun(1,list_param=list_param_pred_conversion)
pysal_pred_mle <- lapply(1:length(lf_shp),FUN=pysal_pred_conversion_fun,
                        list_param=list_param_pred_conversion)
#Use parallel processing on MAC and Linux/Unix systems
#pred_spat_mle <-mclapply(1:n_pred,FUN=predict_spat_reg_fun,list_param=list_param_spat_reg,,mc.preschedule=FALSE,mc.cores = 2)

#pysal_pred_conversion_fun <- function(i,lf_csv,lf_shp,date,out_suffix){
r_stack <- stack(list.files(pattern="*.tif",path=out_dir,full.names=T))

l_names <- lapply(1:length(pysal_pred_mle),FUN=function(i){pysal_pred_mle[[i]]$r_stack_name})
#l_names <- lapply(1:length(pysal_pred_mle),FUN=function(i){pysal_pred_mle[[i]]$r_stack_name})
names(r_stack) <- unlist(l_names)
#writeRaster(r_stack,
#            filename=file.path(out_dir,paste("pysal_mle_",out_suffix,".rst",sep="")),bylayer=T,suffix=names(r_stack))

writeRaster(r_stack,
            filename=file.path(out_dir,paste("pysal_mle",".rst",sep="")),
            bylayer=T,suffix=paste(names(r_stack),"_",out_suffix,sep=""),
            overwrite=T)

r_res_py_stack <- stack(list.files(pattern="*res.*.rst",path=out_dir,full.names=T))
r_pred_py_stack <- stack(list.files(pattern="*pred.*.rst",path=out_dir,full.names=T))
r_rnd_py_stack <- stack(list.files(pattern="*rnd.*.rst",path=out_dir,full.names=T))
#r_res_stack <- stack(list.files(pattern="*res.*.rst",path=out_dir,full.names=T))

l_names <- lapply(1:length(pysal_pred_mle),FUN=function(i){pysal_pred_mle[[i]]$r_stack_name})

r_res_mle_stack <- stack(list.files(path=in_dir1,pattern="r_spat_res_mle.*.rst",full.names=T))
r_pred_mle_stack <- stack(list.files(path=in_dir1,pattern="r_spat_pred_mle.*.rst",full.names=T))

levelplot(r_pred_mle_stack,col.regions=matlab.like(25))
levelplot(r_pred_py_stack,col.regions=matlab.like(25))

levelplot(r_res_stack,col.regions=matlab.like(25))
levelplot(r_rnd_stack,col.regions=matlab.like(25))

df_py <- cbind(as.data.frame(r_res_stack),as.data.frame(r_pred_stack))
df_mle <- cbind(as.data.frame(r_res_mle_stack),as.data.frame(r_pred_mle_stack))

df_py_pred <- as.data.frame(r_pred_stack)
df_py_res <- as.data.frame(r_res_stack)

df_mle_pred <- as.data.frame(r_pred_mle_stack)
df_mle_res <- as.data.frame(r_res_mle_stack)

### Make this a function later...

l_date <- 2001:2012
l_pred_p <- vector("list",length=length(l_date)) #list of plots containing comparison for res
l_res_p <- vector("list",length=length(l_date))
l_df_val <- vector("list",length=length(l_date))

for(i in 1:length(l_date)){
  date <- l_date[i]
  
  data_pred_tmp <- as.data.frame(cbind(df_py_pred[,i],df_mle_pred[,i]))
  names(data_pred_tmp) <- c("y_py","y_mle")
  data_res_tmp <- as.data.frame(cbind(df_py_res[,i],df_mle_res[,i]))
  names(data_res_tmp) <- c("y_py","y_mle")
  
  #plot(values(subset(r_pred_mle_stack,1)),values(subset(r_pred_py_stack,1)))
  l_pred_p[[i]] <- xyplot(y_mle ~ y_py,data=data_pred_tmp, main=paste("Predicted values for R mle and pysal mle for ",date,sep=""))

  mae_fun <- function(x,y){mean(abs(x-y),na.rm=T)}
  #plot(values(subset(r_res_mle_stack,1)),values(subset(r_res_py_stack,1)))
  l_res_p[[i]] <- xyplot(y_mle ~ y_py,data=data_res_tmp, main=paste("Predicted values for R mle and pysal mle for ",date,sep=""))

  cor_val<- cor(data_res_tmp$y_py,data_res_tmp$y_mle,use="complete.obs")
  mae_val_comp <- mae_fun(data_res_tmp$y_py,data_res_tmp$y_mle)
  mean_val_py <- mean(data_res_tmp$y_py,na.rm=T)
  mean_val_mle<- mean(data_res_tmp$y_mle,na.rm=T)
  mae_val_py <- mean(abs(data_res_tmp$y_py),na.rm=T)
  mae_val_mle<- mean(abs(data_res_tmp$y_mle),na.rm=T)
  
  df_val <- data.frame(x1=cor_val,x2=mae_val_comp,x3=mean_val_py,x4=mean_val_mle,x5=mae_val_py,x6=mae_val_mle)
  names(df_val) <- c("cor_val","mae_val_comp","mean_val_py","mean_val_mle","mae_val_py","mae_val_mle")
  l_df_val[[i]]<- df_val
}

df_val <- do.call(rbind,l_df_val)
df_val$date <- l_date
#Now deal with coefficient...
plot(mae_val_py ~ date, data=df_val,type="b",pch=1,ylim=c(0,1500))
lines(mae_val_mle ~ date, data=df_val,type="b",pch=2)
#points(mae_val_mle ~ date, data=df_val)

plot(cor_val ~ date,data=df_val,type="b",ylim=c(0,1))
plot(mae_val_comp ~ date,data=df_val,type="b",ylim=c(0,1500))
lines(mean_val_py ~ date,data=df_val,type="b",ylim=c(0,1500))
plot(mae_val_comp ~ date,data=df_val,type="b",ylim=c(0,1500))


#Plot of residuals compariing pysal and R mle sp package
layout_m <- c(3,4)

png(paste("Figure_residuals_comparison_pysal_R_mle",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
grid.arrange(l_res_p[[1]],l_res_p[[2]],l_res_p[[3]],
             l_res_p[[4]],l_res_p[[5]],l_res_p[[6]],
             l_res_p[[7]],l_res_p[[8]],l_res_p[[9]],
             l_res_p[[10]],l_res_p[[11]],l_res_p[[12]],
             ncol=3)
dev.off()   

layout_m <- c(3,4)

png(paste("Figure_prediction_comparison_pysal_R_mle",out_suffix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
grid.arrange(l_pred_p[[1]],l_pred_p[[2]],l_pred_p[[3]],
             l_pred_p[[4]],l_pred_p[[5]],l_pred_p[[6]],
             l_pred_p[[7]],l_pred_p[[8]],l_pred_p[[9]],
             l_pred_p[[10]],l_pred_p[[11]],l_pred_p[[12]],
             ncol=3)
dev.off()   

#Get the slope for lambda etc and compare to mle...

txt <- readLines(list.files(pattern="*.txt",path=out_dir,full.names=T)[1])

################### END OF SCRIPT ######################