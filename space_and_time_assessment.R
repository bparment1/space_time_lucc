#######################################    Space Time Analyses Project   #######################################
############################ Assessing Space and Time Model for Space Beats Time framework #######################################
#This script contains functions to assess predictions for Space Beats Time Framework.       
#Temporal predictions use OLS with the image of the previous time step or ARIMA.
#Spatial predictions use spatial regression (lag error model) with different estimation methods (e.g. eigen, chebyshev etc.).
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 11/07/2017 
#DATE MODIFIED: 12/09/2017
#Version: 1

#PROJECT: Space beats time Framework
#TO DO:
#
#COMMIT: more changes and testing of code with tiles space and time predictions assessment
#

#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast) #ARIMA forecasting
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
library(sphet) #contains spreg
library(BMS) #contains hex2bin and bin2hex
library(bitops)

#################################
### Functions used in the script:

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

function_analyses_paper <-"MODIS_and_raster_processing_functions_11172017.R"
script_path <- "/home/bparmentier/Google Drive/Space_beats_time/sbt_scripts"  #path to script functions
source(file.path(script_path,function_analyses_paper)) #source all functions used in this script.

################################
###### Parameters and arguments

#Should use the data that is mosaiced!!

args<-commandArgs(TRUE)

args_table <- args[1]

#args_table <- "/home/bparmentier/Google Drive/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_10292017.csv"
#args_table <- "/home/parmentier/Data/Space_beats_time/Data/input_arguments_sbt_script_NDVI_Rita_11072017.csv"
args_table <- "/home/parmentier/Data/Space_beats_time/Data/input_arguments_sbt_assessment_script_NDVI_Rita_11152017.csv"
df_args <- read.table(args_table,sep=",",stringsAsFactors = FALSE)

### use column 2,3,4 etc.
#index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile1
index_val <- 2 #this is set up for parallelization, if we have multiple regions/tiles, tile 2 rita

in_dir <- df_args[1,index_val]
out_dir <- df_args[2,index_val]
proj_str <- df_args[3,index_val]
file_format <- df_args[4,index_val]
NA_flag_val <- df_args[5,index_val]
out_suffix <- df_args[6,index_val]
create_out_dir_param <- df_args[7,index_val] 
data_fname <- df_args[8,index_val] 
coord_names <- df_args[9,index_val]  
zonal_colnames <- df_args[10,index_val] 
var_names <- df_args[11,index_val] 
num_cores <- df_args[12,index_val] 

n_time_event <- df_args[13,index_val]
time_window_selected <- df_args[14,index_val] 
#previous_step <- df_args[15,index_val] 
date_range <- df_args[15,index_val]  
r_ref <- df_args[16,index_val]  
r_temp_pred <- df_args[17,index_val]
r_spat_pred <- df_args[18,index_val]
method_space <- df_args[19,index_val] 
method_time <- df_args[20,index_val] 
pixel_index <- df_args[21,index_val] 

#P1
in_dir <- "/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017"
r_temp_pred <- list.files(path=in_dir,
                          pattern="r_temp_pred_arima_arima_.*._tile_2_NDVI_Rita_11062017.tif",
                          full.names=T)
out_file <- file.path(in_dir,"raster_temp_files_list_tile_2.txt")
write.table(r_temp_pred,out_file)
#r_spat_pred <- list.files(path=in_dir,
#                                         pattern="r_spat_.*._tile_2_NDVI_Rita_11062017.tif",
#                                         full.names=T)
r_spat_pred <- list.files(path=in_dir,
                                         pattern="r_spat_pred_mle_eigen_no_previous_step_.*._tile_2_NDVI_Rita_11062017.tif",
                                         full.names=T)
out_file <- "/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017/raster_spat_files_list_tile_2.txt"
#/home/parmentier/Data/Space_beats_time/outputs/output_tile_2_NDVI_Rita_11062017
write.table(r_spat_pred,out_file)
#r_spat_pred_mle_eigen_no_previous_step__t_113_tile_2_NDVI_Rita_11062017.tif
s_raster <- "/home/parmentier/Data/Space_beats_time/Data/data_Rita_NDVI/rev_project_output/tile_2/raster_files_list_tile_2.txt"
date_range <- "2001.01.01;2010.12.31;16"
time_window_predicted <- "105;114"
r_zonal <- "crop_r_zonal_rev_tile_2"
r_ref <- NULL
method_space <- "mle;eigen" #method for space and time used to predict space and time respectively
method_time <- "arima;arima;TRUE"
out_suffix <- "assessment_tile_2_NDVI_Rita_11062017"
out_dir <- "output_tile_1_2_combined_NDVI_Rita_11062017"
create_out_dir_param <- FALSE  


######################################################
########################  BEGIN SCRIPT  #############

#Create output directory

if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

#####################################
#### STEP 1:  MOSAIC   ####


if(steps_to_run$mosaic==TRUE){
  
  if(is.null(mosaic_dir)){
    
    #out_dir_tmp <- paste0("mosaic_",out_suffix)
    out_dir_tmp <- paste0("mosaic_output")
    
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    out_dir_s <- file.path(out_dir,out_dir_tmp)
  }else{
    out_dir_s <- mosaic_dir
  }
  
  list_m_var <- vector("list",length(list_tiles_modis))  
  l_df_raster_name <- vector("list",length(list_tiles_modis))  
  
  names(list_m_var)<- list_tiles_modis
  list_m_qc <- vector("list",length(list_tiles_modis))  
  names(list_m_qc)<- list_tiles_modis
  
  for (j in 1:length(list_tiles_modis)){
    #out_suffix_s <- paste(list_tiles_modis[j],"_",sprintf( "%03d", product_version),"_",var_modis_name,"_",out_suffix,sep="")
    #file_format_s <-file_format
    #out_suffix_s <- paste(list_tiles_modis[j],"_",sprintf( "%03d", product_version),"_",var_modis_name,"_",out_suffix,file_format_s,sep="")
    #out_dir_s <- file.path(dirname(out_dir),list_tiles_modis)[j]
    #out_dir_s <- file.path(out_dir,list_tiles_modis)[j]
    #list_m_var[[j]]<-list.files(pattern=paste(out_suffix_s,"$",sep=""),
    #                            path=out_dir_s,
    #                            full.names=TRUE) #inputs for moasics
    file_pattern <- paste0(".*.",product_type,".*.",
                           out_suffix,file_format,"$")
    
    in_dir_tmp <- paste0("mask_qc_",list_tiles_modis[j])
    #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
    in_dir_s <- file.path(out_dir,in_dir_tmp) #input dir is import out dir
    
    #list_r_var_s <- list.files(path=out_dir_s,
    #                           pattern=file_pattern,
    #                           full.names=T)
    
    #MOD11A2_A2012353_h08v05_006_LST_Day_1km_arizona_10092017.rst
    #list_var_mosaiced <-list.files(pattern=".*.LST_Day_1km_arizona_10092017.rst$",
    #                               path=out_dir_s,
    #                               full.names=TRUE) #inputs for moasics
    
    list_r_var_s <-list.files(pattern=file_pattern,
                              #patter=".*.NDVI.*.arizona_10182017.rst$",
                              path=in_dir_s,
                              full.names=TRUE) #inputs for moasics
    list_m_var[[j]] <- list_r_var_s
    #list_m_var[[j]] <-create_raster_list_from_file_pat(out_suffix_s,file_pat="",in_dir=out_dir_s,out_prefix="",file_format=file_format_s)
    df_m_var <- lapply(1:length(list_m_var[[j]]),
                       FUN=extract_dates_from_raster_name,
                       list_files=list_m_var[[j]])
    #df_m_var <- lapply(1:length(list_m_var),
    #                   FUN=extract_dates_from_raster_name,
    #                   list_files=list_m_var)
    df_m_var <- do.call(rbind,df_m_var)
    names(df_m_var) <- c(paste("raster_name",j,sep="_"),"date")
    
    df_m_var[,1] <- as.character(df_m_var[,1])
    l_df_raster_name[[j]] <- df_m_var
    
  }
  
  #test <- merge_all(l_df_raster_name,by="date") #,all.x=T,all.y=T) #does not work properly since df have to be order from most to least complete!!!
  #test <- merge(l_df_raster_name[[1]],l_df_raster_name[[2]],by="date",all.x=T,all.y=T)
  df_m_mosaics <- merge_multiple_df(l_df_raster_name,"date")
  x <- subset(df_m_mosaics,select= -c(date))# drop column with name "date"
  
  #report on missing dates:
  #st <- as.Date(start_date,format="%Y.%m.%d")
  #en <- as.Date(end_date,format="%Y.%m.%d")
  #ll <- seq.Date(st, en, by="1 day")
  #dates_queried <- format(ll,"%Y.%m.%d")
  #mosaic_list_var <-mapply(FUN="c",list_m_var,SIMPLIFY=T)
  #x <-mapply(FUN=as.matrix,list_m_var,SIMPLIFY=T)
  
  #x <-mapply(FUN="c",x,SIMPLIFY=T)
  #MODIS_product <- "MOD13A2.005" #NDVI/EVI 1km product (monthly) #param12
  #strsplit(MODIS_product,"[.]")
  #MODIS_product <- "MOD11A1.005"
  MODIS_product_name <- gsub("[.]","_",MODIS_product)
  date_str <- df_m_mosaics$date
  df_m_mosaics$out_rastnames_var <- paste(MODIS_product_name,date_str,"mosaic",product_type,out_suffix,sep="_") #no file format added!
  mosaic_list_var <- lapply(seq_len(nrow(x)), function(i){x[i,]}) #list of tiles by batch to mosaic
  #Prepare list of output names without extension
  out_rastnames_var <- (basename(gsub(list_tiles_modis[1],"",list_m_var[[1]])))
  out_rastnames_var <- gsub(extension(out_rastnames_var),"",out_rastnames_var)
  #out_rastnames_var <- df_m_mosaics$out_rastnames_var
  
  j <- 1
  
  #out_dir_mosaic <-     
  #out_dir_tmp <- paste0("mosaic_output_",out_suffix)
  out_dir_tmp <- paste0("mosaic_output")
  #out_dir_s <- file.path(out_dir,list_tiles_modis[j])
  out_dir_s <- file.path(out_dir,out_dir_tmp)
  if(!file.exists(out_dir_s)){
    dir.create(out_dir_s)
  }
  
  out_dir_mosaic <- out_dir_s
  #list_param_mosaic<-list(j,mosaic_list_var,out_rastnames_var,out_dir,file_format,NA_flag_val)
  list_param_mosaic<-list(j,
                          mosaic_list_var,
                          out_rastnames_var,
                          out_dir_mosaic,
                          file_format,
                          NA_flag_val)
  
  names(list_param_mosaic)<-c("j",
                              "mosaic_list",
                              "out_rastnames",
                              "out_path",
                              "file_format",
                              "NA_flag_val")
  #debug(mosaic_m_raster_list)
  #list_var_mosaiced <- mosaic_m_raster_list(1,list_param_mosaic)
  #list_var_mosaiced <-mclapply(1:11, 
  #                             list_param=list_param_mosaic, 
  #                             mosaic_m_raster_list,
  #                             mc.preschedule=FALSE,
  #                             mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  #Parallelization
  #started at 22:07 -22:11
  list_var_mosaiced <-mclapply(1:length(mosaic_list_var), 
                               list_param=list_param_mosaic, 
                               mosaic_m_raster_list,
                               mc.preschedule=FALSE,
                               mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  #list_var_mosaiced <- lapply(1:length(mosaic_list_var), list_param=list_param_mosaic, mosaic_m_raster_list) #This is the end bracket from mclapply(...) statement
  
  #r_test <- stack(list_var_mosaiced)
  #plot(r_test,y=1:2)
  
  
}

#####################################
#### STEP 2:  ACCURACY ASSESSMENT   ####

debug(accuracy_space_time_calc)

test <- accuracy_space_time_calc(r_temp_pred=r_temp_pred,
                                 r_spat_pred=r_spat_pred,
                                 s_raster=data_fname,
                                 time_window_predicted=time_window_predicted,
                                 r_zonal=zonal_colnames,
                                 method_space=method_space,
                                 method_time=method_time,
                                 r_ref=r_ref,
                                 out_suffix=out_suffix,
                                 var_names=var_names,
                                 NA_flag_val=NA_flag_val,
                                 file_format=file_format, 
                                 date_range=date_range,
                                 out_dir=out_dir,
                                 create_out_dir_param=create_out_dir_param)
  

#################### END OF SCRIPT ################################################################