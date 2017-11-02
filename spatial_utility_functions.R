
#debug(generate_outline_and_grid)
r_val <- subset(s_raster,1)
r_test <- generate_outline_and_grid(r_val)
debug(generate_outline_and_grid)
r_test <- generate_outline_and_grid(r_val,n_grid=c(2,1))
tile_sf <- subset(r_test,1)
tile_sf <- r_test[1,]

generate_raster_tiles(r_stack,tile_sf,num_cores,out_suffix,out_dir){
  tile_sp <- as(tile_sf,"Spatial")

  
  cropping_raster_stack(i,tile_s_sp,r_stack,file_format,out_dir,out_suffix_s){
    ## This functions crops a stack of raster using a set of polygons.
    ## Polygons are mostly likely representing processing tiles
    tile_poly_sp <- tile_sp[i,]
    
    out_dir_tile <- paste("tile_",i,sep="")
    
    if(!dir.exists(out_dir_tile)){
      dir.create(out_dir_tile)
    }
    #out_dir_tile <- dir.create(paste("tile_",i,sep=""))
    r_s_crop <- crop(r_stack,
                     tile_poly_sp,
                     filename=file.path(out_dir_tile,paste0("crop",file_format)),
                     bylayer=T,
                     suffix=paste(names(r_stack),out_suffix_s),
                     overwrite=T)
    
    lf <- list.files(path=file.path(out_dir,out_dir_tile),
                     pattern=paste0(file_format,"$"),
                     full.names=T)
    out_list_filename <- file.path(out_dir,out_dir_tile,
                                   paste0("raster_files_list_",out_dir_tile,".txt"))
    write.table(lf,out_list_filename,row.names=F,sep=",")
    
    return(out_list_file_name)
  }
  
  mclapply(1:nrow(tile_sp),
           
           FUN=function(i){crop(r_stack,tile_sp[i,])})
  
}
##
#To make overlapping grid, buffer teh region outline by % of area 
##and then generate grid.
#Crop/clip from the original using intersect?

#Or select each new grid feature,
#buffer
#dissolve
#reassemle the features
generate_outline_and_grid <- function(reg_layer,n_grid=NULL,out_suffix="",out_dir="."){
  #This function generates a grid from an input layer
  if(class(reg_layer)=="RasterLayer"){
    r <- reg_layer
    r[] <- 1 #change this later
    r <- mask(r,reg_layer)
    reg_sp <- rasterToPolygons(r,dissolve=T)
    reg_layer <- as(reg_sp, "sf") #this makes a class of "sfc_polygon" "sfc"
  }
  if(class(reg_layer)=="SpatialPolygons"){
    reg_layer <- as(reg_layer, "sf")
  }
  ## generate grid: default is 2x2
  if(is.null(n_grid)){
    n_grid <- c(2,2)
  }
  #reg_grid <- st_make_grid(reg_layer,n=n_grid,offset = n_distance)
  reg_grid <- st_make_grid(reg_layer,n=n_grid)

  reg_sf <- st_sf(tile_id=1:length(reg_grid),reg_grid)
  
  ## save grid figures
  plot(reg_layer)
  plot(reg_grid,border="red",add=T)
  
  ### buffer: make overlapping area later!!!
  #n_distance <- 5000
  #reg_buffer <- st_buffer(reg_layer,dist=n_distance)
  
  out_filename <- file.path(out_dir,
                            paste0("reg_sf_",out_suffix,".shp"))
  ## overwrite if existing
  st_write(reg_sf,out_filename,delete_dsn =T)
  
  return(reg_sf)
}


    