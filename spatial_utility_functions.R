
#debug(generate_outline_and_grid)
r_test <- generate_outline_and_grid(r_val)
r_test <- generate_outline_and_grid(r_val,n_grid=c(2,1))
tile_sf <- subset(r_test,1)
tile_sf <- r_test[1,]
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
  
  ## save grid figures
  plot(reg_sf)
  plot(reg_grid,border="red",add=T)
  
  reg_sf <- st_sf(tile_id=1:length(reg_grid),reg_grid)
  
  ### buffer: make overlapping area later!!!
  #n_distance <- 5000
  #reg_buffer <- st_buffer(reg_layer,dist=n_distance)
  
  out_filename <- file.path(out_dir,
                            paste0("reg_grid_",out_suffix,".shp"))
  ## overwrite if existing
  st_write(reg_sf,out_filename,delete_dsn =T)
  
  return(reg_sf)
}


    