generate_outline_and_grid <- function(reg_layer,out_suffix,out_dir){
  
  if(class(reg_layer)=="RasterStack"){
    r <- reg_layer
    r[] <- 1 #change this later
    r <- mask(r,reg_layer)
    reg_sp <- rasterToPolygons(r,dissolve=T)
    reg_layer <- as(reg_sp, "sf")
  }
  if(class(reg_layer)=="SpatialPolygons"){
    reg_layer <- as(reg_layer, "sf")
  }
  ## generate grid: default is 2x2
  reg_grid <- st_make_grid(reg_layer,n=2)
}
    