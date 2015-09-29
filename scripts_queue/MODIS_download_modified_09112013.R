
# title         : MODIS_download.R
# purpose       : Download of MODIS EVI images for the British Colombia;
# reference     : http://spatial-analyst.net/wiki/index.php?title=Download_and_resampling_of_MODIS_images
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, 14 Oct 2010.
# inputs        : Coordinates of the area of interest; proj4 parameters; ftp addresses etc.;
# outputs       : a series of EVI images (geoTIFF) projected in the local coordinate system;
# remarks 1     : To run this script, you need to obtain and install the MODIS resampling tool from [https://lpdaac.usgs.gov/lpdaac/tools/modis_reprojection_tool];
# remarks 2     : You should also obtain the WGET from [http://users.ugent.be/~bpuype/wget/]  --- simply copy the wget exe to windows system folder;
# remarks 3     : make sure you disable your antivirus tools such as Norton or McAfee otherwise it might block wget from running!

library(rgdal)
library(RCurl)
# Obtain the MODIS tool from: http://lpdaac.usgs.gov/landdaac/tools/modis/index.asp
# setwd("E:/PineBeetleBC/MODIS")
# location of the MODIS 1 km monthly blocks:
MOD13A3 <- "ftp://e4ftl01.cr.usgs.gov/MOLT/MOD13A3.005/"
#ftp://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/2000.03.05/
MOD13A3 <- "http://e4ftl01.cr.usgs.gov/MOLT/MOD13A3.005/"
MOD11A1 <- "ftp://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/"
#MOD13A3 <-  "http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/"
#MOD13A3a <- "ftp://anonymous:test@e4ftl01u.ecs.nasa.gov/MOLT/MOD13A3.005/"
MOD13A3a <- "ftp://anonymous:test@e4ftl01.cr.usgs.gov/MOLT/MOD13A3.005/"
# location of the mosiacing tool:
MRT <- 'E:\\MODIS\\MRT\\bin\\'
workd <- 'E:\\PineBeetleBC\\MODIS\\'
options(download.file.method="auto")

in_dir<- "/Users/benoitparmentier/Dropbox/Data/NCEAS/MODIS_processing"
#ok this works for now...
download.file("http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/2000.07.11/MOD11A1.A2000193.h12v04.005.2007200005030.hdf",
              destfile="./test.hdf")
download.file("http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/2000.07.11/MOD11A1.A2000193.h12v04.005.*.hdf",
              destfile="./test.hdf")

list.files(path="http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/")
# get the list of directories (thanks to Barry Rowlingson):
items <- strsplit(getURL(MOD13A3), "\n")[[1]]
items <- strsplit(getURL(MOD11A1), "\n")[[1]]

getURLContent("http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/2000.07.11/") #This works to list content!!!

url<-"http://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.005/"
filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
# items[2]
# [1] "drwxr-xr-x  2 90 118784 Jan  5  2009 2000.02.01\r"
# you get the folders (and files) but the folder names are in the form of a unix directory listing
# get the last word of any lines that start with 'd':
folderLines <- items[substr(items, 1, 1)=='d']
# get the directory names and create a new data frame:
dirs <- unlist(lapply(strsplit(folderLines, " "), function(x){x[length(x)]}))
dates <- data.frame(dirname=unlist(strsplit(dirs, "\r")))

# get the list of *.hdf files:
dates$BLOCK1 <- rep(NA, length(dates$dirname))
dates$BLOCK2 <- rep(NA, length(dates$dirname))
dates$BLOCK3 <- rep(NA, length(dates$dirname))
dates$BLOCK4 <- rep(NA, length(dates$dirname))
dates$BLOCK5 <- rep(NA, length(dates$dirname))
dates$BLOCK6 <- rep(NA, length(dates$dirname))
dates$BLOCK7 <- rep(NA, length(dates$dirname))
dates$BLOCK8 <- rep(NA, length(dates$dirname))
dates$BLOCK9 <- rep(NA, length(dates$dirname))

for (i in 9:length(dates$dirname)){
getlist <- strsplit(getURL(paste(MOD13A3, dates$dirname[[i]], "/", sep=""), .opts=curlOptions(ftplistonly=TRUE)), "\r\n")[[1]]
BLOCK1 <- getlist[grep(getlist, pattern="MOD13A3.*.h09v03.*.hdf")[1]]
BLOCK2 <- getlist[grep(getlist, pattern="MOD13A3.*.h09v04.*.hdf")[1]]
BLOCK3 <- getlist[grep(getlist, pattern="MOD13A3.*.h10v02.*.hdf")[1]]
BLOCK4 <- getlist[grep(getlist, pattern="MOD13A3.*.h10v03.*.hdf")[1]]
BLOCK5 <- getlist[grep(getlist, pattern="MOD13A3.*.h10v04.*.hdf")[1]]
BLOCK6 <- getlist[grep(getlist, pattern="MOD13A3.*.h11v02.*.hdf")[1]]
BLOCK7 <- getlist[grep(getlist, pattern="MOD13A3.*.h11v03.*.hdf")[1]]
BLOCK8 <- getlist[grep(getlist, pattern="MOD13A3.*.h12v03.*.hdf")[1]]
BLOCK9 <- getlist[grep(getlist, pattern="MOD13A3.*.h12v02.*.hdf")[1]]

# write up the file names back to the dates.txt:
for(j in 2:10){
   dates[i,j] <- get(paste("BLOCK", j-1, sep=""))
}

# Download all blocks from the list to a local drive:
# while(!is.na(dates[i,2])&!is.na(dates[i,3])&!is.na(dates[i,4])&!is.na(dates[i,5])&!is.na(dates[i,6])&!is.na(dates[i,7])&!is.na(dates[i,8])&!is.na(dates[i,9])&!is.na(dates[i,10])){
download.file(paste(MOD13A3a, dates$dirname[[i]], "/", BLOCK1,sep=""), destfile=paste(getwd(), "/", BLOCK1, sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK2,sep=""), destfile=paste(getwd(), "/", BLOCK2,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK3,sep=""), destfile=paste(getwd(), "/", BLOCK3,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK4,sep=""), destfile=paste(getwd(), "/", BLOCK4,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK5,sep=""), destfile=paste(getwd(), "/", BLOCK5,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK6,sep=""), destfile=paste(getwd(), "/", BLOCK6,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK7,sep=""), destfile=paste(getwd(), "/", BLOCK7,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK8,sep=""), destfile=paste(getwd(), "/", BLOCK8,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)
download.file(paste(MOD13A3, dates$dirname[[i]], "/", BLOCK9,sep=""), destfile=paste(getwd(), "/", BLOCK9,sep=""), mode='wb', method='wget', quiet=T, cacheOK=FALSE)

# remove "." from the file name:
dirname1 <- sub(sub(pattern="\\.", replacement="_", dates$dirname[[i]]), pattern="\\.", replacement="_", dates$dirname[[i]])
# mosaic the blocks:
mosaicname = file(paste(MRT, "TmpMosaic.prm", sep=""), open="wt")
write(paste(workd, BLOCK1, sep=""), mosaicname)
write(paste(workd, BLOCK2, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK3, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK4, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK5, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK6, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK7, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK8, sep=""), mosaicname, append=T)
write(paste(workd, BLOCK9, sep=""), mosaicname, append=T)
close(mosaicname)
# generate temporary mosaic:
shell(cmd=paste(MRT, 'mrtmosaic -i ', MRT, 'TmpMosaic.prm -s "0 1 0 0 0 0 0 0 0 0 0" -o ', workd, 'TmpMosaic.hdf', sep=""))

# resample to epsg=3005:
filename = file(paste(MRT, "mrt", dirname1, ".prm", sep=""), open="wt")
write(paste('INPUT_FILENAME = ', workd, 'TmpMosaic.hdf', sep=""), filename) 
# write(paste('INPUT_FILENAMES = ( ', workd, BLOCK1, ' ', workd, BLOCK2, ' ', workd, BLOCK3, ' ', workd, BLOCK4, ' ', workd, BLOCK5, ' ', workd, BLOCK6, ' ', workd, BLOCK7, ' ', workd, BLOCK8, ' ', workd, BLOCK9, ' )', sep=""), filename)  # unfortunatelly does not work via command line  :(
write('  ', filename, append=TRUE) 
# write('SPECTRAL_SUBSET = ( 0 1 0 0 0 0 0 0 0 0 0 )', filename, append=TRUE)
write('SPECTRAL_SUBSET = ( 1 )', filename, append=TRUE)
write('  ', filename, append=TRUE)
write('SPATIAL_SUBSET_TYPE = OUTPUT_PROJ_COORDS', filename, append=TRUE)
write('  ', filename, append=TRUE)
write('SPATIAL_SUBSET_UL_CORNER = ( 637278.0 1701350.0 )', filename, append=TRUE)
write('SPATIAL_SUBSET_LR_CORNER = ( 1907278.0 335350.0 )', filename, append=TRUE)
write('  ', filename, append=TRUE)
write(paste('OUTPUT_FILENAME = ', workd, 'tmp', dirname1, '.tif', sep=""), filename, append=TRUE)
write('  ', filename, append=TRUE)
write('RESAMPLING_TYPE = NEAREST_NEIGHBOR', filename, append=TRUE)
write('  ', filename, append=TRUE)
write('OUTPUT_PROJECTION_TYPE = AEA', filename, append=TRUE)
write('  ', filename, append=TRUE)
write('OUTPUT_PROJECTION_PARAMETERS = ( ', filename, append=TRUE)
write(' 0.0 0.0 50.0', filename, append=TRUE)
write(' 58.5 -126.0 45.0', filename, append=TRUE)
write(' 1000000.0 0.0 0.0', filename, append=TRUE)
write(' 0.0 0.0 0.0', filename, append=TRUE)
write(' 0.0 0.0 0.0 )', filename, append=TRUE)
write('  ', filename, append=TRUE)
write('DATUM = NAD83', filename, append=TRUE)
write('  ', filename, append=TRUE)
write('OUTPUT_PIXEL_SIZE = 1000', filename, append=TRUE)
write('  ', filename, append=TRUE)
close(filename)

# Mosaic the images to get the whole area:
shell(cmd=paste(MRT, 'resample -p ', MRT, 'mrt', dirname1, '.prm', sep=""))
# delete all hdf files!
unlink(paste(getwd(), '/', BLOCK1, sep=""))
unlink(paste(getwd(), '/', BLOCK2, sep=""))
unlink(paste(getwd(), '/', BLOCK3, sep=""))
unlink(paste(getwd(), '/', BLOCK4, sep=""))
unlink(paste(getwd(), '/', BLOCK5, sep=""))
unlink(paste(getwd(), '/', BLOCK6, sep=""))
unlink(paste(getwd(), '/', BLOCK7, sep=""))
unlink(paste(getwd(), '/', BLOCK8, sep=""))
unlink(paste(getwd(), '/', BLOCK9, sep=""))
}
#}

# Check the validity:
GDALinfo("tmp2000_02_01.1_km_monthly_EVI.tif")


# end of script!


###### other script for http url

#http://pastebin.com/aiw2qtSj#
#R code to download MODIS via HTTP
#BY: A GUEST ON AUG 3RD, 2013  |  SYNTAX: NONE  |  SIZE: 10.58 KB  |  HITS: 46  |  EXPIRES: NEVER
#DOWNLOAD  |  RAW  |  EMBED  |  REPORT ABUSE  |  PRINT


##### Initialization
cat("Initializing..\n")
rm(list=ls())
test=function() {source("D:\\Jeff\\modis\\script\\downloadMODIS_V1B5.r")}

### change defaults
options(download.file.method="auto")
options(stringsAsFactors=FALSE)

### set user params
# computer params
inpPath="D:\\Jeff\\modis\\inputs\\brisbaneCoord.txt"
wDir="D:\\Jeff\\modis\\temp"
expTifDir="D:\\Jeff\\modis\\outputs\\Geo_TIFF"
expJpgDir="D:\\Jeff\\modis\\outputs\\preview_jpg"
expMetaDir="D:\\Jeff\\modis\\outputs\\metadata"
expMetaPath="D:\\Jeff\\modis\\outputs\\metadata.txt"
mrtPath="C:\\Users\\uqjhans4\\Documents\\MRT\\bin"
mshpPath="D:\\Jeff\\modis\\inputs\\MODIS_GRID_GEO.shp"

## modis params
dataType="MOD13Q1.005"
jpgBrowseName="MOD13A1"
bands="1,0,0,0,0,0,0,0,0,0,0,0" # gets NDVI
esriPrjFile="C:\\Program Files (x86)\\ArcGIS\\Desktop10.0\\Coordinate Systems\\Projected Coordinate Systems\\World\\Mercator (world).prj"

### load in dependencies
require(RCurl)
require(rgdal)
require(rgeos)
require(stringr)
require(XML)

### define built-in functions
wgs84crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

endswith=function(x, char) {
  currSub = substr(x, as.numeric(nchar(x)-nchar(char))+1,nchar(x))
  if (currSub==char) {return(TRUE)}
  return(FALSE)
}

extractFolders=function(urlString) {
  htmlString=getURL(urlString)
  ret=gsub("]", "", str_replace_all(str_extract_all(htmlString, paste('DIR',".([^]]+).", '/\">',sep=""))[[1]], "[a-zA-Z\"= <>/]", ""))
  return(ret[which(nchar(ret)>0)])
}

extractFiles=function(urlString, selHV) {
  # get filename strings
  htmlString=getURL(urlString)
  allVec=gsub('\">', '', gsub('<a href=\"', "", str_extract_all(htmlString, paste('<a href=\"',"([^]]+)", '\">',sep=""))[[1]]))
  ret=c()
  for (currSel in selHV) {
    ret=c(ret, grep(currSel, allVec, value=TRUE))
  }
  # select specific files
  jpg=sapply(ret, FUN=endswith, char=".jpg")
  xml=sapply(ret, FUN=endswith, char=".xml")
  hdf=sapply(ret, FUN=endswith, char=".hdf")
  retList=list(jpg=ret[which(jpg)], xml=ret[which(xml)], hdf=ret[which(hdf)])
  return(retList)
}

generatePRM=function(inFile, outFile) {
  ### generates PRM file to resample to wgs 1984 CRS
  filename = paste(wDir, "\\mrt.prm", sep="")
  cat(paste('INPUT_FILENAME = ', inFile, "\n", sep=""), file=filename)
  cat(paste('SPECTRAL_SUBSET = (',bands,')\n',sep=""), file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  cat('SPATIAL_SUBSET_TYPE = INPUT_LAT_LONG\n', file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  cat('SPATIAL_SUBSET_UL_CORNER = ( -27.25 152.665 )\n', file=filename, append=TRUE)
  cat('SPATIAL_SUBSET_LR_CORNER = ( -27.66 153.200 )\n', file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  cat(paste('OUTPUT_FILENAME = ', outFile, "\n", sep=""), file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  cat('RESAMPLING_TYPE = NN\n', file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  
  # mercator
  cat('OUTPUT_PROJECTION_TYPE = MERCAT\n', file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  cat('DATUM = WGS84\n', file=filename, append=TRUE)
  cat('  \n', file=filename, append=TRUE)
  cat('OUTPUT_PROJECTION_PARAMETERS = ( \n', file=filename, append=TRUE)
  cat(' 0.0 0.0 0.0\n', file=filename, append=TRUE)
  cat(' 0.0 0.0 1.0\n', file=filename, append=TRUE)
  cat(' 0.0 0.0 0.0\n', file=filename, append=TRUE)
  cat(' 0.0 0.0 0.0\n', file=filename, append=TRUE)
  cat(' 0.0 0.0 0.0 )\n', file=filename, append=TRUE)            
}


parseXML2DF=function(inFile) {
  # general formatting for xml
  xmlTree = xmlTreeParse(inFile)
  xmlParse = unlist(xmlTree, recursive=TRUE)
  xmlSub = xmlParse[which(sapply(names(xmlParse), FUN=endswith, char="value"))]
  names(xmlSub) = str_replace_all(names(xmlSub), c("doc"), "")
  names(xmlSub) = str_replace_all(names(xmlSub), c("children"), "")
  names(xmlSub) = str_replace_all(names(xmlSub), c("text"), "")
  names(xmlSub) = str_replace_all(names(xmlSub), c("value"), "")
  names(xmlSub) = gsub(".", "_", names(xmlSub), fixed=TRUE)
  names(xmlSub) = gsub("^_*|(?<=_)_|_*$", "", names(xmlSub), perl=TRUE)
  xmlDF=as.data.frame(matrix(xmlSub, nrow=1))
  colnames(xmlDF)=names(xmlSub)
  xmlDF=xmlDF[,-grep("GranuleURMetaData_InputGranule_InputPointer",names(xmlDF))]
  # properly parse the PSA fields
  PSANameFieldsPos= grep("GranuleMetaDataFile_GranuleURMetaData_PSAs_PSA_PSAName", names(xmlDF))
  PSANameFieldsValue = unlist(xmlDF[1,PSANameFieldsPos], recursive=TRUE)
  PSAValueFieldsPos=grep("GranuleMetaDataFile_GranuleURMetaData_PSAs_PSA_PSAValue", names(xmlDF))
  PSAValueFieldsValue = unlist(xmlDF[1,PSAValueFieldsPos], recursive=TRUE)               
  xmlDF=xmlDF[,-c(PSANameFieldsPos, PSAValueFieldsPos)]
  for (currPos in seq_along(PSANameFieldsPos)) {
    xmlDF[PSANameFieldsValue[currPos]] = PSAValueFieldsValue[currPos]
  }
  # more general field name parsing
  names(xmlDF)=gsub("GranuleMetaDataFile_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub("GranuleURMetaData_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub("SpatialDomainContainer_HorizontalSpatialDomainContainer_GPolygon_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub("SpatialDomainContainer_HorizontalSpatialDomainContainer_GPolygon_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub("MeasuredParameter_MeasuredParameterContainer_QAStats_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub("MeasuredParameter_MeasuredParameterContainer_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub("QAFlags_", "", names(xmlDF), fixed=TRUE)
  names(xmlDF)=gsub(".", "_", names(xmlDF), fixed=TRUE)
  # properly parse the QAStats fields
  PSA_nfp=grep("ParameterName", names(xmlDF))
  PSA_nfv=unlist(xmlDF[,PSA_nfp])
  xmlDF=xmlDF[,-PSA_nfp]
  PSA_pfp1=grep("QAPercentMissingData",names(xmlDF))
  PSA_pfp2=grep("QAPercentOutofBoundsData",names(xmlDF))
  PSA_pfp3=grep("QAPercentInterpolatedData",names(xmlDF))
  PSA_pfp4=grep("QAPercentCloudCover",names(xmlDF))
  PSA_pfp5=grep("AutomaticQualityFlag",names(xmlDF))
  PSA_pfp6=grep("AutomaticQualityFlagExplanation",names(xmlDF))
  PSA_pfp7=grep("ScienceQualityFlag",names(xmlDF))
  PSA_pfp8=grep("ScienceQualityFlagExplanation",names(xmlDF))
  for (currPos in seq_along(PSA_nfp)) {
    currVal=gsub(" ", "_", PSA_nfv[currPos], fixed=TRUE)
    names(xmlDF)[PSA_pfp1[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp1[currPos]],sep="_")
    names(xmlDF)[PSA_pfp2[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp2[currPos]],sep="_")
    names(xmlDF)[PSA_pfp3[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp3[currPos]],sep="_")
    names(xmlDF)[PSA_pfp4[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp4[currPos]],sep="_")
    names(xmlDF)[PSA_pfp5[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp5[currPos]],sep="_")
    names(xmlDF)[PSA_pfp6[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp6[currPos]],sep="_")
    names(xmlDF)[PSA_pfp7[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp7[currPos]],sep="_")
    names(xmlDF)[PSA_pfp8[currPos]] = paste("B",currVal,names(xmlDF)[PSA_pfp8[currPos]],sep="_")
  }
  
  return(xmlDF)
}

rbind_cust=function(x,y) {
  if (paste(names(x), collapse="_") == paste(names(y), collapse="_")) {return(rbind(x,y))}
  xUniCols=names(x)[which(!names(x) %in% names(y))]
  yUniCols=names(y)[which(!names(y) %in% names(x))]
  if (length(xUniCols)>0) {
    for (xi in xUniCols) {y[xi]=NA}
  }
  if (length(yUniCols)>0) {
    for (yi in yUniCols) {x[yi]=NA}
  }
  return(rbind(x,y))
}

#### Preliminary processing
cat("Prpearing data for processing..\n")
### load in data
inpDF=read.table(inpPath, sep=",", header=TRUE, as.is=TRUE)
mrtSHP=suppressWarnings(readOGR(dirname(mshpPath), gsub(".shp", "", basename(mshpPath), fixed=TRUE), verbose=FALSE))

### get list of HV numbers
mrtSHP@data$hv = paste("h", formatC(mrtSHP@data$h, width = 2, format = "d", flag = "0"), "v", formatC(mrtSHP@data$v, width = 2, format = "d", flag = "0"), sep="")
inpSHP=SpatialPoints(coords=data.frame(inpDF$lon, inpDF$lat), proj4string=wgs84crs)
selRows=gContains(mrtSHP, inpSHP, byid=TRUE)
selHVs=unique(mrtSHP@data[which(selRows),"hv"])

### get directories
mainDir=paste("http://e4ftl01.cr.usgs.gov/MOLT/",dataType,"/",sep="")
temporalDirs=paste(mainDir,extractFolders(mainDir),"/",sep="")
dlDir=c(expJpgDir,expMetaDir,wDir)

#### Main processing
### download files
xmlPaths=c()
jpgPaths=c()
for (currDirPos in seq_along(temporalDirs)[1:5]) {
  
  cat(paste("Starting time increment ",currDirPos,"out of",length(temporalDirs),"\n"))
  # get files
  allFiles=extractFiles(temporalDirs[currDirPos], selHVs)
  inFiles=list(jpg=allFiles[[1]], xml=allFiles[[2]], hdf=allFiles[[3]])
  outFiles=list()
  
  # download files
  cat("\tDownloading files..\n")
  for (currFileTypePos in seq_len(3)) {
    outFiles[[currFileTypePos]]=gsub(".", "_", basename(inFiles[[currFileTypePos]]), fixed=TRUE)
    outFiles[[currFileTypePos]]=paste(dlDir[[currFileTypePos]],"\\",gsub(paste("_", names(inFiles)[currFileTypePos], sep=""), paste(".", names(inFiles)[currFileTypePos], sep=""), outFiles[[currFileTypePos]], fixed=TRUE),sep="")
    download.file(url=paste(temporalDirs[currDirPos],inFiles[[currFileTypePos]],sep="")[1], destfile=outFiles[[currFileTypePos]][1], mode="wb", quiet=TRUE)
  }
  jpgPaths=c(jpgPaths,outFiles[[1]][1])
  xmlPaths=c(xmlPaths,outFiles[[2]])
  
  ### reformat and project .hdf files
  cat("\tReformatting and projecting files..\n")
  pb=txtProgressBar(0, length(inFiles[[3]]), 0 ,style=3)
  for (currFilePos in seq_along(inFiles[[3]])) {
    # generate prm file
    inFile=outFiles[[3]][[currFilePos]]
    outFile=paste(expTifDir, "\\", gsub(".hdf", ".tif", basename(inFile)), sep="")
    prjFile=gsub(".tif", ".prj", outFile)
    generatePRM(inFile, outFile)
    
    # execute modis resample tool
    system(paste(mrtPath,"\\resample.exe -p ",wDir,"\\mrt.prm",sep=""), show.output.on.console=FALSE)
    
    # copy projection file
    file.copy (esriPrjFile, prjFile, overwrite=TRUE)
    
    # update progress bar
    setTxtProgressBar(pb, currFilePos)
  }
  close(pb)
}

### merge all metadatafiles into a single csv
cat("Generating metadata file..\n")
mergedMDF=parseXML2DF(xmlPaths[1])
pb2=txtProgressBar(0, length(xmlPaths), 1, style=3)
counter=1
for (currXml in xmlPaths[-1]) {
  # load xml file
  currDF=parseXML2DF(currXml)
  
  # add parsed file to df
  mergedMDF=rbind_cust(mergedMDF, currDF)
  
  # update progress bar
  counter=counter+1
  setTxtProgressBar(pb2, counter)
}
close(pb2)

# add in image preview links for use in excel
mergedMDF$preview=paste('=HYPERLINK("',jpgPaths,'", "LINK")', sep="")

# save complete metadata file
write.table(mergedMDF, expMetaPath, sep="\t", row.names=FALSE, quote=FALSE)

# print finishing message
cat("DONE!\n")
create a new version of this paste RAW Paste Data


