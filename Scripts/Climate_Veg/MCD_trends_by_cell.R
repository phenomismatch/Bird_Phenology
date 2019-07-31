library(dggridR)
library(ggplot2)
library(viridis)
library(gdalUtils)
library(raster)
user <- 'JacobSocolar'
setwd(paste0('/Users/', user, '/Desktop/useful_datasets/Veg_phenology/MCD12Q2/'))
years <- c(2001:2016)
jdayMap <- 365*c(1:14) + 1 + 1*(c(2000:2013) > 2004) + 1*(c(2000:2013) > 2008) + 1*(c(2000:2013) > 2012)


# # Get the start of season 1 layer and save as GEOTIFF
# hdfs <- dir(pattern = ".hdf")
# 
# for(i in 1:length(hdfs)){
#   print(i)
#   sds <- get_subdatasets(hdfs[i])
#   gdal_translate(sds[1], dst_dataset = paste0('onset_inc_', i, '.tif'))
#   gdal_translate(sds[2], dst_dataset = paste0('onset_max_', i, '.tif'))
# }

# Import GEOTIFFs as rasters
tiffs <- dir(pattern = ".tif")
xmls <- dir(pattern = ".xml")
tiffs <- tiffs[-which(tiffs %in% xmls)]
tiffs.inc <- tiffs[grep('_inc_', tiffs)]
tiffs.max <- tiffs[grep('_max_', tiffs)]

inds <- c(1:713)[order(as.character(1:713))]
yrs <- vector()

for(i in 1:713){
  yrs[i] <- as.numeric(substr(hdfs[inds[i]], 10,13))
}

for(y in years){
  TI <- tiffs.inc[which(yrs==y)]
  raster1 <- raster::aggregate(raster(TI[1]), fact=10, na.rm = T)
  for(i in 2:length(TI)){
    print(paste("inc",y,i))
    raster1 <- raster::merge(raster1, raster::aggregate(raster(TI[i]), fact=10, na.rm = T), tolerance=10^-6)
  }
  assign(paste0("inc_",y), raster1)
  df <- as.data.frame(raster1, xy=T)
  df$year <- y
  if(y == 2001){
    inc_data <- df
  }else{
    inc_data <- rbind(inc_data, df)
  }
}

for(y in years){
  TI <- tiffs.max[which(yrs==y)]
  raster1 <- raster::aggregate(raster(TI[1]), fact=10, na.rm = T)
  for(i in 2:length(TI)){
    print(paste("max",y,i))
    raster1 <- raster::merge(raster1, raster::aggregate(raster(TI[i]), fact=10, na.rm = T), tolerance=10^-6)
  }
  assign(paste0("max_",y), raster1)
  df <- as.data.frame(raster1, xy=T)
  df$year <- y
  if(y == 2001){
    max_data <- df
  }else{
    max_data <- rbind(max_data, df)
  }
}

mcd12Q2_bounded <- inc_data[which(!is.na(inc_data$layer)), ]
mcd_sp <- SpatialPointsDataFrame(mcd12Q2_bounded[,c(1,2)], mcd12Q2_bounded[,c(3,4)], proj4string = CRS(proj4string(raster1)))
mcd_sp2 <- spTransform(mcd_sp, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
mcd3 <- as.data.frame(mcd_sp2)
mcd3$jday <- mcd3$layer - jdayMap[mcd3$year - 2000]

# Associate each pixel with the cell of a hexagonal grid
hexgrid7 <- dgconstruct(res=7) # Construct geospatial hexagonal grid
mcd3$cell7 <- dgGEO_to_SEQNUM(hexgrid7, mcd3$x, mcd3$y)$seqnum
names(mcd3) <- c("raw_day","year","x","y","jday","cell")

# Get the mean onset greenness for each cell-year
hex7_mcd12Q2 <- doBy::summaryBy(jday~cell+year,data=mcd3, na.rm=T)
save(hex7_mcd12Q2, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex7_mcd12Q2.Rdata")

# Associate each pixel with the cell of a hexagonal grid
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
mcd3$cell6 <- dgGEO_to_SEQNUM(hexgrid6, mcd3$x, mcd3$y)$seqnum

# Get the mean onset greenness for each cell-year
hex6_mcd12Q2 <- doBy::summaryBy(jday~cell6+year,data=mcd3, na.rm=T)
save(hex6_mcd12Q2, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex6_mcd12Q2.Rdata")




### Now for the max

mcd12Q2_bounded <- max_data[which(!is.na(max_data$layer)), ]
mcd_sp <- SpatialPointsDataFrame(mcd12Q2_bounded[,c(1,2)], mcd12Q2_bounded[,c(3,4)], proj4string = CRS(proj4string(raster1)))
mcd_sp2 <- spTransform(mcd_sp, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
mcd3 <- as.data.frame(mcd_sp2)
mcd3$jday <- mcd3$layer - jdayMap[mcd3$year - 2000]


# Associate each pixel with the cell of a hexagonal grid
hexgrid7 <- dgconstruct(res=7) # Construct geospatial hexagonal grid
mcd3$cell7 <- dgGEO_to_SEQNUM(hexgrid7, mcd3$x, mcd3$y)$seqnum
names(mcd3) <- c("raw_day","year","x","y","jday","cell")

# Get the mean onset greenness for each cell-year
hex7_mcd12Q2_max <- doBy::summaryBy(jday~cell+year,data=mcd3, na.rm=T)
save(hex7_mcd12Q2_max, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex7_mcd12Q2_max.Rdata")

# Associate each pixel with the cell of a hexagonal grid
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
mcd3$cell6 <- dgGEO_to_SEQNUM(hexgrid6, mcd3$x, mcd3$y)$seqnum

# Get the mean onset greenness for each cell-year
hex6_mcd12Q2_max <- doBy::summaryBy(jday~cell6+year,data=mcd3, na.rm=T)
save(hex6_mcd12Q2_max, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex6_mcd12Q2_max.Rdata")
