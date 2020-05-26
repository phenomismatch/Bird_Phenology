library(gdalUtils)
library(raster)
library(dggridR)
library(ggplot2)
library(viridis)
'%ni%' <- Negate('%in%')
setwd('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/USGS_EROS/')

# Import GEOTIFFS as rasters
zips <- dir(pattern = ".ZIP")
all.folders <- list.files()[list.files() %ni% zips]
years <- c(1989:2013)

for(i in 1:length(years)){
  assign(paste0("SOS.",years[i]), raster(paste0('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/USGS_EROS/',
                                                all.folders[i], '/av_SOST', years[i], 'v4.tif')))
}

# Convert rasters to a single dataframe
avhrr_data <- as.data.frame(SOS.1989, xy=T)
for(i in 2:length(years)){
  avhrr_data <- cbind(avhrr_data, as.data.frame(get(paste('SOS.', years[i], sep="")), xy=F))
}



# Get rid of rows with no data and crop to North America
nodatarows <- which(avhrr_data[,3] < -150 | avhrr_data[,3] > 366 | 
                      avhrr_data[,4] < -150 | avhrr_data[,4] > 366 | 
                      avhrr_data[,5] < -150 | avhrr_data[,5] > 366 | 
                      avhrr_data[,6] < -150 | avhrr_data[,6] > 366 | 
                      avhrr_data[,7] < -150 | avhrr_data[,7] > 366 | 
                      avhrr_data[,8] < -150 | avhrr_data[,8] > 366 | 
                      avhrr_data[,9] < -150 | avhrr_data[,9] > 366 | 
                      avhrr_data[,10] < -150 | avhrr_data[,10] > 366 | 
                      avhrr_data[,11] < -150 | avhrr_data[,11] > 366 | 
                      avhrr_data[,12] < -150 | avhrr_data[,12] > 366 | 
                      avhrr_data[,13] < -150 | avhrr_data[,13] > 366 | 
                      avhrr_data[,14] < -150 | avhrr_data[,14] > 366 | 
                      avhrr_data[,15] < -150 | avhrr_data[,15] > 366 | 
                      avhrr_data[,16] < -150 | avhrr_data[,16] > 366 | 
                      avhrr_data[,17] < -150 | avhrr_data[,17] > 366 | 
                      avhrr_data[,18] < -150 | avhrr_data[,18] > 366 | 
                      avhrr_data[,19] < -150 | avhrr_data[,19] > 366 | 
                      avhrr_data[,20] < -150 | avhrr_data[,20] > 366 | 
                      avhrr_data[,21] < -150 | avhrr_data[,21] > 366 | 
                      avhrr_data[,22] < -150 | avhrr_data[,22] > 366 | 
                      avhrr_data[,23] < -150 | avhrr_data[,23] > 366 | 
                      avhrr_data[,24] < -150 | avhrr_data[,24] > 366 | 
                      avhrr_data[,25] < -150 | avhrr_data[,25] > 366 | 
                      avhrr_data[,26] < -150 | avhrr_data[,26] > 366 | 
                      avhrr_data[,27] < -150 | avhrr_data[,27] > 366)

avhrr_reduced <- avhrr_data[-nodatarows,]
avhrr_reduced_sp <- SpatialPointsDataFrame(avhrr_reduced[,c(1,2)], avhrr_reduced[,c(3:27)], proj4string = CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
avhrr_reduced2 <- spTransform(avhrr_reduced_sp, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
avhrr_reduced3 <- as.data.frame(avhrr_reduced2)
avhrr_bounded <- avhrr_reduced3[which(avhrr_reduced3$x < -50 & avhrr_reduced3$y > 14),]

# Associate each pixel with the cell of a hexagonal grid
hexgrid7 <- dgconstruct(res=7) # Construct geospatial hexagonal grid
avhrr_bounded$cell7 <- dgGEO_to_SEQNUM(hexgrid7, avhrr_bounded$x, avhrr_bounded$y)$seqnum
avhrr_cells <- unique(avhrr_bounded$cell7)
avhrr_long <- reshape2::melt(avhrr_bounded, id=c("cell7","x","y"), measure=names(avhrr_bounded)[1:25])
names(avhrr_long) <- c("cell","x","y","SOSyear","SOS")
avhrr_long$year <- as.numeric(stringr::str_extract(avhrr_long$SOSyear, "[0-9]+"))

# Get the mean SOS greenness for each cell-year
hex7_avhrr <- doBy::summaryBy(SOS~cell+SOSyear,data=avhrr_long, na.rm=T)
hex7_avhrr$year <- as.numeric(stringr::str_extract(hex7_avhrr$SOSyear, "[0-9]+"))
hex7_avhrr <- hex7_avhrr[,c(1,4,3)]

save(hex7_avhrr, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex7_avhrr.Rdata")



hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
avhrr_bounded$cell6 <- dgGEO_to_SEQNUM(hexgrid6, avhrr_bounded$x, avhrr_bounded$y)$seqnum
avhrr_cells <- unique(avhrr_bounded$cell6)
avhrr_long <- reshape2::melt(avhrr_bounded, id=c("cell6","x","y"), measure=names(avhrr_bounded)[1:25])
names(avhrr_long) <- c("cell","x","y","SOSyear","SOS")
avhrr_long$year <- as.numeric(stringr::str_extract(avhrr_long$SOSyear, "[0-9]+"))

# Get the mean SOS greenness for each cell-year
hex6_avhrr <- doBy::summaryBy(SOS~cell+SOSyear,data=avhrr_long, na.rm=T)
hex6_avhrr$year <- as.numeric(stringr::str_extract(hex6_avhrr$SOSyear, "[0-9]+"))
hex6_avhrr <- hex6_avhrr[,c(1,4,3)]

save(hex6_avhrr, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex6_avhrr.Rdata")
