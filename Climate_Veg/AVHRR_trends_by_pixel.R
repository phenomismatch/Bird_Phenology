# This is old code and maybe not useful for much; I'm keeping it around in case some of it becomes
# useful in the future. Most of the veg phenology summaries that I've produced are based on aggragating
# data into averages at the scale of a grid scale, and then (for example) computing linear trends
# at the level of the cell.  This script, written for the AVHRR dataset, computes trends by pixel, and
# then aggregates the estimated trends to the level of the grid cell.

library(gdalUtils)
#library(rgdal)
library(raster)
library(dggridR)
library(ggplot2)
library(viridis)
'%ni%' <- Negate('%in%')
setwd('/Users/Tingleylab/Desktop/useful_datasets/AVHRR_phenology/')

# Import GEOTIFFS as rasters
zips <- dir(pattern = ".ZIP")
all.folders <- list.files()[list.files() %ni% zips]
years <- c(1989:2013)

for(i in 1:length(years)){
  assign(paste0("SOS.",years[i]), raster(paste0('/Users/Tingleylab/Desktop/useful_datasets/AVHRR_phenology/',
                                                all.folders[i], '/av_SOST', years[i], 'v4.tif')))
}

# Convert rasters to a single dataframe
avhrr_data <- as.data.frame(SOS.1989, xy=T)
for(i in 2:length(years)){
  avhrr_data <- cbind(avhrr_data, as.data.frame(get(paste('SOS.', years[i], sep="")), xy=F))
}



# Get rid of rows with no data and crop to North America
nodatarows <- which(avhrr_data[,3] < 1 | avhrr_data[,3] > 366 | 
                      avhrr_data[,4] < 1 | avhrr_data[,4] > 366 | 
                      avhrr_data[,5] < 1 | avhrr_data[,5] > 366 | 
                      avhrr_data[,6] < 1 | avhrr_data[,6] > 366 | 
                      avhrr_data[,7] < 1 | avhrr_data[,7] > 366 | 
                      avhrr_data[,8] < 1 | avhrr_data[,8] > 366 | 
                      avhrr_data[,9] < 1 | avhrr_data[,9] > 366 | 
                      avhrr_data[,10] < 1 | avhrr_data[,10] > 366 | 
                      avhrr_data[,11] < 1 | avhrr_data[,11] > 366 | 
                      avhrr_data[,12] < 1 | avhrr_data[,12] > 366 | 
                      avhrr_data[,13] < 1 | avhrr_data[,13] > 366 | 
                      avhrr_data[,14] < 1 | avhrr_data[,14] > 366 | 
                      avhrr_data[,15] < 1 | avhrr_data[,15] > 366 | 
                      avhrr_data[,16] < 1 | avhrr_data[,16] > 366 | 
                      avhrr_data[,17] < 1 | avhrr_data[,17] > 366 | 
                      avhrr_data[,18] < 1 | avhrr_data[,18] > 366 | 
                      avhrr_data[,19] < 1 | avhrr_data[,19] > 366 | 
                      avhrr_data[,20] < 1 | avhrr_data[,20] > 366 | 
                      avhrr_data[,21] < 1 | avhrr_data[,21] > 366 | 
                      avhrr_data[,22] < 1 | avhrr_data[,22] > 366 | 
                      avhrr_data[,23] < 1 | avhrr_data[,23] > 366 | 
                      avhrr_data[,24] < 1 | avhrr_data[,24] > 366 | 
                      avhrr_data[,25] < 1 | avhrr_data[,25] > 366 | 
                      avhrr_data[,26] < 1 | avhrr_data[,26] > 366 | 
                      avhrr_data[,27] < 1 | avhrr_data[,27] > 366)

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
avhrr_long$pixel <- paste(avhrr_long$x, avhrr_long$y, sep="_")

upix <- unique(avhrr_long$pixel)
upix2 <- sample(upix,100000)
avhrr_long2 <- avhrr_long[which(avhrr_long$pixel %in% upix2),]
avhrrtrend <- avhrr_p <- rep(NA, length(upix2))
for(i in 1:length(upix2)){
  print(i)
  pixeldat <- avhrr_long2[which(avhrr_long2$pixel == upix2[i]),]
  fit <- lm(pixeldat$SOS ~ pixeldat$year)
  avhrrtrend[i] <- fit$coefficients[2]
  avhrr_p[i] <- summary(fit)$coefficients[2,4]
}


at_frame <- data.frame(pixel=upix2, trend=avhrrtrend, p_value=avhrr_p, cell=NA)
for(i in 1:dim(at_frame)[1]){
  at_frame$cell[i] <- unique(avhrr_long2$cell[which(avhrr_long2$pixel==at_frame$pixel[i])])
}

#at_frame <- at_frame[which(at_frame$p_value < .05),]

# Get the mean onset greenness for each cell-year
hex7_avhrr <- doBy::summaryBy(trend~cell,data=at_frame, na.rm=T)

# Plot the which trends are positive (spring gets later) versus negative (spring gets earlier)
grid <- dgcellstogrid(hexgrid7, hex7_avhrr$cell, frame=TRUE,wrapcells=TRUE)
grid <- merge(grid, at_frame, by = "cell")
grid$advancement <- as.numeric(grid$trend < 0)
grid$trend2 <- atan(atan(grid$trend))

p<- ggplot() + 
  #  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=trend2), alpha=0.9)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #  geom_point  (aes(x=cellcenters$lon_deg, y=cellcenters$lat_deg)) +
  scale_fill_gradientn(colors=c('red','yellow','blue'))
p+coord_map("ortho", orientation = c(35, -95, 0))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())

