library(gdalUtils)
library(raster)
library(dggridR)
library(ggplot2)
library(viridis)
setwd('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/NAPHEN_SI/')
years <- c(1981:2014)

# Import GEOTIFFs as rasters
tiffs <- dir(pattern = ".tif")
for(i in 1:length(years)){
  assign(paste("onset.",years[i],sep=""), raster(tiffs[i]))
}

# Make sure the rasters all have the same grid
# test1 <- as.data.frame(onset.1981, xy=T)
# test2 <- as.data.frame(onset.2014, xy=T)
# all.equal(test1$x, test2$x)
# # TRUE
# all.equal(test1$y, test2$y)
# # TRUE


# Convert rasters to a single dataframe
naphen_si_data <- as.data.frame(onset.1981, xy=T)
for(i in 2:length(tiffs)){
  naphen_si_data <- cbind(naphen_si_data, as.data.frame(get(paste('onset.', years[i], sep="")), xy=F))
}

# Get rid of rows with no data and crop to North America
nodatarows <- rowSums(!is.na(naphen_si_data[,3:36]))
naphen_si_reduced <- naphen_si_data[which(nodatarows!=0),]
naphen_si_bounded <- naphen_si_reduced[which(naphen_si_reduced$x < -50 & naphen_si_reduced$y > 14),]
names(naphen_si_bounded)[3:36] <- paste0('onset', years)

# Associate each pixel with the cell of a hexagonal grid
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
naphen_si_bounded$cell6 <- dgGEO_to_SEQNUM(hexgrid6, naphen_si_bounded$x, naphen_si_bounded$y)$seqnum
naphen_si_cells <- unique(naphen_si_bounded$cell6)
naphen_si_long <- reshape2::melt(naphen_si_bounded, id=c("cell6","x","y"), measure=names(naphen_si_bounded)[3:36])
names(naphen_si_long) <- c("cell","x","y","onsetyear","onset")
naphen_si_long$year <- as.numeric(stringr::str_extract_all(naphen_si_long$onsetyear, "[0-9]+"))

# Get the mean onset greenness for each cell-year
hex6_naphen_si <- doBy::summaryBy(onset~cell+onsetyear,data=naphen_si_long, na.rm=T)
hex6_naphen_si$year <- as.numeric(stringr::str_extract_all(hex6_naphen_si$onsetyear, "[0-9]+"))
hex6_naphen_si <- hex6_naphen_si[,c(1,4,3)]

save(hex6_naphen_si, file = '/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex6_naphen_si.Rdata')

# Choose time interval under study
hex6_naphen_si2 <- hex6_naphen_si[which(hex6_naphen_si$year>1900 & hex6_naphen_si$year<2020),]

# Get linear trend in onset greenness for each cell and write to a dataframe
naphen_sitrend <- naphen_si_p <- rep(NA, length(unique(hex6_naphen_si2$cell)))
for(i in 1:length(unique(hex6_naphen_si2$cell))){
  celldat <- hex6_naphen_si2[which(hex6_naphen_si2$cell == unique(hex6_naphen_si2$cell)[i]),]
  if(sum(!is.na(celldat$onset.mean)) > 5){
    fit <- lm(celldat$onset.mean ~ celldat$year)
    naphen_sitrend[i] <- fit$coefficients[2]
    naphen_si_p[i] <- summary(fit)$coefficients[2,4]
  }
}
mt_frame <- data.frame(cell=unique(hex6_naphen_si2$cell), trend=naphen_sitrend, p_value=naphen_si_p)

#mt_frame <- mt_frame[which(mt_frame$p_value < .05),]

# Plot the which trends are positive (spring gets later) versus negative (spring gets earlier)
grid <- dgcellstogrid(hexgrid6, mt_frame$cell, frame=TRUE,wrapcells=TRUE)
grid <- merge(grid, mt_frame, by = "cell")
grid$advancement <- as.numeric(grid$trend < 0)
grid$trend2 <- atan(atan(grid$trend))

p<- ggplot() + 
  #  geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_polygon(data=grid,      aes(x=long, y=lat, group=group, fill=trend), alpha=0.9)    +
  geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
  #  geom_point  (aes(x=cellcenters$lon_deg, y=cellcenters$lat_deg)) +
  scale_fill_gradientn(colors=c('red','yellow','blue'))
p+coord_map("ortho", orientation = c(35, -95, 0))+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())
