library(gdalUtils)
library(raster)
library(dggridR)
library(ggplot2)
library(viridis)
setwd('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/eMODIS/')
years <- c(2001:2016)

# Import GEOTIFFs as rasters
tif_files <- dir(pattern = "_TIF")
zips <- dir(pattern = ".ZIP")
tif_folders <- tif_files[-which(tif_files %in% zips)]
tWs <- tif_folders[grep('W', tif_folders)]
tEs <- tif_folders[-which(tif_folders %in% tWs)]
mcd <- as.list(rep(NA, length(years)))
names(mcd) <- paste0("year", as.character(years))
for(i in 1:16){
  print(i)
  if(i < 15.5){
    mcdE <- raster(paste0('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/eMODIS/',
                          tEs[i],'/SOST',years[i],'_eUSTeM250m_v1.tif'))
    mcdW <- raster(paste0('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/eMODIS/',
                          tWs[i],'/SOST',years[i],'_wUSTeM250m_v1.tif'))
  }else{
    mcdE <- raster(paste0('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/eMODIS/',
                          tEs[i],'/SOST',years[i],'_eUSAeM250m_v2.tif'))
    mcdW <- raster(paste0('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/eMODIS/',
                          tWs[i],'/SOST',years[i],'_wUSAeM250m_v2.tif'))
  }
  mcd[[i]] <- raster::merge(mcdE, mcdW, tolerance=0)
}

# Reduce the resolution of the raster; otherwise we get > 200 million cells per year, and errors like
# the following, which I had never seen before:
#   Error in which(is.na(eMODIS_data)) : 
#      long vectors not supported yet: ../../../../R-3.4.4/src/include/Rinlinedfuns.h:138
mcd2 <- mcd
for(i in 1:length(years)){
  print(i)
  mcd2[[i]][mcd2[[i]] == 1000] <- NA
  mcd2[[i]][mcd2[[i]] < 1] <- NA
}

for(i in 1:length(years)){
  print(i)
  assign(paste("onset.",years[i],sep=""), raster::aggregate(mcd2[[i]], fact=10, na.rm = T))
}

# # Make sure the rasters all have the same grid
# test1 <- as.data.frame(onset.2001, xy=T)
# test2 <- as.data.frame(onset.2014, xy=T)
# all.equal(test1$x, test2$x)
# # TRUE
# all.equal(test1$y, test2$y)
# # TRUE

# Convert rasters to a single dataframe
eMODIS_data <- as.data.frame(onset.2001, xy=T)
for(i in 2:length(mcd)){
  eMODIS_data <- cbind(eMODIS_data, as.data.frame(get(paste('onset.', years[i], sep="")), xy=F))
}
names(eMODIS_data)[3:18] <- paste0("onset",years)


# Get rid of rows with no data and crop to North America
nodatarows <- rowSums(!is.na(eMODIS_data[,3:18]))
eMODIS_reduced <- eMODIS_data[which(nodatarows!=0),]
reproj_frame <- SpatialPointsDataFrame(eMODIS_reduced[,1:2], eMODIS_reduced[,3:18], 
                                       proj4string = CRS(proj4string(onset.2001)))
reproj <- spTransform(reproj_frame, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

mcd_reproj <- as.data.frame(reproj)

eMODIS_bounded <- mcd_reproj[which(mcd_reproj$x < -50 & mcd_reproj$y > 14),]

# Associate each pixel with the cell of a hexagonal grid
hexgrid7 <- dgconstruct(res=7) # Construct geospatial hexagonal grid
eMODIS_bounded$cell7 <- dgGEO_to_SEQNUM(hexgrid7, eMODIS_bounded$x, eMODIS_bounded$y)$seqnum
eMODIS_cells <- unique(eMODIS_bounded$cell7)
eMODIS_long <- reshape2::melt(eMODIS_bounded, id=c("cell7","x","y"), measure=names(eMODIS_bounded)[1:16])
names(eMODIS_long) <- c("cell","x","y","onsetyear","onset")
eMODIS_long$year <- as.numeric(stringr::str_extract_all(eMODIS_long$onsetyear, "[0-9]+"))

# Get the mean onset greenness for each cell-year
hex7_eMODIS <- doBy::summaryBy(onset~cell+onsetyear,data=eMODIS_long, na.rm=T)
hex7_eMODIS$year <- as.numeric(stringr::str_extract_all(hex7_eMODIS$onsetyear, "[0-9]+"))
hex7_eMODIS <- hex7_eMODIS[,c(1,4,3)]

save(hex7_eMODIS, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex7_eMODIS.Rdata")


# Associate each pixel with the cell of a hexagonal grid
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
eMODIS_bounded$cell6 <- dgGEO_to_SEQNUM(hexgrid6, eMODIS_bounded$x, eMODIS_bounded$y)$seqnum
eMODIS_cells <- unique(eMODIS_bounded$cell6)
eMODIS_long <- reshape2::melt(eMODIS_bounded, id=c("cell6","x","y"), measure=names(eMODIS_bounded)[1:16])
names(eMODIS_long) <- c("cell","x","y","onsetyear","onset")
eMODIS_long$year <- as.numeric(stringr::str_extract_all(eMODIS_long$onsetyear, "[0-9]+"))

# Get the mean onset greenness for each cell-year
hex6_eMODIS <- doBy::summaryBy(onset~cell+onsetyear,data=eMODIS_long, na.rm=T)
hex6_eMODIS$year <- as.numeric(stringr::str_extract_all(hex6_eMODIS$onsetyear, "[0-9]+"))
hex6_eMODIS <- hex6_eMODIS[,c(1,4,3)]

save(hex6_eMODIS, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex6_eMODIS.Rdata")


