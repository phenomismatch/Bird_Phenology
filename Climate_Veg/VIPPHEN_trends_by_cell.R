library(gdalUtils)
library(raster)
library(dggridR)
library(ggplot2)
library(viridis)
setwd('/Users/Tingleylab/Desktop/useful_datasets/Veg_phenology/VIPPHEN_EVI2/')
years <- c(1981:2014)

# # Get the start of season 1 layer and save as GEOTIFF
# hdfs <- dir(pattern = ".hdf")
# years <- c(1981:2014)
# for(i in 1:length(years)){
#   sds <- get_subdatasets(hdfs[grep(paste(years[i],".004",sep=""),hdfs)])
#   gdal_translate(sds[26], dst_dataset = paste('reliability', as.character(years[i]), '.tif', sep=""))
# }
# 
# for(i in 1:length(years)){
#   sds <- get_subdatasets(hdfs[grep(paste(years[i],".004",sep=""),hdfs)])
#   gdal_translate(sds[1], dst_dataset = paste('onset', as.character(years[i]), '.tif', sep=""))
# }

# Import GEOTIFFs as rasters
tiffs <- dir(pattern = ".tif")
xmls <- dir(pattern = ".xml")
tiffs <- tiffs[-which(tiffs %in% xmls)]
tiffs.onset <- tiffs[grep('onset', tiffs)]
tiffs.reliability <- tiffs[grep('reliability', tiffs)]
for(i in 1:length(years)){
  assign(paste("onset.",years[i],sep=""), raster(tiffs.onset[i]))
  assign(paste("reliability.",years[i],sep=""), raster(tiffs.reliability[i]))
}

## Make sure the rasters all have the same grid
#test1 <- as.data.frame(onset.1981, xy=T)
#test2 <- as.data.frame(onset.2014, xy=T)
#all.equal(test1$x, test2$x)
## TRUE
#all.equal(test1$y, test2$y)
## TRUE
#test3 <- as.data.frame(reliability.1981, xy=T)
#all.equal(test1$x, test3$x)
## TRUE
#all.equal(test1$y, test3$y)
## TRUE

# Convert rasters to a single dataframe
modis_data <- as.data.frame(onset.1981, xy=T)
for(i in 2:length(tiffs.onset)){
  modis_data <- cbind(modis_data, as.data.frame(get(paste('onset.', years[i], sep="")), xy=F))
}

modis_reliability <- as.data.frame(reliability.1981, xy=T)
for(i in 2:length(tiffs.reliability)){
  modis_reliability <- cbind(modis_reliability, as.data.frame(get(paste('reliability.', years[i], sep="")), xy=F))
}


# Get rid of rows with no data and crop to North America
nodatarows <- rowSums(!is.na(modis_data[,3:36]))
modis_reduced <- modis_data[which(nodatarows!=0),]
modis_bounded <- modis_reduced[which(modis_reduced$x < -50 & modis_reduced$y > 14),]

modis.r_reduced <- modis_reliability[which(nodatarows!=0),]
modis.r_bounded <- modis.r_reduced[which(modis_reduced$x < -50 & modis_reduced$y > 14),]

# Get rid of pixels that are unreliable in any year.
attach(modis.r_bounded)
endcutoff <- 2
begincutoff <- 0

for(j in 3:36){
 cj <- modis_bounded[,j]
 cjr <- modis.r_bounded[,j]
 cj[which(cjr < begincutoff | cjr>endcutoff)] <- NA
 modis_bounded[,j] <- cj
}

# Get rid of pixels with values of -1 (mask) or greater than 365 (no idea what these are, but there aren't too many)
attach(modis_bounded)
endcutoff <- 365
begincutoff <- 1
modis_bounded <- modis_bounded[-which(onset1981 < begincutoff | onset1981 > endcutoff |
                                       onset1982 < begincutoff | onset1982 > endcutoff |
                                       onset1983 < begincutoff | onset1983 > endcutoff |
                                       onset1984 < begincutoff | onset1984 > endcutoff |
                                       onset1985 < begincutoff | onset1985 > endcutoff |
                                       onset1986 < begincutoff | onset1986 > endcutoff |
                                       onset1987 < begincutoff | onset1987 > endcutoff |
                                       onset1988 < begincutoff | onset1988 > endcutoff |
                                       onset1989 < begincutoff | onset1989 > endcutoff |
                                       onset1990 < begincutoff | onset1990 > endcutoff |
                                       onset1991 < begincutoff | onset1991 > endcutoff |
                                       onset1992 < begincutoff | onset1992 > endcutoff |
                                       onset1993 < begincutoff | onset1993 > endcutoff |
                                       onset1994 < begincutoff | onset1994 > endcutoff |
                                       onset1995 < begincutoff | onset1995 > endcutoff |
                                       onset1996 < begincutoff | onset1996 > endcutoff |
                                       onset1997 < begincutoff | onset1997 > endcutoff |
                                       onset1998 < begincutoff | onset1998 > endcutoff |
                                       onset1999 < begincutoff | onset1999 > endcutoff |
                                       onset2000 < begincutoff | onset2000 > endcutoff |
                                       onset2001 < begincutoff | onset2001 > endcutoff |
                                       onset2002 < begincutoff | onset2002 > endcutoff |
                                       onset2003 < begincutoff | onset2003 > endcutoff |
                                       onset2004 < begincutoff | onset2004 > endcutoff |
                                       onset2005 < begincutoff | onset2005 > endcutoff | 
                                       onset2006 < begincutoff | onset2006 > endcutoff | 
                                       onset2007 < begincutoff | onset2007 > endcutoff | 
                                       onset2008 < begincutoff | onset2008 > endcutoff | 
                                       onset2009 < begincutoff | onset2009 > endcutoff | 
                                       onset2010 < begincutoff | onset2010 > endcutoff | 
                                       onset2011 < begincutoff | onset2011 > endcutoff | 
                                       onset2012 < begincutoff | onset2012 > endcutoff | 
                                       onset2013 < begincutoff | onset2013 > endcutoff | 
                                       onset2014 < begincutoff | onset2014 > endcutoff),]

# Associate each pixel with the cell of a hexagonal grid
hexgrid7 <- dgconstruct(res=7) # Construct geospatial hexagonal grid
modis_bounded$cell7 <- dgGEO_to_SEQNUM(hexgrid7, modis_bounded$x, modis_bounded$y)$seqnum
modis_cells <- unique(modis_bounded$cell7)
modis_long <- reshape2::melt(modis_bounded, id=c("cell7","x","y"), measure=names(modis_bounded)[3:36])
names(modis_long) <- c("cell","x","y","onsetyear","onset")
modis_long$year <- as.numeric(stringr::str_extract_all(modis_long$onsetyear, "[0-9]+"))

# Get the mean onset greenness for each cell-year
hex7_modis <- doBy::summaryBy(onset~cell+onsetyear,data=modis_long, na.rm=T)
hex7_modis$year <- as.numeric(stringr::str_extract_all(hex7_modis$onsetyear, "[0-9]+"))
hex7_VIPPHEN <- hex7_modis[,c(1,4,3)]

save(hex7_VIPPHEN, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex7_VIPPHEN.Rdata")



# Associate each pixel with the cell of a hexagonal grid
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
modis_bounded$cell6 <- dgGEO_to_SEQNUM(hexgrid6, modis_bounded$x, modis_bounded$y)$seqnum
modis_cells <- unique(modis_bounded$cell6)
modis_long <- reshape2::melt(modis_bounded, id=c("cell6","x","y"), measure=names(modis_bounded)[3:36])
names(modis_long) <- c("cell","x","y","onsetyear","onset")
modis_long$year <- as.numeric(stringr::str_extract_all(modis_long$onsetyear, "[0-9]+"))

# Get the mean onset greenness for each cell-year
hex6_modis <- doBy::summaryBy(onset~cell+onsetyear,data=modis_long, na.rm=T)
hex6_modis$year <- as.numeric(stringr::str_extract_all(hex6_modis$onsetyear, "[0-9]+"))
hex6_VIPPHEN <- hex6_modis[,c(1,4,3)]

save(hex6_VIPPHEN, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/Hex_gridded_products/hex6_VIPPHEN.Rdata")
