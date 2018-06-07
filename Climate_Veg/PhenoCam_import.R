library(stringr)
library(lubridate)
library(dggridR)

hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid
hexgrid7 <- dgconstruct(res=7)

setwd("/Users/TingleyLab/Desktop/useful_datasets/Veg_phenology/PhenoCam/PhenoCam_V1_1511/data/")

files <- list.files()
sites <- unique(str_extract(files, "[^_]+"))
nsites <- length(sites)
coord.files <- files[grep("_1day.csv",files)]

# Extract longitude and latitude from the relevant files, and save to data.frame
site_locs <- data.frame(site = sites, lat = NA, lon = NA)
for(i in 1:nsites){
  myfile <- read.csv(coord.files[grep(sites[i], coord.files)][1], nrows = 2, skip = 5)
  if(length(grep("-", myfile[1,1]))==0){
    site_locs$lat[i] <- as.numeric(str_extract(myfile[1,1], "\\d+\\.*\\d*"))
  }else{
    site_locs$lat[i] <- as.numeric(paste0("-", str_extract(myfile[1,1], "\\d+\\.*\\d*")))
  }
  
  if(length(grep("-", myfile[2,1]))==0){
    site_locs$lon[i] <- as.numeric(str_extract(myfile[2,1], "\\d+\\.*\\d*"))
  }else{
    site_locs$lon[i] <- as.numeric(paste0("-", str_extract(myfile[2,1], "\\d+\\.*\\d*")))
  }
}

# Get the phenophase dates for each ROI at each site, and save to dataframe
data.files <- files[grep("1day_transition_dates.csv",files)]
nseries <- length(data.files)

phenocam_data <- data.frame()
for(i in 1:nseries){
  myfile <- read.csv(data.files[i], skip = 16)
  myfile$lat <- site_locs$lat[which(as.character(site_locs$site) == as.character(myfile$sitename[1]))]
  myfile$lon <- site_locs$lon[which(as.character(site_locs$site) == as.character(myfile$sitename[1]))]
  phenocam_data <- rbind(phenocam_data, myfile)
}
na_phenocam <- phenocam_data[which(phenocam_data$lat > 24 & phenocam_data$lon < -50 & phenocam_data$direction == "rising"), ]

# Put in columns for year, unique site/vegtype/ROI, and hexagonal grid cells, and convert phenophase dates 
# to julian day
na_phenocam$year <- year(as.Date(na_phenocam$transition_10))
na_phenocam$uniqueSVR <- paste(na_phenocam$sitename, na_phenocam$veg_type, na_phenocam$roi_id, sep="_")
na_phenocam$cell6 <- dgGEO_to_SEQNUM(hexgrid6, na_phenocam$lon, na_phenocam$lat)$seqnum
na_phenocam$cell7 <- dgGEO_to_SEQNUM(hexgrid7, na_phenocam$lon, na_phenocam$lat)$seqnum
for(i in 6:14){
  na_phenocam[,i] <- yday(as.Date(na_phenocam[,i]))
}

save(na_phenocam, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/Veg_phenology/na_phenocam.Rdata")
