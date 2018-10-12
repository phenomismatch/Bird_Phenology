######################
# 1 - import ebird data
#
# Import data from eBird and filter
# Script will largely be replaced by DB query
# Separate data file for each species to parallelize more easily
# Formerly Ebird_import.R
######################


# This script imports data from the eBird reference dataset (ERD) for a vector of species, defined by 
# species_list and a vector of years defined by years.  The script is designed to interact with the ERD
# in the format downloaded directly from eBird, with a separate folder for each year from 2002 onwards.

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# runtime -----------------------------------------------------------------

tt <- proc.time()


# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# load packages -----------------------------------------------------------

library(dggridR)
library(dplyr)


# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)

years <- 2002:2016
nyr <- length(years)




# Filter eBird data using auk ---------------------------------------------------

# # run awk script directly from R to get the column indicies of the species of interest, make sure
# # the indicies are the same across all year files:
# species_indices <- as.data.frame(matrix(data=NA,nrow=length(species_list), ncol=length(years)))
# for(j in 1:length(years)){
#   setwd(paste0("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/",years[j],"/"))
#   for(i in 1:length(species_list)){
#     unix_command <- paste0("head -1 checklists.csv | awk -v RS=\",\" '/",species_list[i],
#                            "/{print NR;}'")
#     species_indices[i,j] <- system(unix_command, intern = T)
#   }
# }
# equalTest <- vector()
# for(j in 2:length(years)){
#   equalTest[j-1] <- all.equal(species_indices[,1], species_indices[,j])
# }
# sum(equalTest)
# spInd <- as.vector(species_indices[,1])
# 
# # For each year, run awk script to extract just the columns of interest
# unix_command <- paste0("cut -d \",\" -f1-8,12-19,",paste(spInd,collapse=","), " checklists.csv > checklists_NA_birdPhen.csv")
# for(j in 1:length(years)){
#   setwd(paste0("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/",years[j],"/"))
#   system(unix_command)
# }
# 
# # For each year, run awk script to extract only observations that occurred on or before jday 200:
# setwd("/Users/TingleyLab/Desktop/useful_datasets/eBird/ERD2016SS/")
# for(j in 1:length(years)){
#   unix_command <- paste0("awk -f extract_awk_time.txt ",years[j],"/checklists_NA_birdPhen.csv > ",
#                          years[j],"/checklists_NA_birdPhen_time.csv")
#   system(unix_command)
# }



# Import data to list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Raw'))

# For each year, import the data, and store in a list:
dataList <- as.list(rep(NA, nyr))

for(i in 1:nyr)
{
  assign(paste0("data", years[i]), 
         read.csv(paste0("eBird/ERD2016SS/", years[i], "/checklists_NA_birdPhen_time.csv"),
                  quote = "", stringsAsFactors = FALSE))
  
  dataList[[i]] <- eval(parse(text = paste0("data", years[i])))
  rm(list = paste("data", years[i], sep=""))
  print(years[i])
}



#save to rds object
ebird_NA_phen <- data.table::rbindlist(dataList)
rm(dataList)

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

saveRDS(ebird_NA_phen, file = 'ebird_NA_phen.rds')




# Filter data using R ------------------------------------------------------------

#ebird_NA_phen <- readRDS('ebird_NA_phen.rds')

ebird_NA_phen_proc <- ebird_NA_phen[which(ebird_NA_phen$PRIMARY_CHECKLIST_FLAG==1),]
rm(ebird_NA_phen)

ebird_NA_phen_proc$EFFORT_HRS <- as.numeric(as.character(ebird_NA_phen_proc$EFFORT_HRS))
ebird_NA_phen_proc <- ebird_NA_phen_proc[which(ebird_NA_phen_proc$EFFORT_HRS > 
                                                 0.1 & ebird_NA_phen_proc$EFFORT_HRS<24), ]
ebird_NA_phen_proc <- ebird_NA_phen_proc[which(ebird_NA_phen_proc$YEAR > 2001), ]
ebird_NA_phen_proc <- ebird_NA_phen_proc[-which(ebird_NA_phen_proc$EFFORT_DISTANCE_KM == "?"),]
ebird_NA_phen_proc <- ebird_NA_phen_proc[-which((ebird_NA_phen_proc$TIME + 
                                                   ebird_NA_phen_proc$EFFORT_HRS) < 6), ]
ebird_NA_phen_proc <- ebird_NA_phen_proc[-which(ebird_NA_phen_proc$TIME > 16), ]
ebird_NA_phen_proc <- ebird_NA_phen_proc[which(ebird_NA_phen_proc$LONGITUDE > -100 & 
                                                 ebird_NA_phen_proc$LONGITUDE < -50 & 
                                                 ebird_NA_phen_proc$LATITUDE > 26), ]
ebird_NA_phen_proc$EFFORT_DISTANCE_KM <- as.numeric(as.character(ebird_NA_phen_proc$EFFORT_DISTANCE_KM))
ebird_NA_phen_proc <- ebird_NA_phen_proc[which(ebird_NA_phen_proc$EFFORT_DISTANCE_KM >= 0 & 
                                                 ebird_NA_phen_proc$EFFORT_DISTANCE_KM < 100),]


#convert to data.frame bc data.table is annoying
ebird_NA_phen_proc2 <- as.data.frame(ebird_NA_phen_proc)

#convert counts to presence/absence (1/0)
for(i in 17:130)
{
  print(i)
  ttt <- ebird_NA_phen_proc2[, i]
  ttt[which(ttt == "X")] <- 1
  ttt <- as.numeric(ttt)
  ttt[ttt > 0] <- 1
  if(min(ttt) < 0){stop()}
  ebird_NA_phen_proc[, i] <- ttt
}

rm(ebird_NA_phen_proc2)

#add jday squared and cubed and effort hours
ebird_NA_phen_proc$sjday <- as.vector(scale(as.numeric(as.character(ebird_NA_phen_proc$DAY))))
ebird_NA_phen_proc$sjday2 <- ebird_NA_phen_proc$sjday^2
ebird_NA_phen_proc$sjday3 <- ebird_NA_phen_proc$sjday^3
ebird_NA_phen_proc$shr <- as.vector(scale(ebird_NA_phen_proc$EFFORT_HRS))



# bin to hex grid ---------------------------------------------------------

# Construct geospatial hexagonal grid

hexgrid6 <- dggridR::dgconstruct(res=6) 
ebird_NA_phen_proc$cell6 <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                                     ebird_NA_phen_proc$LONGITUDE, 
                                                     ebird_NA_phen_proc$LATITUDE)[[1]]

#save to rds object
saveRDS(ebird_NA_phen_proc, file = 'ebird_NA_phen_proc.rds')




# Different rds object each species ---------------------------------------

#ebird_NA_phen_proc <- readRDS('ebird_NA_phen_proc.rds')

setwd(paste0(dir, 'Bird_phenology/Data/Processed/ebird_NA_phen_proc_species'))

for(i in 1:nsp)
{
  #i <- 1
  sdata <- dplyr::select(ebird_NA_phen_proc, 
                         YEAR, DAY, sjday, sjday2, 
                         sjday3, shr, cell6, species_list[i])
  names(sdata)[8] <- "detect"
  
  saveRDS(sdata, file = paste0('ebird_NA_phen_proc_',species_list[i], '.rds'))
}


# runtime -----------------------------------------------------------------

time <- proc.time() - tt
rtime <- round(time[3]/60, 2)

sink('1-runtime.txt')
cat(paste0('Runtime (minutes): ', rtime))
sink()
