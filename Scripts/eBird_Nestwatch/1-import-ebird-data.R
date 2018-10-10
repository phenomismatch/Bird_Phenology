######################
# 1 - import ebird data
#
# Import data from eBird
#
# Formerly Ebird_import.R
######################


# This script imports data from the eBird reference dataset (ERD) for a vector of species, defined by 
# species_list and a vector of years defined by years.  The script is designed to interact with the ERD
# in the format downloaded directly from eBird, with a separate folder for each year from 2002 onwards.


cy_dir <- '~/Google_Drive/R/'



# Load packages -----------------------------------------------------------

library(data.table)



# set wd ------------------------------------------------------------------

setwd(paste0(cy_dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list <- read.table('eBird_species_list.txt')



# Set filter parameters ---------------------------------------------------

years <- 2002:2016
nyr <- length(years)

nsp <- NROW(species_list)

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

setwd(paste0(cy_dir, 'Bird_Phenology/Data/Raw'))

# For each year, import the data, and store in a list:
dataList <- as.list(rep(NA, nyr))

for(i in 1:nyr)
{
  assign(paste0("data", years[i]), 
         read.csv(paste0("eBird/ERD2016SS/", years[i], "/checklists_NA_birdPhen_time.csv"),
                  quote = ""))
  
  dataList[[i]] <- eval(parse(text = paste0("data", years[i])))
  rm(list = paste("data", years[i], sep=""))
  print(years[i])
}


# save to rds object ------------------------------------------------------

data_NA_birdPhen <- data.table::rbindlist(dataList)

setwd(paste0(cy_dir, 'Bird_Phenology/Data/Processed'))

saveRDS(data_NA_birdPhen, file = 'ebird_NA_phen.rds')


