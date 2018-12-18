####################
# 6 - arrival date (IAR) ~ nesting date (NW)
#
# *Nestwatch data
####################



# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(dggridR)
library(doParallel)
library(foreach)


# import IAR data ---------------------------------------------------------


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

IAR_data <- readRDS(paste0('master_arrival_', IAR_out_date, '.rds'))



# import eBird species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)

#remove underscore and coerce to vector
species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))

#combine species names into a single string with quotes
SL <- paste0("'", species_list_i2, "'", collapse = ", ")



# create hex grid ---------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6) 



# Nestwatch ---------------------------------------------------------------

#each entry is a different nest in the data, right? Not multiple observers recording the same nest?

setwd(paste0(dir, 'Bird_Phenology/Data/Nestwatch'))


nw_data <- read.csv('Nestwatch_2018_1026.csv')

#add grid cell
nw_data$CELL <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                         in_lon_deg = nw_data$LONGITUDE, 
                                         in_lat_deg = nw_data$LATITUDE)[[1]]

#add year
nw_data$YEAR <- substr(nw_data$PROJ_PERIOD_ID, start = 4, stop = 7)


###########
#add latin species names to nw data (rather than 6 letter codes)
###########


specices <- unique(IAR_data$species)


#combine nest watch data with IAR data
nw_IAR <- data.frame()
for (i in 1:length(species))
{
  #i <- 96 #Vireo olivaceus
  
  #filter by species
  sp <- species[i]
  t_IAR <- dplyr::filter(IAR_data, species == sp)
  t_nw <- dplyr::filter(nw_data, species == sp)
  
  t_yrs <- unique(f_IAR$year)
  
  for (j in 1:t_yrs)
  {
    t2_IAR <- dplyr::filter(t_IAR, year == t_yrs[j])
    t2_nw <- dplyr::filter(t_nw, year == t_yrs[j])
    
    
    #if there are enough data for that species, year, cell
    #by cell
    t_cells <- unique(t2_IAR$cell)
    
    
  }
}

#hierarchical model to relate arrival date (from IAR model) to nesting date (from Nest Watch)

#for species/year/cell

#obs_IAR_arrival ~ N(true_IAR_arrival, IAR_uncertainty)
#obs_NW_nesting ~ N(true_NW_nesting, nw_uncertainty)


#tru_NW_nesting ~ N(mu, sd)
#mu = alpha + beta * true_IAR_nesting











