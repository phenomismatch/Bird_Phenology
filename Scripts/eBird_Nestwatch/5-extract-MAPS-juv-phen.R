######################
# 5 - process young hitting nets MAPS
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# out dir ------------------------------------------------------------

arr_br_dir <- 'arr_br_2019-08-19'


# load packages -----------------------------------------------------------

library(dplyr)
library(plyr)
library(dggridR)


# read in data ------------------------------------------------------------

#read in - RDS create with 1-query-db.R in wing_chord_changes project
setwd(paste0(dir, 'Bird_Phenology/Data'))

data <- readRDS('MAPS-age-filled.rds')

#filter for only age 0 birds (birth year) and exclude NA bands
data_a0 <- dplyr::filter(data, true_age == 0, !is.na(band_id))



# create grid -------------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6)
data_a0$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                       in_lon_deg = data_a0$lng, 
                                       in_lat_deg = data_a0$lat)[[1]]




# filter by species/cell/year ---------------------------------------------

ind_out <- data.frame(sci_name = rep(NA, NROW(data_a0)), 
                      cell = NA, 
                      year = NA, 
                      band_id = NA, 
                      jday = NA)

counter <- 1
#SPECIES
usp <- unique(data_a0$sci_name)
#progress bar
pb <- txtProgressBar(min = 0, max = length(usp), style = 3)

for (s in 1:length(usp))
{
  #s <- 1
  temp <- dplyr::filter(data_a0, sci_name == usp[s])
  
  #CELL
  ucell <- unique(temp$cell)
  for (i in 1:length(ucell))
  {
    #i <- 1
    temp2 <- dplyr::filter(temp, cell == ucell[i])
    
    uyear <- unique(temp2$year)
    
    #YEAR
    for (j in 1:length(uyear))
    {
      #j <- 4
      temp3 <- dplyr::filter(temp2, year == uyear[j])
      
      #get # of captures for each band_id
      bid_cnt <- plyr::count(temp3, 'band_id')
      #which band_id were captured > 1 time
      mc_id <- dplyr::filter(bid_cnt, freq > 1)[,1]
      
      if (NROW(mc_id) > 0)
      {
        #INDIVIDUAL
        for (k in 1:length(mc_id))
        {
          #k <- 1
          temp4 <- dplyr::filter(temp3, band_id == mc_id[k])
          
          ind_out$sci_name[counter] <- usp[s]
          ind_out$cell[counter] <- ucell[i]
          ind_out$year[counter] <- uyear[j]
          ind_out$band_id[counter] <- mc_id[k]
          ind_out$jday[counter] <- min(temp4$day)
          
          counter <- counter + 1
        }
      }
    }
  }
  setTxtProgressBar(pb, s)
}

close(pb)


#trim NA from end of df
sn_na <- min(which(is.na(ind_out$sci_name)))
ind_out2 <- ind_out[-c(sn_na:NROW(ind_out)),]



# save as RDS -------------------------------------------------------------

#create output image dir if it doesn't exist
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_br_dir)),
       dir.create(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_br_dir)),
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_br_dir))

saveRDS(ind_out2, 'arr_br_in.rds')
