################################
#Extract temperature data for NA
#
#Set for 9 cores
################################

#mean min temp for Feb-April

#resources:
#http://rpubs.com/boyerag/297592
#https://daac.ornl.gov/workshops/NetCDF_webinar_08302017.html
#used this one vvv
#https://github.com/ornldaac/thredds_opendap_r_max_temperature/blob/master/opendap_r_v1.Rmd


# Load packages -----------------------------------------------------------

library(ncdf4)
library(dplyr)
library(dggridR)
library(doParallel)
library(foreach)



# top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# args --------------------------------------------------------------------

#args is year - 2002 to 2017
args <- commandArgs(trailingOnly = TRUE)
#args <- 2002



# read in data -----------------------------------------------------------

#only through 2017 at the moment

daymet_url <- paste0('https://thredds.daac.ornl.gov/thredds/dodsC/ornldaac/1328/', 
                     args, '/daymet_v3_tmax_', args, '_na.nc4')
daymet_data <- ncdf4::nc_open(daymet_url)

#get lat/lons
daymet_lats <- ncdf4::ncvar_get(daymet_data, "lat")
daymet_lons <- ncdf4::ncvar_get(daymet_data, "lon")

#pair up lat/lons - filled by columns
LL_daymet <- data.frame(LATS = as.vector(daymet_lats), 
                        LONS = as.vector(daymet_lons),
                        ROW_IND = rep(1:NROW(daymet_lats), times = NCOL(daymet_lats)),
                        COL_IND = rep(1:NCOL(daymet_lats), each = NROW(daymet_lats)))

#filter by relevant location
f_LL_daymet <- dplyr::filter(LL_daymet, 
                             LONS > -100 & LONS < -50 & LATS > 26)

#create grid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#see which cells the lat/lons are associated with
f_LL_daymet$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                             in_lon_deg = f_LL_daymet$LONS, 
                                             in_lat_deg = f_LL_daymet$LATS)[[1]]

#check to make sure things are sorted properly
# mylat <- f_LL_daymet$LATS[500000]
# mylon <- f_LL_daymet$LONS[500000]
# f_LL_daymet[500000,] == which(daymet_lats == mylat & daymet_lons == mylon, arr.ind = TRUE)

#get unique cells in area of interest
cells <- sort(unique(f_LL_daymet$cell))

#get julian days - Feb 1 - April 30 (following Hurlbert and Liang PLOS One)
start_jday <- as.numeric(format(as.Date(paste0('01-02-', YEAR), 
                                        format = '%d-%m-%Y'), format = '%j'))
end_jday <- as.numeric(format(as.Date(paste0('30-04-', YEAR), 
                                      format = '%d-%m-%Y'), format = '%j'))

#spring mean for dm cell
f_LL_daymet$mean_dmcell_val <- NA
#average of all dm cells within that hex cell
f_LL_daymet$mean_hex_val <- NA



# process cells in parallel -----------------------------------------------

tt <- proc.time()

doParallel::registerDoParallel(cores = 3)
OUT <- foreach::foreach(i = 1:10, .combine = 'rbind') %dopar% 
#OUT <- foreach::foreach(i = 1:length(cells), .combine = 'rbind') %dopar% 
{
  #i <- 1
  dm_temp <- dplyr::filter(f_LL_daymet, cell == cells[i])
  
  #extract values over temporal period of interest for each daymet cell that falls within hex grid cell
  dm_cell_tmax <- rep(NA, NROW(dm_temp))
  for (j in 1:NROW(dm_temp))
  {
    #j <- 1
    #cells over water will have NA
    t_tmax <- ncdf4::ncvar_get(daymet_data, "tmax", 
                               start = c(dm_temp$ROW_IND[j], dm_temp$COL_IND[j], start_jday), 
                               count = c(1, 1, (end_jday - start_jday)))
    
    #take mean of values for that period
    mt_tmax <- mean(t_tmax, na.rm = TRUE)
    dm_cell_tmax[j] <- mt_tmax
  }
  
  dm_temp$mean_dmcell_val <- dm_cell_tmax
  dm_temp$mean_hex_val <- mean(dm_cell_tmax, na.rm = TRUE)
  #ind <- which(f_LL_daymet$cell == cells[i])
  #f_LL_daymet[ind,'mean_dmcell_val'] <- dm_cell_tmax
  #f_LL_daymet[ind,'mean_hex_val'] <- mean(dm_cell_tmax, na.rm = TRUE)
  
  return(dm_temp)
}


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

saveRDS(OUT, 'daymet_hex_', args, '.rds')


#3-4 hours?
#submit each year as its own job - 16 jobs


#test - works as expected
# t_df <- data.frame(group = rep(1:5, each = 2), num1 = rnorm(10), num2 = NA)
# 
# ttt_out <- foreach::foreach(i = 1:5, .combine = 'rbind') %dopar% 
# {
#   #i <- 1
#   test_temp <- dplyr::filter(t_df, group == i)
#   test_temp$num2 <- mean(test_temp$num1)
#   return(test_temp)
# }

