################################
#Extract temperature data for NA
#
#Requires ~ 200GB RAM and 4 cores on cluster; ~ ? hours runtime
################################

#mean daily min temp for Feb-April, following Hurlbert and Liang 2012

#resources:
#http://rpubs.com/boyerag/297592
#https://daac.ornl.gov/workshops/NetCDF_webinar_08302017.html
#used this one vvv
#https://github.com/ornldaac/thredds_opendap_r_max_temperature/blob/master/opendap_r_v1.Rmd


# Load packages -----------------------------------------------------------

library(ncdf4)
library(dplyr)
library(dggridR)
library(abind)
library(foreach)
library(doParallel)


# top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# define function -----------------------------------------------------------

#input ncdf file for each year
#returns dataframe of hex sorted 
daymet_fun <- function(input, var = 'tmax', YEAR)
{
  st_time <- proc.time()[3]/60
  setwd('/home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Raw/Daymet_data')
  #setwd('~/Desktop/Bird_Phenology_Offline/Data/Daymet_data/')
  
  daymet_data_temp <- ncdf4::nc_open(input)
  
  #get lat/lons
  daymet_lats <- ncdf4::ncvar_get(daymet_data_temp, "lat")
  daymet_lons <- ncdf4::ncvar_get(daymet_data_temp, "lon")
  
  t_elapsed_1 <- round((proc.time()[3]/60 - st_time), 1)
  
  sink(paste0(dir, 'Bird_Phenology/Data/Processed/daymet/log3.txt'), append = TRUE)
  cat(paste0(var, ' lat/lon read-in complete: ', YEAR, 
             ' - ', t_elapsed_1, ' min \n'))
  
  #Feb, Mar, Apr, and FMA
  #calc julian day
  feb_start <- julian(as.Date(paste0(YEAR, '-02-01')), 
                      origin = as.Date(paste0(YEAR, '-01-01')))[1]
  mar_start <- julian(as.Date(paste0(YEAR, '-03-01')), 
                      origin = as.Date(paste0(YEAR, '-01-01')))[1]
  apr_start <- julian(as.Date(paste0(YEAR, '-04-01')), 
                      origin = as.Date(paste0(YEAR, '-01-01')))[1]
  
  #get data
  feb_d_temp <- ncdf4::ncvar_get(daymet_data_temp, var, start = c(1, 1, feb_start), 
                                 count = c(7814, 8075, 28))
  mar_d_temp <- ncdf4::ncvar_get(daymet_data_temp, var, start = c(1, 1, mar_start), 
                                 count = c(7814, 8075, 31))
  apr_d_temp <- ncdf4::ncvar_get(daymet_data_temp, var, start = c(1, 1, apr_start), 
                                 count = c(7814, 8075, 30))
  
  t_elapsed_2 <- round((proc.time()[3]/60 - t_elapsed_1), 1)
  cat(paste0(var, ' data read-in complete: ', YEAR, 
             ' - ', t_elapsed_2, ' min \n'))
  
  ncdf4::nc_close(daymet_data_temp)
  rm(daymet_data_temp)
  gc()
  
  Feb_temp <- apply(feb_d_temp[,,], c(1,2), mean)
  rm(feb_d_temp)
  Mar_temp <- apply(mar_d_temp[,,], c(1,2), mean)
  rm(mar_d_temp)
  Apr_temp <- apply(apr_d_temp[,,], c(1,2), mean)
  rm(apr_d_temp)
  FMA_temp_array <- abind::abind(Feb_temp, Mar_temp, Apr_temp, along = 3)
  FMA_temp <- apply(FMA_temp_array, c(1,2), mean)
  rm(FMA_temp_array)
  gc()
  
  t_elapsed_3 <- round((proc.time()[3]/60 - t_elapsed_2), 1)
  cat(paste0(var, ' averaging complete: ', YEAR, 
             ' - ', t_elapsed_3, ' min \n'))
  
  LL_daymet <- data.frame(LATS = as.vector(daymet_lats), 
                          LONS = as.vector(daymet_lons),
                          ROW_IND = rep(1:NROW(daymet_lats), times = NCOL(daymet_lats)),
                          COL_IND = rep(1:NCOL(daymet_lats), each = NROW(daymet_lats)),
                          F_temp = as.vector(Feb_temp), 
                          M_temp = as.vector(Mar_temp), 
                          A_temp = as.vector(Apr_temp), 
                          FMA_temp = as.vector(FMA_temp))
  
  rm(daymet_lats)
  rm(daymet_lons)
  rm(Feb_temp)
  rm(Mar_temp)
  rm(Apr_temp)
  rm(FMA_temp)
  
  
  #filter by relevant location
  f_LL_daymet <- dplyr::filter(LL_daymet, 
                               LONS > -100 & LONS < -50 & LATS > 24)
  
  rm(LL_daymet)
  gc()
  
  #create grid
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  
  #see which cells the lat/lons are associated with
  f_LL_daymet$hex_cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                                   in_lon_deg = f_LL_daymet$LONS, 
                                                   in_lat_deg = f_LL_daymet$LATS)[[1]]
  hex_cells <- sort(unique(f_LL_daymet$hex_cell))
  
  #average over daymet cells for each hexcell - 
  HEX_daymet <- data.frame(cell = hex_cells,
                           year = YEAR,
                           F_temp = NA,
                           M_temp = NA,
                           A_temp = NA,
                           FMA_temp = NA)
  colnames(HEX_daymet)[3:6] <- c(paste0('F_', var), paste0('M_', var),
                                 paste0('A_', var), paste0('FMA_', var))
  
  t_elapsed_4 <- round((proc.time()[3]/60 - t_elapsed_3), 1)
  cat(paste0(var, ' starting hex loop: ', YEAR, 
             ' - ', t_elapsed_4, ' min \n'))
  
  #fill df
  for (i in 1:length(hex_cells))
  {
    #i <- 1
    t_daymet <- dplyr::filter(f_LL_daymet, hex_cell == hex_cells[i])
    HEX_daymet[i,3] <- mean(t_daymet$F_temp, na.rm = TRUE)
    HEX_daymet[i,4] <- mean(t_daymet$M_temp, na.rm = TRUE)
    HEX_daymet[i,5] <- mean(t_daymet$A_temp, na.rm = TRUE)
    HEX_daymet[i,6] <- mean(t_daymet$FMA_temp, na.rm = TRUE)
  }
  
  return(HEX_daymet)
}




# run tmax ----------------------------------------------------------------

YEARS <- 2002:2018
#YEARS <- 2002:2003

doParallel::registerDoParallel(cores = 3)
OUT_tmax <- foreach::foreach(k = 1:length(YEARS), .combine = 'rbind') %dopar% 
  {
    tt_tmax <- daymet_fun(input = paste0('daymet_v3_tmax_', YEARS[k], '_na.nc4'), 
                          var = 'tmax', YEAR = YEARS[k])
    return(tt_tmax)
  }

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/daymet'))
saveRDS(OUT_tmax, 'daymet_hex_tmax.rds')
rm(OUT_tmax)
gc()
print('Finished tmax')



# run tmin ----------------------------------------------------------------

OUT_tmin <- foreach::foreach(k = 1:length(YEARS), .combine = 'rbind') %dopar%
  {
    tt_tmin <- daymet_fun(input = paste0('daymet_v3_tmin_', YEARS[k], '_na.nc4'),
                          var = 'tmin', YEAR = YEARS[k])
    return(tt_tmin)
  }

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/daymet'))
saveRDS(OUT_tmin, 'daymet_hex_tmin.rds')
rm(OUT_tmin)
gc()
print('Finished tmin')



# run precip --------------------------------------------------------------

OUT_precip <- foreach::foreach(k = 1:length(YEARS), .combine = 'rbind') %dopar%
  {
    tt_precip <- daymet_fun(input = paste0('daymet_v3_prcp_', YEARS[k], '_na.nc4'),
                     var = 'prcp', YEAR = YEARS[k])
    return(tt_precip)
  }

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/daymet'))
saveRDS(OUT_precip, 'daymet_hex_precip.rds')
print('Finished precip')
print('I completed!')