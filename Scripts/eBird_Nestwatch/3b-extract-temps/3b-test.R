################################
#Extract temperature data for NA
#
#Requires ~ 185GB RAM and 7 cores on cluster; ~ 2 hours runtime; eventually run on desktop
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
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'

setwd('/home/CAM/cyoungflesh/phenomismatch/Bird_Phenology/Data/Raw/Daymet_data')


# read in data -----------------------------------------------------------

tt <- proc.time()

YEARS <- 2002:2018



#doParallel::registerDoParallel(cores = 4)
#OUT <- foreach::foreach(k = 1:length(YEARS), .combine = 'rbind') %dopar% 
#  {
    k <- 1
    
    daymet_data_tmax <- ncdf4::nc_open(paste0('daymet_v3_tmax_', YEARS[k], '_na.nc4'))
    daymet_data_tmin <- ncdf4::nc_open(paste0('daymet_v3_tmin_', YEARS[k], '_na.nc4'))
    daymet_data_precip <- ncdf4::nc_open(paste0('daymet_v3_prcp_', YEARS[k], '_na.nc4'))
    
    #get lat/lons
    daymet_lats <- ncdf4::ncvar_get(daymet_data_tmax, "lat")
    daymet_lons <- ncdf4::ncvar_get(daymet_data_tmax, "lon")
    
    #Feb, Mar, Apr, and FMA
    #calc julian day
    feb_start <- julian(as.Date(paste0(YEARS[k], '-02-01')), 
                        origin = as.Date(paste0(YEARS[k], '-01-01')))[1]
    feb_end <- julian(as.Date(paste0(YEARS[k], '-02-28')), 
                      origin = as.Date(paste0(YEARS[k], '-01-01')))[1]
    mar_start <- julian(as.Date(paste0(YEARS[k], '-03-01')), 
                        origin = as.Date(paste0(YEARS[k], '-01-01')))[1]
    mar_end <- julian(as.Date(paste0(YEARS[k], '-03-31')), 
                      origin = as.Date(paste0(YEARS[k], '-01-01')))[1]
    apr_start <- julian(as.Date(paste0(YEARS[k], '-04-01')), 
                        origin = as.Date(paste0(YEARS[k], '-01-01')))[1]
    apr_end <- julian(as.Date(paste0(YEARS[k], '-04-30')), 
                      origin = as.Date(paste0(YEARS[k], '-01-01')))[1]
    
    #get data
    feb_d_tmax <- ncdf4::ncvar_get(daymet_data_tmax, 'tmax', start = c(1, 1, feb_start), count =c(7814, 8075, 28))
    feb_d_tmin <- ncdf4::ncvar_get(daymet_data_tmin, 'tmin', start = c(1, 1, feb_start), count =c(7814, 8075, 28))
    feb_d_precip <- ncdf4::ncvar_get(daymet_data_precip, 'prcp', start = c(1, 1, feb_start), count =c(7814, 8075, 28))
    
    mar_d_tmax <- ncdf4::ncvar_get(daymet_data_tmax, 'tmax', start = c(1, 1, mar_start), count =c(7814, 8075, 31))
    mar_d_tmin <- ncdf4::ncvar_get(daymet_data_tmin, 'tmin', start = c(1, 1, mar_start), count =c(7814, 8075, 31))
    mar_d_precip <- ncdf4::ncvar_get(daymet_data_precip, 'prcp', start = c(1, 1, mar_start), count =c(7814, 8075, 31))
    
    apr_d_tmax <- ncdf4::ncvar_get(daymet_data_tmax, 'tmax', start = c(1, 1, apr_start), count =c(7814, 8075, 30))
    apr_d_tmin <- ncdf4::ncvar_get(daymet_data_tmin, 'tmin', start = c(1, 1, apr_start), count =c(7814, 8075, 30))
    apr_d_precip <- ncdf4::ncvar_get(daymet_data_precip, 'prcp', start = c(1, 1, apr_start), count =c(7814, 8075, 30))
    
    Feb_tmax <- apply(feb_d_tmax[,,], c(1,2), mean)
    Mar_tmax <- apply(mar_d_tmax[,,], c(1,2), mean)
    Apr_tmax <- apply(apr_d_tmax[,,], c(1,2), mean)
    FMA_tmax_array <- abind::abind(Feb_tmax, Mar_tmax, Apr_tmax, along = 3)
    FMA_tmax <- apply(FMA_tmax_array, c(1,2), mean)
    
    Feb_tmin <- apply(feb_d_tmin[,,], c(1,2), mean)
    Mar_tmin <- apply(mar_d_tmin[,,], c(1,2), mean)
    Apr_tmin <- apply(apr_d_tmin[,,], c(1,2), mean)
    FMA_tmin_array <- abind::abind(Feb_tmin, Mar_tmin, Apr_tmin, along = 3)
    FMA_tmin <- apply(FMA_tmin_array, c(1,2), mean)
    
    Feb_precip <- apply(feb_d_precip[,,], c(1,2), mean)
    Mar_precip <- apply(mar_d_precip[,,], c(1,2), mean)
    Apr_precip <- apply(apr_d_precip[,,], c(1,2), mean)
    FMA_precip_array <- abind::abind(Feb_precip, Mar_precip, Apr_precip, along = 3)
    FMA_precip <- apply(FMA_precip_array, c(1,2), mean)
    
    
    #lat/lons, temps in data.frame
    LL_daymet <- data.frame(LATS = as.vector(daymet_lats), 
                            LONS = as.vector(daymet_lons),
                            ROW_IND = rep(1:NROW(daymet_lats), times = NCOL(daymet_lats)),
                            COL_IND = rep(1:NCOL(daymet_lats), each = NROW(daymet_lats)),
                            F_tmax = as.vector(Feb_tmax), 
                            M_tmax = as.vector(Mar_tmax), 
                            A_tmax = as.vector(Apr_tmax), 
                            FMA_tmax = as.vector(FMA_tmax),
                            F_tmin = as.vector(Feb_tmin), 
                            M_tmin = as.vector(Mar_tmin), 
                            A_tmin = as.vector(Apr_tmin), 
                            FMA_tmin = as.vector(FMA_tmin),
                            F_precip = as.vector(Feb_precip), 
                            M_precip = as.vector(Mar_precip), 
                            A_precip = as.vector(Apr_precip), 
                            FMA_precip = as.vector(FMA_precip))
    
    #filter by relevant location
    f_LL_daymet <- dplyr::filter(LL_daymet, 
                                 LONS > -100 & LONS < -50 & LATS > 24)
    
    #create grid
    hexgrid6 <- dggridR::dgconstruct(res = 6)
    
    #see which cells the lat/lons are associated with
    f_LL_daymet$hex_cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                                     in_lon_deg = f_LL_daymet$LONS, 
                                                     in_lat_deg = f_LL_daymet$LATS)[[1]]
    
    #check to make sure things are sorted properly
    # mylat <- f_LL_daymet$LATS[500000]
    # mylon <- f_LL_daymet$LONS[500000]
    # f_LL_daymet[500000,] == which(daymet_lats == mylat & daymet_lons == mylon, arr.ind = TRUE)
    
    hex_cells <- sort(unique(f_LL_daymet$hex_cell))
    
    #average over daymet cells for each hexcell - 
    HEX_daymet <- data.frame(cell = hex_cells,
                             year = YEARS[k],
                             F_tmax = NA,
                             M_tmax = NA,
                             A_tmax = NA,
                             FMA_tmax = NA,
                             F_tmin = NA,
                             M_tmin = NA,
                             A_tmin = NA,
                             FMA_tmin = NA,
                             F_precip = NA,
                             M_precip = NA,
                             A_precip = NA,
                             FMA_precip = NA)
    
    
    pb <- txtProgressBar(min = 0, max = length(hex_cells), style = 3)
    #fill df
    for (i in 1:length(hex_cells))
    {
      #i <- 1
      t_daymet <- dplyr::filter(f_LL_daymet, hex_cell == hex_cells[i])
      HEX_daymet$F_tmax[i] <- mean(t_daymet$F_tmax, na.rm = TRUE)
      HEX_daymet$M_tmax[i] <- mean(t_daymet$M_tmax, na.rm = TRUE)
      HEX_daymet$A_tmax[i] <- mean(t_daymet$A_tmax, na.rm = TRUE)
      HEX_daymet$FMA_tmax[i] <- mean(t_daymet$FMA_tmax, na.rm = TRUE)
      
      HEX_daymet$F_tmin[i] <- mean(t_daymet$F_tmin, na.rm = TRUE)
      HEX_daymet$M_tmin[i] <- mean(t_daymet$M_tmin, na.rm = TRUE)
      HEX_daymet$A_tmin[i] <- mean(t_daymet$A_tmin, na.rm = TRUE)
      HEX_daymet$FMA_tmin[i] <- mean(t_daymet$FMA_tmin, na.rm = TRUE)
      
      HEX_daymet$F_precip[i] <- mean(t_daymet$F_precip, na.rm = TRUE)
      HEX_daymet$M_precip[i] <- mean(t_daymet$M_precip, na.rm = TRUE)
      HEX_daymet$A_precip[i] <- mean(t_daymet$A_precip, na.rm = TRUE)
      HEX_daymet$FMA_precip[i] <- mean(t_daymet$FMA_precip, na.rm = TRUE)
      
      print(paste0(YEARS[k]))
      setTxtProgressBar(pb, i)
    }
    
#    return(HEX_daymet)
#  }
#proc.time() - tt


# write to rds file -------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/daymet'))

saveRDS(HEX_daymet, 'test_daymet_hex.rds')

print('I completed!')
