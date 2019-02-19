################################
#Extract temperature data for NA
#
#Requires ~ 185GB RAM and 7 cores; ~ 2 hours runtime
################################

#mean min temp for Feb-April, following Hurlbert and Liang 2012

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



# read in data -----------------------------------------------------------

tt <- proc.time()

#only through 2017 at the moment
YEARS <- 2002:2017

#for (k in 1:length(YEARS))

doParallel::registerDoParallel(cores = 6)
OUT <- foreach::foreach(k = 1:length(YEARS), .combine = 'rbind') %dopar% 
{
  #k <- 1

  daymet_url <- paste0('https://thredds.daac.ornl.gov/thredds/dodsC/ornldaac/1345/daymet_v3_tmax_monavg_', YEARS[k], '_na.nc4')

  daymet_data <- ncdf4::nc_open(daymet_url)

  #get lat/lons
  daymet_lats <- ncdf4::ncvar_get(daymet_data, "lat")
  daymet_lons <- ncdf4::ncvar_get(daymet_data, "lon")

  #get temp data
  daymet_tmax <- ncdf4::ncvar_get(daymet_data, 'tmax')

  #Feb, Mar, Apr, and FMA
  Feb_tmax <- daymet_tmax[,,2]
  Mar_tmax <- daymet_tmax[,,3]
  Apr_tmax <- daymet_tmax[,,4]
  FMA_array <- abind::abind(Feb_tmax, Mar_tmax, Apr_tmax, along = 3)
  FMA_tmax <- apply(FMA_array, c(1,2), mean)


  #lat/lons, temps in data.frame
  LL_daymet <- data.frame(LATS = as.vector(daymet_lats), 
                          LONS = as.vector(daymet_lons),
                          ROW_IND = rep(1:NROW(daymet_lats), times = NCOL(daymet_lats)),
                          COL_IND = rep(1:NCOL(daymet_lats), each = NROW(daymet_lats)),
                          F_tmax = as.vector(Feb_tmax), 
                          M_tmax = as.vector(Mar_tmax), 
                          A_tmax = as.vector(Apr_tmax), 
                          FMA_tmax = as.vector(FMA_tmax))

  #filter by relevant location
  f_LL_daymet <- dplyr::filter(LL_daymet, 
                               LONS > -100 & LONS < -50 & LATS > 26)

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
                           HC_F_tmax = NA,
                           HC_M_tmax = NA,
                           HC_A_tmax = NA,
                           HC_FMA_tmax = NA)

  
  pb <- txtProgressBar(min = 0, max = length(hex_cells), style = 3)
  #fill df
  for (i in 1:length(hex_cells))
  {
    #i <- 1
    t_daymet <- dplyr::filter(f_LL_daymet, hex_cell == hex_cells[i])
    HEX_daymet$HC_F_tmax[i] <- mean(t_daymet$F_tmax, na.rm = TRUE)
    HEX_daymet$HC_M_tmax[i] <- mean(t_daymet$M_tmax, na.rm = TRUE)
    HEX_daymet$HC_A_tmax[i] <- mean(t_daymet$A_tmax, na.rm = TRUE)
    HEX_daymet$HC_FMA_tmax[i] <- mean(t_daymet$FMA_tmax, na.rm = TRUE)

    print(paste0(YEARS[k]))
    setTxtProgressBar(pb, i)
  }
  
  return(HEX_daymet)
}
proc.time() - tt

  
  

# write to rds file -------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/daymet'))

saveRDS(OUT, 'daymet_hex_tmax.rds')

