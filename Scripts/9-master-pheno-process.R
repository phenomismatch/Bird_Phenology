######################
# 12 - Compile master pheno file
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'


# other dir ---------------------------------------------------------------

arr_master_date <- '2020-07-21'
br_master_date <- '2021-03-29'
NW_RUN_DATE <- '2020-12-03'
pro_date <- '2021-03-30'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(rstan)
library(dplyr)
library(dggridR)
library(rgdal)
library(sp)
library(rgeos)


# Load data ------------------------------------------------------

#arrival
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/arrival_master_', arr_master_date))
arr_master <- readRDS(paste0('arrival_master_', arr_master_date, '.rds'))

#breeding
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/breeding_master_', br_master_date))
br_master <- readRDS(paste0('breeding_master_', br_master_date, '.rds'))

#NW intervals
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
nw_pi <- readRDS(paste0('Nestwatch_pheno_dates-', NW_RUN_DATE, '.rds'))


# merge -------------------------------------------------------------------

arr_master_2 <- dplyr::select(arr_master, species, cell, mig_cell, breed_cell, 
                              cell_lat, cell_lng, per_ovr, year, arr_GAM_mean, 
                              arr_GAM_sd, VALID_arr_GAM = VALID_GAM, arr_IAR_mean, 
                              arr_IAR_sd)

br_master_E <- dplyr::filter(br_master, metric == 'E')
br_master_Y <- dplyr::filter(br_master, metric == 'Y')
br_master_F <- dplyr::filter(br_master, metric == 'F')

br_master_E_2 <- dplyr::select(br_master_E, species, cell, 
                             br_E_GAM_mean = br_GAM_mean, br_E_GAM_sd = br_GAM_sd, 
                             VALID_br_E_GAM = VALID)
br_master_Y_2 <- dplyr::select(br_master_Y, species, cell, 
                               br_Y_GAM_mean = br_GAM_mean, br_Y_GAM_sd = br_GAM_sd, 
                               VALID_br_Y_GAM = VALID)
br_master_F_2 <- dplyr::select(br_master_F, species, cell, 
                               br_F_GAM_mean = br_GAM_mean, br_F_GAM_sd = br_GAM_sd, 
                               VALID_br_F_GAM = VALID)

nw_pi_2 <- dplyr::select(nw_pi, species, med_LH_imp, med_LF_imp)

#join
mrg <- dplyr::full_join(arr_master_2, br_master_E_2, by = c('species', 'cell')) %>%
  dplyr::full_join(br_master_Y_2, by = c('species', 'cell')) %>%
  dplyr::full_join(br_master_F_2, by = c('species', 'cell')) %>%
  dplyr::left_join(nw_pi_2, by = 'species')


# write to RDS ------------------------------------------------------------

#create dir if doesn't exist
ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Data/processed/pheno_master_', pro_date)), 
       dir.create(paste0(dir, 'Bird_Phenology/Data/processed/pheno_master_', pro_date)), 
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Data/processed/pheno_master_', pro_date))
saveRDS(mrg, paste0('pheno_master_', pro_date, '.rds'))


# pheno windows -----------------------------------------------------------

# PERIODS
#========
# ARRIVAL: IAR arr model
# FLEDGE: IAR br model
# LAY: FLEDGE - LF_imp (from NW)
# HATCH: FLEDGE - (LF_imp - LH_imp) (from NW)


# #filter by species
# sp <- 'Vireo_olivaceus'
# 
# f1 <- dplyr::filter(mrg2, species == sp)
# br_years <- unique(dplyr::filter(f1, !is.na(br_IAR_mean))$year)
# f2 <- dplyr::filter(f1, year %in% br_years)
# 
# #function to calculate cross-year means for ind cells
# c_fun <- function(CELL)
# {
#   c1 <- dplyr::filter(f2, cell == CELL)
#   
#   mean_ARR <- mean(c1$arr_IAR_mean, na.rm = TRUE)
#   sd_ARR <- sd(c1$arr_IAR_mean, na.rm = TRUE)
#   
#   mean_FLEDGE <- mean(c1$br_IAR_mean, na.rm = TRUE)
#   sd_FLEDGE <- sd(c1$br_IAR_mean, na.rm = TRUE)
#   
#   mean_LAY <- mean_FLEDGE - c1$med_LF_imp[1]
#   
#   mean_HATCH <- mean_FLEDGE - (c1$med_LF_imp[1] - c1$med_LH_imp[1])
#  
#   out <- data.frame(cell = CELL, 
#                     mean_ARR = round(mean_ARR, 0),
#                     sd_ARR = round(sd_ARR, 1),
#                     mean_LAY = round(mean_LAY, 0),
#                     sd_LAY = round(sd_FLEDGE, 1),
#                     mean_HATCH = round(mean_HATCH, 0),
#                     sd_HATCH = round(sd_FLEDGE, 1),
#                     mean_FLEDGE = round(mean_FLEDGE, 0),
#                     sd_FLEDGE = round(sd_FLEDGE, 1))
#   
#   return(out)
# }
# 
# #run function
# 
# #filter by cell
# # 595-purple (boston)
# # 618-blue (chicago)
# # 676 green (DC area)
# 
# dd <- rbind(c_fun(CELL = 595),
#             c_fun(CELL = 618),
#             c_fun(CELL = 676))
# sp_pheno_dates <- data.frame(species = sp, dd)
# 
# setwd('~/Desktop')
# write.csv(sp_pheno_dates, paste0('V_olivaceus_pheno_dates_', pro_date, '.csv'))

