######################
# 11 - Extract breeding dates from joint breeding IAR output
#
# Compiles dataframe with pre and post IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'


# other dir ---------------------------------------------------------------

br_in_dir <- 'breeding_master_2020-06-04'
juv_in_dir <- 'juv_master_2020-06-04'

IAR_out_dir <- 'bj_IAR_hm_2020-08-28'
master_out_dir <- 'bj_master_2020-08-28'

br_date <- substr(br_in_dir, start = 17, stop = 26)
juv_date <- substr(juv_in_dir, start = 12, stop = 21)
IAR_out_date <- substr(IAR_out_dir, start = 11, stop = 20)


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(rstan)
library(dplyr)
library(dggridR)
library(rgdal)
library(sp)
library(rgeos)


# Filter data by species/years ------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', br_in_dir))
br_master <- readRDS(paste0('breeding_master_', br_date, '.rds'))

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', juv_in_dir))
juv_master <- readRDS(paste0('juv_master_', juv_date, '.rds'))

br_f <- dplyr::select(br_master, species, year, cell, br_GAM_mean, 
                       br_GAM_sd, VALID)

juv_f <- dplyr::select(juv_master, species, year, cell, juv_GAM_mean, 
                        juv_GAM_sd, breed_cell, other_cell, VALID, 
                        per_ovr, cell_lat, cell_lng)

#join
mrg1 <- dplyr::full_join(br_f, juv_f, by = c('species', 'year', 'cell'))


species <- as.character(read.table(paste0(dir, 'Bird_Phenology/Data/IAR_species_list.txt'))[,1])
#species <- 'Vireo_olivaceus'


# create empty dataframes to fill -----------------------------------------

#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)

out <- data.frame(species = rep(NA, NROW(mrg1)), cell = NA, 
                  breed_cell = NA, other_cell = NA, 
                  cell_lat = NA, cell_lng = NA, per_ovr = NA, year = NA, 
                  br_GAM_mean = NA, br_GAM_sd = NA, VALID_br_GAM = NA, 
                  juv_GAM_mean = NA, juv_GAM_sd = NA, VALID_juv_GAM = NA, 
                  bj_IAR_mean = NA, bj_IAR_sd = NA, 
                  sigma_beta0_mean = NA, sigma_beta0_sd = NA,
                  beta_gamma_mean = NA, beta_gamma_sd = NA,
                  alpha_mean = NA, alpha_sd = NA,
                  num_diverge = NA, max_Rhat = NA, min_neff = NA)


# run loop to fill empty dfs ----------------------------------------------------------------

counter <- 1
for (i in 1:length(species))
{
  #i <- 35
  #i <- 45
  
  #filter by species
  sp <- species[i]
  print(sp)
  
  #switch to out dir
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
  
  #if that species RDS object exists in dir
  if (length(grep(paste0(sp, '-bj-iar-hm-stan_output-', IAR_out_date, '.rds'), list.files())) > 0)
  {
    
    # filter and read in data -------------------------------------------------
    mrg2 <- dplyr::filter(mrg1, species == sp, year >= 2002, year <= 2017, 
                            per_ovr >= 0.05, breed_cell == TRUE, other_cell == FALSE)
    
    mrg2_temp <- mrg2
    br_na <- which(mrg2_temp$VALID.x == FALSE)
    mrg2_temp$br_GAM_mean[br_na] <- NA
    mrg2_temp$br_GAM_sd[br_na] <- NA
    juv_na <- which(mrg2_temp$VALID.y == FALSE)
    mrg2_temp$juv_GAM_mean[juv_na] <- NA
    mrg2_temp$juv_GAM_sd[juv_na] <- NA
    
    
    agg_br <- aggregate(br_GAM_mean ~ year, data = mrg2_temp, function(x) sum(!is.na(x)))
    agg_juv <- aggregate(juv_GAM_mean ~ year, data = mrg2_temp, function(x) sum(!is.na(x)))
    
    #filter for valid years
    agg_mrg <- dplyr::full_join(agg_br, agg_juv, by = 'year')
    agg_mrg$j <- apply(agg_mrg[,2:3], 1, function(x) sum(x, na.rm = TRUE))
    vyrs <- dplyr::filter(agg_mrg, j >=3)$year
    
    mrg3 <- dplyr::filter(mrg2, year %in% vyrs)
    
    #read in IAR model output and input
    t_fit <- readRDS(paste0(sp, '-bj-iar-hm-stan_output-', IAR_out_date, '.rds'))
    t_data <- readRDS(paste0(sp, '-bj-iar-hm-stan_input-', IAR_out_date, '.rds'))
    
    #only cells and years that were modeled (to account for any lone cells that were dropped in 4-IAR-arr.R)
    f_in <- dplyr::filter(mrg3, cell %in% t_data$cells)
    
    t_cells <- unique(f_in$cell)
    t_years <- unique(f_in$year)
    
    # extract posteriors ------------------------------------------------------
    
    #extract median and sd for IAR arrival dates
    fit_mean <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = mean)[[1]]
    fit_med <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = median)[[1]]
    fit_sd <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = sd)[[1]]
    
    #extract variance year effect (sigma_beta0)
    sigma_beta0_mean <- MCMCvis::MCMCpstr(t_fit, params = 'sigma_beta0', func = mean)[[1]]
    sigma_beta0_sd <- MCMCvis::MCMCpstr(t_fit, params = 'sigma_beta0', func = sd)[[1]]
    
    #extract lat gradient (beta_gamma)
    beta_gamma_mean <- MCMCvis::MCMCpstr(t_fit, params = 'beta_gamma',
                                         func = mean)[[1]]
    beta_gamma_sd <- MCMCvis::MCMCpstr(t_fit, params = 'beta_gamma',
                                       func = sd)[[1]]
    
    #extract phenological offset (alpha)
    alpha_mean <- MCMCvis::MCMCpstr(t_fit, params = 'alpha',
                                         func = mean)[[1]]
    alpha_sd <- MCMCvis::MCMCpstr(t_fit, params = 'alpha',
                                       func = sd)[[1]]
    
    #diagnostics
    num_diverge <- rstan::get_num_divergent(t_fit)
    model_summary <- MCMCvis::MCMCsummary(t_fit, excl = 'y_rep', round = 3)
    max_Rhat <- max(as.vector(model_summary[, grep('Rhat', colnames(model_summary))]))
    min_neff <- min(as.vector(model_summary[, grep('n.eff', colnames(model_summary))]))
    
    #loop through years
    for (j in 1:length(t_years))
    {
      #j <- 4
      print(paste0('species: ', sp, ', ',
                   'year: ', t_years[j]))
      
      t_f_in <- dplyr::filter(f_in, year == t_years[j])
      
      t_full <- data.frame(t_f_in[,c('species', 'cell', 'breed_cell', 'other_cell')],
                           cell_lat = t_f_in$cell_lat,
                           cell_lng = t_f_in$cell_lng,
                           per_ovr = t_f_in$per_ovr,
                           year = t_f_in$year,
                           br_GAM_mean = t_f_in$br_GAM_mean,
                           br_GAM_sd = t_f_in$br_GAM_sd,
                           VALID_br_GAM = t_f_in$VALID.x,
                           juv_GAM_mean = t_f_in$juv_GAM_mean,
                           juv_GAM_sd = t_f_in$juv_GAM_sd,
                           VALID_juv_GAM = t_f_in$VALID.y,
                           bj_IAR_mean = fit_mean[j,],
                           bj_IAR_sd = fit_sd[j,],
                           sigma_beta0_mean,
                           sigma_beta0_sd,
                           beta_gamma_mean,
                           beta_gamma_sd,
                           alpha_mean,
                           alpha_sd,
                           num_diverge,
                           max_Rhat,
                           min_neff)
      
      #fill empty df
      out[counter:(counter + NROW(t_full) - 1),] <- t_full
      
      #advance counter
      counter <- counter + NROW(t_full)
      
    } #end year loop
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop

#remove zeros
if (sum(is.na(out$species)) > 0)
{
  out2 <- out[-c(min(which(is.na(out$species))):NROW(out)),]
} else {
  out2 <- out
}


# write to file -----------------------------------------------------------

#create dir if it does not exist

ifelse(!dir.exists(paste0(dir, 'Bird_Phenology/Data/Processed/', master_out_dir)),
       dir.create(paste0(dir, 'Bird_Phenology/Data/Processed/', master_out_dir)),
       FALSE)

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', master_out_dir))

saveRDS(out2, file = paste0('bj_master_', IAR_out_date, '.rds'))

print('I completed!')

