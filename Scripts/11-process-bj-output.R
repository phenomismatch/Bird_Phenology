######################
# 11 - Extract breeding dates from joint breeding IAR output
#
# Compiles dataframe with pre and post IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'


# other dir ---------------------------------------------------------------

br_in_dir <- 'br_IAR_input_2020-12-03'
juv_in_dir <- 'juv_IAR_input_2021-01-11'

IAR_out_dir <- 'bj_IAR_hm_2021-01-11'
master_out_dir <- 'bj_master_2021-01-11'

br_date <- substr(br_in_dir, start = 14, stop = 23)
juv_date <- substr(juv_in_dir, start = 15, stop = 24)
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
br_master <- readRDS(paste0('br_IAR_input_', br_date, '.rds'))

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', juv_in_dir))
juv_master <- readRDS(paste0('juv_IAR_input_', juv_date, '.rds'))


#cell info
cell_f1 <- dplyr::select(br_master, species, cell,
                         per_ovr, cell_lat, cell_lng,
                         breed_cell, other_cell)

cell_f2 <- dplyr::select(juv_master, species, cell,
                         per_ovr, cell_lat, cell_lng,
                         breed_cell, other_cell)

cell_c1 <- unique(rbind(cell_f1, cell_f2))

#pheno info
br_f <- dplyr::select(br_master, species, year, cell, br_GAM_mean, 
                      br_GAM_sd, VALID_br = VALID)

juv_f <- dplyr::select(juv_master, species, year, cell, juv_logis_mean, 
                       juv_logis_sd, VALID_juv = VALID)

#join
mrg1 <- dplyr::full_join(br_f, juv_f, by = c('species', 'year', 'cell'))
mrg2 <- dplyr::left_join(mrg1, cell_c1, by = c('species', 'cell'))

species <- as.character(read.table(paste0(dir, 'Bird_Phenology/Data/arr_species_list.txt'))[,1])
#species <- 'Vireo_olivaceus'


# create empty dataframes to fill -----------------------------------------

#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)

out <- data.frame(species = rep(NA, NROW(mrg2)), cell = NA, 
                  breed_cell = NA, other_cell = NA, 
                  cell_lat = NA, cell_lng = NA, per_ovr = NA, year = NA, 
                  br_GAM_mean = NA, br_GAM_sd = NA, VALID_br_GAM = NA, 
                  juv_logis_mean = NA, juv_logis_sd = NA, VALID_juv_logis = NA, 
                  bj_IAR_mean = NA, bj_IAR_sd = NA, 
                  sigma_beta0_mean = NA, sigma_beta0_sd = NA,
                  beta_gamma_mean = NA, beta_gamma_sd = NA,
                  alpha_mean = NA, alpha_sd = NA,
                  num_diverge = NA, max_Rhat = NA, min_neff = NA, 
                  minutes = NA, ITER = NA, novr = NA)


# run loop to fill empty dfs ----------------------------------------------------------------

counter <- 1
for (i in 1:length(species))
{
  #i <- 1
  #i <- 45
  
  #filter by species
  sp <- species[i]
  print(sp)
  
  #switch to out dir
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
  
  #if that species RDS object exists in dir
  if (length(grep(paste0(sp, '-bj-iar-hm-stan_output-', IAR_out_date, '.rds'), list.files())) > 0)
  {
    
    #NA for species/year/cells with non-valid GAM results
    br_na <- which(mrg2$VALID_br == FALSE)
    mrg2$br_GAM_mean[br_na] <- NA
    mrg2$br_GAM_sd[br_na] <- NA
    juv_na <- which(mrg2$VALID_juv == FALSE)
    mrg2$juv_logis_mean[juv_na] <- NA
    mrg2$juv_logis_sd[juv_na] <- NA
    
    #filter by year and species
    mrg3 <- dplyr::filter(mrg2, year >= 2002, year <= 2017, per_ovr >= 0.05, 
                          species == sp, breed_cell == TRUE, other_cell == FALSE)
    
    if (NROW(dplyr::filter(mrg3, !is.na(br_GAM_mean))) > 0)
    {
      agg_br <- aggregate(br_GAM_mean ~ year, data = mrg3, function(x) sum(!is.na(x)))
    } else {
      agg_br <- data.frame(year = NA, br_GAM_mean = NA)
    }
    if (NROW(dplyr::filter(mrg3, !is.na(juv_logis_mean))) > 0)
    {
      agg_juv <- aggregate(juv_logis_mean ~ year, data = mrg3, function(x) sum(!is.na(x)))
    } else {
      agg_juv <- data.frame(year = NA, juv_GAM_mean = NA)
    }
    
    #filter for valid years
    agg_mrg <- dplyr::full_join(agg_br, agg_juv, by = 'year')
    agg_mrg$j <- apply(agg_mrg[,2:3], 1, function(x) sum(x, na.rm = TRUE))
    vyrs <- dplyr::filter(agg_mrg, j >= 3)$year
    
    mrg4 <- dplyr::filter(mrg3, year %in% vyrs)
    
    #ensure that cells are ordered
    mrg5 <- dplyr::arrange(mrg4, species, cell, year)
    
    # #number of cell/years with overlapping data
    novr <- NROW(dplyr::filter(mrg5, !is.na(br_GAM_mean), !is.na(juv_logis_mean)))
    
    #read in IAR model output and input
    t_fit <- readRDS(paste0(sp, '-bj-iar-hm-stan_output-', IAR_out_date, '.rds'))
    t_data <- readRDS(paste0(sp, '-bj-iar-hm-stan_input-', IAR_out_date, '.rds'))
    
    #only cells and years that were modeled (to account for any lone cells that were dropped in 4-IAR-arr.R)
    f_in <- dplyr::filter(mrg5, cell %in% t_data$cells)
    
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
    
    #extract offset (alpha)
    alpha_mean <- MCMCvis::MCMCpstr(t_fit, params = 'alpha',
                                         func = mean)[[1]]
    alpha_sd <- MCMCvis::MCMCpstr(t_fit, params = 'alpha',
                                       func = sd)[[1]]
    
    #diagnostics
    num_diverge <- rstan::get_num_divergent(t_fit)
    model_summary <- MCMCvis::MCMCsummary(t_fit, excl = c('br_rep', 'juv_rep'), round = 3)
    max_Rhat <- max(as.vector(model_summary[, grep('Rhat', colnames(model_summary))]))
    min_neff <- min(as.vector(model_summary[, grep('n.eff', colnames(model_summary))]))
    
    #get runtime and ITER from results file
    ff <- paste0(sp, '-bj-iar-hm-stan_results-', IAR_out_date, '.txt')
    if (length(grep(ff, list.files())) > 0)
    {
      con <- file(ff, 'r')
      f3 <- readLines(con, n = 3)
      close(con)
      minutes <- round(as.numeric(strsplit(f3[2], ' ')[[1]][3]), 0)
      ITER <- as.numeric(strsplit(f3[3], ' ')[[1]][2])
    } else {
      minutes <- NA
      ITER <- NA
    }
    
    #loop through years
    for (j in 1:length(t_years))
    {
      #j <- 1
      print(paste0('species: ', sp, ', ',
                   'year: ', t_years[j]))
      
      t_f_in <- dplyr::filter(f_in, year == t_years[j])
      
      #insert FALSE where VALID is NA (where cells were not modeled by respective IAR)
      v_br_na <- which(is.na(t_f_in$VALID_br))
      v_juv_na <- which(is.na(t_f_in$VALID_juv))
      if (length(v_br_na) > 0)
      {
        t_f_in$VALID_br[v_br_na] <- FALSE
      }
      if (length(v_juv_na) > 0)
      {
        t_f_in$VALID_juv[v_juv_na] <- FALSE
      }
      
      
      t_full <- data.frame(t_f_in[,c('species', 'cell', 'breed_cell', 'other_cell')],
                           cell_lat = t_f_in$cell_lat,
                           cell_lng = t_f_in$cell_lng,
                           per_ovr = t_f_in$per_ovr,
                           year = t_f_in$year,
                           br_GAM_mean = t_f_in$br_GAM_mean,
                           br_GAM_sd = t_f_in$br_GAM_sd,
                           VALID_br_GAM = t_f_in$VALID_br,
                           juv_logis_mean = t_f_in$juv_logis_mean,
                           juv_logis_sd = t_f_in$juv_logis_sd,
                           VALID_juv_GAM = t_f_in$VALID_juv,
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
                           min_neff,
                           minutes,
                           ITER,
                           novr)
      
      #fill empty df
      out[counter:(counter + NROW(t_full) - 1),] <- t_full
      
      #advance counter
      counter <- counter + NROW(t_full)
      
    } #end year loop
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory. Species may not have met data requirements for IAR (e.g., 3 valid years of data).'))
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

#unique(out2[,c('species', 'minutes', 'ITER', 'max_Rhat', 'min_neff')])

