######################
# Extract AOS pheno trends model output
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

in_dir <- '/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/trends_output_2019-06-16'


# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(dggridR)


# setwd -------------------------------------------------------------------

setwd(in_dir)

nc_in_dir <- nchar(in_dir)
in_date <- substr(in_dir, start = (nc_in_dir - 9), stop = nc_in_dir)


# Filter data by species/years ------------------------------------------------------

files <- list.files()
species <- unique(sapply(strsplit(files, split = '-'), head , 1))

out <- data.frame()
for (i in 1:length(species))
{
  #i <- 1
  
  #filter by species
  sp <- species[i]
  
  #if that species RDS object exists in dir
  #if (length(grep(paste0(sp, '-', in_date, '-pheno_trends_stan_output.rds'), list.files())) > 0)
  if (length(grep(paste0(sp, '-2019-05-26-pheno_trends_stan_output.rds'), list.files())) > 0)
  {
    #read in model output
    #t_fit <- readRDS(paste0(sp, '-', in_date, '-pheno_trends_stan_output.rds'))
    t_fit <- readRDS(paste0(sp, '-2019-05-26-pheno_trends_stan_output.rds'))
    
    #add year range and cell range (maybe range area calculated from range maps)?
    
    mn_mu_alpha <- MCMCvis::MCMCpstr(t_fit, params = 'mu_alpha', func = mean)[[1]]
    sd_mu_alpha <- MCMCvis::MCMCpstr(t_fit, params = 'mu_alpha', func = sd)[[1]]
    
    mn_sigma <- MCMCvis::MCMCpstr(t_fit, params = 'sigma', func = mean)[[1]]
    sd_sigma <- MCMCvis::MCMCpstr(t_fit, params = 'sigma', func = sd)[[1]]
    
    mn_alpha_beta <- MCMCvis::MCMCpstr(t_fit, params = 'alpha_beta', func = mean)[[1]]
    sd_alpha_beta <- MCMCvis::MCMCpstr(t_fit, params = 'alpha_beta', func = sd)[[1]]
    
    mn_beta_beta <- MCMCvis::MCMCpstr(t_fit, params = 'beta_beta', func = mean)[[1]]
    sd_beta_beta <- MCMCvis::MCMCpstr(t_fit, params = 'beta_beta', func = sd)[[1]]
    
    #diagnostics
    num_diverge <- rstan::get_num_divergent(t_fit)
    model_summary <- MCMCvis::MCMCsummary(t_fit, excl = 'y_rep', round = 2)
    max_rhat <- max(as.vector(model_summary[, grep('Rhat', colnames(model_summary))]))
    min_neff <- min(as.vector(model_summary[, grep('n.eff', colnames(model_summary))]))
    
    t_full <- data.frame(species = sp,
                         mn_mu_alpha,
                         sd_mu_alpha,
                         mn_sigma,
                         sd_sigma,
                         mn_alpha_beta,
                         sd_alpha_beta,
                         mn_beta_beta,
                         sd_beta_beta,
                         num_diverge,
                         max_rhat,
                         min_neff)
    
    out <- rbind(out, t_full)
    
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

saveRDS(out, file = paste0('arrival_master_', IAR_out_date, '.rds'))

print('I completed!')
