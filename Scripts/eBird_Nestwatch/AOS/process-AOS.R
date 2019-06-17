######################
# Extract AOS pheno trends model output
#
# mu_alpha - absolute arrival date for species (intercepts drawn from this mean)
# sigma - interannual variability in arrival date after accounting for trend
# alpha_beta - magnitude of phenological change over time
# beta_beta - effect of lat on magnitude of phenological change over time
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

in_dir <- '/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/trends_output_2019-06-16'
out_dir <- '/Users/caseyyoungflesh/Desktop/Bird_Phenology_Offline/Data/Processed/'

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
  print(sp)
  
  #if that species RDS object exists in dir
  #if (length(grep(paste0(sp, '-', in_date, '-pheno_trends_stan_output.rds'), list.files())) > 0)
  if (length(grep(paste0(sp, '-2019-05-26-pheno_trends_stan_output.rds'), list.files())) > 0)
  {
    #read in model output
    #t_fit <- readRDS(paste0(sp, '-', in_date, '-pheno_trends_stan_output.rds'))
    #t_in <- readRDS(paste0(sp, '-', in_date, ' -pheno_trends_stan_input.rds'))
    t_fit <- readRDS(paste0(sp, '-2019-05-26-pheno_trends_stan_output.rds'))
    #t_in <- readRDS(paste0(sp, '-2019-05-26-pheno_trends_stan_input.rds'))
    
    #make sure there are samples
    if (length(t_fit@par_dims) > 0)
    {
      #add year range and cell range (maybe range area calculated from range maps)?
      
      mn_mu_alpha <- MCMCvis::MCMCpstr(t_fit, params = 'mu_alpha', func = mean)[[1]]
      sd_mu_alpha <- MCMCvis::MCMCpstr(t_fit, params = 'mu_alpha', func = sd)[[1]]
      mu_alpha_ch <- MCMCvis::MCMCchains(t_fit, params = 'mu_alpha')
      colnames(mu_alpha_ch) <- sp
      
      mn_sigma <- MCMCvis::MCMCpstr(t_fit, params = 'sigma', func = mean)[[1]]
      sd_sigma <- MCMCvis::MCMCpstr(t_fit, params = 'sigma', func = sd)[[1]]
      sigma_ch <- MCMCvis::MCMCchains(t_fit, params = 'sigma')
      colnames(sigma_ch) <- sp
      
      mn_alpha_beta <- MCMCvis::MCMCpstr(t_fit, params = 'alpha_beta', func = mean)[[1]]
      sd_alpha_beta <- MCMCvis::MCMCpstr(t_fit, params = 'alpha_beta', func = sd)[[1]]
      alpha_beta_ch <- MCMCvis::MCMCchains(t_fit, params = 'alpha_beta')
      colnames(alpha_beta_ch) <- sp
      
      mn_beta_beta <- MCMCvis::MCMCpstr(t_fit, params = 'beta_beta', func = mean)[[1]]
      sd_beta_beta <- MCMCvis::MCMCpstr(t_fit, params = 'beta_beta', func = sd)[[1]]
      beta_beta_ch <- MCMCvis::MCMCchains(t_fit, params = 'beta_beta')
      colnames(beta_beta_ch) <- sp
      
      #n_cells <- length(t_in$NC)
      #n_years <- length(t_in$years)
      #rng_lat <- range(t_in$lat)[2] - range(t_in$lat)[1]
      
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
                           #n_cells,
                           #n_years,
                           #rng_lat,
                           num_diverge,
                           max_rhat,
                           min_neff)
      
      out <- rbind(out, t_full)
      
      if (i == 1)
      {
        mu_alpha_post <- mu_alpha_ch
        sigma_post <- sigma_ch
        alpha_beta_post <- alpha_beta_ch
        beta_beta_post <- beta_beta_ch
      } else {
        mu_alpha_post <- cbind(mu_alpha_post, mu_alpha_ch)
        sigma_post <- cbind(sigma_post, sigma_ch)
        alpha_beta_post <- cbind(alpha_beta_post, alpha_beta_ch)
        beta_beta_post <- cbind(beta_beta_post, beta_beta_ch)
      }
    } else {
      print('No samples for species')
    }
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop


setwd(out_dir)

#saveRDS(out, file = paste0('pheno_trends_master_', in_date, '.rds'))
#saveRDS(mu_alpha_post, file = paste0('pheno_trends_mu_alpha_post_', in_date, '.rds'))
#saveRDS(sigma_post, file = paste0('pheno_trends_sigma_post_', in_date, '.rds'))
#saveRDS(alpha_beta_post, file = paste0('pheno_trends_alpha_beta_post_', in_date, '.rds'))
#saveRDS(beta_beta_post, file = paste0('pheno_trends_beta_beta_post_', in_date, '.rds'))

print('I completed!')
