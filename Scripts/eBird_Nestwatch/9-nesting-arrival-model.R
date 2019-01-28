####################
# 9 - nesting date (from Nestwatch and others) ~ nesting date (from IAR model)
#
# *How well do arrival dates (as determined from the IAR model) predict nesting dates (as determined from Nestwatch)
####################


# data fusion model for arrival ~ breeding --------------------------------


#obs_EB_breeding_date is halfmax estimate for EB breeding date for that species/cell/year
#known_EB_uncertainty is uncertainty in halfmax estimate for EB breeding date for that species/cell/year
#obs_NW_breeding_date is mean NW breeding date for that species/cell/year
#known_NW_uncertainty is sd of NW breeding date for that species/cell/year
#obs_MAPS_breeding_date is mean of breeding date period for that species/cell/year
#known_MAPS_uncertainty is length of period (uniform distribution)

###observation models
#obs_IAR_arrival_date ~ N(true_IAR_arrival_date, known_IAR_uncertainty)
#obs_EB_breeding_date ~ N(true_EB_breeding_date, known_EB_uncertainty)
#obs_NW_breeding_date ~ N(true_NW_breeding_date, known_NW_uncertainty)
#obs_MAPS_breeding_date ~ N(true_MAPS_breeding_date, known_MAPS_uncertainty)


#ONE WAY (arrival ~ breeding):
###process model for explanatory var
#true_EB_breeding_date ~ N(mu_EB, sigma_EB)
#mu_EB = alpha_EB + beta1_EB * true_NW_breeding_date + beta2_EB * true_MAPS_breeding date

###process model for response var
#true_IAR_arrival_date ~ N(mu, sigma)
#mu = alpha + beta * true_EB_breeding_date


#ALTERNATIVELY (breeding ~ arrival):
#true_EB_breeding_date ~ N(master_breeding_date, sigma_EB)
#mu_EB = alpha_EB + beta1_EB * true_NW_breeding_date + beta2_EB * true_MAPS_breeding date

#master_breeding_date ~ N(mu, sigma)
#mu = alpha + beta * true_IAR_arrival_date


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)
library(ggplot2)
library(rstan)
library(MCMCvis)


# import IAR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

IAR_out_dir <- 'IAR_output_2018-11-12'
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)

IAR_data <- readRDS(paste0('master_arrival_', IAR_out_date, '.rds'))



# import IAR species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)

#remove underscore and coerce to vector
species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))


# create hex grid ---------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6) 



#hierarchical model to relate arrival date (from IAR model) to nesting date (from Nest Watch)

#one model per species






#############################
#Process data for model input
#############################






#create data list for Stan
DATA <- list(J = nyr,
             N = ncel, 
             N_obs = len_y_obs_in,
             N_mis = len_y_mis_in,
             N_edges = nrow(ninds), 
             node1 = ninds[,1],
             node2 = ninds[,2],
             y_obs = y_obs_in,
             sigma_y = sigma_y_in,
             scaling_factor = scaling_factor,
             ii_obs = ii_obs_in,
             ii_mis = ii_mis_in)



# Stan model --------------------------------------------------------------

#Might be alright - data needs to be formatted properly though
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

nest_arrival <- '
data {
int<lower = 0> N;                                     // number of cell/years
int<lower = 0> J;                                     // max number of obs in a cell/year for Nestwatch
real<lower = 0, upper = 200> obs_IAR[N];              // mean output IAR model
real sigma_IAR[N];                                    // sd output IAR model
real p_obs_NW[J,N];                                   // nestwatch data
int<lower = 0> ii_obs[J, N];                          // indices of observed data
int<lower = 0> ii_mis[J, N];                          // indices of missing data
int<lower = 0> N_obs[N];                              // number of non-missing for each cell/year
int<lower = 0> N_mis[N];                              // number missing for each cell/year
}


parameters {
real<lower = 0, upper = 200> p_obs_NW_mis[J, N];      // missing nestwatch data
real<lower = 0> true_IAR[N];                          //true arrival date
real<lower = 0> true_NW[N];                           //true nesting date
real<lower = 0> sigma_NW;                             //sd for cell/year averaging
real mu[N];
real<lower = 0> sigma;
real alpha;
real beta;
}

transformed parameters {
real<lower = 0, upper = 200> obs_NW[J, N];            // nestwatch data to be modeled

// indexing to avoid NAs
for (n in 1:N)
{
  obs_NW[ii_obs[1:N_obs[n], n], n] = p_obs_NW[1:N_obs[n], n];
  obs_NW[ii_mis[1:N_mis[n], n], n] = p_obs_NW_mis[1:N_mis[n], n];
}
}

model {

// observation model - modeling true state as a function of some observed state

obs_IAR ~ normal(true_IAR, sigma_IAR); //sigma_IAR is given (one for each observation)

for (n in 1:N) //n is cell/year number
{
  obs_NW[,n] ~ normal(true_NW[n], sigma_NW); //sigma_NW is estimated (one per species)
}

  true_NW ~ normal(mu, sigma); //one sigma, one mu for each cell/year
  mu = alpha + beta * true_IAR; //one alpha, one beta
}'




# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tt <- proc.time()
fit <- stan(model_code = nest_arrival,
            data = DATA,
            chains = 4,
            iter = 2,
            cores = 4,
            pars = c('true_IAR', 'true_NW', 'sigma_NW', 'alpha', 'beta', 'mu', 'sigma'),
            control = list(max_treedepth = 25, adapt_delta = 0.95, stepsize = 0.005)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60







#save to RDS
# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
# saveRDS(fit, file = paste0('IAR_stan_', args, '-', IAR_out_date, '.rds'))
# fit <- readRDS('IAR_stan_Catharus_minimus-2018-11-12.rds')




# diagnostics -------------------------------------------------------------


# pairs(fit, pars = c('sigma', 'rho'))

# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
# max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
# get_elapsed_time(fit)

# MCMCtrace(fit)
# MCMCsummary(fit, params = c('sigma', 'rho', 'beta0'), n.eff = TRUE)
# MCMCsummary(fit, params = c('theta', 'phi'), n.eff = TRUE)

# print(fit, pars = c('sigma', 'rho'))

# #shiny stan
# library(shinystan)
# launch_shinystan(fit)




# write model results to file ---------------------------------------------

# options(max.print = 50000)
# sink(paste0('IAR_results_', args, '.txt'))
# cat(paste0('IAR results ', args, ' \n'))
# cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
# print(fit)
# sink()








