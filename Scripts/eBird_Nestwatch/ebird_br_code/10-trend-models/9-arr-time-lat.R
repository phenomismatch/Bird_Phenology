####################
# 9 - arrival data ~ year + lat
#
####################


# top-level dir --------------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)


# import ARR/BR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))


MODEL_DATE <- '2019-02-07'

#arrival data
IAR_out <- 'IAR_output_2019-01-16'
IAR_out_date <- substr(IAR_out, start = 12, stop = 21)

IAR_data <- readRDS(paste0('arrival_master_', IAR_out_date, '.rds'))



mdf <- dplyr::select(IAR_data, species, cell, year, cell_lat, cell_lon, 
                     mean_post_IAR, sd_post_IAR)



#OVERVIEW
#=======#
#arrival date ~ year + lat


#DETAILS
#======#
# #observation model - where obs and sigma are known
# x_obs[i] ~ normal(x_true[i], sigma_x[i])
# 
# #arrival as a function of year and latitude
# x_true[i] ~ normal(mu_arr[i], sigma_arr)
# mu_arr[i] = alpha_arr_j + beta1_arr_j * year[i] + beta2_arr_j * lat[i] + beta3_j * year[i] * lat[i]



#do not merge with breeding data
mdf3 <- mdf

#number codes for species
sp_num <- as.numeric(factor(mdf3$species))

#number codes for years
yr_num <- as.numeric(factor(mdf3$year))

#create data list for Stan
DATA <- list(x_obs = mdf3$mean_post_IAR,
             sigma_x = mdf3$sd_post_IAR,
             sp_id = sp_num,
             year = yr_num,
             lat = mdf3$cell_lat,
             US = length(unique(sp_num)),
             N = NROW(mdf3))



# Stan model --------------------------------------------------------------

arr_time_lat <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
real<lower = 0, upper = 200> x_obs[N];                // mean halfmax IAR
real<lower = 0> sigma_x[N];                           // sd halfmax IAR
int<lower = 1, upper = US> sp_id[N];                  // species ids
vector<lower = 1, upper = 17>[N] year;
vector<lower = 26, upper = 90>[N] lat;
}

parameters {
real<lower = 0, upper = 200> x_true[N];                           //true arrival
real mu_alpha_raw;
real mu_beta1_raw;
real mu_beta2_raw;
real mu_beta3_raw;
real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_beta1_raw;
real<lower = 0> sigma_beta2_raw;
real<lower = 0> sigma_beta3_raw;
real<lower = 0> sigma_raw;
real alpha_raw[US];
real beta1_raw[US];
real beta2_raw[US];
real beta3_raw[US];
}

transformed parameters {
real mu_alpha;
real mu_beta1;
real mu_beta2;
real mu_beta3;
real sigma_alpha;
real sigma_beta1;
real sigma_beta2;
real sigma_beta3;
real sigma;
real alpha[US];
real beta1[US];
real beta2[US];
real beta3[US];
real mu[N];

// non-centered parameterization
mu_alpha = mu_alpha_raw * 20 + 70;                       // implies mu_alpha ~ normal(70, 20)
mu_beta1 = mu_beta1_raw * 2 + 1;                           // implies mu_beta ~ normal(1, 2)
mu_beta2 = mu_beta2_raw * 2 + 1;
mu_beta3 = mu_beta3_raw * 2 + 1;
sigma_alpha = sigma_alpha_raw * 10;                      // implies sigma_alpha ~ halfnormal(0, 10)
sigma_beta1 = sigma_beta1_raw * 3;                         // implies sigma_beta ~ halfnormal(0, 3)
sigma_beta2 = sigma_beta2_raw * 3;
sigma_beta3 = sigma_beta3_raw * 3;
sigma = sigma_raw * 10;                                  // implies sigma ~ halfnormal(0, 10)

for (j in 1:US)
{
  alpha[j] = alpha_raw[j] * sigma_alpha + mu_alpha;      // implies alpha[j] ~ normal(mu_alpha, sigma_alpha)
  beta1[j] = beta1_raw[j] * sigma_beta1 + mu_beta1;          // implies beta[j] ~ normal(mu_beta, sigma_beta)
  beta2[j] = beta2_raw[j] * sigma_beta2 + mu_beta2;
  beta3[j] = beta3_raw[j] * sigma_beta3 + mu_beta3;
}

for (i in 1:N)
{
  mu[i] = alpha[sp_id[i]] + beta1[sp_id[i]] * year[i] + beta2[sp_id[i]] * lat[i] + beta3[sp_id[i]] * year[i] * lat[i];
}
}

model {

// observation model - modeling true state as a function of some observed state
x_obs ~ normal(x_true, sigma_x);

// non-centered parameterization
mu_alpha_raw ~ normal(0, 1);
mu_beta1_raw ~ normal(0, 1);
mu_beta2_raw ~ normal(0, 1);
mu_beta3_raw ~ normal(0, 1);
sigma_alpha_raw ~ normal(0, 1);
sigma_beta1_raw ~ normal(0, 1);
sigma_beta2_raw ~ normal(0, 1);
sigma_beta3_raw ~ normal(0, 1);
sigma_raw ~ normal(0, 1);

for (j in 1:US)
{
  alpha_raw[j] ~ normal(0, 1);
  beta1_raw[j] ~ normal(0, 1);
  beta2_raw[j] ~ normal(0, 1);
  beta3_raw[j] ~ normal(0, 1);
}

x_true ~ normal(mu, sigma);
}

generated quantities {
// real y_rep[N];
// real errors[N];
// real PPC_mean;
// real BR2;
// real arr_br[N];
// #traditional R^2
// RSS = dot_self(y - mu);
// TSS = dot_self(y - mean(y));
// R2 = 1 - RSS/TSS;
// #new Bayes R^2 - http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
// errors = y - mu;
// BR2 = var(mu)/(var(mu) + var(errors));
// #PPC
// y_rep = normal_rng(mu, sigma);
// PPC_mean = mean(y_rep)
// #arrival date - breeding date
// arr_br = x_true - y_true
}
'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


DELTA <- 0.97
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit <- rstan::stan(model_code = arr_time_lat,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 'beta1', 'beta2', 'beta3',
                            'mu_alpha', 'mu_beta1', 'mu_beta2', 'mu_beta3',
                            'sigma_alpha', 'sigma_beta1', 'sigma_beta2', 'sigma_beta3',
                            'sigma', 'x_true'),
                   control = list(max_treedepth = TREE_DEPTH, adapt_delta = DELTA, stepsize = STEP_SIZE)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit, file = paste0('temp_ARR_YEAR_LAT_stan_', MODEL_DATE, '.rds'))
#fit <- readRDS(paste0('temp_ARR_YEAR_LAT_stan_', MODEL_DATE, '.rds'))




# diagnostics -------------------------------------------------------------


# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = c('alpha', 'beta'), ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'sigma', ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'mu', ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'x_true')
# MCMCvis::MCMCplot(fit, params = 'beta1', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'beta2', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'beta3', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'mu', ISB = FALSE, rank = TRUE)
# #MCMCtrace(fit)
# 
# (num_diverge <- rstan::get_num_divergent(fit))
# (num_tree <- rstan::get_num_max_treedepth(fit))
# (num_BFMI <- rstan::get_low_bfmi_chains(fit))


#shiny stan
# library(shinystan)
# launch_shinystan(fit)


#plot results - true states with sd error bars


# write model results to file ---------------------------------------------

# options(max.print = 50000)
# sink(paste0('BR_ARR_results_', MODEL_DATE, '.txt'))
# cat(paste0('BR_ARR results ', MODEL_DATE, ' \n'))
# cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
# cat(paste0('Adapt delta: ', DELTA, ' \n'))
# cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
# #cat(paste0('Step size: ', STEP_SIZE, ' \n'))
# cat(paste0('Number of divergences: ', num_diverge, ' \n'))
# cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
# cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
# print(fit)
# sink()


