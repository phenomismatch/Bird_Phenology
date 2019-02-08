####################
# 9 - br data ~ year + lat
# 
# Individual species
#
####################


# top-level dir --------------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'

MODEL_DATE <- '2019-02-08'
BR_TIME_LAT_IND_DIR <- paste0('trend_models_', MODEL_DATE)


# species arg -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Vireo_olivaceus')


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)


# import ARR/BR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))


#arrival data - for cell lat/lon
IAR_out <- 'IAR_output_2019-01-16'
IAR_out_date <- substr(IAR_out, start = 12, stop = 21)

IAR_data <- readRDS(paste0('arrival_master_', IAR_out_date, '.rds'))


#breeding data
BR_out <- 'halfmax_breeding_2019-01-30'
BR_out_date <- substr(BR_out, start = 18, stop = 27)

BR_data <- readRDS(paste0('temp_breeding_master_', BR_out_date, '.rds'))


# merge datasets ----------------------------------------------------------

#left_join for now - could use full_join to include all BR codes where ARR is missing

#same naems for cols
colnames(BR_data)[1:3] <- c('species', 'cell', 'year')

#merge data sets - keeping all rows from IAR, matching BR rows
master_data <- dplyr::left_join(IAR_data, BR_data, by = c('species', 'cell', 'year'))


mdf <- dplyr::select(master_data, species, cell, year, cell_lat, cell_lon, 
                     EB_HM_mean, EB_HM_sd)


#OVERVIEW
#=======#
#arrival ~ N(alpha_j + beta_j * year)
#beta_j = alpha2 + beta2 * lat


#DETAILS
#======#
# #observation model - where obs and sigma are known
# x_obs[i] ~ normal(x_true[i], sigma_x[i])
# 
# #arrival as a function of year and latitude
# x_true[i] ~ normal(mu_arr[i], sigma_arr)
# mu_arr[i] = alpha[cn[i]] + beta[cn[i]] * year[i]
# for (j in J)
#   beta1[j] = alpha2 + beta2 * lat[j]


#do not merge with breeding data
mdf3 <- filter(mdf, species == args)

#number codes for species
#sp_num <- as.numeric(factor(mdf3$species))

#number code for cells (not actual cell number)
cell_num <- as.numeric(factor(mdf3$cell))

#lats of cells
cell_mrg <- data.frame(cell = mdf3$cell, cell_num = cell_num, 
                       cell_lat = mdf3$cell_lat)

u_cell_mrg <- unique(cell_mrg)

#number codes for years
yr_num <- as.numeric(factor(mdf3$year))



# create Stan data object -------------------------------------------------

ncell <- NROW(u_cell_mrg)
nyr <- length(unique(yr_num))

#indices for observed data
ii_obs_in <- which(!is.na(mdf3$EB_HM_mean))
ii_mis_in <- which(is.na(mdf3$EB_HM_mean))

#number of observation and NAs for each year
len_y_obs_in <- length(ii_obs_in)
len_y_mis_in <- length(ii_mis_in)

y_obs_in <- mdf3$EB_HM_mean[ii_obs_in]
sigma_y_in <- mdf3$EB_HM_sd

#fill 0.1 for sd - does not effect value for y_true - just need a place holder
sigma_y_in[which(is.na(sigma_y_in), arr.ind = TRUE)] <- 0.1


#create data list for Stan
DATA <- list(y_obs = y_obs_in,
             sigma_y = sigma_y_in,
             cn_id = cell_num,
             year = yr_num,
             lat = u_cell_mrg$cell_lat,
             US = ncell,
             N = NROW(mdf3),
             N_obs = len_y_obs_in,
             N_mis = len_y_mis_in,
             ii_obs = ii_obs_in,
             ii_mis = ii_mis_in)



# Stan model --------------------------------------------------------------

br_time_lat_ind <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
int<lower = 0> N_obs;
int<lower = 0> N_mis;
int<lower = 1> ii_obs[N_obs];
int<lower = 1> ii_mis[N_mis];
real<lower = 0, upper = 300> y_obs[N_obs];                // mean halfmax IAR
real<lower = 0> sigma_y[N];                           // sd halfmax IAR
int<lower = 1, upper = US> cn_id[N];                  // species ids
vector<lower = 1, upper = 17>[N] year;
vector<lower = 26, upper = 90>[US] lat;
}

parameters {
real<lower = 0, upper = 300> y_mis[N_mis];
real<lower = 0, upper = 300> y_true[N];                           //true arrival
real mu_alpha_raw;
real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_raw;
real alpha_raw[US];
real alpha2_raw;
real beta2_raw;
}

transformed parameters {
real<lower = 0, upper = 300> y[N];
real mu_alpha;
real sigma_alpha;
real sigma;
real alpha[US];
real beta[US];
real alpha2;
real beta2;
real mu[N];

// non-centered parameterization
mu_alpha = mu_alpha_raw * 20 + 70;                       // implies mu_alpha ~ normal(70, 20)
sigma_alpha = sigma_alpha_raw * 10;                      // implies sigma_alpha ~ halfnormal(0, 10)
sigma = sigma_raw * 10;                                  // implies sigma ~ halfnormal(0, 10)
alpha2 = alpha2_raw * 10;                                 // implies alpha2 ~ normal(0, 10)
beta2 = beta2_raw * 2 + 1;                               // implies beta2 ~ normal(1, 2)

for (j in 1:US)
{
  alpha[j] = alpha_raw[j] * sigma_alpha + mu_alpha;      // implies alpha[j] ~ normal(mu_alpha, sigma_alpha)
  beta[j] = alpha2 + beta2 * lat[j];
}

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
}

// indexing to avoid NAs
y[ii_obs] = y_obs;
y[ii_mis] = y_mis;

}

model {

// observation model - modeling true state as a function of some observed state
y ~ normal(y_true, sigma_y);

// non-centered parameterization
mu_alpha_raw ~ normal(0, 1);
sigma_alpha_raw ~ normal(0, 1);
sigma_raw ~ normal(0, 1);
alpha2_raw ~ normal(0, 1);
beta2_raw ~ normal(0, 1);

for (j in 1:US)
{
  alpha_raw[j] ~ normal(0, 1);
}

y_true ~ normal(mu, sigma);
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
CHAINS <- 1
ITER <- 30

tt <- proc.time()
fit <- rstan::stan(model_code = br_time_lat_ind,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 'beta', 'mu_alpha', 
                            'alpha2', 'beta2',
                            'sigma_alpha', 'sigma', 'y_true'),
                   control = list(adapt_delta = DELTA, 
                                  max_treedepth = TREE_DEPTH, 
                                  stepsize = STEP_SIZE))
run_time <- (proc.time() - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', BR_TIME_LAT_IND_DIR))
saveRDS(fit, file = paste0('temp_BR_YEAR_LAT_IND_stan_', MODEL_DATE, '_', args, '.rds'))
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

options(max.print = 50000)
sink(paste0('BR_YEAR_LAT_IND_results_', MODEL_DATE, '_', args, '.txt'))
cat(paste0('BR_YEAR_LAT_IND_results_', MODEL_DATE, '_', args, ' \n'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
print(fit)
sink()



