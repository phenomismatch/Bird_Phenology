######################
# Interannual variance (with uncertainty) for each species/cell
# 
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/'



# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
juv_date <- '2019-11-07'

#IAR data
arr_date <- '2019-05-26'

#run date
run_date <- '2019-11-08'


# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS
#from Bird_Phenology repo
setwd(paste0(dir, 'pheno_trends/Data/juv_master_', juv_date))

#read in 
juvs_master <- readRDS(paste0('juv_master_', juv_date, '.rds'))

#master arrival data (from IAR output)
setwd(paste0(dir, 'pheno_trends/Data/arrival_master_', arr_date))

arr_master <- readRDS(paste0('arrival_master_', arr_date, '.rds'))

p_juv <- dplyr::filter(juvs_master, !is.na(juv_mean))
p_arr <- dplyr::filter(arr_master, !is.na(mean_pre_IAR))


sp_idx <- as.numeric(factor(p_juv$species))
cn_id <- as.numeric(p_juv$cell)

# Stan model --------------------------------------------------------------

DATA <- list(y = p_juv$juv_mean,
             sd_y = p_juv$juv_sd,
             N = NROW(p_juv),
             sp = sp_idx,
             cn_id = cn_id, 
             Nsp = length(unique(sp_idx)))


#interannual variance in arrival and fledging dates for each species/cell

stanmodel1 <- "
data {
int<lower=0> N;                      // number of data points
vector<lower=0>[N] y;                 // response
vector<lower=0>[N] sd_y;               // uncertainty in response
int<lower=0> sp[N];              
int<lower=0> Nsp;
real<lower=0> cell_lat[N];
real sc_tmax[N];
}

parameters {
real<lower = 0> sigma_raw;
vector[N] mu_y_raw;
real mu_alpha_raw;
real mu_beta_raw;
real mu_gamma_raw;
}

transformed parameters {
real<lower = 0> sigma;
vector[N] mu_y;
vector[N] mu;
vector<lower = 0>[3] sigma_sp;
real mu_alpha;
real mu_beta;
real mu_gamma;

mu_alpha = mu_alpha_raw * 200;
mu_beta = mu_beta_raw * 2;
mu_gamma = mu_gamma_raw * 2;
sigma = sigma_raw * 10;
mu_y = mu_y_raw * 40 + 60;


alpha = mu_alpha + abg[,1];
beta = mu_beta + abg[,2];
gamma = mu_gamma + abg[,3];

// implies mu_y[i] ~ normal(mu[i], sigma)
mu_y = mu_y_raw * sigma + mu; 

}

model {

// observation model for difference
y ~ normal(mu_y, sd_y);

}

"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.999
TREE_DEPTH <- 16
STEP_SIZE <- 0.0001
CHAINS <- 3
ITER <- 4000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma',
                            'abg',
                            'mu_alpha',
                            'mu_beta',
                            'mu_gamma',
                            'sigma_sp',
                            'Rho',
                            'sigma',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


