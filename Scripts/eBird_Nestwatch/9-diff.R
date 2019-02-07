
# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)


# import ARR/BR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))


MODEL_DATE <- '2019-02-04'

#arrival data
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
                     mean_post_IAR, sd_post_IAR, EB_HM_mean, EB_HM_sd, 
                     NW_mean_cid, NW_sd_cid, MAPS_l_bounds, MAPS_u_bounds, d_avail)




# nesting date ~ arrival date ---------------------------------------------

#OVERVIEW
#=======#
# diff ~ arrival + year + lat


#DETAILS
#======#





#remove na vals for BR cods
to.rm <- which(is.na(mdf$EB_HM_mean))
mdf2 <- mdf[-to.rm, ]

#remove species where there are less than 3 obs
obs_per_sp <- plyr::count(mdf2, 'species')
sp.to.rm <- obs_per_sp[which(obs_per_sp[,2] < 3), 1]
mdf3 <- mdf2[-which(mdf2$species %in% sp.to.rm),]

#number codes for species
sp_num <- as.numeric(factor(mdf3$species))
#number codes for years
yr_num <- as.numeric(factor(mdf3$year))
range(yr_num)



arr_br_mean <- MCMCpstr(fit, params = 'arr_br', func = mean)[[1]]
arr_br_sd <- MCMCpstr(fit, params = 'arr_br', func = sd)[[1]]
x_true_mean <- MCMCpstr(fit, params = 'x_true', func = mean)[[1]]
x_true_sd <- MCMCpstr(fit, params = 'x_true', func = sd)[[1]]

#create data list for Stan
DATA <- list(y_obs = arr_br_mean,
             sigma_y = arr_br_sd,
             x_obs = x_true_mean,
             sigma_x = x_true_sd,
             year = mdf3$year,
             sp_id = sp_num,
             US = length(unique(sp_num)),
             N = NROW(arr_br_mean))



# Stan model --------------------------------------------------------------

diff_arr_br <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
real y_obs[N];                // mean halfmax BR codes
real<lower = 0> sigma_y[N];                           // sd halfmax BR codes
real<lower = 0, upper = 200> x_obs[N];                // mean halfmax IAR
real<lower = 0> sigma_x[N];                           // sd halfmax IAR
int<lower = 1, upper = US> sp_id[N];                  // species ids
int<lower = 2002, upper = 2017> year[N];
}

parameters {
real y_true[N];                           //true arrival date
real<lower = 0, upper = 200> x_true[N];                           //true nesting date
real mu_alpha_raw;
real mu_beta1_raw;
real mu_beta2_raw;
real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_beta1_raw;
real<lower = 0> sigma_beta2_raw;
real<lower = 0> sigma_raw;
real alpha_raw[US];
real beta1_raw[US];
real beta2_raw[US];
}

transformed parameters {

real mu_alpha;
real mu_beta1;
real mu_beta2;
real sigma_alpha;
real sigma_beta1;
real sigma_beta2;
real sigma;
real alpha[US];
real beta1[US];
real beta2[US];
real mu[N];

// non-centered parameterization

mu_alpha = mu_alpha_raw * 20 + 70;                       // implies mu_alpha ~ normal(70, 20)
mu_beta1 = mu_beta1_raw * 2 + 1;                           // implies mu_beta ~ normal(1, 2)
mu_beta2 = mu_beta2_raw * 2 + 1;                           // implies mu_beta ~ normal(1, 2)
sigma_alpha = sigma_alpha_raw * 10;                      // implies sigma_alpha ~ halfnormal(0, 10)
sigma_beta1 = sigma_beta1_raw * 3;                         // implies sigma_beta ~ halfnormal(0, 3)
sigma_beta2 = sigma_beta2_raw * 3;                         // implies sigma_beta ~ halfnormal(0, 3)
sigma = sigma_raw * 10;                                  // implies sigma ~ halfnormal(0, 10)

for (j in 1:US)
{
  alpha[j] = alpha_raw[j] * sigma_alpha + mu_alpha;      // implies alpha[j] ~ normal(mu_alpha, sigma_alpha)
  beta1[j] = beta1_raw[j] * sigma_beta1 + mu_beta1;          // implies beta[j] ~ normal(mu_beta, sigma_beta)
  beta2[j] = beta2_raw[j] * sigma_beta2 + mu_beta2;          // implies beta[j] ~ normal(mu_beta, sigma_beta)
}

for (i in 1:N)
{
  mu[i] = alpha[sp_id[i]] + beta1[sp_id[i]] * x_true[i] + beta2[sp_id[i]] * year[i];
}

}

model {

// observation model - modeling true state as a function of some observed state

y_obs ~ normal(y_true, sigma_y);
x_obs ~ normal(x_true, sigma_x);

// non-centered parameterization

mu_alpha_raw ~ normal(0, 1);
mu_beta1_raw ~ normal(0, 1);
mu_beta2_raw ~ normal(0, 1);
sigma_alpha_raw ~ normal(0, 1);
sigma_beta1_raw ~ normal(0, 1);
sigma_beta2_raw ~ normal(0, 1);
sigma_raw ~ normal(0, 1);

for (j in 1:US)
{
  alpha_raw[j] ~ normal(0, 1);
  beta1_raw[j] ~ normal(0, 1);
  beta2_raw[j] ~ normal(0, 1);
}

y_true ~ normal(mu, sigma);

}

generated quantities {

//vector[N] resids;
//real BR2;
//vector[N] y_rep;
//real PPC_mean;
//real var_mu;
//real var_resids;

// #traditional R^2
// RSS = dot_self(y - mu);
// TSS = dot_self(y - mean(y));
// R2 = 1 - RSS/TSS;

// new Bayes R^2 - http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
// calculate residuals variance for mu and errors to use in R^2
//resids = y_true - mu;
//var_mu = (dot_self(mu - mean(mu))) / (N - 1);
//var_resids = (dot_self(resids - mean(resids))) / (N - 1);
//BR2 = var_mu/(var_mu + var_resids);

// PPC
//for (i in 1:N)
//{
//  y_rep[N] = normal_rng(mu[N], sigma);
//}
//PPC_mean = mean(y_rep);
}
'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


DELTA <- 0.97
TREE_DEPTH <- 16
STEP_SIZE <- 0.0005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit2 <- rstan::stan(model_code = diff_arr_br,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha', 'beta1', 'beta2', 'mu_alpha', 'mu_beta1', 'mu_beta2', 'sigma_alpha', 
                             'sigma_beta1', 'sigma_beta2', 'sigma', 'y_true', 'x_true'),
                    control = list(max_treedepth = TREE_DEPTH, adapt_delta = DELTA, stepsize = STEP_SIZE)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit2, file = paste0('temp_ARRBR_YEAR_stan_', MODEL_DATE, '.rds'))
#fit2 <- readRDS(paste0('temp_ARRBR_YEAR_stan_', MODEL_DATE, '.rds'))




# diagnostics -------------------------------------------------------------


MCMCvis::MCMCsummary(fit2, n.eff = TRUE, params = c('alpha', 'beta1', 'beta2'))
MCMCvis::MCMCsummary(fit2, n.eff = TRUE, params = 'sigma')
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'y_true')
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'x_true')
#MCMCtrace(fit)

(num_diverge <- rstan::get_num_divergent(fit2))
(num_tree <- rstan::get_num_max_treedepth(fit2))
(num_BFMI <- rstan::get_low_bfmi_chains(fit2))


#shiny stan
# library(shinystan)
# launch_shinystan(fit2)
