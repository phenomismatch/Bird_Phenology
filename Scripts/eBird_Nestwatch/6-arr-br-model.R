######################
# 6 - arrival ~ young hitting nets
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# model dir ------------------------------------------------------------

#juveniles MAPS
arr_br_dir <- 'arr_br_2019-08-19'

#IAR data
arr_master_dir <- 'arrival_master_2019-05-26'


# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Process data ------------------------------------------------------------------

#juveniles hitting nets - MAPS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_br_dir))

juv_MAPS <- readRDS('arr_br_in.rds')

#IAR arrival data
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_master_dir))

arr_master <- readRDS(paste0(arr_master_dir, '.rds'))

#save only species/cell/years that match arr_master





# Stan model --------------------------------------------------------------

DATA <- list(N = NROW(data),
             Nsp = length(usp),
             y = response,
             sp = as.numeric(stan_data$sp_f), #species
             year = stan_data$year_f,
             lat = stan_data$lat,
             elev = stan_data$elev)


stanmodel3 <- "
data {
int<lower=0> N;                     // number of obs
int<lower=0> Nsp;                   // number of species
real<lower=0> y[N];                 // response
int<lower=1, upper=Nsp> sp[N];      // species ids
vector<lower=0>[N] year;
vector<lower=0>[N] lat;
vector<lower=0>[N] elev;
}

parameters {
real mu_alpha_raw;
real mu_beta_raw;
real mu_gamma_raw;
real mu_theta_raw;
vector<lower = 0>[Nsp] sigma_raw;
vector<lower = 0>[4] sigma_sp_raw;
cholesky_factor_corr[4] L_Rho;             // correlation matrix
matrix[4, Nsp] z;                          // z-scores
}

transformed parameters {
vector[N] mu;
matrix[Nsp, 4] abgt;                              // matrix for alpha, beta, gamma, and theta
matrix[4, 4] Rho;                                 // covariance matrix
vector[Nsp] alpha;
vector[Nsp] beta;
vector[Nsp] gamma;
vector[Nsp] theta;
vector[Nsp] alpha_c;
vector[Nsp] beta_c;
vector[Nsp] gamma_c;
vector[Nsp] theta_c;

vector<lower = 0>[Nsp] sigma;
vector<lower = 0>[4] sigma_sp;
real mu_alpha;
real mu_beta;
real mu_gamma;
real mu_theta;

sigma = sigma_raw * 5;
mu_alpha = mu_alpha_raw * 10 + 70;
mu_beta = mu_beta_raw * 0.1;
mu_gamma = mu_gamma_raw * 1;
mu_theta = mu_theta_raw * 0.01;
sigma_sp[1] = sigma_sp_raw[1] * 20;             // variance alpha
sigma_sp[2] = sigma_sp_raw[2] * 0.1;              // variance beta
sigma_sp[3] = sigma_sp_raw[3] * 1;              // variance gamma
sigma_sp[4] = sigma_sp_raw[4] * 0.01;              // variance gamma


// cholesky factor of covariance matrix multiplied by z score
abgt = (diag_pre_multiply(sigma_sp, L_Rho) * z)';
alpha = abgt[,1];
beta = abgt[,2];
gamma = abgt[,3];
theta = abgt[,4];
Rho = L_Rho * L_Rho';

for (i in 1:N)
{
  mu[i] = (mu_alpha + abgt[sp[i], 1]) + 
  (mu_beta + abgt[sp[i], 2]) * year[i] + 
  (mu_gamma + abgt[sp[i], 3]) * lat[i] +
  (mu_theta + abgt[sp[i], 4]) * elev[i];
}

alpha_c = mu_alpha + alpha;
beta_c = mu_beta + beta;
gamma_c = mu_gamma + gamma;
theta_c = mu_theta + theta;
}

model {
sigma_raw ~ std_normal();
mu_alpha_raw ~ std_normal();
mu_beta_raw ~ std_normal();
mu_gamma_raw ~ std_normal();
mu_theta_raw ~ std_normal();
sigma_sp_raw ~ std_normal();

to_vector(z) ~ std_normal();         // faster than normal(0, 1);
L_Rho ~ lkj_corr_cholesky(2);

y ~ normal(mu, sigma[sp]);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu, sigma[sp]);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.85
TREE_DEPTH <- 14
STEP_SIZE <- 0.0025
CHAINS <- 4
ITER <- 5000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel3,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma',
                            'theta',
                            'alpha_c',
                            'beta_c',
                            'gamma_c',
                            'theta_c',
                            'mu_alpha',
                            'mu_beta',
                            'mu_gamma',
                            'mu_theta',
                            'Rho',
                            'sigma_sp',
                            'sigma',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'wing_chord_changes/Results'))
saveRDS(fit, file = paste0('wc-tle-stan-output-', DATE, '.rds'))
#fit <- readRDS('wc-tle-stan-output-2019-05-02.rds')


num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))



# Calc diagnostics ---------------------------------------------------

#fit <- readRDS('Ictinia_mississippiensis-2019-05-26-pheno_trends_stan_output.rds')
# library(shinystan)
# launch_shinystan(fit)

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mn_stepsize <- sapply(sampler_params, 
                      function(x) mean(x[, 'stepsize__']))
mn_treedepth <- sapply(sampler_params, 
                       function(x) mean(x[, 'treedepth__']))
accept_stat <- sapply(sampler_params, 
                      function(x) mean(x[, 'accept_stat__']))



# Summaries ---------------------------------------------------------------

#get summary of model output
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, round = 2, excl = 'y_rep')

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])

#y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')

# bayesplot::ppc_stat(DATA$y_obs, y_rep, stat = 'mean')
# bayesplot::ppc_dens_overlay(DATA$y_obs, y_rep[1:500,])



# write model results to file ---------------------------------------------

options(max.print = 50000)
sink(paste0('wc-tle-stan-results-', DATE, '.txt'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
cat(paste0('Mean stepsize: ', round(mean(mn_stepsize), 5), ' \n'))
cat(paste0('Mean treedepth: ', round(mean(mn_treedepth), 1), ' \n'))
cat(paste0('Mean accept stat: ', round(mean(accept_stat), 2), ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output, na.rm = TRUE), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output, na.rm = TRUE), ' \n'))
print(model_summary)
sink()




# check output ------------------------------------------------------------


# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'alpha')
# MCMCvis::MCMCplot(fit, params = c('alpha'), labels = usp, main = 'alpha', rank = TRUE)
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'beta')
# MCMCvis::MCMCplot(fit, params = c('beta_c'), labels = usp,
#                   main = 'beta - change over time')
#  
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'gamma')
# MCMCvis::MCMCplot(fit, params = c('gamma_c'),
#                   labels = usp, main = 'gamma - change over latitude')
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'theta')
# MCMCvis::MCMCplot(fit, params = c('theta_c'),
#                   labels = usp, main = 'theta - change over elevation')
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'sigma')
# 
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 5, params = c('mu_beta', 'mu_gamma', 'mu_theta'), ISB = FALSE)
# MCMCvis::MCMCplot(fit, params = 'mu', ISB = FALSE)

# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# sapply(sampler_params, function(x) mean(x[, 'stepsize__']))
# sapply(sampler_params, function(x) mean(x[, 'divergent__']))
# sapply(sampler_params, function(x) mean(x[, 'treedepth__']))
# sapply(sampler_params, function(x) mean(x[, 'accept_stat__']))
# sapply(sampler_params, function(x) mean(x[, 'n_leapfrog__']))
# sapply(sampler_params, function(x) mean(x[, 'energy__']))

# library(shinystan)
# launch_shinystan(fit)


# # PPC ---------------------------------------------------------------------
# 
# y_val <- DATA$y
# y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
# bayesplot::ppc_dens_overlay(DATA$y, y_rep[1:100,])
# plot(DATA$y, y_rep[1,], pch = '.', xlim = c(0, 200), ylim = c(0, 200))
# abline(a = 0, b = 1, col = 'red', lty = 2)
# 
# 
# 
# # PPO ---------------------------------------------------------------------
# 
# sigma = sigma_raw * 5;
# mu_alpha = mu_alpha_raw * 10 + 70;
# mu_beta = mu_beta_raw * 0.1;
# mu_gamma = mu_gamma_raw * 1;
# mu_theta = mu_theta_raw * 0.01;
# sigma_sp[1] = sigma_sp_raw[1] * 20;             // variance alpha
# sigma_sp[2] = sigma_sp_raw[2] * 0.1;              // variance beta
# sigma_sp[3] = sigma_sp_raw[3] * 1;              // variance gamma
# sigma_sp[4] = sigma_sp_raw[4] * 0.01;              // variance gamma
# 
# 
# 
# #mu_alpha ~ N(70, 10)
# PR <- rnorm(10000, 70, 10)
# MCMCvis::MCMCtrace(fit,
#                    params = 'mu_alpha',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #mu_beta ~ N(0, 2)
# PR <- rnorm(10000, 0, 0.1)
# MCMCvis::MCMCtrace(fit,
#                    params = 'mu_beta',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #mu_gamma ~ N(0, 2)
# PR <- rnorm(10000, 0, 1)
# MCMCvis::MCMCtrace(fit,
#                    params = 'mu_gamma',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #mu_theta ~ N(0, 1)
# PR <- rnorm(10000, 0, 0.01)
# MCMCvis::MCMCtrace(fit,
#                    params = 'mu_theta',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #sigma_sp[1] ~ HN(0, 20)
# PR_p <- rnorm(10000, 0, 20)
# PR <- PR_p[which(PR_p > 0)]
# MCMCvis::MCMCtrace(fit,
#                    params = 'sigma_sp\\[1',
#                    ISB = 'FALSE',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #sigma_sp[2] ~ HN(0, 0.1)
# PR_p <- rnorm(10000, 0, 0.1)
# PR <- PR_p[which(PR_p > 0)]
# MCMCvis::MCMCtrace(fit,
#                    params = 'sigma_sp\\[2',
#                    ISB = 'FALSE',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #sigma_sp[3] ~ HN(0, 1)
# PR_p <- rnorm(10000, 0, 1)
# PR <- PR_p[which(PR_p > 0)]
# MCMCvis::MCMCtrace(fit,
#                    params = 'sigma_sp\\[3',
#                    ISB = 'FALSE',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #sigma_sp[4] ~ HN(0, 0.01)
# PR_p <- rnorm(10000, 0, 0.01)
# PR <- PR_p[which(PR_p > 0)]
# MCMCvis::MCMCtrace(fit,
#                    params = 'sigma_sp\\[4',
#                    ISB = 'FALSE',
#                    priors = PR,
#                    pdf = FALSE)
# 
# #sigma ~ HN(0, 5)
# PR_p <- rnorm(10000, 0, 5)
# PR <- PR_p[which(PR_p > 0)]
# MCMCvis::MCMCtrace(fit,
#                    params = 'sigma',
#                    priors = PR,
#                    pdf = FALSE)
 




