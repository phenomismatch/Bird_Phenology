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



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_br_dir))

juv_MAPS <- readRDS('arr_br_in.rds')

#add undercore to species names to match format of arrival master
juv_MAPS$species <- gsub(' ', '_', juv_MAPS$species)

#master arrival data (from IAR output)
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_master_dir))

arr_master <- readRDS(paste0(arr_master_dir, '.rds'))

#filter for cells that we have data for pre-IAR
arr_master2 <- dplyr::filter(arr_master, !is.na(mean_pre_IAR))


#save only species/cell/years that match arr_master - to be used for juv
juv_data <- dplyr::inner_join(juv_MAPS, arr_master2, by = c('species', 'cell', 'year'))

#filter species based on number of cells of data like in AOS analyses?
#dplyr::filter()

#arrival data for unique species/cells/years - to be used for y and sd_y
arr_data <- unique(arr_data[,c('species', 'cell', 'year', 'mean_post_IAR', 'sd_post_IAR')])



# Reformat data -----------------------------------------------------------

#unqiue species in processed data
usp_f <- unique(arr_data$species)
nusp_f <- length(usp_f)
nuyr_f <- length(unique(arr_data$year))
nuc_f <- length(unique(arr_data$cell))
nind_f <- length(unique(juv_data$band_id))

#s,i,j
mn_arr_array <- array(NA, dim = c(nusp_f, nuc_f, nuyr_f))
sd_arr_array <- array(NA, dim = c(nusp_f, nuc_f, nuyr_f))
#s,i,j,k
juv_array <- array(NA, dim = c(nusp_f, nuc_f, nuyr_f, nind_f))

for (s in 1:length(usp_f))
{
  #SPECIES
  #s <- 1
  temp <- dplyr::filter(juv_data, species == usp_f[s])

  #CELL  
  uc_t <- unique(temp$cell)
  for (i in 1:length(uc_t))
  {
    #i <- 1
    temp2 <- dplyr::filter(temp, cell == uc_t[i])
    
    cell_counter <- 1
    #YEAR
    uyr_t <- unique(temp2$year)
    for (j in 1:length(uyr_t))
    {
      #j <- 1
      temp3 <- dplyr::filter(temp2, year == uyr_t[j])
      
      mn_arr_array[s,cell_counter,j] <- temp3$mean_post_IAR[1]
      sd_arr_array[s,i,j] <- temp3$sd_post_IAR[1]
      
      #INDIVIDUAL
      uind_t <- unique(temp3$band_id)
      for (k in 1:length(uind_t))
      {
        #k <- 1
        temp4 <- dplyr::filter(temp3, band_id == uind_t[k])
        
        juv_array[s,i,j,k] <- temp4$jday
      }
    }
  }
}

arr_data$mean_post_IAR #[s,i,j]
arr_data$sd_post_IAR #[s,i,j]


# Stan model --------------------------------------------------------------

DATA <- list(y = arr_data$mean_post_IAR, #[s,i,j]
             sd_y = arr_data$sd_post_IAR, #[s,i,j]
             juv = juv_data, #[s,i,j,k]
             NS = XX,
             NY = XX, #[s]
             NJ = XX, #[s,i]
             NK = XX, #[s,i,j]
             MY = XX, #max # years
             MJ = XX, #max # cells
             MK = XX) #max #ind


stanmodel1 <- "
data {
int<lower=0> NS;                      // number of species
int<lower=0> MY;                      // max # years
int<lower=0> MJ;                      // max # cells
int<lower=0> MK;                      // max # ind
int<lower=0> NY[NS];                  // array - number of years for each species
int<lower=0> NJ[NS,MY];               // array - number of cells for each species/year
int<lower=0> NK[NS,MY,MJ];            // array - number of ind for each species/year/cell
real<lower=0> y[NS, MY, MJ, MK];      // response
real<lower=0> sd_y[NS, MY, MJ, MK];   // uncertainty in response
real<lower=0> juv[NS, MY, MJ, MK];    // response
}

parameters {
real<lower = 0> sigma_raw;
real<lower = 0> sigma_juv_raw;
real alpha_raw;
real beta_raw;
mu_raw[NS,MY,MJ];
mu_juv_raw[NS,MY,MJ];
}

transformed parameters {
real sigma;
real sigma_juv;
real alpha;
real beta;
real mu[NS,MY,MJ];
real mu_juv[NS,MY,MJ];

sigma = sigma_raw * 10;
sigma_juv = sigma_juv_raw * 10;
alpha = alpha_raw * 10;
beta = beta_raw * 10;

for (s in 1:NS)
{
  for (i in 1:NY[s])
  {
    for (j in 1:NJ[s,i])
    {
      // non-centered prior for mu_juv
      mu_juv[s,i,j] = mu_juv_raw[s,i,j] * 10;
      
      // process model
      mu2[s,i,j] = alpha + beta * mu_juv[s,i,j];
      
      // implies mu[s,i,j] ~ normal(mu2[s,i,j], sigma)
      mu[s,i,j] = mu_raw[s,i,j] * sigma + mu2[s,i,j]; 
    }
  }
}

}

model {
sigma_raw ~ std_normal();
sigma_juv_raw ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();


for (s in 1:NS)
{
  for (i in 1:NY[s])
  {
    for (j in 1:NJ[s,i])
    {
      // non centered for mu_juv
      mu_juv_raw[s,i,j] ~ std_normal();
      
      for (k in 1:NK[s,i,j])
      {
        // mean juv date for each species/cell/year
        juv[s,i,j,k] ~ normal(mu_juv[s,i,j], sigma_juv);
      }
      
      // non-centered for mu
      mu_raw[s,i,j] ~ std_normal();
      
      // observation model for arrival
      y[s,i,j] ~ normal(mu[s,i,j], sd_y[s,i,j]);
    }
  }
}
}

// generated quantities {
// real y_rep;

// y_rep = normal_rng(mu, sigma);
// }
"

# stanmodel3 <- "
# data {
# int<lower=0> N;                     // number of obs
# int<lower=0> Nsp;                   // number of species
# real<lower=0> y[N];                 // response
# int<lower=1, upper=Nsp> sp[N];      // species ids
# vector<lower=0>[N] year;
# vector<lower=0>[N] lat;
# vector<lower=0>[N] elev;
# }
# 
# parameters {
# real mu_alpha_raw;
# real mu_beta_raw;
# real mu_gamma_raw;
# real mu_theta_raw;
# vector<lower = 0>[Nsp] sigma_raw;
# vector<lower = 0>[2] sigma_sp_raw;
# cholesky_factor_corr[2] L_Rho;             // correlation matrix
# matrix[2, Nsp] z;                          // z-scores
# }
# 
# transformed parameters {
# vector[N] mu;
# matrix[Nsp, 2] abgt;                              // matrix for alpha, beta, gamma, and theta
# matrix[2, 2] Rho;                                 // covariance matrix
# vector[Nsp] alpha;
# vector[Nsp] beta;
# vector[Nsp] alpha_c;
# vector[Nsp] beta_c;
# 
# vector<lower = 0>[Nsp] sigma;
# vector<lower = 0>[2] sigma_sp;
# real mu_alpha;
# real mu_beta;
# 
# sigma = sigma_raw * 5;
# mu_alpha = mu_alpha_raw * 10 + 70;
# mu_beta = mu_beta_raw * 0.1;
# sigma_sp[1] = sigma_sp_raw[1] * 20;             // variance alpha
# sigma_sp[2] = sigma_sp_raw[2] * 0.1;              // variance beta
# 
# // cholesky factor of covariance matrix multiplied by z score
# ab = (diag_pre_multiply(sigma_sp, L_Rho) * z)';
# alpha = ab[,1];
# beta = ab[,2];
# Rho = L_Rho * L_Rho';
# 
# for (i in 1:N)
# {
#   mu[i] = (mu_alpha + ab[sp[i], 1]) + 
#   (mu_beta + ab[sp[i], 2]) * year[i];
# }
# 
# alpha_c = mu_alpha + alpha;
# beta_c = mu_beta + beta;
# }
# 
# model {
# sigma_raw ~ std_normal();
# mu_alpha_raw ~ std_normal();
# mu_beta_raw ~ std_normal();
# sigma_sp_raw ~ std_normal();
# 
# to_vector(z) ~ std_normal();         // faster than normal(0, 1);
# L_Rho ~ lkj_corr_cholesky(2);
# 
# y ~ normal(mu, sigma[sp]);
# }
# 
# generated quantities {
# real y_rep[N];
# 
# y_rep = normal_rng(mu, sigma[sp]);
# }
# "

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
                            'sigma_juv',
                            'sigma',
                            'mu_juv'), 
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
 




