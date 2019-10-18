######################
#NO MODEL ALPHA
#
# y_{i} ~ N(mu_{i}, sigma)
# mu_{i} = alpha_{j} + beta_{j} * year_{i}
# alpha_{j} ~ N(mu_alpha_{k}, sigma_alpha)
# beta_{j} ~ N(mu_beta_{j}, sigma_beta)
# mu_alpha_{k} ~ N(mm, sm)
# mu_beta_{j} = pi_{k} + nu_{k} * lat_{j}
# [pi_k, nu_k] ~ MVN(mu_pi, mu_nu, Sigma_pn)
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
juv_date <- '2019-10-15'
run_date <- '2019-10-17'


# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/juv_master_', juv_date))

#read in 
juvs_master <- readRDS(paste0('juv-output-', juv_date, '.rds'))

#only species/cells/years with data for juvs
j1 <- dplyr::filter(juvs_master, !is.na(juv_mean))

#species/cells with at least 3 years
cnt_arr <- plyr::count(j1, c('species', 'cell'))
ff <- filter(cnt_arr, freq >= 3)
fcr <- data.frame(species = ff[,'species'], cell = ff[,'cell'])
j2 <- dplyr::inner_join(fcr, j1)

#add cell lat to df
hexgrid6 <- dggridR::dgconstruct(res = 6)
j2$cell_lat <- dggridR::dgSEQNUM_to_GEO(hexgrid6, 
                                        in_seqnum = j2$cell)$lat_deg
j2$cell_lng <- dggridR::dgSEQNUM_to_GEO(hexgrid6, 
                                        in_seqnum = j2$cell)$lon_deg

#more than 5 degrees lat difference between cells with data
cnt_j2 <- aggregate(cell_lat ~ species, data = j2, FUN = function(x) max(x) - min(x))
g_idx <- cnt_j2[,2] > 5
sp_k <- cnt_j2[g_idx,1]

j3 <- dplyr::filter(j2, species %in% sp_k)
j3$sp_idx <- as.numeric(factor(j3$species))
usp <- unique(j3$sp_idx)


#add cell id
j3$cn_id <- NA
cn_id <- c()
counter <- 0
for (i in 1:length(usp))
{
  #i <- 1
  t_idx <- which(j3$sp_id == usp[i])
  temp <- j3[t_idx,]
  tj <- as.numeric(factor(temp$cell))
  j3$cn_id[t_idx] <- counter + tj
  
  counter <- counter + max(tj)
}

#species index and cell_lat
u_cn_id <- sort(unique(j3$cn_id))
sp_id <- c()
cell_lat <- c()
for (i in 1:length(u_cn_id))
{
  #i <- 1
  s_idx <- which(j3$cn_id == u_cn_id[i])
  sn <- j3[s_idx,'sp_idx'][1]
  sp_id <- c(sp_id, sn)
  
  cl <- j3[s_idx,'cell_lat'][1]
  cell_lat <- c(cell_lat, cl)
}



# quick freq check --------------------------------------------------------

# plyr::count(j2, 'species')
# 
# tt <- dplyr::filter(j2, species == 'Geothlypis_trichas')
# tt_h <- dplyr::filter(tt,  cell_lat > 45)
# tt_l <- dplyr::filter(tt,  cell_lat < 40)
# 
# summary(lm(juv_mean ~ year, data = tt_h))
# summary(lm(juv_mean ~ year, data = tt_l))



# setwd('~/Desktop')
# sink('output.txt')
# for (i in 1:length(sp_f))
# {
#   #i <- 3
#   temp <- dplyr::filter(j2, species == sp_f[i])
# 
#   t_rng <- quantile(temp$cell_lat, probs = c(0.1, 0.90))
#   tt_h <- dplyr::filter(temp,  cell_lat >= t_rng[2])
#   tt_l <- dplyr::filter(temp,  cell_lat <= t_rng[1])
#   if (NROW(tt_h) > 3 & NROW(tt_l) > 3)
#   {
#     f1 <- summary(lm(juv_mean ~ year, data = tt_h))
#     f2 <- summary(lm(juv_mean ~ year, data = tt_l))
#     plot(tt_h$year, tt_h$juv_mean)
#     plot(tt_l$year, tt_l$juv_mean)
#     print(sp_f[i])
#     print(paste0('slope high lat: ', round(f1$coefficients[2,1], 2)))
#     print(paste0('pval high lat: ', round(f1$coefficients[2,4], 2)))
#     print(paste0('slope low lat: ', round(f2$coefficients[2,1], 2)))
#     print(paste0('pval low lat: ', round(f2$coefficients[2,4], 2)))
#     print('')
#   }
# }
# sink()


# Stan model --------------------------------------------------------------

DATA <- list(N = NROW(j3),
             NC = length(u_cn_id),
             Nsp = length(unique(sp_id)),
             y = j3$juv_mean,
             sd_y = j3$juv_sd,
             sp_id = sp_id,
             cn_id = j3$cn_id,
             year = as.numeric(factor(j3$year)),
             cell_lat = cell_lat,
             P = 2)


stanmodel1 <- "
data {
int<lower=0> N;                      // number of data points
int<lower=0> NC;                      // number of species/cells
int<lower=0> Nsp;                     // number of species
vector<lower=0>[N] y;                  // response
vector<lower=0>[N] sd_y;               // uncertainty in response
int<lower=0> sp_id[NC];                 // species ids
int<lower=0> cn_id[N];                // species/cell ids
vector<lower=0>[N] year;
vector<lower=0>[NC] cell_lat;
int<lower=0> P;
}

parameters {
real<lower = 0> sigma_raw;
vector[N] mu_y_raw;
// real mu_gamma_raw;
// real mu_theta_raw;
real mu_pi_raw;
real mu_nu_raw;
// vector<lower = 0>[P] sigma_ab_raw;
// vector<lower = 0>[P] sigma_gt_raw;
vector<lower = 0>[P] sigma_pn_raw;
// cholesky_factor_corr[P] L_Rho_ab;                 // cholesky factor of corr matrix
// cholesky_factor_corr[P] L_Rho_gt;
cholesky_factor_corr[P] L_Rho_pn;
// matrix[P, NC] z_ab;
// matrix[P, Nsp] z_gt;
matrix[P, Nsp] z_pn;
real mu_mu_alpha_raw;
real<lower = 0> sigma_mu_alpha_raw;
vector[Nsp] mu_alpha_raw;
vector[NC] alpha_raw;
vector[NC] beta_raw;
real <lower = 0> sigma_alpha_raw;
real <lower = 0> sigma_beta_raw;
}

transformed parameters {
real<lower = 0> sigma;
vector[N] mu_y;
vector[N] mu;
// matrix[NC, P] ab;
// matrix[Nsp, P] gt;
matrix[Nsp, P] pn;
// matrix[P, P] Rho_ab;    // correlation matrix
// matrix[P, P] Rho_gt;
matrix[P, P] Rho_pn;
vector[Nsp] mu_alpha;
vector[NC] mu_beta;
// real mu_gamma;
// real mu_theta;
real mu_pi;
real mu_nu;
vector[NC] alpha;
vector[NC] beta;
// vector[Nsp] gamma;
// vector[Nsp] theta;
vector[Nsp] pi;
vector[Nsp] nu;
// vector<lower = 0>[P] sigma_ab;
// vector<lower = 0>[P] sigma_gt;
vector<lower = 0>[P] sigma_pn;
real mu_mu_alpha;
real <lower = 0> sigma_mu_alpha;
real <lower = 0> sigma_alpha;
real <lower = 0> sigma_beta;

// mu_gamma = mu_gamma_raw * 100;
// mu_theta = mu_theta_raw * 2;
mu_pi = mu_pi_raw * 2;
mu_nu = mu_nu_raw * 2;

sigma = sigma_raw * 10;
// sigma_ab[1] = sigma_ab_raw[1] * 10;
// sigma_ab[2] = sigma_ab_raw[2] * 10;
// sigma_gt[1] = sigma_gt_raw[1] * 40;
// sigma_gt[2] = sigma_gt_raw[2] * 1;
sigma_pn[1] = sigma_pn_raw[1] * 1;
sigma_pn[2] = sigma_pn_raw[2] * 1;
mu_mu_alpha = mu_mu_alpha_raw * 40 + 200;
sigma_mu_alpha = sigma_mu_alpha_raw * 20;
sigma_alpha = sigma_alpha_raw * 30;
sigma_beta = sigma_beta_raw * 30;

// cholesky factor of covariance matrix (i.e., diagonal matrix of scale times cholesky factor of correlation matrix) multiplied by z score
// cholesky factor transforms uncorrelated variables (z scores) into variables whose variances and covariances are given by Sigma (i.e., sigma[diag of scale] * Rho[corr matrix] * sigma) and are centered on 0 (since corr_xy = cov_xy / (sigma_x * simga_y)) - remember to use matrix multiplication
// implies gt ~ MVN(0, Sigma)

// ab = (diag_pre_multiply(sigma_ab, L_Rho_ab) * z_ab)';
// Rho_ab = L_Rho_ab * L_Rho_ab';

// gt = (diag_pre_multiply(sigma_gt, L_Rho_gt) * z_gt)';
// Rho_gt = L_Rho_gt * L_Rho_gt';

pn = (diag_pre_multiply(sigma_pn, L_Rho_pn) * z_pn)';
Rho_pn = L_Rho_pn * L_Rho_pn';

// gamma = mu_gamma + gt[,1];
// theta = mu_theta + gt[,2];
pi = mu_pi + pn[,1];
nu = mu_nu + pn[,2];

mu_alpha = mu_alpha_raw * sigma_mu_alpha + mu_mu_alpha;

for (j in 1:NC)
{
  // mu_alpha[j] = gamma[sp_id[j]] + theta[sp_id[j]] * cell_lat[j];
  mu_beta[j] = pi[sp_id[j]] + nu[sp_id[j]] * cell_lat[j];
}

// alpha = mu_alpha + ab[,1];
// beta = mu_beta + ab[,2];

alpha = alpha_raw * sigma_alpha + mu_alpha[sp_id];
beta = beta_raw * sigma_beta + mu_beta;       // beta_j ~ N(mu_beta_j, sigma_beta)

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
}

// implies mu_y[i] ~ normal(mu[i], sigma)
mu_y = mu_y_raw * sigma + mu; 

}

model {
// sigma_ab_raw ~ std_normal();
// sigma_gt_raw ~ std_normal();
sigma_pn_raw ~ std_normal();
// mu_gamma_raw ~ std_normal();
// mu_theta_raw ~ std_normal();
mu_mu_alpha_raw ~ std_normal();
sigma_mu_alpha_raw ~ std_normal();
mu_alpha_raw ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();
sigma_alpha_raw ~ std_normal();
sigma_beta_raw ~ std_normal();
mu_pi_raw ~ std_normal();
mu_nu_raw ~ std_normal();
sigma_raw ~ std_normal();
mu_y_raw ~ std_normal();

// to_vector(z_ab) ~ std_normal();
// L_Rho_ab ~ lkj_corr_cholesky(1);
// to_vector(z_gt) ~ std_normal();
// L_Rho_gt ~ lkj_corr_cholesky(1);
to_vector(z_pn) ~ std_normal();
L_Rho_pn ~ lkj_corr_cholesky(1);

// observation model for juveniles
y ~ normal(mu_y, sd_y);

}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu_y, sd_y);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.97
TREE_DEPTH <- 15
STEP_SIZE <- 0.003
CHAINS <- 4
ITER <- 5000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'mu_alpha',
                            'mu_beta',
                            'sigma_alpha',
                            'sigma_beta',
                            'mu_mu_alpha',
                            'sigma_mu_alpha',
                            #'gamma',
                            #'theta',
                            'pi',
                            'nu',
                            #'ab',
                            #'gt',
                            'pn',
                            #'mu_gamma',
                            #'mu_theta',
                            'mu_pi',
                            'mu_nu',
                            #'sigma_ab',
                            #'sigma_gt',
                            'sigma_pn',
                            #'Rho_ab',
                            #'Rho_gt',
                            'Rho_pn',
                            'sigma',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd('~/Desktop/')
#setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_', juv_date))
saveRDS(fit, file = paste0('juv-joint-cell-no-alpha-model-stan-output-', run_date, '.rds'))
#fit <- readRDS('juv-time-bad-stan-output-2019-08-31.rds')

MCMCvis::MCMCsummary(fit, 
                     params = c('pi', 'nu'), 
                     round = 3)
MCMCvis::MCMCsummary(fit, 
                     params = c('mu_pi', 'mu_nu'), 
                     round = 4)
MCMCvis::MCMCsummary(fit, 
                     params = c('mu_alpha', 'mu_beta'), 
                     round = 4)
MCMCvis::MCMCsummary(fit, 
                     params = 'sigma',
                     ISB = FALSE,
                     round = 4)
MCMCvis::MCMCplot(fit, params = 'beta', ref_ovl = TRUE)
MCMCvis::MCMCplot(fit, params = 'beta', ref_ovl = TRUE)
MCMCvis::MCMCplot(fit, params = 'mu_beta', rank = TRUE)
MCMCvis::MCMCplot(fit, params = 'pi')
MCMCvis::MCMCplot(fit, params = 'nu')

beta_s <- MCMCvis::MCMCsummary(fit, params = 'beta')
which(beta_s[,'97.5%'] < 0)
num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))



# Calc diagnostics ---------------------------------------------------

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
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, 
                                      n.eff = TRUE, 
                                      round = 2, 
                                      excl = 'mu_y')

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])



# write model results to file ---------------------------------------------

options(max.print = 50000)
sink(paste0('juv-time-stan-results-', run_date, '.txt'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Iterations: ', ITER, ' \n'))
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



# # PPC ---------------------------------------------------------------------
# 
y_val <- DATA$y
y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
bayesplot::ppc_stat(DATA$y, y_rep, stat = 'mean')
bayesplot::ppc_dens_overlay(DATA$y, y_rep[1:100,])
PPC_fun <- function(FUN, YR = y_rep, D = DATA$y)
{
  out <- sum(apply(YR, 1, FUN) > FUN(D)) / NROW(YR)
  print(out)
}
PPC_fun(mean)
PPC_fun(min)
PPC_fun(max)
# 
# 
# 
# # PPO ---------------------------------------------------------------------

#mu_gamma ~ N(0, 100)
PR <- rnorm(10000, 0, 100)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_gamma',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_mu_gamma.pdf'))

#mu_theta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_theta',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_mu_theta.pdf'))

#mu_pi ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_pi',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_mu_pi.pdf'))

#mu_nu ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_nu',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_mu_nu.pdf'))

#sigma_ab[1] ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_ab\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma_ab[1].pdf'))

#sigma_ab[2] ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_ab\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma_ab[2].pdf'))


#sigma_gt[1] ~ HN(0, 40)
PR_p <- rnorm(10000, 0, 40)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gt\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma_gt[1].pdf'))

#sigma_gt[2] ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_gt\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma_gt[2].pdf'))

#sigma_pn[1] ~ HN(0, 40)
PR_p <- rnorm(10000, 0, 40)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_pn\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma_pn[1].pdf'))

#sigma_pn[2] ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_pn\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma_pn[2].pdf'))

#sigma ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', run_date, '-trace_sigma.pdf'))




# Plot fit ----------------------------------------------------------------


data_vis_fun <- function(SPECIES = 'all')
{
  #SPECIES <- 'Geothlypis_trichas'
  
  #extract posterior estimates for true states for y and x
  y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                   func = mean)[[1]]
  y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.025)))[[1]]
  y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  DATA_PLOT <- data.frame(mean_y = y_true_mean,
                          mean_y_l = y_true_LCI,
                          mean_y_u = y_true_UCI,
                          mean_x = DATA$year, 
                          y_obs = DATA$y,
                          y_obs_l = DATA$y - (1.96 * DATA$sd_y),
                          y_obs_u = DATA$y + (1.96 * DATA$sd_y),
                          x_obs = DATA$year,
                          sp_id = factor(j2$species))
  
  
  if (SPECIES == 'all')
  {
    #model fit for mu_beta and mu_alpha
    a_ch <- MCMCchains(fit, params = 'mu_alpha')[,1]
    b_ch <- MCMCchains(fit, params = 'mu_beta')[,1]
    
    DATA_PLOT2 <- DATA_PLOT
  } else {
    
    idx <- which(unique(j2$species) == SPECIES)
    if (length(idx) > 0)
    {
      a_ch <- MCMCvis::MCMCchains(fit, ISB = FALSE, params = paste0('alpha_c\\[', idx, '\\]'))[,1]
      b_ch <- MCMCvis::MCMCchains(fit, ISB = FALSE, params = paste0('beta_c\\[', idx, '\\]'))[,1]
      
      DATA_PLOT2 <- dplyr::filter(DATA_PLOT, sp_id == SPECIES)
    } else {
      stop(paste0('Species: ', SPECIES, ' not found!'))
    }
  }
  
  sim_x <- seq(min(DATA_PLOT2$mean_x) - 1, max(DATA_PLOT2$mean_x) + 1, length = 100)
  
  mf <- matrix(nrow = length(a_ch), ncol = 100)
  for (i in 1:length(sim_x))
  {
    mf[,i] <- a_ch + b_ch * sim_x[i]
  }
  
  med_mf <- apply(mf, 2, median)
  LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
  UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
  
  FIT_PLOT <- data.frame(MN = med_mf,
                         MN_X = sim_x,
                         LCI = LCI_mf,
                         UCI = UCI_mf)
  
  p <- ggplot(data = DATA_PLOT2, aes(mean_x, mean_y)) +
    #model fit
    geom_ribbon(data = FIT_PLOT,
                aes(x = MN_X, ymin = LCI, ymax = UCI),
                fill = 'grey', alpha = 0.7,
                inherit.aes = FALSE) +
    geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
              alpha = 0.9,
              inherit.aes = FALSE,
              size = 1.4) +
    #latent state
    geom_errorbar(data = DATA_PLOT2, 
                  aes(ymin = mean_y_l, ymax = mean_y_u), width = 0.3,
                  color = 'black', alpha = 0.2) +
    geom_point(data = DATA_PLOT2, aes(mean_x, mean_y), color = 'black',
               inherit.aes = FALSE, size = 3, alpha = 0.3) +
    # geom_point(data = DATA_PLOT2, aes(mean_x, mean_y, color = sp_id),
    #            inherit.aes = FALSE, size = 3, alpha = 0.3) +
    #observed data
    # geom_errorbar(data = DATA_PLOT2,
    #               aes(x = x_obs, ymin = y_obs_l, ymax = y_obs_u), width = 0.3,
    #               color = 'red', alpha = 0.2) +
    # geom_errorbarh(data = DATA_PLOT2,
    #                aes(y = y_obs, xmin = x_obs_l, xmax = x_obs_u), height = 0.005,
    #                color = 'red', alpha = 0.2) +
    # geom_point(data = DATA_PLOT2, aes(x_obs, y_obs), color = 'red',
    #            inherit.aes = FALSE, size = 3, alpha = 0.3) +
  theme_bw() +
    #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
    ylab('BR halfmax') +
    xlab('Year') +
    ggtitle(paste0('Species: ', SPECIES)) +
    theme(
      plot.title = element_text(size = 22),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
      axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
  
  print(p)
}



#all species together
pdf('juv-time-ALL.pdf')
data_vis_fun(SPECIES = 'all')
dev.off()

#each species individually
sps <- unique(j2$species)
for (i in 1:length(sps))
{
  #i <- 3
  pdf(paste0('juv-time-', sps[i], '.pdf'))
  data_vis_fun(SPECIES = sps[i])
  dev.off()
}

