######################
# 7 - arrival ~ juv hitting nets
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
juv_date <- '2019-08-21'

#IAR data
arr_master_dir <- 'arrival_master_2019-05-26'



# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/arr_br_', juv_date))

#read in 
juvs_master <- readRDS(paste0('juv-output-', juv_date, '.rds'))

#master arrival data (from IAR output)
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_master_dir))

arr_master <- readRDS(paste0(arr_master_dir, '.rds'))

#merge data sets
mrg <- dplyr::inner_join(juvs_master, arr_master, by = c('species', 'cell', 'year'))

#only species/cells/years with data for juvs and input data for IAR
mrg_f <- dplyr::filter(mrg, !is.na(juv_mean), !is.na(mean_pre_IAR))


# #only species that have more than 3 data points
# cnt_arr <- plyr::count(mrg_f, 'species')
# sp_f <- filter(cnt_arr, freq > 3)$species

sp_idx <- as.numeric(factor(mrg_f$species))


# Stan model --------------------------------------------------------------

DATA <- list(y = mrg_f$mean_post_IAR,
             sd_y = mrg_f$sd_post_IAR,
             juv = mrg_f$juv_mean,
             sd_juv = mrg_f$juv_sd,
             N = NROW(mrg_f),
             sp = sp_idx,
             Nsp = length(unique(sp_idx)))

stanmodel1 <- "
data {
int<lower=0> N;                      // number of data points
vector<lower=0>[N] y;                  // response
vector<lower=0>[N] sd_y;               // uncertainty in response
vector<lower=0>[N] juv;
vector<lower=0>[N] sd_juv;
int<lower=0> sp[N];              
int<lower=0> Nsp;
}

parameters {
real<lower = 0> sigma_raw;
vector[N] mu_y_raw;
vector[N] mu_juv_raw;
real mu_alpha_raw;
real mu_beta_raw;
vector<lower = 0>[2] sigma_sp_raw;
cholesky_factor_corr[2] L_Rho;             // correlation matrix
matrix[2, Nsp] z;                          // z-scores
}

transformed parameters {
real<lower = 0> sigma;
vector[N] mu_y;
vector[N] mu_juv;
vector[N] mu;
matrix[Nsp, 2] ab;                              // matrix for alpha, beta, gamma, and theta
matrix[2, 2] Rho;                                 // covariance matrix
vector[Nsp] alpha;
vector[Nsp] beta;
vector[Nsp] alpha_c;
vector[Nsp] beta_c;
vector<lower = 0>[2] sigma_sp;
real mu_alpha;
real mu_beta;

sigma = sigma_raw * 10;
mu_juv = mu_juv_raw * 100 + 200;
mu_alpha = mu_alpha_raw * 200;
mu_beta = mu_beta_raw * 5;
sigma_sp[1] = sigma_sp_raw[1] * 200;             // variance alpha
sigma_sp[2] = sigma_sp_raw[2] * 5;              // variance beta

// cholesky factor of covariance matrix multiplied by z score
ab = (diag_pre_multiply(sigma_sp, L_Rho) * z)';
alpha = ab[,1];
beta = ab[,2];
Rho = L_Rho * L_Rho';

for (i in 1:N)
{
  mu[i] = (mu_alpha + ab[sp[i], 1]) +
  (mu_beta + ab[sp[i], 2]) * mu_juv[i];
}

alpha_c = mu_alpha + alpha;
beta_c = mu_beta + beta;

// implies mu_y[i] ~ normal(mu[i], sigma)
mu_y = mu_y_raw * sigma + mu; 

}

model {
sigma_raw ~ std_normal();
sigma_sp_raw ~ std_normal();
mu_alpha_raw ~ std_normal();
mu_beta_raw ~ std_normal();
mu_y_raw ~ std_normal();
mu_juv_raw ~ std_normal();

to_vector(z) ~ std_normal();
L_Rho ~ lkj_corr_cholesky(2);

// observation model for juvs
juv ~ normal(mu_juv, sd_juv);

// observation model for arrival
y ~ normal(mu_y, sd_y);

}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu_y, sd_y);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.90
TREE_DEPTH <- 15
STEP_SIZE <- 0.001
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha_c',
                            'beta_c',
                            'alpha',
                            'beta',
                            'mu_alpha',
                            'mu_beta',
                            'sigma_sp',
                            'Rho',
                            'sigma',
                            'mu_juv',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/arr_br_', juv_date))
saveRDS(fit, file = paste0('arr-br-stan-output', juv_date, '.rds'))


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
model_summary <- MCMCvis::MCMCsummary(fit, Rhat = TRUE, n.eff = TRUE, round = 2, excl = 'y_rep')

#extract Rhat and neff values
rhat_output <- as.vector(model_summary[, grep('Rhat', colnames(model_summary))])
neff_output <- as.vector(model_summary[, grep('n.eff', colnames(model_summary))])

#y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')

# bayesplot::ppc_stat(DATA$y_obs, y_rep, stat = 'mean')
# bayesplot::ppc_dens_overlay(DATA$y_obs, y_rep[1:500,])



# write model results to file ---------------------------------------------

options(max.print = 50000)
sink(paste0('arr-br-stan-results-', DATE, '.txt'))
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
 




# new plot ----------------------------------------------------------------



data_vis_fun <- function(SPECIES = 'all')
{
  #SPECIES <- 'Vireo_olivaceus'
  
  #extract posterior estimates for true states for y and x
  y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                   func = mean)[[1]]
  y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.025)))[[1]]
  y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  x_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_juv', type = 'summary', 
                                   func = mean)[[1]]
  x_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_juv', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.025)))[[1]]
  x_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_juv', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  #need true latent states
  DATA_PLOT <- data.frame(mean_y = y_true_mean,
                          mean_y_l = y_true_LCI,
                          mean_y_u = y_true_UCI,
                          mean_x = x_true_mean, 
                          mean_x_l = x_true_LCI,
                          mean_x_u = x_true_UCI,
                          sp_id = factor(mrg_f$species))
  
  if (SPECIES == 'all')
  {
    #model fit for mu_beta and mu_alpha
    alpha_ch <- MCMCchains(fit, params = 'mu_alpha')[,1]
    beta_ch <- MCMCchains(fit, params = 'mu_beta')[,1]
    
    DATA_PLOT2 <- DATA_PLOT
  } else {
    idx <- which(unique(mrg_f$species) == SPECIES)
    if (length(idx) > 0)
    {
      alpha_ch <- MCMCchains(fit, params = paste0('alpha\\[', idx, '\\]'), ISB = FALSE)[,1]
      beta_ch <- MCMCchains(fit, params = paste0('beta\\[', idx, '\\]'), ISB = FALSE)[,1]
      
      DATA_PLOT2 <- dplyr::filter(DATA_PLOT, sp_id == SPECIES)
    } else {
      stop(paste0('Species: ', SPECIES, ' not found!'))
    }
  }
  
  sim_x <- seq(min(DATA_PLOT2$mean_x_l) - 1, max(DATA_PLOT2$mean_x_u) + 1, length = 100)
  
  mf <- matrix(nrow = length(beta_ch), ncol = 100)
  for (i in 1:length(sim_x))
  {
    mf[,i] <- alpha_ch + beta_ch * sim_x[i]
  }
  
  med_mf <- apply(mf, 2, median)
  LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
  UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
  
  FIT_PLOT <- data.frame(MN = med_mf,
                         MN_X = sim_x,
                         LCI = LCI_mf,
                         UCI = UCI_mf)
  
  p <- ggplot(data = DATA_PLOT2, aes(mean_x, mean_y)) +
    geom_ribbon(data = FIT_PLOT,
                aes(x = MN_X, ymin = LCI, ymax = UCI),
                fill = 'grey', alpha = 0.7,
                inherit.aes = FALSE) +
    geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
              alpha = 0.9,
              inherit.aes = FALSE,
              size = 1.4) +
    geom_errorbar(data = DATA_PLOT2, 
                  aes(ymin = mean_y_l, ymax = mean_y_u), width = 0.3,
                  color = 'black', alpha = 0.2) +
    geom_errorbarh(data = DATA_PLOT2, 
                   aes(xmin = mean_x_l, xmax = mean_x_u), height = 0.005,
                   color = 'black', alpha = 0.2) +
    # geom_point(data = DATA_PLOT2, aes(mean_x, mean_y), color = 'black',
    #            inherit.aes = FALSE, size = 1, alpha = 0.3) +
    geom_point(data = DATA_PLOT2, aes(mean_x, mean_y, color = sp_id),
               inherit.aes = FALSE, size = 1, alpha = 0.3) +
    theme_bw() +
    #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
    ylab('True ARR halfmax') +
    xlab('True BR halfmax') +
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
data_vis_fun(SPECIES = 'all')

#each species individually
sps <- unique(mrg_f$species)
for (i in 1:length(sps))
{
  #i <- 3
  data_vis_fun(SPECIES = sps[i])
}




################################
#vvvvvvvvv DEPRECATED? vvvvvvvvv
################################


# plot --------------------------------------------------------------------

#extract posterior estimates for true states for y and x
y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                 func = mean)[[1]]
y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                func = function(x) quantile(x, probs = c(0.025)))[[1]]
y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_y', type = 'summary', 
                                func = function(x) quantile(x, probs = c(0.975)))[[1]]

x_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_juv', type = 'summary', func = mean)[[1]]
x_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_juv', type = 'summary', func = function(x) quantile(x, probs = c(0.025)))[[1]]
x_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_juv', type = 'summary', func = function(x) quantile(x, probs = c(0.975)))[[1]]

#need true latent states
DATA_PLOT <- data.frame(mean_y = y_true_mean,
                        mean_y_l = y_true_LCI,
                        mean_y_u = y_true_UCI,
                        mean_x = x_true_mean, 
                        mean_x_l = x_true_LCI,
                        mean_x_u = x_true_UCI,
                        sp_id = factor(arr_data$species))
  
#model fit for mu_beta and mu_alpha
alpha_ch <- MCMCchains(fit, params = 'mu_alpha')[,1]
beta_ch <- MCMCchains(fit, params = 'mu_beta')[,1]

sim_x <- seq(min(DATA_PLOT$mean_x_l) - 1, max(DATA_PLOT$mean_x_u) + 1, length = 100)
  
mf <- matrix(nrow = length(beta_ch), ncol = 100)
for (i in 1:length(sim_x))
{
  mf[,i] <- alpha_ch + beta_ch * sim_x[i]
}

med_mf <- apply(mf, 2, median)
LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))

FIT_PLOT <- data.frame(MN = med_mf,
                       MN_X = sim_x,
                       LCI = LCI_mf,
                       UCI = UCI_mf)

p <- ggplot(data = DATA_PLOT, aes(mean_x, mean_y)) +
  # geom_ribbon(data = FIT_PLOT,
  #             aes(x = MN_X, ymin = LCI, ymax = UCI),
  #             fill = 'grey', alpha = 0.7,
  #             inherit.aes = FALSE) +
  # geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
  #           alpha = 0.9,
  #           inherit.aes = FALSE,
  #           size = 1.4) +
  geom_errorbar(data = DATA_PLOT, 
                aes(ymin = mean_y_l, ymax = mean_y_u), width = 0.3,
                color = 'black', alpha = 0.2) +
  geom_errorbarh(data = DATA_PLOT, 
                 aes(xmin = mean_x_l, xmax = mean_x_u), height = 0.005,
                 color = 'black', alpha = 0.2) +
  geom_point(data = DATA_PLOT, aes(mean_x, mean_y, color = sp_id),
             inherit.aes = FALSE, size = 3, alpha = 0.3) +
  theme_bw() +
  #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
  xlab('BR') +
  ylab('ARR') +
  #ggtitle(paste0('Species: ', SPECIES)) +
  theme(
    plot.title = element_text(size = 22),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick

print(p)



arr <- dplyr::filter(arr_data2, species == 'Geothlypis_trichas')

means <- apply(juv_array, 1, function(x) mean(x[which(x > 0)]))
da <- apply(juv_array, 1, function(x) sum(x > 0))
GT <- which(arr_data2$species == 'Geothlypis_trichas')

#number of birds for each cell/year
da[GT]

plot(means[GT], tt$mean_post_IAR)
