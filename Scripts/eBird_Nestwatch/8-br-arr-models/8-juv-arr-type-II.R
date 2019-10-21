######################
# 8 - juv hitting nets ~ arrival
#
# y_{obs_{i}} \sim N(y_{true_{i}}, \sigma_{y_{i}})
# 
# x_{obs_{i}} \sim N(x_{true_{i}}, \sigma_{x_{i]})
# 
# \begin{bmatrix} y_{true_{i}} \\ x_{true_{i}} \end{bmatrix} \sim MVN \left(\begin{bmatrix} 0 \\ 0 \end{bmatrix}, \Sigma \right)
# 
# \beta = \frac{\sigma_{x} - \sigma_{y} + \sqrt{(\sigma_{x} - \sigma_{y})^{2} + 4 \sigma_{xy}}}{2 \sqrt{\sigma_{xy}}}
# 
# \alpha = \overline{y_{obs}} - \beta * \overline{x_{obs}}
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
juv_master_dir <- 'juv_master_2019-10-15'

#IAR data
arr_master_dir <- 'arrival_master_2019-05-26'



# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS

#setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', juv_master_dir))
setwd(paste0('~/Google_Drive/', juv_master_dir))

#read in 
juvs_master <- readRDS(paste0(juv_master_dir, '.rds'))

#master arrival data (from IAR output)
#setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', arr_master_dir))
setwd(paste0('~/Google_Drive/', arr_master_dir))

arr_master <- readRDS(paste0(arr_master_dir, '.rds'))

#merge data sets
mrg <- dplyr::inner_join(juvs_master, 
                         arr_master, by = c('species', 'cell', 'year'))

#only species/cells/years with data for juvs
mrg_f <- dplyr::filter(mrg, !is.na(juv_mean), !is.na(mean_pre_IAR))

#only species that have at least 5 data points
cnt_arr <- plyr::count(mrg_f, 'species')
sp_f <- filter(cnt_arr, freq >= 5)$species

mrg_f2 <- dplyr::filter(mrg_f, species %in% sp_f)
sp_idx <- as.numeric(factor(mrg_f2$species))


# Stan model --------------------------------------------------------------

DATA <- list(y = mrg_f2$juv_mean,
             sd_y = mrg_f2$juv_sd,
             arr = mrg_f2$mean_post_IAR,
             sd_arr = mrg_f2$sd_post_IAR,
             N = NROW(mrg_f2),
             sp = sp_idx,
             Nsp = length(unique(sp_idx)))


# Latex for orthogonal slope (specific case of Deming regression where errors are known, I believe):
# \beta = \frac{\sigma_{x} - \sigma_{y} + \sqrt{(\sigma_{x} - \sigma_{y})^{2} + 4 \sigma_{xy}}}{2 \sqrt{\sigma_{xy}}}
# helpful: 
# https://www.wikiwand.com/en/Deming_regression
# https://davegiles.blogspot.com/2014/11/orthogonal-regression-first-steps.html
# https://bayes.wustl.edu/etj/articles/leapz.pdf

# less helpful:
# https://stats.stackexchange.com/questions/6163/what-is-the-prediction-error-while-using-deming-regression-weighted-total-least
# https://stats.stackexchange.com/questions/143378/bayesian-estimates-for-deming-regression-coinciding-with-least-squares-estimates

# hierarchical prior on covariance matrix: 
# https://discourse.mc-stan.org/t/covariance-matrix-with-more-hierarchy-levels/696
# and
# http://modernstatisticalworkflow.blogspot.com/2017/04/hierarchical-vector-autoregression.html


stanmodel1 <- "
data {
int<lower=0> N;                        // number of data points
vector<lower=0>[N] y;                  // response
vector<lower=0>[N] sd_y;               // uncertainty in response
vector<lower=0>[N] arr;
vector<lower=0>[N] sd_arr;
// int<lower=0> sp[N];              
// int<lower=0> Nsp;
}

parameters {
vector[N] mu_y;
vector[N] mu_arr;
vector<lower=0>[2] sigma_aj_raw;
corr_matrix[2] Omega;                       // corr matrix
vector[2] gamma_raw;
}

transformed parameters {
vector[2] aj[N];                                  // matrix for arr and juv
cov_matrix[2] Sigma;                              // covariance matrix
real beta;                                        // orthogonal slope
vector<lower=0>[2] sigma_aj;                      // sd for arr and juv
vector[2] gamma;                                  // means for MVN

// sd of arrival and juv respectively
sigma_aj[1] = sigma_aj_raw[1] * 20;
sigma_aj[2] = sigma_aj_raw[2] * 20;

gamma[1] = gamma_raw[1] * 30 + 200;      // prior mean juv
gamma[2] = gamma_raw[2] * 30 + 125;      // prior mean arr

// fill aj matrix with arr and juv values - 2nd dim is 'vector' in mixed object
for (i in 1:N)
{
  aj[i, 1] = mu_y[i];
  aj[i, 2] = mu_arr[i];
}

// derived covariance matrix - diag(sigma_aj) * Omega * diag(sigma_aj)
Sigma = quad_form_diag(Omega, sigma_aj);


// Orthogonal regression slope from Dave Miller in Slack General Channel ~ Sep 3, 2019
// see links about for derivation of intercept parameter
// Sigma[2,1] is cov_xy
beta = (sigma_aj[1] - sigma_aj[2] + sqrt((sigma_aj[1] - sigma_aj[2])^2 + 4 * Sigma[2,1])) / 2 * sqrt(Sigma[2,1]);
}

model {
sigma_aj_raw ~ std_normal();
Omega ~ lkj_corr(2);
gamma_raw  ~ std_normal();


// observation model for arr
arr ~ normal(mu_arr, sd_arr);

// observation model for juveniles
y ~ normal(mu_y, sd_y);

// estimate covariance matrix
aj ~ multi_normal(gamma, Sigma);

}

generated quantities {
real y_rep[N];
real y_bar;
real x_bar;
real alpha;

y_rep = normal_rng(mu_y, sd_y);
y_bar = mean(y);
x_bar = mean(arr);
alpha = y_bar - beta * x_bar;
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 13
STEP_SIZE <- 0.001
CHAINS <- 3
ITER <- 1000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha',
                            'beta',
                            'gamma',
                            'aj',
                            'sigma_aj',
                            'Omega',
                            'Sigma',
                            'mu_arr',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60



MCMCsummary(fit, params = c('alpha', 'beta', 'gamma'))

aj <- MCMCvis::MCMCpstr(fit, params = 'aj')[[1]]

mu_arr <- MCMCvis::MCMCpstr(fit, params = 'mu_arr')[[1]]
mu_y <- MCMCvis::MCMCpstr(fit, params = 'mu_y')[[1]]



alpha <- MCMCvis::MCMCpstr(fit, params = 'alpha')[[1]]
beta <- MCMCvis::MCMCpstr(fit, params = 'beta')[[1]]
tt <- cbind(DATA$y,DATA$arr)

tt - aj


plot(0,0, xlim = c(-1000, 1000), ylim = c(-1000, 1000))
abline(a = alpha, b = beta)
points(tt)




setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_', juv_date))
saveRDS(fit, file = paste0('juv-arr-stan-output-', juv_date, '.rds'))
#saveRDS(mrg_f2, file = paste0('mrg-data-juv-', juv_date, '.rds'))
#fit <- readRDS(paste0('juv-arr-stan-output-', juv_date, '.rds'))

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



# write model results to file ---------------------------------------------

options(max.print = 50000)
sink(paste0('juv-arr-stan-results-', juv_date, '.txt'))
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
# 
# mu_alpha = mu_alpha_raw * 200;
# mu_beta = mu_beta_raw * 2;
# sigma = sigma_raw * 10;
# mu_arr = mu_arr_raw * 40 + 180;
# sigma_sp[1] = sigma_sp_raw[1] * 20;
# sigma_sp[2] = sigma_sp_raw[2] * 1;
# 
# 
# 
#mu_alpha ~ N(0, 200)
PR <- rnorm(10000, 0, 200)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_alpha',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-arr-', juv_date, '-trace_mu_alpha.pdf'))

#mu_beta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_beta',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-arr-', juv_date, '-trace_mu_beta.pdf'))

#sigma_sp[1] ~ HN(0, 40)
PR_p <- rnorm(10000, 0, 40)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sp\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-arr-', juv_date, '-trace_sigma_sp[1].pdf'))

#sigma_sp[2] ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sp\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-arr-', juv_date, '-trace_sigma_sp[2].pdf'))


#sigma ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-arr-', juv_date, '-trace_sigma.pdf'))

#mu_arr ~ N(180, 40)
PR <- rnorm(10000, 180, 40)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_arr',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-arr-', juv_date, '-trace_mu_juv.pdf'))



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
  
  x_true_mean <- MCMCvis::MCMCpstr(fit, params = 'mu_arr', type = 'summary', 
                                   func = mean)[[1]]
  x_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_arr', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.025)))[[1]]
  x_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_arr', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  DATA_PLOT <- data.frame(mean_y = y_true_mean,
                          mean_y_l = y_true_LCI,
                          mean_y_u = y_true_UCI,
                          mean_x = x_true_mean, 
                          mean_x_l = x_true_LCI,
                          mean_x_u = x_true_UCI,
                          y_obs = DATA$y,
                          y_obs_l = DATA$y - (1.96 * DATA$sd_y),
                          y_obs_u = DATA$y + (1.96 * DATA$sd_y),
                          x_obs = DATA$arr,
                          x_obs_l = DATA$arr - (1.96 * DATA$sd_arr),
                          x_obs_u = DATA$arr + (1.96 * DATA$sd_arr),
                          sp_id = factor(mrg_f2$species))
  
  
  if (SPECIES == 'all')
  {
    #model fit for mu_beta and mu_alpha
    a_ch <- MCMCchains(fit, params = 'mu_alpha')[,1]
    b_ch <- MCMCchains(fit, params = 'mu_beta')[,1]
    
    DATA_PLOT2 <- DATA_PLOT
  } else {
    
    idx <- which(unique(mrg_f2$species) == SPECIES)
    if (length(idx) > 0)
    {
      a_ch <- MCMCvis::MCMCchains(fit, ISB = FALSE, params = paste0('alpha_c\\[', idx, '\\]'))[,1]
      b_ch <- MCMCvis::MCMCchains(fit, ISB = FALSE, params = paste0('beta_c\\[', idx, '\\]'))[,1]
      
      DATA_PLOT2 <- dplyr::filter(DATA_PLOT, sp_id == SPECIES)
    } else {
      stop(paste0('Species: ', SPECIES, ' not found!'))
    }
  }
  
  sim_x <- seq(min(DATA_PLOT2$mean_x_l) - 1, max(DATA_PLOT2$mean_x_u) + 1, length = 100)
  
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
    geom_errorbarh(data = DATA_PLOT2, 
                   aes(xmin = mean_x_l, xmax = mean_x_u), height = 0.005,
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
    xlab('ARR halfmax') +
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
pdf('juv-arr-ALL.pdf')
data_vis_fun(SPECIES = 'all')
dev.off()

#each species individually
sps <- unique(mrg_f2$species)
for (i in 1:length(sps))
{
  #i <- 3
  pdf(paste0('juv-arr-', sps[i], '.pdf'))
  data_vis_fun(SPECIES = sps[i])
  dev.off()
}

