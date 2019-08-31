######################
# ? - juv hitting nets ~ time
#
# juv_{i} \sim N(\mu_{i}, \sigma)
# \mu_{i} = \alpha_{i} + \beta_{i} \times year
# \alpha_{i} \sim N(\mu_{\alpha_{i}}, \sigma_{\alpha})
# \beta_{i} \sim N(\mu_{\beta_{i}}, \sigma_{\beta})
# \mu_{\alpha_{i}} = \gamma_{j} + \theta_{j} \times lat_{i}
# \mu_{\beta_{i}} = \pi_{j} + \nu_{j} \times lat_{i}
# \begin{bmatrix} \gamma_{j} \\ \theta_{j} \end{bmatrix} \sim MVN \left(\begin{bmatrix} \mu_{\gamma} \\ \mu_{\theta} \end{bmatrix}, \Sigma_{\gamma\theta} \right)
# \begin{bmatrix} \pi_{j} \\ \nu_{j} \end{bmatrix} \sim MVN \left(\begin{bmatrix} \mu_{\pi} \\ \mu_{\nu} \end{bmatrix}, \Sigma_{\pi\nu} \right)
#
# see McElreath p. 393
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
juv_date <- '2019-08-26'

#IAR data
arr_master_dir <- 'arrival_master_2019-05-26'



# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_', juv_date))

#read in 
juvs_master <- readRDS(paste0('juv-output-', juv_date, '.rds'))

#only species/cells/years with data for juvs
j1 <- dplyr::filter(juvs_master, !is.na(juv_mean))

#only species that have at least 5 data points
cnt_arr <- plyr::count(j1, 'species')
sp_f <- filter(cnt_arr, freq >= 5)$species

j2 <- dplyr::filter(j1, species %in% sp_f)
sp_idx <- as.numeric(factor(j2$species))


#add cell lat to df
hexgrid6 <- dggridR::dgconstruct(res = 6)
j2$cell_lat <- dggridR::dgSEQNUM_to_GEO(hexgrid6, 
                                        in_seqnum = j2$cell)$lat_deg
j2$cell_lng <- dggridR::dgSEQNUM_to_GEO(hexgrid6, 
                                        in_seqnum = j2$cell)$lon_deg



# Stan model --------------------------------------------------------------

DATA <- list(y = j2$juv_mean,
             sd_y = j2$juv_sd,
             year = j2$year,
             cell_lat = j2$cell_lat,
             N = NROW(j2),
             sp = sp_idx,
             Nsp = length(unique(sp_idx)),
             P = 2)

stanmodel1 <- "
data {
int<lower=0> N;                      // number of data points
vector<lower=0>[N] y;                  // response
vector<lower=0>[N] sd_y;               // uncertainty in response
vector<lower=0>[N] year;
int<lower=0> sp[N];              
int<lower=0> Nsp;
int<lower=0> P;
vector<lower=0>[N] cell_lat;
}

parameters {
real<lower = 0> sigma_raw;
vector[N] mu_y_raw;
real mu_gamma_raw;
real mu_theta_raw;
real mu_pi_raw;
real mu_nu_raw;
vector<lower = 0>[P] sigma_ab_raw;
vector<lower = 0>[P] sigma_gt_raw;
vector<lower = 0>[P] sigma_pn_raw;
cholesky_factor_corr[P] L_Rho_ab;
cholesky_factor_corr[P] L_Rho_gt;
cholesky_factor_corr[P] L_Rho_pn;
matrix[P, N] z_ab;
matrix[P, Nsp] z_gt;
matrix[P, Nsp] z_pn;
}

transformed parameters {
real<lower = 0> sigma;
vector[N] mu_y;
vector[N] mu;
matrix[N, P] gt;
matrix[Nsp, P] gt;
matrix[Nsp, P] pn;
matrix[P, P] Rho_ab;
matrix[P, P] Rho_gt;
matrix[P, P] Rho_pn;
real mu_gamma;
real mu_theta;
real mu_pi;
real mu_nu;
vector[N] alpha;
vector[N] beta;
vector[Nsp] gamma;
vector[Nsp] theta;
vector[Nsp] pi;
vector[Nsp] nu;
vector[Nsp] gamma_c;
vector[Nsp] theta_c;
vector[Nsp] pi_c;
vector[Nsp] nu_c;
vector<lower = 0>[P] sigma_ab;
vector<lower = 0>[P] sigma_gt;
vector<lower = 0>[P] sigma_pn;

mu_gamma = mu_gamma_raw * 200;
mu_theta = mu_theta_raw * 2;
mu_pi = mu_pi_raw * 200;
mu_nu = mu_nu_raw * 2;

sigma = sigma_raw * 10;
sigma_ab[1] = sigma_ab_raw[1] * 40;
sigma_ab[2] = sigma_ab_raw[2] * 10;
sigma_gt[1] = sigma_gt_raw[1] * 40;
sigma_gt[2] = sigma_gt_raw[2] * 1;
sigma_pn[1] = sigma_pn_raw[1] * 40;
sigma_pn[2] = sigma_pn_raw[2] * 1;

// cholesky factor of covariance matrix multiplied by z score
gt = (diag_pre_multiply(sigma_gt, L_Rho_gt) * z_gt)';
gamma = gt[,1];
theta = gt[,2];
Rho_gt = L_Rho_gt * L_Rho_gt';

pn = (diag_pre_multiply(sigma_pn, L_Rho_pn) * z_pn)';
pi = pn[,1];
nu = pn[,2];
Rho_pn = L_Rho_pn * L_Rho_pn';


ab = (diag_pre_multiply(sigma_ab, L_Rho_ab) * z_ab)';
for (i in 1:N)
{
  ab[,1] = (mu_gamma + gt[sp[i], 1]) +
  (mu_theta + gt[sp[i], 2]) * cell_lat[i];
  
  ab[i,2] = (mu_pi + pn[sp[i], 1]) +
  (mu_nu + pn[sp[i], 2]) * cell_lat[i];
}

alpha = ab[,1];
beta = ab[,2];
Rho_ab = L_Rho_ab * L_Rho_ab';


// alpha_c = mu_alpha + alpha;
// beta_c = mu_beta + beta;
gamma_c = mu_gamma + gamma;
theta_c = mu_theta + theta;
pi_c = mu_pi + pi;
nu_c = mu_nu + nu;

for (i in 1:N)
{
  mu[i] = alpha[i] + beta[i] * year[i];
}

// implies mu_y[i] ~ normal(mu[i], sigma)
mu_y = mu_y_raw * sigma + mu; 

}

model {
sigma_ab_raw ~ std_normal();
sigma_gt_raw ~ std_normal();
sigma_pn_raw ~ std_normal();
mu_gamma_raw ~ std_normal();
mu_theta_raw ~ std_normal();
mu_pi_raw ~ std_normal();
mu_nu_raw ~ std_normal();
sigma_raw ~ std_normal();
mu_y_raw ~ std_normal();

to_vector(z_ab) ~ std_normal();
L_Rho_ab ~ lkj_corr_cholesky(1);
to_vector(z_gt) ~ std_normal();
L_Rho_gt ~ lkj_corr_cholesky(1);
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

DELTA <- 0.90
TREE_DEPTH <- 16
STEP_SIZE <- 0.0001
CHAINS <- 4
ITER <- 4000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel1,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('gamma_c',
                            'theta_c',
                            'pi_c',
                            'nu_c',
                            'gamma',
                            'theta',
                            'pi',
                            'nu',
                            'mu_gamma',
                            'mu_theta',
                            'mu_pi',
                            'mu_nu',
                            'sigma_gt',
                            'sigma_pn',
                            'Rho_gt',
                            'Rho_pn',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_', juv_date))
saveRDS(fit, file = paste0('juv-time-bad-stan-output-', juv_date, '.rds'))
#fit <- readRDS(paste0('arr-br-juv-stan-output-', juv_date, '.rds'))

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
sink(paste0('juv-time-stan-results-', juv_date, '.txt'))
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
                   filename = paste0('juv-time-', juv_date, '-trace_mu_alpha.pdf'))

#mu_beta ~ N(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_beta',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', juv_date, '-trace_mu_beta.pdf'))

#sigma_sp[1] ~ HN(0, 40)
PR_p <- rnorm(10000, 0, 40)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sp\\[1',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', juv_date, '-trace_sigma_sp[1].pdf'))

#sigma_sp[2] ~ HN(0, 1)
PR_p <- rnorm(10000, 0, 1)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_sp\\[2',
                   ISB = 'FALSE',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', juv_date, '-trace_sigma_sp[2].pdf'))


#sigma ~ HN(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', juv_date, '-trace_sigma.pdf'))

#mu_arr ~ N(180, 40)
PR <- rnorm(10000, 180, 40)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_arr',
                   priors = PR,
                   pdf = TRUE,
                   filename = paste0('juv-time-', juv_date, '-trace_mu_juv.pdf'))



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

