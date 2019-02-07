####################
# 9 - breeding data ~ year + lat
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




#OVERVIEW
#=======#
#breeding date ~ year + lat


#DETAILS
#======#
# #observation model - where obs and sigma are known
# y_obs[i] ~ normal(y_true[i], sigma_y[i])
# 
# #arrival as a function of year and latitude
# y_true[i] ~ normal(mu_br[i], sigma_br)
# mu_br[i] = alpha_br_j + beta1_br_j * year[i] + beta2_br_j * lat[i]



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

#create data list for Stan
DATA <- list(y_obs = mdf3$EB_HM_mean,
             sigma_y = mdf3$EB_HM_sd,
             sp_id = sp_num,
             year = yr_num,
             lat = mdf3$cell_lat,
             US = length(unique(sp_num)),
             N = NROW(mdf3))



# Stan model --------------------------------------------------------------

br_time_lat <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
vector<lower = 0, upper = 200>[N] y_obs;                // mean halfmax IAR
vector<lower = 0>[N] sigma_y;                           // sd halfmax IAR
int<lower = 1, upper = US> sp_id[N];                  // species ids
vector<lower = 1, upper = 17>[N] year;
vector<lower = 26, upper = 90>[N] lat;
}

parameters {
vector<lower = 0, upper = 300>[N] y_true;                           //true arrival
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
  mu[i] = alpha[sp_id[i]] + beta1[sp_id[i]] * year[i] + beta2[sp_id[i]] * lat[i];
}
}

model {

// observation model - modeling true state as a function of some observed state
y_obs ~ normal(y_true, sigma_y);

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
fit <- rstan::stan(model_code = br_time_lat,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 'beta1', 'beta2', 
                            'mu_alpha', 'mu_beta1', 'mu_beta2',
                            'sigma_alpha', 'sigma_beta1', 'sigma_beta2', 
                            'sigma', 'y_true'),
                   control = list(max_treedepth = TREE_DEPTH, adapt_delta = DELTA, stepsize = STEP_SIZE)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit, file = paste0('temp_BR_YEAR_LAT_stan_', MODEL_DATE, '.rds'))
#fit <- readRDS(paste0('temp_BR_ARR_stan_', MODEL_DATE, '.rds'))




# diagnostics -------------------------------------------------------------


# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = c('alpha', 'beta'), ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'sigma')
# # MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'y_true')
# # MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'x_true')
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





# # Plot results ------------------------------------------------------------
# 
# 
# data_vis_fun <- function(SPECIES = 'all')
# {
#   #SPECIES <- 'Vireo_olivaceus'
#   
#   #extract posterior estimates for true states for y and x
#   y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'summary', 
#                                    func = mean)[[1]]
#   y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'summary', 
#                                   func = function(x) quantile(x, probs = c(0.025)))[[1]]
#   y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'summary', 
#                                   func = function(x) quantile(x, probs = c(0.975)))[[1]]
#   
#   x_true_mean <- MCMCvis::MCMCpstr(fit, params = 'x_true', type = 'summary', func = mean)[[1]]
#   x_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'x_true', type = 'summary', func = function(x) quantile(x, probs = c(0.025)))[[1]]
#   x_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'x_true', type = 'summary', func = function(x) quantile(x, probs = c(0.975)))[[1]]
#   
#   #need true latent states
#   DATA_PLOT <- data.frame(mean_y = y_true_mean,
#                           mean_y_l = y_true_LCI,
#                           mean_y_u = y_true_UCI,
#                           mean_x = x_true_mean, 
#                           mean_x_l = x_true_LCI,
#                           mean_x_u = x_true_UCI,
#                           sp_id = sp_num)
#   
#   if (SPECIES == 'all')
#   {
#     #model fit for mu_beta and mu_alpha
#     alpha_ch <- MCMCchains(fit, params = 'mu_alpha')[,1]
#     beta_ch <- MCMCchains(fit, params = 'mu_beta')[,1]
#     
#     DATA_PLOT2 <- DATA_PLOT
#   } else {
#     idx <- which(unique(mdf3$species) == SPECIES)
#     if (length(idx) > 0)
#     {
#       alpha_ch <- MCMCchains(fit, params = paste0('alpha\\[', idx, '\\]'), ISB = FALSE)[,1]
#       beta_ch <- MCMCchains(fit, params = paste0('beta\\[', idx, '\\]'), ISB = FALSE)[,1]
#       
#       DATA_PLOT2 <- dplyr::filter(DATA_PLOT, sp_id == idx)
#     } else {
#       stop(paste0('Species: ', SPECIES, ' not found!'))
#     }
#   }
#   
#   sim_x <- seq(min(DATA_PLOT2$mean_x_l) - 1, max(DATA_PLOT2$mean_x_u) + 1, length = 100)
#   
#   mf <- matrix(nrow = length(beta_ch), ncol = 100)
#   for (i in 1:length(sim_x))
#   {
#     mf[,i] <- alpha_ch + beta_ch * sim_x[i]
#   }
#   
#   med_mf <- apply(mf, 2, median)
#   LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
#   UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
#   
#   FIT_PLOT <- data.frame(MN = med_mf,
#                          MN_X = sim_x,
#                          LCI = LCI_mf,
#                          UCI = UCI_mf)
#   
#   p <- ggplot(data = DATA_PLOT2, aes(mean_x, mean_y)) +
#     geom_ribbon(data = FIT_PLOT,
#                 aes(x = MN_X, ymin = LCI, ymax = UCI),
#                 fill = 'grey', alpha = 0.7,
#                 inherit.aes = FALSE) +
#     geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
#               alpha = 0.9,
#               inherit.aes = FALSE,
#               size = 1.4) +
#     geom_errorbar(data = DATA_PLOT2, 
#                   aes(ymin = mean_y_l, ymax = mean_y_u), width = 0.3,
#                   color = 'black', alpha = 0.2) +
#     geom_errorbarh(data = DATA_PLOT2, 
#                    aes(xmin = mean_x_l, xmax = mean_x_u), height = 0.005,
#                    color = 'black', alpha = 0.2) +
#     geom_point(data = DATA_PLOT2, aes(mean_x, mean_y), color = 'black',
#                inherit.aes = FALSE, size = 1, alpha = 0.3) +
#     theme_bw() +
#     #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
#     xlab('True ARR halfmax') +
#     ylab('True BR halfmax') +
#     ggtitle(paste0('Species: ', SPECIES)) +
#     theme(
#       plot.title = element_text(size = 22),
#       axis.text = element_text(size = 16),
#       axis.title = element_text(size = 18),
#       axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
#       axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
#       axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
#   
#   print(p)
# }
# 
# 
# 
# #all species together
# data_vis_fun(SPECIES = 'all')
# 
# #each species individually
# sps <- unique(mdf3$species)
# for (i in 1:length(sps))
# {
#   #i <- 3
#   data_vis_fun(SPECIES = sps[i])
# }
