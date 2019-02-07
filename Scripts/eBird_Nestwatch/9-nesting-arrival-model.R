####################
# 9 - nesting date (from Nestwatch and others) ~ nesting date (from IAR model)
#
# *How well do arrival dates (as determined from the IAR model) predict nesting dates (as determined from Nestwatch)
####################


# data fusion model for arrival ~ breeding --------------------------------


#each species different model? Each species own intercept/slope? How to accomodate missing data from entire species (MAPS data might not be available for everyone)?

#need to determine how to assign uncertainty for NW data (maybe have one sigma for NW, where that is informed by sigma values for each cell obs - sd tend to be large for this dataset, though some cells only have one obs so no sd [e.g., perfect detection])

#obs_EB_breeding_date is halfmax estimate for EB breeding date for that species/cell/year
#known_EB_uncertainty is uncertainty in halfmax estimate for EB breeding date for that species/cell/year
#obs_NW_breeding_date is mean NW breeding date for that species/cell/year
#known_NW_uncertainty is sd of NW breeding date for that species/cell/year
#obs_MAPS_breeding_date is mean of breeding date period for that species/cell/year
#known_MAPS_uncertainty is length of period (uniform distribution)

###observation models
#obs_IAR_arrival_date[i] ~ N(true_IAR_arrival_date[i], known_IAR_uncertainty[i])
#obs_EB_breeding_date[i] ~ N(true_EB_breeding_date[i], known_EB_uncertainty[i])
#obs_NW_breeding_date[i] ~ N(true_NW_breeding_date[i], known_NW_uncertainty[i]) #index so that positions without data aren't modeled here
#known_NW_uncertainty[i] ~ N(mu_NW_sigma, sigma_NW_sigma) #to accomodate missing values in covariate - See Statisitical Rethinking Ch. 14; can accomodate relationship of this with other predictor as well, using: mu_NW_sigma[i] = alpha + beta*true_MAPS_breeding_date[i]
#true_MAPS_breeding_date[i] ~ U(lower_bounds_period[i], upper_bounds_period[i]) #informed by obs_MAPS_breeding_date; true in as NAs; assuming that MAPS obs for most EB obs - this true?


#ONE WAY (arrival ~ breeding):
###process model for explanatory var
#true_EB_breeding_date ~ N(master_breeding_date, sigma_EB)
#master_breeding_date = alpha_EB + beta1_EB * true_NW_breeding_date + beta2_EB * true_MAPS_breeding date

#OR

#true_EB_breeding_date ~ N(master_breeding_date, sigma_EB)
#true_NW_breeding_date ~ N(master_breeding_date, sigma_NW) #NAs in places where there are EB data but not NW data
#master_breeding_date = alpha_EB + beta2_EB * true_MAPS_breeding date

###process model for response var
#true_IAR_arrival_date ~ N(mu, sigma)
#mu = alpha + beta * master_breeding_date


#ALTERNATIVELY (breeding ~ arrival):
#true_EB_breeding_date ~ N(master_breeding_date, sigma_EB)
#master_breeding_date = alpha_EB + beta1_EB * true_NW_breeding_date + beta2_EB * true_MAPS_breeding date

#OBS_EB_breeding_date ~ N(master_breeding_date, known_EB_uncertainty)
#OBS_NW_breeding_date ~ N(master_breeding_date, known_NW_uncertainty)
#OBS_MAPS_breeding_date ~ N(master_breeding_date, known_MAPS_uncertainty)
#master_breeding_date = alpha_EB + beta1_EB * true_NW_breeding_date + beta2_EB * true_MAPS_breeding date

#OR (though may not be possible bc MAPS is a uniform distribution)

#OBS_EB_breeding_date ~ N(master_breeding_date, known_EB_uncertainty)
#OBS_NW_breeding_date ~ N(master_breeding_date, known_NW_uncertainty)
#OBS_MAPS_breeding_date ~ N(master_breeding_date, known_MAPS_uncertainty)

#master_breeding_date ~ N(mu, sigma)
#mu = alpha + beta * true_IAR_arrival_date


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



# ARR - BR ----------------------------------------------------------------

# #should calculate this in Stan as a derived qty
# mdf$diff_arr_br <- mdf$mean_post_IAR - mdf$EB_HM_mean
# 
# #this could then be used in a separate model to look at how this differences has changed over time, and what enviromental factors lead to years with larger differences between arrival and breeding
# 
# plot(mdf$year, mdf$diff_arr_br, pch = 19, col = rgb(0,0,0,0.5))
# summary(lm(mdf$diff_arr_br ~ mdf$year))


# plot data ---------------------------------------------------------------

# head(mdf)
# plot(mdf$mean_post_IAR, mdf$EB_HM_mean, pch = 19, col = rgb(0,0,0,0.5),
#      xlim = c(50, 220), ylim = c(50, 220))
# abline(a = 0, b = 1, lty = 2, col = 'red')
# summary(lm(mdf$EB_HM_mean ~ mdf$mean_post_IAR))
# 
# 
# temp <- dplyr::filter(mdf, species == 'Vireo_olivaceus')
# plot(temp$mean_post_IAR, temp$EB_HM_mean, pch = 19, col = rgb(0,0,0,0.5),
#      xlim = c(50, 220), ylim = c(50, 220))
# abline(a = 0, b = 1, lty = 2, col = 'red')
# summary(lm(temp$EB_HM_mean ~ temp$mean_post_IAR))



#OVERVIEW
#=======#
#breeding date ~ arrival date
#breeding date ~ year + lat
#arrival date ~ year + lat


#DETAILS
#======#
# #observation model - where obs and sigma are known
# y_obs[i] ~ normal(y_true[i], sigma_y[i])
# x_obs[i] ~ normal(x_true[i], sigma_x[i])
# 
# #breeding as a function of arrival
# y_true[i] ~ normal(mu[i], sigma)
# mu[i] = alpha_j + beta_j * true_arr[i]
# 
# #breeding as a function of year and latitude
# y_true[i] ~ normal(mu_br[i], sigma_br)
# mu_br[i] = alpha_br_j + beta1_br_j * year[i] + beta2_br_j * lat[i]
# 
# #arrival as a function of year and latitude
# x_true[i] ~ normal(mu_arr[i], sigma_arr)
# mu_arr[i] = alpha_arr_j + beta1_arr_j * year[i] + beta2_arr_j * lat[i]
# 
# #derived qty - change in difference over time and space for each species
# diff[i] = true_arr[i] - true_br[i]


#FIGS
#======#
#[Fig 1a - for each species: plot year on x, DOY on y - points for arrival and breeding date each year - slope for change over year] - might not be able to plot slope here for multiple regression - see McElreath
#[Fig 1b - for each species: plot year on x, diff on y - points for diff each year]
#[Fig 2 - for each species: plot lat on x, DOY on y - points for arrival and breeding date each lat - slope for change over lat]



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

#create data list for Stan
DATA <- list(y_obs = mdf3$EB_HM_mean,
             sigma_y = mdf3$EB_HM_sd,
             x_obs = mdf3$mean_post_IAR,
             sigma_x = mdf3$sd_post_IAR,
             year = yr_num,
             lat = mdf3$cell_lat,
             sp_id = sp_num,
             US = length(unique(sp_num)),
             N = NROW(mdf3))



# Stan model --------------------------------------------------------------

br_arr <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
vector<lower = 0, upper = 300>[N] y_obs;              // mean halfmax BR codes
vector<lower = 0>[N] sigma_y;                         // sd halfmax BR codes
vector<lower = 0, upper = 200>[N] x_obs;              // mean halfmax IAR
vector<lower = 0>[N] sigma_x;                         // sd halfmax IAR
int<lower = 1, upper = US> sp_id[N];                  // species ids
vector<lower = 1, upper = 7>[N] year;                 // year number
vector<lower = 26, upper = 90>[N] lat;                // lat of hex cell
}

parameters {
vector<lower = 0, upper = 300>[N] y_true;             // true arrival date
vector<lower = 0, upper = 200>[N] x_true;             // true nesting date

real mu_alpha_raw;
real mu_beta_raw;
real mu_alpha_br_raw;
real mu_beta1_br_raw;
real mu_beta2_br_raw;
real mu_alpha_arr_raw;
real mu_beta1_arr_raw;
real mu_beta2_arr_raw;

real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_beta_raw;
real<lower = 0> sigma_alpha_br_raw;
real<lower = 0> sigma_beta1_br_raw;
real<lower = 0> sigma_beta2_br_raw;
real<lower = 0> sigma_alpha_arr_raw;
real<lower = 0> sigma_beta1_arr_raw;
real<lower = 0> sigma_beta2_arr_raw;

real<lower = 0> sigma_raw;
real<lower = 0> sigma_br_raw;
real<lower = 0> sigma_arr_raw;

real alpha_raw[US];
real beta_raw[US];
real alpha_br_raw[US];
real beta1_br_raw[US];
real beta2_br_raw[US];
real alpha_arr_raw[US];
real beta1_arr_raw[US];
real beta2_arr_raw[US];
}

transformed parameters {

real mu_alpha;
real mu_beta;
real sigma_alpha;
real sigma_beta;
real mu_alpha_br;
real mu_beta1_br;
real mu_beta2_br;
real sigma_alpha_br;
real sigma_beta1_br;
real sigma_beta2_br;
real mu_alpha_arr;
real mu_beta1_arr;
real mu_beta2_arr;
real sigma_alpha_arr;
real sigma_beta1_arr;
real sigma_beta2_arr;

real sigma;
real sigma_br;
real sigma_arr;

real alpha[US];
real beta[US];
real alpha_br[US];
real beta1_br[US];
real beta2_br[US];
real alpha_arr[US];
real beta1_arr[US];
real beta2_arr[US];

real mu[N];
real mu_br[N];
real mu_arr[N];

// non-centered parameterization

mu_alpha = mu_alpha_raw * 20 + 70;                       // implies mu_alpha ~ normal(70, 20)
mu_beta = mu_beta_raw * 2 + 1;                           // implies mu_beta ~ normal(1, 2)
sigma_alpha = sigma_alpha_raw * 10;                      // implies sigma_alpha ~ halfnormal(0, 10)
sigma_beta = sigma_beta_raw * 3;                         // implies sigma_beta ~ halfnormal(0, 3)

mu_alpha_br = mu_alpha_br_raw * 20 + 70;
mu_beta1_br = mu_beta1_br_raw * 2 + 1;
mu_beta2_br = mu_beta2_br_raw * 2 + 1;
sigma_alpha_br = sigma_alpha_br_raw * 10;
sigma_beta1_br = sigma_beta1_br_raw * 3;
sigma_beta2_br = sigma_beta2_br_raw * 3;

mu_alpha_arr = mu_alpha_arr_raw * 20 + 70;
mu_beta1_arr = mu_beta1_arr_raw * 2 + 1;
mu_beta2_arr = mu_beta2_arr_raw * 2 + 1;
sigma_alpha_arr = sigma_alpha_arr_raw * 10;
sigma_beta1_arr = sigma_beta1_arr_raw * 3;
sigma_beta2_arr = sigma_beta2_arr_raw * 3;


sigma = sigma_raw * 10;                                  // implies sigma ~ halfnormal(0, 10)
sigma_br = sigma_br_raw * 10;
sigma_arr = sigma_arr_raw * 10;


for (j in 1:US)
{
  alpha[j] = alpha_raw[j] * sigma_alpha + mu_alpha;      // implies alpha[j] ~ normal(mu_alpha, sigma_alpha)
  beta[j] = beta_raw[j] * sigma_beta + mu_beta;          // implies beta[j] ~ normal(mu_beta, sigma_beta)

  alpha_br[j] = alpha_br_raw[j] * sigma_alpha_br + mu_alpha_br;
  beta1_br[j] = beta1_br_raw[j] * sigma_beta1_br + mu_beta1_br;
  beta2_br[j] = beta2_br_raw[j] * sigma_beta2_br + mu_beta2_br;

  alpha_arr[j] = alpha_arr_raw[j] * sigma_alpha_arr + mu_alpha_arr;
  beta1_arr[j] = beta1_arr_raw[j] * sigma_beta1_arr + mu_beta1_arr;
  beta2_arr[j] = beta2_arr_raw[j] * sigma_beta2_arr + mu_beta2_arr;
}

for (i in 1:N)
{
  mu[i] = alpha[sp_id[i]] + beta[sp_id[i]] * x_true[i];
  mu_br[i] = alpha_br[sp_id[i]] + beta1_br[sp_id[i]] * year[i] + beta2_br[sp_id[i]] * lat[i];
  mu_arr[i] = alpha_arr[sp_id[i]] + beta1_arr[sp_id[i]] * year[i] + beta2_arr[sp_id[i]] * lat[i];
}

}

model {

// observation model - modeling true state as a function of some observed state

y_obs ~ normal(y_true, sigma_y);
x_obs ~ normal(x_true, sigma_x);

// non-centered parameterization

mu_alpha_raw ~ normal(0, 1);
mu_beta_raw ~ normal(0, 1);
sigma_alpha_raw ~ normal(0, 1);
sigma_beta_raw ~ normal(0, 1);

mu_alpha_br_raw ~ normal(0, 1);
mu_beta1_br_raw ~ normal(0, 1);
mu_beta2_br_raw ~ normal(0, 1);
sigma_alpha_br_raw ~ normal(0, 1);
sigma_beta1_br_raw ~ normal(0, 1);
sigma_beta2_br_raw ~ normal(0, 1);

mu_alpha_arr_raw ~ normal(0, 1);
mu_beta1_arr_raw ~ normal(0, 1);
mu_beta2_arr_raw ~ normal(0, 1);
sigma_alpha_arr_raw ~ normal(0, 1);
sigma_beta1_arr_raw ~ normal(0, 1);
sigma_beta2_arr_raw ~ normal(0, 1);

sigma_raw ~ normal(0, 1);
sigma_br_raw ~ normal(0, 1);
sigma_arr_raw ~ normal(0, 1);

for (j in 1:US)
{
  alpha_raw[j] ~ normal(0, 1);
  beta_raw[j] ~ normal(0, 1);

  alpha_br_raw[j] ~ normal(0, 1);
  beta1_br_raw[j] ~ normal(0, 1);
  beta2_br_raw[j] ~ normal(0, 1);

  alpha_arr_raw[j] ~ normal(0, 1);
  beta1_arr_raw[j] ~ normal(0, 1);
  beta2_arr_raw[j] ~ normal(0, 1);
}

y_true ~ normal(mu, sigma);
y_true ~ normal(mu_br, sigma_br);
x_true ~ normal(mu_arr, sigma_arr);

}

generated quantities {

// vector[N] resids;
// real BR2;
// vector[N] y_rep;
// real PPC_mean;
vector[N] arr_br;
// real var_mu;
// real var_resids;

// #traditional R^2
// RSS = dot_self(y - mu);
// TSS = dot_self(y - mean(y));
// R2 = 1 - RSS/TSS;

// new Bayes R^2 - http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf
// calculate residuals variance for mu and errors to use in R^2
// resids = y_true - mu;
// var_mu = (dot_self(mu - mean(mu))) / (N - 1);
// var_resids = (dot_self(resids - mean(resids))) / (N - 1);
// BR2 = var_mu/(var_mu + var_resids);

// PPC
//for (i in 1:N)
//{
//  y_rep[N] = normal_rng(mu[N], sigma);
//}
//PPC_mean = mean(y_rep);

// arrival date - breeding date
arr_br = x_true - y_true;
}
'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


DELTA <- 0.97
TREE_DEPTH <- 18
STEP_SIZE <- 0.0005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit <- rstan::stan(model_code = br_arr,
            data = DATA,
            chains = CHAINS,
            iter = ITER,
            cores = CHAINS,
            pars = c('alpha', 'beta', 
                     'alpha_br', 'beta1_br', 'beta2_br',
                     'alpha_arr', 'beta1_arr', 'beta2_arr',
                     'mu_alpha', 'mu_beta', 
                     'sigma_alpha', 'sigma_beta', 
                     'mu_alpha_br', 'mu_beta1_br', 'mu_beta2_br',
                     'sigma_alpha_br', 'sigma_beta1_br', 'sigma_beta2_br', 
                     'mu_alpha_arr', 'mu_beta1_arr', 'mu_beta2_arr',
                     'sigma_alpha_arr', 'sigma_beta1_arr', 'sigma_beta2_arr', 
                     'sigma', 'sigma_br', 'sigma_arr',
                     'y_true', 'x_true', 'arr_br'),
            control = list(max_treedepth = TREE_DEPTH, adapt_delta = DELTA, stepsize = STEP_SIZE)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit, file = paste0('temp_BR_ARR_stan_', MODEL_DATE, '.rds'))
#fit <- readRDS(paste0('temp_BR_ARR_stan_', MODEL_DATE, '.rds'))




# diagnostics -------------------------------------------------------------


MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = c('alpha', 'beta'), ISB = FALSE)
MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'sigma')
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'y_true')
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'x_true')
#MCMCtrace(fit)

(num_diverge <- rstan::get_num_divergent(fit))
(num_tree <- rstan::get_num_max_treedepth(fit))
(num_BFMI <- rstan::get_low_bfmi_chains(fit))


#shiny stan
# library(shinystan)
# launch_shinystan(fit)



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





# Plot results ------------------------------------------------------------


data_vis_fun <- function(SPECIES = 'all')
{
  #SPECIES <- 'Vireo_olivaceus'
  
  #extract posterior estimates for true states for y and x
  y_true_mean <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'summary', 
                                   func = mean)[[1]]
  y_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.025)))[[1]]
  y_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'summary', 
                                  func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  x_true_mean <- MCMCvis::MCMCpstr(fit, params = 'x_true', type = 'summary', func = mean)[[1]]
  x_true_LCI <- MCMCvis::MCMCpstr(fit, params = 'x_true', type = 'summary', func = function(x) quantile(x, probs = c(0.025)))[[1]]
  x_true_UCI <- MCMCvis::MCMCpstr(fit, params = 'x_true', type = 'summary', func = function(x) quantile(x, probs = c(0.975)))[[1]]
  
  #need true latent states
  DATA_PLOT <- data.frame(mean_y = y_true_mean,
                           mean_y_l = y_true_LCI,
                           mean_y_u = y_true_UCI,
                           mean_x = x_true_mean, 
                           mean_x_l = x_true_LCI,
                           mean_x_u = x_true_UCI,
                           sp_id = sp_num)
  
  if (SPECIES == 'all')
  {
    #model fit for mu_beta and mu_alpha
    alpha_ch <- MCMCchains(fit, params = 'mu_alpha')[,1]
    beta_ch <- MCMCchains(fit, params = 'mu_beta')[,1]
    
    DATA_PLOT2 <- DATA_PLOT
  } else {
    idx <- which(unique(mdf3$species) == SPECIES)
    if (length(idx) > 0)
    {
      alpha_ch <- MCMCchains(fit, params = paste0('alpha\\[', idx, '\\]'), ISB = FALSE)[,1]
      beta_ch <- MCMCchains(fit, params = paste0('beta\\[', idx, '\\]'), ISB = FALSE)[,1]
    
      DATA_PLOT2 <- dplyr::filter(DATA_PLOT, sp_id == idx)
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
    geom_point(data = DATA_PLOT2, aes(mean_x, mean_y), color = 'black',
               inherit.aes = FALSE, size = 1, alpha = 0.3) +
    theme_bw() +
    #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
    xlab('True ARR halfmax') +
    ylab('True BR halfmax') +
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







#vvvv HERE DOWN NOT IN OTHER SCRIPT THAT ACCOUNTS FOR CORR BETWEEN ALPHA AND BETA vvvv
# change is arr_br over time ----------------------------------------------

#Diff_arr_br ~ alpha_j + beta1_j * year + beta2_j * arr #where j is species

#change to model correction between alpha and betas - see wishart model stan


#read in RDS
#setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
#fit <- readRDS(paste0('temp_BR_ARR_stan_', MODEL_DATE, '.rds'))


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


