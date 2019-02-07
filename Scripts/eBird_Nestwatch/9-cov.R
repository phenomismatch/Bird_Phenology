####################
# 9 - nesting date (from Nestwatch and others) ~ nesting date (from IAR model)
#
# *How well do arrival dates (as determined from the IAR model) predict nesting dates (as determined from Nestwatch)
####################


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



# nesting date ~ arrival date ---------------------------------------------

#arrival date (from IAR model) to nesting date (from BR codes and possible NW and MAPS)

#remove na vals for BR cods
to.rm <- which(is.na(mdf$EB_HM_mean))
mdf2 <- mdf[-to.rm, ]

#remove species where there are less than 3 obs
obs_per_sp <- plyr::count(mdf2, 'species')
sp.to.rm <- obs_per_sp[which(obs_per_sp[,2] < 3), 1]
mdf3 <- mdf2[-which(mdf2$species %in% sp.to.rm),]

#number codes for species
sp_num <- as.numeric(factor(mdf3$species))

#create data list for Stan
DATA <- list(y_obs = mdf3$EB_HM_mean,
             sigma_y = mdf3$EB_HM_sd,
             x_obs = mdf3$mean_post_IAR,
             sigma_x = mdf3$sd_post_IAR,
             sp_id = sp_num,
             US = length(unique(sp_num)),
             N = NROW(mdf3))




# Stan model --------------------------------------------------------------

br_arr_corr <- "
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
vector<lower = 0, upper = 300>[N] y_obs;                // mean halfmax BR codes
vector<lower = 0>[N] sigma_y;                           // sd halfmax BR codes
vector<lower = 0, upper = 200>[N] x_obs;                // mean halfmax IAR
vector<lower = 0>[N] sigma_x;                           // sd halfmax IAR
int<lower = 1, upper = US> sp_id[N];                  // species ids
}

parameters {
vector<lower = 0, upper = 300>[N] y_true;                           //true arrival date
vector<lower = 0, upper = 200>[N] x_true;                           //true nesting date
real mu_alpha_raw;
real mu_beta_raw;
real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_beta_raw;
real<lower = 0> sigma_raw;
real<lower = -1, upper = 1> rho;
vector[2] B_temp;
}


model {

vector[N] mu;
matrix[2, 2] Sigma_b;
real mu_alpha;
real mu_beta;
real sigma_alpha;
real sigma_beta;
real sigma;
matrix[2, US] B_hat;
vector[2] B_hat_temp;
matrix[2, US] B;
vector[US] alpha;
vector[US] beta;


// observation model - modeling true state as a function of some observed state

y_obs ~ normal(y_true, sigma_y);
x_obs ~ normal(x_true, sigma_x);

// priors - non-centered parameterization

mu_alpha_raw ~ normal(0, 1);
mu_beta_raw ~ normal(0, 1);
sigma_alpha_raw ~ normal(0, 1);
sigma_beta_raw ~ normal(0, 1);
sigma_raw ~ normal(0, 1);
rho ~ uniform(-1, 1);


// non-centered parameterization

mu_alpha = mu_alpha_raw * 20 + 70;                       // implies mu_alpha ~ normal(70, 20)
mu_beta = mu_beta_raw * 2 + 1;                           // implies mu_beta ~ normal(1, 2)
sigma_alpha = sigma_alpha_raw * 10;                      // implies sigma_alpha ~ halfnormal(0, 10)
sigma_beta = sigma_beta_raw * 3;                         // implies sigma_beta ~ halfnormal(0, 3)
sigma = sigma_raw * 10;                                  // implies sigma ~ halfnormal(0, 10)

// fill Sigma_b

Sigma_b[1, 1] = pow(sigma_alpha, 2);
Sigma_b[2, 2] = pow(sigma_beta, 2);
Sigma_b[1, 2] = rho * sigma_alpha * sigma_beta;
Sigma_b[2, 1] = Sigma_b[1, 2];


// draw alpha and beta from multivariate N

for (j in 1:US)
{
  B_hat[1,j] = mu_alpha;
  B_hat[2,j] = mu_beta;
  B_hat_temp[1] = mu_alpha;
  B_hat_temp[2] = mu_beta;
  B_temp ~ multi_normal(B_hat_temp, Sigma_b);
  B[1, j] = B_temp[1];
  B[2, j] = B_temp[2];
}

for (j in 1:US)
{
  alpha[j] = B[1, j];
  beta[j] = B[2, j];
}

for (i in 1:N)
{
  mu[i] = alpha[sp_id[i]] + beta[sp_id[i]] * x_true[i];
}

y_true ~ normal(mu, sigma);

}

generated quantities {

// vector[N] resids;
// real BR2;
// vector[N] y_rep;
// real PPC_mean;
// vector[N] arr_br;
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
// arr_br = x_true - y_true;
}
"



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


DELTA <- 0.95
TREE_DEPTH <- 15
STEP_SIZE <- 0.005
CHAINS <- 1
ITER <- 30

tt <- proc.time()
fit <- rstan::stan(model_code = br_arr_corr,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 'beta', 'mu_alpha', 'mu_beta', 'sigma_alpha', 'sigma_beta', 
                            'rho', 'sigma', 'y_true', 'x_true'),
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



#all species together
data_vis_fun(SPECIES = 'all')

#each species individually
sps <- unique(mdf3$species)
for (i in 1:length(sps))
{
  #i <- 3
  data_vis_fun(SPECIES = sps[i])
}

