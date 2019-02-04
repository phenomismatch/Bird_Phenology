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
library(dggridR)
library(ggplot2)
library(rstan)
library(MCMCvis)


# import ARR/BR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

#arrival data
IAR_out <- 'IAR_output_2019-01-16'
IAR_out_date <- substr(IAR_out, start = 12, stop = 21)

IAR_data <- readRDS(paste0('arrival_master_', IAR_out_date, '.rds'))


#breeding data
BR_out <- 'halfmax_breeding_2019-01-30'
BR_out_date <- substr(BR_out, start = 18, stop = 27)

BR_data <- readRDS(paste0('temp_breeding_master_', BR_out_date, '.rds'))




# # import IAR species list -----------------------------------------------------
# 
# setwd(paste0(dir, 'Bird_Phenology/Data/'))
# 
# species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)
# 
# #remove underscore and coerce to vector
# species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))



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

#number codes for species
sp_num <- as.numeric(factor(mdf2$species))

#create data list for Stan
DATA <- list(y_obs = mdf2$EB_HM_mean,
             sigma_y = mdf2$EB_HM_sd,
             x_obs = mdf2$mean_post_IAR,
             sigma_x = mdf2$sd_post_IAR,
             sp_id = sp_num,
             US = length(unique(sp_num)),
             N = NROW(mdf2))



# Stan model --------------------------------------------------------------

br_arr <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
real<lower = 0, upper = 200> y_obs[N];                // mean halfmax BR codes
real<lower = 0> sigma_y[N];                           // sd halfmax BR codes
real<lower = 0, upper = 200> x_obs[N];                // mean halfmax IAR
real<lower = 0> sigma_x[N];                           // sd halfmax IAR
int<lower = 1, upper = US> sp_id[N];                  // species ids
}

parameters {
real<lower = 0> y_true[N];                           //true arrival date
real<lower = 0> x_true[N];                           //true nesting date
real<lower = 0> sigma;
real alpha[US];
real beta[US];
real mu_alpha;
real mu_beta;
real<lower = 0> sigma_alpha;
real<lower = 0> sigma_beta;
}

transformed parameters {

//one alpha, one beta for each species
real mu[N];

for (i in 1:N)
{
  mu[i] = alpha[sp_id[i]] + beta[sp_id[i]] * x_true[i];
}

}

model {

// observation model - modeling true state as a function of some observed state

y_obs ~ normal(y_true, sigma_y);
x_obs ~ normal(x_true, sigma_x);

for (j in 1:US)
{
  alpha[j] ~ normal(mu_alpha, sigma_alpha);
  beta[j] ~ normal(mu_beta, sigma_beta);
}

y_true ~ normal(mu, sigma);

}'



# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

tt <- proc.time()
fit <- stan(model_code = br_arr,
            data = DATA,
            chains = 4,
            iter = 1000,
            cores = 4,
            pars = c('y_true', 'x_true', 'sigma', 'alpha', 'beta', 'mu_alpha',
                     'mu_beta', 'sigma_alpha', 'sigma_beta'),
            control = list(max_treedepth = 18, adapt_delta = 0.95))#, stepsize = 0.005)) # modified control parameters based on warnings
run_time <- (proc.time() - tt[3]) / 60


MCMCsummary(fit)




#save to RDS
# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
# saveRDS(fit, file = paste0('IAR_stan_', args, '-', IAR_out_date, '.rds'))
# fit <- readRDS('IAR_stan_Catharus_minimus-2018-11-12.rds')




# diagnostics -------------------------------------------------------------


# pairs(fit, pars = c('sigma', 'rho'))

# sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
# mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
# max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
# get_elapsed_time(fit)

# MCMCtrace(fit)
# MCMCsummary(fit, params = c('sigma', 'rho', 'beta0'), n.eff = TRUE)
# MCMCsummary(fit, params = c('theta', 'phi'), n.eff = TRUE)

# print(fit, pars = c('sigma', 'rho'))

# #shiny stan
# library(shinystan)
# launch_shinystan(fit)




# write model results to file ---------------------------------------------

# options(max.print = 50000)
# sink(paste0('IAR_results_', args, '.txt'))
# cat(paste0('IAR results ', args, ' \n'))
# cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
# print(fit)
# sink()



# Plot results ------------------------------------------------------------

DATA <- list(y_obs = mdf2$EB_HM_mean,
             sigma_y = mdf2$EB_HM_sd,
             x_obs = mdf2$mean_post_IAR,
             sigma_x = mdf2$sd_post_IAR,
             sp_id = sp_num,
             US = length(unique(sp_num)),
             N = NROW(mdf2))


DATA_PLOT2 <- data.frame(mean_y = DATA$y_obs,
                         mean_y_l = (DATA$y_obs - DATA$simga_y),
                         mean_y_u = (DATA$y_obs + DATA$simga_y),
                         mean_x = DATA$x_obs, 
                         mean_x_l = (DATA$x_obs - DATA$simga_x),
                         mean_x_u = (DATA$x_obs - DATA$simga_x))

#model fit to plot - for one alpha/beta only (could plot top level beta parameter)
# alpha_ch <- MCMCchains(fit, params = 'alpha')[,1]
# beta_ch <- MCMCchains(fit, params = 'beta')[,1]
# 
# sim_med_c <- seq(min(med_c)-1, max(med_c)+1, length = 100)
# 
# mf <- matrix(nrow = length(beta_ch), ncol = 100)
# for (i in 1:length(sim_med_c))
# {
#   mf[,i] <- alpha_ch + beta_ch * sim_med_c[i]
# }
# 
# med_mf <- apply(mf, 2, median)
# LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
# UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
# 
# FIT_PLOT <- data.frame(MN = med_mf, 
#                        MN_X = sim_med_c,
#                        LCI = LCI_mf,
#                        UCI = UCI_mf)


ggplot(data = DATA_PLOT2, aes(mean_x, mean_y), color = 'black', alpha = 0.6) +
  # geom_ribbon(data = FIT_PLOT, 
  #             aes(x = MN_X, ymin = LCI, ymax = UCI),
  #             fill = 'grey', alpha = 0.7,
  #             inherit.aes = FALSE) +
  # geom_line(data = FIT_PLOT, aes(MN_X, MN), color = 'red',
  #           alpha = 0.9,
  #           inherit.aes = FALSE,
  #           size = 1.4) +
  geom_errorbar(data = DATA_PLOT2, 
                aes(ymin = mean_y_l, ymax = mean_y_u), width = 0.3,
                color = 'black', alpha = 0.2) +
  geom_errorbarh(data = DATA_PLOT2, 
                 aes(xmin = mean_x_l, xmax = mean_x_u), height = 0.005,
                 color = 'black', alpha = 0.2) +
  geom_point(data = DATA_PLOT2, aes(mean_x, mean_y), color = 'black',
             inherit.aes = FALSE, size = 3, alpha = 0.7) +
  theme_bw() +
  #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
  xlab('True BR halfmax') +
  ylab('True ARR halfmax') +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick







