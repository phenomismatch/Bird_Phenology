#######
# phenology ~ time using post IAR
# Vireo olivaceus
#######


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)
library(dggridR)



# read in data produced by 5 ----------------------------------------------------

setwd('~/Desktop/')

data <- readRDS('Vireo_olivaceus_pro_IAR.rds')



# Process data ------------------------------------------------------------

#cell years with input data
data_f <- data[which(!is.na(data$mean_pre_IAR)),]

#cells with more than three years of data
cnts <- plyr::count(data_f, 'cell')
u_cells <- cnts[which(cnts[,2] > 3),1]

data_f2 <- dplyr::filter(data_f, cell %in% u_cells)

t_cl <- unique(data_f2[,c('cell', 'cell_lat')])
ot_cl <- t_cl[order(t_cl[,1]),]

#create data list for Stan
DATA <- list(N = NROW(data_f2),
             y_obs = data_f2$mean_post_IAR,
             y_sd = data_f2$sd_post_IAR,
             cn_id = as.numeric(factor(data_f2$cell)),
             NC = length(u_cells),
             year = (data_f2$year - 2001),
             lat = ot_cl$cell_lat)


# ggplot(data_f2, aes(x = year, y = mean_post_IAR, col = factor(cell))) + 
#   geom_point() +
#   stat_smooth(method = 'lm')

# Stan model --------------------------------------------------------------

model <- "
data {
int<lower = 0> N;                                     // number of obs
vector<lower = 0, upper = 200>[N] y_obs;
vector<lower = 0>[N] y_sd;
int<lower = 1> cn_id[N];                  // species ids
int<lower = 0> NC;                                    // number of cells
vector<lower = 1, upper = 17>[N] year;
vector[NC] lat;
}

parameters {
vector[NC] alpha_raw;
vector[NC] beta_raw;
vector[N] y_true_raw;
// real<lower = 0> sigma_raw;
// vector<lower = 0>[NC] sigma_raw;
real<lower = 0> sigma_beta_raw;
real mu_beta_raw;
real alpha_beta_raw;
real beta_beta_raw;
real<lower = 0> sigma_sigma_raw;
real alpha_sigma_raw;
real beta_sigma_raw;
vector<lower = 0> [NC] sigma;
// real mu_alpha_raw;
// real<lower = 0> sigma_alpha_raw;
}

transformed parameters {
vector[N] mu;
vector[N] y_true;
vector[NC] alpha;
vector[NC] beta;
// real<lower = 0> sigma;
vector<lower = 0>[NC] mu_sigma;
real<lower = 0> sigma_sigma;
real<lower = 0> sigma_beta;
vector[NC] mu_beta;
real alpha_beta;
real beta_beta;
real<lower = 0> alpha_sigma;
real<lower = 0> beta_sigma;
// real mu_alpha;
// real<lower = 0> sigma_alpha;

// sigma = sigma_raw * 3;
alpha_sigma = alpha_sigma_raw * 3;
beta_sigma = beta_sigma_raw * 1;
mu_sigma = alpha_sigma + beta_sigma * lat;
sigma_sigma = sigma_sigma_raw * 1;

// alpha = alpha_raw * sigma_alpha + mu_alpha;
// mu_alpha = mu_alpha_raw * 1;
// sigma_alpha = sigma_alpha_raw * 1;
alpha = alpha_raw * 30 + 120;

alpha_beta = alpha_beta_raw * 3;
beta_beta = beta_beta_raw * 3;
mu_beta = alpha_beta + beta_beta * lat;
sigma_beta = sigma_beta_raw * 3;
beta = beta_raw * sigma_beta + mu_beta;

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
  // mu[i] = alpha[cn_id[i]] + beta * year[i]; // varying intercept, one slope
  y_true[i] = y_true_raw[i] * sigma[cn_id[i]] + mu[i];
}
}

model {

y_true_raw ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();

mu_beta_raw ~ std_normal();
sigma_beta_raw ~ std_normal();

// mu_alpha_raw ~ std_normal();
// sigma_alpha_raw ~ std_normal();

alpha_beta_raw ~ std_normal();
beta_beta_raw ~ std_normal();
alpha_sigma_raw ~ std_normal();
beta_sigma_raw ~ std_normal();

// sigma_raw ~ std_normal();

sigma ~ normal(mu_sigma, sigma_sigma);
y_obs ~ normal(y_true, y_sd);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(y_true, y_sd);
}
"


# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.97
TREE_DEPTH <- 15
STEP_SIZE <- 0.001
CHAINS <- 4
ITER <- 6000

tt <- proc.time()
fit <- rstan::stan(model_code = model,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 'beta', 'mu_beta', 
                            'sigma_beta', 'alpha_beta', 'beta_beta',
                            'sigma', 'sigma_sigma', 'alpha_sigma',
                            'beta_sigma', 'y_true', 'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time() - tt[3]) / 60





# Calc diagnostics ---------------------------------------------------

# library(shinystan)
# launch_shinystan(fit)

num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
mn_stepsize <- round(sapply(sampler_params, 
                            function(x) mean(x[, 'stepsize__'])), 5)
mn_treedepth <- round(sapply(sampler_params, 
                             function(x) mean(x[, 'treedepth__'])), 1)
accept_stat <- round(sapply(sampler_params, 
                            function(x) mean(x[, 'accept_stat__'])), 2)






# old ---------------------------------------------------------------------


y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
bayesplot::ppc_dens_overlay(DATA$y_obs, y_rep[1:50,])

MCMCvis::MCMCsummary(fit, excl = c('y_true', 'y_rep'), round = 2)
MCMCvis::MCMCplot(fit, params = 'beta', rank = TRUE)
#MCMCvis::MCMCplot(fit, params = 'sigma', rank = TRUE)




# plot results ------------------------------------------------------------

alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')

x_sim <- seq(0, 18, length = 50)

FIT_PLOT <- data.frame(cell = rep(1:DATA$NC, each = length(x_sim)), 
                       x_sim = rep(x_sim, times = DATA$NC),
                       mu_rep_mn = NA, 
                       mu_rep_LCI = NA, 
                       mu_rep_UCI = NA)

counter <- 1
for (i in 1:NCOL(alpha_ch))
{
  #i <- 1
  for (j in 1:length(x_sim))
  {
    #j <- 1
    temp <- alpha_ch[,i] + beta_ch[,i] * x_sim[j]
    tmn <- mean(temp)
    tsd <- sd(temp)
    tLCI <- tmn - tsd
    tUCI <- tmn + tsd
    FIT_PLOT$mu_rep_mn[counter] <- tmn
    FIT_PLOT$mu_rep_LCI[counter] <- tLCI
    FIT_PLOT$mu_rep_UCI[counter] <- tUCI
    counter <- counter + 1
  }
}


#mean and +- 1 sd
y_true_mn <- MCMCvis::MCMCpstr(fit, params = 'y_true', func = mean)[[1]]
y_true_sd <- MCMCvis::MCMCpstr(fit, params = 'y_true', func = sd)[[1]]
y_true_LCI <- y_true_mn - y_true_sd
y_true_UCI <- y_true_mn + y_true_sd

DATA_PLOT <- data.frame(y_true_mn, y_true_LCI, y_true_UCI,
                        cell = DATA$cn_id,
                        year = DATA$year,
                        y_obs = DATA$y_obs,
                        lat = data_f2$cell_lat)



#pdf('Figure_6.pdf', height = 11, width = 9, useDingbats = FALSE)
ggplot(data = DATA_PLOT, aes(DATA$year, y_true_mn), color = 'black', alpha = 0.6) +
  # geom_ribbon(data = FIT_PLOT,
  #             aes(x = x_sim, ymin = mu_rep_LCI, ymax = mu_rep_UCI, 
  #                 group = cell),#, fill = cell),
  #             #fill = 'grey',
  #             alpha = 0.6) +
  # geom_line(data = FIT_PLOT, aes(x_sim, mu_rep_mn, group = cell, col = factor(cell)),
  #           alpha = 0.9,
  #           inherit.aes = FALSE,
  #           size = 0.8) +
  geom_point(data = DATA_PLOT,
             aes(year, y_true_mn), color = 'black',
             inherit.aes = FALSE, size = 2, alpha = 0.7) +
  geom_point(data = DATA_PLOT,
             aes(year, y_true_mn, color = factor(cell)),
             inherit.aes = FALSE, size = 1.5, alpha = 1) +
  # geom_point(data = DATA_PLOT,
  #            aes(year, y_obs), color = 'black',
  #            inherit.aes = FALSE, size = 1.5, alpha = 0.5) +
  # geom_errorbar(data = DATA_PLOT,
  #               aes(ymin = y_true_LCI, ymax = y_true_UCI), width = 0.4,
  #               color = 'black', alpha = 0.4) +
  theme_bw() +
  theme(legend.position='none') +
  xlab('Year') +
  ylab('Arrival date') +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
#dev.off()



# plot anomalies ----------------------------------------------------------

out <- data.frame()
for (i in 1:DATA$NC)
{
  #i <- 1
  temp <- dplyr::filter(DATA_PLOT, cell == i)
  tsc <- scale(temp$y_true_mn, scale = FALSE)
  tsc_obs <- scale(temp$y_obs, scale = FALSE)
  temp$y_true_sc <- tsc
  temp$y_obs_sc <- tsc_obs
  temp$year <- temp$year + 2001
  out <- rbind(out, temp)
}



pdf('Vireo_olivaceus_anomaly.pdf', height = 5, width = 7, useDingbats = FALSE)
ggplot(data = out, aes(year, y_true_sc), color = 'black', alpha = 0.6) +
  # geom_ribbon(data = FIT_PLOT,
  #             aes(x = x_sim, ymin = mu_rep_LCI, ymax = mu_rep_UCI),
  #             fill = 'grey',
  #             inherit.aes = FALSE,
  #             alpha = 0.6) +
  # geom_line(data = FIT_PLOT, aes(x_sim, mu_rep_mn),
  #           col = 'red',
  #           alpha = 0.9,
  #           inherit.aes = FALSE,
  #           size = 0.8) +
  # geom_point(data = out,
  #            aes(year, y_true_sc), color = 'black',
  #            inherit.aes = FALSE, size = 2, alpha = 0.7) +
  geom_point(data = out,
             aes(year, y_true_sc, color = lat), #color = factor(cell)),
             inherit.aes = FALSE, size = 3, alpha = 0.7) +
  scale_color_gradient(low = 'red', high = 'white') +
  # geom_point(data = DATA_PLOT,
  #            aes(year, y_obs), color = 'black',
  #            inherit.aes = FALSE, size = 1.5, alpha = 0.5) +
  # geom_errorbar(data = DATA_PLOT,
  #               aes(ymin = y_true_LCI, ymax = y_true_UCI), width = 0.4,
  #               color = 'black', alpha = 0.4) +
  theme_bw() +
  #theme(legend.position='none') +
  xlab('Year') +
  ylab('Arrival anomaly') +
  ylim(c(-5, 5)) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm'),) #length of axis tick
dev.off()


