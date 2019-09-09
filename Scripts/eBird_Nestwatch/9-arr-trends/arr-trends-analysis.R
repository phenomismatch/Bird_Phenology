#############################
# Analysis of pheno trends output
#
#############################

# PARAMETER DESCRIPTIONS
# ======================

# alpha_gamma - arrival date of that species at 0 degrees lat
# beta_gamma - speed of migration (days / degrees lat)
# mu_alpha - absolute arrival date for species (intercepts drawn from this mean)
# sigma - interannual variability in arrival date after accounting for trend
# alpha_beta - magnitude of phenological change over time
# beta_beta - effect of lat on magnitude of phenological change over time


# set dirs ----------------------------------------------------------------

trends_date <- '2019-09-04'
IAR_date <- '2019-05-26'

trends_summary_dir <- paste0('~/Desktop/Bird_Phenology_Offline/Data/Processed/trends_summary_', trends_date)

arrival_master_dir <- paste0('~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_master_', IAR_date)


# load packages -----------------------------------------------------------

library(dplyr)
library(rstan)
library(MCMCvis)


# read in data ------------------------------------------------------------

setwd(trends_summary_dir)

pt_out <- readRDS(paste0('pheno_trends_master_', trends_date, '.rds'))
sp <- pt_out$species

pt2 <- dplyr::filter(pt_out, per_cd > 0.2, rng_lat > 5)

sigma_post <- readRDS(paste0('pheno_trends_sigma_post_', trends_date, '.rds'))
alpha_beta_post <- readRDS(paste0('pheno_trends_alpha_beta_post_', trends_date, '.rds'))
beta_beta_post <- readRDS(paste0('pheno_trends_beta_beta_post_', trends_date, '.rds'))

setwd(arrival_master_dir)

alpha_gamma_post <- readRDS('arrival_alpha_gamma_post_2019-05-26.rds')
beta_gamma_post <- readRDS('arrival_beta_gamma_post_2019-05-26.rds')



# plot beta_beta ----------------------------------------------------------

#simulate beta_beta (gamma) posteriors
temp <- matrix(NA, nrow = 1000, ncol = NROW(pt2))
for (i in 1:NROW(pt2))
{
  #i <- 1
  temp[,i] <- rnorm(1000, pt2$mn_beta_beta[i], pt2$sd_beta_beta[i])
}
colnames(temp) <- unique(pt2$species)
labs <- gsub('_', ' ', unique(pt2$species))


setwd('~/Desktop')
pdf('gamma.pdf', height = 14, width = 8.5, useDingbats = FALSE)
MCMCvis::MCMCplot(object = temp, 
                  labels = labs,
                  main = 'gamma - lat effect on slope',
                  guide_lines = TRUE, sz_labels = 0.8,
                  rank = TRUE)
dev.off()




# size doesn't matter for halfmax -----------------------------------------

# tt <- rnorm(100, out$mn_sigma[1], out$sd_sigma[1])
# d <- density(tt)
# plot(d)
# xmax <- d$x[d$y == max(d$y)]
# 
# x1 <- d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax] - max(d$y)/2))]
# x2 <- d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax] - max(d$y)/2))]
# points(c(x1, x2), c(d$y[d$x==x1], d$y[d$x==x2]), col="red")




# alpha_beta (and beta_beta) ~ sigma - lm in loop --------------------------------------

ab_slope <- c()
bb_slope <- c()
for (i in 1:NROW(sigma_post))
{
  #i <- 1

  tv_sigma <- sigma_post[i,]
  tv_alpha_beta <- alpha_beta_post[i,]
  tv_beta_beta <- beta_beta_post[i,]

  tfit <- lm(abs(tv_alpha_beta) ~ tv_sigma)
  tfit2 <- lm(abs(tv_beta_beta) ~ tv_sigma)

  #slope
  tslope <- coefficients(tfit)[2]
  tslope2 <- coefficients(tfit2)[2]

  #output
  ab_slope <- c(ab_slope, tslope)
  bb_slope <- c(bb_slope, tslope2)
}


setwd(trends_summary_dir)
pdf('alpha-beta_sigma_hist.pdf')
hist(ab_slope, xlim = c(0, 0.5), xlab = 'Slope estimate', main = '|alpha_beta| ~ sigma')
dev.off()

pdf('beta-beta_sigma_hist.pdf')
hist(bb_slope, xlim = c(0, 0.15), xlab = 'Slope estimate', main = '|beta_beta| ~ sigma')
dev.off()




# alpha_gamma ~ sigma ------------------------------------------------------

rng_sigma <- range(pt2$mn_sigma)

#data for stan model
DATA <- list(N = NROW(pt_out),
             y_obs = pt_out$mn_alpha_gamma,
             sd_y = pt_out$sd_alpha_gamma,
             x_obs = pt_out$mn_sigma,
             sd_x = pt_out$sd_sigma,
             pred_sim =  seq(from = rng_sigma[1], to = rng_sigma[2], length.out = 100),
             N_pred_sim = 100)

#lower limit is 0 for predictor
stanmodel <- "
data {
int<lower=0> N;
vector[N] y_obs;
vector<lower=0>[N] sd_y;
vector<lower=0>[N] x_obs;
vector<lower=0>[N] sd_x;
int<lower=0> N_pred_sim;
vector[N_pred_sim] pred_sim;
}

parameters {
vector[N] y_true_raw;
vector<lower=0>[N] x_true_raw;
real<lower=0> sigma_raw;
real alpha_raw;
real beta_raw;
}

transformed parameters {
vector[N] mu;
vector[N_pred_sim] mu_rep;
real<lower=0> sigma;
real alpha;
real beta;
vector[N] y_true;
vector<lower=0>[N] x_true;

sigma = sigma_raw * 10;
alpha = alpha_raw * 10;
beta = beta_raw * 10;
x_true = x_true_raw * 10;

mu = alpha + beta * x_true;
mu_rep = alpha + beta * pred_sim;

y_true = y_true_raw * sigma + mu;
}

model {

alpha_raw ~ std_normal();
beta_raw ~ std_normal();
sigma_raw ~ std_normal();
x_true_raw ~ std_normal();

// account for observation error in y
y_obs ~ normal(y_true, sd_y);

// account for observation error in x
x_obs ~ normal(x_true, sd_x);

y_true_raw ~ std_normal();      // implies y_true ~ normal(mu, sigma);
}

generated quantities {
real y_rep[N];
real x_rep[N];

y_rep = normal_rng(y_true, sd_y);
x_rep = normal_rng(x_true, sd_x);
}
"

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000


tt <- proc.time()
fit1 <- rstan::stan(model_code = stanmodel,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'sigma',
                             'y_true',
                             'x_true',
                             'y_rep',
                             'x_rep',
                             'mu_rep'),
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(trends_summary_dir)
saveRDS(fit1, file = 'alpha-gamma_sigma.rds')




# beta_gamma ~ sigma ------------------------------------------------------

DATA <- list(N = NROW(pt_out),
             y_obs = pt_out$mn_beta_gamma,
             sd_y = pt_out$sd_beta_gamma,
             x_obs = pt_out$mn_sigma,
             sd_x = pt_out$sd_sigma,
             pred_sim =  seq(from = rng_sigma[1], to = rng_sigma[2], length.out = 100),
             N_pred_sim = 100)

DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit2 <- rstan::stan(model_code = stanmodel,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'sigma',
                             'y_true',
                             'x_true',
                             'y_rep',
                             'x_rep',
                             'mu_rep'),
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(trends_summary_dir)
saveRDS(fit2, file = 'beta-gamma_sigma.rds')



# alpha_beta ~ sigma ------------------------------------------------------

# mn_alpha_beta <- apply(alpha_beta_post, 2, function(x) mean(abs(x)))
# sd_alpha_beta <- apply(alpha_beta_post, 2, function(x) sd(abs(x)))

DATA <- list(N = NROW(pt2),
             # y_obs = mn_alpha_beta,
             # sd_y = sd_alpha_beta,
             y_obs = pt2$mn_alpha_beta,
             sd_y = pt2$sd_alpha_beta,
             x_obs = pt2$mn_sigma,
             sd_x = pt2$sd_sigma,
             pred_sim =  seq(from = rng_sigma[1], to = rng_sigma[2], length.out = 100),
             N_pred_sim = 100)

DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit3 <- rstan::stan(model_code = stanmodel,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'sigma',
                             'y_true',
                             'x_true',
                             'y_rep',
                             'x_rep',
                             'mu_rep'),
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(trends_summary_dir)
saveRDS(fit3, file = 'alpha-beta_sigma.rds')



# beta_beta ~ sigma -------------------------------------------------------

# mn_beta_beta <- apply(beta_beta_post, 2, function(x) mean(abs(x)))
# sd_beta_beta <- apply(beta_beta_post, 2, function(x) sd(abs(x)))

DATA <- list(N = NROW(pt2),
             # y_obs = mn_beta_beta,
             # sd_y = sd_beta_beta,
             y_obs = pt2$mn_beta_beta,
             sd_y = pt2$sd_beta_beta,
             x_obs = pt2$mn_sigma,
             sd_x = pt2$sd_sigma,
             pred_sim =  seq(from = rng_sigma[1], to = rng_sigma[2], length.out = 100),
             N_pred_sim = 100)

DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit4 <- rstan::stan(model_code = stanmodel,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'sigma',
                             'y_true',
                             'x_true',
                             'y_rep',
                             'x_rep',
                             'mu_rep'),
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(trends_summary_dir)
saveRDS(fit4, file = 'beta-beta_sigma.rds')


# alpha_gamma ~ beta_gamma -------------------------------------------------------

#rng_beta_gamma <- range(pt_out$mn_beta_gamma)
rng_beta_gamma <- c(1, 4.5)

DATA <- list(N = NROW(pt_out),
             y_obs = pt_out$mn_alpha_gamma,
             sd_y = pt_out$sd_alpha_gamma,
             x_obs = pt_out$mn_beta_gamma,
             sd_x = pt_out$sd_beta_gamma,
             pred_sim =  seq(from = rng_beta_gamma[1], to = rng_beta_gamma[2], 
                             length.out = 100),
             N_pred_sim = 100)


#lower limit is 0 for predictor
stanmodel_agbg <- "
data {
int<lower=0> N;
vector[N] y_obs;
vector<lower=0>[N] sd_y;
vector<lower=0>[N] x_obs;
vector<lower=0>[N] sd_x;
int<lower=0> N_pred_sim;
vector[N_pred_sim] pred_sim;
}

parameters {
vector[N] y_true_raw;
vector<lower=0>[N] x_true_raw;
real<lower=0> sigma_raw;
real alpha_raw;
real beta_raw;
}

transformed parameters {
vector[N] mu;
vector[N_pred_sim] mu_rep;
real<lower=0> sigma;
real alpha;
real beta;
vector[N] y_true;
vector<lower=0>[N] x_true;

sigma = sigma_raw * 50;
alpha = alpha_raw * 75;
beta = beta_raw * 50;
x_true = x_true_raw * 10;

mu = alpha + beta * x_true;
mu_rep = alpha + beta * pred_sim;

y_true = y_true_raw * sigma + mu;
}

model {

alpha_raw ~ std_normal();
beta_raw ~ std_normal();
sigma_raw ~ std_normal();
x_true_raw ~ std_normal();

// account for observation error in y
y_obs ~ normal(y_true, sd_y);

// account for observation error in x
x_obs ~ normal(x_true, sd_x);

y_true_raw ~ std_normal();      // implies y_true ~ normal(mu, sigma);
}

generated quantities {
real y_rep[N];
real x_rep[N];

y_rep = normal_rng(y_true, sd_y);
x_rep = normal_rng(x_true, sd_x);
}
"


DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.005
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit5 <- rstan::stan(model_code = stanmodel_agbg,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('alpha',
                             'beta',
                             'sigma',
                             'y_true',
                             'x_true',
                             'y_rep',
                             'x_rep',
                             'mu_rep'),
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(trends_summary_dir)
saveRDS(fit5, file = 'alpha-gamma_beta-gamma.rds')



# plot plasticity ---------------------------------------------------------


plt_fun <- function(INPUT, RNG, MAIN, XLAB, YLAB, TITLE)
{
  #model fit
  mu_rep_mn <- MCMCvis::MCMCpstr(INPUT, params = 'mu_rep', 
                                 func = mean)[[1]]
  mu_rep_LCI <- MCMCvis::MCMCpstr(INPUT, params = 'mu_rep', 
                                  func = function(x) quantile(x, probs = 0.025))[[1]]
  mu_rep_UCI <- MCMCvis::MCMCpstr(INPUT, params = 'mu_rep', 
                                  func = function(x) quantile(x, probs = 0.975))[[1]]
  
  #data to plot model fit
  FIT_PLOT <- data.frame(x_sim = seq(from = RNG[1], 
                                     to = RNG[2], length.out = 100), 
                         mu_rep_mn, mu_rep_LCI, mu_rep_UCI)
  
  #mean and +- 1 sd for data
  y_true_mn <- MCMCvis::MCMCpstr(INPUT, params = 'y_true', func = mean)[[1]]
  y_true_sd <- MCMCvis::MCMCpstr(INPUT, params = 'y_true', func = sd)[[1]]
  y_true_LCI <- y_true_mn - y_true_sd
  y_true_UCI <- y_true_mn + y_true_sd
  
  x_true_mn <- MCMCvis::MCMCpstr(INPUT, params = 'x_true', func = mean)[[1]]
  x_true_sd <- MCMCvis::MCMCpstr(INPUT, params = 'x_true', func = sd)[[1]]
  x_true_LCI <- x_true_mn - x_true_sd
  x_true_UCI <- x_true_mn + x_true_sd
  
  
  DATA_PLOT <- data.frame(y_true_mn, y_true_LCI, y_true_UCI,
                          x_true_mn, x_true_LCI, x_true_UCI)
  
  plt <- ggplot(data = DATA_PLOT, aes(x_true_mn, y_true_mn), color = 'black', alpha = 0.6) +
    geom_ribbon(data = FIT_PLOT, 
                aes(x = x_sim, ymin = mu_rep_LCI, ymax = mu_rep_UCI),
                fill = 'grey', alpha = 0.6,
                inherit.aes = FALSE) +
    geom_line(data = FIT_PLOT, aes(x_sim, mu_rep_mn), color = 'red',
              alpha = 0.9,
              inherit.aes = FALSE,
              size = 1.4) +
    geom_point(data = DATA_PLOT,
               aes(x_true_mn, y_true_mn), color = 'black',
               inherit.aes = FALSE, size = 3, alpha = 0.4) +
    geom_errorbar(data = DATA_PLOT,
                  aes(ymin = y_true_LCI, ymax = y_true_UCI), #width = 0.05,
                  color = 'black', alpha = 0.2) +
    geom_errorbarh(data = DATA_PLOT,
                   aes(xmin = x_true_LCI, xmax = x_true_UCI), #height = 0.05,
                   color = 'black', alpha = 0.2) +
    theme_bw() +
    ggtitle(MAIN) +
    xlab(XLAB) +
    ylab(YLAB) +
    theme(
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
      axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
  
  ggsave(paste0(TITLE, '.pdf'), plt)
}


setwd(trends_summary_dir)
plt_fun(fit1, rng_sigma, MAIN = 'alpha_gamma ~ sigma', 
        XLAB = 'Interannual variability', 
        YLAB = 'Relative arrival date',
        TITLE = 'alpha-gamma_sigma')
plt_fun(fit2, rng_sigma, MAIN = 'beta_gamma ~ sigma', 
        XLAB = 'Interannual variability', 
        YLAB = 'Migration speed (days/degree lat)',
        TITLE = 'beta-gamma_sigma')
plt_fun(fit3, rng_sigma, MAIN = 'alpha_beta ~ sigma', 
        XLAB = 'Interannual variability', 
        YLAB = 'Phenological change',
        TITLE = 'alpha-beta_sigma')
plt_fun(fit4, rng_sigma, MAIN = 'beta_beta ~ sigma', 
        XLAB = 'Interannual variability', 
        YLAB = 'Effect of lat on rate of phenological change',
        TITLE = 'beta-beta_sigma')
plt_fun(fit5, rng_beta_gamma, MAIN = 'alpha_gamma ~ beta_gamma', 
        XLAB = 'Migration speed (days/degree lat)', 
        YLAB = 'Relative arrival date',
        TITLE = 'alpha-gamma_beta_gamma')



MCMCsummary(fit1, round = 2, excl = c('y_true', 'x_true', 'y_rep', 'x_rep', 'mu_rep'))
MCMCsummary(fit2, round = 2, excl = c('y_true', 'x_true', 'y_rep', 'x_rep', 'mu_rep'))
MCMCsummary(fit3, round = 2, excl = c('y_true', 'x_true', 'y_rep', 'x_rep', 'mu_rep'))
MCMCsummary(fit4, round = 2, excl = c('y_true', 'x_true', 'y_rep', 'x_rep', 'mu_rep'))
MCMCsummary(fit5, round = 2, excl = c('y_true', 'x_true', 'y_rep', 'x_rep', 'mu_rep'))










# plot parameter estimates ------------------------------------------------


MCMCvis::MCMCplot(sigma_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(0, 10), main = 'sigma')
MCMCvis::MCMCplot(alpha_beta_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(-1, 1.5), main = 'alpha_beta')
MCMCvis::MCMCplot(beta_beta_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(-0.25, 0.5), main = 'beta_beta')
MCMCvis::MCMCplot(alpha_gamma_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, #xlim = c(-1, 1.5), 
                  main = 'alpha_gamma')
MCMCvis::MCMCplot(beta_gamma_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, #xlim = c(-0.25, 0.5), 
                  main = 'beta_gamma')


# traits ------------------------------------------------------------------

setwd('~/Google_Drive/R/Bird_Phenology/Data/Traits')

traits <- read.csv('Trait_database-2019-06-17.csv', stringsAsFactors = FALSE)
#add underscore
traits$species <- gsub(' ', '_', traits$SCI_NAME)

out_j <- dplyr::left_join(pt_out, traits, by = 'species')

# out_jf <- dplyr::select(out_j, species, mn_mu_alpha, sd_mu_alpha, mn_sigma, 
#                         sd_sigma, mn_alpha_beta, sd_alpha_beta, mn_beta_beta,
#                         sd_beta_beta, IUCN_STATUS, BODY_MASS_ELTON, MIGRATION_DISTANCE_LASORTE, 
#                         SPRING_MIGRATION_SPEED_LASORTE, CLUTCH_SIZE_BONA, BROODS_PER_YEAR_BONA, 
#                         EGG_MASS_BONA, MAX_LIFESPAN_BONA, INCUBATION_PERIOD_LENGTH_BONA, 
#                         FLEDGING_AGE_BONA, DIET_INV_ELTON, DIET_VUNK_ELTON, DIET_SCAV_ELTON, 
#                         DIET_FRUIT_ELTON, DIET_SEED_ELTON, rng_lat, n_cells, n_years, num_diverge, max_rhat, min_neff)
# 
# 
# 
# idx <- which(colnames(out_jf) %in% c('species', 'IUCN_STATUS', 'sd_mu_alpha', 'sd_sigma', 
#                                      'sd_alpha_beta', 'sd_beta_beta', 
#                                      'num_diverge', 'max_rhat', 'min_neff'))
# 
# pairs(out_jf[,-idx])

summary(lm(out_j$mn_beta_gamma ~ out_j$SPRING_MIGRATION_SPEED_LASORTE))
summary(lm(out_j$SPRING_MIGRATION_SPEED_LASORTE ~ out_j$BODY_MASS_LASORTE))
summary(lm(out_j$mn_beta_gamma ~ out_j$BODY_MASS_LASORTE))


out_jf2 <- dplyr::filter(out_jf, n_cells > 10, n_years > 10, rng_lat > 10)

# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_mu_alpha)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_mu_alpha)
# 
# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_sigma)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_sigma)
# 
# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_alpha_beta)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_alpha_beta)
# 
# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_beta_beta)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_beta_beta)

fit_fun(out_jf$SPRING_MIGRATION_SPEED_LASORTE, out_jf$MIGRATION_DISTANCE_LASORTE)
fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$CLUTCH_SIZE_BONA)


