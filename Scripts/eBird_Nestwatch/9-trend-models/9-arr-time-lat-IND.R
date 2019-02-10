####################
# 9 - arrival data ~ year + lat
# 
# Individual species
#
####################


# top-level dir --------------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'

MODEL_DATE <- '2019-02-08'
ARR_TIME_LAT_IND_DIR <- paste0('trend_models_', MODEL_DATE)


# species arg -----------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Vireo_olivaceus')


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)
library(dggridR)
library(gridExtra)


# import ARR/BR data ---------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

#arrival data
IAR_out <- 'IAR_output_2019-01-16'
IAR_out_date <- substr(IAR_out, start = 12, stop = 21)

IAR_data <- readRDS(paste0('arrival_master_', IAR_out_date, '.rds'))


mdf <- dplyr::select(IAR_data, species, cell, year, cell_lat, cell_lon, 
                     mean_post_IAR, sd_post_IAR)



#OVERVIEW
#=======#
#arrival ~ N(alpha_j + beta_j * year)
#beta_j = alpha2 + beta2 * lat


#DETAILS
#======#
# #observation model - where obs and sigma are known
# x_obs[i] ~ normal(x_true[i], sigma_x[i])
# 
# #arrival as a function of year and latitude
# x_true[i] ~ normal(mu_arr[i], sigma_arr)
# mu_arr[i] = alpha[cn[i]] + beta[cn[i]] * year[i]
# for (j in J)
#   beta1[j] = alpha2 + beta2 * lat[j]


#do not merge with breeding data
mdf3 <- filter(mdf, species == args)

#number codes for species
#sp_num <- as.numeric(factor(mdf3$species))

#number code for cells (not actual cell number)
cell_num <- as.numeric(factor(mdf3$cell))

#lats of cells
cell_mrg <- data.frame(cell = mdf3$cell, cell_num = cell_num, 
                       cell_lat = mdf3$cell_lat)

u_cell_mrg <- unique(cell_mrg)

#number codes for years
yr_num <- as.numeric(factor(mdf3$year))


#create data list for Stan
DATA <- list(x_obs = mdf3$mean_post_IAR,
             sigma_x = mdf3$sd_post_IAR,
             cn_id = cell_num,
             year = yr_num,
             lat = u_cell_mrg$cell_lat,
             US = NROW(u_cell_mrg),
             N = NROW(mdf3))



# Stan model --------------------------------------------------------------

arr_time_lat_ind <- '
data {
int<lower = 0> N;                                     // number of obs
int<lower = 0> US;                                    // number of species
vector<lower = 0, upper = 200>[N] x_obs;                // mean halfmax IAR
vector<lower = 0>[N] sigma_x;                           // sd halfmax IAR
int<lower = 1, upper = US> cn_id[N];                  // species ids
vector<lower = 1, upper = 17>[N] year;
vector<lower = 26, upper = 90>[US] lat;
}

parameters {
vector<lower = 0, upper = 200>[N] x_true;                           //true arrival
real mu_alpha_raw;
real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_raw;
real alpha_raw[US];
real alpha2_raw;
real beta2_raw;
}

transformed parameters {
real mu_alpha;
real sigma_alpha;
real sigma;
real alpha[US];
real beta[US];
real alpha2;
real beta2;
real mu[N];

// non-centered parameterization
mu_alpha = mu_alpha_raw * 20 + 70;                       // implies mu_alpha ~ normal(70, 20)
sigma_alpha = sigma_alpha_raw * 10;                      // implies sigma_alpha ~ halfnormal(0, 10)
sigma = sigma_raw * 10;                                  // implies sigma ~ halfnormal(0, 10)
alpha2 = alpha2_raw * 10;                                 // implies alpha2 ~ normal(0, 10)
beta2 = beta2_raw * 2 + 1;                               // implies beta2 ~ normal(1, 2)

for (j in 1:US)
{
  alpha[j] = alpha_raw[j] * sigma_alpha + mu_alpha;      // implies alpha[j] ~ normal(mu_alpha, sigma_alpha)
  beta[j] = alpha2 + beta2 * lat[j];
}

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
}
}

model {

// observation model - modeling true state as a function of some observed state
x_obs ~ normal(x_true, sigma_x);

// non-centered parameterization
mu_alpha_raw ~ normal(0, 1);
sigma_alpha_raw ~ normal(0, 1);
sigma_raw ~ normal(0, 1);
alpha2_raw ~ normal(0, 1);
beta2_raw ~ normal(0, 1);

for (j in 1:US)
{
  alpha_raw[j] ~ normal(0, 1);
}

x_true ~ normal(mu, sigma);
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


DELTA <- 0.85
TREE_DEPTH <- 15
STEP_SIZE <- 0.05
CHAINS <- 3
ITER <- 3000

tt <- proc.time()
fit <- rstan::stan(model_code = arr_time_lat_ind,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 'beta', 'mu_alpha', 
                            'alpha2', 'beta2',
                            'sigma_alpha', 'sigma', 'x_true'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time() - tt[3]) / 60



#save to RDS
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', ARR_TIME_LAT_IND_DIR))
saveRDS(fit, file = paste0('temp_ARR_YEAR_LAT_IND_stan_', MODEL_DATE, '_', args, '.rds'))
#fit <- readRDS(paste0('temp_ARR_YEAR_LAT_IND_stan_', MODEL_DATE, '_', args, '.rds'))




# diagnostics -------------------------------------------------------------


# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'alpha', ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'beta', ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'sigma', ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'mu', ISB = FALSE)
# MCMCvis::MCMCsummary(fit, n.eff = TRUE, params = 'x_true')
# MCMCvis::MCMCplot(fit, params = 'beta1', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'beta2', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'beta3', rank = TRUE)
# MCMCvis::MCMCplot(fit, params = 'mu', ISB = FALSE, rank = TRUE)
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

options(max.print = 50000)
sink(paste0('ARR_YEAR_LAT_IND_results_', MODEL_DATE, '_', args, '.txt'))
cat(paste0('ARR_YEAR_LAT_IND_results_', MODEL_DATE, '_', args, ' \n'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
print(fit)
sink()


# Extract posterior estimates ---------------------------------------------

#extract median and sd estimates for x_true
x_true_mean <- MCMCvis::MCMCpstr(fit, params = 'x_true', func = mean)[[1]]
x_true_sd <- MCMCvis::MCMCpstr(fit, params = 'x_true', func = sd)[[1]]

#extract slope estimates and CI
beta_med <- MCMCvis::MCMCpstr(fit, params = 'beta', func = median)[[1]]
beta_LCI <- MCMCvis::MCMCpstr(fit, params = 'beta', func = function(x) quantile(x, probs = c(0.025)))[[1]]
beta_UCI <- MCMCvis::MCMCpstr(fit, params = 'beta', func = function(x) quantile(x, probs = c(0.975)))[[1]]
beta_sd <- MCMCvis::MCMCpstr(fit, params = 'beta', func = sd)[[1]]


# plots arr over time -----------------------------------------------------

#for each cell:
#X-AXIS = year
#Y-AXIS = julian day
#x_true with error bars
#regression line

#need true latent states
DATA_PLOT <- data.frame(mean_x = x_true_mean,
                        mean_x_l = x_true_mean - x_true_sd,
                        mean_x_u = x_true_mean + x_true_sd,
                        IAR_x = mdf3$mean_post_IAR,
                        IAR_x_l = mdf3$mean_post_IAR - mdf3$sd_post_IAR,
                        IAR_x_u = mdf3$mean_post_IAR + mdf3$sd_post_IAR,
                        year = mdf3$year,
                        cell = cell_mrg$cell,
                        beta_med = beta_med,
                        beta_LCI = beta_LCI,
                        beta_UCI = beta_UCI,
                        beta_sd = beta_sd)

#plot for each cell
for (i in 1:NROW(u_cell_mrg))
{
  #i <- 1
  print(i)
  
  #model fit to plot
  alpha_ch <- MCMCchains(fit, params = paste0('alpha\\[', 
                                              u_cell_mrg$cell_num[i], '\\]'), ISB = FALSE)[,1]
  beta_ch <- MCMCchains(fit, params = paste0('beta\\[', 
                                             u_cell_mrg$cell_num[i], '\\]'), ISB = FALSE)[,1]

  sim_year <- seq(1, length(unique(mdf3$year)), length = 100)
  sim_year_actual <- seq(min(mdf3$year), max(mdf3$year), length = 100)
  
  mf <- matrix(nrow = length(beta_ch), ncol = 100)
  for (j in 1:length(sim_year))
  {
    mf[,j] <- alpha_ch + beta_ch * sim_year[j]
  }
  
  med_mf <- apply(mf, 2, median)
  LCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.025))
  UCI_mf <- apply(mf, 2, function(x) quantile(x, probs = 0.975))
  
  FIT_PLOT <- data.frame(MN = med_mf, 
                         MN_YR = sim_year_actual,
                         LCI = LCI_mf,
                         UCI = UCI_mf)

  DATA_PLOT2 <- dplyr::filter(DATA_PLOT, cell == u_cell_mrg$cell[i])
  
  p <- ggplot(data = DATA_PLOT2, aes(year, mean_x)) +
    geom_ribbon(data = FIT_PLOT,
                aes(x = MN_YR, ymin = LCI, ymax = UCI),
                fill = 'grey', alpha = 0.7,
                inherit.aes = FALSE) +
    geom_line(data = FIT_PLOT, aes(MN_YR, MN), color = 'red',
              alpha = 0.9,
              inherit.aes = FALSE,
              size = 1.4) +
    geom_point(data = DATA_PLOT2, aes(year - 0.1, mean_x), color = 'black',
               inherit.aes = FALSE, size = 1.75, alpha = 0.5) +
    geom_errorbar(data = DATA_PLOT2,
                  aes(ymin = mean_x_l, ymax = mean_x_u, x = year - 0.1), width = 0.5,
                  color = 'black', alpha = 0.4) +
    geom_point(data = DATA_PLOT2, aes(year + 0.1, IAR_x), color = 'red',
               inherit.aes = FALSE, size = 1.75, alpha = 0.3) +
    geom_errorbar(data = DATA_PLOT2,
                  aes(ymin = IAR_x_l, ymax = IAR_x_u, x = year + 0.1), width = 0.5,
                  color = 'red', alpha = 0.3) +
    theme_bw() +
    #scale_x_discrete(limits = c(seq(18,30, by = 2))) +
    xlab('Year') +
    ylab('True ARR') +
    ylim(low = 80, high = 170) +
    ggtitle(paste0(args, '; Cell: ', u_cell_mrg$cell[i])) +
    theme(
      plot.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.title.x = element_text(margin = margin(t = 5, r = 15, b = 0, l = 0)),
      axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
  
  assign(paste0('p_', i), p)
}


#save figs to pdfs
setwd(paste0(dir, 'Bird_Phenology/Figures/arrival_trends'))
pdf(paste0(args, '-plots-ARR-time.pdf'), height = 6, width = 9, useDingbats = FALSE)

counter <- 1
for (i in 1:ceiling(NROW(u_cell_mrg)/4))
{
  if ((NROW(u_cell_mrg) - counter) >= 4)
  {
    gridExtra::grid.arrange(eval(as.name(paste0('p_', counter))), 
                            eval(as.name(paste0('p_', counter + 1))), 
                            eval(as.name(paste0('p_', counter + 2))), 
                            eval(as.name(paste0('p_', counter + 3)))) 
  } else {
    if ((NROW(u_cell_mrg) - counter) == 3)
    {
      gridExtra::grid.arrange(eval(as.name(paste0('p_', counter))), 
                              eval(as.name(paste0('p_', counter + 1))), 
                              eval(as.name(paste0('p_', counter + 2)))) 
    }
    if ((NROW(u_cell_mrg) - counter) == 2)
    {
      gridExtra::grid.arrange(eval(as.name(paste0('p_', counter))), 
                              eval(as.name(paste0('p_', counter + 1))))
    }
    if ((NROW(u_cell_mrg) - counter) == 1)
    {
      gridExtra::grid.arrange(eval(as.name(paste0('p_', counter))))
    }
  }
  counter <- counter + 4
}
dev.off()




# plot slope estimates on each cell on map --------------------------------------

#slope in grey

#make hexgrid
hexgrid6 <- dggridR::dgconstruct(res = 6)

#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, u_cell_mrg$cell)
cell_grid$cell <- as.numeric(cell_grid$cell)
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, u_cell_mrg$cell)
ll_df <- data.frame(cell = u_cell_mrg$cell, 
                    lon_deg = cell_centers$lon_deg, 
                    lat_deg = cell_centers$lat_deg)

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])


#determine min/max for plotting
f_rng <- range(beta_med)
MIN <- round(min(f_rng), 3)
MAX <- round(max(f_rng), 3)



#merge hex spatial data with HM data
to_plt <- dplyr::inner_join(DATA_PLOT, cell_grid, by = 'cell')
to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')

#plot
fp <- ggplot() +
  geom_path(data = usamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = canadamap, 
            aes(x = x, y = y), color = 'black') + 
  geom_path(data = mexicomap, 
            aes(x = x, y = y), color = 'black') + 
  # geom_polygon(data = nrng.df, 
  #           aes(x = long, y = lat, group=group), fill = 'green', alpha = 0.4) + 
  # geom_polygon(data = nrng_rm.df, 
  #              aes(x = long, y = lat, group=group), fill = 'orange', alpha = 0.4) + 
  coord_map("ortho", orientation = c(35, -80, 0), 
            xlim = c(-100, -55), ylim = c(25, 66)) + 
  geom_polygon(data = to_plt2, aes(x = long, y = lat, group = group, fill = beta_med), 
               alpha = 0.4) +
  geom_path(data = to_plt2, aes(x = long, y = lat, group = group), 
            alpha = 0.4, color = 'black') + 
  scale_fill_gradientn(colors = c('red', 'blue'),
                       limits = c(MIN, MAX)) +
  labs(fill = 'Slope estimate') +
  annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5, 
           label = round(to_plt2$beta_med, digits = 2), col = 'black', alpha = 0.2,
           size = 3) +
  annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5, 
           label = round(to_plt2$beta_sd, digits = 2), 
                          col = 'white', alpha = 0.3,
           size = 2.5) +
  ggtitle(paste0(args, '-slope-ARR-time')) +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')

setwd(paste0(dir, 'Bird_Phenology/Figures/arrival_trends'))
ggsave(plot = fp, filename = paste0(args, '-slope-map-ARR-time.pdf'))

