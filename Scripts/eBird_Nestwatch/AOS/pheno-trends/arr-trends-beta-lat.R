#######
# phenology ~ time using post IAR
#
# i = obs
# j = cell
#
# y_{obs_{i}} \sim N(\mu_{y_{i}}, \sigma_{y})
# \mu_{y_{i}} \sim N(\mu_{i}, \sigma)
# \mu_{i} = \alpha_{j} + \beta_{j} \times year_{i}
# \alpha_{j} \sim N(\mu_{\alpha}, \sigma_{\alpha})
# \beta_{j} \sim N(\mu_{\beta_{j}}, \sigma_{\beta})
# \mu_{\beta_{j}} = \alpha_{\beta} + \beta_{\beta} \times lat_{j}
#######


# top-level dir --------------------------------------------------------------

#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# species -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
#args <- as.character('Vireo_olivaceus')
#args <- as.character('Ammospiza_nelsoni')
#args <- as.character('Agelaius_phoeniceus')
#args <- as.character('Turdus_migratorius')

#args <- 'Empidonax_virescens'
#args <- 'Hirundo_rustica'
#args <- 'Icterus_spurius'

# other dir ---------------------------------------------------------------

IAR_in_date <- '2019-05-03'
IAR_out_date <- '2019-05-26'
run_date <- '2019-09-03'

IAR_in_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/IAR_input_', IAR_in_date)

IAR_out_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/IAR_output_', IAR_out_date)

#IAR_out_dir <- paste0('~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_', IAR_out_date)

trends_out_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/trends_output_', run_date)

#trends_out_dir <- paste0('~/Desktop/Bird_Phenology_Offline/Data/Processed/trends_output_', run_date)


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(rstan)
library(MCMCvis)
library(dggridR)



# Filter data by species/years ------------------------------------------------------

setwd(IAR_in_dir)

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

#switch to out dir
setwd(IAR_out_dir)


#if that species RDS object exists in dir
pro_data <- data.frame()
if (length(grep(paste0(args, '-', IAR_out_date, '-iar-stan_output.rds'), list.files())) > 0)
{
  f_in_p <- dplyr::filter(df_master, species == args & MODEL == TRUE)
  
  #read in IAR model output and input
  t_fit <- readRDS(paste0(args, '-', IAR_out_date, '-iar-stan_output.rds'))
  t_data <- readRDS(paste0(args, '-', IAR_out_date, '-iar-stan_input.rds'))
  
  #only cells and years that were modeled (to account for any lone cells that were dropped in 4-IAR-arr.R)
  f_in <- dplyr::filter(f_in_p, cell %in% t_data$cells)
  
  t_cells <- unique(f_in$cell)
  t_years <- unique(f_in$year)
  
  #make hexgrid
  hexgrid6 <- dggridR::dgconstruct(res = 6)
  
  #get hexgrid cell centers
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)
  
  #extract median and sd for IAR arrival dates
  mean_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = mean)[[1]]
  med_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = median)[[1]]
  sd_fit <- MCMCvis::MCMCpstr(t_fit, params = 'y_true', func = sd)[[1]]
  
  #loop through years
  for (j in 1:length(t_years))
  {
    #j <- 1
    print(paste0('species: ', args, ', ', 
                 'year: ', t_years[j]))
    
    t_f_in <- dplyr::filter(f_in, year == t_years[j])
    
    t_full <- data.frame(t_f_in[,c('species','cell')], 
                         cell_lat = round(cellcenters$lat_deg, digits = 2), 
                         cell_lon = round(cellcenters$lon_deg, digits = 2),
                         t_f_in[,c('year', 'HM_mean', 'HM_sd')],
                         mean_post_IAR = mean_fit[,j], sd_post_IAR = sd_fit[,j])
    
    colnames(t_full)[6:7] <- c('mean_pre_IAR', 'sd_pre_IAR')
    
    pro_data <- rbind(pro_data, t_full)
  } #end year loop
} else {
  stop(paste0('.rds file for ', args, ' not found in directory'))
}



# Process data ------------------------------------------------------------

#cell years with input data
#data_f <- pro_data[which(!is.na(pro_data$mean_pre_IAR)),]

#all cell years
data_f <- pro_data

#cells with at least 5 years of data
cnts <- plyr::count(data_f, 'cell')
u_cells <- cnts[which(cnts[,2] >= 5),1]

data_f2 <- dplyr::filter(data_f, cell %in% u_cells)

t_cl <- unique(data_f2[,c('cell', 'cell_lat')])
ot_cl <- t_cl[order(t_cl[,1]),]

#years to generate low, mid, high lat phenology
sim_year <- min(data_f2$year - 2001):max(data_f2$year - 2001)

#create data list for Stan
DATA <- list(N = NROW(data_f2),
             y_obs = data_f2$mean_post_IAR,
             y_sd = data_f2$sd_post_IAR,
             cn_id = as.numeric(factor(data_f2$cell)),
             NC = NROW(ot_cl$cell),
             year = (data_f2$year - 2001),
             lat = scale(ot_cl$cell_lat, scale = FALSE)[,1],
             lat_usc = ot_cl$cell_lat)


# ggplot(data_f2, aes(x = year, y = mean_post_IAR, col = factor(cell))) +
#   geom_point() +
#   stat_smooth(method = 'lm')

# Stan model --------------------------------------------------------------

model <- "
data {
int<lower = 0> N;                                 // number of obs
vector<lower = 0, upper = 200>[N] y_obs;
vector<lower = 0>[N] y_sd;
int<lower = 1> cn_id[N];                          // cell ids
int<lower = 0> NC;                                // number of cells
vector<lower = 1, upper = 17>[N] year;
vector[NC] lat;
}

parameters {
vector[N] y_true_raw;
vector[NC] alpha_raw;
vector[NC] beta_raw;
vector[NC] mu_beta_raw;
real<lower = 0> sigma_beta_raw;
real alpha_beta_raw;
real beta_beta_raw;
real mu_alpha_raw;
real<lower = 0> sigma_alpha_raw;
real<lower = 0> sigma_raw;
}

transformed parameters {
vector[N] mu;
vector[N] y_true;
vector[NC] alpha;
vector[NC] beta;
vector[NC] mu_beta;
real<lower = 0> sigma_beta;
real alpha_beta;
real beta_beta;
real mu_alpha;
real<lower = 0> sigma_alpha;
real<lower = 0> sigma;

sigma = sigma_raw * 5;
mu_alpha = mu_alpha_raw * 50 + 120;
sigma_alpha = sigma_alpha_raw * 20;
alpha = alpha_raw * sigma_alpha + mu_alpha;

alpha_beta = alpha_beta_raw * 1;
beta_beta = beta_beta_raw * 1;
mu_beta = alpha_beta + beta_beta * lat;
sigma_beta = sigma_beta_raw * 3;
beta = beta_raw * sigma_beta + mu_beta;

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
  //y_true[i] = y_true_raw[i] * sigma[cn_id[i]] + mu[i];
  y_true[i] = y_true_raw[i] * sigma + mu[i];
}
}

model {

y_true_raw ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();

sigma_beta_raw ~ std_normal();

mu_alpha_raw ~ std_normal();
sigma_alpha_raw ~ std_normal();

alpha_beta_raw ~ std_normal();
beta_beta_raw ~ std_normal();
sigma_raw ~ std_normal();
// mu_sigma_raw ~ std_normal();
// sigma_sigma_raw ~ std_normal();

y_obs ~ normal(y_true, y_sd);
}

generated quantities {

real y_rep[N];

// PPC
y_rep = normal_rng(y_true, y_sd);
}
"


# Run model ---------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.99
TREE_DEPTH <- 15
STEP_SIZE <- 0.001
CHAINS <- 4
ITER <- 8000

tt <- proc.time()
fit <- rstan::stan(model_code = model,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('alpha', 
                            'mu_alpha', 
                            'sigma_alpha', 
                            'beta',
                            'sigma_beta', 
                            'alpha_beta',
                            'beta_beta',
                            'sigma', 
                            'y_true', 
                            'y_rep'),
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))


#rerun model with higher target acceptance if divergences exist
if (num_diverge > 0)
{
  DELTA <- 0.9999
  
  tt <- proc.time()
  fit <- rstan::stan(model_code = model,
                     data = DATA,
                     chains = CHAINS,
                     iter = ITER,
                     cores = CHAINS,
                     pars = c('alpha', 'mu_alpha', 'sigma_alpha', 'beta',
                              'sigma_beta', 'alpha_beta', 'beta_beta',
                              'sigma', 'y_true', 'y_rep'),
                     control = list(adapt_delta = DELTA,
                                    max_treedepth = TREE_DEPTH,
                                    stepsize = STEP_SIZE))
  run_time <- (proc.time()[3] - tt[3]) / 60
  
  num_diverge <- rstan::get_num_divergent(fit)
  num_tree <- rstan::get_num_max_treedepth(fit)
  num_BFMI <- length(rstan::get_low_bfmi_chains(fit))
}


#create dir if doesn't exist
ifelse(!dir.exists(trends_out_dir),
       dir.create(trends_out_dir),
       FALSE)

#save to RDS
setwd(trends_out_dir)

saveRDS(fit, file = paste0(args, '-', run_date, '-pheno_trends_stan_output.rds'))
saveRDS(DATA, file = paste0(args, '-', run_date, '-pheno_trends_stan_input.rds'))


# Calc diagnostics ---------------------------------------------------

#fit <- readRDS('Ictinia_mississippiensis-2019-05-26-pheno_trends_stan_output.rds')
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

y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')

# bayesplot::ppc_stat(DATA$y_obs, y_rep, stat = 'mean')
# bayesplot::ppc_dens_overlay(DATA$y_obs, y_rep[1:500,])


# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0(args, '-', run_date, '-pheno_trends_results.txt'))
cat(paste0('Pheno trends results ', args, ' \n'))
cat(paste0('Total minutes: ', round(run_time, digits = 2), ' \n'))
cat(paste0('Iterations: ', ITER, ' \n'))
cat(paste0('Adapt delta: ', DELTA, ' \n'))
cat(paste0('Max tree depth: ', TREE_DEPTH, ' \n'))
cat(paste0('Step size: ', STEP_SIZE, ' \n'))
cat(paste0('Number of divergences: ', num_diverge, ' \n'))
cat(paste0('Number of tree exceeds: ', num_tree, ' \n'))
cat(paste0('Number chains low BFMI: ', num_BFMI, ' \n'))
cat(paste0('Mean stepsize: ', round(mean(mn_stepsize), 5), ' \n'))
cat(paste0('Mean treedepth: ', round(mean(mn_treedepth), 1), ' \n'))
cat(paste0('Mean accept stat: ', round(mean(accept_stat), 2), ' \n'))
cat(paste0('Max Rhat: ', max(rhat_output), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output), ' \n'))
print(model_summary)
sink()



# density overlay plot ----------------------------------------------------

#modified bayesplot::ppc_dens_overlay function
tdata <- bayesplot::ppc_data(DATA$y_obs, y_rep[1:100,])

p <- ggplot(tdata) + 
  aes_(x = ~value) + 
  stat_density(aes_(group = ~rep_id, color = "yrep"), 
               data = function(x) dplyr::filter(x, !tdata$is_y), 
               geom = "line", position = "identity", size = 0.25, 
               alpha = 0.3, trim = FALSE, bw = 'nrd0', adjust = 1, 
               kernel = 'gaussian', n = 1024) + 
  stat_density(aes_(color = "y"), data = function(x) dplyr::filter(x, tdata$is_y), 
               geom = "line", position = "identity", lineend = "round", size = 1, trim = FALSE, 
               bw = 'nrd0', adjust = 1, kernel = 'gaussian', n = 1024) + 
  theme_classic() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 24)) +
  labs(colour = '') +
  scale_color_manual(values = c('red', 'black')) +
  ggtitle(paste0(args))

ggsave(paste0(args, '-', run_date, '-dens_overlay.pdf'), p)




# Map of trends ------------------------------------------

#estimated slope in grey, sd in white

#extract median and sd estimates for mu params
med_fit <- MCMCvis::MCMCpstr(fit, params = 'beta', func = median)[[1]]
sd_fit <- MCMCvis::MCMCpstr(fit, params = 'beta', func = sd)[[1]]

#transform cells to grid
cell_grid <- dggridR::dgcellstogrid(hexgrid6, ot_cl$cell)
cell_grid$cell <- as.numeric(cell_grid$cell)
cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, ot_cl$cell)
ll_df <- data.frame(cell = ot_cl$cell,
                    lon_deg = cell_centers$lon_deg,
                    lat_deg = cell_centers$lat_deg)

#load maps
usamap <- data.frame(maps::map("world", "USA", plot = FALSE)[c("x", "y")])
canadamap <- data.frame(maps::map("world", "Canada", plot = FALSE)[c("x", "y")])
mexicomap <- data.frame(maps::map("world", "Mexico", plot = FALSE)[c("x", "y")])

#min/max for plotting using output data
MIN <- min(med_fit)
MAX <- max(med_fit)

# #read in breeding/migration range shp file
# setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))
# sp_rng <- rgdal::readOGR(f_out$shp_fname[1], verbose = FALSE)
# 
# #filter by breeding (2) and migration (4) range - need to convert spdf to sp
# nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
# nrng_sp <- sp::SpatialPolygons(nrng@polygons)
# 
# #filter by resident (1) and over winter (3) range - need to convert spdf to sp
# nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]
# nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
# 
# #plotting species range
# nrng@data$id <- rownames(nrng@data)
# nrng.points <- ggplot2::fortify(nrng, region = "id")
# nrng.df <- plyr::join(nrng.points, nrng@data, by = "id")
# 
# nrng_rm@data$id <- rownames(nrng_rm@data)
# nrng_rm.points <- ggplot2::fortify(nrng_rm, region = "id")
# nrng_rm.df <- plyr::join(nrng_rm.points, nrng_rm@data, by = "id")

#median of mu and sd of mu
m_fit <- data.frame(med_beta = med_fit, sd_beta = sd_fit, cell = ot_cl$cell)

#merge hex spatial data with HM data
to_plt <- dplyr::inner_join(m_fit, cell_grid, by = 'cell')
to_plt2 <- dplyr::inner_join(to_plt, ll_df, by = 'cell')


#plot
p_beta <- ggplot() +
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
            xlim = c(-100, -65), ylim = c(23, 55)) +
  geom_polygon(data = to_plt2, aes(x = long, y = lat, group = group, fill = med_beta),
               alpha = 0.5) +
  geom_path(data = to_plt2, aes(x = long, y = lat, group = group),
            alpha = 0.4, color = 'black') +
  scale_fill_gradient2(low = 'indianred', high = 'royalblue', mid = 'lightgoldenrod',
                       limits = c(MIN, MAX), midpoint = 0) +
  # scale_fill_gradientn(colors = c('orange', 'grey', 'light blue'),
  #                        #c(hcl(h = 240, c = 35, l = 35), hcl(h = 180, c = 15, l = 92)),
  #                      limits = c(MIN, MAX)) +
  # scale_fill_gradient2(low = 'indianred2', mid = 'grey60', high = 'lightblue2', 
  #                      limits = c(MIN, MAX), midpoint = 0) +
  labs(fill = 'Slope') +
  # annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg + 0.5,
  #          label = round(to_plt2$med_beta, digits = 2), col = 'black', alpha = 0.2,
  #          size = 3) +
  # annotate('text', x = to_plt2$lon_deg, y = to_plt2$lat_deg - 0.5,
  #          label = round(to_plt2$sd_beta, digits = 2), col = 'white', alpha = 1,
  #          size = 2) +
  ggtitle(paste0(args, ' - Arrival ~ time')) +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')


ggsave(plot = p_beta,
       filename = paste0(args, '-', run_date, '-pheno_trends_map.pdf'))




# plot trends ------------------------------------------------------------

alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')

x_sim <- seq(0, 18, length = 500)

FIT_PLOT <- data.frame(cell = rep(1:DATA$NC, each = length(x_sim)), 
                       x_sim = rep(x_sim, times = DATA$NC),
                       mu_rep_mn = NA, 
                       mu_rep_LCI = NA, 
                       mu_rep_UCI = NA)

counter <- 1
for (i in 1:NCOL(alpha_ch))
{
  #i <- 1
  c_idx <- which(data_f2$cell == ot_cl$cell[i])
  
  yrs <- DATA$year[c_idx]
  min_yr <- min(yrs)
  max_yr <- max(yrs)
  
  for (j in 1:length(x_sim))
  {
    #j <- 29
    if (x_sim[j] >= (min_yr - 1) & x_sim[j] <= (max_yr + 1))
    {
      temp <- alpha_ch[,i] + beta_ch[,i] * x_sim[j]
      tmn <- mean(temp)
      tsd <- sd(temp)
      tLCI <- tmn - tsd
      tUCI <- tmn + tsd
      FIT_PLOT$mu_rep_mn[counter] <- tmn
      FIT_PLOT$mu_rep_LCI[counter] <- tLCI
      FIT_PLOT$mu_rep_UCI[counter] <- tUCI
    } else {
      FIT_PLOT$mu_rep_mn[counter] <- NA
      FIT_PLOT$mu_rep_LCI[counter] <- NA
      FIT_PLOT$mu_rep_UCI[counter] <- NA
    }
    counter <- counter + 1
  }
}


#mean and +- 1 sd - y_out
y_true_mn <- MCMCvis::MCMCpstr(fit, params = 'y_true', func = mean)[[1]]
y_true_sd <- MCMCvis::MCMCpstr(fit, params = 'y_true', func = sd)[[1]]
y_true_LCI <- y_true_mn - y_true_sd
y_true_UCI <- y_true_mn + y_true_sd

DATA_PLOT <- data.frame(y_true_mn, y_true_LCI, y_true_UCI,
                        cell = DATA$cn_id,
                        year = DATA$year,
                        y_obs = DATA$y_obs,
                        y_obs_LCI = DATA$y_obs - DATA$y_sd,
                        y_obs_UCI = DATA$y_obs + DATA$y_sd,
                        lat = data_f2$cell_lat)

pdf(paste0(args, '-', run_date, '-pheno_trends_fig.pdf'), 
    height = 11, width = 9, useDingbats = FALSE)
ggplot(data = DATA_PLOT, aes(DATA$year, y_true_mn), color = 'black', alpha = 0.6) +
  # geom_ribbon(data = FIT_PLOT,
  #             inherit.aes = FALSE,
  #             aes(x = x_sim, ymin = mu_rep_LCI, ymax = mu_rep_UCI,
  #                 group = cell),# color = cell),
  #             #fill = cell,
  #             alpha = 0.1) +
  geom_line(data = FIT_PLOT, aes(x_sim, mu_rep_mn, group = cell, col = factor(cell)),
            alpha = 0.9,
            inherit.aes = FALSE,
            size = 0.8) +
  geom_point(data = DATA_PLOT,
             aes(year, y_true_mn, col = factor(cell)),
             inherit.aes = FALSE, size = 1.5, alpha = 0.8) +
  geom_errorbar(data = DATA_PLOT,
                aes(ymin = y_true_LCI, ymax = y_true_UCI), width = 0.2,
                color = 'black', alpha = 0.4) +
  # geom_point(data = DATA_PLOT,
  #            aes(year, y_obs), col = 'black',
  #            inherit.aes = FALSE, size = 2, alpha = 0.8) +
  # geom_point(data = DATA_PLOT,
  #            aes(year, y_obs, col = factor(cell)),
  #            inherit.aes = FALSE, size = 1.5, alpha = 0.8) +
  # geom_errorbar(data = DATA_PLOT,
  #               aes(ymin = y_obs_LCI, ymax = y_obs_UCI), width = 0.2,
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
dev.off()


pdf(paste0(args, '-', IAR_out_date, '-betas.pdf'), 
    height = 11, width = 9, useDingbats = FALSE)
MCMCvis::MCMCplot(fit, params = 'beta', main = args)
dev.off()



# plot anomalies ----------------------------------------------------------

# out <- data.frame()
# for (i in 1:DATA$NC)
# {
#   #i <- 1
#   temp <- dplyr::filter(DATA_PLOT, cell == i)
#   tsc <- scale(temp$y_true_mn, scale = FALSE)
#   tsc_obs <- scale(temp$y_obs, scale = FALSE)
#   temp$y_true_sc <- tsc
#   temp$y_obs_sc <- tsc_obs
#   temp$year <- temp$year + 2001
#   out <- rbind(out, temp)
# }
# 
# pdf('Vireo_olivaceus_anomaly.pdf', height = 5, width = 7, useDingbats = FALSE)
# ggplot(data = out, aes(year, y_true_sc), color = 'black', alpha = 0.6) +
#   # geom_ribbon(data = FIT_PLOT,
#   #             aes(x = x_sim, ymin = mu_rep_LCI, ymax = mu_rep_UCI),
#   #             fill = 'grey',
#   #             inherit.aes = FALSE,
#   #             alpha = 0.6) +
#   # geom_line(data = FIT_PLOT, aes(x_sim, mu_rep_mn),
#   #           col = 'red',
#   #           alpha = 0.9,
#   #           inherit.aes = FALSE,
#   #           size = 0.8) +
#   # geom_point(data = out,
# #            aes(year, y_true_sc), color = 'black',
# #            inherit.aes = FALSE, size = 2, alpha = 0.7) +
# geom_point(data = out,
#            aes(year, y_true_sc, color = lat), #color = factor(cell)),
#            inherit.aes = FALSE, size = 3, alpha = 0.7) +
#   scale_color_gradient(low = 'red', high = 'white') +
#   # geom_point(data = DATA_PLOT,
#   #            aes(year, y_obs), color = 'black',
#   #            inherit.aes = FALSE, size = 1.5, alpha = 0.5) +
#   # geom_errorbar(data = DATA_PLOT,
#   #               aes(ymin = y_true_LCI, ymax = y_true_UCI), width = 0.4,
#   #               color = 'black', alpha = 0.4) +
#   theme_bw() +
#   #theme(legend.position='none') +
#   xlab('Year') +
#   ylab('Arrival anomaly') +
#   ylim(c(-5, 5)) +
#   theme(
#     axis.text = element_text(size = 16),
#     axis.title = element_text(size = 18),
#     axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
#     axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
#     axis.ticks.length= unit(0.2, 'cm'),) #length of axis tick
# dev.off()



# Trace plots with PPO ----------------------------------------------------

#mu_alpha ~ normal(120, 50)
PR <- rnorm(10000, 120, 50)
MCMCvis::MCMCtrace(fit,
                   params = 'mu_alpha',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_mu_alpha.pdf'))

#sigma_alpha ~ halfnormal(0, 10)
PR_p <- rnorm(10000, 0, 10)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_alpha',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_sigma_alpha.pdf'))

#sigma ~ halfnormal(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_sigma.pdf'))

#alpha_beta ~ normal(0, 10)
PR <- rnorm(10000, 0, 10)
MCMCvis::MCMCtrace(fit,
                   params = 'alpha_beta',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_alpha_beta.pdf'))

#beta_beta ~ normal(0, 1)
PR <- rnorm(10000, 0, 1)
MCMCvis::MCMCtrace(fit,
                   params = 'beta_beta',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_beta_beta.pdf'))

#sigma_beta ~ halfnormal(0, 3)
PR_p <- rnorm(10000, 0, 3)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_beta',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_sigma_beta.pdf'))


if ('Rplots.pdf' %in% list.files())
{
  file.remove('Rplots.pdf')
}


print('I completed!')
