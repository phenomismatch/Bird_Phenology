######################
# 10 - juv hitting nets ~ time - models for each species
#
# i = obs
# j = cell
#
# y_{i} \sim N(\mu_{y_{i}}, \sigma_{y_{i}})
# \mu_{y_{i}} \sim N(mu_{i}, \sigma)
# \mu_{i} = \alpha_{i} + \beta_{i} \times year_{i}
# \alpha_{i} \sim N(\mu_{\alpha_{i}}, \sigma_{\alpha})
# \mu_{\alpha_{i}} = \gamma + \theta \times lat_{i}
# \beta_{i} \sim N(\mu_{\beta_{i}}, \sigma_{\beta})
# \mu_{\beta_{i}} = \nu + \eta \times lat_{i}
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# species -----------------------------------------------------------------

#args <- commandArgs(trailingOnly = TRUE)
#args <- as.character('Catharus_minimus')
#args <- as.character('Empidonax_virescens')
#args <- as.character('Vireo_olivaceus')
#args <- as.character('Ammospiza_nelsoni')
#args <- as.character('Agelaius_phoeniceus')
#args <- as.character('Turdus_migratorius')
#args <- as.character('Setophaga_pinus')
args <- as.character('Seiurus_aurocapilla')




# model dir ------------------------------------------------------------

#juveniles MAPS - date input data processed
juv_date <- '2019-10-15'
run_date <- '2019-10-15'

#dirs
juv_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/juv_master_', juv_date)
juv_trends_dir <- paste0(dir, 'Bird_Phenology/Data/Processed/juv_trends_output_', run_date)



# Load packages -----------------------------------------------------------

library(rstan)
library(ggplot2)
library(dplyr)
library(MCMCvis)



# Filter data ------------------------------------------------------------------

#juveniles hitting nets - MAPS

setwd(juv_dir)

#read in 
juvs_master <- readRDS(paste0('juv-output-', juv_date, '.rds'))

#only species/cells/years with data for juvs
j1 <- dplyr::filter(juvs_master, !is.na(juv_mean), species == args)

#cells with at least 3 years of data
cnts <- plyr::count(j1, 'cell')
u_cells <- cnts[which(cnts[,2] >= 3),1]
j2 <- dplyr::filter(j1, cell %in% u_cells)
# j2 <- j1


#add cell lat to df
hexgrid6 <- dggridR::dgconstruct(res = 6)
j2$cell_lat <- dggridR::dgSEQNUM_to_GEO(hexgrid6, 
                                        in_seqnum = j2$cell)$lat_deg
j2$cell_lng <- dggridR::dgSEQNUM_to_GEO(hexgrid6, 
                                        in_seqnum = j2$cell)$lon_deg

#filter by study region
j3 <- dplyr::filter(j2, cell_lng > -95, cell_lat > 24)

#only species that have at least 10 data points
if (NROW(j3) < 10)
{
  stop('Species has fewer than 10 data points')
}

#order cells (and corresponding cell lats) so they match factor
t_cl <- unique(j3[,c('cell', 'cell_lat')])
ot_cl <- t_cl[order(t_cl[,1]),]


DATA <- list(y = j3$juv_mean,
             sd_y = j3$juv_sd,
             year = as.numeric(factor(j3$year)),
             cn_id = as.numeric(factor(j3$cell)),
             NC = NROW(ot_cl),
             # lat = scale(ot_cl$cell_lat, scale = FALSE)[,1],
             lat = ot_cl$cell_lat,
             #lat = scale(j3$cell_lat, scale = FALSE)[,1],
             N = NROW(j3))



# Stan model --------------------------------------------------------------

model <- "
data {
int<lower = 0> N;                                 // number of obs
vector<lower = 0>[N] y;
vector<lower = 0>[N] sd_y;
int<lower = 1> cn_id[N];                          // cell ids
int<lower = 0> NC;                                // number of cells
vector<lower = 1>[N] year;
vector[NC] lat;
}

parameters {
vector[N] mu_y_raw;
vector[NC] alpha_raw;
vector[NC] beta_raw;
vector[NC] mu_alpha_raw;
real<lower = 0> sigma_alpha_raw;
vector[NC] mu_beta_raw;
real<lower = 0> sigma_beta_raw;
real gamma_raw;
real theta_raw;
real pi_raw;
real nu_raw;
real<lower = 0> sigma_raw;
}

transformed parameters {
vector[N] mu;
vector[N] mu_y;
vector[NC] alpha;
vector[NC] beta;
vector[NC] mu_alpha;
real<lower = 0> sigma_alpha;
vector[NC] mu_beta;
real<lower = 0> sigma_beta;
real gamma;
real theta;
real pi;
real nu;
real<lower = 0> sigma;

sigma = sigma_raw * 5;

gamma = gamma_raw * 20 + 200;           // implies gamma ~ N(200, 20)
theta = theta_raw * 3;
mu_alpha = gamma + theta * lat;
sigma_alpha = sigma_alpha_raw * 5;
alpha = alpha_raw * sigma_alpha + mu_alpha;

pi = pi_raw * 1;
nu = nu_raw * 2;
mu_beta = pi + nu * lat;
sigma_beta = sigma_beta_raw * 5;
beta = beta_raw * sigma_beta + mu_beta;

for (i in 1:N)
{
  mu[i] = alpha[cn_id[i]] + beta[cn_id[i]] * year[i];
  mu_y[i] = mu_y_raw[i] * sigma + mu[i];
}
}

model {

mu_y_raw ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();

sigma_beta_raw ~ std_normal();
sigma_alpha_raw ~ std_normal();

gamma_raw ~ std_normal();
theta_raw ~ std_normal();
pi_raw ~ std_normal();
nu_raw ~ std_normal();
sigma_raw ~ std_normal();

y ~ normal(mu_y, sd_y);
}

generated quantities {

real y_rep[N];

// PPC
y_rep = normal_rng(mu_y, sd_y);
}
"

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.99
TREE_DEPTH <- 16
STEP_SIZE <- 0.003
CHAINS <- 4
ITER <- 5000

tt <- proc.time()
fit <- rstan::stan(model_code = model,
                   data = DATA,
                   chains = CHAINS,
                   iter = ITER,
                   cores = CHAINS,
                   pars = c('gamma',
                            'theta',
                            'pi',
                            'nu',
                            'alpha',
                            'beta',
                            'sigma_alpha',
                            'sigma_beta',
                            'sigma',
                            'mu_y',
                            'y_rep'), 
                   control = list(adapt_delta = DELTA,
                                  max_treedepth = TREE_DEPTH,
                                  stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60


num_diverge <- rstan::get_num_divergent(fit)
num_tree <- rstan::get_num_max_treedepth(fit)
num_BFMI <- length(rstan::get_low_bfmi_chains(fit))


# #rerun model with higher target acceptance if divergences exist
# if (num_diverge > 0)
# {
#   DELTA <- 0.99
#   
#   tt <- proc.time()
#   fit <- rstan::stan(model_code = model,
#                      data = DATA,
#                      chains = CHAINS,
#                      iter = ITER,
#                      cores = CHAINS,
#                      pars = c('gamma',
#                               'theta',
#                               'pi',
#                               'nu',
#                               'alpha',
#                               'beta',
#                               'sigma_alpha',
#                               'sigma_beta',
#                               'sigma',
#                               'mu_y',
#                               'y_rep'), 
#                      control = list(adapt_delta = DELTA,
#                                     max_treedepth = TREE_DEPTH,
#                                     stepsize = STEP_SIZE))
#   run_time <- (proc.time()[3] - tt[3]) / 60
#   
#   num_diverge <- rstan::get_num_divergent(fit)
#   num_tree <- rstan::get_num_max_treedepth(fit)
#   num_BFMI <- length(rstan::get_low_bfmi_chains(fit))
# }



#save to RDS
setwd(juv_trends_dir)

saveRDS(fit, file = paste0(args, '-', run_date, '-juv_trends_stan_output.rds'))
saveRDS(DATA, file = paste0(args, '-', run_date, '-juv_trends_stan_input.rds'))



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



# PPC ---------------------------------------------------------------------

y_rep <- MCMCvis::MCMCchains(fit, params = 'y_rep')
y_val <- DATA$y

# bayesplot::ppc_stat(y_val, y_rep, stat = 'mean')
# bayesplot::ppc_dens_overlay(y_val, y_rep[1:100,])

PPC_fun <- function(FUN, YR = y_rep, D = y_val)
{
  out <- sum(apply(YR, 1, FUN) > FUN(D)) / NROW(YR)
  print(out)
}

ppc_mn <- round(PPC_fun(mean), 3)
ppc_min <- round(PPC_fun(min), 3)
ppc_max <- round(PPC_fun(max), 3)




# write model results to file ---------------------------------------------

options(max.print = 5e6)
sink(paste0(args, '-', run_date, '-juv_trends_results.txt'))
cat(paste0('Juv trends results ', args, ' \n'))
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
cat(paste0('Max Rhat: ', max(rhat_output, na.rm = TRUE), ' \n'))
cat(paste0('Min n.eff: ', min(neff_output, na.rm = TRUE), ' \n'))
cat(paste0('PPC mean: ', ppc_mn, ' \n'))
cat(paste0('PPC min: ', ppc_min, ' \n'))
cat(paste0('PPC man: ', ppc_max, ' \n'))
print(model_summary)
sink()




# density overlay plot ----------------------------------------------------

#modified bayesplot::ppc_dens_overlay function
tdata <- bayesplot::ppc_data(DATA$y, y_rep[1:100,])

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
  ggtitle(paste0(args, ' - Fledging ~ time')) +
  theme_bw() +
  xlab('Longitude') +
  ylab('Latitude')


ggsave(plot = p_beta,
       filename = paste0(args, '-', run_date, '-juv_trends_map.pdf'))




# plot trends ------------------------------------------------------------

alpha_ch <- MCMCvis::MCMCchains(fit, params = 'alpha')
beta_ch <- MCMCvis::MCMCchains(fit, params = 'beta')

x_sim <- seq(0, max(DATA$year), length = 500)

FIT_PLOT <- data.frame(cell = rep(1:DATA$NC, each = length(x_sim)), 
                       x_sim = rep(x_sim, times = DATA$NC),
                       mu_rep_mn = NA, 
                       mu_rep_LCI = NA, 
                       mu_rep_UCI = NA)


counter <- 1
for (i in 1:NCOL(alpha_ch))
{
  #i <- 1
  c_idx <- which(j2$cell == ot_cl$cell[i])
  
  yrs <- DATA$year[c_idx]
  min_yr <- min(yrs)
  max_yr <- max(yrs)
  
  for (j in 1:length(x_sim))
  {
    #j <- 1
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
y_true_mn <- MCMCvis::MCMCpstr(fit, params = 'mu_y', func = mean)[[1]]
y_true_sd <- MCMCvis::MCMCpstr(fit, params = 'mu_y', func = sd)[[1]]
y_true_LCI <- y_true_mn - y_true_sd
y_true_UCI <- y_true_mn + y_true_sd

DATA_PLOT <- data.frame(y_true_mn, 
                        y_true_LCI, 
                        y_true_UCI,
                        cell = DATA$cn_id,
                        year = DATA$year,
                        y_obs = DATA$y,
                        y_obs_LCI = DATA$y - DATA$sd_y,
                        y_obs_UCI = DATA$y + DATA$sd_y)

pdf(paste0(args, '-', run_date, '-juv_trends_fig.pdf'), 
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
  ylab('Fledge date') +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 15, r = 15, b = 0, l = 0)),
    axis.ticks.length= unit(0.2, 'cm')) #length of axis tick
dev.off()


pdf(paste0(args, '-', run_date, '-betas.pdf'), 
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

#sigma ~ halfnormal(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_sigma.pdf'))

#gamma ~ normal(200, 20)
PR <- rnorm(10000, 200, 20)
MCMCvis::MCMCtrace(fit,
                   params = 'gamma',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_gamma.pdf'))

#theta ~ normal(0, 3)
PR <- rnorm(10000, 0, 3)
MCMCvis::MCMCtrace(fit,
                   params = 'theta',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_theta.pdf'))

#sigma_alpha ~ halfnormal(0, 5)
PR_p <- rnorm(10000, 0, 5)
PR <- PR_p[which(PR_p > 0)]
MCMCvis::MCMCtrace(fit,
                   params = 'sigma_alpha',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_sigma_alpha.pdf'))

#pi ~ normal(0, 1)
PR <- rnorm(10000, 0, 1)
MCMCvis::MCMCtrace(fit,
                   params = 'pi',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_pi.pdf'))

#nu ~ normal(0, 2)
PR <- rnorm(10000, 0, 2)
MCMCvis::MCMCtrace(fit,
                   params = 'nu',
                   priors = PR,
                   open_pdf = FALSE,
                   filename = paste0(args, '-', run_date, '-trace_nu.pdf'))

#sigma_beta ~ halfnormal(0, 5)
PR_p <- rnorm(10000, 0, 5)
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

