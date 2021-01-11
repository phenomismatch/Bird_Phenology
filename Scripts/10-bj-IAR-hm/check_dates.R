##########################################
# look at estimates for cells with MAPS vs. eBird data
#
# MAPS estimates seem to be much later than eBird estimates
# seems to be due to lack of clear 'hump' in data
##########################################

library(ggplot2)
library(dplyr)

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/bj_master_2020-12-07")
bj_master <- readRDS('bj_master_2020-12-07.rds')
f_idx <- which(bj_master$VALID_br_GAM == FALSE & bj_master$VALID_juv_GAM == FALSE)
bj_master$bj_IAR_mean[f_idx] <- NA
bj_master$bj_IAR_sd[f_idx] <- NA

env_date <- '2020-08-06'
setwd(paste0('~/Google_Drive/R/pheno_trends/Data/environment/processed/', env_date))
gr <- readRDS(paste0('MidGreenup-', env_date, '-forest.rds'))

#merge with greenup data
df_gr <- which(colnames(gr) %in% c('cell_lat', 'cell_lng'))
mrg_f <- dplyr::left_join(bj_master, gr[,-df_gr], by = c('cell', 'year'))

#only breeding cells (not migratory, resident, or non-breeding)
#only rows that have greenup data
#only cells that had > 10000 cells with valid greenup data (from forest pixels) - 2500 km^2 greenup data
mrg_f2 <- dplyr::filter(mrg_f, !is.na(gr_mn), 
                        gr_ncell > 10000, other_cell == FALSE)
cn_mrg_f2 <- colnames(mrg_f2)
sp_k <- unique(mrg_f2$species)

#what proportion of species range are there data available for
fd <- dplyr::filter(mrg_f2, !is.na(bj_IAR_mean))
#number total cells
ntc <- aggregate(cell ~ species, mrg_f2, function(x) length(unique(x)))
#total lat range
ntl <- aggregate(cell_lat ~ species, mrg_f2, function(x) max(x) - min(x))
#number obs each cell/species
cnt_csp <- plyr::count(fd, c('species', 'cell', 'cell_lat'))
ge3_csp <- dplyr::filter(cnt_csp, freq >= 3)
#number cells greater than or equal to 3
ncd <- aggregate(cell ~ species, ge3_csp, function(x) length(unique(x)))
#range of cells greater than or equal to 3
ncl <- aggregate(cell_lat ~ species, ge3_csp, function(x) max(x) - min(x))
#merge
jj <- dplyr::left_join(ntc, ncd, by = 'species')
jj2 <- dplyr::left_join(jj, ntl, by = 'species')
jj3 <- dplyr::left_join(jj2, ncl, by = 'species')
#proportion cell data
jj3$pcd <- jj3[,3] /  jj3[,2]
#proportion lat range
jj3$lrng <- jj3[,5] /  jj3[,4]
colnames(jj3) <- c('species', 'n_cells', 'n_cells_ge_3', 'rng_lat', 
                   'rng_lat_ge3', 'prop_cell', 'prop_lat')
tsid <- which(jj3$n_cells_ge_3 >= 3)
sp_k2 <- sp_k[tsid]

#filter species
mrg_f3 <- dplyr::filter(mrg_f2, species %in% sp_k2)



#filter data
bj_data <- dplyr::filter(mrg_f3, !is.na(bj_IAR_mean))
br_data <- dplyr::filter(mrg_f3, !is.na(br_GAM_mean))

#fit models
bj_fit <- lm(bj_IAR_mean ~ gr_mn, data = bj_data)
bj_data$resid <- residuals(bj_fit)
bj_fit2 <- lm(bj_IAR_mean ~ gr_mn, data = br_data)


br_fit <- lm(br_GAM_mean ~ gr_mn, data = br_data)
br_data$resid <- residuals(br_fit)


#plots (br ~ GR)
#BJ data
#GREEN = cells with MAPS data
bj_IAR <- dplyr::filter(bj_data, VALID_juv_GAM != TRUE)
# ggplot(bj_IAR, aes(gr_mn, bj_IAR_mean, color = factor(species))) +
#         geom_point() +
#         xlim(c(85, 155)) +
#         ylim(c(115, 205)) +
#         theme(legend.position = 'none')

plot(bj_IAR$gr_mn, bj_IAR$bj_IAR_mean, 
     xlab = 'Greenup', ylab = 'IAR-derived fledge date',
     main = 'IAR-derived data',
     pch = 19, col = rgb(0,0,0,0.2),
     xlim = c(85, 155), ylim = c(115, 205))
abline(bj_fit, col = 'red')
abline(bj_fit2, col = 'green')
bj_juv <- dplyr::filter(bj_data, VALID_juv_GAM == TRUE)
points(bj_juv$gr_mn, bj_juv$bj_IAR_mean, col = rgb(0,1,0,0.4), pch = 19)

#BR data
plot(br_data$gr_mn, br_data$br_GAM_mean, 
     xlab = 'Greenup', ylab = 'GAM-derived fledge date (from eBird data)',
     main = 'GAM-derived data',
     pch = 19, col = rgb(0,0,0,0.2),
     xlim = c(85, 155), ylim = c(115, 205))
abline(br_fit, col = 'red')

# dplyr::filter(mrg_f3, VALID_br_GAM == TRUE, VALID_juv_GAM == TRUE)
# plot(bj_data$br_GAM_mean, bj_data$juv_GAM_mean, 
#      pch = 19)

# number data points for br GAM per species
plyr::count(br_data, 'species')
# number data points for juv GAM per species - V griseus and G formosa are only juv
plyr::count(bj_juv, 'species')
#check per species basis
#tsp <- dplyr::filter(bj_data, species == 'Empidonax_virescens')
#tsp <- dplyr::filter(bj_data, species == 'Troglodytes_aedon')
tsp <- dplyr::filter(bj_data, species == 'Vireo_griseus')
plot(tsp$gr_mn, tsp$bj_IAR_mean, pch = 19)
tsp_juv <- dplyr::filter(tsp, VALID_juv_GAM)
points(tsp_juv$gr_mn, tsp_juv$bj_IAR_mean, pch = 19, col = 'green')


#read in juv GAM data (to get plmax)
setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/juv_master_2020-12-04")
juv_master <- readRDS('juv_master_2020-12-04.rds')
dplyr::filter(juv_master, year >= 2002, year <= 2017, per_ovr >= 0.05, 
              plmax > 0.99, breed_cell == TRUE, other_cell == FALSE)
tt2 <- dplyr::select(juv_master, species, year, cell, plmax)
bj_juv2 <- dplyr::left_join(bj_juv, tt2)
#look at residuals in relation to plmax
plot(bj_juv2$plmax, bj_juv2$resid, 
     xlab = 'Proportion of GAM iterations with local maximum',
     ylab = 'Residual from breeding ~ greenup model above',
     pch = 19, col = rgb(0,0,0,0.5))

#only three MAPS with plmax > 0.99
dplyr::filter(bj_juv2, plmax > 0.99)
