#calc var for each cell (interannual variability)
#more variation at bottom edge of range, it seems


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



# args --------------------------------------------------------------------

species <- 'Vireo_olivaceus'


# read in data produced by 5 ----------------------------------------------------

setwd('~/Desktop/')

data <- readRDS('Vireo_olivaceus_pro_IAR.rds')

yrs <- unique(data$year)
cells <- unique(data[,c('cell', 'cell_lat')])

#which cell/years have data
d_ind <- data.frame()
for (i in 1:length(yrs))
{
  for (j in 1:length(cells[,1]))
  {
    #i <- 1
    #j <- 1
    temp <- dplyr::filter(data, year == yrs[i], cell == cells[j,1])
    
    if (!is.na(temp$mean_pre_IAR))
    {
      tt <- data.frame(row = j, col = i)
      d_ind <- rbind(d_ind, tt)
    }
  }
}

#indicices to feed MCMCvis
ext_ind <- paste0('y_true\\[', d_ind[,1], ',', d_ind[,2], '\\]')


# read in IAR fit data ----------------------------------------------------

setwd('~/Desktop/Bird_Phenology_Offline/Data/Processed/IAR_output_2019-05-26/')

fit <- readRDS(paste0(species, '-2019-05-26-iar-stan_output.rds'))

yt_ch_p <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'chains')[[1]]

#filter cell lats
cell_num <- sort(unique(d_ind2$row))
cl <- cells[cell_num,]



# read in GAM fit data ----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_species_2019-05-03'))
GAM_halfmax <- readRDS(paste0('halfmax_df_arrival_', species, '.rds'))



# calc sd among years -----------------------------------------------------

#for each iter, sd among years for each cell
sd_iter <- c()
dsd_out <- c()
GAM_sd <- c()
sd_iter_yr <- c()
for (i in 1:length(cell_num))
{
  #i <- 11
  #get value for each iter for every year at specified cell
  tt <- yt_ch_p[cell_num[i], , ]
  #only years which there were data
  yr_idx <- dplyr::filter(d_ind, row == cell_num[i])[,2]
  tt2 <- yt_ch_p[cell_num[i], yr_idx, ]
  
  #IAR input sd
  # dsd <- sd(dplyr::filter(data, cell == cl$cell[i])$mean_pre_IAR, na.rm = TRUE)
  # dsd_out <- c(dsd_out, dsd)
  
  #GAM sd
  tt_GAM <- dplyr::filter(GAM_halfmax, cell == cl$cell[i])
  halfmax <- as.matrix(tt_GAM[,grep('iter', colnames(tt_GAM))])
  hm_mns <- apply(halfmax, 1, mean)
  t_GAM_sd <- sd(hm_mns, na.rm = TRUE)
  GAM_sd <- c(GAM_sd, t_GAM_sd)
  
  #for each iter calculate sd among years for specified cell
  if (NCOL(tt) > 1)
  {
    t_sd <- apply(tt, 2, sd)
    sd_iter <- cbind(sd_iter, t_sd)
  } else {
    sd_iter <- cbind(sd_iter, rep(NA, 24000))
  }
  if (NCOL(tt2) > 2)
  {
    t2_sd <- apply(tt2, 2, sd)
    sd_iter_yr <- cbind(sd_iter_yr, t2_sd)
  } else {
    sd_iter_yr <- cbind(sd_iter_yr, rep(NA, 24000))
  }

  #hist(t_sd, xlim = c(0, 7), breaks = 20)
}

colnames(sd_iter) <- cl$cell
colnames(sd_iter_yr) <- cl$cell

#cell numbers and lat sorted by lat
cl_srt <- cl[order(cl[,2]),]



# slope estimate for each iter --------------------------------------------

sl_out <- c()
sl_out2 <- c()
for (j in 1:NROW(sd_iter))
{
  #j <- 1
  sl <- summary(lm(sd_iter[j,] ~ cl$cell_lat))$coef[2,1]
  sl_out <- c(sl_out, sl)
  
  sl2 <- summary(lm(sd_iter_yr[j,] ~ cl$cell_lat))$coef[2,1]
  sl_out2 <- c(sl_out2, sl2)
}


# #histogram fig
# pdf(paste0(species, '_sd_lat_hist.pdf'))
# hist(sl_out, xlab = 'Slope for sd ~ lat', main = paste0(species))
# abline(v = 0, lty = 2, lwd = 4)
# dev.off()
# 
# hist(sl_out2, xlab = 'Slope for sd ~ lat', main = paste0(species))
# abline(v = 0, lty = 2, lwd = 4)



# plot interannual variabilty ---------------------------------------------

#black points are values at each iteration (only years with data)
#red points are mean
pdf(paste0(species, '_sd_lat_plot.pdf'))
plot(rep(cl$cell_lat[1], 24000), sd_iter_yr[,1], pch = '.',
     col = rgb(0,0,0,0.2), xlim = range(cl$cell_lat), ylim = c(0,8),
     xlab = 'Latitude', ylab = 'Inerannual variability (sd)',
     main = species)
points(cl$cell_lat[1], mean(sd_iter_yr[,1]), cex = 0.5, pch = 19, col = rgb(1,0,0,0.3))
for (i in 2:NCOL(sd_iter_yr))
{
  #i <- 2
  points(rep(cl$cell_lat[i], 24000), sd_iter_yr[,i], pch = '.', 
         col = rgb(0,0,0,0.2))
  points(cl$cell_lat[i], mean(sd_iter_yr[,i]), cex = 0.5, pch = 19, col = rgb(1,0,0,0.3))
}
dev.off()

#same with MCMCplot
MCMCplot(sd_iter_yr, params = cl_srt$cell,
         horiz = FALSE, labels = cl_srt$cell_lat, ref = NULL)

