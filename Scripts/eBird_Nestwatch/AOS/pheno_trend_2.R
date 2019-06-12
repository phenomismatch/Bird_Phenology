#calc var for each cell (interannual variability)


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

fit <- readRDS('Vireo_olivaceus-2019-05-26-iar-stan_output.rds')
fit_in <- readRDS('Vireo_olivaceus-2019-05-26-iar-stan_input.rds')

yt_ch_p <- MCMCvis::MCMCpstr(fit, params = 'y_true', type = 'chains')[[1]]

cell_num <- sort(unique(d_ind2$row))

#filter cell lats
cl <- cells[cell_num,]



# calc sd among years -----------------------------------------------------

#for each iter, sd among years for each cell
sd_iter <- c()
for (i in 1:length(cell_num))
{
  #i <- 2
  #get value for each iter for every year at specified cell
  tt <- yt_ch_p[cell_num[i], , ]
  #for each iter calculate sd among years for specified cell
  t_sd <- apply(tt, 2, sd)

  #hist(t_sd, xlim = c(0, 7), breaks = 20)
  sd_iter <- cbind(sd_iter, t_sd)
}

colnames(sd_iter) <- cl$cell

#cell numbers and lat sorted by lat
cl_srt <- cl[order(cl[,2]),]



# slope estimate for each iter --------------------------------------------

sl_out <- c()
for (j in 1:NROW(sd_iter))
{
  #j <- 1
  sl <- summary(lm(sd_iter[j,] ~ cl$cell_lat))$coef[2,1]
  sl_out <- c(sl_out, sl)
}

hist(sl_out)


# plot interannual variabilty ---------------------------------------------

plot(rep(cl$cell_lat[1], 24000), sd_iter[,1], pch = '.',
     col = rgb(0,0,0,0.2), xlim = range(cl$cell_lat), ylim = c(0,8),
     xlab = 'Latitude', ylab = 'Inerannual variability (sd)')
for (i in 2:NCOL(sd_iter))
{
  #i <- 2
  points(rep(cl$cell_lat[i], 24000), sd_iter[,i], pch = '.', 
         col = rgb(0,0,0,0.2))
}

#same with MCMCplot
MCMCplot(sd_iter, params = cl_srt$cell, 
         horiz = FALSE, labels = cl_srt$cell_lat, ref = NULL)

