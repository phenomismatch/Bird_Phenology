#interannual variation for cells to check

library(ggplot2)
library(dplyr)

setwd("~/Google_Drive/R/pheno_trends/Data/arrival_master_2019-05-26")
arr_master <- readRDS('arrival_master_2019-05-26.rds')

setwd("~/Google_Drive/R/pheno_trends/Data/juv_master_2019-10-21")
juv_master <- readRDS('juv_master_2019-10-21.rds')


usp <- unique(arr_master$species)
for (i in 1:length(usp))
{
  #i <- 94
  #temp <- dplyr::filter(arr_master, species == usp[i])
  temp <- dplyr::filter(arr_master, species == usp[i], !is.na(mean_pre_IAR))
  hist(temp$sd_pre_IAR)
  hist(aggregate(mean_post_IAR ~ cell, data = temp, FUN = sd)[,2])
  ggplot(temp, aes(year, mean_post_IAR, color = factor(cell))) +
    geom_line() +
    theme(legend.position="none") +
    ylim(c(40, 240))
  
  temp_j <- dplyr::filter(juv_master, species == usp[i], !is.na(juv_mean))
  hist(temp_j$juv_sd)
  hist(aggregate(juv_mean ~ cell, data = temp_j, FUN = sd)[,2])
  ggplot(temp_j, aes(year, juv_mean, color = factor(cell))) +
    geom_line() +
    theme(legend.position="none") +
    ylim(c(40, 240))
}
