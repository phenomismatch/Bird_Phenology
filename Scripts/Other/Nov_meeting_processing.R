#processing for Nov Phenomismatch meeting

#2001-2016 for your hexes: 595, 676, 729 for any birds you think might be particularly informative


# arrival -----------------------------------------------------------------

setwd('~/Google_Drive/R/pheno_trends/Data/arrival_master_2019-05-26/')

tdata <- readRDS('arrival_master_2019-05-26.rds')
str(tdata)

tcell <- c(595, 676, 729)
fdata <- dplyr::filter(tdata, cell %in% tcell, year <= 2016, year >= 2001)
fdata2 <- dplyr::filter(fdata, min_neff > 200, max_rhat < 1.1)

fdata3 <- dplyr::select(fdata2, species, cell, year, 
                       mig_cell, breed_cell, cell_lat, 
                       cell_lon, mean_pre_IAR, sd_pre_IAR, 
                       mean_post_IAR, sd_post_IAR, F_prcp, 
                       M_prcp, A_prcp, FMA_prcp, F_tmax, 
                       M_tmax, A_tmax, FMA_tmax, F_tmin, 
                       M_tmin, A_tmin, FMA_tmin)

setwd('~/Desktop/')
write.csv(fdata3, 'Bird_arrival_2019-11-06.csv')


# breeding ----------------------------------------------------------------

setwd('~/Google_Drive/R/pheno_trends/Data/juv_master_2019-10-21/')

tdata <- readRDS('juv_master_2019-10-21.rds')
str(tdata)

tcell <- c(595, 676, 729)
fdata <- dplyr::filter(tdata, cell %in% tcell, year <= 2016, year >= 2001)
fdata2 <- dplyr::filter(fdata, min_neff > 200, max_Rhat < 1.1, !is.na(juv_mean))

fdata3 <- dplyr::select(fdata2, species, cell, year, 
                        cell, juv_mean, juv_sd)

setwd('~/Desktop/')
write.csv(fdata3, 'Bird_fledge_2019-11-06.csv')


