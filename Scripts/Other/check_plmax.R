
# load datra --------------------------------------------------------------

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/arrival_master_2020-04-15")
arr_master <- readRDS('arrival_master_2020-04-15.rds')

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/IAR_input_2020-02-26")
GAM_master <- readRDS('IAR_input-2020-02-26.rds')

setwd("~/Google_Drive/R/Bird_Phenology/Data/Processed/juv_master_2020-03-25")
juv_master <- readRDS('juv_master_2020-03-25.rds')

# process -----------------------------------------------------------------

head(arr_master)

arr2 <- dplyr::filter(arr_master, !is.na(arr_GAM_mean))
arr2$diff <- arr2$arr_IAR_mean - arr2$arr_GAM_mean

mrg2 <- dplyr::left_join(arr2, GAM_master, 
                         by = c('species', 'cell', 'year'))

mrg3 <- dplyr::select(mrg2, species, cell, year, cell_lat, 
                      arr_GAM_mean, arr_GAM_sd, 
                      arr_IAR_mean, arr_IAR_sd,
                      PPC_mn_bias, PPC_mn_pval,
                      diff, mlmax, plmax)




# explore -----------------------------------------------------------------

sp_diff <- aggregate(diff ~ species, mrg3, mean)
sp_PPC <- aggregate(PPC_mn_pval ~ species, mrg3, mean)
sp_mrg <- dplyr::left_join(sp_diff, sp_PPC)

plot(sp_mrg$diff, sp_mrg$PPC_mn_pval, col = rgb(0,0,0,0.3), pch = 19)
abline(h = 0.05, col = 'red', lty = 2)
abline(h = 0.95, col = 'red', lty = 2)


sp_mrg
sp <- 'Ictinia_mississippiensis'
sp <- 'Cardellina_pusilla'
sp <- 'Catharus_minimus'
sp <- 'Coccyzus_americanus'
sp <- 'Contopus_virens'
sp <- 'Dumetella_carolinensis'
sp <- 'Empidonax_virescens'
sp <- 'Geothlypis_trichas'
sp <- 'Hirundo_rustica'
sp <- 'Icterus_spurius'
sp <- 'Oreothlypis_celata'
sp <- 'Oreothlypis_peregrina'
sp <- 'Parkesia_motacilla'
sp <- 'Passerina_caerulea'
sp <- 'Passerina_cyanea'
sp <- 'Petrochelidon_pyrrhonota'
sp <- 'Pooecetes_gramineus'
sp <- 'Protonotaria_citrea'
sp <- 'Regulus_calendula'
sp <- 'Sayornis_phoebe'
sp <- 'Setophaga_americana'
sp <- 'Setophaga_castanea'
sp <- 'Setophaga_discolor'
sp <- 'Setophaga_dominica'
sp <- 'Setophaga_palmarum'
sp <- 'Setophaga_striata'
sp <- 'Setophaga_tigrina'
sp <- 'Spizella_passerina'
sp <- 'Stelgidopteryx_serripennis'
sp <- 'Tachycineta_bicolor'
sp <- 'Troglodytes_aedon'
sp <- 'Tyrannus_tyrannus'
sp <- 'Vermivora_cyanoptera'
sp <- 'Vireo_flavifrons'
sp <- 'Vireo_griseus'
sp <- 'Vireo_olivaceus'

tt <- dplyr::filter(mrg3, species == sp)
dplyr::filter(mrg3, species == sp, diff > 8)

plot(tt$arr_GAM_sd, tt$diff)
#plot(tt$arr_GAM_sd, tt$plmax)
plot(tt$plmax, tt$diff)
abline(h = 0, col = 'red', lty = 2)

library(ggplot2)
ggplot(mrg3, aes(arr_GAM_sd, diff)) +
  geom_point() +
  geom_smooth()

plot(mrg3$arr_GAM_sd, mrg3$diff, #pch = '.', 
     col = rgb(0,0,0,0.1))
summary(lm(diff ~ arr_GAM_sd + plmax, data = mrg3))
abline(h = 0, col = 'red', lty = 2)
plot(mrg3$plmax, mrg3$diff, #pch = '.', 
     col = rgb(0,0,0,0.1))#, 
     #xlim = c(0.8, 1))
abline(h = 0, col = 'red', lty = 2)
abline(v = 0.9, col = 'blue', lty = 2)


dplyr::filter(mrg3, diff < -30, plmax > 0.8)
NROW(dplyr::filter(mrg3, plmax > 0.9))
NROW(mrg3)




jmf <- dplyr::filter(juv_master, !is.na(juv_mean))
NROW(jmf)
jmf2 <- dplyr::filter(jmf, plmax > 0.8)
NROW(jmf2)

dplyr::filter(plyr::count(jmf2, c('species', 'cell')), freq >= 3)
NROW(dplyr::filter(juv_master, !is.na(juv_mean), plmax > 0.8))

dplyr::inner_join(mrg3, jmf)


mrg4 <- dplyr::filter(mrg3, plmax > 0.9)
plot(mrg4$arr_GAM_sd, mrg4$diff, #pch = '.', 
     col = rgb(0,0,0,0.1))
abline(h = 0, col = 'red', lty = 2)
plot(mrg4$plmax, mrg4$diff, #pch = '.', 
     col = rgb(0,0,0,0.1), ylim = c(-20, 20))#, 


#conditions
plmax > 0.8
arr_GAM_sd < 10
