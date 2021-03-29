######################
# Analyze Nestwatch
#
# multimodal distribution of pheno events
#
# pheno period interval length
######################


# load packages -----------------------------------------------------------

library(tidyverse)
library(Rphylopars)
library(ape)


# date --------------------------------------------------------------------

NW_QUERY_DATE <- '2020-10-28'
RUN_DATE <- '2020-12-03'


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/Bird_Phenology'

setwd(paste0(dir, "/Data/Processed"))
nw_data <- readRDS(paste0('Nestwatch-pro-', NW_QUERY_DATE, '.rds'))

#load arr data
ARR_DATE <- '2020-07-21'
setwd(paste0(dir, "/Data/Processed/arrival_master_", ARR_DATE))
arr_data <- readRDS(paste0('arrival_master_', ARR_DATE, '.rds'))

#fill NA where GAM not valid
arr_data$arr_IAR_mean[which(arr_data$VALID_GAM == FALSE)] <- NA
arr_data$arr_IAR_sd[which(arr_data$VALID_GAM == FALSE)] <- NA

arr_data2 <- dplyr::select(arr_data, species, cell, mig_cell, breed_cell, 
                           cell_lat, cell_lng, year, arr_IAR_mean)


# count num obs --------------------------------------------------------------

#usp_arr <- unique(arr_data2$species)
#nw_data2 <- dplyr::filter(nw_data, species %in% usp_arr)
nw_data2 <- nw_data

usp <- sort(unique(nw_data2$species))

#filter by only valid lay dates
sp_cnt <- dplyr::filter(nw_data, !is.na(lay)) %>%
  dplyr::count(species)

#hist(sp_cnt[,2], breaks = 10000, xlim = c(0, 100))

#only species with at least 40 observations
usp2 <- sp_cnt[which(sp_cnt[,2] > 40),1]


# distribution of pheno dates (hist) -----------------------------------------------

#LAY
#dev.off()
par(mfrow = c(3,3))
for (i in 1:length(usp2))
{
  #i <- 36
  temp <- dplyr::filter(nw_data2, species == usp2[i])
  if (sum(!is.na(temp$lay)) > 0)
  {
    #hist(temp$lay, main = paste0('LAY - ', usp2[i]))
    plot(temp$lat, temp$lay, pch = 19, 
         col = rgb(0,0,0,0.015), 
         main = paste0('LAY ~ LAT - ', usp2[i]),
         xlab = 'Latitude', ylab = 'Lay date')
  }
}

for (i in 1:length(usp2))
{
  #i <- 1
  temp <- dplyr::filter(nw_data, species == usp2[i])
  if (sum(!is.na(temp$hatch)) > 0)
  {
    #hist(temp$hatch, main = paste0('HATCH - ', usp2[i]))
    plot(temp$lat, temp$lay, pch = 19, 
         col = rgb(0,0,0,0.05), 
         main = paste0('HATCH ~ LAT - ', usp2[i]))
  }
}

for (i in 1:length(usp2))
{
  #i <- 39
  temp <- dplyr::filter(nw_data, species == usp2[i])
  if (sum(!is.na(temp$fledge)) > 0)
  {
    #hist(temp$fledge, main = paste0('FLEDGE - ', usp2[i]))
    plot(temp$lat, temp$lay, pch = 19, 
         col = rgb(0,0,0,0.05), 
         main = paste0('FLEDGE ~ LAT - ', usp2[i]))
  }
}

#can test multimodality using likelihood ratio test
#https://stats.stackexchange.com/questions/138223/how-to-test-if-my-distribution-is-multimodal


# interval length ---------------------------------------------------------

dd <- data.frame(species = rep(NA, length(usp2)),
                 med_LH = NA,
                 med_LF = NA)

par(mfrow = c(2,2))
for (i in 1:length(usp2))
{
  #i <- 1
  temp <- dplyr::filter(nw_data2, species == usp2[i])
  temp_hl <- temp$hatch - temp$lay
  temp_fl <- temp$fledge - temp$lay
  
  dd$species[i] <- usp2[i]
  dd$med_LH[i] <- median(temp_hl, na.rm = TRUE)
  dd$med_LF[i] <- median(temp_fl, na.rm = TRUE)
  
  # if (sum(!is.na(temp_hl)) > 0)
  # {
  #   hist(temp_hl, main = paste0('Lay -> Hatch - ', usp2[i]))
  #   abline(v = median(temp_hl, na.rm = TRUE), col = 'red', lwd = 2)
  # } else {
  #   hist(0, main = paste0('Lay -> Hatch - ', usp2[i]))
  # }
  # 
  # if (sum(!is.na(temp_fl)) > 0)
  # {
  #   hist(temp_fl, main = paste0('Lay -> Fledge - ', usp2[i]))
  #   abline(v = median(temp_fl, na.rm = TRUE), col = 'red', lwd = 2)
  # } else {
  #   hist(0, main = paste0('Lay -> Fledge - ', usp2[i]))
  # }
}


sp_med_LH <- median(dd$med_LH, na.rm = TRUE)
sp_med_LF <- median(dd$med_LF, na.rm = TRUE)

hist(dd$med_LH)
abline(v = sp_med_LH, col = 'red', lwd = 2)
hist(dd$med_LF)
abline(v = median(sp_med_LF, na.rm = TRUE), col = 'red', lwd = 2)


# pheno periods -----------------------------------------------------------------

#on a cell/year basis, bc of possible compression between arr and lay over lat

# PERIODS
#========
# impute missing periods with phylogenetic relationship
# ARRIVAL: IAR arr model
# FLEDGE: IAR br model
# LAY: FLEDGE - interval from NW
# HATCH: FLEDGE - interval from NW

'%ni%' <- Negate('%in%')

#arr species
arr_sp <- unique(arr_data2$species)
dd_sp <- unique(dd$species)

#missing sp in dd
miss_sp <- arr_sp[which(arr_sp %ni% dd_sp)]

#new df with all species
dd2 <- data.frame(species = c(dd$species, miss_sp), 
                  med_LH = c(dd$med_LH, rep(NA, length.out = length(miss_sp))),
                  med_LF = c(dd$med_LF, rep(NA, length.out = length(miss_sp))))


# impute using phylogentic relationship -----------------------------------

#exclude the following species (non-near passerines)
to.rm <- c('Aix_sponsa', 'Branta_canadensis', 'Buteo_jamaicensis', 'Falco_sparverius', 
           'Hirundo_rustica_erythrogaster', 'Larus_argentatus', 'Lophodytes_cucullatus', 
           'Megascops_asio', 'Poecile_carolinensis/atricapillus', 
           'Troglodytes_aedon_[aedon_Group]', 'Charadrius_vociferus')

to.rm.idx <- which(dd2$species %in% to.rm)
dd3 <- dd2[-to.rm.idx, ]

setwd("~/Google_Drive/R/Bird_Phenology/Data/NW_bird_phylo")
bt_names <- read.csv('nw_birdtree_names.txt', header = FALSE)
dd3$bt_species <- gsub(' ', '_', bt_names[,1])

#read in tree (100 trees)
#tree downlaoded from: http://birdtree.org
setwd("~/Google_Drive/R/Bird_Phenology/Data/NW_bird_phylo/tree-pruner-4daa14d3-e2e8-4730-bcf0-af6ced22cae0")
tree <- ape::read.nexus('output.nex')

#impute at each tree
#https://github.com/ericgoolsby/Rphylopars/wiki/Example-1:-Getting-Started
LH_out <- data.frame(species = dd3$species, bt_species = dd3$bt_species)
LF_out <- data.frame(species = dd3$species, bt_species = dd3$bt_species)
for (i in 1:length(tree))
{
  #i <- 2
  tt <- tree[[i]]
  tnames <- tt$tip.label
  
  #order names in df
  tt2 <- left_join(data.frame(species = tnames), dd3[,-1], 
                   by = c('species' = 'bt_species'))
  
  #fit multivariate Brown motion model using traits (both) and tree
  pBM <- Rphylopars::phylopars(trait_data = tt2, tree = tt)
  
  #imputed means
  i_mns <- pBM$anc_recon[1:NROW(tt2),]
  #i_v <- sqrt(pBM$anc_var[1:NROW(tt2),])*1.96
  
  # #just to check: differences between true and imputed are 0
  # pBM$trait_data[,2] - i_mns[,1]
  # pBM$trait_data[,3] - i_mns[,2]
  LH_out <- dplyr::left_join(LH_out, data.frame(bt_species = row.names(i_mns), 
                                                med_LH = i_mns[,1]), by = 'bt_species')
  LF_out <- dplyr::left_join(LF_out, data.frame(bt_species = row.names(i_mns), 
                                                med_LF = i_mns[,2]), by = 'bt_species')
}

colnames(LH_out) <- c('species', 'bt_species', paste0('med_LH_t', 1:length(tree)))
colnames(LF_out) <- c('species', 'bt_species', paste0('med_LF_t', 1:length(tree)))

#take mean across all trees
dd3$med_LH_imp <- round(apply(LH_out[,-c(1:2)], 1, mean))
dd3$med_LF_imp <- round(apply(LF_out[,-c(1:2)], 1, mean))

dd4 <- dplyr::select(dd3, species, bt_species, med_LH, med_LF, med_LH_imp, med_LF_imp)

setwd(paste0(dir, '/Data/Processed/'))
saveRDS(dd4, paste0('Nestwatch_pheno_dates-', RUN_DATE, '.rds'))

#just species from arr_master
setwd(paste0(dir, '/Data/'))
write.table(unique(arr_data$species), 'arr_species_list.txt', 
            row.names = FALSE, col.names = FALSE)

#use this in 8-br-GAM model

#impute missing values with medians - alternative to phylogenetic imputation
# na_med_LH <- which(is.na(dd2$med_LH))
# na_med_LF <- which(is.na(dd2$med_LF))
# dd2$med_LH[na_med_LH] <- sp_med_LH
# dd2$med_LF[na_med_LF] <- sp_med_LF


# pheno periods per nest basis --------------------------------------------

nw_data$hl <- NA
nw_data$fh <- NA
for (i in 1:NROW(nw_data))
{
  #i <- 1
  print(paste0('Processing : ', i, ' of ', NROW(nw_data)))
  
  nw_data$hl[i] <- nw_data[i,'hatch'] - nw_data[i,'lay']
  nw_data$fh[i] <- nw_data[i,'fledge'] - nw_data[i,'hatch']
}

usp <- unique(nw_data$species)
td <- dplyr::filter(nw_data, hl > 0, fh > 0)
cnt <- plyr::count(td, 'species')
nsp <- cnt[which(cnt$freq > 400), 'species']

nsp2 <- nsp[-which(nsp == 'Aix_sponsa')]
td2 <- dplyr::filter(td, species %in% nsp2)
td2$bs <- td2$n_fledged / td2$n_eggs
td2$nc <- td2$n_chicks / td2$n_eggs

ggplot(td2, aes(lat, hl, color = factor(species))) +
#ggplot(td2, aes(lat, fh, color = factor(species))) +
#ggplot(td2, aes(hl, bs, color = factor(species))) +
#ggplot(td2, aes(fh, bs, color = factor(species))) +
#ggplot(td2, aes(hl, nc, color = factor(species))) +
  geom_point(alpha = 0.05) +
  geom_line(stat = 'smooth', method = 'lm',
            se = FALSE, alpha = 0.5, size = 1.5) +
  theme_bw() +
  theme(legend.position = 'none') +
  ylim(c(0, 1))
  #ylim(c(0, 50))


out <- data.frame(species = rep(NA, length(nsp2)),
           hl_sl = NA,
           hl_pv = NA,
           fh_sl = NA,
           fh_pv = NA,
           hl_bs_sl = NA,
           hl_bs_pv = NA,
           fh_bs_sl = NA,
           fh_bs_pv = NA,
           rng_lat = NA,
           N_hl = NA,
           N_fh = NA)
for (i in 1:length(nsp2))
{
  #i <- 1
  tt <- dplyr::filter(td2, species == nsp2[i])
  hl_tfit <- summary(lm(hl ~ lat, data = tt))
  fh_tfit <- summary(lm(fh ~ lat, data = tt))
  
  # hl_bs_tfit <- summary(lm(bs ~ hl, data = tt))
  # fh_bs_tfit <- summary(lm(bs ~ fh, data = tt))
  
  out$species[i] <- nsp2[i]
  out$hl_sl[i] <- hl_tfit$coefficients[2,1]
  out$hl_pv[i] <- hl_tfit$coefficients[2,4]
  out$fh_sl[i] <- fh_tfit$coefficients[2,1]
  out$fh_pv[i] <- fh_tfit$coefficients[2,4]
  
  out$hl_bs_sl[i] <- hl_bs_tfit$coefficients[2,1]
  out$hl_bs_pv[i] <- hl_bs_tfit$coefficients[2,4]
  out$fh_bs_sl[i] <- fh_bs_tfit$coefficients[2,1]
  out$fh_bs_pv[i] <- fh_bs_tfit$coefficients[2,4]
  
  out$N_hl[i] <- NROW(dplyr::filter(tt, !is.na(hl)))
  out$N_fh[i] <- NROW(dplyr::filter(tt, !is.na(fh)))
  out$rng_lat[i] <- range(tt$lat)[2] - range(tt$lat)[1]
}

plot(out$fh_sl, out$hl_sl)
hist(out$hl_sl)
hist(dplyr::filter(out, hl_pv < 0.05)$hl_sl, main = 'hatch - lay ~ lat')
hist(out$fh_sl)
hist(dplyr::filter(out, fh_pv < 0.05)$fh_sl, main = 'fledge - hatch ~ lat')
