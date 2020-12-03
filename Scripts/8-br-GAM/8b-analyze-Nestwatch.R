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

# #LAY
# #dev.off()
# par(mfrow = c(3,3))
# for (i in 1:length(usp2))
# {
#   #i <- 1
#   temp <- dplyr::filter(nw_data2, species == usp2[i])
#   if (sum(!is.na(temp$lay)) > 0)
#   {
#     hist(temp$lay, main = paste0('LAY - ', usp2[i]))
#   }
# }
# 
# for (i in 1:length(usp2))
# {
#   #i <- 1
#   temp <- dplyr::filter(nw_data, species == usp2[i])
#   if (sum(!is.na(temp$hatch)) > 0)
#   {
#     hist(temp$hatch, main = paste0('HATCH - ', usp2[i]))
#   }
# }
# 
# for (i in 1:length(usp2))
# {
#   #i <- 39
#   temp <- dplyr::filter(nw_data, species == usp2[i])
#   if (sum(!is.na(temp$fledge)) > 0)
#   {
#     hist(temp$fledge, main = paste0('FLEDGE - ', usp2[i]))
#   }
# }

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


# pheno windows -----------------------------------------------------------------

#on a cell/year basis, bc of compression between arr and lay over lat

# PERIODS (ideal)
#================
# impute missing periods with phylogenetic relationship
# ARRIVAL: IAR arr model
# FLEDGE: IAR bj model
# LAY: FLEDGE - median from NW
# HATCH: FLEDGE - median from NW

# PERIODS (current):
#================
# impute missing periods with cross-species median
# ARRIVAL: IAR arr model
# FLEDGE: LAY + NW interval
# LAY: ARRIVAL + 30
# HATCH: LAY + NW interval

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


#use this in 8-br-GAM model

#impute missing values with medians - alternative to phylogenetic imputation
# na_med_LH <- which(is.na(dd2$med_LH))
# na_med_LF <- which(is.na(dd2$med_LF))
# dd2$med_LH[na_med_LH] <- sp_med_LH
# dd2$med_LF[na_med_LF] <- sp_med_LF



# get windows -------------------------------------------------------------

# #species in arr master data
# usp_arr <- unique(arr_data2$species)
# 
# #empty df
# pw_data <- data.frame(species = rep(NA, NROW(arr_data2)),
#                       year = NA,
#                       cell = NA,
#                       mig_cell = NA,
#                       breed_cell = NA,
#                       cell_lat = NA,
#                       cell_lng = NA,
#                       arr = NA,
#                       lay = NA,
#                       hatch = NA,
#                       fledge = NA)
#                       
# #fill df
# counter <- 1
# for (i in 1:length(usp_arr))
# {
#   #i <- 1
#   temp <- dplyr::filter(arr_data2, species == usp_arr[i])
#   nrt <- NROW(temp)
#   tdd <- dplyr::filter(dd4, species == usp_arr[i])
#   
#   pw_data$species[counter:(counter + nrt - 1)] <- temp$species
#   pw_data$year[counter:(counter + nrt - 1)] <- temp$year
#   pw_data$cell[counter:(counter + nrt - 1)] <- temp$cell
#   pw_data$mig_cell[counter:(counter + nrt - 1)] <- temp$mig_cell
#   pw_data$breed_cell[counter:(counter + nrt - 1)] <- temp$breed_cell
#   pw_data$cell_lat[counter:(counter + nrt - 1)] <- temp$cell_lat
#   pw_data$cell_lng[counter:(counter + nrt - 1)] <- temp$cell_lng
#   
#   #filter dd4
#   #what were we going to do about lay date?
#   #backcalculate from estimated fledge...
#   
#   #arr
#   t_arr <- round(temp$arr_IAR_mean, 0)
#   pw_data$arr[counter:(counter + nrt - 1)] <- t_arr
#   #lay
#   pw_data$lay[counter:(counter + nrt - 1)] <- t_arr + 30
#   #hatch
#   pw_data$hatch[counter:(counter + nrt - 1)] <- t_arr + 30 + dd4$
#   #fledge
#   pw_data$fledge[counter:(counter + nrt - 1)] <- t_arr + 30 + tdd$med_LF
#   
#   # windows
#   # #arr -> lay
#   # pw_data$arr_lay_start[counter:(counter + nrt - 1)] <- temp$arr_IAR_mean
#   # arr_lay_end <- temp$arr_IAR_mean + 29
#   # pw_data$arr_lay_end[counter:(counter + nrt - 1)] <- arr_lay_end
#   # #lay -> hatch
#   # pw_data$lay_hatch_start[counter:(counter + nrt - 1)] <- arr_lay_end + 1
#   # lay_hatch_end <- arr_lay_end + 1 + tdd$med_LH
#   # pw_data$lay_hatch_end[counter:(counter + nrt - 1)] <- lay_hatch_end
#   # #hatch -> fledge
#   # pw_data$hatch_fledge_start[counter:(counter + nrt - 1)] <- lay_hatch_end + 1
#   # hatch_fledge_end <- arr_lay_end + 1 + tdd$med_LF
#   # pw_data$hatch_fledge_end[counter:(counter + nrt - 1)] <- lay_hatch_end + 1 + hatch_fledge_end
# 
#   counter <- counter + nrt
# }
# 
# tail(pw_data)


# write to rds ------------------------------------------------------------

# setwd(paste0(dir, '/Data/Processed/'))
# #saveRDS(pw_data, paste0('bird-pheno-periods-', RUN_DATE, '.rds'))
# write.csv(pw_data, paste0('bird-pheno-periods-', RUN_DATE, '.csv'), row.names = FALSE)
# 
# round(median(pw_data$arr, na.rm = TRUE), 0)
# round(median(pw_data$lay, na.rm = TRUE), 0)
# round(median(pw_data$hatch, na.rm = TRUE), 0)
# round(median(pw_data$fledge, na.rm = TRUE), 0)
