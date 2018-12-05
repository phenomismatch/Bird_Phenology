######################
# 5 - Extract arrival dates
#
# Compiles dataframe with pre (2-logit_cubic.R) and post (4-IAR-model.R) IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
#dir <- '~/Google_Drive/R/'

#Xanadu
dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# other dir ---------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2018-11-12'
IAR_out_dir <- 'IAR_output_2018-11-12'



# Load packages -----------------------------------------------------------

library(MCMCvis)
library(dplyr)
library(dggridR)


# setwd -------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_in_dir))

IAR_in_date <- substr(IAR_in_dir, start = 11, stop = 20)
IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)



# Filter data by species/years ------------------------------------------------------

#read in master df
df_master <- readRDS(paste0('IAR_input-', IAR_in_date, '.rds'))

species <- as.character(read.table('../../IAR_species_list.txt')[,1])

#switch to out dir
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))



#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)
out <- data.frame()
for (i in 1:length(species))
{
  #i <- 96 #Vireo olivaceus
  
  #filter by species
  sp <- species[i]
  f_in <- dplyr::filter(df_master, species == sp)
  
  #define cells and years to be modeled
  t_cells <- unique(f_in$cell)
  t_years <- unique(f_in$year)

  #make hexgrid
  hexgrid6 <- dggridR::dgconstruct(res = 6)

  #get hexgrid cell centers
  cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)

  #read in IAR model output
  t_fit <- readRDS(paste0('IAR_stan_', sp, '-', IAR_out_date, '.rds'))

  #extract median and sd for IAR arrival dates
  med_fit <- round(MCMCpstr(t_fit, params = 'mu', func = median)[[1]], digits = 2)
  sd_fit <- round(MCMCpstr(t_fit, params = 'mu', func = sd)[[1]], digits = 2)
  
  #loop through years
  for (j in 1:length(t_years))
  {
    #j <- 1
    
    t_f_in <- dplyr::filter(f_in, year == t_years[j])
    
    t_full <- data.frame(t_f_in[,c('species','cell')], 
                       cell_lat = round(cellcenters$lat_deg, digits = 2), 
                       cell_lon = round(cellcenters$lon_deg, digits = 2),
                       t_f_in[,c('year', 'HM_mean', 'HM_sd')],
                       med_post_IAR = med_fit[,j], sd_post_IAR = sd_fit[,j])
    
    colnames(t_full)[6:7] <- c('med_pre_IAR', 'sd_pre_IAR')
    
    t_full$med_pre_IAR <- round(t_full$med_pre_IAR, digits = 2)
    t_full$sd_pre_IAR <- round(t_full$sd_pre_IAR, digits = 2)
    
    out <- rbind(out, t_full)
  }
}


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

#write.csv(out, paste0(species[i], '_arrival_', IAR_out_date, '.csv'), row.names = FALSE)
write.csv(out, paste0('master_arrival_', IAR_out_date, '.csv'), row.names = FALSE)









# vvvvvv OLD vvvvvv --------------------------------------------------------------



# This script produces an index of the phenology of a given year in a given cell, based on bird arrival
# (see NA_birdPhen1.R and NA_birdPhen2.R). For now, it just uses the raw (unsmoothed) estimated half-max 
# dates.  For each cell, the script determines the list of species with a half-max estimate in all years 
# from 2009-present and from 2004-present.  It then computes a phenological index for that cell by taking
# the inverse-variance-weighted average of the posterior means for those species. By focusing only on
# species that have valid phenological estimates in every year, we ensure that interannual fluctuations 
# in the index do not result from variation in which species are included in the index. Note, however,
# that differences between cells could indeed be due to differences in which species are included.

# In the future, when we have spatial models available for most species-years, we will be able to 
# predict the half-max even in species-cell-years where we couldn't fit it directly, and we will be
# able to include a lot more species in these metrics.

library(dggridR)
library(coda)
library(ggplot2)
library(viridis)
library(rstan)
'%ni%' <- Negate('%in%')
hexgrid6 <- dgconstruct(res=6) # Construct geospatial hexagonal grid

years <- 2002:2016
nyr <- length(years)

species_list <- c("Empidonax_virescens", "Myiarchus_crinitus", "Contopus_virens", "Vireo_olivaceus",
                  "Vireo_solitarius", "Vireo_gilvus", "Vireo_flavifrons", "Catharus_fuscescens",
                  "Dumetella_carolinensis", "Setophaga_dominica", "Limnothlypis_swainsonii",
                  "Setophaga_citrina", "Geothlypis_formosa", "Parkesia_motacilla", "Parkesia_noveboracensis",
                  "Mniotilta_varia", "Setophaga_americana", "Setophaga_ruticilla",
                  "Setophaga_virens", "Setophaga_caerulescens", "Protonotaria_citrea", "Setophaga_cerulea",
                  "Seiurus_aurocapilla", "Cardellina_canadensis", "Piranga_olivacea", "Piranga_rubra",
                  "Pheucticus_ludovicianus", "Icterus_galbula", "Empidonax_traillii", "Empidonax_alnorum", #30
                  "Empidonax_minimus", "Tyrannus_tyrannus", "Vireo_bellii", "Vireo_griseus", "Tachycineta_bicolor",
                  "Stelgidopteryx_serripennis","Hirundo_rustica","Riparia_riparia","Petrochelidon_pyrrhonota",
                  "Progne_subis", "Vermivora_cyanoptera","Vermivora_chrysoptera","Oreothlypis_ruficapilla", #43
                  "Setophaga_pensylvanica", "Setophaga_petechia", "Setophaga_discolor", "Geothlypis_philadelphia",
                  "Pooecetes_gramineus", "Ammodramus_nelsoni", "Passerina_cyanea", #50
                  "Passerina_caerulea","Spiza_americana","Icterus_spurius","Dolichonyx_oryzivorus","Contopus_cooperi", #55
                  "Empidonax_flaviventris","Regulus_satrapa","Regulus_calendula","Vireo_philadelphicus", #59
                  "Troglodytes_hiemalis","Catharus_guttatus","Catharus_ustulatus","Catharus_bicknelli", #63
                  "Catharus_minimus","Setophaga_fusca","Setophaga_striata","Setophaga_tigrina","Oreothlypis_peregrina", #68
                  "Setophaga_castanea","Setophaga_palmarum","Oreothlypis_celata","Cardellina_pusilla", #72
                  "Oporornis_agilis","Setophaga_magnolia","Setophaga_coronata","Passerella_iliaca",
                  "Melospiza_lincolnii","Spizelloides_arborea","Junco_hyemalis","Zonotrichia_leucophrys", #80
                  "Zonotrichia_albicollis","Euphagus_carolinus","Ictinia_mississippiensis",
                  "Elanoides_forficatus","Archilochus_colubris","Coccyzus_erythropthalmus","Coccyzus_americanus",
                  "Antrostomus_vociferus","Antrostomus_carolinensis","Sphyrapicus_varius","Sayornis_phoebe",
                  "Bombycilla_cedrorum","Cyanocitta_cristata","Turdus_migratorius","Polioptila_caerulea",
                  "Setophaga_pinus","Peucaea_aestivalis","Spinus_tristis",
                  "Geothlypis_trichas","Icteria_virens","Mimus_polyglottos","Lanius_ludovicianus",
                  "Melospiza_georgiana","Quiscalus_quiscula","Passerculus_sandwichensis","Ammodramus_caudacutus",
                  "Spizella_passerina","Pipilo_erythrophthalmus","Troglodytes_aedon","Sturnella_magna",
                  "Sialia_sialis","Cistothorus_palustris","Corvus_ossifragus","Agelaius_phoeniceus")
nsp <- length(species_list)



load("/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/birdPhenRaw.Rdata")
diagnostics_frame <- birdPhenRaw

diag_frame <- diagnostics_frame[which(diagnostics_frame$n1 > 29 & diagnostics_frame$n0 > 29 &
                                        diagnostics_frame$n1W < diagnostics_frame$n1/50 &
                                        diagnostics_frame$njd0i > 29 & diagnostics_frame$njd1 > 19), ]

winter_spcels <- unique(diagnostics_frame$spCel[which(diagnostics_frame$n1W > diagnostics_frame$n1/50)])

diag_frame <- diag_frame[which(diag_frame$spCel %ni% winter_spcels), ]

spcs <- unique(diag_frame$spCel)
nspc <- length(spcs)

cells <- unique(diag_frame$cell)
ncel <- length(cells)

for(i in 2002:2016){print(length(unique(diag_frame$spCel[which(diag_frame$year==i)])))}
hist(diag_frame$phen_mean)

bp2009 <- data.frame()
for(i in 1:nspc){
  print(i)
  B <- diag_frame[which(diag_frame$spCel==spcs[i]), ]
  YR <- B$year
  if(sum(c(2009:2016) %in% YR) == 8){
    bp2009 <- rbind(bp2009, B)
  }
}
bp2009 <- bp2009[which(bp2009$year > 2008),]
dim(bp2009)

bp2009$anom <- NA

for(i in 1:dim(bp2009)[1]){
  B <- bp2009[which(bp2009$spCel == bp2009$spCel[i]),]
  bp2009$anom[i] <- bp2009$phen_mean[i] - mean(B$phen_mean)
  bp2009$anomSE[i] <- sqrt(((7/8)*bp2009$phen_sd[i])^2 + sum(B$phen_sd[which(B$year != bp2009$year[i])]^2)/7)
}

hex6_birdAll2009 <- data.frame(cell=rep(cells,each=8), year=c(2009:2016), mean=NA, SE=NA)
for(i in 1:ncel){
  B <- bp2009[which(bp2009$cell == cells[i]), ]
  if(dim(B)[1] > 0){
    for(k in 1:8){
      B2 <- B[which(B$year == (2008+k)), ]
      mus <- vector()
      sds <- vector()
      counter <- 0
      for(j in 1:nsp){
        if(species_list[j] %in% B2$species){
          counter <- counter+1
          B3 <- B2[which(B2$species == species_list[j]), ]
          mus[counter] <- B3$anom[which(B3$species == species_list[j])]
          sds[counter] <- B3$anomSE[which(B3$species == species_list[j])]
        }
      }
      hex6_birdAll2009$mean[which(hex6_birdAll2009$cell == cells[i] & hex6_birdAll2009$year == (2008+k))] <- 
        weighted.mean(mus, 1/(sds^2))
      hex6_birdAll2009$SE[which(hex6_birdAll2009$cell == cells[i] & hex6_birdAll2009$year == (2008+k))] <-
        sqrt(1/sum(1/(sds^2)))
    }
  }
}

save(hex6_birdAll2009, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/hex6_birdAll2009.Rdata")
