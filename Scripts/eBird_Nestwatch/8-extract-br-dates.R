######################
# 5 - Extract arrival dates from IAR output
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

saveRDS(out, file = paste0('master_arrival_', IAR_out_date, '.rds'))



