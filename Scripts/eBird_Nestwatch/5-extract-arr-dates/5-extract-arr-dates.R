######################
# 5 - Extract arrival dates from IAR output
#
# Compiles dataframe with pre (2-halfmax-arr.R) and post (4-IAR-arr.R) IAR data for every year/cell
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# other dir ---------------------------------------------------------------

IAR_in_dir <- 'IAR_input_2019-02-02'
IAR_out_dir <- 'IAR_output_2019-02-13'



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
#setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', IAR_out_dir))
#setwd(paste0('~/Desktop/Bird_Phenology_Offline/Bird_Phenology/Data/Processed/', IAR_out_dir))


#combine pre and post IAR data for every year/cell that was modeled (including cell lat/lon)
out <- data.frame()
for (i in 1:length(species))
{
  #i <- 96 #Vireo olivaceus
  #i <- 24 #Empidonax virescens
  
  #filter by species
  sp <- species[i]
  
  #if that species RDS object exists in dir
  if (length(grep(paste0(sp, '-', IAR_out_date, '-iar-stan_output.rds'), list.files())) > 0)
  {
    f_in <- dplyr::filter(df_master, species == sp & MODEL == TRUE)
    
    #cells and years that were modeled
    t_cells <- unique(f_in$cell)
    t_years <- unique(f_in$year)
  
    #make hexgrid
    hexgrid6 <- dggridR::dgconstruct(res = 6)
  
    #get hexgrid cell centers
    cellcenters <- dggridR::dgSEQNUM_to_GEO(hexgrid6, t_cells)
  
    #read in IAR model output
    t_fit <- readRDS(paste0(sp, '-', IAR_out_date, '-iar-stan_output.rds'))
  
    #extract median and sd for IAR arrival dates
    mean_fit <- round(MCMCpstr(t_fit, params = 'y_true', func = mean)[[1]], digits = 2)
    med_fit <- round(MCMCpstr(t_fit, params = 'y_true', func = median)[[1]], digits = 2)
    sd_fit <- round(MCMCpstr(t_fit, params = 'y_true', func = sd)[[1]], digits = 2)
    
    #loop through years
    for (j in 1:length(t_years))
    {
      #j <- 1
      print(paste0('species: ', sp, ', ', 
                   'year: ', t_years[j]))
      
      t_f_in <- dplyr::filter(f_in, year == t_years[j])
      
      t_full <- data.frame(t_f_in[,c('species','cell')], 
                         cell_lat = round(cellcenters$lat_deg, digits = 2), 
                         cell_lon = round(cellcenters$lon_deg, digits = 2),
                         t_f_in[,c('year', 'HM_mean', 'HM_sd')],
                         mean_post_IAR = mean_fit[,j], sd_post_IAR = sd_fit[,j])
      
      colnames(t_full)[6:7] <- c('mean_pre_IAR', 'sd_pre_IAR')
      
      t_full$mean_pre_IAR <- round(t_full$mean_pre_IAR, digits = 2)
      t_full$sd_pre_IAR <- round(t_full$sd_pre_IAR, digits = 2)
     
      out <- rbind(out, t_full)
    } #end year loop
  } else {
    print(paste0('.rds file for ', sp, ' not found in directory'))
  }
} #end species loop


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

saveRDS(out, file = paste0('arrival_master_', IAR_out_date, '.rds'))

print('I completed!')

