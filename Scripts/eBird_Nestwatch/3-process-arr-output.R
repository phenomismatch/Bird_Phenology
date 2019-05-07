######################
# 3 - proces arrival cubic output
#
# Aggeregate posterior info and diagnostic info from 2-logit-cubic.R to be used in IAR model
# Determine which cells should be used in IAR model (cells that overlap breeding and migratory ranges, but do not overlap resident or non-breeding ranges) - 18 hour run time?
#
# Parts formerly in ICAR_parallel.R
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

hm_dir <- 'halfmax_species_2019-03-29'
hm_date <- substr(hm_dir, start = 17, stop = 26)


# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)
library(sp)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]

#temporary filter out these species
# rm_sp <- which(species_list %in% c('Empidonax_virescens', 
#                                    'Empidonax_alnorum', 'Setophaga_americana'))
# species_list <- species_list[-rm_sp]
nsp <- length(species_list)

years <- 2002:2018
nyr <- length(years)



# Proces halfmax results -----------------------------------------------------------------

#get number of cell/years
cell_years <- 0
for (i in 1:nsp)
{
  #i <- 110
  
  #import halfmax estimates and diagnostics from logit cubic model
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
  temp_halfmax <- readRDS(paste0('halfmax_df_arrival_', species_list[i], '.rds'))

  if (!is.na(temp_halfmax$year[1]))
  {
    cell_years <- cell_years + NROW(temp_halfmax)
  }
}

counter <- 1
for (i in 1:nsp)
{
  #i <- 1
  
  #import halfmax estimates and diagnostics from logit cubic model
  #setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
  temp_halfmax <- readRDS(paste0('halfmax_df_arrival_', species_list[i], '.rds'))
  
  if (i == 1)
  {
    na_reps <- rep(NA, cell_years)
    
    diagnostics_frame <- data.frame(species = na_reps,
                                    year = na_reps,
                                    cell = na_reps,
                                    HM_mean = na_reps,
                                    HM_sd = na_reps,
                                    n1 = na_reps,
                                    n1W = na_reps,
                                    n0 = na_reps,
                                    n0i = na_reps,
                                    njd1 = na_reps,
                                    njd0 = na_reps,
                                    njd0i = na_reps,
                                    max_Rhat = na_reps,
                                    min_neff = na_reps,
                                    nphen_bad = na_reps,
                                    delta = na_reps,
                                    tree_depth = na_reps,
                                    num_diverge = na_reps,
                                    num_tree = na_reps,
                                    num_BFMI = na_reps)
  }
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 1
    tt_halfmax1 <- dplyr::filter(temp_halfmax, 
                                 year == years[j])
    
    cells <- unique(tt_halfmax1$cell)
    ncell <- length(cells)
    
    if (ncell > 0)
    {
      for (k in 1:ncell)
      {
        #k <- 1
        
        diagnostics_frame$species[counter] <- species_list[i]
        diagnostics_frame$year[counter] <- years[j]
        diagnostics_frame$cell[counter] <- cells[k]
        
        #get model fits
        tt_halfmax2 <- dplyr::filter(tt_halfmax1, 
                             cell == cells[k])
        
        print(paste0('species: ', species_list[i], ', ',
                     'year: ', years[j], ', ', 
                     'cell: ', cells[k]))

        #number of surveys where species was detected
        diagnostics_frame$n1[counter] <- tt_halfmax2$n1
        #number of surveys where species was not detected
        diagnostics_frame$n0[counter] <- tt_halfmax2$n0
        #number of detections that came before jday 60
        diagnostics_frame$n1W[counter] <- tt_halfmax2$n1W
        #number of non-detections before first detection
        diagnostics_frame$n0i[counter] <- tt_halfmax2$n0i
        #number of unique days with detections
        diagnostics_frame$njd1[counter] <- tt_halfmax2$njd1
        #number of unique days with non-detection
        diagnostics_frame$njd0[counter] <- tt_halfmax2$njd0
        #number of unique days of non-detections before first detection
        diagnostics_frame$njd0i[counter] <- tt_halfmax2$njd0i
          
        diagnostics_frame$min_neff[counter] <- tt_halfmax2$min_neff
        diagnostics_frame$max_Rhat[counter] <- tt_halfmax2$max_Rhat
        
        diagnostics_frame$delta[counter] <- tt_halfmax2$delta
        diagnostics_frame$tree_depth[counter] <- tt_halfmax2$tree_depth
        diagnostics_frame$num_diverge[counter] <- tt_halfmax2$num_diverge
        diagnostics_frame$num_tree[counter] <- tt_halfmax2$num_tree
        diagnostics_frame$num_BFMI[counter] <- tt_halfmax2$num_BFMI
        
        iter_ind <- grep('iter', colnames(tt_halfmax2))
        halfmax_posterior <- as.numeric(tt_halfmax2[,iter_ind])
        
        #determine how many estimates are 1 and not 1 (estimates of 1 are bogus)
        diagnostics_frame$nphen_bad[counter] <- sum(halfmax_posterior == 1)
        #halfmax_posterior2 <- halfmax_posterior[which(halfmax_posterior != 1)]
        
        #calculate posterior mean and sd
        if (sum(!is.na(halfmax_posterior)) > 0)
        {
          diagnostics_frame$HM_mean[counter] <- mean(halfmax_posterior)
          diagnostics_frame$HM_sd[counter] <- sd(halfmax_posterior)
        }
        
        counter <- counter + 1
      } # if loop for at least one cell - species without sufficient ranges have 0 cells
    } # k -cell
  } # j - year
} # i - species


#add 'meets criteria' column
diagnostics_frame$MODEL <- NA
diagnostics_frame$shp_fname <- NA

#no longer happens?
# #remove rows that were unfilled - happens due to species not having sufficient range
# to.rm <- which(is.na(diagnostics_frame$n1))
# diagnostics_frame2 <- diagnostics_frame[-to.rm,]
diagnostics_frame2 <- diagnostics_frame


### add NA for both HM_mean and HM_sd if any of the following conditions are met
#num_diverge > 0
#max_Rhat >= 1.1
#min_neff < 250
#num_BFMI > 0

to.NA <- which(diagnostics_frame2$num_diverge > 0 | 
                 diagnostics_frame2$max_Rhat >= 1.1 |
                 diagnostics_frame2$min_neff < 200 |
                 diagnostics_frame2$num_BFMI > 0 |
                 diagnostics_frame2$HM_sd > 20)

# #1.4% of cells are bad
# length(to.NA)/sum(!is.na(diagnostics_frame2$HM_mean))
# diagnostics_frame2[to.NA,c('species', 'cell', 'year', 
#                            'HM_mean', 'HM_sd', 'min_neff', 'num_diverge')]

if (length(to.NA) > 0)
{
  diagnostics_frame2[to.NA,'HM_mean'] <- NA
  diagnostics_frame2[to.NA,'HM_sd'] <- NA
}



# Filter data based on criteria -----------------------------------------------------------

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 'NC' cells in all three years from 2016-2018
#     Species-years with at least 'NC' cells for those species

NC <- 3

'%ni%' <- Negate('%in%')

#reference key for species synonyms
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')

#change dir to shp files
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))

df_out <- data.frame()
#which species/years meet criteria for model
for (i in 1:length(species_list))
{
    #i <- 1
    
    print(i)
    #filter by species
    t_sp <- dplyr::filter(diagnostics_frame2, species == species_list[i])
    
    if (NROW(t_sp) > 0)
    {
      #match species name to shp file name
      g_ind <- grep(species_list[i], sp_key$file_names_2016)
    
      #check for synonyms if there are no matches
      if (length(g_ind) == 0)
      {
        g_ind2 <- grep(species_list[i], sp_key$BL_Checklist_name)
      } else {
        g_ind2 <- g_ind
      }
    
      #get filename to put in df
      fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
      t_sp$shp_fname <- fname
    
      #number of cells with good data in each year from 2016-2018
      nobs_yr <- c()
      for (j in 2016:2018)
      {
        #j <- 2018
        ty_sp3 <- dplyr::filter(t_sp, year == j)
        ind <- which(!is.na(ty_sp3$HM_mean))
        nobs_yr <- c(nobs_yr, length(ind))
        #ty_sp[ind,]
      }

      #if all three years have greater than or = to 'NC' cells of data, figure 
      #...out which years have at least 'NC' cells
      yrs_kp <- c()
      if (sum(nobs_yr >= NC) == 3)
      {
        #see which years have more than 3 cells of data
        nobs_yr2 <- c()
        for (j in min(years):max(years))
        {
          #j <- 2012
          ty2_sp3 <- dplyr::filter(t_sp, year == j)
          ind2 <- which(!is.na(ty2_sp3$HM_mean))
          nobs_yr2 <- c(nobs_yr2, length(ind2))
        }
      
        #years to keep (more than three cells of data)
        yrs_kp <- years[which(nobs_yr2 >= NC)]
      }
    
      if (length(yrs_kp) > 0)
      {
        t_sp[which(t_sp$year %in% yrs_kp),]$MODEL <- TRUE
      }
    
      df_out <- rbind(df_out, t_sp)
    }
}
run_time <- (proc.time()[3] - tt[3]) / 60



# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- df_out[with(df_out, order(species, year, cell)),]


# write to RDS --------------------------------------------------

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_input_', hm_date)

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(df_master, paste0('IAR_input-', hm_date, '.rds'))



# create list of species to run through IAR model -------------------------

species_tm <- aggregate(MODEL ~ species, data = df_master, FUN = function(x) sum(x, na.rm = TRUE))$species

setwd(paste0(dir, 'Bird_Phenology/Data/'))
write.table(species_tm, file = paste0('IAR_species_list.txt'), row.names = FALSE, col.names = FALSE)




# # create dfs that show # cells with data in each year/species, and # years with data in each cell/species -----------------
# 
# #create df with species/cells/n_yrs per sp_cell
# cells_frame <- data.frame(species = rep(NA, cell_years), 
#                           cell = rep(NA, cell_years), 
#                           n_yrs = rep(NA, cell_years))
# 
# yrs_frame <- data.frame(species = rep(NA, cell_years), 
#                         year = rep(NA, cell_years), 
#                         n_cells = rep(NA, cell_years))
# 
# counter_cell <- 1
# counter_year <- 1
# #fill cells_frame and yrs_frame
# for (i in 1:nsp)
# {
#   #i <- 101
#   print(i)
#   
#   tspf <- dplyr::filter(df_master, species == species_list[i])
#   tcells <- unique(tspf$cell)
#   tyears <- unique(tspf$year)
#   
#   if (NROW(tspf) > 0)
#   {
#     for (k in 1:length(tcells))
#     {
#       #k <- 1
#       t_cell <- dplyr::filter(tspf, cell == tcells[k])
#     
#       cells_frame[counter_cell, 'species'] <- species_list[i]
#       cells_frame[counter_cell, 'cell'] <- tcells[k]
# 
#       #insert number of yrs with data
#       yrs_d <- t_cell$year[which(!is.na(t_cell$HM_mean))]
#       cells_frame[counter_cell,'n_yrs'] <- length(yrs_d)
#       counter_cell <- counter_cell + 1
#     }
#   
#     for (j in 1:length(tyears))
#     {
#       #j <- 1
#       t_year <- dplyr::filter(tspf, year == tyears[j])
#     
#       yrs_frame[counter_year, 'species'] <- species_list[i]
#       yrs_frame[counter_year, 'year'] <- tyears[j]
#     
#       #insert number of yrs with data
#       yrs_d <- t_cell$year[which(!is.na(t_year$HM_mean))]
#       yrs_frame[counter_year,'n_cells'] <- length(yrs_d)
#       counter_year <- counter_year + 1
#     }
#   }
# }
# 
# 
# #remove NA from end of df
# to.rm.cell <- min(which(is.na(cells_frame$species))):NROW(cells_frame)
# to.rm.year <- min(which(is.na(yrs_frame$species))):NROW(yrs_frame)
# 
# cells_frame2 <- cells_frame[-to.rm.cell,]
# yrs_frame2 <- yrs_frame[-to.rm.year,]
# 
# #write to RDS
# setwd(IAR_dir_path)
# saveRDS(cells_frame2, paste0('cells_frame-', hm_date, '.rds'))
# saveRDS(yrs_frame2, paste0('yrs_frame-', hm_date, '.rds'))



# explore data ------------------------------------------------------------

# aggregate(n_cells ~ species, data = yrs_frame, FUN = max)
# aggregate(n_cells ~ species, data = yrs_frame, FUN = mean)


