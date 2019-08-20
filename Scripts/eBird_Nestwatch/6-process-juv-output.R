######################
# 6 - proces juv GAM output
#
# Aggeregate posterior info and diagnostic info from 5-halfmax-juvs.R to be used in arr ~ br model
#
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'


# db/juv query dir ------------------------------------------------------------

juv_dir <- 'halfmax_juvs_2019-08-20'
juv_date <- substr(juv_dir, start = 14, stop = 23)

arr_br_dir <- 'arr_br_2019-08-19'


# runtime -----------------------------------------------------------------

tt <- proc.time()



# Load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)
library(sp)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)
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
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', juv_dir))
  temp_halfmax <- readRDS(paste0('halfmax_juvs_', species_list[i], '.rds'))

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
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', juv_dir))
  temp_halfmax <- readRDS(paste0('halfmax_juvs_', species_list[i], '.rds'))
  
  if (i == 1)
  {
    na_reps <- rep(NA, cell_years)
    
    diagnostics_frame <- data.frame(species = na_reps,
                                    year = na_reps,
                                    cell = na_reps,
                                    juv_mean = na_reps,
                                    juv_sd = na_reps,
                                    n1 = na_reps,
                                    #n1W = na_reps,
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
        # #number of detections that came before jday 60
        # diagnostics_frame$n1W[counter] <- tt_halfmax2$n1W
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
          diagnostics_frame$juv_mean[counter] <- mean(halfmax_posterior)
          diagnostics_frame$juv_sd[counter] <- sd(halfmax_posterior)
        }
        
        counter <- counter + 1
      } # if loop for at least one cell - species without sufficient ranges have 0 cells
    } # k -cell
  } # j - year
} # i - species



### add NA for both juv_mean and juv_sd if any of the following conditions are met

to.NA <- which(diagnostics_frame$num_diverge > 0 | 
                 diagnostics_frame$max_Rhat >= 1.1 |
                 diagnostics_frame$min_neff < 200 |
                 diagnostics_frame$num_BFMI > 0 |
                 diagnostics_frame$juv_sd > 10)

# #1.4% of cells are bad
# length(to.NA)/sum(!is.na(diagnostics_frame2$juv_mean))
# diagnostics_frame2[to.NA,c('species', 'cell', 'year',
#                            'juv_mean', 'juv_sd', 'min_neff', 'num_diverge')]

if (length(to.NA) > 0)
{
  diagnostics_frame2[to.NA,'juv_mean'] <- NA
  diagnostics_frame2[to.NA,'juv_sd'] <- NA
}




# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- diagnostics_frame[with(diagnostics_frame, order(species, year, cell)),]


# write to RDS --------------------------------------------------

ab_dir_path <- paste0(dir, 'Bird_Phenology/Data/Processed/', arr_br_dir)

dir.create(ab_dir_path)
setwd(ab_dir_path)

saveRDS(df_master, paste0('juv-output-', juv_date, '.rds'))
