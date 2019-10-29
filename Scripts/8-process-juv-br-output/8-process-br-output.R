######################
# 8 - Process breeding code logit cubic output
#
# Compiles dataframe with breeding data from Nestwatch (XXX), eBird breeding codes (7-logit-cubic-breeding.R), and MAPS observations (6-query-nesting-data.R) for every year/cell

#Run time: ~ 15 hours
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/labs/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

halfmax_breeding_date <- '2019-02-13'
NW_date <- '2019-01-28'
MAPS_date <- '2019-01-31'
eBird_query_date <- '2019-02-02'


# Load packages -----------------------------------------------------------

library(dplyr)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('IAR_species_list-2019-02-02.txt', stringsAsFactors = FALSE)
species_list <- species_list_i[,1]
nsp <- length(species_list)


#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)
years <- 2002:2017
nyr <- length(years)


# read in NW and MAPS data ------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))

NW_data <- readRDS(paste0('breeding_NW_', NW_date, '.rds'))
MAPS_data <- readRDS(paste0('breeding_MAPS_obs_', MAPS_date, '.rds'))


# create and fill df ------------------------------------------------------

#read in ebird surveys to get full list of cells
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/eBird_query_', eBird_query_date))
t_data <- readRDS(paste0('ebird_NA_phen_proc_', species_list[1], '.rds'))

na_reps <- rep(NA, (nsp*nyr*length(unique(t_data$cell))))
rm(t_data)


m_breeding_df <- data.frame(species = na_reps,
                            cell = na_reps,
                            year = na_reps,
                            EB_HM_mean = na_reps,
                            EB_HM_sd = na_reps,
                            EB_n1 = na_reps,
                            EB_n1W = na_reps,
                            EB_n0 = na_reps,
                            EB_n0i = na_reps,
                            EB_njd1 = na_reps,
                            EB_njd0 = na_reps,
                            EB_njd0i = na_reps,
                            EB_max_Rhat = na_reps,
                            EB_min_neff = na_reps,
                            EB_nphen_bad = na_reps,
                            EB_delta = na_reps,
                            EB_tree_depth = na_reps,
                            EB_num_diverge = na_reps,
                            EB_num_tree = na_reps,
                            EB_num_BFMI = na_reps,
                            NW_mean_cid = na_reps,
                            NW_sd_cid = na_reps,
                            NW_num_obs = na_reps,
                            MAPS_midpoint = na_reps,
                            MAPS_l_bounds = na_reps,
                            MAPS_u_bounds = na_reps,
                            MAPS_n_stations = na_reps)

counter <- 1
for (i in 1:nsp)
{
  #i <- 87
  
  #readin halfmax data
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_breeding_', halfmax_breeding_date))
  temp_halfmax <- readRDS(paste0('halfmax_df_breeding_', species_list[i], '.rds'))
  
  cells <- sort(unique(temp_halfmax$cell))
  ncell <- length(cells)
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 2
    tt_halfmax1 <- dplyr::filter(temp_halfmax, 
                                 year == years[j])
    
    cells <- unique(tt_halfmax1$cell)
    ncell <- length(cells)
    
    for (k in 1:ncell)
    {
      #k <- 75
      
      #write species/cell/year data
      m_breeding_df$species[counter] <- species_list[i]
      m_breeding_df$year[counter] <- years[j]
      m_breeding_df$cell[counter] <- cells[k]

      print(paste0('species: ', species_list[i], ', ',
                   'year: ', years[j], ', ',
                   'cell: ', cells[k]))
      
      #####################
      #ebird breeding codes
      #####################

      #get model fits
      tt_halfmax2 <- dplyr::filter(tt_halfmax1, 
                                   cell == cells[k])
      
      #number of surveys where breeding was detected (confirmed or probable)
      m_breeding_df$EB_n1[counter] <- tt_halfmax2$n1
      #number of surveys where breeding was not detected (bird not seen breeding or not seen)
      m_breeding_df$EB_n0[counter] <- tt_halfmax2$n0
      #number of detections that came before jday 60
      m_breeding_df$EB_n1W[counter] <- tt_halfmax2$n1W
      #number of non-detections before first detection
      m_breeding_df$EB_n0i[counter] <- tt_halfmax2$n0i
      #number of unique days with detections
      m_breeding_df$EB_njd1[counter] <- tt_halfmax2$njd1
      #number of unique days with non-detection
      m_breeding_df$EB_njd0[counter] <- tt_halfmax2$njd0
      #number of unique days of non-detections before first detection
      m_breeding_df$EB_njd0i[counter] <- tt_halfmax2$njd0i

      m_breeding_df$EB_min_neff[counter] <- tt_halfmax2$min_neff
      m_breeding_df$EB_max_Rhat[counter] <- tt_halfmax2$max_Rhat
      
      m_breeding_df$delta[counter] <- tt_halfmax2$delta
      m_breeding_df$tree_depth[counter] <- tt_halfmax2$tree_depth
      m_breeding_df$num_diverge[counter] <- tt_halfmax2$num_diverge
      m_breeding_df$num_tree[counter] <- tt_halfmax2$num_tree
      m_breeding_df$num_BFMI[counter] <- tt_halfmax2$num_BFMI
      
      iter_ind <- grep('iter', colnames(tt_halfmax2))
      halfmax_posterior <- as.numeric(tt_halfmax2[,iter_ind])
      
      
      #determine how many estimates are 1 and not 1 (estimates of 1 are bogus)
      m_breeding_df$EB_nphen_bad[counter] <- sum(halfmax_posterior == 1)
      #halfmax_posterior2 <- halfmax_posterior[which(halfmax_posterior != 1)]

      #calculate posterior mean and sd
      if (sum(!is.na(halfmax_posterior)) > 0)
      {
        m_breeding_df$EB_HM_mean[counter] <- mean(halfmax_posterior)
        m_breeding_df$EB_HM_sd[counter] <- sd(halfmax_posterior)
      }
      
      ########
      #NW data
      ########
      
      t_NW <- dplyr::filter(NW_data, 
                            SCI_NAME == species_list[i],
                            CELL == cells[k],
                            YEAR == years[j])
      
      if (NROW(t_NW) > 0)
      {
        m_breeding_df$NW_mean_cid[counter] <- as.numeric(t_NW$MEAN_FIRST_LAY)
        m_breeding_df$NW_sd_cid[counter] <- as.numeric(t_NW$SD_FIRST_LAY)
        m_breeding_df$NW_num_obs[counter] <- as.numeric(t_NW$NUM_OBS)
      }
      
      #########
      #MAPS obs
      #########

      t_MAPS <- dplyr::filter(MAPS_data, 
                              SCINAME == species_list[i], 
                              cell == cells[k],
                              YR == years[j])
      
      if (NROW(t_MAPS) > 0)
      {
        #julian day used if breeding was observed
        #values of 0 for midpoint indicate observed, but not observed breeding - NA indicated not observed
        if (!is.na(t_MAPS$midpoint))
        {
          if (t_MAPS$midpoint > 0)
          {
            m_breeding_df$MAPS_midpoint[counter] <- as.numeric(t_MAPS$midpoint)
            m_breeding_df$MAPS_l_bounds[counter] <- as.numeric(t_MAPS$l_bounds)
            m_breeding_df$MAPS_u_bounds[counter] <- as.numeric(t_MAPS$u_bounds)
            m_breeding_df$MAPS_n_stations[counter] <- as.numeric(t_MAPS$n_stations)
          }
        }
      }
      
      counter <- counter + 1
    } # k -cell
  } # j - year
} # i - species

#remove NA padding (find first year row to have NA and subtract one from that index)
to.rm <- which(is.na(m_breeding_df$EB_n1))
m_breeding_df2 <- m_breeding_df[-to.rm,]


### add NA for both EB_HM_mean and EB_HM_sd if any of the following conditions are met
#num_diverge > 0
#max_Rhat >= 1.1
#min_neff < 250
#num_BFMI > 0

to.NA <- which(m_breeding_df2$num_diverge > 0 | 
                 m_breeding_df2$max_Rhat >= 1.1 |
                 m_breeding_df2$min_neff < 200 |
                 m_breeding_df2$num_BFMI > 0)

if (length(to.NA) > 0)
{
  m_breeding_df2[to.NA,'EB_HM_mean'] <- NA
  m_breeding_df2[to.NA,'EB_HM_sd'] <- NA
}

# #Look at plots that have a number of nphen_bad...why is this hapenning
# range(df_out$nphen_bad, na.rm = TRUE)
# hist(df_out$nphen_bad)
# df_out[which(df_out$nphen_bad > 10),]
# 
# 
# #remove if sd is greater than 15? That should filter out higih nphen_bad values
# NROW(df_out[which(df_out$nphen_bad > 10),])
# NROW(df_out[which(df_out$nphen_bad > 0),])
# hist(df_out$nphen_bad[-(which(df_out$nphen_bad > 5))])
# 
# NROW(dplyr::filter(df_out, HM_sd > 15))
# hist(df_out$HM_sd, breaks = 20)

#dplyr::filter(df_out, species == 'Agelaius_phoeniceus', year == 2014)




#Add column with codes:
#Find which rows have data for EB, NW, and MAPS
#CODE = 000 [EB|NW|MAPS]
#where 1 means dataset is available
#for example: 001 = no EB, no NW, yes MAPS
#111 = yes EB, yes NW, yes MAPS

EB <- as.numeric(!is.na(m_breeding_df2$EB_HM_mean))
NW <- as.numeric(!is.na(m_breeding_df2$NW_mean_cid))
MAPS <- as.numeric(!is.na(m_breeding_df2$MAPS_midpoint))

m_breeding_df2$d_avail <- paste0(EB, NW, MAPS)

# sum(m_breeding_df2$d_avail == '111')
# sum(m_breeding_df2$d_avail == '110')
# sum(m_breeding_df2$d_avail == '101')
# sum(m_breeding_df2$d_avail == '100')


#add 'meets criteria' column
m_breeding_df2$MODEL <- NA
m_breeding_df2$shp_fname <- NA



# save RDS for later use --------------------------------------


saveRDS(m_breeding_df2, 'temp_br_processing.rds')




# Filter data based on criteria -----------------------------------------------------------

# Which species-years are worth modeling? At a bare minimum:
#     Species with at least 'NC' cells in all three years from 2015-2017
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
  t_sp <- dplyr::filter(m_breeding_df2, species == species_list[i])
  
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
    
    #number of cells with good data in each year from 2015-2017
    nobs_yr <- c()
    for (j in 2015:2017)
    {
      #j <- 2017
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

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_br_input_', halfmax_breeding_date)

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(df_master, paste0('IAR_br_input-', halfmax_breeding_date, '.rds'))



# create list of species to run through IAR model -------------------------

species_tm <- aggregate(MODEL ~ species, data = df_master, FUN = function(x) sum(x, na.rm = TRUE))$species

setwd(paste0(dir, 'Bird_Phenology/Data/'))
write.table(species_tm, file = paste0('IAR_br_species_list-', halfmax_breeding_date, '.txt'), 
            row.names = FALSE, col.names = FALSE)


