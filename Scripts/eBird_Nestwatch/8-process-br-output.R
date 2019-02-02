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
#dir <- '/home/CAM/cyoungflesh/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

halfmax_breeding_date <- '2019-01-30'
bc_query_date <- '2019-01-30'
NW_date <- '2019-01-28'
MAPS_date <- '2019-01-31'
eBird_query_dir <- 'eBird_query_2018-10-15'


# Load packages -----------------------------------------------------------

library(dplyr)


# Set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)
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
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', eBird_query_dir))
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
  #i <- 1
  
  #readin halfmax data
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/halfmax_breeding_', halfmax_breeding_date))
  temp_halfmax <- readRDS(paste0('halfmax_df_breeding_', species_list[i], '.rds'))
  
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/breeding_cat_query_', bc_query_date))
  temp_bc <- readRDS(paste0('ebird_NA_breeding_cat_', species_list[i], '.rds'))
  
  cells <- sort(unique(temp_halfmax$cell))
  years <- sort(unique(temp_halfmax$year))
  nyr <- length(years)
  ncell <- length(cells)
  
  #loop through years
  for (j in 1:nyr)
  {
    #j <- 1
    
    for (k in 1:ncell)
    {
      #k <- 1
      print(paste0('species: ', species_list[i], ', ',
                   'year: ', years[j], ', ',
                   'cell: ', cells[k]))
      
      #####################
      #ebird breeding codes
      #####################
      
      cysdata <- dplyr::filter(temp_bc, 
                               year == years[j] & 
                               cell == cells[k])
      
      tt_halfmax <- filter(temp_halfmax, 
                           year == years[j], 
                           cell == cells[k])
      
      #write species/cell/year data
      m_breeding_df$SPECIES[counter] <- species_list[i]
      m_breeding_df$CELL[counter] <- cells[k]
      m_breeding_df$YEAR[counter] <- years[j]
        
      #number of surveys where breeding was detected (confirmed or probable)
      m_breeding_df$EB_n1[counter] <- tt_halfmax$n1
      #number of surveys where breeding was not detected (bird not seen breeding or not seen)
      m_breeding_df$EB_n0[counter] <- tt_halfmax$n0
      #number of detections that came before jday 60
      m_breeding_df$EB_n1W[counter] <- tt_halfmax$n1W
      #number of unique days with detections
      m_breeding_df$EB_njd1[counter] <- tt_halfmax$njd1
      #number of unique days with non-detection
      m_breeding_df$EB_njd0[counter] <- tt_halfmax$njd0
      #number of unique days of non-detections before first detection
      m_breeding_df$EB_njd0i[counter] <- tt_halfmax$njd0i
      #number of non-detections before first detection
      m_breeding_df$EB_n0i[counter] <- tt_halfmax$n0i
      
      m_breeding_df$EB_min_neff[counter] <- tt_halfmax$min_neff
      m_breeding_df$EB_max_Rhat[counter] <- tt_halfmax$max_Rhat
      
      iter_ind <- grep('iter', colnames(tt_halfmax))
      halfmax_posterior <- as.numeric(tt_halfmax[,iter_ind])
      
      #determine how many estimates are 1 and not 1 (estimates of 1 are bogus)
      m_breeding_df$EB_nphen_bad[counter] <- sum(halfmax_posterior == 1)

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
        m_breeding_df$NW_mean_cid[counter] <- t_NW$MEAN_FIRST_LAY
        m_breeding_df$NW_sd_cid[counter] <- t_NW$SD_FIRST_LAY
        m_breeding_df$NW_num_obs[counter] <- t_NW$NUM_OBS
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
            m_breeding_df$MAPS_midpoint[counter] <- t_MAPS$midpoint
            m_breeding_df$MAPS_l_bounds[counter] <- t_MAPS$l_bounds
            m_breeding_df$MAPS_u_bounds[counter] <- t_MAPS$u_bounds
            m_breeding_df$MAPS_n_stations[counter] <- t_MAPS$n_stations
          }
        }
      }
      
      counter <- counter + 1
    } # k -cell
  } # j - year
} # i - species


#remove NA padding (find first year row to have NA and subtract one from that index)
fin_ind <- min(which(is.na(m_breeding_df$YEAR)))
m_breeding_df2 <- m_breeding_df[1:(fin_ind - 1),]


#check diags - look okay
ct <- m_breeding_df2[which(!is.na(m_breeding_df$EB_HM_mean)),]
min(ct$EB_min_neff)
max(ct$EB_max_Rhat)



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

sum(m_breeding_df2$d_avail == '111')
sum(m_breeding_df2$d_avail == '110')
sum(m_breeding_df2$d_avail == '101')
sum(m_breeding_df2$d_avail == '100')


#ALL YEARS/CELLS MODELED IN IAR MODEL
#write to rds
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
saveRDS(m_breeding_df2, paste0('temp_breeding_master_', Sys.Date(), '.rds'))









#m_breeding_df2 <- readRDS('temp_breeding_master_2019-01-30.rds')

#data where ebird is present
tt <- m_breeding_df2[which(m_breeding_df2$d_avail == '100' |
                       m_breeding_df2$d_avail == '110' |
                       m_breeding_df2$d_avail == '111' |
                       m_breeding_df2$d_avail == '101'),]

tt2 <- dplyr::select(tt, SPECIES, CELL, YEAR, EB_HM_mean, EB_HM_sd, 
                     NW_mean_cid, NW_sd_cid, NW_num_obs,
                     MAPS_midpoint, MAPS_n_stations, d_avail)

tt2$MAPS_midpoint[which(tt2$MAPS_midpoint == 0)] <- NA
tt2[which(tt2$NW_num_obs > 5),]
tt2[which(tt2$MAPS_midpoint > 180),]
tt2[which(abs(tt2$MAPS_midpoint - tt2$EB_HM_mean) > 20),]

#quick plots - so much uncertainty, hard to say anything from this
plot(tt2$NW_mean_cid, tt2$EB_HM_mean, pch = 19, col = rgb(0,0,0,0.5))
plot(tt2$MAPS_midpoint, tt2$EB_HM_mean, pch = 19, col = rgb(0,0,0,0.5))

tt2[which(tt2$MAPS_midpoint == 0),]

aa <- readRDS('Processed/arrival_master_2019-01-16.rds')
head(tt2)
colnames(tt2)[1:3] <- c('species', 'cell', 'year')
m2 <- dplyr::left_join(aa, tt2, by = c('species', 'cell', 'year'))

head(m2)
plot(m2$mean_post_IAR, m2$EB_HM_mean, pch = 19, col = rgb(0,0,0,0.5),
     xlim = c(50, 220), ylim = c(50, 220))
abline(a = 0, b = 1, lty = 2)
summary(lm(m2$EB_HM_mean ~ m2$mean_post_IAR))

#vvv OLD OLD vvv



#add 'meets criteria' column
diagnostics_frame$MODEL <- NA
diagnostics_frame$shp_fname <- NA


# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- df_out[with(df_out, order(species, year, cell)),]

#create scaling factor column
df_master$scaling_factor <- NA



# write to RDS --------------------------------------------------

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_', Sys.Date())

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(df_master, paste0('IAR_input-', Sys.Date(), '.rds'))


# explore data ------------------------------------------------------------

aggregate(n_cells ~ species, data = yrs_frame, FUN = max)
aggregate(n_cells ~ species, data = yrs_frame, FUN = mean)


