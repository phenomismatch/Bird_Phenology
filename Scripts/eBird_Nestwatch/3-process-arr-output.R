######################
# 3 - proces arrival cubic output
#
# Aggeregate posterior info and diagnostic info from 2-logit-cubic.R to be used in IAR model
# Determine which cells should be used in IAR model (cells that overlap breeding and migratory ranges, but do not overlap resident or non-breeding ranges) 
#
# Parts formerly in ICAR_parallel.R
######################


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'

#Xanadu
#dir <- '/UCHC/LABS/Tingley/phenomismatch/'



# db/hm query dir ------------------------------------------------------------

hm_dir <- 'halfmax_species_2019-02-02'
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
nsp <- length(species_list)


#DATA ONLY VALID THROUGH 2017 (2018 data only goes to ~ jday 60 as of 2018-10-15 query)
years <- 2002:2017
nyr <- length(years)


# get all potentially relevant cells --------------------------------------

#generate points over relevant coordinates
set.seed(1)

N <- 100000
u     <- runif(N)
v     <- runif(N)
theta <- 2*pi*u      * 180/pi
phi   <- acos(2*v-1) * 180/pi
lon   <- theta-180
lat   <- phi-90
sdf    <- data.frame(lat = lat, lon = lon)

sdf2 <- dplyr::filter(sdf , lat > 26 & lon > -100 & lon < -50)


#construct hex cells using points in US
hexgrid6 <- dggridR::dgconstruct(res = 6)
tcells <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                  in_lon_deg = sdf2$lon, in_lat_deg = sdf2$lat)[[1]]
cell_grid <- dggridR::dgcellstogrid(hexgrid6, tcells)


#plot of cells in relevant region
# ggplot() +
#   geom_path(data = usamap,
#             aes(x = x, y = y), color = 'black') +
#   geom_path(data = canadamap,
#             aes(x = x, y = y), color = 'black') +
#   geom_path(data = mexicomap,
#             aes(x = x, y = y), color = 'black') +
#   coord_map("ortho", orientation = c(35, -80, 0),
#             xlim = c(-100, -50), ylim = c(25, 90)) +
#   geom_path(data = cell_grid, aes(x = long, y = lat, group = group),
#             alpha = 0.5, color = 'red') +
#   xlab('Longitude') +
#   ylab('Latitude') +
#   theme_bw()


# Proces logit cubic results -----------------------------------------------------------------

#get number of cell/years
cell_years <- 0
for (i in 1:nsp)
{
  #i <- 113
  
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
  #i <- 113
  
  #import halfmax estimates and diagnostics from logit cubic model
  setwd(paste0(dir, 'Bird_Phenology/Data/Processed/', hm_dir))
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
                                    tree_depth = na_reps)
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

#remove rows that were unfilled
to.rm <- which(is.na(diagnostics_frame$n1))
diagnostics_frame2 <- diagnostics_frame[-to.rm,]



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
    t_sp <- dplyr::filter(diagnostics_frame2, species == species_list[i])
    
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
run_time <- (proc.time()[3] - tt[3]) / 60



# order -------------------------------------------------------------------

#order diagnostics frame by species, year, and cell #
df_master <- df_out[with(df_out, order(species, year, cell)),]


# write to RDS --------------------------------------------------

IAR_dir_path <- paste0(dir, 'Bird_phenology/Data/Processed/IAR_input_', hm_date)

dir.create(IAR_dir_path)
setwd(IAR_dir_path)

saveRDS(df_master, paste0('IAR_input-', hm_date, '.rds'))





# vvv HERE DOWN vvv





# create list of species to run through IAR model -------------------------

species_tm <- aggregate(MODEL ~ species, data = df_master, FUN = function(x) sum(x, na.rm = TRUE))$species

setwd(paste0(dir, 'Bird_Phenology/Data/'))
write.table(species_tm, file = 'IAR_species_list.txt', row.names = FALSE, col.names = FALSE)




# create dfs that show # cells with data in each year/species, and # years with data in each cell/species -----------------

#create vector of year names
yrs_vec <- c()
for (i in min(years):max(years))
{
  yrs_vec <- c(yrs_vec, paste0('yr_', i))
}

#create vector of cell names
cell_vec <- c()
for (i in 1:ncell)
{
  cell_vec <- c(cell_vec, paste0('cell_', cells[i]))
}


#create df with species/cells/n_yrs per sp_cell
cells_frame <- data.frame(species = rep(species_list, each = length(cells)), 
                          cell = rep(cells, length(species_list)), 
                          n_yrs = NA)

yrs_frame <- data.frame(species = rep(species_list, each = length(years)), 
                        year = rep(years, length(species_list)), 
                        n_cells = NA)


#add years and cell cols
cells_frame[yrs_vec] <- NA
yrs_frame[cell_vec] <- NA

#fill cells_frame and yrs_frame
for (i in 1:nsp)
{
  #i <- 80
  
  for (k in 1:ncell)
  {
    #k <- 30
    t_cell <- dplyr::filter(diagnostics_frame, species == species_list[i], cell == cells[k])
    
    #get index for species/cell row
    idx <- which(cells_frame$species == species_list[i] & 
                   cells_frame$cell == cells[k])
    
    #insert number of yrs with data
    yrs_d <- t_cell$year[which(!is.na(t_cell$HM_mean))]
    cells_frame[idx,'n_yrs'] <- length(yrs_d)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    temp_yrs <- paste0('yr_', yrs_d)
    cells_frame[idx, which(colnames(cells_frame) %in% temp_yrs)] <- TRUE
  }
  
  for (j in min(years):max(years))
  {
    #j <- 2015
    t_yr <- dplyr::filter(diagnostics_frame, species == species_list[i], year == j)
    
    #get index for species/year row
    idx_yr <- which(yrs_frame$species == species_list[i] & 
                      yrs_frame$year == j)
    
    #insert number of cells with data
    cells_d <- t_yr$cell[which(!is.na(t_yr$HM_mean))]
    yrs_frame[idx_yr,'n_cells'] <- length(cells_d)
    
    #find out which years have data and insert TRUE into df - leave NA otherwise
    temp_cells <- paste0('cell_', cells_d)
    yrs_frame[idx_yr, which(colnames(yrs_frame) %in% temp_cells)] <- TRUE
  }
}

#write to RDS
setwd(IAR_dir_path)
saveRDS(cells_frame, paste0('cells_frame-', Sys.Date(), '.rds'))
saveRDS(yrs_frame, paste0('yrs_frame-', Sys.Date(), '.rds'))



# explore data ------------------------------------------------------------

# aggregate(n_cells ~ species, data = yrs_frame, FUN = max)
# aggregate(n_cells ~ species, data = yrs_frame, FUN = mean)


