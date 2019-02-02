####################
# 1 - query ebird data
#
# Filters eBird data from DB
# Zero fills
# Creates directory of processed data (rds file for each species) and copy of this 
# ... script (for reproducability) in /Data/Processed/eBird_query_<DATE>
#
# Runtime: 8-9 hours on 7 core machine - could be sped if run on on Xanadu but need to add IP to DB whitelist
#
# Replaces 1-process-ebird-data.R, which processed data based on local copy of eBird reference dataset
####################


# #can access DB from command line with:
# psql "sslmode=disable dbname=sightings user=cyoungflesh hostaddr=35.221.16.125"



tt <- proc.time()


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(dggridR)
library(doParallel)
library(foreach)

# set wd ------------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list_i <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)

#remove underscore and coerce to vector
species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))

#combine species names into a single string with quotes
SL <- paste0("'", species_list_i2, "'", collapse = ", ")
#SL <- paste0("'Empidonax virescens'")


# access DB ---------------------------------------------------------------

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")



# create query dir and navigate there -------------------------------------------


query_dir_path <- paste0('Processed/eBird_query_', Sys.Date())

dir.create(query_dir_path)
setwd(query_dir_path)




# filter dataset notes ----------------------------------------------------------

#SQL FILTER BY:
# * only ebird data
# * species of interest (import from text file)
# * only complete checklists
# * duration minutes (6 min to 24 hours)
# * lat (> 26) and lon (-100 to -50)
# * year > 2001
# * day of year (< julian day 200)
# * time started before 16:00
# * radius < 100km


#FILTERED IN R:
# * only unique groups when not null


#EBIRD DB FIELDS FROM FROM RAFE
# place_json contains:
# COUNTRY_CODE
# STATE_CODE
# COUNTY_CODE
# IBA_CODE
# BCR_CODE
# USFWS_CODE
# ATLAS_BLOCK
# LOCALITY_ID
# LOCALITY_TYPE
# EFFORT_AREA_HA
# 
# event_json contains:
# SAMPLING_EVENT_IDENTIFIER
# EFFORT_AREA_HA
# APPROVED
# REVIEWED
# NUMBER_OBSERVERS
# ALL_SPECIES_REPORTED
# OBSERVATION_DATE
# GROUP_IDENTIFIER
# DURATION_MINUTES
# 
# counts_json contains:
# SCIENTIFIC_NAME
# GLOBAL_UNIQUE_IDENTIFIER
# LAST_EDITED_DATE
# TAXONOMIC_ORDER
# CATEGORY
# SUBSPECIES_SCIENTIFIC_NAME
# BREEDING_BIRD_ATLAS_CODE
# BREEDING_BIRD_ATLAS_CATEGORY
# AGE_SEX
# OBSERVER_ID
# HAS_MEDIA




# Query and filter - event_id ----------------------------------------------------

#*get info for each unique event
#*create species columns with NA
#*in a loop, query each species individually and fill obs with 1s, no obs with 0s
#*merge with hex cells
#*create rds objects for each species


#filter all unique event_ids that meet criteria - about 38 min to complete query

data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT DISTINCT ON (event_id) event_id, year, day, place_id, lat, lng, started, 
                                    radius,
                                    (event_json ->> 'ALL_SPECIES_REPORTED')::int AS all_species_reported,
                                    (event_json ->> 'DURATION_MINUTES')::int AS duration_minutes,
                                    count_json ->> 'OBSERVER_ID' AS observer_id,
                                    count_json ->> 'BREEDING_BIRD_ATLAS_CODE' AS BBA_code,
                                    (event_json ->> 'NUMBER_OBSERVERS')::int AS number_observers,
                                    event_json ->> 'GROUP_IDENTIFIER' AS group_identifier
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxons USING (taxon_id)
                                    WHERE dataset_id = 'ebird'
                                    AND year > 2001
                                    AND day < 200
                                    AND lng BETWEEN -100 AND -50
                                    AND lat > 26
                                    AND (sci_name IN (", SL,"))
                                    AND (event_json ->> 'ALL_SPECIES_REPORTED')::int = 1
                                    AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
                                    AND LEFT(started, 2)::int < 16
                                    AND RADIUS < 100000;
                                    "))



#only values with unique group identifiers (or no group identifier)
data2 <- data[!duplicated(data[,'group_identifier'], 
                          incomparables = NA),]

rm(data)



# add sjday, sjday^2, sjday^3 and shr ---------------------------------------------------

#calculate polynomial then center data
#julian day, julian day^2, and julian day^3
SJDAY  <- as.vector(data2$day)
SJDAY2 <- as.vector(data2$day^2)
SJDAY3 <- as.vector(data2$day^3)

data2$sjday <- SJDAY
data2$sjday2 <- SJDAY2
data2$sjday3 <- SJDAY3

#scaled effort hours
SHR <- as.vector(scale((data2$duration_minutes/60), scale = FALSE))

data2$shr <- SHR



# bin lat/lon to hex grid and add to data ---------------------------------------------------------

# Construct geospatial hexagonal grid

hexgrid6 <- dggridR::dgconstruct(res = 6) 
data2$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                        in_lon_deg = data2$lng, 
                                        in_lat_deg = data2$lat)[[1]]




# query individual species, zero fill, and create RDS objects ----------------------------------


#create species columns
data2[species_list_i[,1]] <- NA

nsp <- NROW(species_list_i)

#run in parallel with 7 logical cores
doParallel::registerDoParallel(cores = 7)

tt <- proc.time()
foreach::foreach(i = 1:nsp) %dopar%
{
  #i <- 89
  print(i)
  
  pg <- DBI::dbDriver("PostgreSQL")
  
  cxn <- DBI::dbConnect(pg, 
                        user = "cyoungflesh", 
                        password = pass, 
                        host = "35.221.16.125", 
                        port = 5432, 
                        dbname = "sightings")
  
  temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id
                                      FROM events
                                      JOIN places USING (place_id)
                                      JOIN counts USING (event_id)
                                      JOIN taxons USING (taxon_id)
                                      WHERE dataset_id = 'ebird'
                                      AND year > 2001
                                      AND day < 200
                                      AND lng BETWEEN -100 AND -50
                                      AND lat > 26
                                      AND (sci_name IN ('", species_list_i2[i],"'));
                                      "))
  
  #indices to fill with 1s (observations of species made for these events)
  ind <- which(data2$event_id %in% temp$event_id)
  
  #indices to fill with 0s (observations of species not made for these events)
  n_ind <- (1:NROW(data2))[-ind]
  
  #fill observersations for species i with 1s
  data2[ind, species_list_i[i,1]] <- 1
  
  #fill no observations for species i with 0s
  data2[n_ind, species_list_i[i,1]] <- 0
  
  sdata <- dplyr::select(data2, 
                         year, day, sjday, sjday2, 
                         sjday3, shr, cell, species_list_i[i,1])
  
  names(sdata)[8] <- "detect"
  sdata['species'] <- species_list_i[i,1]
  
  saveRDS(sdata, file = paste0('ebird_NA_phen_proc_', species_list_i[i,1], '.rds'))
  DBI::dbDisconnect(cxn)
}
proc.time() - tt



# Find files that weren’t created -----------------------------------------

fls <- list.files()
nms <- c()
for (i in 1:length(fls))
{
  #i <- 1
  nms <- c(nms, substr(fls[i], 20, (nchar(fls[i]) - 4)))
}

#missed species
m_sp <- species_list_i[which(!species_list_i[,1] %in% nms),1]
#remove underscore
m_sp2 <- gsub("_", " ", m_sp)

#create those files
for (i in 1:length(m_sp2))
{
  #i <- 1
  print(i)
  
  pg <- DBI::dbDriver("PostgreSQL")
  
  cxn <- DBI::dbConnect(pg, 
                        user = "cyoungflesh", 
                        password = pass, 
                        host = "35.221.16.125", 
                        port = 5432, 
                        dbname = "sightings")
  
  temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id
                                      FROM events
                                      JOIN places USING (place_id)
                                      JOIN counts USING (event_id)
                                      JOIN taxons USING (taxon_id)
                                      WHERE dataset_id = 'ebird'
                                      AND year > 2001
                                      AND day < 200
                                      AND lng BETWEEN -100 AND -50
                                      AND lat > 26
                                      AND (sci_name IN ('", m_sp2[i],"'));
                                      "))
  
  #indices to fill with 1s (observations of species made for these events)
  ind <- which(data2$event_id %in% temp$event_id)
  
  #indices to fill with 0s (observations of species not made for these events)
  n_ind <- (1:NROW(data2))[-ind]
  
  #fill observersations for species i with 1s
  data2[ind, m_sp2[i]] <- 1
  
  #fill no observations for species i with 0s
  data2[n_ind, m_sp2[i]] <- 0
  
  sdata <- dplyr::select(data2, 
                         year, day, sjday, sjday2, 
                         sjday3, shr, cell, m_sp2[i])
  
  names(sdata)[8] <- "detect"
  sdata['species'] <- m_sp2[i]
  
  saveRDS(sdata, file = paste0('ebird_NA_phen_proc_', m_sp[i], '.rds'))
}



# filter by cell ----------------------------------------------------------

read

#reference key for species synonyms
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/metadata/'))
sp_key <- read.csv('species_filenames_key.csv')

#change dir to shp files
setwd(paste0(dir, 'Bird_Phenology/Data/BirdLife_range_maps/shapefiles/'))

df_out <- data.frame()
#which species/years meet criteria for model
for (i in 1:length(species_list))
{
  #i <- 16
  
  print(i)
  #filter by species
  t_sp <- dplyr::filter(diagnostics_frame, species == species_list[i])
  
  #filter by breeding/migration cells
  #match species name to shp file name
  g_ind <- grep(species_list[i], sp_key$file_names_2016)
  
  #check for synonyms if there are no matches
  if (length(g_ind) == 0)
  {
    g_ind2 <- grep(species_list[i], sp_key$BL_Checklist_name)
  } else {
    g_ind2 <- g_ind
  }
  
  #get filename and read in
  fname <- as.character(sp_key[g_ind2,]$filenames[grep('.shp', sp_key[g_ind2, 'filenames'])])
  t_sp$shp_fname <- fname
  sp_rng <- rgdal::readOGR(fname, verbose = FALSE)
  
  #filter by breeding (2) and migration (4) range - need to convert spdf to sp
  nrng <- sp_rng[which(sp_rng$SEASONAL == 2 | sp_rng$SEASONAL == 4),]
  
  #filter by resident (1) and non-breeding (3) to exclude hex cells that contain 2/4 and 1/3
  nrng_rm <- sp_rng[which(sp_rng$SEASONAL == 1 | sp_rng$SEASONAL == 3),]
  
  #only process if there is a seasonal range
  if (NROW(nrng@data) > 0)
  {
    #good cells
    nrng_sp <- sp::SpatialPolygons(nrng@polygons)
    sp::proj4string(nrng_sp) <- sp::CRS(sp::proj4string(nrng))
    ptsreg <- sp::spsample(nrng, 50000, type = "regular")
    br_mig_cells <- as.numeric(which(!is.na(sp::over(hge, ptsreg))))
    
    #bad cells
    nrng_rm_sp <- sp::SpatialPolygons(nrng_rm@polygons)
    sp::proj4string(nrng_rm_sp) <- sp::CRS(sp::proj4string(nrng_rm))
    ptsreg_rm <- sp::spsample(nrng_rm_sp, 50000, type = "regular")
    res_ovr_cells <- as.numeric(which(!is.na(sp::over(hge, ptsreg_rm))))
    
    #remove cells that appear in resident and overwinter range that also appear in breeding range
    cell_mrg <- c(br_mig_cells, res_ovr_cells)
    to_rm <- cell_mrg[duplicated(cell_mrg)]
    
    if (length(to_rm) > 0)
    {
      overlap_cells <- br_mig_cells[-which(br_mig_cells %in% to_rm)]
    } else {
      overlap_cells <- br_mig_cells
    }
    
    #get cell centers
    cell_centers <- dggridR::dgSEQNUM_to_GEO(hexgrid6, overlap_cells)
    cc_df <- data.frame(cell = overlap_cells, lon = cell_centers$lon_deg, 
                        lat = cell_centers$lat_deg)
    
    #cells only within the range that ebird surveys were filtered to
    n_cc_df <- cc_df[which(cc_df$lon > -100 & cc_df$lon < -50 & cc_df$lat > 26),]
    cells <- n_cc_df$cell
    
    #retain rows that match selected cells
    t_sp2 <- t_sp[which(t_sp$cell %in% cells),]
    
    #create rows for cells that were missing in ebird data
    missing_cells <- cells[which(cells %ni% t_sp2$cell)]
    
    temp_dff <- t_sp2[1,]
    temp_dff[,2:20] <- NA
    
    nmc <- length(missing_cells)
    nreps <- nmc * nyr
    
    temp_dff2 <- temp_dff[rep(row.names(temp_dff), nreps),]
    rownames(temp_dff2) <- NULL
    
    temp_dff2$year <- rep(years, nmc)
    temp_dff2$cell <- rep(missing_cells, each = nyr)
    
    #combine filtered data with missing cells
    t_sp3 <- rbind(t_sp2, temp_dff2)
    
    
    #number of cells with good data in each year from 2015-2017
    nobs_yr <- c()
    for (j in 2015:2017)
    {
      #j <- 2017
      ty_sp3 <- dplyr::filter(t_sp3, year == j)
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
        ty2_sp3 <- dplyr::filter(t_sp3, year == j)
        ind2 <- which(!is.na(ty2_sp3$HM_mean))
        nobs_yr2 <- c(nobs_yr2, length(ind2))
      }
      
      #years to keep (more than three cells of data)
      yrs_kp <- years[which(nobs_yr2 >= NC)]
    }
    
    if (length(yrs_kp) > 0)
    {
      t_sp3[which(t_sp3$year %in% yrs_kp),]$MODEL <- TRUE
    }
    
    df_out <- rbind(df_out, t_sp3)
  } else {
    #merge unchanged data if there isn't a seasonal range
    df_out <- rbind(df_out, t_sp)
  }
}



# copy script to query folder for records ---------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/1-query-ebird.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path, '/1-query-ebird-', Sys.Date(), '.R'))


proc.time() - tt
