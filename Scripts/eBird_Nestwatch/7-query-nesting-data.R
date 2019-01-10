####################
# 7 - query other nesting phenology data
#
# *eBird breeding codes
#   -0 if not observed in that survey at all (independent of breeding code)
#   -NA if observed but no breeding code recorded
#   -letter code if observed and breeding code recorded
# *MAPS observations
# *MAPS data
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(dggridR)
library(doParallel)
library(foreach)


# import IAR data ---------------------------------------------------------

# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
# 
# IAR_out_dir <- 'IAR_output_2018-11-12'
# IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)

#IAR_data <- readRDS(paste0('master_arrival_', IAR_out_date, '.rds'))



# import IAR species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)

#remove underscore and coerce to vector
species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))


#combine species names into a single string with quotes
SL <- paste0("'", species_list_i2, "'", collapse = ", ")


# eBird breeding codes ----------------------------------------------------


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


query_dir_path <- paste0('Processed/breeding_cat_query_', Sys.Date())

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

#filter all unique event_ids that meet criteria - about 38 min to complete query

data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT DISTINCT ON (event_id) event_id, year, day, place_id, lat, lng, started, 
                                    radius,
                                    (event_json ->> 'ALL_SPECIES_REPORTED')::int AS all_species_reported,
                                    (event_json ->> 'DURATION_MINUTES')::int AS duration_minutes,
                                    count_json ->> 'OBSERVER_ID' AS observer_id,
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


#calculate polynomial then center data
#scaled julian day, scaled julian day^2, and scaled julian day^3
SJDAY  <- as.vector(scale(data2$day, scale = FALSE))
SJDAY2 <- as.vector(scale(data2$day^2, scale = FALSE))
SJDAY3 <- as.vector(scale(data2$day^3, scale = FALSE))

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

#run in parallel with 6 logical cores
doParallel::registerDoParallel(cores = 4)

tt <- proc.time()
foreach::foreach(i = 1:nsp) %dopar%
{
  #i <- 96
  print(i)
  
  pg <- DBI::dbDriver("PostgreSQL")
  
  cxn <- DBI::dbConnect(pg, 
                        user = "cyoungflesh", 
                        password = pass, 
                        host = "35.221.16.125", 
                        port = 5432, 
                        dbname = "sightings")
  
  temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, sci_name,
                                      count_json ->> 'BREEDING_BIRD_ATLAS_CODE' AS bba_code,
                                      count_json ->> 'BREEDING_BIRD_ATLAS_CATEGORY' AS bba_category
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
  
  #observations of species made for these events) - indices
  ind <- which(data2$event_id %in% temp$event_id)
  
  #observations of species not made for these events - indices
  n_ind <- (1:NROW(data2))[-ind]
  
  
  #0 if not observed in that survey at all (independent of breeding code)
  #NA if observed but no breeding code recorded
  #letter code if observed and breeding code recorded
  data2[ind, species_list_i[i,1]] <- temp$bba_category[ind]
  
  #fill no observations for species i with 0s
  data2[n_ind, species_list_i[i,1]] <- 0
  
  sdata <- dplyr::select(data2, 
                         year, day, cell, sjday, sjday2, 
                         sjday3, shr, species_list_i[i,1])
  
  names(sdata)[8] <- "bba_breeding_category"
  sdata['species'] <- species_list_i[i,1]
  
  saveRDS(sdata, file = paste0('ebird_NA_breeding_cat_', species_list_i[i,1], '.rds'))
  
  DBI::dbDisconnect(cxn)
}
proc.time() - tt



# Find files that weren’t created -----------------------------------------

fls <- list.files()
nms <- c()
for (i in 1:length(fls))
{
  #i <- 1
  nms <- c(nms, substr(fls[i], 23, (nchar(fls[i]) - 4)))
}

#missed species
m_sp <- species_list_i[which(!species_list_i[,1] %in% nms),1]
#remove underscore
m_sp2 <- gsub("_", " ", m_sp)

if (length(m_sp2) > 0)
{
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
    
    temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, sci_name,
                                        count_json ->> 'BREEDING_BIRD_ATLAS_CODE' AS bba_code,
                                        count_json ->> 'BREEDING_BIRD_ATLAS_CATEGORY' AS bba_category
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
    
    
    #observations of species made for these events) - indices
    ind <- which(data2$event_id %in% temp$event_id)
    
    #observations of species not made for these events - indices
    n_ind <- (1:NROW(data2))[-ind]
    
    
    #0 if not observed in that survey at all (independent of breeding code)
    #NA if observed but no breeding code recorded
    #letter code if observed and breeding code recorded
    data2[ind, m_sp2[i]] <- temp$bba_category[ind]
    
    #fill no observations for species i with 0s
    data2[n_ind, m_sp2[i]] <- 0
    
    sdata <- dplyr::select(data2, 
                           year, day, sjday, sjday2, 
                           sjday3, shr, cell, m_sp2[i])
    
    names(sdata)[8] <- "bba_breeding_category"
    sdata['species'] <- m_sp2[i]
    
    saveRDS(sdata, file = paste0('ebird_NA_breeding_cat_', msp[i,1], '.rds'))
    DBI::dbDisconnect(cxn)
  }
}


# copy script to query folder for records ---------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/7-query-nesting-data.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path, '/7-query-nesting-data-', Sys.Date(), '.R'))


proc.time() - tt









# MAPS obs ----------------------------------------------------------------


setwd(paste0(dir, 'Bird_phenology/Data/MAPS_Obs'))

MAPS_obs <- read.csv('1117M.csv')

#0 - observed but not breeding
#P - probably breeding
#C - confirmed breeding

#PS1 - May 1-May 10
#PS2 - May 11-May 20
#PS3 - May 21-May 30
#PS4 - May 31-June 9
#PS5 - June 10-June 19
#PS6 - June 20-June 29
#PS7 - June 30-July 9
#PS8 - July 10-July 19
#PS9 - July 20-July 29
#PS10 - July 30-August 8
#PS11 - August 9-August 18


#get grid cell of each station by lat/lon
setwd(paste0(dir, 'Bird_phenology/Data/MAPS_Obs/CntrlStations'))

MAPS_stations <- read.csv('STATIONS.csv', skipNul = TRUE)

MAPS_mrg <- dplyr::left_join(MAPS_obs, MAPS_stations, by = 'LOC')

MAPS_mrg$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                       in_lon_deg = MAPS_mrg$DECLNG, 
                                       in_lat_deg = MAPS_mrg$DECLAT)[[1]]

#species codes - merge with MAPS_mrg

#filter by species, get breeding date (period) for each year/cell



