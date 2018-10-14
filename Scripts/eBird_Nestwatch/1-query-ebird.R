####################
# 1b - query ebird data
#
# Filters eBird data from DB
# Zero fills
# Creates directory of processed data (rds file for each species) and copy of this 
# ... script (for reproducability) in /Data/Processed/db_query_<DATE>
#
# Replaces 1-process-ebird-data.R, which processed data based on local copy of eBird reference dataset
####################


# #can access DB from command line with:
# psql "sslmode=disable dbname=sightings user=cyoungflesh hostaddr=35.221.16.125"


dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(dggridR)

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




# Query and filter ----------------------------------------------------

#*get info for each unique event
#*create 153 species columns with NA
#*in a loop, query each species individually and fill obs with 1s, no obs with 0s
#*merge with hex cells
#*create rds objects for each species


#filter all unique event_ids that meet criteria - about 42 min to complete query
tt <- proc.time()
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
proc.time() - tt



#only values with unique group identifiers (or no group identifier)
data2 <- data[!duplicated(data[,'group_identifier'], 
                          incomparables = NA),]

#create species columns
data2[species_list_i[,1]] <- NA

#zero fill
nsp <- NROW(species_list_i)

#~2-3 minutes for each species query
tt <- proc.time()
for (i in 1:nsp)
{
  #i <- 1
  print(i)
  
  temp <- DBI::dbGetQuery(cxn, paste0("
                                      SELECT event_id, year, day, place_id, lat, lng, started, radius, 
                                      sci_name, common_name, count,
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
                                      AND (sci_name IN ('", species_list_i2[i],"'))
                                      AND (event_json ->> 'ALL_SPECIES_REPORTED')::int = 1
                                      AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
                                      AND LEFT(started, 2)::int < 16
                                      AND RADIUS < 100000;
                                      "))
  
  #indices to fill with 1s (observations of species made for these events)
  ind <- which(data2$event_id %in% temp$event_id)
  
  #indices to fill with 0s (observations of species not made for these events)
  n_ind <- (1:NROW(data2))[-ind]
  
  #fill observersations for species i with 1s
  data2[ind, species_list_i[i,1]] <- 1
  
  #fill no observations for species i with 0s
  data2[n_ind, species_list_i[i,1]] <- 0
}
proc.time() - tt


# add jday^2, jday^3 and shr ---------------------------------------------------

#scaled julian day, scaled julian day^2, and scaled julian day^3
SJDAY  <- as.vector(scale(data2$day, scale = FALSE))
SJDAY2 <- as.vector(scale(data2$day^2, scale = FALSE))
SJDAY3 <- as.vector(scale(data2$day^3, scale = FALSE))

data2$sjday <- SJDAY
data2$sjday2 <- SJDAY2
data2$sjday3 <- SJDAY3

#scaled effort hours
SHR <- as.vector(scale((data2$duration_minutes/60)))

data2$shr <- SHR



# bin to hex grid ---------------------------------------------------------

# Construct geospatial hexagonal grid

hexgrid6 <- dggridR::dgconstruct(res = 6) 
data2$cell6 <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                        in_lon_deg = data2$lng, 
                                        in_lat_deg = data2$lat)[[1]]

query_dir_path <- paste0('Processed/db_query_', Sys.Date())

dir.create(query_dir_path)
setwd(query_dir_path)

#save to rds object
saveRDS(data2, file = 'ebird_NA_phen_proc_ALL_SPECIES.rds')




# Different rds object each species ---------------------------------------

for(i in 1:nsp)
{
  #i <- 1
  sdata <- dplyr::select(data2, 
                         year, day, sjday, sjday2, 
                         sjday3, shr, cell6, species_list_i[i,1])
  
  names(sdata)[8] <- "detect"
  
  saveRDS(sdata, file = paste0('ebird_NA_phen_proc_', species_list_i[i,1], '.rds'))
}

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/1b-query-ebird.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path))

