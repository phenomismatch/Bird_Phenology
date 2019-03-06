##############################
# query ebird database for Val
##############################

#Lat/longs: -74.0267, -70.82259, 40.7584, 42.27485  (xmin, xmax, ymin, ymax)
#Years: 2013-2018
#May through July (05-01 to 07-31)
#Yes, zero-filled. Complete checklists only (or retain that field so I can filter later). Yes, surveys as rows with species in columns. It's probably safest to retain all the eBird fields (BREEDING.BIRD.ATLAS.CODE, STATE_PROVINCE, ATLAS.BLOCK, LOCALITY, LOCALITY.ID, LOCALITY.TYPE, LATITUDE, LONGITUDE, OBSERVATION.DATE, TIME.OBSERVATIONS.STARTED, PROTOCOL.TYPE, DURATION.MINUTES, EFFORT.DISTANCE.KM, .. etc. )

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


query_dir_path <- paste0('Processed/eBird_query_Val')

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
# * time started before 18:00
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

# #some years will be one day later earlier to start and one day earlier to end (due to Leap Years)
# #regular year - 121-212
# format(as.Date('2019-05-01'), '%j')
# format(as.Date('2019-07-31'), '%j')
# #leap year - 122-213
# format(as.Date('2020-05-01'), '%j')
# format(as.Date('2020-07-31'), '%j')


#get list of all species
data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT event_id, year, day, place_id, lat, lng,
                                    event_json, place_json, count_json
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxa USING (taxon_id)
                                    WHERE year BETWEEN 2013 AND 2018
                                    AND events.dataset_id = 'ebird'
                                    AND day BETWEEN 121 AND 213
                                    AND lng BETWEEN -74.5 AND -70
                                    AND lat BETWEEN 40 AND 42.5
                                    LIMIT 1000;
                                    "))



#get list of all event_ids
data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT event_id, year, day, place_id, lat, lng,
                                    event_json, place_json, count_json
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxa USING (taxon_id)
                                    WHERE day BETWEEN 121 AND 213
                                    AND dataset_id = 'ebird'
                                    AND year BETWEEN 2013 AND 2018
                                    AND lng BETWEEN -75 AND -70
                                    AND lat BETWEEN 40 AND 43
                                    LIMIT 1000;
"))






#only values with unique group identifiers (or no group identifier)
data2 <- data[!duplicated(data[,'group_identifier'], 
                          incomparables = NA),]