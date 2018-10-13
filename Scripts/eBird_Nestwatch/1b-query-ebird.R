####################
#Sample DB access code from Rafe
####################

#https://db.rstudio.com/dplyr/



# #can access DB from command line with:
# psql "sslmode=disable dbname=sightings user=cyoungflesh hostaddr=35.221.16.125"


cy_dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)


# set wd ------------------------------------------------------------------

setwd(paste0(cy_dir, 'Bird_Phenology/Data/'))



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

#FILTER BY:
# * only ebird data
# * species of interest (import from text file)
# * only complete checklists
# * duration minutes (6 min to 24 hours)
# * lat (> 26) and lon (-100 to -50)
# * year > 2001
# * day of year (< julian day 200)
# * time started before 16:00
# * radius < 100km


#NOT INCORPORATED:
# * only unique groups when not null


#FROM RAFE
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

data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT event_id, year, day, place_id, lat, lng, started, 
                                    count, radius, sci_name, common_name,
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
                                    AND RADIUS < 100000
                                    LIMIT 500;
                                    "))

str(data)
head(data)


#*get info for each unique event
#*create 153 species columns with NA
#*in a loop, query each species individually and fill NAs
#*change all remaining NA values to 0
#*merge with hex cells
#*create rds objects for each species


#all unique event_ids that meet criteria
tt <- proc.time()
data2 <- DBI::dbGetQuery(cxn, paste0("
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

str(data2)
head(data2)

#select unique group IDs when not NULL
unique((data2$event_id))
#zero fill


