######################
# X - check ebird breeding code data
#
# filter ebird basic dataset (on desktop) for one species, only breeding codes C3/C4 (probable/confirmed breeding)
# check availability of breeding codes over time and compare with results from DB query
#
######################

#Run using bash - issue running on cluster
#awk -F'\t' '$6=="Contopus virens" {print $0}' ebd_relFeb-2018.txt > Contopus_virens.txt
#awk -F'\t' '$11=="C3" || $11=="C4" {print $0}' Contopus_virens.txt > C34_Contopus_virens.txt


# load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)


# read in sorted br code data ---------------------------------------------


setwd('~/Desktop/ebd_relFeb-2018/')

raw_data <- read.csv('Contopus_virens.txt', 
                 header = FALSE,
                 sep = '\t',
                 quote= '',
                 stringsAsFactors = FALSE)


ebd_data <- read.csv('C34_Contopus_virens.txt', 
                 header = FALSE,
                 sep = '\t',
                 quote= '',
                 stringsAsFactors = FALSE)



# colnames from EBD -------------------------------------------------------

cn <- read.csv('ebd_relFeb-2018.txt', 
                 header = FALSE,
                 sep = '\t',
                 quote= '',
                 nrows = 1,
                 stringsAsFactors = FALSE)

  
colnames(raw_data) <- gsub(' ', '_', cn)
colnames(ebd_data) <- gsub(' ', '_', cn)

rd_filt <- dplyr::select(raw_data[,-50], c('SCIENTIFIC_NAME', 'BREEDING_BIRD_ATLAS_CATEGORY', 'LATITUDE', 'LONGITUDE', 
                                           'OBSERVATION_DATE', 'TIME_OBSERVATIONS_STARTED', 'PROJECT_CODE', 'DURATION_MINUTES', 
                                           'EFFORT_DISTANCE_KM', 'ALL_SPECIES_REPORTED', 'APPROVED'))


d_filt <- dplyr::select(ebd_data, c('SCIENTIFIC_NAME', 'BREEDING_BIRD_ATLAS_CATEGORY', 'LATITUDE', 'LONGITUDE', 
                            'OBSERVATION_DATE', 'TIME_OBSERVATIONS_STARTED', 'PROJECT_CODE', 'DURATION_MINUTES', 
                            'EFFORT_DISTANCE_KM', 'ALL_SPECIES_REPORTED', 'APPROVED'))


# new cols for year and day -----------------------------------------------


rd_filt$YEAR <- as.numeric(format(as.Date(rd_filt[,'OBSERVATION_DATE']), '%Y'))
rd_filt$DAY <- as.numeric(format(as.Date(rd_filt[,'OBSERVATION_DATE']), '%j'))
rd_filt$TIME <- as.numeric(substr(rd_filt[,'TIME_OBSERVATIONS_STARTED'], start = 1, stop = 2))

d_filt$YEAR <- as.numeric(format(as.Date(d_filt[,'OBSERVATION_DATE']), '%Y'))
d_filt$DAY <- as.numeric(format(as.Date(d_filt[,'OBSERVATION_DATE']), '%j'))
d_filt$TIME <- as.numeric(substr(d_filt[,'TIME_OBSERVATIONS_STARTED'], start = 1, stop = 2))



# filter based on relevant fields -----------------------------------------


d_filt2 <- dplyr::filter(d_filt, 
              YEAR > 2001,
              DAY < 200,
              LONGITUDE > -100,
              LONGITUDE < -50,
              LATITUDE > 26)#,
              #ALL_SPECIES_REPORTED == 1,
              #TIME < 16,
              #DURATION_MINUTES > 6,
              #DURATION_MINUTES < 1440,
              #EFFORT_DISTANCE_KM < 100)

plyr::count(d_filt2, vars = 'YEAR')


rd_filt2 <- dplyr::filter(rd_filt, 
                         YEAR > 2001,
                         DAY < 200,
                         LONGITUDE > -100,
                         LONGITUDE < -50,
                         LATITUDE > 26)#,
#ALL_SPECIES_REPORTED == 1,
#TIME < 16,
#DURATION_MINUTES > 6,
#DURATION_MINUTES < 1440,
#EFFORT_DISTANCE_KM < 100)


head(raw_data)
#check subspecies common names




# db query ---------------------------------------------------

setwd('~/Google_Drive/R/Bird_Phenology/Data')

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")


temp_db <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT event_id, year, day, place_id, lat, lng, started, 
                                       radius, sci_name, 
                                       (event_json ->> 'ALL_SPECIES_REPORTED')::int AS all_species_reported,
                                       (event_json ->> 'DURATION_MINUTES')::int AS duration_minutes,
                                       count_json ->> 'BREEDING_BIRD_ATLAS_CODE' AS bba_code,
                                       count_json ->> 'BREEDING_BIRD_ATLAS_CATEGORY' AS bba_category,
                                       count_json ->> 'GLOBAL_UNIQUE_IDENTIFIER' AS global_unique_identifier,
                                       (event_json ->> 'NUMBER_OBSERVERS')::int AS number_observers
                                       FROM places
                                       JOIN events USING (place_id)
                                       JOIN counts USING (event_id)
                                       JOIN taxons USING (taxon_id)
                                       WHERE dataset_id = 'ebird'
                                       AND (sci_name IN ('Contopus virens'));
                                       "))

#compare this query to raw_data
head(temp_db)
head(raw_data)

NROW(temp_db)
NROW(raw_data)

#just to check - same as the number of rows
#length(unique(temp_db$global_unique_identifier))
#length(unique(raw_data[,1]))


#which IDs (ebird generateed ids) are missing in db?

db_id <- temp_db$global_unique_identifier
rd_id <- raw_data$GLOBAL_UNIQUE_IDENTIFIER

#are all the ids from db query in the raw data?
sum(db_id %in% rd_id) #YES

'%ni%' <- Negate('%in%')

#which ids are missing from db query
rd_id[which(rd_id %ni% db_id)] #THESE ARE THE MISSING EBIRD IDS


#sample missing ids - which raw data ARE NOT in query
raw_data[head(which(rd_id %ni% db_id)),]
#sample non-missing ids - which raw data ARE in query
raw_data[head(which(rd_id %in% db_id)),]




#filter for C3/C4
raw_data$YEAR <- as.numeric(format(as.Date(raw_data[,'V28']), '%Y'))
rd_filt <- dplyr::filter(raw_data, V11 == 'C3' | V11 == 'C4')

db_filt <- dplyr::filter(temp_db, bba_category == 'C3' | bba_category == 'C4')


plyr::count(db_filt, vars = 'year')
plyr::count(rd_filt, vars = 'YEAR')




#then compare C3/C4 to data

NROW(temp_db_filt)
NROW(data)



# temp_db <- DBI::dbGetQuery(cxn, paste0("
#                                     SELECT year, day, place_id, lat, lng, started, 
#                                     radius,
#                                     (event_json ->> 'ALL_SPECIES_REPORTED')::int AS all_species_reported,
#                                     (event_json ->> 'DURATION_MINUTES')::int AS duration_minutes,
#                                     count_json ->> 'BREEDING_BIRD_ATLAS_CODE' AS bba_code,
#                                     count_json ->> 'BREEDING_BIRD_ATLAS_CATEGORY' AS bba_category,
#                                     (event_json ->> 'NUMBER_OBSERVERS')::int AS number_observers
#                                     FROM places
#                                     JOIN events USING (place_id)
#                                     JOIN counts USING (event_id)
#                                     JOIN taxons USING (taxon_id)
#                                     WHERE dataset_id = 'ebird'
#                                     AND year > 2001
#                                     AND day < 200
#                                     AND lng BETWEEN -100 AND -50
#                                     AND lat > 26
#                                     AND (sci_name IN ('Contopus virens'))
#                                     AND (event_json ->> 'ALL_SPECIES_REPORTED')::int = 1
#                                     AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
#                                     AND LEFT(started, 2)::int < 16
#                                     AND RADIUS < 100000;
#                                     "))

NROW(temp_db)

db_filt <- dplyr::filter(temp_db, bba_category == 'C3' | bba_category == 'C4')

plyr::count(db_filt, vars = 'year')


# DATE_BC <- '2019-01-09'
# 
# setwd(paste0('~/Google_Drive/R/Bird_phenology/Data/Processed/breeding_cat_query_', DATE_BC))
# C_virens <- readRDS('ebird_NA_breeding_cat_Contopus_virens.rds')
# 
# tt <- dplyr::filter(C_virens, bba_breeding_category == 'C3' | bba_breeding_category == 'C4')
# 
# plyr::count(tt, vars = 'year')
# 
