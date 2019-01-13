######################
# X - check ebird breeding code data
#
# filter ebird basic dataset (on desktop) for one species, only breeding codes C3/C4 (probable/confirmed breeding)
# check availability of breeding codes over time and compare with results from DB query
#
######################

#Run using bash - issue runningon cluster
#awk -F'\t' '$6=="Contopus virens" {print $0}' ebd_relFeb-2018.txt > Contopus_virens.txt
#awk -F'\t' '$11=="C3" || $11=="C4" {print $0}' Contopus_virens.txt > C34_Contopus_virens.txt


# load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)


# read in sorted br code data ---------------------------------------------


setwd('~/Desktop/ebd_relFeb-2018/')


data <- read.csv('C34_Contopus_virens.txt', 
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

  

colnames(data) <- gsub(' ', '_', cn)

d_filt <- dplyr::select(data, c('SCIENTIFIC_NAME', 'BREEDING_BIRD_ATLAS_CATEGORY', 'LATITUDE', 'LONGITUDE', 
                            'OBSERVATION_DATE', 'TIME_OBSERVATIONS_STARTED', 'PROJECT_CODE', 'DURATION_MINUTES', 
                            'EFFORT_DISTANCE_KM', 'ALL_SPECIES_REPORTED', 'APPROVED'))




# new cols for year and day -----------------------------------------------


d_filt$YEAR <- as.numeric(format(as.Date(d_filt[,'OBSERVATION_DATE']), '%Y'))
d_filt$DAY <- as.numeric(format(as.Date(d_filt[,'OBSERVATION_DATE']), '%j'))
d_filt$TIME <- as.numeric(substr(d_filt[,'TIME_OBSERVATIONS_STARTED'], start = 1, stop = 2))



# filter based on relevant fields -----------------------------------------


d_filt2 <- dplyr::filter(d_filt, 
              YEAR > 2001,
              DAY < 200,
              LONGITUDE > -100,
              LONGITUDE < -50,
              LATITUDE > 26,
              ALL_SPECIES_REPORTED == 1,
              TIME < 16,
              DURATION_MINUTES > 6,
              DURATION_MINUTES < 1440,
              EFFORT_DISTANCE_KM < 100)

plyr::count(d_filt2, vars = 'YEAR')


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

temp_db <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, sci_name,
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
                                    AND (sci_name IN ('Contopus virens'));
                                    "))


DATE_BC <- '2019-01-09'

setwd(paste0('~/Google_Drive/R/Bird_phenology/Data/Processed/breeding_cat_query_', DATE_BC))
C_virens <- readRDS('ebird_NA_breeding_cat_Contopus_virens.rds')

tt <- dplyr::filter(C_virens, bba_breeding_category == 'C3' | bba_breeding_category == 'C4')

plyr::count(tt, vars = 'year')

