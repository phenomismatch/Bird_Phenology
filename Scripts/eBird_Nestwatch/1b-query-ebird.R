####################
#Sample DB access code from Rafe
####################

#https://db.rstudio.com/dplyr/



cy_dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(dbplyr)
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



# access DB ---------------------------------------------------------------

pass <- as.character(read.lines('db_pass.txt'))
pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")

# filter dataset ----------------------------------------------------------

data <- dbGetQuery(cxn, paste0("
                   SELECT event_id, year, day, place_id, lat, lng, started, count, radius, sci_name, common_name,
                   event_json ->> 'GROUP_IDENTIFIER' AS group_identifier,
                   event_json -> 'ALL_SPECIES_REPORTED' AS all_species_reported,
                   event_json -> 'DURATION_MINUTES' AS duration_minutes,
                   event_json -> 'OBSERVER_ID' AS observer_id,
                   event_json -> 'NUMBER_OBSERVERS' AS number_observers
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
                   LIMIT 1000;
                   "))

#AND radius < 100



#filter by only all_species_report == 1
#filter by started
#filter by duration_minutes
#select unqiue group IDs when not NULL


