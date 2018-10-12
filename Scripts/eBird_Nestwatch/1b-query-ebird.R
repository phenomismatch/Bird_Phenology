####################
#Sample DB access code from Rafe
####################

#https://db.rstudio.com/dplyr/


#NEED TO KNOW HOW TO:

#filter by only all_species_report == 1
#filter by started
#filter by duration_minutes
#select unique group IDs when not NULL


# #access DB from command line
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
SL <- 'Empidonax virescens'


# access DB ---------------------------------------------------------------

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")

# filter dataset ----------------------------------------------------------

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


#NEED RAFE:
#Bird atlas breeding code?
#Protocol type?



# Filter and zero fill ----------------------------------------------------

#*get event_ids
#*create blank df with that many rows
#*query each species individually and fill blank df (THROUGH LOOP)
#*change all NA values to 0
#*merge hex cells
#*create rds objects for each species


data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT event_id, year, day, place_id, lat, lng, started, 
                                    count, radius, sci_name, common_name,
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
                                    AND (sci_name IN ('", SL,"'))
                                    AND (event_json ->> 'ALL_SPECIES_REPORTED')::int = 1
                                    AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
                                    AND LEFT(started, 2)::int < 16
                                    AND RADIUS < 100000;
                                    "))

str(data)
head(data)


