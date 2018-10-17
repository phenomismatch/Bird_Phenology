###################
#Investigate BBS issue
#
###################

# Set dir -----------------------------------------------------------------

dir <- '~/Google_Drive/R/'

setwd(paste0(dir, 'Bird_Phenology/Data/'))



# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)



# species -----------------------------------------------------------------

SL <- c('Cardinalis_cardinalis', 'Tachycineta_bicolor')
SL_sp <- gsub("_", " ", SL)


# db connection ----------------------------------------------------------------


pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")

# query bbs database' -----------------------------------------------------


i <- 1

temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, year, day, 
                                    started, ended, lng, lat, common_name, count,
                                    (place_json ->> 'countrynum')::int AS country_number,
                                    (place_json ->> 'statenum')::int AS state_number,
                                    (place_json ->> 'route')::int AS route_number,
                                    (place_json ->> 'routename') AS route_name,
                                    (place_json ->> 'stratum') AS straum,
                                    (place_json ->> 'bcr') AS bcr,
                                    (place_json ->> 'routetypeid') AS routetypeid,
                                    (place_json ->> 'routetypedetailid') AS routetypedetailid,
                                    (event_json ->> 'routedataid') AS routedataid,
                                    (event_json ->> 'rpid') AS rpid,
                                    (event_json ->> 'obsn') AS obsn,
                                    (event_json ->> 'totalspp') AS totalspp,
                                    (count_json ->> 'record_id') AS record_id
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxons USING (taxon_id)
                                    WHERE dataset_id = 'bbs'
                                    AND year = 2016
                                    AND sci_name = '", SL_sp[i],"';
                                    "))


# different nnmbers of rows -----------------------------------------------

NROW(temp)
length(unique(temp$event_id))
NROW(unique(temp[c('country_number', 'state_number', 'route_number', 'day')]))



# find duplicates ---------------------------------------------------------

dups <- which(duplicated(temp[c('country_number', 'state_number', 'route_number', 'day')]))

dup_ind <- sort(c(dups, dups - 1))

data_dup <- temp[dup_ind,]

setwd('~/Desktop/')
write.csv(data_dup, 'bbs_data_issues_2018-10-17.csv', row.names = FALSE)


