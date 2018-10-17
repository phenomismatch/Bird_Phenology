###################
# Query DB for BBS data 
#
# Queries and zero fills BBS data for Northern Cardinal and Tree Swallow
# Some species counts are entered twice for some routes. Counts are summed across species for each route
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



# Get all routes run in 2016 ----------------------------------------------------

#each route is its own event_id

data <- DBI::dbGetQuery(cxn, paste0("SELECT DISTINCT ON (event_id) event_id, year, day, 
                                    started, ended, lng, lat,
                                    (place_json ->> 'countrynum')::int AS country_number,
                                    (place_json ->> 'statenum')::int AS state_number,
                                    (place_json ->> 'route')::int AS route_number,
                                    (place_json ->> 'routename') AS route_name
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxons USING (taxon_id)
                                    WHERE dataset_id = 'bbs'
                                    AND year = 2016;
                                    "))

#add species columns
data[SL] <- NA

for (i in 1:2)
{
  #i <- 1
  #multiple species records in events -> sum count across event
  temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id,
                                      SUM (count)
                                      FROM places
                                      JOIN events USING (place_id)
                                      JOIN counts USING (event_id)
                                      JOIN taxons USING (taxon_id)
                                      WHERE dataset_id = 'bbs'
                                      AND year = 2016
                                      AND sci_name = '", SL_sp[i],"'
                                      GROUP BY event_id;
                                      "))
  
  #indices to fill with 1s (observations of species made for these events)
  ind <- which(data$event_id %in% temp$event_id)
  
  #indices to fill with 0s (observations of species not made for these events)
  n_ind <- (1:NROW(data))[-ind]
  
  #fill observersations for species i with 1s
  data[ind, SL[i]] <- temp[,'sum']
  
  #fill no observations for species i with 0s
  data[n_ind, SL[i]] <- 0
}


#only TRSW
TRSW <- data[, -which(names(data) %in% SL[1])]
#only NOCA
NOCA <- data[, -which(names(data) %in% SL[2])]



# write to file -----------------------------------------------------------


setwd('~/Desktop/')

write.csv(TRSW, 'TRSW_2016_bbs_ZF.csv', row.names = FALSE)
write.csv(NOCA, 'NOCA_2016_bbs_ZF.csv', row.names = FALSE)


