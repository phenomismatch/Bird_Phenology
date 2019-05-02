###############################
# Find issues in MAPS database
#
# ids that have multiple species names
# Sting from 199172901 - 199173000 shoudl not be used
###############################



# top-level dir -----------------------------------------------------------

dir <- '~/Google_Drive/R/'



# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(ggplot2)
library(rstan)


# query DB ----------------------------------------------------------------


setwd(paste0(dir, 'Bird_Phenology/Data/'))

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg,
                      user = "cyoungflesh",
                      password = pass,
                      host = "35.221.16.125",
                      port = 5432,
                      dbname = "sightings")


maps_data_full <- DBI::dbGetQuery(cxn, paste0("SELECT lng, lat, year, day, common_name, sci_name, event_id, started, ended,
                                         (count_json ->> 'ANET') AS anet,
                                         (event_json ->> 'STATION') AS station,
                                         (place_json ->> 'HABITAT') AS habitat,
                                         (place_json ->> 'ELEV')::int AS elev,
                                         (count_json ->> 'C') capture_code,
                                         (count_json ->> 'BAND') AS band_id,
                                         (count_json ->> 'AGE')::int AS age,
                                         (count_json ->> 'SEX') AS sex,
                                         (count_json ->> 'FW') AS feather_wear,
                                         (count_json ->> 'CP')::int AS cloacal_pro,
                                         (count_json ->> 'BP')::int AS brood_patch,
                                         (count_json ->> 'F')::int AS fat_content,
                                         (count_json ->> 'WNG')::float AS wing_chord,
                                         (count_json ->> 'WEIGHT')::float AS weight,
                                         (count_json ->> 'N') AS standard_effort
                                         FROM places
                                         JOIN events USING (place_id)
                                         JOIN counts USING (event_id)
                                         JOIN taxa USING (taxon_id)
                                         WHERE events.dataset_id = 'maps';
                                         "))



# find band_ids with multiple species -------------------------------------

un_bid <- unique(maps_data_full$band_id)

bad_ids <- data.frame()
pb <- txtProgressBar(min = 0, max = length(un_bid), style = 3)
for (i in 1:length(un_bid))
{
  #i <- 1
  temp <- dplyr::filter(maps_data_full, band_id == un_bid[i])
  
  if (length(unique(temp$sci_name)) > 1)
  {
    bad_ids <- rbind(bad_ids, temp)
  }
  setTxtProgressBar(pb, i)
}
close(pb)


#setwd('~/Desktop')

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# saveRDS(bad_ids, 'MAPS_bad_ids.rds')
bad_ids <- readRDS('MAPS-bad-ids.rds')

head(bad_ids)
bad_ids2 <- dplyr::filter(bad_ids, standard_effort == '?', band_id < 199172901, band_id > 199173000)
NROW(bad_ids2)


output <- dplyr::select(bad_ids, year, day, band_id, common_name, sci_name, station)

write.csv(output, 'MAPS-bad-ids.csv', row.names = FALSE)

# find sites with lat/lon 0,0 ---------------------------------------------


#EASA has lat/lon 0,0
llz <- dplyr::filter(maps_data_full, lat == 0, lng == 0)
colnames(llz)[24] <- 'N'

setwd('~/Desktop')
write.csv(llz, 'MAPS-EASA.csv', row.names = FALSE)

