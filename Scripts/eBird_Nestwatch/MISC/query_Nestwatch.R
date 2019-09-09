######################
# ? - Query Nestwatch data
#
######################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'

DATE <- '2019-08-27'



# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(dggridR)



# Nestwatch data - query db -----------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")


nestwatch_data <- DBI::dbGetQuery(cxn, paste0("SELECT lng, lat, year, day, common_name, sci_name,
                                              event_id, count_id, count,
                                              (event_json ->> 'FIRST_LAY_DT') AS lay,
                                              (event_json ->> 'HATCH_DT') AS hatch,
                                              (event_json ->> 'FLEDGE_DT') AS fledge,
                                              (event_json ->> 'event_type') AS event_type,
                                              (count_json ->> 'count_type') AS count_type,
                                              (event_json ->> 'ATTEMPT_ID') AS attempt_id
                                              FROM places
                                              JOIN events USING (place_id)
                                              JOIN counts USING (event_id)
                                              JOIN taxa USING (taxon_id)
                                              WHERE events.dataset_id = 'nestwatch';
                                              "))



# Process Nestwatch data ---------------------------------------------------------------

nestwatch_data$sci_name <- gsub(' ', '_', nestwatch_data$sci_name)
uaid <- unique(nestwatch_data$attempt_id)

nw_pro <- data.frame(species = rep(NA, length(uaid)), year = NA, lat = NA, lng = NA, 
                     attempt_id = NA, lay = NA, hatch = NA, fledge = NA,
                     n_eggs = NA, n_chicks = NA, n_unh_eggs = NA, n_fledged = NA,
                     n_chicks_dead = NA)

pb <- txtProgressBar(min = 0, max = length(uaid), style = 3)
for (i in 1:length(uaid))
{
  #i <- 1
  temp <- dplyr::filter(nestwatch_data, attempt_id == uaid[i])
  
  nw_pro$species[i] <- temp$sci_name[1]
  nw_pro$year[i] <- temp$year[1]
  nw_pro$lat[i] <- temp$lat[1]
  nw_pro$lng[i] <- temp$lng[1]
  nw_pro$attempt_id[i] <- temp$attempt_id[1]
  
  nw_pro$lay[i] <- as.numeric(format(as.Date(temp$lay, format = '%d-%b-%Y'), '%j'))[1]
  nw_pro$hatch[i] <- as.numeric(format(as.Date(temp$hatch, format = '%d-%b-%Y'), '%j'))[1]
  nw_pro$fledge[i] <- as.numeric(format(as.Date(temp$fledge, format = '%d-%b-%Y'), '%j'))[1]
  
  n_eggs <- dplyr::filter(temp, count_type == 'CLUTCH_SIZE_HOST_ATLEAST')
  n_chicks <- dplyr::filter(temp, count_type == 'YOUNG_HOST_TOTAL_ATLEAST')
  n_unh_eggs <- dplyr::filter(temp, count_type == 'EGGS_HOST_UNH_ATLEAST')
  n_fledged <- dplyr::filter(temp, count_type == 'YOUNG_HOST_FLEDGED_ATLEAST')
  n_chicks_dead <- dplyr::filter(temp, count_type == 'YOUNG_HOST_DEAD_ATLEAST')
  
  if (NROW(n_eggs) > 0) nw_pro$n_eggs[i] <- n_eggs$count
  if (NROW(n_chicks) > 0) nw_pro$n_chicks[i] <- n_chicks$count
  if (NROW(n_unh_eggs) > 0) nw_pro$n_unh_eggs[i] <- n_unh_eggs$count
  if (NROW(n_fledged) > 0) nw_pro$n_fledged[i] <- n_fledged$count
  if (NROW(n_chicks_dead) > 0) nw_pro$n_chicks_dead[i] <- n_chicks_dead$count
  setTxtProgressBar(pb, i)
}
close(pb)



#create grid and add to df
hexgrid6 <- dggridR::dgconstruct(res = 6) 
nw_pro$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                       in_lon_deg = nw_pro$lng, 
                                       in_lat_deg = nw_pro$lat)[[1]]




# write to RDS ------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
saveRDS(nw_pro, paste0('Nestwatch-pro-', DATE, '.rds'))

