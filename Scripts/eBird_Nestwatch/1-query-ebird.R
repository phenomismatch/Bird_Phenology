####################
# 1 - query ebird data
#
# Filters eBird data from DB
# Zero fills
# Creates directory of processed data (rds file for each species) and copy of this 
# ... script (for reproducability) in /Data/Processed/eBird_query_<DATE>
#
# Runtime: 8-9 hours on 7 core machine - could be sped if run on on Xanadu?
####################


# #can access DB from command line with:
# psql "sslmode=disable dbname=sightings user=cyoungflesh hostaddr=35.221.16.125"



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



# create query dir and navigate there -------------------------------------------


query_dir_path <- paste0('Processed/eBird_query_', Sys.Date())

dir.create(query_dir_path)
setwd(query_dir_path)




# filter dataset notes ----------------------------------------------------------

#From Rafe:
#all ebird records must meet the following for inclusion in DB:
#*ALL_SPECIES_REPORTED = 1 (complete surveys)
#*APPROVED = 1
#*Must have date
#*Must have count
#*Lng between -95 and -50
#*Lat between 20 and 90
#*Sci Name on approved ~120 passerine list


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

#*get info for each unique event
#*create species columns with NA
#*in a loop, query each species individually and fill obs with 1s, no obs with 0s
#*merge with hex cells
#*create rds objects for each species


#filter all unique event_ids that meet criteria - about 38 min to complete query
data <- DBI::dbGetQuery(cxn, paste0("
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
                                    JOIN taxa USING (taxon_id)
                                    WHERE events.dataset_id = 'ebird'
                                    AND year > 2001
                                    AND day < 200
                                    AND lng BETWEEN -95 AND -50
                                    AND lat > 26
                                    AND (sci_name IN (", SL,"))
                                    AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
                                    AND LEFT(started, 2)::int < 18
                                    AND RADIUS < 100000;
                                    "))


#only values with unique group identifiers (or no group identifier)
data2 <- data[!duplicated(data[,'group_identifier'], 
                          incomparables = NA),]

rm(data)



# add sjday, sjday^2, sjday^3 and shr ---------------------------------------------------

#julian day
data2$jday <- as.vector(data2$day)

#scaled effort hours
data2$shr <- as.vector(scale((data2$duration_minutes/60), scale = FALSE))



# bin lat/lon to hex grid and add to data ---------------------------------------------------------

# Construct geospatial hexagonal grid

hexgrid6 <- dggridR::dgconstruct(res = 6) 
data2$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                        in_lon_deg = data2$lng, 
                                        in_lat_deg = data2$lat)[[1]]




# query individual species, zero fill, and create RDS objects ----------------------------------


#create species columns
data2[species_list_i[,1]] <- NA

nsp <- NROW(species_list_i)

#run in parallel with 7 logical cores
doParallel::registerDoParallel(cores = 4)

tt <- proc.time()
foreach::foreach(i = 1:nsp) %dopar%
{
  #i <- 89
  print(i)
  
  pg <- DBI::dbDriver("PostgreSQL")
  
  cxn <- DBI::dbConnect(pg, 
                        user = "cyoungflesh", 
                        password = pass, 
                        host = "35.221.16.125", 
                        port = 5432, 
                        dbname = "sightings")
  
  temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id
                                      FROM events
                                      JOIN places USING (place_id)
                                      JOIN counts USING (event_id)
                                      JOIN taxa USING (taxon_id)
                                      WHERE events.dataset_id = 'ebird'
                                      AND year > 2001
                                      AND day < 200
                                      AND lng BETWEEN -95 AND -50
                                      AND lat > 26
                                      AND (sci_name IN ('", species_list_i2[i],"'));
                                      "))
  
  #indices to fill with 1s (observations of species made for these events)
  ind <- which(data2$event_id %in% temp$event_id)
  
  #indices to fill with 0s (observations of species not made for these events)
  n_ind <- (1:NROW(data2))[-ind]
  
  #fill observersations for species i with 1s
  data2[ind, species_list_i[i,1]] <- 1
  
  #fill no observations for species i with 0s
  data2[n_ind, species_list_i[i,1]] <- 0
  
  sdata <- dplyr::select(data2, 
                         year, day, jday,
                         shr, cell, species_list_i[i,1])
  
  names(sdata)[8] <- "detect"
  sdata['species'] <- species_list_i[i,1]
  
  saveRDS(sdata, file = paste0('ebird_NA_phen_proc_', species_list_i[i,1], '.rds'))
  DBI::dbDisconnect(cxn)
}
proc.time() - tt



# Find files that weren’t created -----------------------------------------

fls <- list.files()
nms <- c()
for (i in 1:length(fls))
{
  #i <- 1
  nms <- c(nms, substr(fls[i], 20, (nchar(fls[i]) - 4)))
}

#missed species
m_sp <- species_list_i[which(!species_list_i[,1] %in% nms),1]
#remove underscore
m_sp2 <- gsub("_", " ", m_sp)

#create those files
if (length(m_sp2) > 0)
{
  for (i in 1:length(m_sp2))
  {
    #i <- 1
    print(i)
    
    pg <- DBI::dbDriver("PostgreSQL")
    
    cxn <- DBI::dbConnect(pg, 
                          user = "cyoungflesh", 
                          password = pass, 
                          host = "35.221.16.125", 
                          port = 5432, 
                          dbname = "sightings")
    
    temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id
                                        FROM events
                                        JOIN places USING (place_id)
                                        JOIN counts USING (event_id)
                                        JOIN taxa USING (taxon_id)
                                        WHERE events.dataset_id = 'ebird'
                                        AND year > 2001
                                        AND day < 200
                                        AND lng BETWEEN -95 AND -50
                                        AND lat > 26
                                        AND (sci_name IN ('", m_sp2[i],"'));
                                        "))
    
    #indices to fill with 1s (observations of species made for these events)
    ind <- which(data2$event_id %in% temp$event_id)
    
    #indices to fill with 0s (observations of species not made for these events)
    n_ind <- (1:NROW(data2))[-ind]
    
    #fill observersations for species i with 1s
    data2[ind, m_sp2[i]] <- 1
    
    #fill no observations for species i with 0s
    data2[n_ind, m_sp2[i]] <- 0
    
    sdata <- dplyr::select(data2, 
                           year, day, jday,
                           shr, cell, m_sp2[i])
    
    names(sdata)[8] <- "detect"
    sdata['species'] <- m_sp2[i]
    
    saveRDS(sdata, file = paste0('ebird_NA_phen_proc_', m_sp[i], '.rds'))
  }
}



# copy script to query folder for records ---------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/1-query-ebird.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path, '/1-query-ebird-', Sys.Date(), '.R'))


proc.time() - tt
