####################
# 6 - query nesting phenology data (eBird breeding codes, MAPS, and Nestwatch)
#
# *eBird breeding codes
#   -0 if not observed in that survey at all (independent of breeding code)
#   -NA if observed but no breeding code recorded
#   -letter code if observed and breeding code recorded
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(dggridR)
library(doParallel)
library(foreach)


# import IAR data ---------------------------------------------------------

# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
# 
# IAR_out_dir <- 'IAR_output_2018-11-12'
# IAR_out_date <- substr(IAR_out_dir, start = 12, stop = 21)

#IAR_data <- readRDS(paste0('master_arrival_', IAR_out_date, '.rds'))



# import IAR species list -----------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/'))

species_list_i <- read.table(paste0('IAR_species_list.txt'), stringsAsFactors = FALSE)

#remove underscore and coerce to vector
species_list_i2 <- as.vector(apply(species_list_i, 2, function(x) gsub("_", " ", x)))



# eBird breeding codes ----------------------------------------------------


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


query_dir_path <- paste0('Processed/breeding_cat_query_', Sys.Date())

dir.create(query_dir_path)
setwd(query_dir_path)




# filter dataset notes ----------------------------------------------------------

#SQL FILTER BY:
# * only ebird data
# * species of interest (import from text file)
# * only complete checklists
# * duration minutes (6 min to 24 hours)
# * lat (> 26) and lon (-100 to -50)
# * year > 2001
# * day of year (< julian day 300)
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

#filter all unique event_ids that meet criteria - about 38 min to complete query
#same conditions for query as for arrival data

data <- DBI::dbGetQuery(cxn, paste0("
                                    SELECT DISTINCT ON (event_id) event_id, year, day, place_id, lat, lng, started, 
                                    radius,
                                    (event_json ->> 'ALL_SPECIES_REPORTED')::int AS all_species_reported,
                                    (event_json ->> 'DURATION_MINUTES')::int AS duration_minutes,
                                    count_json ->> 'OBSERVER_ID' AS observer_id,
                                    count_json ->> 'GLOBAL_UNIQUE_IDENTIFIER' AS global_unique_id,
                                    (event_json ->> 'NUMBER_OBSERVERS')::int AS number_observers,
                                    event_json ->> 'GROUP_IDENTIFIER' AS group_identifier
                                    FROM events
                                    INNER JOIN places USING (place_id)
                                    LEFT JOIN counts USING (event_id)                                    
                                    WHERE events.dataset_id = 'ebird'
                                    AND year > 2001
                                    AND day < 300
                                    AND lng BETWEEN -95 AND -50
                                    AND lat > 24
                                    AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
                                    AND LEFT(started, 2)::int < 18
                                    AND RADIUS < 100000;
                                    "))


#only values with unique group identifiers (or no group identifier)
data2 <- data[!duplicated(data[,'group_identifier'], 
                          incomparables = NA),]

rm(data)

#add jday and shr
cn_id <- grep('day', colnames(data2))
colnames(data2)[cn_id] <- 'jday'


#scaled effort hours
data2$shr <- as.vector(scale((data2$duration_minutes/60), scale = FALSE))



# bin lat/lon to hex grid and add to data ---------------------------------------------------------

# Construct geospatial hexagonal grid

hexgrid6 <- dggridR::dgconstruct(res = 6) 
data2$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                       in_lon_deg = data2$lng, 
                                       in_lat_deg = data2$lat)[[1]]




# query individual species, zero fill, and create RDS objects ----------------------------------

data2[species_list_i[,1]] <- NA

nsp <- NROW(species_list_i)

#run in parallel with 6 logical cores
doParallel::registerDoParallel(cores = 2)

tt <- proc.time()
foreach::foreach(i = 1:nsp) %dopar%
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
  
  temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, sci_name,
                                      count_json ->> 'BREEDING_BIRD_ATLAS_CODE' AS bba_code,
                                      count_json ->> 'BREEDING_BIRD_ATLAS_CATEGORY' AS bba_category
                                      FROM events
                                      INNER JOIN places USING (place_id)
                                      INNER JOIN counts USING (event_id)
                                      INNER JOIN taxa USING (taxon_id)
                                      WHERE events.dataset_id = 'ebird'
                                      AND year > 2001
                                      AND day < 300
                                      AND lng BETWEEN -95 AND -50
                                      AND lat > 24
                                      AND (sci_name IN ('", species_list_i2[i],"'));
                                      "))
  
  #cannot use global ID, as this is per species seen - there may be two species for a given survey (due to sub species)
  
  #find 'highest' breeding code and insert that for a given event_id
  '%ni%' <- Negate('%in%')
  #event ids of dups
  tn <- table(temp$event_id)
  ev2 <- as.numeric(names(tn[tn > 1]))
  #data without dups (neither dup)
  temp2 <- temp[which(temp[,'event_id'] %ni% ev2),]
  
  
  #take 'furthest along' observation from each survey
  
  #could use Nestwatch to get at more accurate species-specific numbers
  #FL = Recently fledged young - F
  #NY = Nest with young - Y
  #FY = Feeding young - Y
  #CS = Carrying fecal sac - Y
  #CF = Carrying food - Y
  #DD = Distraction display -Y
  #NE = Nest with egg - E
  #ON = Occupied nest - E
  #PE = Brood patch - E
  
  FLEDGE <- c('FL')
  YOUNG <- c('NY', 'FY', 'CS', 'CF', 'DD')
  EGG <- c('NE', 'ON', 'PE')
  
  #add dups in (only one) with corresponding br code
  if (length(ev2) > 0)
  {
    for (j in 1:length(ev2))
    {
      #j <- 5
      lt <- dplyr::filter(temp, event_id == ev2[j])
      
      #if there is a F, insert that
      #if there is a Y, insert that
      #if there is a E, insert that
      #if none, insert NA
      if (sum(lt$bba_code %in% FLEDGE, na.rm = TRUE) > 0)
      {
        tsind <- min(which(lt$bba_code %in% FLEDGE))
        temp2 <- rbind(temp2, lt[tsind,])
      } else {
        if (sum(lt$bba_code %in% YOUNG, na.rm = TRUE) > 0)
        {
          tsind <- min(which(lt$bba_code %in% YOUNG))
          temp2 <- rbind(temp2, lt[tsind,])
        } else {
          if (sum(lt$bba_code %in% EGG, na.rm = TRUE) > 0)
          {
            tsind <- min(which(lt$bba_code %in% EGG))
            temp2 <- rbind(temp2, lt[tsind,])
          } else { 
            tsind <- min(which(is.na(lt$bba_category)))
            temp2 <- rbind(temp2, lt[tsind,])
            
          }
        }
      }
    }
  }
  
  #NAs for codes we're not interested in
  NONE <- c('NB', 'CN', 'T', 'C', 'N', 'A', 'P', 'S', 'H', 'F', NA)
  temp2$bba_code[which(temp2$bba_code %in% NONE)] <- NA
  
  FLEDGE_uid <- temp2$event_id[which(temp2$bba_code %in% FLEDGE)]
  YOUNG_uid <- temp2$event_id[which(temp2$bba_code %in% YOUNG)]
  EGG_uid <- temp2$event_id[which(temp2$bba_code %in% EGG)]
  NONE_uid <- temp2$event_id[which(is.na(temp2$bba_code))]
  
  #check combined length matches number of rows
  #length(c(FLEDGE_uid, YOUNG_uid, EGG_uid, NONE_uid)) == NROW(temp2)
  
  #indices for data2 where the event ids match
  FLEDGE_ind <- which(data2$event_id %in% FLEDGE_uid)
  YOUNG_ind <- which(data2$event_id %in% YOUNG_uid)
  EGG_ind <- which(data2$event_id %in% EGG_uid)
  NONE_ind <- which(data2$event_id %in% NONE_uid)
  
  #0 if not observed in that survey at all (independent of breeding code)
  #NA if observed but no applicable breeding code recorded
  #all entries zero to start
  data2[, species_list_i[i,1]] <- 0
  data2[FLEDGE_ind, species_list_i[i,1]] <- 'F'
  data2[YOUNG_ind, species_list_i[i,1]] <- 'Y'
  data2[EGG_ind, species_list_i[i,1]] <- 'E'
  data2[NONE_ind, species_list_i[i,1]] <- NA
  
  sdata <- dplyr::select(data2, 
                         event_id, year, jday,
                         shr, cell, species_list_i[i,1])
  
  names(sdata)[6] <- "breeding_code"
  sdata['species'] <- species_list_i[i,1]
  
  saveRDS(sdata, file = paste0('ebird_NA_breeding_code_', species_list_i[i,1], '.rds'))
  
  DBI::dbDisconnect(cxn)
}
proc.time() - tt



# Find files that weren’t created -----------------------------------------

fls <- list.files()
nms <- c()
for (i in 1:length(fls))
{
  #i <- 1
  nms <- c(nms, substr(fls[i], 23, (nchar(fls[i]) - 4)))
}

#missed species
m_sp <- species_list_i[which(!species_list_i[,1] %in% nms),1]
#remove underscore
m_sp2 <- gsub("_", " ", m_sp)

if (length(m_sp2) > 0)
{
  #create those files
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
    
    temp <- DBI::dbGetQuery(cxn, paste0("SELECT event_id, sci_name,
                                        count_json ->> 'BREEDING_BIRD_ATLAS_CATEGORY' AS bba_category
                                        FROM events
                                        INNER JOIN places USING (place_id)
                                        INNER JOIN counts USING (event_id)
                                        INNER JOIN taxa USING (taxon_id)
                                        WHERE events.dataset_id = 'ebird'
                                        AND year > 2001
                                        AND day < 300
                                        AND lng BETWEEN -95 AND -50
                                        AND lat > 24
                                        AND (sci_name IN ('", m_sp2[i],"'));
                                        "))
    
    #find 'highest' breeding code and insert that for a given event_id
    '%ni%' <- Negate('%in%')
    #event ids of dups
    tn <- table(temp$event_id)
    ev2 <- as.numeric(names(tn[tn > 1]))
    #data without dups (neither dup)
    temp2 <- temp[which(temp[,'event_id'] %ni% ev2),]
    
    
    #take 'furthest along' observation from each survey
    
    #could use Nestwatch to get at more accurate species-specific numbers
    #FL = Recently fledged young - F
    #NY = Nest with young - Y
    #FY = Feeding young - Y
    #CS = Carrying fecal sac - Y
    #CF = Carrying food - Y
    #DD = Distraction display - Y
    #NE = Nest with egg - E
    #ON = Occupied nest - E
    #PE = Brood patch - E
    
    FLEDGE <- c('FL')
    YOUNG <- c('NY', 'FY', 'CS', 'CF', 'DD')
    EGG <- c('NE', 'ON', 'PE')
    
    #add dups in (only one) with corresponding br code
    if (length(ev2) > 0)
    {
      for (j in 1:length(ev2))
      {
        #j <- 5
        lt <- dplyr::filter(temp, event_id == ev2[j])
        
        #if there is a F, insert that
        #if there is a Y, insert that
        #if there is a E, insert that
        #if none, insert NA
        if (sum(lt$bba_code %in% FLEDGE, na.rm = TRUE) > 0)
        {
          tsind <- min(which(lt$bba_code %in% FLEDGE))
          temp2 <- rbind(temp2, lt[tsind,])
        } else {
          if (sum(lt$bba_code %in% YOUNG, na.rm = TRUE) > 0)
          {
            tsind <- min(which(lt$bba_code %in% YOUNG))
            temp2 <- rbind(temp2, lt[tsind,])
          } else {
            if (sum(lt$bba_code %in% EGG, na.rm = TRUE) > 0)
            {
              tsind <- min(which(lt$bba_code %in% EGG))
              temp2 <- rbind(temp2, lt[tsind,])
            } else { 
              tsind <- min(which(is.na(lt$bba_category)))
              temp2 <- rbind(temp2, lt[tsind,])
              
            }
          }
        }
      }
    }
    
    #NAs for codes we're not interested in
    NONE <- c('NB', 'CN', 'T', 'C', 'N', 'A', 'P', 'S', 'H', 'F', NA)
    temp2$bba_code[which(temp2$bba_code %in% NONE)] <- NA
    
    FLEDGE_uid <- temp2$event_id[which(temp2$bba_code %in% FLEDGE)]
    YOUNG_uid <- temp2$event_id[which(temp2$bba_code %in% YOUNG)]
    EGG_uid <- temp2$event_id[which(temp2$bba_code %in% EGG)]
    NONE_uid <- temp2$event_id[which(is.na(temp2$bba_code))]
    
    #check combined length matches number of rows
    #length(c(FLEDGE_uid, YOUNG_uid, EGG_uid, NONE_uid)) == NROW(temp2)
    
    #indices for data2 where the event ids match
    FLEDGE_ind <- which(data2$event_id %in% FLEDGE_uid)
    YOUNG_ind <- which(data2$event_id %in% YOUNG_uid)
    EGG_ind <- which(data2$event_id %in% EGG_uid)
    NONE_ind <- which(data2$event_id %in% NONE_uid)
    
    #0 if not observed in that survey at all (independent of breeding code)
    #NA if observed but no applicable breeding code recorded
    #all entries zero to start
    data2[, m_sp[i]] <- 0
    data2[FLEDGE_ind, m_sp[i]] <- 'F'
    data2[YOUNG_ind, m_sp[i]] <- 'Y'
    data2[EGG_ind, m_sp[i]] <- 'E'
    data2[NONE_ind, m_sp[i]] <- NA
    
    sdata <- dplyr::select(data2, 
                           event_id, year, jday,
                           shr, cell, m_sp[i])
    
    names(sdata)[6] <- "breeding_code"
    sdata['species'] <- m_sp[i]
    
    saveRDS(sdata, file = paste0('ebird_NA_breeding_code_', 
                                 m_sp[i], '.rds'))
    
    DBI::dbDisconnect(cxn)
  }
}


# copy script to query folder for records ---------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/ebird_br_code/6-query-nesting-data.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path, 
              '/6-query-nesting-data-', Sys.Date(), '.R'))


proc.time() - tt





# MAPS obs - DB -----------------------------------------------------------


dir <- '~/Google_Drive/R/'

run_date <- '2019-10-18'

# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)


# query DB ----------------------------------------------------------------


setwd(paste0(dir, 'wing_chord_changes/Data/'))

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg,
                      user = "cyoungflesh",
                      password = pass,
                      host = "35.221.16.125",
                      port = 5432,
                      dbname = "sightings")

#Entire North America - only days 100-300

maps_data <- DBI::dbGetQuery(cxn, paste0("SELECT lng, lat, year, day, common_name, sci_name, event_id, started, ended,
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
                                         WHERE events.dataset_id = 'maps'
                                         AND day BETWEEN 100 AND 300
                                         AND (count_json ->> 'C') in ('N', 'R', 'U')
                                         AND (count_json ->> 'N') in ('-', 'D', 'S', 'T', 'O', '+');
                                         "))

#C (Capture code): Only codes N,R,U are used for analyses according to MAPS docs
#N = newly banded bird
#R = recaptured bird
#U = unbanded bird

#N (Indicator to include in productivity and survivorship analyses): Okay to use these for analyses that do not require to control for effort
#0 = not caught at MAPS station
#T = outside normaps MAPS operation for that station
#S = caught within MAPS station boundary but not in a MAPS net
#D = date outside of MAPS periods
#- = record examined with current MAPS analytical procedure (taken to be standard capture methods)
#+ = record examined with preliminary MAPS analytical procedure

#Age:
#0 - unknown
#4 - local (young bird incapable of flight)
#2 - hatching-year bird
#1 - after hatching-year bird
#5 - second-year bird
#6 - after second-year bird
#7 - third-year bird
#8 - after third-year bird

#Sex:
#M - male
#F - female
#U - unknown
#X - unattempted

#Cloacal protuberance:
#0 - none
#1 - small
#2 - medium
#3 - large

#Brood patch:
#0 - none
#1 - smooth
#2 - vascularized
#3 - heavy
#4 - wrinkled
#5 - molting

#Fat content:
#0 - none
#1 - trace
#2 - light
#3 - half
#4 - full
#5 - bulging
#6 - greatly bulging
#7 - very excessive


# fill in true ages ------------------------------------------------------------

#true_age 1 means 1 year old (born in previous season)

#add col for true age (actual age, not 'year' of bird)
maps_data$true_age <- NA

#fill 0 for all birds banded in birth year
#2 = hatching year bird
#4 = local (young bird incapable of flight))
#5 = second year bird
#7 = third year bird - unreliable due to molting in second year (id as third year bird late second year)
maps_data$true_age[which(maps_data$age %in% c(2, 4))] <- 0

#ids for birds banded in birth year
ids <- unique(as.numeric(maps_data$band_id))

# enter in true age for all birds banded in birth year (or entered as second year or third year)
# create progress bar
pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
for (i in 1:length(ids))
{
  #i <- 61
  #filter for band_id (all years)
  temp <- dplyr::filter(maps_data, band_id == ids[i])
  t_yrs <- sort(as.numeric(unique(temp$year)))
  
  birth_year <- NA
  
  #define birth year based on assigned age class (band in birth year trumps other ages)
  if (sum(!is.na(temp$true_age)) > 0 & length(t_yrs) > 1)
  {
    birth_year <- t_yrs[1]
  } else {
    if (sum(temp$age == 5) > 0)
    {
      birth_year <- (t_yrs[1] - 1)
    }
  }
  
  #fill in age values for subsequent year if birth year was defined
  if (!is.na(birth_year))
  {
    for (j in 1:length(t_yrs))
    {
      #j <- 1
      #insert age
      maps_data[which(maps_data$band_id == ids[i] & maps_data$year == t_yrs[j]),
                'true_age'] <- (t_yrs[j] - birth_year)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

#one time fill with feather wear to avoid refilling age data
# qq <- readRDS('MAPS-age-filled.rds')
# qq$feather_wear <- maps_data$feather_wear
# setwd(paste0(dir, 'wing_chord_time_lat/Data/'))
# saveRDS(qq, 'MAPS-age-filled.rds')


#save RDS file
setwd(paste0(dir, 'wing_chord_changes/Data/'))
saveRDS(maps_data, paste0('MAPS-age-filled-', run_date, '.rds'))







# Nestwatch  ---------------------------------------------------------------

#to derive time period between egg -> hatch -> fledge
#query db

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
                                              (count_json ->> 'count_type') AS count_type
                                              FROM places
                                              JOIN events USING (place_id)
                                              JOIN counts USING (event_id)
                                              JOIN taxa USING (taxon_id)
                                              WHERE events.dataset_id = 'nestwatch'
                                              AND lng BETWEEN -95 AND -50
                                              AND lat > 24;
                                              "))



#only keep entries that have dates - convert to julian day
nestwatch_data$lay_j <- format(as.Date(nestwatch_data$lay, 
                                        format = '%d-%b-%y'), '%j')
nestwatch_data$hatch_j <- format(as.Date(nestwatch_data$hatch, 
                                       format = '%d-%b-%y'), '%j')
nestwatch_data$fledge_j <- format(as.Date(nestwatch_data$fledge, 
                                       format = '%d-%b-%y'), '%j')


nestwatch_data$sci_name <- gsub(' ', '_', nestwatch_data$sci_name)


out <- data.frame(species = rep(NA, length(species_list_i[,1])),
                  m_FH = NA,
                  m_FL = NA,
                  m_HL = NA,
                  n_FH = NA,
                  n_FL = NA,
                  n_HL = NA)
counter <- 1
#get metrics from nestwatch data
for (i in 1:length(species_list_i[,1]))
{
  #i <- 1
  
  #filter by species
  sp <- species_list_i[i,1]
  t_nw <- dplyr::filter(nestwatch_data, sci_name == sp)
  
  #fledge - hatch
  FH <- as.numeric(t_nw$fledge_j) - as.numeric(t_nw$hatch_j)
  FL <- as.numeric(t_nw$fledge_j) - as.numeric(t_nw$lay_j)
  HL <- as.numeric(t_nw$hatch_j) - as.numeric(t_nw$lay_j)

  out$species[counter] <- sp
  out$m_FH[counter] <- median(FH, na.rm = TRUE)
  out$m_FL[counter] <- median(FL, na.rm = TRUE)
  out$m_HL[counter] <- median(HL, na.rm = TRUE)
  out$n_FH[counter] <- sum(!is.na(FH))
  out$n_FL[counter] <- sum(!is.na(FL))
  out$n_HL[counter] <- sum(!is.na(HL))
  
  counter <- counter + 1
}

out2 <- dplyr::filter(out, n_FL > 10, n_HL > 10)
median(out2$m_FH)
median(out2$m_FL)
median(out2$m_HL)


