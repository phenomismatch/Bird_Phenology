####################
# 6 - query other nesting phenology data
#
# *eBird breeding codes - approx 7 hours
#   -0 if not observed in that survey at all (independent of breeding code)
#   -NA if observed but no breeding code recorded
#   -letter code if observed and breeding code recorded
# *MAPS observations
# *Nestwatch
# *MAPS data
####################


# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'

hm_date <- '2019-02-02'

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


#combine species names into a single string with quotes
SL <- paste0("'", species_list_i2, "'", collapse = ", ")





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
                                    (event_json ->> 'NUMBER_OBSERVERS')::int AS number_observers,
                                    event_json ->> 'GROUP_IDENTIFIER' AS group_identifier
                                    FROM places
                                    JOIN events USING (place_id)
                                    JOIN counts USING (event_id)
                                    JOIN taxa USING (taxon_id)
                                    WHERE events.dataset_id = 'ebird'
                                    AND year > 2001
                                    AND day < 300
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


#calculate polynomial then center data
#scaled julian day, scaled julian day^2, and scaled julian day^3
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

data2[species_list_i[,1]] <- NA

nsp <- NROW(species_list_i)

#run in parallel with 6 logical cores
doParallel::registerDoParallel(cores = 4)

tt <- proc.time()
foreach::foreach(i = 1:nsp) %dopar%
{
  #i <- 100
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
                                      JOIN places USING (place_id)
                                      JOIN counts USING (event_id)
                                      JOIN taxa USING (taxon_id)
                                      WHERE events.dataset_id = 'ebird'
                                      AND year > 2001
                                      AND day < 300
                                      AND lng BETWEEN -100 AND -50
                                      AND lat > 26
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
  
  #add dups in (only one) with corresponding br code
  if (length(ev2) > 0)
  {
    for (j in 1:length(ev2))
    {
      #i <- 2
      lt <- dplyr::filter(temp, event_id == ev2[j])
      
      #if there is a C4, insert that
      #if there is a C3, insert that
      #if neither, insert NA
      if (sum(lt$bba_category == 'C4', na.rm = TRUE) > 0)
      {
        tsind <- min(which(lt$bba_category == 'C4'))
        temp2 <- rbind(temp2, lt[tsind,])
      } else {
        if (sum(lt$bba_category == 'C3', na.rm = TRUE) > 0)
        {
          tsind <- min(which(lt$bba_category == 'C3'))
          temp2 <- rbind(temp2, lt[tsind,])
        } else {
          tsind <- min(which(is.na(lt$bba_category)))
          temp2 <- rbind(temp2, lt[tsind,])
        }
      }
    }
  }
  
  #NAs for C1 and C2 (not/possible breeding) as we're not interested in these metrics
  temp2$bba_category[which(temp2$bba_category == 'C1' | temp2$bba_category == 'C2')] <- NA
  
  #just c3/4
  c3_uid <- temp2$event_id[which(temp2$bba_category == 'C3')]
  c4_uid <- temp2$event_id[which(temp2$bba_category == 'C4')]
  z_uid <- temp2$event_id[which(is.na(temp2$bba_category))]
  #check combined length matches number of rows
  #length(c(c3_uid, c4_uid, z_uid)) == NROW(temp2)
  
  #indices for data2 where the event ids match
  c3_ind <- which(data2$event_id %in% c3_uid)
  c4_ind <- which(data2$event_id %in% c4_uid)
  z_ind <- which(data2$event_id %in% z_uid)
  
  #0 if not observed in that survey at all (independent of breeding code)
  #NA if observed but no breeding code recorded
  #letter code if observed and breeding code recorded
  
  data2[c3_ind, species_list_i[i,1]] <- 'C3'
  data2[c4_ind, species_list_i[i,1]] <- 'C4'
  data2[z_ind, species_list_i[i,1]] <- 0
  
  sdata <- dplyr::select(data2, 
                         year, day, cell, jday,
                         shr, species_list_i[i,1])
  
  names(sdata)[8] <- "bba_category"
  sdata['species'] <- species_list_i[i,1]
  
  saveRDS(sdata, file = paste0('ebird_NA_breeding_cat_', species_list_i[i,1], '.rds'))
  
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
                                        JOIN places USING (place_id)
                                        JOIN counts USING (event_id)
                                        JOIN taxa USING (taxon_id)
                                        WHERE events.dataset_id = 'ebird'
                                        AND year > 2001
                                        AND day < 300
                                        AND lng BETWEEN -95 AND -50
                                        AND lat > 26
                                        AND (sci_name IN ('", m_sp2[i],"'));
                                        "))
    
    #find 'highest' breeding code and insert that for a given event_id
    '%ni%' <- Negate('%in%')
    #event ids of dups
    tn <- table(temp$event_id)
    ev2 <- as.numeric(names(tn[tn > 1]))
    #data without dups (neither dup)
    temp2 <- temp[which(temp[,'event_id'] %ni% ev2),]
    
    if (length(ev2) > 0)
    {
      #add dups in (only one) with corresponding br code
      for (j in 1:length(ev2))
      {
        #i <- 2
        lt <- dplyr::filter(temp, event_id == ev2[j])
        
        #if there is a C4, insert that
        #if there is a C3, insert that
        #if neither, insert NA
        if (sum(lt$bba_category == 'C4', na.rm = TRUE) > 0)
        {
          tsind <- min(which(lt$bba_category == 'C4'))
          temp2 <- rbind(temp2, lt[tsind,])
        } else {
          if (sum(lt$bba_category == 'C3', na.rm = TRUE) > 0)
          {
            tsind <- min(which(lt$bba_category == 'C3'))
            temp2 <- rbind(temp2, lt[tsind,])
          } else {
            tsind <- min(which(is.na(lt$bba_category)))
            temp2 <- rbind(temp2, lt[tsind,])
          }
        }
      }
    }
    
    #NAs for C1 and C2 (not/possible breeding) as we're not interested in these metrics
    temp2$bba_category[which(temp2$bba_category == 'C1' | temp2$bba_category == 'C2')] <- NA
    
    #just c3/4
    c3_uid <- temp2$event_id[which(temp2$bba_category == 'C3')]
    c4_uid <- temp2$event_id[which(temp2$bba_category == 'C4')]
    z_uid <- temp2$event_id[which(is.na(temp2$bba_category))]
    #check combined length matches number of rows
    #length(c(c3_uid, c4_uid, z_uid)) == NROW(temp2)
    
    #indices for data2 where the event ids match
    c3_ind <- which(data2$event_id %in% c3_uid)
    c4_ind <- which(data2$event_id %in% c4_uid)
    z_ind <- which(data2$event_id %in% z_uid)
    
    #0 if not observed in that survey at all (independent of breeding code)
    #NA if observed but no breeding code recorded
    #letter code if observed and breeding code recorded
    
    data2[c3_ind, m_sp[i]] <- 'C3'
    data2[c4_ind, m_sp[i]] <- 'C4'
    data2[z_ind, m_sp[i]] <- 0
    
    sdata <- dplyr::select(data2, 
                           year, day, cell, jday,
                           shr, m_sp[i])
    
    names(sdata)[8] <- "bba_category"
    sdata['species'] <- m_sp[i]
    
    saveRDS(sdata, file = paste0('ebird_NA_breeding_cat_', m_sp[i], '.rds'))
    DBI::dbDisconnect(cxn)
  }
}


# copy script to query folder for records ---------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/6-query-nesting-data.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path, 
              '/6-query-nesting-data-', Sys.Date(), '.R'))


proc.time() - tt





# counts for breeding codes -----------------------------------------------

# #IAR input data (to get relevant cells)
# DATE_MA <- '2018-11-12'
# 
# setwd(paste0(dir, 'Bird_phenology/Data/Processed/IAR_input_', DATE_MA))
# df_master <- readRDS(paste0('IAR_input-', DATE_MA, '.rds'))
# 
# nsp <- NROW(species_list_i)
# 
# DATE_BC <- '2019-01-28'
# 
# output_df <- data.frame()
# for (i in 1:nsp)
# {
#   #i <- 17
#   #read in ebird breeding code data
#   
#   setwd(paste0(dir, 'Bird_phenology/Data/Processed/breeding_cat_query_', DATE_BC))
#   temp_bc <- readRDS(paste0('ebird_NA_breeding_cat_', species_list_i[i,1], '.rds'))
#   temp_master <- dplyr::filter(df_master, species == species_list_i[i,1])
#   
#   #only cells that are in IAR input
#   kp_cells <- unique(temp_master$cell)
#   temp_bc_f <- dplyr::filter(temp_bc, cell %in% kp_cells)
#   
#   #probable/confirmed
#   t_C34 <- dplyr::filter(temp_bc_f,
#                          bba_category == 'C3' |
#                            bba_category == 'C4')
#   
#   #how many C3/C4 obs for a given year
#   tt <- plyr::count(t_C34, vars = c('year'))
#   #how many C3/C4 obs for a given year/cell
#   cy <- plyr::count(t_C34, vars = c('cell', 'year'))
#   #filter by > 20 obs in cell
#   gr20 <- dplyr::filter(cy, freq > 20)
#   #how cells in each year have # C3/C4 obs > 20
#   nc_gr20 <- plyr::count(gr20[,1:2], vars = c('year'))
#   colnames(nc_gr20)[2] <- 'num_usable_cells'
#   #merge with number of C3/C4 obs
#   mrg <- dplyr::left_join(tt, nc_gr20, by = 'year')
#   colnames(mrg)[2] <- 'num_br_obs'
#   #merge with years
#   nn <- data.frame(year = 2002:2017)
#   mrg2 <- dplyr::left_join(nn, mrg, by = 'year')
#   #insert zeros for NA vals
#   to.z.use <- which(is.na(mrg2[,'num_usable_cells']))
#   mrg2[to.z.use, 'num_usable_cells'] <- 0
#   to.z.obs <- which(is.na(mrg2[,'num_br_obs']))
#   mrg2[to.z.obs, 'num_br_obs'] <- 0
#   
#   t_out <- data.frame(species = species_list_i[i,1], mrg2)
#   output_df <- rbind(output_df, t_out)
# }
# 
# 
# num_years <- aggregate(num_usable_cells ~ species, data = output_df, FUN = function(x) sum(x > 0))
# colnames(num_years)[2] <- 'num_years'
# 
# summary_output_df <- aggregate(num_usable_cells ~ species, data = output_df, FUN = sum)
# colnames(summary_output_df)[2] <- 'num_cell_years'
# summary_output_df$num_years <- num_years[,'num_years']
# 
# 
# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
# write.csv(output_df, paste0('br_code_data_avail-', 
#                             DATE_BC, '.csv'), row.names = FALSE)
# write.csv(summary_output_df, paste0('summary_br_code_data_avail-', 
#                                     DATE_BC, '.csv'), row.names = FALSE)

# #read in data
# setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
# summary_output_df <- read.csv('summary_br_code_data_avail.csv')
# sum(summary_output_df[,2] > 20)
# hist(summary_output_df[,2], col = 'grey',
#      xlab = 'Number cell/years of data',
#      main = 'eBird breeding code availability')
# 
# output_df <- read.csv('br_code_data_avail.csv')





# MAPS obs ----------------------------------------------------------------


setwd(paste0(dir, 'Bird_phenology/Data/MAPS_Obs/Export_2018_11'))

MAPS_obs <- read.csv('1106Band18.csv', stringsAsFactors = FALSE)

#Capture_code: only codes N,R,U are used for analyses
#Species: 4-letter codes
#Age: 
#0 - unknown
#4 - local (young bird incapable of flight)
#2 - hatching-year bird
#1 - after hartching-year bird
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

#get grid cell of each station by lat/lon
setwd(paste0(dir, 'Bird_phenology/Data/MAPS_Obs/CntrlStations'))
MAPS_stations <- read.csv('STATIONS.csv', skipNul = TRUE, stringsAsFactors = FALSE)

MAPS_mrg <- dplyr::left_join(MAPS_obs, MAPS_stations, by = 'STA')

#check to make sure STATION names are in same
#sum(MAPS_mrg$STATION.x != MAPS_mrg$STATION.y)

#LAT/LON PROBLEM STATIONS: 
#DFME (Harrison, OH) - NOT IN OBS
#SSPT (Ashford, WA) - NOT IN OBS
#HDQR (Weldon, CA) - NOT IN OBS
#MAMA (Weldon, CA) - NOT IN OBS
#PLMR (Weldon, CA) - NOT IN OBS
#SUNY (Old Westbury, NY) - NOT IN OBS
#RKCK (REMOVE) - NOT IN OBS
#EASA (New Waverly, TX) - (30.54, -95.48)
#ANWR (Kotz Springs, LA) - (30.54, -91.75)
#REGR (Columbus, AR) - (33.78, -93.82)
#CEFB (Clemson, SC) - (34.68, -82.81)
#HRAE (REMOVE) - NOT IN OBS

#These stations do not appear in the data
# #insert lat/lons for these stations
# st_rep <- data.frame(STATION = c('EASA', 'ANWR', 'REGR', 'CEFB'), 
#                      LAT = c(30.54, 30.54, 33.78, 34.68),
#                      LNG = c(-95.48, -91.75, -93.82, -82.81))
# 
# 
# for (i in 1:NROW(st_rep))
# {
#   #i <- 4
#   tt <- which(MAPS_mrg$STATION.x == st_rep$STATION[i])
#   MAPS_mrg[tt,'DECLAT'] <- st_rep$LAT[i]
#   MAPS_mrg[tt,'DECLNG'] <- st_rep$LNG[i]
# }


hexgrid6 <- dggridR::dgconstruct(res = 6) 
MAPS_mrg$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                          in_lon_deg = MAPS_mrg$DECLNG, 
                                          in_lat_deg = MAPS_mrg$DECLAT)[[1]]

#merge MAPS data with species codes
setwd(paste0(dir, 'Bird_Phenology/Data/MAPS_Obs'))
sp_codes <- read.csv('sp_4_letter_codes.csv', stringsAsFactors = FALSE)
MAPS_mrg2 <- dplyr::left_join(MAPS_mrg, sp_codes, by = 'SPEC')

#split date into separate cols
dsplit <- strsplit(MAPS_mrg2$DATE, '/')
MAPS_mrg2$month <- unlist(dsplit)[3*(1:length(MAPS_mrg2$DATE)) - 2]
MAPS_mrg2$day <- unlist(dsplit)[3*(1:length(MAPS_mrg2$DATE)) - 1]
MAPS_mrg2$year <- unlist(dsplit)[3*(1:length(MAPS_mrg2$DATE))]

post_2000 <- which(as.numeric(MAPS_mrg2$year) < 70)
pre_2000 <- which(as.numeric(MAPS_mrg2$year) > 70)

MAPS_mrg2$year[post_2000] <- paste0('20', MAPS_mrg2$year[post_2000])
MAPS_mrg2$year[pre_2000] <- paste0('19', MAPS_mrg2$year[pre_2000])

MAPS_mrg2$jday <- as.numeric(format(as.Date(paste0(MAPS_mrg2$year, '-', 
                                                  MAPS_mrg2$month, '-', MAPS_mrg2$day)), '%j'))

MAPS_mrg3 <- dplyr::select(MAPS_mrg2, STA, LOC.x, cell, DATE, year, month, day, jday, C, BAND, 
                           SPEC, SCINAME, AGE, SEX, CP, BP, F, WEIGHT)

colnames(MAPS_mrg3) <- c('STA', 'Location_code', 'Cell', 'Date', 'Year', 'Month', 'Day', 
                         'Jday', 'Capture_code', 'Band_id', 'Species', 'Sci_name', 'Age_code', 'Sex', 
                         'Cloacal_prot', 'Brood_patch', 'Fat_content', 'Weight')

#filter by species and capture code (only codes N,R,U are used for analyses)
MAPS_mrg4 <- dplyr::filter(MAPS_mrg3, Sci_name %in% species_list_i2 & 
                             Capture_code %in% c('N', 'R', 'U') &
                             Year > 2001)

#Age split up into AFTER... codes - these may or may not have been advanced in database for recaptures
#What characterizes breeding?
#Does effort need to be controlled for?

#2 - zero year
#4 - zero year
#5 - second year
#7 - third year
#...

#add col for true age (actual age, not 'year' of bird)
MAPS_mrg5 <- MAPS_mrg4
MAPS_mrg5$True_age <- NA

#fill 0 for all birds banded in birth year
MAPS_mrg5$True_age[which(MAPS_mrg5$Age_code %in% c(2,4))] <- 0 

#ids for birds banded in birth year
ids <- unique(as.numeric(MAPS_mrg5$Band_id))


#enter in true age for all birds banded in birth year
# create progress bar
pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
for (i in 1:length(ids))
{
  #i <- 61
  #filter for band_id (all years)
  temp <- dplyr::filter(MAPS_mrg5, Band_id == ids[i])
  t_yrs <- sort(as.numeric(unique(temp$Year)))
  
  birth_year <- NA
  
  #define birth year based on assigned age class (band in birth year trumps other ages)
  if (sum(!is.na(temp$True_age)) > 0 & length(t_yrs) > 1)
  {
    birth_year <- t_yrs[1]
  } else {
    if (sum(temp$Age_code == 5) > 0 & length(t_yrs) > 1)
    {
      birth_year <- (t_yrs[1] - 1)
    } else {
      if (sum(temp$Age_code == 7) > 0 & length(t_yrs) > 1)
      {
        birth_year <- (t_yrs[1] - 2)
      }
    }
  }
  
  #fill in age values for subsequent year if birth year was defined
  if (!is.na(birth_year))
  {
    for (j in 2:length(t_yrs))
    {
      #j <- 2
      #insert age
      MAPS_mrg5[which(MAPS_mrg5$Band_id == ids[i] & MAPS_mrg5$Year == t_yrs[j]), 
                'True_age'] <- (t_yrs[j] - birth_year)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

#proportion of captures that have a True_age
sum(!is.na(MAPS_mrg5$True_age))/NROW(MAPS_mrg5)


# setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# saveRDS(MAPS_mrg5, 'MAPS_age_filled.rds')
# MAPS_mrg5 <- readRDS('MAPS_age_filled.rds')


#--------------------------------#
#logistic model for presence of brood patch

#just hatch year birds
MAPS_f <- dplyr::filter(MAPS_mrg5, True_age == 0)
plyr::count(MAPS_f, 'Sci_name')

oc_0 <- filter(MAPS_f, Sci_name == 'Oreothlypis celata')
plyr::count(oc_0, c('Cell', 'Year'))
tt <- dplyr::filter(oc_0, Cell == 444, Year == 2011)

plot(density(tt$Jday))


#--------------------------------#


#--------------------------------#
#logistic model for presence of brood patch

MAPS_age <- dplyr::filter(MAPS_mrg5, True_age > 0)
plyr::count(MAPS_age, 'Sci_name')

C_ustulatus <- filter(MAPS_age, Sci_name == 'Catharus ustulatus')

#breeding defined as:
#BP 2-4
#CP 1-3

C_ustulatus <- dplyr::filter(MAPS_age, Sci_name == 'Catharus ustulatus')

cu_f <- dplyr::filter(C_ustulatus, Sex == 'F')
cu_m <- dplyr::filter(C_ustulatus, Sex == 'M')

#add breeder col
C_ustulatus$BR <- NA
cu_f$BR <- NA

#male female
plot(density(C_ustulatus[which(C_ustulatus$Sex == 'M'),]$Jday))
lines(density(C_ustulatus[which(C_ustulatus$Sex == 'F'),]$Jday), col = 'red')

#cloacal protuberance probably not a good metric
plot(density(cu_m[which(cu_m$Cloacal_prot <= 2),]$Jday))
lines(density(cu_m[which(cu_m$Cloacal_prot > 2),]$Jday), col = 'red')

#brood patch maybe better
plot(density(cu_f[which(cu_f$Brood_patch <= 2 | cu_f$Brood_patch == 5),]$Jday))
lines(density(cu_f[which(cu_f$Brood_patch > 2 & cu_f$Brood_patch < 5),]$Jday), 
      col = 'red')

cu_f[which(cu_f$Brood_patch > 2 & cu_f$Brood_patch < 5),]$BR <- 1
cu_f[which(cu_f$Brood_patch <= 2 | cu_f$Brood_patch == 5),]$BR <- 0

#DATA <- dplyr::filter(cu_f, Jday < DAY, Jday > 50)
DATA <- cu_f

DATA$sjday <- DATA$Jday
DATA$sjday2 <- DATA$Jday^2
DATA$sjday3 <- DATA$Jday^3

plyr::count(cu_f, c('Cell', 'Year'))


library(rstanarm)
ITER <- 1000
#ITER <- 10
CHAINS <- 3

fit2 <- rstanarm::stan_glm(BR ~ sjday + sjday2 + sjday3,
                           data = DATA,
                           family = binomial(link = "logit"),
                           algorithm = 'sampling',
                           iter = ITER,
                           chains = CHAINS,
                           cores = CHAINS)

summary(fit2)

predictDays <- range(DATA$sjday)[1]:range(DATA$sjday)[2]
predictDays2 <- predictDays^2
predictDays3 <- predictDays^3

newdata <- data.frame(sjday = predictDays, sjday2 = predictDays2, sjday3 = predictDays3)

dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))

for (L in 1:((ITER/2)*CHAINS))
{
  #L <- 1
  rowL <- as.vector(dfit[L,])
  halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
}


mn_dfit <- apply(dfit, 2, mean)
LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
mn_hm <- mean(halfmax_fit)
LCI_hm <- quantile(halfmax_fit, probs = 0.025)
UCI_hm <- quantile(halfmax_fit, probs = 0.975)

plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
     ylim = c(0, max(UCI_dfit)),
     xlab = 'Julian Day', ylab = 'Detection Probability')
lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
lines(predictDays, mn_dfit, lwd = 2)
DATA$BR_pl <- NA
DATA[which(DATA$BR == 1),]$BR_pl <- max(UCI_dfit)
DATA[which(DATA$BR == 0),]$BR_pl <- 0
points(DATA$Jday, DATA$BR_pl, col = rgb(0,0,0,0.25))
abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
mean(halfmax_fit)
sd(halfmax_fit)



#--------------------------------#


# 
# jfit <- lm(MAPS_age$Jday ~ MAPS_age$True_age)
# summary(jfit)
# 
# plot(MAPS_age$True_age, MAPS_age$Jday, pch = '.', col = rgb(0,0,0, 0.5))
# abline(jfit, col = 'red')
# a1 <- dplyr::filter(MAPS_age, True_age == 1)
# a2 <- dplyr::filter(MAPS_age, True_age == 2)
# a3 <- dplyr::filter(MAPS_age, True_age == 3)
# a4 <- dplyr::filter(MAPS_age, True_age == 4)
# a5 <- dplyr::filter(MAPS_age, True_age == 5)
# plot(density(a1$Jday))
# lines(density(a2$Jday))
# lines(density(a3$Jday))
# lines(density(a4$Jday))
# lines(density(a5$Jday))




#*breeding timing
#Each species/cell/year:
#maybe each capture different row, 1 for brood patch/CP, 0 if not -> fit a logistic to this for each cell/year?
#y[i] ~ bern(p[i])
#logit(p[i]) = alpha + beta * DAY[i]

#*phenology and age
#could fit logistic regression with age as random effect -> each age class would have different beta param
#for each cell/year: predict start of breeding for each age class (age 1, age 2, age 3+)
#calculate derived qty - difference in halfmax estimates
#y[i] ~ bern(p[i])
#logit(p[i]) = alpha[id[i]] + beta[id[i]] * DAY[i]

#*difference in arrival between males and females - how this maps onto phlogeny/traits


#*fat content and age - Fat[i] ~ alpha[id[i]] + beta1[id[i]] * Age[i] + beta2[id[i]] * Day[i] + beta3[id[i]] * Day[i]^2
plot(MAPS_age$True_age, MAPS_age$Fat_content, col = rgb(0,0,0, 0.5))
plot(MAPS_age$Jday, MAPS_age$Fat_content, col = rgb(0,0,0, 0.5))
summary(lm(Fat_content ~ True_age + poly(Jday, 2, raw = TRUE), data = MAPS_age))

#*Weight and age - Weight[i] ~ alpha[id[i]] + beta1[id[i]] * Age[i] + beta2[id[i]] * Day[i] + beta3[id[i]] * Day[i]^2
plot(MAPS_age$True_age, MAPS_age$Weight, col = rgb(0,0,0, 0.5))
plot(MAPS_age$Jday, MAPS_age$Weight, col = rgb(0,0,0, 0.5))
summary(lm(Weight ~ True_age + poly(Jday, 2, raw = TRUE), data = MAPS_age))



#filter by species, get breeding date (breeding period in this case) for each year/cell
na_reps <- rep(NA, NROW(unique(MAPS_mrg2[c('YR', 'cell', 'SCINAME')])))

MAPS_out <- data.frame(YR = na_reps,
                       cell = na_reps,
                       SCINAME = na_reps,
                       midpoint = na_reps,
                       l_bounds = na_reps,
                       u_bounds = na_reps,
                       n_stations = na_reps)

counter <- 1
for (i in 1:length(species_list_i2))
{
  #i <- 19
  temp <- dplyr::filter(MAPS_mrg2, SCINAME == species_list_i2[i])
  if (NROW(temp) > 0)
  {
    
    t_cell <- unique(temp$cell)
    
    for (j in 1:length(t_cell))
    {
      #j <- 4
      temp2 <- dplyr::filter(temp, cell == t_cell[j])
      t_yr <- unique(temp2$YR)
      
      for (k in 1:length(t_yr))
      {
        
        print(paste0('species: ', species_list_i2[i], ', ',
                     'cell: ', t_cell[j], ', ',
                     'year: ', t_yr[k]))
        
        #k <- 12
        temp3 <- dplyr::filter(temp2, YR == t_yr[k])
        
        #C = confirmed breeder
        #P = probably breeder
        #O = observed
        #- = not observed
        
        #input is 'MM-DD' (05-01)
        #returns julian day
        dfun <- function(input)
        {
          t_date_lb <- paste0(t_yr[k], '-', input)
          t_j_lb <- format(as.Date(t_date_lb), '%j')
          return(t_j_lb)
        }
        
        PS1 <- dfun('05-01')
        PS2 <- dfun('05-11')
        PS3 <- dfun('05-21')
        PS4 <- dfun('05-31')
        PS5 <- dfun('06-10')
        PS6 <- dfun('06-20')
        PS7 <- dfun('06-30')
        PS8 <- dfun('07-10')
        PS9 <- dfun('07-20')
        PS10 <- dfun('07-30')
        PS11 <- dfun('08-09')
        PS12 <- dfun('08-18')
        
        periods <- c(PS1, PS2, PS3, PS4, PS5, PS6, PS7, PS8, PS9, PS10, PS11, PS12)
        
        #from first observaton observed breeding in hex cell:
        #input julian day if breeding
        #input 0 if observed in the season but not breeding
        #input NA if not observed at all
        
        cind <- grep('PS', colnames(temp3))
        
        #confirmed or possbile breeding for each station
        cp_ind <- c()
        ob_ind <- c()
        no_ind <- c()
        for (m in 1:NROW(temp3))
        {
          #m <- 1
          cp_ind <- c(cp_ind, which(temp3[m,cind] == 'C' | temp3[m,cind] == 'P'))
          ob_ind <- c(ob_ind, which(temp3[m,cind] == 'O'))
          no_ind <- c(no_ind, which(temp3[m,cind] == '-'))
        }
        
        if (length(cp_ind) > 0)
        {
          #lower bounds of period observed
          first_br_obs_lb <- as.numeric(periods[min(cp_ind)])
          #upper bounds of period observed
          first_br_obs_ub <- first_br_obs_lb + 9
          midpoint_br <- mean(c(first_br_obs_lb, first_br_obs_ub))
        } else {
          if (length(ob_ind) > 0)
          {
            first_br_obs_lb <- 0
            first_br_obs_ub <- 0
            midpoint_br <- 0
          } else {
            first_br_obs_lb <- NA
            first_br_obs_ub <- NA
            midpoint_br <- NA
          }
        }
        
        #column for midpoint breeding date
        #column for lower bounds breeding date
        #column for upper bounds breeding date
        #column for cell
        #column for year
        #column for species SCI NAME (with underscore)
        
        t_info <- dplyr::select(temp3, YR, cell, SCINAME)[1,]
        
        t_info$SCINAME <- gsub(' ', '_', t_info$SCINAME)
        
        t_info$midpoint <- midpoint_br
        t_info$l_bounds <- first_br_obs_lb
        t_info$u_bounds <- first_br_obs_ub
        #n_stations is number of stations that recorded probably or confirmed breeding
        t_info$n_stations <- length(cp_ind)
        
        MAPS_out[counter,] <- t_info
        counter <- counter + 1
      }
    }
  }
}



#remove NA padding (find first year row to have NA and subtract one from that index)
fin_ind <- min(which(is.na(MAPS_out$YR)))
MAPS_out2 <- MAPS_out[1:(fin_ind - 1),]


#NA = station operating in that year/cell, but bird not observed
#0 = station operating in that year/cell, bird observed, but not observed breeding
#JDAY = station operating in that year/cell, bird observed breeding

#write to rds
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
saveRDS(MAPS_out2, paste0('breeding_MAPS_obs_', Sys.Date(), '.rds'))





# Process Nestwatch data ---------------------------------------------------------------

#ebird specices codes
setwd(paste0(dir, 'Bird_Phenology/Data/'))

ebird_tax <- read.csv('eBird_Taxonomy_v2018_14Aug2018.csv')
ebird_tax2 <- dplyr::select(ebird_tax, SPECIES_CODE, PRIMARY_COM_NAME, SCI_NAME)
#replace space with underscore
ebird_tax2$SCI_NAME <- gsub(" ", "_", ebird_tax2$SCI_NAME)


#nestwatch data
setwd(paste0(dir, 'Bird_Phenology/Data/Nestwatch'))

nw_data <- read.csv('Nestwatch_2018_1026.csv')

#add grid cells
nw_data$CELL <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                         in_lon_deg = nw_data$LONGITUDE, 
                                         in_lat_deg = nw_data$LATITUDE)[[1]]

#add year
nw_data$YEAR <- substr(nw_data$PROJ_PERIOD_ID, start = 4, stop = 7)

#merge with ebird taxonomy
nw_data2 <- dplyr::left_join(nw_data, ebird_tax2, by = 'SPECIES_CODE')

#only keep entries that have dates - convert to julian day
nw_data2$FIRST_LAY_DT <- format(as.Date(substr(nw_data2$FIRST_LAY_DT, 1, 9), format = '%d%B%Y'), '%j')
nw_data2$HATCH_DT <- format(as.Date(substr(nw_data2$HATCH_DT, 1, 9), format = '%d%B%Y'), '%j')
nw_data2$FLEDGE_DT <- format(as.Date(substr(nw_data2$FLEDGE_DT, 1, 9), format = '%d%B%Y'), '%j')
to.rm <- which(is.na(nw_data2$FIRST_LAY_DT))
nw_data3 <- nw_data2[-to.rm, ]

#filter by species used in IAR
nw_data4 <- dplyr::filter(nw_data3, SCI_NAME %in% species_list_i[,1])

#number of Nestwatch observations for each species
#plyr::count(nw_data4, 'SCI_NAME')


species <- sort(unique(nw_data4$SCI_NAME))

na_reps <- rep(NA, NROW(unique(nw_data4[,c('SCI_NAME', 'YEAR', 'CELL')])))

NW_df <- data.frame(SCI_NAME = na_reps, 
                    YEAR = na_reps,
                    CELL = na_reps,
                    NUM_OBS = na_reps,
                    MEAN_FIRST_LAY = na_reps,
                    SD_FIRST_LAY = na_reps)

counter <- 1
#get metrics from nestwatch data
for (i in 1:length(species))
{
  #i <- 1
  
  #filter by species
  sp <- species[i]
  t_nw <- dplyr::filter(nw_data4, SCI_NAME == sp)
  
  t_cell <- unique(t_nw$CELL)
  for (j in 1:length(t_cell))
  {
    #j <- 1
    t_nw2 <- dplyr::filter(t_nw, CELL == t_cell[j])
    
    t_yr <- unique(t_nw2$YEAR)
    
    for (k in 1:length(t_yr))
    {
      #k <- 2
      t_nw3 <- dplyr::filter(t_nw2, YEAR == t_yr[k])
      
      NW_df[counter,] <- sp
      NW_df[counter,'YEAR'] <- t_yr[k]
      NW_df[counter,'CELL'] <- t_cell[j]
      NW_df[counter,'NUM_OBS'] <- NROW(t_nw3)
      NW_df[counter,'MEAN_FIRST_LAY'] <- round(mean(as.numeric(t_nw3$FIRST_LAY_DT)), 3)
      NW_df[counter,'SD_FIRST_LAY'] <- round(sd(as.numeric(t_nw3$FIRST_LAY_DT)), 3)
      counter <- counter + 1
    }
  }
}



#NA for sd = only one obs
#very large sd in observed lay date for a given cell/year

#write to rds
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
saveRDS(NW_df, paste0('breeding_NW_', Sys.Date(), '.rds'))

