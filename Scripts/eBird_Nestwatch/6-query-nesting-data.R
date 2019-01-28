####################
# 7 - query other nesting phenology data
#
# *eBird breeding codes - approx 7 hours
#   -0 if not observed in that survey at all (independent of breeding code)
#   -NA if observed but no breeding code recorded
#   -letter code if observed and breeding code recorded
# *MAPS observations
# *MAPS data
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

species_list_i <- read.table('IAR_species_list.txt', stringsAsFactors = FALSE)

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
# * day of year (< julian day 200)
# * time started before 16:00
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
                                    JOIN taxons USING (taxon_id)
                                    WHERE dataset_id = 'ebird'
                                    AND year > 2001
                                    AND day < 200
                                    AND lng BETWEEN -100 AND -50
                                    AND lat > 26
                                    AND (sci_name IN (", SL,"))
                                    AND (event_json ->> 'ALL_SPECIES_REPORTED')::int = 1
                                    AND (event_json ->> 'DURATION_MINUTES')::int BETWEEN 6 AND 1440
                                    AND LEFT(started, 2)::int < 16
                                    AND RADIUS < 100000;
                                    "))


#only values with unique group identifiers (or no group identifier)
data2 <- data[!duplicated(data[,'group_identifier'], 
                          incomparables = NA),]

rm(data)


#calculate polynomial then center data
#scaled julian day, scaled julian day^2, and scaled julian day^3
SJDAY  <- as.vector(scale(data2$day, scale = FALSE))
SJDAY2 <- as.vector(scale(data2$day^2, scale = FALSE))
SJDAY3 <- as.vector(scale(data2$day^3, scale = FALSE))

data2$sjday <- SJDAY
data2$sjday2 <- SJDAY2
data2$sjday3 <- SJDAY3

#scaled effort hours
SHR <- as.vector(scale((data2$duration_minutes/60), scale = FALSE))

data2$shr <- SHR




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
                                      JOIN taxons USING (taxon_id)
                                      WHERE dataset_id = 'ebird'
                                      AND year > 2001
                                      AND day < 200
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
                         year, day, cell, sjday, sjday2, 
                         sjday3, shr, species_list_i[i,1])
  
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
                                        JOIN taxons USING (taxon_id)
                                        WHERE dataset_id = 'ebird'
                                        AND year > 2001
                                        AND day < 200
                                        AND lng BETWEEN -100 AND -50
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
                           year, day, cell, sjday, sjday2, 
                           sjday3, shr, m_sp[i])
    
    names(sdata)[8] <- "bba_category"
    sdata['species'] <- m_sp[i]
    
    saveRDS(sdata, file = paste0('ebird_NA_breeding_cat_', m_sp[i], '.rds'))
    DBI::dbDisconnect(cxn)
  }
}


# copy script to query folder for records ---------------------------------

system(paste0('cp ', dir, 'Bird_Phenology/Scripts/ebird_Nestwatch/7-query-nesting-data.R ', 
              dir, 'Bird_Phenology/Data/', query_dir_path, '/7-query-nesting-data-', Sys.Date(), '.R'))


proc.time() - tt





# counts for breeding codes -----------------------------------------------

#IAR input data (to get relevant cells)
DATE_MA <- '2018-11-12'

setwd(paste0(dir, 'Bird_phenology/Data/Processed/IAR_input_', DATE_MA))
df_master <- readRDS(paste0('IAR_input-', DATE_MA, '.rds'))

nsp <- NROW(species_list_i)

output_df <- data.frame()
for (i in 1:nsp)
{
  #i <- 17
  #read in ebird breeding code data
  DATE_BC <- '2019-01-15'
  
  setwd(paste0(dir, 'Bird_phenology/Data/Processed/breeding_cat_query_', DATE_BC))
  temp_bc <- readRDS(paste0('ebird_NA_breeding_cat_', species_list_i[i,1], '.rds'))
  temp_master <- dplyr::filter(df_master, species == species_list_i[i,1])
  
  #only cells that are in IAR input
  kp_cells <- unique(temp_master$cell)
  temp_bc_f <- dplyr::filter(temp_bc, cell %in% kp_cells)
  
  #probable/confirmed
  t_C34 <- dplyr::filter(temp_bc_f,
                         bba_category == 'C3' |
                           bba_category == 'C4')
  
  #how many C3/C4 obs for a given year
  tt <- plyr::count(t_C34, vars = c('year'))
  #how many C3/C4 obs for a given year/cell
  cy <- plyr::count(t_C34, vars = c('cell', 'year'))
  #filter by > 20 obs in cell
  gr20 <- dplyr::filter(cy, freq > 20)
  #how cells in each year have # C3/C4 obs > 20
  nc_gr20 <- plyr::count(gr20[,1:2], vars = c('year'))
  colnames(nc_gr20)[2] <- 'num_usable_cells'
  #merge with number of C3/C4 obs
  mrg <- dplyr::left_join(tt, nc_gr20, by = 'year')
  colnames(mrg)[2] <- 'num_br_obs'
  #merge with years
  nn <- data.frame(year = 2002:2017)
  mrg2 <- dplyr::left_join(nn, mrg, by = 'year')
  #insert zeros for NA vals
  to.z.use <- which(is.na(mrg2[,'num_usable_cells']))
  mrg2[to.z.use, 'num_usable_cells'] <- 0
  to.z.obs <- which(is.na(mrg2[,'num_br_obs']))
  mrg2[to.z.obs, 'num_br_obs'] <- 0
  
  t_out <- data.frame(species = species_list_i[i,1], mrg2)
  output_df <- rbind(output_df, t_out)
}


num_years <- aggregate(num_usable_cells ~ species, data = output_df, FUN = function(x) sum(x > 0))
colnames(num_years)[2] <- 'num_years'

summary_output_df <- aggregate(num_usable_cells ~ species, data = output_df, FUN = sum)
colnames(summary_output_df)[2] <- 'num_cell_years'
summary_output_df$num_years <- num_years[,'num_years']


setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
write.csv(output_df, 'br_code_data_avail.csv', row.names = FALSE)
write.csv(summary_output_df, 'summary_br_code_data_avail.csv', row.names = FALSE)

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


setwd(paste0(dir, 'Bird_phenology/Data/MAPS_Obs'))

MAPS_obs <- read.csv('1117M.csv', stringsAsFactors = FALSE)

#0 - observed but not breeding
#P - probably breeding
#C - confirmed breeding

#PS1 - May 1-May 10
#PS2 - May 11-May 20
#PS3 - May 21-May 30
#PS4 - May 31-June 9
#PS5 - June 10-June 19
#PS6 - June 20-June 29
#PS7 - June 30-July 9
#PS8 - July 10-July 19
#PS9 - July 20-July 29
#PS10 - July 30-August 8
#PS11 - August 9-August 18


#get grid cell of each station by lat/lon
setwd(paste0(dir, 'Bird_phenology/Data/MAPS_Obs/CntrlStations'))

MAPS_stations <- read.csv('STATIONS.csv', skipNul = TRUE, stringsAsFactors = FALSE)

#STA2 appears to be unqiue station code
# head(MAPS_obs)
# head(MAPS_stations)
# unique(MAPS_obs$STA)
# unique(MAPS_stations$STA2)

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

#insert lat/lons for these stations
st_rep <- data.frame(STATION = c('EASA', 'ANWR', 'REGR', 'CEFB'), 
                     LAT = c(30.54, 30.54, 33.78, 34.68),
                     LNG = c(-95.48, -91.75, -93.82, -82.81))

for (i in 1:NROW(st_rep))
{
  #i <- 2
  tt <- which(MAPS_mrg$STATION.x == st_rep$STATION[i])
  MAPS_mrg[tt,'DECLAT'] <- st_rep$LAT[i]
  MAPS_mrg[tt,'DECLNG'] <- st_rep$LNG[i]
}

MAPS_mrg$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                          in_lon_deg = MAPS_mrg$DECLNG, 
                                          in_lat_deg = MAPS_mrg$DECLAT)[[1]]

#filter by species, get breeding date (breeding period in this case) for each year/cell

setwd(paste0(dir, 'Bird_Phenology/Data/MAPS_Obs'))
sp_codes <- read.csv('sp_4_letter_codes.csv')





