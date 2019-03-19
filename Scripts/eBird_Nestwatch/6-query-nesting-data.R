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

SL <- "'Vireo olivaceus'"



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
                                      AND lng BETWEEN -95 AND -50
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
  
  names(sdata)[6] <- "bba_category"
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
    
    names(sdata)[6] <- "bba_category"
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




# MAPS obs - DB -----------------------------------------------------------

pass <- readLines('db_pass.txt')

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = pass, 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")

maps_data <- DBI::dbGetQuery(cxn, paste0("SELECT lng, lat, year, day,
                                      (event_json ->> 'STATION') AS station,
                                      (place_json ->> 'HABITAT') AS habitat,
                                      (place_json ->> 'ELEV')::int AS elev,
                                      (count_json ->> 'C') capture_code,
                                      (count_json ->> 'BAND') AS band_id,
                                      (count_json ->> 'AGE')::int AS age,
                                      (count_json ->> 'SEX') AS sex,
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
                                      AND day < 300
                                      AND lng BETWEEN -95 AND -50
                                      AND lat > 26
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

head(maps_data)



# add hex cells -----------------------------------------------------------

hexgrid6 <- dggridR::dgconstruct(res = 6) 
maps_data$cell <- dggridR::dgGEO_to_SEQNUM(hexgrid6, 
                                          in_lon_deg = maps_data$lng, 
                                          in_lat_deg = maps_data$lat)[[1]]



# fill in true ages ------------------------------------------------------------

#add col for true age (actual age, not 'year' of bird)
maps_data$true_age <- NA

#fill 0 for all birds banded in birth year
maps_data$true_age[which(maps_data$age %in% c(2, 4))] <- 0 

#ids for birds banded in birth year
ids <- unique(as.numeric(maps_data$band_id))


#enter in true age for all birds banded in birth year
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
    if (sum(temp$age == 5) > 0 & length(t_yrs) > 1)
    {
      birth_year <- (t_yrs[1] - 1)
    } else {
      if (sum(temp$age == 7) > 0 & length(t_yrs) > 1)
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
      maps_data[which(maps_data$band_id == ids[i] & maps_data$year == t_yrs[j]), 
                'true_age'] <- (t_yrs[j] - birth_year)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

#proportion of captures that have a True_age
sum(!is.na(maps_data$true_age))/NROW(maps_data)


# setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# saveRDS(MAPS_mrg5, 'MAPS_age_filled.rds')
# MAPS_mrg5 <- readRDS('MAPS_age_filled.rds')






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

