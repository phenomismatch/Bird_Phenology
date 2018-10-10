####################
#Sample DB access code from Rafe
####################

#https://db.rstudio.com/dplyr/



cy_dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(dbplyr)
library(RPostgreSQL)


# set wd ------------------------------------------------------------------

setwd(paste0(cy_dir, 'Bird_Phenology/Data/'))



# import eBird species list -----------------------------------------------------

species_list <- read.table('eBird_species_list.txt', stringsAsFactors = FALSE)

#remove underscore and coerce to vector
SL <- as.vector(apply(species_list, 2, function(x) gsub("_", " ", x)))


# access DB ---------------------------------------------------------------

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = "change_me", 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")

#IGNORE WARNINGS (Rob said no big deal abotu JSON warnings)
options(warn=-1)
places <- dplyr::tbl(cxn, "places")
events <- dplyr::tbl(cxn, "events")
counts <- dplyr::tbl(cxn, "counts")
taxons <- dplyr::tbl(cxn, "taxons")

codes <- dplyr::tbl(cxn, "codes")
countries <- dplyr::tbl(cxn, "countries")
datasets <- dplyr::tbl(cxn, "datasets")
version <- dplyr::tbl(cxn, "version")
options(warn=0)


# filter dataset ----------------------------------------------------------

test <- counts %>%
  inner_join(events, by = 'event_id') %>%
  inner_join(places, by = 'place_id') %>%
  inner_join(taxons, by = 'taxon_id') %>%
  filter(dataset_id == 'ebird') %>%
  filter(year > 2001) %>%
  filter(day < 200) %>%
  filter(lat > 26) %>%
  filter(lng > -100) %>%
  filter(lng < -50) %>%
  #filter(radius > 100000) %>%
  filter(sci_name %in% SL) %>%
  head(1000) %>%
  select(c(event_id, year, day, 
           place_id, lat, lng, 
           started, ended, count, radius, 
           sci_name, common_name)) 

t2 <- as.data.frame(test)

