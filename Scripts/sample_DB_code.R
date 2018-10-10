####################
#Sample DB access code from Rafe
####################

#https://db.rstudio.com/dplyr/


library(dplyr)
library(dbplyr)
library(RPostgreSQL)

pg <- DBI::dbDriver("PostgreSQL")

cxn <- DBI::dbConnect(pg, 
                      user = "cyoungflesh", 
                      password = "change_me", 
                      host = "35.221.16.125", 
                      port = 5432, 
                      dbname = "sightings")

#IGNORE WARNINGS

places <- dplyr::tbl(cxn, "places")
events <- dplyr::tbl(cxn, "events")
counts <- dplyr::tbl(cxn, "counts")
taxons <- dplyr::tbl(cxn, "taxons")


codes <- dplyr::tbl(cxn, "codes")
countries <- dplyr::tbl(cxn, "countries")
datasets <- dplyr::tbl(cxn, "datasets")
version <- dplyr::tbl(cxn, "version")


test <- counts %>%
  inner_join(events, by = 'event_id') %>%
  inner_join(places, by = 'place_id') %>%
  inner_join(taxons, by = 'taxon_id') %>%
  filter(dataset_id == 'ebird') %>%  #only ebird data
  filter(year > 2001) %>% #filter based on year
  filter(day < 200) %>% #filter based on julian day
  filter(common_name == 'Yellow Warbler')
#filter based on effort hours
#filter based on effort distance
#filter based on whether or not it was a group checklist (need only one of each of these)


#NEED:
#GROUP IDENTIFIER in eBird basic dataset
#EFFORT DISTANCE KM




t2 <- as.data.frame(head(test))


t2

# "dplyr" is lazy and only evaluates as little as possible to create the desired output.

# This helps with efficiency. So, you may want to avoid the as.data.frame() function,
# it is there only to show how to convert the queries into a dataframe.


all_tables <- places %>% 
  filter(between(lng, -73, -72)) %>%
  filter(between(lat, 40, 44)) %>%
  inner_join(events, by = "place_id") %>%
  filter(between(year, 2010, 2014)) %>%
  filter(between(day, 100, 200)) %>%
  head(100) %>%
  inner_join(counts, by = "event_id") %>%
  inner_join(taxons, by = "taxon_id") %>%
  head(100) %>%
  as.data.frame()

