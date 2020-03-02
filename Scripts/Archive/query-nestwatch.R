# Query nestwatch 


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
