####################
# Investigate changes in wing_chord over age
#
# Database query:
# ---------------
# Day between 100 and 300
# Long between -95 and -50
# Lat > 26
# Only capture code 'N', 'R', 'U' (new, recaptures, unbanded birds)
# 
# 
# Data processing:
# ----------------
# Fill true age for each bird for each year. Done using birth year derived from birds either caught during the hatch year or birds listed as Second Year Bird or Third Year. Remove all 0 values for weight and wing_chord. Filter to only known age adults (greater than 1 year old AKA Second Year or later).
# 
# For each individual band_id, take mean wing_chord for Age 1 (Second Year) and mean wing chord for Age 2+ (Third Year +); also calculate wing_chord coefficient of variation (sd/mean) Age 1 and Age 2+. Remove any individuals that have Age 1 CV of greater than 0.03 and Age 2+ CV of greater than 0.05 (sds that are greater than 3% and 5% of mean, respectively). 
# 
# 
# Analysis:
# ---------
# Run paired two-tailed t-test for each species that had measurements for greater than 10 individuals (to test difference in wing_chord between Age 1 and Age 2+ age classes). Estimated difference and p-values returned for each species.
####################


# top-level dir -----------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)



# read in data ------------------------------------------------------------

#if after query
# setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# maps_data <- readRDS('MAPS_wing_chord_age_filled.rds')


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

#Entire North America

maps_data <- DBI::dbGetQuery(cxn, paste0("SELECT lng, lat, year, day, common_name, sci_name, event_id, started, ended,
                                         (count_json ->> 'ANET') AS anet,
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

# setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# saveRDS(maps_data, 'MAPS_wing_chord_age_filled.rds')







# Process wing chord data --------------------------------------------


#remove records that have weight 0 and wing chord 0
to.rm <- which(maps_data$weight == 0 | maps_data$wing_chord == 0 | 
                 is.na(maps_data$weight) | is.na(maps_data$fat_content) | is.na(maps_data$wing_chord))
maps_data_qc <- maps_data[-to.rm, ]

#only known-age adults
maps_adults_qc <- dplyr::filter(maps_data_qc, true_age > 0)

#find unique band ids
bid_cnt <- plyr::count(maps_adults_qc$band_id)

#individuals that have been captured more than 2 times
bid_c <- dplyr::filter(bid_cnt, freq > 2)
c_birds <- bid_c[,1]

#only band_ids of interest
maps_c <- dplyr::filter(maps_data_qc, band_id %in% c_birds)

#band ids
nid <- unique(maps_c$band_id)



# plots metrics over age ---------------------------------------------------

#weight
ggplot(maps_c, aes(true_age, weight, col = band_id)) +
  geom_line() +
  theme(legend.position="none")

#sweight
ggplot(maps_c, aes(true_age, (weight/wing_chord), col = band_id)) +
  geom_line() +
  theme(legend.position="none")

#fat content
ggplot(maps_c, aes(true_age, fat_content, col = band_id)) +
  geom_line() +
  theme(legend.position="none")

#wing chord
ggplot(maps_c, aes(true_age, wing_chord, col = band_id)) +
  geom_line() +
  theme(legend.position="none")



# calc mean and cv wing chord for each bird -----------------------------------------------------------

#For each individual, calculate mean and CV (sd/mean) for Age 1 wing chord and Age 2+ wing chord
ttd <- data.frame()
for (i in 1:length(nid))
{
  #i <- 6
  temp <- dplyr::filter(maps_c, band_id == nid[i])
  
  #only if have at least one year age > 1 and one year age = 1
  if (sum(temp$true_age > 1) > 0 & sum(temp$true_age == 1) > 0)
  {
    #age 1
    a1 <- mean(dplyr::filter(temp, true_age == 1)$wing_chord)
    a1_sd <- sd(dplyr::filter(temp, true_age == 1)$wing_chord)
    
    #age 2+
    a2p <- mean(dplyr::filter(temp, true_age > 1)$wing_chord)
    a2p_sd <- sd(dplyr::filter(temp, true_age > 1)$wing_chord)
    
    temp <- data.frame(sp = temp$sci_name[1], band_id = nid[i], 
                       a1_wc = a1, a1_cv = a1_sd/a1, 
                       a2p_wc = a2p, a2p_cv = a2p_sd/a2p)
    ttd <- rbind(ttd, temp)
  }
}


#if sd within year 1 is greater than 3% of mean - exclude
#if sd year 2+ is greater than 5% of mean - exclude
to.rm.a1 <- which(ttd$a1_cv > 0.03)
to.rm.a2p <- which(ttd$a2p_cv > 0.05)
ctm <- c(to.rm.a1, to.rm.a2p)
ctm2 <- ctm[!(duplicated(ctm) | duplicated(ctm, fromLast = TRUE))]
#remove values
ttd2 <- ttd[-ctm2,]

#only species that have > 10 individuals
csp <- plyr::count(ttd2, 'sp')
nsp <- csp[which(csp$freq > 10), 1]

#filter relevant species
ttd3 <- dplyr::filter(ttd2, sp %in% nsp) 


# paired t-test for each species ------------------------------------------

#run paired t-test and calc mean wing chord and sd wing chord for each species
out <- data.frame()
for (i in 1:length(nsp))
{
  #i <- 1
  temp <- dplyr::filter(ttd3, sp == nsp[i])
  tfit <- t.test(temp$a1_wc, temp$a2p_wc, paired = TRUE)
  tt <- data.frame(sp = nsp[i],
                   a1_mn = mean(temp$a1_wc),
                   a1_sd = sd(temp$a1_wc, na.rm = TRUE),
                   a2p_mn = mean(temp$a2p_wc),
                   a2p_sd = sd(temp$a2p_wc, na.rm = TRUE),
                   est = round(tfit$estimate, 2), 
                   p = round(tfit$p.value, 2))
  out <- rbind(out, tt)
}
row.names(out) <- NULL



# plot results ------------------------------------------------------------


#histograms of estimated differences and p-values
hist(out$est)
hist(out$p)




#boxplot
tplt <- reshape2::melt(ttd3[,c('sp', 'a1_wc', 'a2p_wc')], id = 'sp')
tplt$spvar <- interaction(tplt$sp, tplt$var)

ggplot(aes(y = value, x = sp, fill = variable), data = tplt) + 
  geom_boxplot() +
  theme_bw() +
  theme(legend.position="none")



#plots of all individuals used in final analysis
nmc <- dplyr::filter(maps_c, band_id %in% ttd3$nid)
#wing chord
ggplot(nmc, aes(true_age, wing_chord, col = band_id)) +
  geom_line(alpha = 0.5) +
  theme(legend.position="none") +
  theme_bw()


