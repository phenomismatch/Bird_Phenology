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
# Fill true age for each bird for each year. Done using birth year derived from birds either caught during the hatch year or birds listed as Second Year Bird. Remove all 0 values for weight and wing_chord. Filter to only known age adults (greater than 1 year old AKA Second Year or later).
# 
# For each individual band_id, take mean wing_chord for Age 1 (Second Year) and mean wing chord for Age 2+ (Third Year +); also calculate wing_chord coefficient of variation (sd/mean) Age 1 and Age 2+. Remove any individuals that have Age 1 CV of greater than 0.03 and Age 2+ CV of greater than 0.05 (sds that are greater than 3% and 5% of mean, respectively).
# 
# 
# Analysis:
# ---------
# Run paired two-tailed t-test for each species that had measurements for greater than 20 individuals (to test difference in wing_chord between Age 1 and Age 2+ age classes). Estimated difference and p-values returned for each species.
####################


# top-level dir -----------------------------------------------------------

dir <- '~/Google_Drive/R/'


# Load packages -----------------------------------------------------------

library(dplyr)
library(RPostgreSQL)
library(DBI)
library(ggplot2)
library(rstan)
library(MCMCvis)


# # query DB ----------------------------------------------------------------
# 
# 
# setwd(paste0(dir, 'Bird_Phenology/Data/'))
# 
# pass <- readLines('db_pass.txt')
# 
# pg <- DBI::dbDriver("PostgreSQL")
# 
# cxn <- DBI::dbConnect(pg, 
#                       user = "cyoungflesh", 
#                       password = pass, 
#                       host = "35.221.16.125", 
#                       port = 5432, 
#                       dbname = "sightings")
# 
# #Entire North America - only days 100-300
# 
# maps_data <- DBI::dbGetQuery(cxn, paste0("SELECT lng, lat, year, day, common_name, sci_name, event_id, started, ended,
#                                          (count_json ->> 'ANET') AS anet,
#                                          (event_json ->> 'STATION') AS station,
#                                          (place_json ->> 'HABITAT') AS habitat,
#                                          (place_json ->> 'ELEV')::int AS elev,
#                                          (count_json ->> 'C') capture_code,
#                                          (count_json ->> 'BAND') AS band_id,
#                                          (count_json ->> 'AGE')::int AS age,
#                                          (count_json ->> 'SEX') AS sex,
#                                          (count_json ->> 'FW') AS feather_wear,
#                                          (count_json ->> 'CP')::int AS cloacal_pro,
#                                          (count_json ->> 'BP')::int AS brood_patch,
#                                          (count_json ->> 'F')::int AS fat_content,
#                                          (count_json ->> 'WNG')::float AS wing_chord,
#                                          (count_json ->> 'WEIGHT')::float AS weight,
#                                          (count_json ->> 'N') AS standard_effort
#                                          FROM places
#                                          JOIN events USING (place_id)
#                                          JOIN counts USING (event_id)
#                                          JOIN taxa USING (taxon_id)
#                                          WHERE events.dataset_id = 'maps'
#                                          AND day BETWEEN 100 AND 300
#                                          AND (count_json ->> 'C') in ('N', 'R', 'U')
#                                          AND (count_json ->> 'N') in ('-', 'D', 'S', 'T', 'O', '+');
#                                          "))
# 
# #C (Capture code): Only codes N,R,U are used for analyses according to MAPS docs
# #N = newly banded bird
# #R = recaptured bird
# #U = unbanded bird
# 
# #N (Indicator to include in productivity and survivorship analyses): Okay to use these for analyses that do not require to control for effort
# #0 = not caught at MAPS station
# #T = outside normaps MAPS operation for that station
# #S = caught within MAPS station boundary but not in a MAPS net
# #D = date outside of MAPS periods
# #- = record examined with current MAPS analytical procedure (taken to be standard capture methods)
# #+ = record examined with preliminary MAPS analytical procedure
# 
# #Age: 
# #0 - unknown
# #4 - local (young bird incapable of flight)
# #2 - hatching-year bird
# #1 - after hatching-year bird
# #5 - second-year bird
# #6 - after second-year bird
# #7 - third-year bird
# #8 - after third-year bird
# 
# #Sex: 
# #M - male
# #F - female
# #U - unknown
# #X - unattempted
# 
# #Cloacal protuberance:
# #0 - none
# #1 - small
# #2 - medium
# #3 - large
# 
# #Brood patch:
# #0 - none
# #1 - smooth
# #2 - vascularized
# #3 - heavy
# #4 - wrinkled
# #5 - molting
# 
# #Fat content:
# #0 - none
# #1 - trace
# #2 - light
# #3 - half
# #4 - full
# #5 - bulging
# #6 - greatly bulging
# #7 - very excessive
# 
# 
# # fill in true ages ------------------------------------------------------------
# 
# #true_age 1 means 1 year old (born in previous season)
# 
# #add col for true age (actual age, not 'year' of bird)
# maps_data$true_age <- NA
# 
# #fill 0 for all birds banded in birth year
# #2 = hatching year bird
# #4 = local (young bird incapable of flight))
# #5 = second year bird
# #7 = third year bird - unreliable due to molting in second year (id as third year bird late second year)
# maps_data$true_age[which(maps_data$age %in% c(2, 4))] <- 0 
# 
# #ids for birds banded in birth year
# ids <- unique(as.numeric(maps_data$band_id))
# 
# # enter in true age for all birds banded in birth year (or entered as second year or third year)
# # create progress bar
# pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
# for (i in 1:length(ids))
# {
#   #i <- 61
#   #filter for band_id (all years)
#   temp <- dplyr::filter(maps_data, band_id == ids[i])
#   t_yrs <- sort(as.numeric(unique(temp$year)))
#   
#   birth_year <- NA
#   
#   #define birth year based on assigned age class (band in birth year trumps other ages)
#   if (sum(!is.na(temp$true_age)) > 0 & length(t_yrs) > 1)
#   {
#     birth_year <- t_yrs[1]
#   } else {
#     if (sum(temp$age == 5) > 0)
#     {
#       birth_year <- (t_yrs[1] - 1)
#     }
#   }
#   
#   #fill in age values for subsequent year if birth year was defined
#   if (!is.na(birth_year))
#   {
#     for (j in 1:length(t_yrs))
#     {
#       #j <- 1
#       #insert age
#       maps_data[which(maps_data$band_id == ids[i] & maps_data$year == t_yrs[j]), 
#                 'true_age'] <- (t_yrs[j] - birth_year)
#     }
#   }
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# 
# #one time fill with feather wear to avoid refilling age data
# # qq <- readRDS('MAPS-age-filled.rds')
# # qq$feather_wear <- maps_data$feather_wear
# # setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# # saveRDS(qq, 'MAPS-age-filled.rds')
# 
# 
# #save RDS file
# # setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# # saveRDS(maps_data, 'MAPS-age-filled.rds')



# Process wing chord data --------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
#setwd(paste0(dir, '../'))
maps_data <- readRDS('MAPS-age-filled.rds')


#remove records that have wing chord 0/NA nad feather_wear NA
to.rm <- which(maps_data$wing_chord == 0 | is.na(maps_data$wing_chord) | is.na(maps_data$feather_wear))
maps_data_qc <- maps_data[-to.rm, ]

#only known-age adults
maps_adults_qc <- dplyr::filter(maps_data_qc, true_age > 0)


#remove WC values outside 3 sd of mean
usp_maps_ad <- unique(maps_adults_qc$sci_name)
maps_adults_qc2 <- data.frame()
trm <- 0
for (i in 1:length(usp_maps_ad))
{
  #i <- 7
  temp <- dplyr::filter(maps_adults_qc, sci_name == usp_maps_ad[i])
  sd_wc <- sd(temp$wing_chord)
  mn_wc <- mean(temp$wing_chord)
  #outside 3 sds
  low <- mn_wc - 3*sd_wc
  high <-  mn_wc + 3*sd_wc
  to.rm <- which(temp$wing_chord < low | temp$wing_chord > high)
  if (length(to.rm) > 0)
  {
    tt <- temp[-to.rm,]
  } else {
    tt <- temp
  }
  
  trm <- trm + length(to.rm)
  maps_adults_qc2 <- rbind(maps_adults_qc2, tt)
}


#find unique band ids
bid_cnt <- plyr::count(maps_adults_qc2$band_id)

#individuals that have been captured at least 2 times
bid_c <- dplyr::filter(bid_cnt, freq >= 2)
c_birds <- bid_c[,1]

#only band_ids of interest
maps_c <- dplyr::filter(maps_adults_qc2, band_id %in% c_birds)

#band ids
nid <- unique(maps_c$band_id)



# plots metrics over age ---------------------------------------------------

# #weight
# ggplot(maps_c, aes(true_age, weight, col = band_id)) +
#   geom_line() +
#   theme(legend.position="none")
# 
# #sweight
# ggplot(maps_c, aes(true_age, (weight/wing_chord), col = band_id)) +
#   geom_line() +
#   theme(legend.position="none")
# 
# #fat content
# ggplot(maps_c, aes(true_age, fat_content, col = band_id)) +
#   geom_line() +
#   theme(legend.position="none")
# 
# #wing chord
# ggplot(maps_c, aes(true_age, wing_chord, col = band_id)) +
#   geom_line() +
#   theme(legend.position="none")
# 
# plot(maps_c$true_age, maps_c$feather_wear)


# calc mean and cv wing chord for each bird -----------------------------------------------------------

#For each individual, calculate mean and CV (sd/mean) for Age 1 wing chord and Age 2+ wing chord
# ttd <- data.frame()
# for (i in 1:length(nid))
# {
#   #i <- 2
#   temp <- dplyr::filter(maps_c, band_id == nid[i])
#   
#   #only if have at least one year age > 1 and one year age = 1
#   if (sum(temp$true_age > 1) > 0 & sum(temp$true_age == 1) > 0)
#   {
#     #age 1
#     a1 <- mean(dplyr::filter(temp, true_age == 1)$wing_chord)
#     a1_sd <- sd(dplyr::filter(temp, true_age == 1)$wing_chord)
#     
#     #age 2+
#     a2p <- mean(dplyr::filter(temp, true_age > 1)$wing_chord)
#     a2p_sd <- sd(dplyr::filter(temp, true_age > 1)$wing_chord)
#     
#     temp <- data.frame(sp = temp$sci_name[1], cn = temp$common_name[1],
#                        band_id = nid[i], 
#                        a1_wc = a1, a1_cv = a1_sd/a1, 
#                        a2p_wc = a2p, a2p_cv = a2p_sd/a2p)
#     ttd <- rbind(ttd, temp)
#   }
# }
# 
# 
# #if sd within year 1 is greater than 3% of mean - exclude
# #if sd year 2+ is greater than 5% of mean - exclude
# to.rm.a1 <- which(ttd$a1_cv > 0.03)
# to.rm.a2p <- which(ttd$a2p_cv > 0.05)
# ctm <- c(to.rm.a1, to.rm.a2p)
# ctm2 <- ctm[!(duplicated(ctm) | duplicated(ctm, fromLast = TRUE))]
# #remove values
# ttd2 <- ttd[-ctm2,]
# 
# #only species that have > 20 individuals
# csp <- plyr::count(ttd2, 'sp')
# nsp <- csp[which(csp$freq > 20), 1]
# 
# #filter relevant species
# ttd3 <- dplyr::filter(ttd2, sp %in% nsp) 


# paired t-test for each species ------------------------------------------

#run paired t-test and calc mean wing chord and sd wing chord for each species
# out <- data.frame()
# for (i in 1:length(nsp))
# {
#   #i <- 1
#   temp <- dplyr::filter(ttd3, sp == nsp[i])
#   tfit <- t.test(temp$a1_wc, temp$a2p_wc, paired = TRUE)
#   tt <- data.frame(sp = nsp[i],
#                    cn = temp$cn[1],
#                    a1_mn = round(mean(temp$a1_wc), 2),
#                    a1_sd = round(sd(temp$a1_wc, na.rm = TRUE), 2),
#                    a2p_mn = round(mean(temp$a2p_wc), 2),
#                    a2p_sd = round(sd(temp$a2p_wc, na.rm = TRUE), 2),
#                    est = round(tfit$estimate, 3), 
#                    p = round(tfit$p.value, 5))
#   out <- rbind(out, tt)
# }
# row.names(out) <- NULL
# 
# #Benjamini & Hochberg p-value correction
# out$p_adj <-  round(p.adjust(out$p, method="BH"), 3)



# plot results ------------------------------------------------------------


# #histograms of estimated differences and p-values
# hist(out$est)
# hist(out$p)
# 
# 
# 
# 
# #boxplot
# tplt <- reshape2::melt(ttd3[,c('sp', 'a1_wc', 'a2p_wc')], id = 'sp')
# tplt$spvar <- interaction(tplt$sp, tplt$var)
# 
# ggplot(aes(y = value, x = sp, fill = variable), data = tplt) + 
#   geom_boxplot() +
#   theme_bw() +
#   theme(legend.position="none")



#plots of all individuals used in final analysis
# nmc <- dplyr::filter(maps_c, band_id %in% ttd3$band_id)
# #wing chord
# ggplot(nmc, aes(true_age, wing_chord, col = band_id)) +
#   geom_line() +
#   theme_bw() +
#   theme(legend.position="none")
  





#plot of before and after lines for each species



# setwd(paste0(dir, '/Bird_Phenology/Data/Traits'))
# traits <- read.csv('Trait_database-2019-03-28.csv')
# 
# head(out)
# head(traits)
# 
# j_out <- dplyr::left_join(out, traits, by = c('sp' = 'SCI_NAME'))
# head(j_out)
# 
# 
# plot(j_out$a1_mn, j_out$est, xlab = 'SY Wing Chord', 
#      ylab = 'Diff')
# fit <- lm(est ~ a2p_mn, data = j_out)
# abline(fit, col = 'red')
# 
# res <- residuals(fit)
# plot(res, j_out$est/j_out$a2p_mn)
# plot(res, j_out$MIGRATION_DISTANCE_LASORTE, 
#      xlab = 'Wing Chord Residuals', ylab = 'Migration Distance')
# plot(res, j_out$SPRING_MIGRATION_SPEED_LASORTE,
#      xlab = 'Wing Chord Residuals', ylab = 'Migration Speed')
# 
# 
# plot((j_out$est/j_out$a2p_mn), j_out$MIGRATION_DISTANCE_LASORTE)
# plot((j_out$est/j_out$a2p_mn), j_out$SPRING_MIGRATION_SPEED_LASORTE)
# plot((j_out$est), j_out$MIGRATION_DISTANCE_LASORTE)
# plot((j_out$est), j_out$SPRING_MIGRATION_SPEED_LASORTE)
# 
# 
# plot(j_out$BODY_MASS_ELTON, j_out$est)
# plot(j_out$MIGRATION_DISTANCE_LASORTE, j_out$est)
# plot(j_out$SPRING_MIGRATION_SPEED_LASORTE, j_out$est)
# plot(j_out$FALL_MIGRATION_SPEED_LASORTE, j_out$est)
# 
# plot(j_out$BODY_MASS_ELTON, j_out$MIGRATION_DISTANCE_LASORTE)
# plot(j_out$BODY_MASS_ELTON, j_out$SPRING_MIGRATION_SPEED_LASORTE)
# plot(j_out$MIGRATION_DISTANCE_LASORTE, j_out$SPRING_MIGRATION_SPEED_LASORTE)
# plot(j_out$MIGRATION_DISTANCE_LASORTE, j_out$FALL_MIGRATION_SPEED_LASORTE)
# plot(j_out$SPRING_MIGRATION_SPEED_LASORTE, j_out$FALL_MIGRATION_SPEED_LASORTE)
# plot(j_out$EGG_MASS_BONA, j_out$BODY_MASS_ELTON)




# stan model --------------------------------------------------------------

#// nu prior
#// https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/

#modified BEST
#https://vuorre.netlify.com/post/2017/01/02/how-to-compare-two-groups-with-robust-bayesian-estimation-using-r-stan-and-brms/#bayesian-estimation-of-the-t-test
#https://stats.stackexchange.com/questions/130389/bayesian-equivalent-of-two-sample-t-test

m_band_id <- data.frame()
for (i in 1:length(nid))
{
  #i <- 21003
  temp <- dplyr::filter(maps_c, band_id == nid[i])

  if (length(unique(temp$sci_name)) == 1)
  {
    #only if have at least one year age > 1 and one year age = 1
    if (sum(temp$true_age > 1) > 0 & sum(temp$true_age == 1) > 0)
    {
      #age 1
      temp <- data.frame(sp = temp$sci_name[1], 
                         cn = temp$common_name[1],
                         band_id = temp$band_id[1])
      m_band_id <- rbind(m_band_id, temp)
    }
  }
}


#only species that have > 20 individuals
csp <- plyr::count(m_band_id, 'sp')
nsp <- as.character(csp[which(csp$freq > 20), 1])

#filter relevant species
m_band_id2 <- dplyr::filter(m_band_id, sp %in% nsp)

u_bid <- unique(m_band_id2$band_id)

#filter based on new band_ids
maps_c2 <- dplyr::filter(maps_c, band_id %in% u_bid)


#NROW(maps_data) #1644585 - ALL
#NROW(maps_data_qc) #1472931 - REMOVE WING CHORD NA/0
#NROW(maps_adults_qc) #199256 - ONLY KNOWN AGE ADULTS WITH FEATHER WEAR RECORDS
#NROW(maps_adults_qc2) #198252 - REMOVE RECORDING ERRORS
#NROW(bid_c) #30649 (number of ind) - ONLY BIRDS CAPTURES > 2 times
#NROW(maps_c) #88027 - ONLY BAND_IDS OF INTEREST
#NROW(m_band_id) #13175 (numer of ind) - ONLY BIRDS CAPTURED AS SY AND ASY
#NROW(m_band_id2) #12416 (number of ind) - ONLY RELEVANT SPECIES
#NROW(maps_c2) #43733 - data for relevant individuals


#age ids
#when true_age == 1, x = 1
#when true_age > 1, x = 2
maps_c2$x <- NA
maps_c2$x[which(maps_c2$true_age == 1)] <- 1
maps_c2$x[which(maps_c2$true_age > 1)] <- 2

#sci names as numbers
maps_c2$sci_name_f <- as.numeric(factor(maps_c2$sci_name))

#individual ids as numbers
maps_c2$band_id_f <- as.numeric(factor(maps_c2$band_id))

usp <- unique(maps_c2$sci_name)
ucn <- unique(maps_c2$common_name)
usp_f <- unique(maps_c2$sci_name_f)
u_bid2 <- unique(maps_c2$band_id)
u_bid2_f <- unique(maps_c2$band_id_f)


# #species used in analysis
# sp_n <- cbind(sci_name = usp, common_name = ucn)
# sp_n2 <- sp_n[order(sp_n[,1]),]
# #setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
# setwd(paste0(dir, '../../Desktop'))
# write.csv(sp_n2, 'wing_cord_age_species.csv', row.names = FALSE)


#species id for each individual
sp_j <- rep(NA, length(u_bid2_f))
for (i in 1:length(u_bid2_f))
{
  #i <- 1
  temp <- dplyr::filter(maps_c2, band_id_f == i)
  sp_j[i] <- temp$sci_name_f[1]
}

#mean and sd for each species
mean_sp <- c()
sd_sp <- c()
for (i in 1:length(usp_f))
{
  #i <- 25
  temp <- dplyr::filter(maps_c2, sci_name_f == i)
  temp_mn <- mean(temp$wing_chord)
  temp_sd <- sd(temp$wing_chord)
  
  mean_sp <- c(mean_sp, temp_mn)
  sd_sp <- c(sd_sp, temp_sd)
}

#mean each individual (i) for SY/ASY (k)
#y[i] = wing chord obs
#ind[i] = id # j at obs i
#sp[i] = species # s for measure i
#sp_j[j] = species # s for ind j
#x[i] = stage k
#FW[i] = feather wear for that obs
#alpha[j,k] = intercept for ind j at stage k
#beta[s] = effect of feather wear on wing chord for species s
#sigma[s] = unaccounted for variance in mean as a function of feather wear
#mu_sp[s,k] = mean for species s, stage k
#sigma_sp[s,k] = sd for species s, stage k
#nu = degrees of freedom


#data for stan model
DATA <- list(N = NROW(maps_c2),
             NS = length(usp),
             NJ = length(u_bid2),
             y = maps_c2$wing_chord,
             sp = maps_c2$sci_name_f, #species
             ind = maps_c2$band_id_f, #individuals
             x = maps_c2$x,
             FW = as.numeric(maps_c2$feather_wear) + 1,
             sp_j = sp_j,
             mu_mu_sp = mean_sp,
             mu_sigma_sp = sd_sp)

stanmodel <- "
data {
int<lower=0> N;                     // number of obs
int<lower=0> NS;                    // number of species
int<lower=0> NJ;                    // number of individuals
real<lower=0> y[N];                 // response
int<lower=0> sp[N];                 // species id
int<lower=0> ind[N];                // individual id
int<lower=1> x[N];                  // stage id
int<lower=1> FW[N];                 // feather wear index
int<lower=1> sp_j[NJ];              // species id for individual
real<lower = 0> mu_mu_sp[NS];       // mean species wing chord (both stages)
real<lower = 0> mu_sigma_sp[NS];    // sd species wing chord (both stages)
}

parameters {
real<lower = 0> sigma[NS];
matrix[NJ, 2] alpha;                // individual mean wing chord
real beta[NS];                      // effect of feather wear on wing chord
matrix[NS, 2] mu_sp;                // mean species wing chord (SY/ASY)
matrix<lower = 0>[NS, 2] sigma_sp;  // sd species wing chord (SY/ASY)
real<lower = 0> nu[NS];             // degrees of freedom for t-dist
real mu_beta_raw;
real<lower = 0> sigma_beta_raw;
}

transformed parameters {
real mu_y[N];
real mu_beta;
real<lower = 0> sigma_beta;

mu_beta = mu_beta_raw * 3;
sigma_beta = sigma_beta_raw * 2;

for (i in 1:N)
{
  mu_y[i] = alpha[ind[i], x[i]] + beta[sp[i]] * FW[i];
}
}

model {

for (i in 1:N)
{
  y[i] ~ normal(mu_y[i], sigma[sp[i]]);
}

for (j in 1:NJ)
{
  for (k in 1:2)
  {
    alpha[j,k] ~ student_t(nu[sp_j[j]], mu_sp[sp_j[j], k], sigma_sp[sp_j[j], k]);
  }
}

for (s in 1:NS)
{
  beta[s] ~ normal(mu_beta, sigma_beta);
  nu[s] ~ gamma(2, 0.1);
  
  for (k in 1:2)
  {
    mu_sp[s,k] ~ normal(mu_mu_sp[s], 100);
    sigma_sp[s,k] ~ normal(mu_sigma_sp[s], 100);
  }
}

mu_beta_raw ~ normal(0, 1);
sigma_beta_raw ~ normal(0, 1);
}

generated quantities {
real mu_sp_d[NS];

for (s in 1:NS)
{
  mu_sp_d[s] = mu_sp[s,2] - mu_sp[s,1];
}
}
"

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

DELTA <- 0.95
TREE_DEPTH <- 16
STEP_SIZE <- 0.05
CHAINS <- 4
ITER <- 3000

tt <- proc.time()
fit <- rstan::stan(model_code = stanmodel,
                    data = DATA,
                    chains = CHAINS,
                    iter = ITER,
                    cores = CHAINS,
                    pars = c('sigma',
                             'alpha',
                             'beta',
                             'mu_sp',
                             'sigma_sp',
                             'nu',
                             'mu_beta',
                             'sigma_beta',
                             'mu_sp_d'), 
                    control = list(adapt_delta = DELTA,
                                   max_treedepth = TREE_DEPTH,
                                   stepsize = STEP_SIZE))
run_time <- (proc.time()[3] - tt[3]) / 60

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
saveRDS(fit, file = 'MAPS-wc-age-BEST-stan-output.rds')

# MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2,
#                      excl = c('mu_sp_d','alpha', 'beta', 'mu_sp', 'sigma_sp'))
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'beta')
MCMCvis::MCMCplot(fit, params = 'beta')
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'mu_sp_d')
MCMCvis::MCMCplot(fit, params = 'mu_sp_d')
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'nu')
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'sigma_sp')
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'sigma')
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'mu_beta')
MCMCvis::MCMCsummary(fit, n.eff = TRUE, round = 2, params = 'sigma_beta')

# library(shinystan)
# launch_shinystan(fit)

mu_sp_d <- MCMCvis::MCMCpstr(fit, params = 'mu_sp_d')[[1]]
mu_sp_d_UCI <- MCMCvis::MCMCpstr(fit, params = 'mu_sp_d', 
                                 func = function(x) quantile(x, probs = 0.975))[[1]]
mu_sp_d_LCI <- MCMCvis::MCMCpstr(fit, params = 'mu_sp_d',
                                 func = function(x) quantile(x, probs = 0.025))[[1]]
WC_b <- MCMCvis::MCMCpstr(fit, params = 'mu_sp')[[1]]

mu_sp_d_ch <- MCMCvis::MCMCchains(fit, params = 'mu_sp_d')

sp_df <- data.frame(sp_id = sort(unique(maps_c2$sci_name_f)), 
                    sci_name = sort(unique(maps_c2$sci_name)),
                    wing_chord_ASY = WC_b[,2],
                    mu_sp_d = mu_sp_d,
                    mu_sp_d_LCI = mu_sp_d_LCI,
                    mu_sp_d_UCI = mu_sp_d_UCI)

WC_age <- read.csv('WC-age.csv')

wc_b <- dplyr::full_join(sp_df, WC_age, by = c('sci_name' = 'SCINAME'))

wc_b$sdiff <- wc_b$mu_sp_d / wc_b$wing_chord_ASY
  
#FROM PETER PYLE
# JA - outer primaries juvenile in 'all' SYs 
# JM - outer primaries juvenile in >75% of SYs (eccentric patterns) 
# JS - outer primaries juvenile in <75% of SYs 
# (eccentric patterns except for one species, Northern Cardinal) 
# FE - outer primaries formative in 'all' birds (eccentric patterns) 
# FC - outer primaries formative in all birds (not eccentric patterns) 
# 
# So we will expect to see greater differences 
# between SY and ASY wing chords for JA, a little 
# bit less for JM, about half as much for JS, and 
# less of a difference for the two F categories. 
# Any difference in these will be noteworthy, 
# though, and indicate formative primaries are 
# shorter than basic primaries, which would be a novel thing to report. 



# MCMCvis plots -----------------------------------------------------------

#order posteriors by PSCORE group
JA_wc <- dplyr::arrange(dplyr::filter(wc_b, PSCORE == 'JA'), desc(mu_sp_d))
JM_wc <- dplyr::arrange(dplyr::filter(wc_b, PSCORE == 'JM'), desc(mu_sp_d))
JS_wc <- dplyr::arrange(dplyr::filter(wc_b, PSCORE == 'JS'), desc(mu_sp_d))
FE_wc <- dplyr::arrange(dplyr::filter(wc_b, PSCORE == 'FE'), desc(mu_sp_d))
FC_wc <- dplyr::arrange(dplyr::filter(wc_b, PSCORE == 'FC'), desc(mu_sp_d))

beta_JA <- paste0('beta\\[', JA_wc$sp_id, '\\]')
beta_JM <- paste0('beta\\[', JM_wc$sp_id, '\\]')
beta_JS <- paste0('beta\\[', JS_wc$sp_id, '\\]')
beta_FE <- paste0('beta\\[', FE_wc$sp_id, '\\]')
beta_FC <- paste0('beta\\[', FC_wc$sp_id, '\\]')

msd_JA <- paste0('mu_sp_d\\[', JA_wc$sp_id, '\\]')
msd_JM <- paste0('mu_sp_d\\[', JM_wc$sp_id, '\\]')
msd_JS <- paste0('mu_sp_d\\[', JS_wc$sp_id, '\\]')
msd_FE <- paste0('mu_sp_d\\[', FE_wc$sp_id, '\\]')
msd_FC <- paste0('mu_sp_d\\[', FC_wc$sp_id, '\\]')


beta_all <- c(beta_JA, beta_JM, beta_JS, beta_FE, beta_FC)
names <- c(as.character(JA_wc$COMMONNAME), as.character(JM_wc$COMMONNAME), 
           as.character(JS_wc$COMMONNAME), as.character(FE_wc$COMMONNAME), 
           as.character(FC_wc$COMMONNAME))

# MCMCvis::MCMCplot(fit, params = beta_all, 
#                   main = 'Effect of feather wear',
#                   ISB = FALSE,
#                   labels = names)

setwd('~/Desktop')
pdf('Figure_1.pdf', height = 11, width = 9, useDingbats = FALSE)
MCMCvis::MCMCplot(fit, params = 'beta', 
                  main = 'Effect of feather wear on wing chord',
                  rank = TRUE,
                  labels = wc_b$COMMONNAME,
                  sz_labels = 0.7)
dev.off()



msd_all <- c(msd_JA, msd_JM, msd_JS, msd_FE, msd_FC)

pdf('Figure_2.pdf', height = 11, width = 9, useDingbats = FALSE)
MCMCvis::MCMCplot(fit, params = msd_all, 
                  main = 'Difference in wing chord across age classes',
                  ISB = FALSE,
                  labels = names,
                  xlim = c(-3, 5),
                  sz_labels = 0.7)
dev.off()


NROW(JA_wc)
NROW(JM_wc)
NROW(JS_wc)
NROW(FE_wc)
NROW(FC_wc)




# mean diff density plots -------------------------------------------------

wc_JA <- dplyr::filter(wc_b, PSCORE == 'JA')
wc_JM <- dplyr::filter(wc_b, PSCORE == 'JM')
wc_JS <- dplyr::filter(wc_b, PSCORE == 'JS')
wc_FE <- dplyr::filter(wc_b, PSCORE == 'FE')
wc_FC <- dplyr::filter(wc_b, PSCORE == 'FC')


# plot(density(wc_JA$mu_sp_d), ylim = c(0, 2), lwd = 2)
# abline(v = mean(wc_JA$mu_sp_d), lwd = 3, lty = 2)
# 
# lines(density(wc_JM$mu_sp_d), col = 'red', lwd = 2)
# abline(v = mean(wc_JM$mu_sp_d), lwd = 3, col = 'red', lty = 2)
# 
# lines(density(wc_JS$mu_sp_d), col = 'green', lwd = 2)
# abline(v = mean(wc_JS$mu_sp_d), lwd = 3, col = 'green', lty = 2)
# 
# lines(density(wc_FE$mu_sp_d), col = 'blue', lwd = 2)
# abline(v = mean(wc_FE$mu_sp_d), lwd = 3, col = 'blue', lty = 2)
# 
# lines(density(wc_FC$mu_sp_d), col = 'pink', lwd = 2)
# abline(v = mean(wc_FC$mu_sp_d), lwd = 3, col = 'pink', lty = 2)


#chains for each PSCORE group
JA_gr <- mu_sp_d_ch[,wc_JA$sp_id]
JM_gr <- mu_sp_d_ch[,wc_JM$sp_id]
JS_gr <- mu_sp_d_ch[,wc_JS$sp_id]
FE_gr <- mu_sp_d_ch[,wc_FE$sp_id]
FC_gr <- mu_sp_d_ch[,wc_FC$sp_id][,-1] #remove red-wing blackbird - may be strange things going on here


#mean difference for PSCORE group
mn_JA <- apply(JA_gr, 1, mean)
mn_JM <- apply(JM_gr, 1, mean)
mn_JS <- apply(JS_gr, 1, mean)
mn_FE <- apply(FE_gr, 1, mean)
mn_FC <- apply(FC_gr, 1, mean)


# plot(density(mn_JA), ylim = c(0, 7), 
#      xlim = c(0.25, 1.75), lwd = 2,
#      main = 'Mean diff by PSCORE group',
#      xlab = 'Wing chord difference (mm)',
#      ylab = '',
#      yaxt = 'n')
# abline(v = mean(mn_JA), lwd = 3, lty = 2)
# 
# lines(density(mn_JM), col = 'red', lwd = 2)
# abline(v = mean(mn_JM), lwd = 3, col = 'red', lty = 2)
# 
# lines(density(mn_JS), col = 'green', lwd = 2)
# abline(v = mean(mn_JS), lwd = 3, col = 'green', lty = 2)
# 
# lines(density(mn_FE), col = 'blue', lwd = 2)
# abline(v = mean(mn_FE), lwd = 3, col = 'blue', lty = 2)
# 
# lines(density(mn_FC), col = 'brown', lwd = 2)
# abline(v = mean(mn_FC), lwd = 3, col = 'brown', lty = 2)

den_df <- rbind(data.frame(val = mn_JA, PSCORE = 'JA'), 
                data.frame(val = mn_JM, PSCORE = 'JM'),
                data.frame(val = mn_JS, PSCORE = 'JS'),
                data.frame(val = mn_FE, PSCORE = 'FE'),
                data.frame(val = mn_FC, PSCORE = 'FC'))

pdf('Figure_3.pdf', height = 5, width = 5, useDingbats = FALSE)
ggplot(den_df, aes(x = val, fill = PSCORE)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 12)) +
  ylab('') +
  ggtitle('Mean diff by PSCORE group') +
  xlab('Wing chord difference (mm)')
dev.off()



#scaled by abs wing chord
mean(wc_JA$sdiff)
sd(wc_JA$sdiff)
mean(wc_JM$sdiff)
sd(wc_JM$sdiff)
mean(wc_JS$sdiff)
sd(wc_JS$sdiff)
mean(wc_FE$sdiff)
sd(wc_FE$sdiff)
mean(wc_FC$sdiff)
sd(wc_FC$sdiff)



Ap <- dplyr::filter(maps_c2, sci_name == 'Agelaius phoeniceus')

uap <- unique(Ap$band_id)
for (i in 1:length(uap))
{
  #i <- 3
  temp <- dplyr::filter(Ap, band_id == uap[i])
  print(sd(temp$wing_chord))
}


