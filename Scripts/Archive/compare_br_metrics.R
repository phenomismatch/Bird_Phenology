#compare breeding phenology metrics


# Top-level dir -----------------------------------------------------------

#desktop/laptop
dir <- '~/Google_Drive/R/'


# load MAPS data ------------------------------------------------------------

#read in - RDS create with 1-query-db.R in wing_chord_changes project
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps <- readRDS('MAPS-age-filled.rds')
colnames(maps)[grep('sci_name', colnames(maps))] <- 'species'
#add underscore to species naems
maps$species <- gsub(' ', '_', maps$species)



# load logistic and GAM ---------------------------------------------------

#GAM
setwd(paste0(dir, 'Bird_Phenology/Data/Processed/br_arr_2019-08-26'))
d26 <- readRDS('juv-output-2019-08-26.rds')


# load nestwatch ----------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
nw <- readRDS('Nestwatch-pro-2019-08-27.rds')
NROW(nw)
agg_nw <- aggregate(fledge ~ species + cell + year, FUN = mean, data = nw)
agg_N <- aggregate(fledge ~ species + cell + year, FUN = length, data = nw)
agg_nw_sd <- aggregate(fledge ~ species + cell + year, FUN = sd, data = nw)

#add sd
agg_nw2 <- data.frame(agg_nw, sd = agg_nw_sd[,4], N = agg_N[,4])
agg_nw3 <- dplyr::filter(agg_nw2, N > 10)


# merge with Nestwatch ----------------------------------------------------

mrg <- dplyr::inner_join(d26, agg_nw3, c('species', 'year', 'cell'))
mrg_f <- dplyr::filter(mrg, !is.na(juv_mean), !is.na(fledge))

plot(mrg_f$juv_mean, mrg_f$fledge)
summary(lm(fledge ~ juv_mean, data = mrg_f))


#####BELOW NOT SETUP FOR NESTWATCH DATA#######

# analyze -----------------------------------------------------------------

df <- data.frame(species = rep(NA, NROW(d21)), cell = NA, year = NA, 
                 first = NA, first_sd = NA, GAM = NA, GAM_sd = NA, 
                 LOG = NA, LOG_sd = NA)
counter <- 1
usp <- unique(d21$species)
for (i in 1:length(usp))
{
  #i <- 1
  d21_f <- dplyr::filter(d21, species == usp[i])
  d22_f <- dplyr::filter(d22, species == usp[i])
  maps_f <- dplyr::filter(maps, species == usp[i])  
  
  uyr <- unique(d21_f$year)
  for (j in 1:length(uyr))
  {
    #j <- 1
    d21_f2 <- dplyr::filter(d21_f, year == uyr[j])
    d22_f2 <- dplyr::filter(d22_f, year == uyr[j])
    maps_f2 <- dplyr::filter(maps_f, year == uyr[j])  
    
    uc <- unique(d21_f2$cell)
    for (k in 1:length(uc))
    {
      #k <- 1
      d21_f3 <- dplyr::filter(d21_f2, cell == uc[k])
      d22_f3 <- dplyr::filter(d22_f2, cell == uc[k])
      maps_f3 <- dplyr::filter(maps_f2, cell == uc[k])
      
      if (NROW(maps_f3) > 0 & !is.na(d21_f3$juv_mean))
      {
        ubid <- unique(maps_f3$band_id)
        out <- rep(NA, length(ubid))
        for (n in 1:length(ubid))
        {
          maps_f4 <- dplyr::filter(maps_f3, band_id == ubid[n])
          out[n] <- min(maps_f4$day)
        }
        
        df$species[counter] <- usp[i]
        df$year[counter] <- uyr[j]
        df$cell[counter] <- uc[k]
        df$first[counter] <- mean(out)
        df$first_sd[counter] <- sd(out)
        df$GAM[counter] <- d21_f3$juv_mean
        df$GAM_sd[counter] <- d21_f3$juv_sd
        df$LOG[counter] <- d22_f3$juv_mean
        df$LOG_sd[counter] <- d22_f3$juv_sd
        counter <- counter + 1
      }
    }
  }
}

#mean first to GAM
plot(df$first, df$GAM)
df2 <- df[is.finite(df$GAM) & is.finite(df$first),]
summary(lm(df2$GAM ~ df2$first))
#mean first to logistic
plot(df$first, df$LOG)
#GAM to logistic
plot(df$GAM, df$LOG)
abline(0, 1, col = 'red', lty = 2, lwd = 2)

