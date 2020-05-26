
#check phenological intervals and change in BS over time for each species

#seems to be a large amount of variability in phenological intervals - how accurate are these data?


# interval ~ lat ----------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed/'))
nw <- readRDS('Nestwatch-pro-2019-08-27.rds')

cnt <- plyr::count(nw, 'species')
idx <- which(cnt$freq > 1000)
sp2 <- cnt[idx, 1]

nw2 <- dplyr::filter(nw, species %in% sp2)

#look at interval between lay, hatch, and fledge

out <- data.frame(species = rep(NA, length(sp2)), 
                  hf_slope = NA, 
                  hf_pval = NA,
                  hl_slope = NA,
                  hl_pval = NA,
                  lat_slope = NA, 
                  lat_pval = NA,
                  bs_slope = NA,
                  bs_pval = NA)

for (i in 1:length(sp2))
{
  print(round((i / length(sp2)),2))
  #i <- 18
  temp <- dplyr::filter(nw2, species == sp2[i])
  
  lf <- temp$fledge - temp$lay
  lh <- temp$hatch - temp$lay
  hf <- temp$fledge - temp$hatch
  
  #find which obs are not within 3 sd of mean and remove from df
  mn_f <- mean(temp$fledge, na.rm = TRUE)
  sd_f <- sd(temp$fledge, na.rm = TRUE)
  mn_l <- mean(temp$lay, na.rm = TRUE)
  sd_l <- sd(temp$lay, na.rm = TRUE)
  
  to.rm1 <- which(temp$fledge > (mn_f + (3 * sd_f)) | temp$fledge < (mn_f - (3 * sd_f)))
  to.rm2 <- which(temp$lay > (mn_l + (3 * sd_l)) | temp$lay < (mn_l - (3 * sd_l)))
  
  crm <- c(to.rm1, to.rm2)
  dups <- which(duplicated(crm))
  if (length(dups) > 0)
  {
    crm2 <- crm[-dups]
  } else {
    crm2 <- crm
  }
  
  if (length(crm2) > 0)
  {
    lf2 <- lf[-crm2]
    lh2 <- lh[-crm2]
    hf2 <- hf[-crm2]
    LAT <- temp$lat[-crm2]
    BS <- temp$n_fledged[-crm2] / temp$n_eggs[-crm2]
  } else {
    lf2 <- lf
    lh2 <- lh
    hf2 <- hf
    LAT <- temp$lat
    BS <- temp$n_fledged / temp$n_eggs
  }

  #change in lay hatch interval over lat
  if (sum(!is.na(lh2)) > 5 & 
      sum(diff(lh2), na.rm = TRUE) != 0 &
      sum(diff(LAT), na.rm = TRUE) != 0)
  {
    sfit <- summary(lm(lh2 ~ LAT))
    out$species[i] <- sp2[i]
    out$lh_slope[i] <- sfit$coefficients[2,1]
    out$lh_pval[i] <- sfit$coefficients[2,4]
    plot(LAT, lh2, main = paste0(sp2[i], ' - Lay -> Hatch'), pch = '.')
    abline(sfit, col = 'red', lwd = 2)
  }
  
  #change in hatch fledge interval over lat
  if (sum(!is.na(hf2)) > 5 & 
      sum(diff(hf2), na.rm = TRUE) != 0 &
      sum(diff(LAT), na.rm = TRUE) != 0)
  {
    sfit <- summary(lm(hf2 ~ LAT))
    out$species[i] <- sp2[i]
    out$hf_slope[i] <- sfit$coefficients[2,1]
    out$hf_pval[i] <- sfit$coefficients[2,4]
    plot(LAT, hf2, main = paste0(sp2[i], ' - Hatch -> Fledge'), pch = '.')
    abline(sfit, col = 'red', lwd = 2)
  }
}
hist(out$hf_slope)
hist(out$hf_pval)
hist(out$lh_slope)
hist(out$lh_pval)



plot(LAT, hf2, main = paste0(sp2[i], ' - Hatch -> Fledge'), pch = 19, 
     col = rgb(0,0,0,0.3), ylim = c(0, 40), xlim = c(30, 55),
     xlab = 'Latitude', ylab = 'Hatch -> Fledge Interval')

mn_hf2 <- mean(hf2, na.rm = TRUE)
sd_hf2 <- sd(hf2, na.rm = TRUE)

mn_hf2 + (1.96 * sd_hf2)
mn_hf2 - (1.96 * sd_hf2)



plot(LAT, lh2, main = paste0(sp2[i], ' - Lay -> Hatch'), pch = 19, 
     col = rgb(0,0,0,0.3), ylim = c(0, 40), xlim = c(30, 55),
     xlab = 'Latitude', ylab = 'Lay -> Hatch Interval')

mn_hl2 <- mean(hl2, na.rm = TRUE)
sd_hl2 <- sd(hl2, na.rm = TRUE)

mn_hl2 + (1.96 * sd_hl2)
mn_hl2 - (1.96 * sd_hl2)


# check nestwatch against what Jacob was using ----------------------------


setwd('~/Desktop/Bird_Phenology_Offline/Raw_data/Nestwatch/')

load('nest_data_formatted_ALL.Rdata')

j_nw <- cnd5@data
d1 <- dplyr::filter(j_nw, SPECIES_CODE == 'amerob')

#convert dates to jday
lay_jday <- as.numeric(format(as.Date(sapply(strsplit(as.character.factor(d1$FIRST_LAY_DT), 
                                                      split = ':'), '[', 1), format = '%d%b%Y'), '%j'))

hatch_jday <- as.numeric(format(as.Date(sapply(strsplit(as.character.factor(d1$HATCH_DT), 
                                                        split = ':'), '[', 1), format = '%d%b%Y'), '%j'))

fledge_jday <- as.numeric(format(as.Date(sapply(strsplit(as.character.factor(d1$FLEDGE_DT), 
                                                         split = ':'), '[', 1), format = '%d%b%Y'), '%j'))

hl <- hatch_jday - lay_jday
fh <- fledge_jday - hatch_jday

mean(hl, na.rm = TRUE)
sd(hl, na.rm = TRUE)
mean(fh, na.rm = TRUE)
sd(fh, na.rm = TRUE)


d2 <- dplyr::filter(nw, species == 'Turdus_migratorius')
hl2 <- d2$hatch - d2$lay
fh2 <- d2$fledge - d2$hatch

mean(hl2, na.rm = TRUE)
sd(hl2, na.rm = TRUE)
mean(fh2, na.rm = TRUE)
sd(fh2, na.rm = TRUE)

