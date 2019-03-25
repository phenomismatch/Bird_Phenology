#############################
# MAPS exploration
#
# *breeding onset with brood patch - see halfmax-bp (may need to aggregate to larger cells)
# *phenology as it relates to sex (difficult to do with MAPS data) - maybe with a few very late species?
# *phenology as it relates to age (difficult to do with MAPS data) - maybe with a few very late species?
# *weight/fat changes with age - CHANGING! Getting larger (wing chord/weight) but losing fat content and relative weight as they age
# *weight/fat changes with age of individual birds - CHANGING! Getting larger (wing chord/weight) but losing fat content and maybe relative weight as they age
# *weight/fat changes over time of population - CHANGING! Declines in weight and fat over time
# *weight/fat changes over time of year - weight complicated, fat yes, but complicated
# *males live longer than females - yes, supported by lit 
# *age distribution over time of population - yes, though due to birds of known age starting when project started
# *NEED TO EXPLORE: proportion (or frequency) of young hiting nets as metric of breeding phenology
#############################



# top-level dir --------------------------------------------------------------

dir <- '~/Google_Drive/R/'


# load packages -----------------------------------------------------------

library(dplyr)
library(dggridR)



# load data ---------------------------------------------------------------

setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps_data <- readRDS('MAPS_age_filled.rds')





# Only MAPS adults --------------------------------------------------------

maps_adults <- dplyr::filter(maps_data, age %in% c('1', '5', '6', '7', '8'))



# Onset of breeding with brood patch --------------------------------------

#breeding defined as:
#BP 3-4
#CP 1-3

#only captures that recorded brood patch/cloacal protuberance
maps_bp_f <- dplyr::filter(maps_adults, sex == 'F', brood_patch %in% c(0:5))
#maps_cp_m <- dplyr::filter(maps_adults, sex == 'M', cloacal_pro %in% c(0:3))

#add breeder col
maps_bp_f$bp <- NA

#brood patch maybe better than CP
plot(density(maps_bp_f[which(maps_bp_f$brood_patch <= 2 | maps_bp_f$brood_patch == 5),]$day))
lines(density(maps_bp_f[which(maps_bp_f$brood_patch > 2 & maps_bp_f$brood_patch < 5),]$day), 
      col = 'red')

maps_bp_f[which(maps_bp_f$brood_patch <= 2 | maps_bp_f$brood_patch == 5),]$bp <- 0
maps_bp_f[which(maps_bp_f$brood_patch > 2 & maps_bp_f$brood_patch < 5),]$bp <- 1


#how many obs for each species - top 10
cnt_bpf <- plyr::count(maps_bp_f, 'sci_name')

sp <- head(cnt_bpf[order(cnt_bpf$freq, decreasing = TRUE),], n = 10)$sci_name

cnt_spf <- plyr::count(maps_bp_f, c('cell', 'year'))

#see halfmax-bp for models




# phenology as it relates to sex ------------------------------------------

#DIFFICULT TO DO WITH MAPS DATA
#can't reliably measure because MAPS starts so late - after bird arrival, even though females are being captured in nets later (may be due to the fact that females do most of the incubation)

#difference between male and female arrival time

plot(density(maps_adults[which(maps_adults$sex == 'M'),]$day))
lines(density(maps_adults[which(maps_adults$sex == 'F'),]$day), col = 'red')

sci_names <- sort(unique(maps_adults$sci_name))

#species with some data that are relatively late arrivers:
#Contopus virens
#Empidonax alnorum
#Empidonax traillii
#Empidonax virescens
#Icteria virens
#Passerina cyanea
#Vireo olivaceus

for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 178
  temp <- dplyr::filter(maps_adults, sci_name == sci_names[i])
  
  dev.off()
  par(mfrow = c(3,3))
  t_yrs <- sort(unique(temp$year))
  for (j in 1:length(t_yrs))
  {
    #j <- 1
    temp2 <- dplyr::filter(temp, year == t_yrs[j])
    #plyr::count(temp2, c('cell'))
    if ((NROW(dplyr::filter(temp2, sex == 'F')) > 2) & (NROW(dplyr::filter(temp2, sex == 'M')) > 2))
    {
      plot(density(temp2[which(temp2$sex == 'M'),]$day), main = paste0(temp2$sci_name[1], ' - ', t_yrs[j]))
      lines(density(temp2[which(temp2$sex == 'F'),]$day), col = 'red')
    }
  }
}





# phenology as it relates to age ------------------------------------------

#DIFFICULT TO DO WITH MAPS DATA

#known age males
maps_age <- dplyr::filter(maps_adults, !is.na(true_age))

for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 178
  temp <- dplyr::filter(maps_age, sci_name == sci_names[i])
  
  # if (NROW(temp) > 30)
  # {
  #   sfit <- summary(lm(temp$day ~ temp$true_age))
  #   print(sfit$coef[2,4])
  # }

  dev.off()
  par(mfrow = c(3,3))
  t_yrs <- sort(unique(temp$year))
  for (j in 1:length(t_yrs))
  {
    #j <- 1
    temp2 <- dplyr::filter(temp, year == t_yrs[j])
    
    a1 <- dplyr::filter(temp2, true_age == 1)
    a2 <- dplyr::filter(temp2, true_age == 2)
    a3 <- dplyr::filter(temp2, true_age == 3)
    a4 <- dplyr::filter(temp2, true_age == 4)
    a5p <- dplyr::filter(temp2, true_age > 4)
    
    if (NROW(a1) > 1)
    {
      plot(density(a1$day), lwd = 2)
      
      if (NROW(a2) > 1)
      {
        lines(density(a2$day), col = 'red', lwd = 2)
      }
      if (NROW(a3) > 1)
      {
        lines(density(a3$day), col = 'orange', lwd = 2)
      }
      if (NROW(a4) > 1)
      {
        lines(density(a4$day), col = 'green', lwd = 2)
      }
      if (NROW(a5p) > 1)
      {
        lines(density(a5p$day), col = 'blue', lwd = 2)
      }
    }
  }
}



# how does weight/fat change with age -------------------------------------------

#getting larger but losing fat content and relative weight as they age

to.rm <- which(maps_adults$weight == 0 | maps_adults$wing_chord == 0 | is.na(maps_adults$weight) | is.na(maps_adults$fat_content))
maps_adults_qc <- maps_adults[-to.rm, ]

ma_qc <-  dplyr::filter(maps_adults_qc, !is.na(true_age))

#quantile regression with fat/weight and age

df_tt3 <- data.frame(sci_name = rep(NA, length(sci_names)),
                     c_name = NA,
                     N = NA,
                     slope_sweight = NA,
                     pv_sweight = NA,
                     slope_weight = NA,
                     pv_weight = NA,
                     slope_fat = NA,
                     pv_fat = NA,
                     slope_wc = NA,
                     pv_wc = NA)
dev.off()
par(mfrow = c(4,4))
for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 7
  temp <- dplyr::filter(ma_qc, sci_name == sci_names[i])
  
  df_tt3$sci_name[i] <- temp$sci_name[1]
  df_tt3$c_name[i] <- temp$common_name[1]
  df_tt3$N[i] <- NROW(temp)
  
  if (length(unique(temp$true_age)) > 2)
  {
    tfit <- summary(lm((weight/wing_chord) ~ true_age, data = temp))
    #plot(temp$true_age, (temp$weight/temp$wing_chord), main = paste0('Weight ', temp$sci_name[1]))
    df_tt3$slope_sweight[i] <- round(tfit$coef[2,1], 3)
    df_tt3$pv_sweight[i] <- round(tfit$coef[2,4], 3)
    
    tfit15 <- summary(lm(weight ~ true_age, data = temp))
    #plot(temp$true_age, (temp$weight), main = paste0('Weight ', temp$sci_name[1]))
    df_tt3$slope_weight[i] <- round(tfit15$coef[2,1], 3)
    df_tt3$pv_weight[i] <- round(tfit15$coef[2,4], 3)
    
    tfit2 <- summary(lm(fat_content ~ true_age, data = temp))
    #plot(temp$true_age, temp$fat_content, main = paste0('Fat ', temp$sci_name[1]))
    df_tt3$slope_fat[i] <- round(tfit2$coef[2,1], 3)
    df_tt3$pv_fat[i] <- round(tfit2$coef[2,4], 3)
    
    tfit3 <- summary(lm(wing_chord ~ true_age, data = temp))
    plot(temp$true_age, temp$wing_chord, main = paste0('Wing Chord ', temp$sci_name[1]))
    abline(tfit3, col = 'red')
    df_tt3$slope_wc[i] <- round(tfit3$coef[2,1], 3)
    df_tt3$pv_wc[i] <- round(tfit3$coef[2,4], 3)
  }
}

ntt3 <- dplyr::filter(df_tt3, N > 200)

#losing standardized weight slightly as they age
hist(ntt3$slope_sweight)
hist(ntt3$pv_sweight)

#overall gaining weight
hist(ntt3$slope_weight)
hist(ntt3$pv_weight)

#losing fat content
hist(ntt3$slope_fat)
hist(ntt3$pv_fat)

#and wing chords growing
hist(ntt3$slope_wc)
hist(ntt3$pv_wc)





# how does weight/fat of poplation change over time --------------------

#*is weight (std by wing chord) changing over time? DECLINING
#*is fat changing over time? DECLINING
#*is wing chord changing over time? relatively stable


df_tt2 <- data.frame(sci_name = rep(NA, length(sci_names)),
                    c_name = NA,
                    N = NA,
                    slope_sweight = NA,
                    pv_sweight = NA,
                    slope_weight = NA,
                    pv_weight = NA,
                    slope_fat = NA,
                    pv_fat = NA,
                    slope_wc = NA,
                    pv_wc = NA)
dev.off()
par(mfrow = c(4,4))
for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 55
  temp <- dplyr::filter(maps_adults_qc, sci_name == sci_names[i])
  
  df_tt2$sci_name[i] <- temp$sci_name[1]
  df_tt2$c_name[i] <- temp$common_name[1]
  df_tt2$N[i] <- NROW(temp)
  
  if (length(unique(temp$year)) > 3)
  {
    tfit <- summary(lm((weight/wing_chord) ~ year, data = temp))
    plot(temp$year, (temp$weight/temp$wing_chord), main = paste0('Weight ', temp$sci_name[1]))
    abline(tfit, col = 'red')
    df_tt2$slope_sweight[i] <- round(tfit$coef[2,1], 3)
    df_tt2$pv_sweight[i] <- round(tfit$coef[2,4], 3)
    
    tfit15 <- summary(lm(weight ~ year, data = temp))
    plot(temp$year, temp$weight, main = paste0('Weight ', temp$sci_name[1]))
    abline(tfit15, col = 'red')
    df_tt2$slope_weight[i] <- round(tfit15$coef[2,1], 3)
    df_tt2$pv_weight[i] <- round(tfit15$coef[2,4], 3)
    
    tfit2 <- summary(lm(fat_content ~ year, data = temp))
    plot(temp$year, temp$fat_content, main = paste0('Fat ', temp$sci_name[1]))
    abline(tfit2, col = 'red')
    df_tt2$slope_fat[i] <- round(tfit2$coef[2,1], 3)
    df_tt2$pv_fat[i] <- round(tfit2$coef[2,4], 3)

    tfit3 <- summary(lm(wing_chord ~ year, data = temp))
    plot(temp$year, temp$wing_chord, main = paste0('Wing Chord ', temp$sci_name[1]))
    abline(tfit3, col = 'red')
    df_tt2$slope_wc[i] <- round(tfit3$coef[2,1], 3)
    df_tt2$pv_wc[i] <- round(tfit3$coef[2,4], 3)
  }
}

ntt2 <- dplyr::filter(df_tt2, N > 200)

#declining overall weight
hist(ntt2$slope_sweight)
hist(ntt2$pv_sweight)

#slight loss in weight
hist(ntt2$slope_weight)
hist(ntt2$pv_weight)

#losing fat over time
hist(ntt2$slope_fat)
hist(ntt2$pv_fat)

#some change in wc, but both ways
hist(ntt2$slope_wc)
hist(ntt2$pv_wc)


# how does weight/fat of poplation change over lat --------------------

#*is weight (std by wing chord) changing over lat? slight decrease in weight with higher lat
#*is weight (non std) changing over lat? slight increse at high lat for most species (decrease for some)
#*is fat changing over lat? decrease in fat with higher lat
#*is wing chord changing over lat? increse in wing chord with higher lat

df_tt5 <- data.frame(sci_name = rep(NA, length(sci_names)),
                     c_name = NA,
                     N = NA,
                     slope_sweight = NA,
                     pv_sweight = NA,
                     slope_weight = NA,
                     pv_weight = NA,
                     slope_fat = NA,
                     pv_fat = NA,
                     slope_wc = NA,
                     pv_wc = NA)
#dev.off()
#par(mfrow = c(4,4))
for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 55
  temp <- dplyr::filter(maps_adults_qc, sci_name == sci_names[i])
  
  df_tt5$sci_name[i] <- temp$sci_name[1]
  df_tt5$c_name[i] <- temp$common_name[1]
  df_tt5$N[i] <- NROW(temp)
  
  if (length(unique(temp$lat)) > 3)
  {
    tfit <- summary(lm((weight/wing_chord) ~ lat, data = temp))
    #plot(temp$lat, (temp$weight/temp$wing_chord), main = paste0('Weight ', temp$sci_name[1]))
    #abline(tfit, col = 'red')
    df_tt5$slope_sweight[i] <- round(tfit$coef[2,1], 3)
    df_tt5$pv_sweight[i] <- round(tfit$coef[2,4], 3)
    
    tfit15 <- summary(lm(weight ~ lat, data = temp))
    #plot(temp$lat, (temp$weight), main = paste0('Weight ', temp$sci_name[1]))
    #abline(tfit, col = 'red')
    df_tt5$slope_weight[i] <- round(tfit15$coef[2,1], 3)
    df_tt5$pv_weight[i] <- round(tfit15$coef[2,4], 3)
    
    tfit2 <- summary(lm(fat_content ~ lat, data = temp))
    #plot(temp$lat, temp$fat_content, main = paste0('Fat ', temp$sci_name[1]))
    #abline(tfit2, col = 'red')
    df_tt5$slope_fat[i] <- round(tfit2$coef[2,1], 3)
    df_tt5$pv_fat[i] <- round(tfit2$coef[2,4], 3)
    
    tfit3 <- summary(lm(wing_chord ~ lat, data = temp))
    #plot(temp$lat, temp$wing_chord, main = paste0('Wing Chord ', temp$sci_name[1]))
    #abline(tfit3, col = 'red')
    df_tt5$slope_wc[i] <- round(tfit3$coef[2,1], 3)
    df_tt5$pv_wc[i] <- round(tfit3$coef[2,4], 3)
  }
}

ntt5 <- dplyr::filter(df_tt5, N > 200)

#no pronounced changes
hist(ntt5$slope_sweight)
hist(ntt5$pv_sweight)

hist(ntt5$slope_weight)
hist(ntt5$pv_weight)


#deccrease in fat as moving to higher lats
hist(ntt5$slope_fat)
hist(ntt5$pv_fat)

#larger wing chord moving north
hist(ntt5$slope_wc)
hist(ntt5$pv_wc)



# how does weight/fat of population change over a year --------------------

#COMPLICATED
#fat scores change over season - differs and is complex, though
#weight more complicated


df_tt <- data.frame(sci_name = rep(NA, length(sci_names)),
                       c_name = NA,
                       N = NA,
                       slope_weight = NA,
                       pv_weight = NA,
                       slope_fat = NA,
                       pv_fat = NA,
                       slope_wf = NA,
                       pv_wf = NA)
dev.off()
par(mfrow = c(4,4))
for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 2
  temp <- dplyr::filter(maps_adults_qc, sci_name == sci_names[i])
  
  df_tt$sci_name[i] <- temp$sci_name[1]
  df_tt$c_name[i] <- temp$common_name[1]
  df_tt$N[i] <- NROW(temp)
  
  if (NROW(temp) > 3)
  {
    tfit <- summary(lm((weight/wing_chord) ~ day, data = temp))
    plot(temp$day, (temp$weight/temp$wing_chord), main = paste0('Weight ', temp$sci_name[1]))
    df_tt$slope_weight[i] <- round(tfit$coef[2,1], 3)
    df_tt$pv_weight[i] <- round(tfit$coef[2,4], 3)
    
    tfit2 <- summary(lm(fat_content ~ day, data = temp))
    plot(temp$day, temp$fat_content, main = paste0('Fat ', temp$sci_name[1]))
    df_tt$slope_fat[i] <- round(tfit2$coef[2,1], 3)
    df_tt$pv_fat[i] <- round(tfit2$coef[2,4], 3)
    
    tfit3 <- summary(lm(fat_content ~ (weight/wing_chord), data = temp))
    plot((temp$weight/temp$wing_chord), temp$fat_content, main = paste0('WF ', temp$sci_name[1]))
    df_tt$slope_wf[i] <- round(tfit3$coef[2,1], 3)
    df_tt$pv_wf[i] <- round(tfit3$coef[2,4], 3)
  }
}
  




# how does weight/fat/wing_chord of individual birds change over time -----------


#DECREASE in fat, INCREASE wing chord, INCREASE weight

#find unique band ids
bid_cnt <- plyr::count(maps_adults_qc$band_id)

#individuals that have been captured more than 5 times
bid_c <- dplyr::filter(bid_cnt, freq > 3)
c_birds <- bid_c[,1]


maps_c <- dplyr::filter(maps_adults_qc, band_id %in% c_birds, !is.na(true_age))

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


df_temp <- data.frame(band_id = rep(NA, length(c_birds)),
                      sci_name = rep(NA, length(c_birds)),
                      slope_weight = rep(NA, length(c_birds)), 
                      pv_weight = rep(NA, length(c_birds)),
                      slope_sweight = rep(NA, length(c_birds)), 
                      pv_sweight = rep(NA, length(c_birds)),
                      slope_fat = rep(NA, length(c_birds)), 
                      pv_fat = rep(NA, length(c_birds)),
                      slope_wc = rep(NA, length(c_birds)), 
                      pv_wc = rep(NA, length(c_birds)),
                      slope_wfc = rep(NA, length(c_birds)), 
                      pv_wfc = rep(NA, length(c_birds)))

for (i in 1:length(c_birds))
{
  #i <- 3100
  temp <- dplyr::filter(maps_c, band_id == c_birds[i])
  
  df_temp$band_id[i] <- temp$band_id[1]
  df_temp$sci_name[i] <- temp$sci_name[1]
  
  if (length(unique(temp$true_age)) > 1)
  {
    #weight
    tfit <- summary(lm(weight ~ true_age, data = temp))
    df_temp$slope_weight[i] <- round(tfit$coef[2,1], 2)
    df_temp$pv_weight[i] <- round(tfit$coef[2,4], 2)
    
    tfit15 <- summary(lm((weight/wing_chord) ~ true_age, data = temp))
    df_temp$slope_sweight[i] <- round(tfit15$coef[2,1], 2)
    df_temp$pv_sweight[i] <- round(tfit15$coef[2,4], 2)
    
    #fat
    tfit2 <- summary(lm(fat_content ~ true_age, data = temp))
    df_temp$slope_fat[i] <- round(tfit2$coef[2,1], 2)
    df_temp$pv_fat[i] <- round(tfit2$coef[2,4], 2)
    
    #wing chord
    tfit3 <- summary(lm(wing_chord ~ true_age, data = temp))
    df_temp$slope_wc[i] <- round(tfit3$coef[2,1], 2)
    df_temp$pv_wc[i] <- round(tfit3$coef[2,4], 2)
    
    #weight ~ fat content
    if (length(unique(temp$fat_content)) > 1)
    {
      tfit4 <- summary(lm(weight ~ fat_content, data = temp))
      df_temp$slope_wfc[i] <- round(tfit4$coef[2,1], 2)
      df_temp$pv_wfc[i] <- round(tfit4$coef[2,4], 2)
    }
  }
}

hist(df_temp$slope_weight)
hist(df_temp$slope_sweight)
hist(df_temp$pv_weight)
hist(df_temp$slope_fat)
hist(df_temp$slope_wc)
hist(df_temp$slope_wfc)
hist(df_temp$pv_wfc)

hist(dplyr::filter(df_temp, pv_weight < 0.05)$slope_weight)
hist(dplyr::filter(df_temp, pv_weight < 0.05)$slope_sweight)


# do females live longer than males? max life span each species -----------------------------------------

#MALES LIVE LONGER - known in lit
#max ages of each species calculated

maps_as <- dplyr::filter(maps_age, sex %in% c('M', 'F'))
maps_m <- dplyr::filter(maps_age, sex %in% c('M'))
maps_f <- dplyr::filter(maps_age, sex %in% c('F'))

#30% more males captures than females
#(NROW(maps_m) - NROW(maps_f)) / NROW(maps_f)

hist(maps_m$true_age)
hist(maps_f$true_age)

mean(maps_m$true_age)
sd(maps_m$true_age)
mean(maps_f$true_age)
sd(maps_f$true_age)



df_temp2 <- data.frame(sci_name = rep(NA, length(sci_names)),
                      c_name = NA,
                      N = NA,
                      max_age = NA,
                      mn_male_age = NA,
                      sd_male_age = NA,
                      mn_female_age = NA,
                      sd_female_age = NA)

for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 142
  
  temp <- dplyr::filter(maps_age, sci_name == sci_names[i])
  df_temp2$sci_name[i] <- temp$sci_name[1]
  df_temp2$c_name[i] <- temp$common_name[1]
  df_temp2$N[i] <- NROW(temp)
  df_temp2$max_age[i] <- max(temp$true_age)
  
  temp_m <- dplyr::filter(temp, sex %in% c('M'))
  temp_f <- dplyr::filter(temp, sex %in% c('F'))
  
  if (NROW(temp) > 1)
  {
    df_temp2$mn_male_age[i] <- round(mean(temp_m$true_age, na.rm = TRUE), 2)
    df_temp2$sd_male_age[i] <- round(sd(temp_m$true_age, na.rm = TRUE), 2)
    df_temp2$mn_female_age[i] <- round(mean(temp_f$true_age, na.rm = TRUE), 2)
    df_temp2$sd_female_age[i] <- round(sd(temp_f$true_age, na.rm = TRUE), 2)
  }
}

to.rm <- which(is.na(df_temp2$sci_name))
df_temp3 <- df_temp2[-to.rm,]

#Blackpol warbler
df_temp2[grep('Shrike', df_temp2$c_name),]





# has age distribution of birds changed over time? ------------------------

#AGE OF KNOWN AGE BIRDS INCREASING - likely due to true age being calculated from birds that were tagged as juveniles

df_temp4 <- data.frame(sci_name = rep(NA, length(sci_names)),
                       c_name = NA,
                       N = NA,
                       slope = NA,
                       pv = NA)

for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 3
  
  temp <- dplyr::filter(maps_age, sci_name == sci_names[i])
  df_temp4$sci_name[i] <- temp$sci_name[1]
  df_temp4$c_name[i] <- temp$common_name[1]
  df_temp4$N[i] <- NROW(temp)
  
  t_yrs <- unique(temp$year)
  
  if (NROW(temp) > 3)
  {
    tdf <- data.frame(year = rep(NA, length(t_yrs)), mn_age = NA, sd_age = NA)
    for (j in 1:length(t_yrs))
    {
      #j <- 1
      temp2 <- dplyr::filter(temp, year == t_yrs[j])
  
      tdf$year[j] <- temp2$year[1]
      tdf$mn_age[j] <- mean(temp2$true_age)
      tdf$sd_age[j] <- sd(temp2$true_age)
    }
    
    tfit <- summary(lm(mn_age ~ year, data = tdf))
    df_temp4$slope[i] <- round(tfit$coef[2,1], 2)
    df_temp4$pv[i] <- round(tfit$coef[2,4], 2)
  }
}

df_temp5 <- dplyr::filter(df_temp4, N > 50)

NROW(df_temp5)


# How does age change over lat --------------------------------------------

#a number of sig trends, but could be due to sampling bias (perhaps more southerly sites sampling less over time?)
df_age_lat <- data.frame(sci_name = rep(NA, length(sci_names)),
                       c_name = NA,
                       N = NA,
                       slope = NA,
                       pv = NA)

for (i in 1:length(sci_names))
{
  #which(sci_names == 'Vireo olivaceus')
  #i <- 3
  
  temp <- dplyr::filter(maps_age, sci_name == sci_names[i])
  df_age_lat$sci_name[i] <- temp$sci_name[1]
  df_age_lat$c_name[i] <- temp$common_name[1]
  df_age_lat$N[i] <- NROW(temp)
  
  t_lat <- unique(temp$lat)
  
  if (NROW(temp) > 50)
  {
    tfit <- summary(lm(true_age ~ lat, data = temp))
    df_age_lat$slope[i] <- round(tfit$coef[2,1], 2)
    df_age_lat$pv[i] <- round(tfit$coef[2,4], 2)
  }
}

hist(df_age_lat$slope)
hist(df_age_lat$pv)

  

