##### Code to run simple models on posteror mean (or median?) parameters from arrival_master_2019-05-26.rds
setwd("/Users/JacobSocolar/Dropbox/Work/Phenomismatch/IAR_output")

IARout <- readRDS('arrival_master_2019-05-26.rds')
IARout$SC <- paste(IARout$species, IARout$cell, sep = '_')
IARout[IARout$species == 'Agelaius_phoeniceus' & IARout$cell == 425, ]


migD <- read.csv('migD.csv')  # gives whether long or short distance migrant, and a summary habitat code for each species
migD[is.na(migD)] <- 0

speed.df <- migD[, 1:4]
speed.df$mag <- speed.df$mbg <- speed.df$lspan <- NA
for(i in 1:nrow(speed.df)){
  speed.df$mag[i] <- unique(IARout$mean_alpha_gamma[IARout$species == speed.df$species[i]])
  speed.df$mbg[i] <- unique(IARout$mean_beta_gamma[IARout$species == speed.df$species[i]])
  speed.df$lspan[i] <- max(IARout$cell_lat[IARout$species == speed.df$species[i]]) - min(IARout$cell_lat[IARout$species == speed.df$species[i]])
}
speed.df[is.na(speed.df)] <- 0

# "naive intercept" is correlated with speed.  This is basically guaranteed by extrapolation back to zero.
plot(mbg ~ mag, data = speed.df)

# long-distance migrants tend to have estimates over a wider span of cells, making it tough to pick a single
# reference latitude to treat as the "intercept"
plot(lspan ~ LDM, data = speed.df)

# migration speed is strongly related to "naive intercept"
speed.model1 <- lm(mbg ~ mag, data = speed.df)
summary(speed.model1)

# speed is strongly related to whether a species is a long-distance migrant
speed.model2 <- lm(mbg ~ LDM, data = speed.df)
summary(speed.model2)

# as expected, migration distance doesn't fully account for the relationship between "naive intercept" and speed
speed.model3 <- lm(mbg ~ mag + LDM, data = speed.df)
summary(speed.model3)


speed2.df <- IARout[!duplicated(IARout$SC), ]
speed2.df$LDM <- NA
for(i in 1:nrow(speed2.df)){
  speed2.df$LDM[i] <- migD$LDM[which(migD$species == speed2.df$species[i])]
}

# a "non-naive" model, where arrival date in cells is modeled as a function of migration speed and random effects for
# both species and cell, still shows a relationship between arrival date and migration speed.
speed.model4 <- lmerTest::lmer(mean_gamma ~ mean_beta_gamma + (1 | species) + (1 | cell), data = speed2.df)
summary(speed.model4)

# this effect is essentially entirely due to who is a long- versus short-distance migrant
speed.model5 <- lmerTest::lmer(mean_gamma ~ mean_beta_gamma + LDM + (1 | species) + (1 | cell), data = speed2.df)
summary(speed.model5)

# below is still under construction
var.df <- migD[1:4, ]
var.df$mag <- var.df$mbg <- NA
for(i in 1:nrow(speed.df)){
  speed.df$mag[i] <- unique(IARout$mean_alpha_gamma[IARout$species == speed.df$species[i]])
  speed.df$mbg[i] <- unique(IARout$mean_beta_gamma[IARout$species == speed.df$species[i]])
  speed.df$lspan[i] <- max(IARout$cell_lat[IARout$species == speed.df$species[i]]) - min(IARout$cell_lat[IARout$species == speed.df$species[i]])
}






