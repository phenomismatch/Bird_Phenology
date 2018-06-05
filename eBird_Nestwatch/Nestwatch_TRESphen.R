library(dggridR)
library(lubridate)
library(doBy)

hexgrid6 <- dgconstruct(res=6)

#### 1) Nestwatch data import, inspection, and manipulation ####
nestwatch <- read.csv("/Users/TingleyLab/Dropbox/Work/Nestwatch/AllLocations.csv")

# Extract only records with known first lay date, then extract jday and year
knest <- nestwatch[nestwatch$FIRST_LAY_DT != nestwatch$FIRST_LAY_DT[1],]
knest$date <- as.Date(substr(knest$FIRST_LAY_DT,1,9),format="%d%B%Y")
knest$jday <- yday(knest$date)
knest$year <- year(knest$date)

# hexagonal grid
knest$cell6 <- dgGEO_to_SEQNUM(hexgrid6, knest$LONGITUDE, nestwatch$LATITUDE)$seqnum

TRES <- knest[which(knest$SPECIES_CODE == "treswa"), ]

hex6_nestwatch <- summaryBy(jday ~ year + cell6, data = TRES)
hex6_nestwatch$n <- NA
for(i in 1:dim(hex6_nestwatch)[1]){
  hex6_nestwatch$n[i] <- length(which(TRES$year == hex6_nestwatch$year[i] & TRES$cell6 == hex6_nestwatch$cell6[i]))
}
names(hex6_nestwatch) <- c("year", "cell", "jday.mean", "n")
save(hex6_nestwatch, file = "/Users/TingleyLab/Dropbox/Work/Phenomismatch/NA_birdPhen/hex6_nestwatch.Rdata")
