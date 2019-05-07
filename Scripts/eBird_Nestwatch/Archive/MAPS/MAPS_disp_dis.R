#####################
#Dispersal distance of birds (caught in more than one site)
#
# Take away: substantial dispersal is VERY RARE
#####################


library(dplyr)
library(geosphere)

dir <- '~/Google_Drive/R/'

#setwd(paste0(dir, '../'))
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))

maps_data <- readRDS('MAPS-age-filled.rds')

#find unique band ids
bid_cnt <- plyr::count(maps_data$band_id)

#individuals that have been captured at least 2 times
bid_c <- dplyr::filter(bid_cnt, freq >= 2)
c_birds <- bid_c[,1]

#only band_ids of interest
maps_c <- dplyr::filter(maps_data, band_id %in% c_birds)

#band ids
nid <- unique(maps_c$band_id)

N <- length(nid)
#N <- 500

mdf <- data.frame(common_name = rep(NA, N), 
                  sci_name = NA,
                  nyrs = NA,
                  band_id = NA,
                  sex = NA,
                  wing_chord = NA,
                  weight = NA,
                  dis = NA,
                  del = NA,
                  ka = NA) #ka = caught at birth grounds (first year)

pb <- txtProgressBar(min = 0, max = N, style = 3)
for (i in 1:N)
{
  #i <- 6177
  temp <- dplyr::filter(maps_c, band_id == nid[i])
  
  nyrs <- length(unique(temp$year))
  lls <- unique(temp[,c(1:2)])
  
  if (NROW(lls) > 1)
  {
    t_dis_mat <- matrix(data = NA, nrow = NROW(lls), ncol = NROW(lls))
    
    for (j in 1:NROW(lls))
    {
      #j <- 1
      for (k in j:NROW(lls))
      {
        #k <- 1
        t_dis_mat[j,k] <- geosphere::distm(c(lls$lng[j], lls$lat[j]), 
                                           c(lls$lng[k], lls$lat[k]))
      }
    }
    d_elev <- max(temp$elev) - min(temp$elev)
    #to display distance in km
    tm <- max(t_dis_mat, na.rm = TRUE) / 1000
  } else {
    tm <- 0
    d_elev <- 0
  }
  
  mdf$common_name[i] <- temp$common_name[1]
  mdf$sci_name[i] <- temp$sci_name[1]
  mdf$nyrs[i] <- nyrs
  mdf$band_id[i] <- temp$band_id[1]
  mdf$sex[i] <- temp$sex[1]
  mdf$wing_chord[i] <- median(temp$wing_chord, na.rm = TRUE)
  mdf$weight[i] <- median(temp$weight, na.rm = TRUE)
  mdf$dis[i] <- tm
  mdf$del[i] <- d_elev
  if (sum(temp$age %in% c(2,4)) > 0)
  {
    ka <- TRUE
  } else {
    ka <- FALSE
  }
  mdf$ka[i] <- ka 
  
  setTxtProgressBar(pb, i)
}

#only birds caught in more than one year
mdf_g1 <- dplyr::filter(mdf, nyrs > 1)

#remove problematic band string (according to Dani)
ndf <- mdf_g1[-which(mdf_g1$band_id >= 199172901, mdf_g1$band_id <= 199173000),]
range(ndf$dis, na.rm = TRUE)

ndf2 <- ndf[-which(ndf$dis == 0),]
hist(ndf2$dis)

#only species with more than 20 ind
sp_cnt <- plyr::count(ndf, 'sci_name')
sp20 <- sp_cnt[which(sp_cnt$freq >= 20), 1]
ndf3 <- dplyr::filter(ndf, sci_name %in% sp20)

mn_agg <- stats::aggregate(ndf3[, c('dis', 'del', 'wing_chord', 'weight')], list(ndf3$sci_name), mean)
colnames(mn_agg) <- c('sci_name', 'mn_dis', 'mn_del', 'mn_wing_chord', 'mn_weight')
med_agg <- stats::aggregate(ndf3[, c('dis', 'del', 'wing_chord', 'weight')], list(ndf3$sci_name), median)
colnames(med_agg) <- c('sci_name', 'med_dis', 'med_del', 'med_wing_chord', 'med_weight')
sd_agg <- stats::aggregate(ndf3[, c('dis', 'del', 'wing_chord', 'weight')], list(ndf3$sci_name), sd)
colnames(sd_agg) <- c('sci_name', 'sd_dis', 'sd_del', 'sd_wing_chord', 'sd_weight')
max_agg <- stats::aggregate(ndf3[, c('dis', 'del', 'wing_chord', 'weight')], list(ndf3$sci_name), max)
colnames(max_agg) <- c('sci_name', 'max_dis', 'max_del', 'max_wing_chord', 'max_weight')


j1 <- dplyr::left_join(mn_agg, med_agg, by = 'sci_name')
j2 <- dplyr::left_join(j1, sd_agg, by = 'sci_name')
j3 <- dplyr::left_join(j2, max_agg, by = 'sci_name')

plot(j3$mn_wing_chord, j3$max_dis)
plot(j3$mn_wing_chord, j3$sd_dis)
plot(j3$mn_weight, j3$max_dis)
plot(j3$mn_weight, j3$sd_dis)


hist(dplyr::filter(ndf3, ka == TRUE)$dis)
hist(dplyr::filter(ndf3, ka == FALSE)$dis)

par(mfrow = c(3,3))
for (i in 1:length(sp20))
{
  #i <- 1
  temp <- dplyr::filter(ndf3, sci_name == sp20[i])
  hist(temp$dis, main = sp20[i])
}
