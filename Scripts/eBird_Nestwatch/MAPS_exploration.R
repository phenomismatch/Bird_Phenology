#############################
#MAPS exploration
#
#
#############################




#--------------------------------#
#logistic model for presence of brood patch


#What characterizes breeding? Brood patch?


#just hatch year birds
MAPS_f <- dplyr::filter(MAPS_mrg5, True_age == 0)
plyr::count(MAPS_f, 'Sci_name')

oc_0 <- filter(MAPS_f, Sci_name == 'Oreothlypis celata')
plyr::count(oc_0, c('Cell', 'Year'))
tt <- dplyr::filter(oc_0, Cell == 444, Year == 2011)

plot(density(tt$Jday))


#--------------------------------#


#--------------------------------#
#logistic model for presence of brood patch

MAPS_age <- dplyr::filter(MAPS_mrg5, True_age > 0)
plyr::count(MAPS_age, 'Sci_name')

C_ustulatus <- filter(MAPS_age, Sci_name == 'Catharus ustulatus')

#breeding defined as:
#BP 2-4
#CP 1-3

C_ustulatus <- dplyr::filter(MAPS_age, Sci_name == 'Catharus ustulatus')

cu_f <- dplyr::filter(C_ustulatus, Sex == 'F')
cu_m <- dplyr::filter(C_ustulatus, Sex == 'M')

#add breeder col
C_ustulatus$BR <- NA
cu_f$BR <- NA

#male female
plot(density(C_ustulatus[which(C_ustulatus$Sex == 'M'),]$Jday))
lines(density(C_ustulatus[which(C_ustulatus$Sex == 'F'),]$Jday), col = 'red')

#cloacal protuberance probably not a good metric
plot(density(cu_m[which(cu_m$Cloacal_prot <= 2),]$Jday))
lines(density(cu_m[which(cu_m$Cloacal_prot > 2),]$Jday), col = 'red')

#brood patch maybe better
plot(density(cu_f[which(cu_f$Brood_patch <= 2 | cu_f$Brood_patch == 5),]$Jday))
lines(density(cu_f[which(cu_f$Brood_patch > 2 & cu_f$Brood_patch < 5),]$Jday), 
      col = 'red')

cu_f[which(cu_f$Brood_patch > 2 & cu_f$Brood_patch < 5),]$BR <- 1
cu_f[which(cu_f$Brood_patch <= 2 | cu_f$Brood_patch == 5),]$BR <- 0

#DATA <- dplyr::filter(cu_f, Jday < DAY, Jday > 50)
DATA <- cu_f

DATA$sjday <- DATA$Jday
DATA$sjday2 <- DATA$Jday^2
DATA$sjday3 <- DATA$Jday^3

plyr::count(cu_f, c('Cell', 'Year'))


library(rstanarm)
ITER <- 1000
#ITER <- 10
CHAINS <- 3

fit2 <- rstanarm::stan_glm(BR ~ sjday + sjday2 + sjday3,
                           data = DATA,
                           family = binomial(link = "logit"),
                           algorithm = 'sampling',
                           iter = ITER,
                           chains = CHAINS,
                           cores = CHAINS)

summary(fit2)

predictDays <- range(DATA$sjday)[1]:range(DATA$sjday)[2]
predictDays2 <- predictDays^2
predictDays3 <- predictDays^3

newdata <- data.frame(sjday = predictDays, sjday2 = predictDays2, sjday3 = predictDays3)

dfit <- rstanarm::posterior_linpred(fit2, newdata = newdata, transform = T)
halfmax_fit <- rep(NA, ((ITER/2)*CHAINS))

for (L in 1:((ITER/2)*CHAINS))
{
  #L <- 1
  rowL <- as.vector(dfit[L,])
  halfmax_fit[L] <- predictDays[min(which(rowL > (max(rowL)/2)))]
}


mn_dfit <- apply(dfit, 2, mean)
LCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.025))
UCI_dfit <- apply(dfit, 2, function(x) quantile(x, probs = 0.975))
mn_hm <- mean(halfmax_fit)
LCI_hm <- quantile(halfmax_fit, probs = 0.025)
UCI_hm <- quantile(halfmax_fit, probs = 0.975)

plot(predictDays, UCI_dfit, type = 'l', col = 'red', lty = 2, lwd = 2,
     ylim = c(0, max(UCI_dfit)),
     xlab = 'Julian Day', ylab = 'Detection Probability')
lines(predictDays, LCI_dfit, col = 'red', lty = 2, lwd = 2)
lines(predictDays, mn_dfit, lwd = 2)
DATA$BR_pl <- NA
DATA[which(DATA$BR == 1),]$BR_pl <- max(UCI_dfit)
DATA[which(DATA$BR == 0),]$BR_pl <- 0
points(DATA$Jday, DATA$BR_pl, col = rgb(0,0,0,0.25))
abline(v = mn_hm, col = rgb(0,0,1,0.5), lwd = 2)
abline(v = LCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
abline(v = UCI_hm, col = rgb(0,0,1,0.5), lwd = 2, lty = 2)
mean(halfmax_fit)
sd(halfmax_fit)



#--------------------------------#


# 
# jfit <- lm(MAPS_age$Jday ~ MAPS_age$True_age)
# summary(jfit)
# 
# plot(MAPS_age$True_age, MAPS_age$Jday, pch = '.', col = rgb(0,0,0, 0.5))
# abline(jfit, col = 'red')
# a1 <- dplyr::filter(MAPS_age, True_age == 1)
# a2 <- dplyr::filter(MAPS_age, True_age == 2)
# a3 <- dplyr::filter(MAPS_age, True_age == 3)
# a4 <- dplyr::filter(MAPS_age, True_age == 4)
# a5 <- dplyr::filter(MAPS_age, True_age == 5)
# plot(density(a1$Jday))
# lines(density(a2$Jday))
# lines(density(a3$Jday))
# lines(density(a4$Jday))
# lines(density(a5$Jday))




#*breeding timing
#Each species/cell/year:
#maybe each capture different row, 1 for brood patch/CP, 0 if not -> fit a logistic to this for each cell/year?
#y[i] ~ bern(p[i])
#logit(p[i]) = alpha + beta * DAY[i]

#*phenology and age
#could fit logistic regression with age as random effect -> each age class would have different beta param
#for each cell/year: predict start of breeding for each age class (age 1, age 2, age 3+)
#calculate derived qty - difference in halfmax estimates
#y[i] ~ bern(p[i])
#logit(p[i]) = alpha[id[i]] + beta[id[i]] * DAY[i]

#*difference in arrival between males and females - how this maps onto phlogeny/traits


#*fat content and age - Fat[i] ~ alpha[id[i]] + beta1[id[i]] * Age[i] + beta2[id[i]] * Day[i] + beta3[id[i]] * Day[i]^2
plot(MAPS_age$True_age, MAPS_age$Fat_content, col = rgb(0,0,0, 0.5))
plot(MAPS_age$Jday, MAPS_age$Fat_content, col = rgb(0,0,0, 0.5))
summary(lm(Fat_content ~ True_age + poly(Jday, 2, raw = TRUE), data = MAPS_age))

#*Weight and age - Weight[i] ~ alpha[id[i]] + beta1[id[i]] * Age[i] + beta2[id[i]] * Day[i] + beta3[id[i]] * Day[i]^2
plot(MAPS_age$True_age, MAPS_age$Weight, col = rgb(0,0,0, 0.5))
plot(MAPS_age$Jday, MAPS_age$Weight, col = rgb(0,0,0, 0.5))
summary(lm(Weight ~ True_age + poly(Jday, 2, raw = TRUE), data = MAPS_age))



#filter by species, get breeding date (breeding period in this case) for each year/cell
na_reps <- rep(NA, NROW(unique(MAPS_mrg2[c('YR', 'cell', 'SCINAME')])))

MAPS_out <- data.frame(YR = na_reps,
                       cell = na_reps,
                       SCINAME = na_reps,
                       midpoint = na_reps,
                       l_bounds = na_reps,
                       u_bounds = na_reps,
                       n_stations = na_reps)

counter <- 1
for (i in 1:length(species_list_i2))
{
  #i <- 19
  temp <- dplyr::filter(MAPS_mrg2, SCINAME == species_list_i2[i])
  if (NROW(temp) > 0)
  {
    
    t_cell <- unique(temp$cell)
    
    for (j in 1:length(t_cell))
    {
      #j <- 4
      temp2 <- dplyr::filter(temp, cell == t_cell[j])
      t_yr <- unique(temp2$YR)
      
      for (k in 1:length(t_yr))
      {
        
        print(paste0('species: ', species_list_i2[i], ', ',
                     'cell: ', t_cell[j], ', ',
                     'year: ', t_yr[k]))
        
        #k <- 12
        temp3 <- dplyr::filter(temp2, YR == t_yr[k])
        
        #C = confirmed breeder
        #P = probably breeder
        #O = observed
        #- = not observed
        
        #input is 'MM-DD' (05-01)
        #returns julian day
        dfun <- function(input)
        {
          t_date_lb <- paste0(t_yr[k], '-', input)
          t_j_lb <- format(as.Date(t_date_lb), '%j')
          return(t_j_lb)
        }
        
        PS1 <- dfun('05-01')
        PS2 <- dfun('05-11')
        PS3 <- dfun('05-21')
        PS4 <- dfun('05-31')
        PS5 <- dfun('06-10')
        PS6 <- dfun('06-20')
        PS7 <- dfun('06-30')
        PS8 <- dfun('07-10')
        PS9 <- dfun('07-20')
        PS10 <- dfun('07-30')
        PS11 <- dfun('08-09')
        PS12 <- dfun('08-18')
        
        periods <- c(PS1, PS2, PS3, PS4, PS5, PS6, PS7, PS8, PS9, PS10, PS11, PS12)
        
        #from first observaton observed breeding in hex cell:
        #input julian day if breeding
        #input 0 if observed in the season but not breeding
        #input NA if not observed at all
        
        cind <- grep('PS', colnames(temp3))
        
        #confirmed or possbile breeding for each station
        cp_ind <- c()
        ob_ind <- c()
        no_ind <- c()
        for (m in 1:NROW(temp3))
        {
          #m <- 1
          cp_ind <- c(cp_ind, which(temp3[m,cind] == 'C' | temp3[m,cind] == 'P'))
          ob_ind <- c(ob_ind, which(temp3[m,cind] == 'O'))
          no_ind <- c(no_ind, which(temp3[m,cind] == '-'))
        }
        
        if (length(cp_ind) > 0)
        {
          #lower bounds of period observed
          first_br_obs_lb <- as.numeric(periods[min(cp_ind)])
          #upper bounds of period observed
          first_br_obs_ub <- first_br_obs_lb + 9
          midpoint_br <- mean(c(first_br_obs_lb, first_br_obs_ub))
        } else {
          if (length(ob_ind) > 0)
          {
            first_br_obs_lb <- 0
            first_br_obs_ub <- 0
            midpoint_br <- 0
          } else {
            first_br_obs_lb <- NA
            first_br_obs_ub <- NA
            midpoint_br <- NA
          }
        }
        
        #column for midpoint breeding date
        #column for lower bounds breeding date
        #column for upper bounds breeding date
        #column for cell
        #column for year
        #column for species SCI NAME (with underscore)
        
        t_info <- dplyr::select(temp3, YR, cell, SCINAME)[1,]
        
        t_info$SCINAME <- gsub(' ', '_', t_info$SCINAME)
        
        t_info$midpoint <- midpoint_br
        t_info$l_bounds <- first_br_obs_lb
        t_info$u_bounds <- first_br_obs_ub
        #n_stations is number of stations that recorded probably or confirmed breeding
        t_info$n_stations <- length(cp_ind)
        
        MAPS_out[counter,] <- t_info
        counter <- counter + 1
      }
    }
  }
}



#remove NA padding (find first year row to have NA and subtract one from that index)
fin_ind <- min(which(is.na(MAPS_out$YR)))
MAPS_out2 <- MAPS_out[1:(fin_ind - 1),]


#NA = station operating in that year/cell, but bird not observed
#0 = station operating in that year/cell, bird observed, but not observed breeding
#JDAY = station operating in that year/cell, bird observed breeding

#write to rds
setwd(paste0(dir, 'Bird_Phenology/Data/Processed'))
saveRDS(MAPS_out2, paste0('breeding_MAPS_obs_', Sys.Date(), '.rds'))