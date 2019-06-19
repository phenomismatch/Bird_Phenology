#############################
# Analysis of pheno trends output
#
#############################

# PARAMETER DESCRIPTIONS
# ======================
# alpha_gamma - arrival date of that species at 0 degrees lat
# beta_gamma - speed of migration (days / degrees lat)

# mu_alpha - absolute arrival date for species (intercepts drawn from this mean)
# sigma - interannual variability in arrival date after accounting for trend
# alpha_beta - magnitude of phenological change over time
# beta_beta - effect of lat on magnitude of phenological change over time



# read in data ------------------------------------------------------------

setwd('~/Desktop/Bird_Phenology_Offline/Data/Processed/trends_summary_2019-06-17/')

out <- readRDS('pheno_trends_master_2019-06-17.rds')
sp <- out$species



# size doesn't matter for halfmax -----------------------------------------

# tt <- rnorm(100, out$mn_sigma[1], out$sd_sigma[1])
# d <- density(tt)
# plot(d)
# xmax <- d$x[d$y == max(d$y)]
# 
# x1 <- d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax] - max(d$y)/2))]
# x2 <- d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax] - max(d$y)/2))]
# points(c(x1, x2), c(d$y[d$x==x1], d$y[d$x==x2]), col="red")




# alpha_beta (and beta_beta) ~ sigma - lm in loop --------------------------------------

oslope <- c()
for (i in 1:1000)
{
  #i <- 1
  
  dt <- c()
  for (j in 1:length(sp))
  {
    #j <- 1
    tmn_sigma <- out[j,'mn_sigma']
    tsd_sigma <- out[j,'sd_sigma']
  
    tvals_sigma <- rnorm(1, tmn_sigma, tsd_sigma)
    while (tvals_sigma <= 0)
    {
      tvals_sigma <- rnorm(1, tmn_sigma, tsd_sigma)
    }
  
    tmn_alpha_beta <- out[j,'mn_alpha_beta']
    tsd_alpha_beta <- out[j,'sd_alpha_beta']
    tvals_alpha_beta <- rnorm(1, tmn_alpha_beta, tsd_alpha_beta)
    
    # tmn_beta_beta <- out[j,'mn_beta_beta']
    # tsd_beta_beta <- out[j,'sd_beta_beta']
    # tvals_beta_beta <- rnorm(1, tmn_beta_beta, tsd_beta_beta)
    
    temp <- c(tvals_sigma, tvals_alpha_beta)
    # temp <- c(tvals_sigma, tvals_beta_beta)
    dt <- rbind(dt, temp)
  }
  colnames(dt) <- c('sigma', 'alpha_beta')
  #colnames(dt) <- c('sigma', 'beta_beta')
  
  tfit <- lm(abs(dt[,'alpha_beta']) ~ dt[,'sigma'])
  #tfit <- lm(abs(dt[,'beta_beta']) ~ dt[,'sigma'])
  
  #slope
  tslope <- coefficients(tfit)[2]
  
  #output
  oslope <- c(oslope, tslope)
}
hist(oslope)





# raw relationship between estimates ------------------------------------------

fit_fun <- function(first, second)
{
  tf <- lm(second ~ first)
  plot(first, second)
  abline(tf, col = 'red')
  #hist(residuals(tf))
  print(summary(tf))
}

fit_fun(out[, 'mn_sigma'], abs(out[, 'mn_mu_alpha']))
fit_fun(out[, 'mn_sigma'], abs(out[, 'mn_alpha_beta']))
fit_fun(out[, 'mn_sigma'], abs(out[, 'mn_beta_beta']))

#more interannual variability = more change over time


# plot parameter estimates ------------------------------------------------

MCMCvis::MCMCplot(mu_alpha_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(50, 175), main = 'mu_alpha')
MCMCvis::MCMCplot(sigma_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(0, 10), main = 'sigma')
MCMCvis::MCMCplot(alpha_beta_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(-1, 1.5), main = 'alpha_beta')
MCMCvis::MCMCplot(beta_beta_post, guide_lines = TRUE, params = sp,
                  sz_labels = 0.6, xlim = c(-0.25, 0.5), main = 'beta_beta')



# traits ------------------------------------------------------------------

setwd("~/Google_Drive/R/Bird_Phenology/Data/Traits")

traits <- read.csv('Trait_database-2019-06-17.csv', stringsAsFactors = FALSE)
#add underscore
traits$species <- gsub(' ', '_', traits$SCI_NAME)

out_j <- dplyr::left_join(out, traits, by = 'species')

out_jf <- dplyr::select(out_j, species, mn_mu_alpha, sd_mu_alpha, mn_sigma, 
                        sd_sigma, mn_alpha_beta, sd_alpha_beta, mn_beta_beta,
                        sd_beta_beta, IUCN_STATUS, BODY_MASS_ELTON, MIGRATION_DISTANCE_LASORTE, 
                        SPRING_MIGRATION_SPEED_LASORTE, CLUTCH_SIZE_BONA, BROODS_PER_YEAR_BONA, 
                        EGG_MASS_BONA, MAX_LIFESPAN_BONA, INCUBATION_PERIOD_LENGTH_BONA, 
                        FLEDGING_AGE_BONA, DIET_INV_ELTON, DIET_VUNK_ELTON, DIET_SCAV_ELTON, 
                        DIET_FRUIT_ELTON, DIET_SEED_ELTON, rng_lat, n_cells, n_years, num_diverge, max_rhat, min_neff)



idx <- which(colnames(out_jf) %in% c('species', 'IUCN_STATUS', 'sd_mu_alpha', 'sd_sigma', 
                                     'sd_alpha_beta', 'sd_beta_beta', 
                                     'num_diverge', 'max_rhat', 'min_neff'))

pairs(out_jf[,-idx])


out_jf2 <- dplyr::filter(out_jf, n_cells > 10, n_years > 10, rng_lat > 10)


fit_fun(out_jf$MIGRATION_DISTANCE_LASORTE, out_jf$mn_mu_alpha)
fit_fun(out_jf$SPRING_MIGRATION_SPEED_LASORTE, out_jf$mn_mu_alpha)

fit_fun(out_jf$MIGRATION_DISTANCE_LASORTE, out_jf$mn_sigma)
fit_fun(out_jf$SPRING_MIGRATION_SPEED_LASORTE, out_jf$mn_sigma)

fit_fun(out_jf2$MIGRATION_DISTANCE_LASORTE, out_jf2$mn_alpha_beta)
fit_fun(out_jf2$SPRING_MIGRATION_SPEED_LASORTE, out_jf2$mn_alpha_beta)

fit_fun(out_jf2$MIGRATION_DISTANCE_LASORTE, out_jf2$mn_beta_beta)
fit_fun(out_jf2$SPRING_MIGRATION_SPEED_LASORTE, out_jf2$mn_beta_beta)

plot(out_jf2$rng_lat, out_jf2$mn_alpha_beta)
plot(out_jf2$rng_lat, out_jf2$mn_beta_beta)

# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_mu_alpha)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_mu_alpha)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_mu_alpha)
# 
# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_sigma)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_sigma)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_sigma)
# 
# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_alpha_beta)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_alpha_beta)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_alpha_beta)
# 
# fit_fun(out_jf$CLUTCH_SIZE_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$BROODS_PER_YEAR_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$EGG_MASS_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$BODY_MASS_ELTON, out_jf$mn_beta_beta)
# fit_fun(out_jf$MAX_LIFESPAN_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$INCUBATION_PERIOD_LENGTH_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$mn_beta_beta)
# fit_fun(out_jf$DIET_INV_ELTON, out_jf$mn_beta_beta)

fit_fun(out_jf$SPRING_MIGRATION_SPEED_LASORTE, out_jf$MIGRATION_DISTANCE_LASORTE)
fit_fun(out_jf$FLEDGING_AGE_BONA, out_jf$CLUTCH_SIZE_BONA)


