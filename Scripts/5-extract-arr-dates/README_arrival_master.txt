==============================
README arrival_master_DATE.rds
==============================

Summary
-------
.rds object containing output from IAR arrival models (one for each species). Model output is produced from 4-IAR-arr.R.

.rds object can be read into R using the following:
obj <- readRDS('arrival_master_DATE.rds')



Dataframe columns
-----------------

species - scientific name for species

cell - hex cell number (hex cells generated with dggrid package in R)

mig_cell - logical specifying whether any portion of that partical hex cell is classified as in that species' migratory range (according to the BirdLife International range maps)

breed_cell - logical specifying whether any portion of that partical hex cell is classified as in that species' breeding range (according to the BirdLife International range maps)

cell_lat - latitude for center of hex cell

cell_lon - longitude for center of hex cell

year - year

arr_GAM_mean - posterior mean for arrival date for that species/cell/year (derived from GAM model - 2-halfmax-arr.R). This is used as an input for the IAR model (4-IAR-arr.R)

arr_GAM_sd - posterior standard deviation for arrival date for that species/cell/year (derived from GAM model - 2-halfmax-arr.R). This is used as an input for the IAR model (4-IAR-arr.R)

arr_IAR_mean - posterior mean for arrival date for that species/cell/year (derived from IAR model - 4-IAR-arr.R)

arr_IAR_sd - posterior standard deviation for arrival date for that species/cell/year (derived from IAR model - 4-IAR-arr.R)

gamma_mean - posterior mean for gamma parameter (derived from IAR model - 4-IAR-arr.R). Gamma is a cell random effect and represents the mean arrival date at a given cell. Each hex cell for a given species will have an associated gamma parameter.

gamma_sd - posterior standard deviation for gamma parameter (derived from IAR model - 4-IAR-arr.R). Gamma is a cell random effect and represents the mean arrival date at a given cell. Each hex cell for a given species will have an associated gamma parameter.

beta0_mean - posterior mean for beta0 parameter (derived from IAR model - 4-IAR-arr.R). Beta0 is a year random effect and represents the arrival date intercept for a given year. Each year for a given species will have an associated beta0 parameter.

beta0_sd - posterior standard deviation for beta0 parameter (derived from IAR model - 4-IAR-arr.R). Beta0 is a year random effect and represents the arrival date intercept for a given year. Each year for a given species will have an associated beta0 parameter.

alpha_gamma_mean - posterior mean for alpha_gamma (derived from IAR model - 4-IAR-arr.R). Alpha_gamma is the intercept for the modeled gamma parameter and represents the expected cell random effect (estimated mean arrival date at a given cell) at 0 degrees latitude.

alpha_gamma_sd - posterior standard deviation for alpha_gamma (derived from IAR model - 4-IAR-arr.R). Alpha_gamma is the intercept for the modeled gamma parameter and represents the expected cell random effect (estimated mean arrival date at a given cell) at 0 degrees latitude.

beta_gamma_mean - posterior mean for beta_gamma (derived from IAR model - 4-IAR-arr.R). Beta_gamma is the slope (as a function of latitude) for the modeled gamma parameter and represents the degree to which mean arrival changes over latitude (days per degree latitude). This can also be thought of as mean northward migration speed across the modeled range of the species. 

beta_gamma_sd - posterior standard deviation for beta_gamma (derived from IAR model - 4-IAR-arr.R). Beta_gamma is the slope (as a function of latitude) for the modeled gamma parameter and represents the degree to which mean arrival changes over latitude (days per degree latitude). This can also be thought of as mean northward migration speed across the modeled range of the species. 

num_diverge - number of diverges for the IAR model for this species. Divergences occur when the NUTS sampler encounters areas of high curvature in the joint probability density. Divergences should be eliminated by reparameterizing or increasing the Adapt_delta argument for Stan.

max_rhat - maximum Rhat value for any parameter from the IAR model.

min_neff - minimum number of effective samples for any parameter from the IAR model.

