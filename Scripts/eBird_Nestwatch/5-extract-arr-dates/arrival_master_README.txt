==============================
arrival_master_DATE.rds README
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

cell_lat - latitude for center of hex cell

cell_lon - longitude for center of hex cell

year - year

mean_pre_IAR - posterior mean for arrival date for that species/cell/year (derived from GAM model - 2-halfmax-arr.R). This is used as an input for the IAR model (4-IAR-arr.R)

sd_pre_IAR - posterior standard deviation for arrival date for that species/cell/year (derived from GAM model - 2-halfmax-arr.R). This is used as an input for the IAR model (4-IAR-arr.R)

mean_post_IAR - posterior mean for arrival date for that species/cell/year (derived from IAR model - 4-IAR-arr.R)

sd_post_IAR - posterior standard deviation for arrival date for that species/cell/year (derived from IAR model - 4-IAR-arr.R)

mean_gamma - posterior mean for gamma parameter (derived from IAR model - 4-IAR-arr.R). Gamma is a cell random effect and represents the mean arrival date at a given cell. Each hex cell for a given species will have an associated gamma parameter.

sd_gamma - posterior standard deviation for gamma parameter (derived from IAR model - 4-IAR-arr.R). Gamma is a cell random effect and represents the mean arrival date at a given cell. Each hex cell for a given species will have an associated gamma parameter.

mean_beta0 - posterior mean for beta0 parameter (derived from IAR model - 4-IAR-arr.R). Beta0 is a year random effect and represents the arrival date intercept for a given year. Each year for a given species will have an associated beta0 parameter.

sd_beta0 - posterior standard deviation for beta0 parameter (derived from IAR model - 4-IAR-arr.R). Beta0 is a year random effect and represents the arrival date intercept for a given year. Each year for a given species will have an associated beta0 parameter.

mean_alpha_gamma - posterior mean for alpha_gamma (derived from IAR model - 4-IAR-arr.R). Alpha_gamma is the intercept for the modeled gamma parameter and represents the expected cell random effect (estimated mean arrival date at a given cell) at 0 degrees latitude.

sd_alpha_gamma - posterior standard deviation for alpha_gamma (derived from IAR model - 4-IAR-arr.R). Alpha_gamma is the intercept for the modeled gamma parameter and represents the expected cell random effect (estimated mean arrival date at a given cell) at 0 degrees latitude.

mean_beta_gamma - posterior mean for beta_gamma (derived from IAR model - 4-IAR-arr.R). Beta_gamma is the slope (as a function of latitude) for the modeled gamma parameter and represents the degree to which mean arrival changes over latitude (days per degree latitude). This can also be thought of as mean northward migration speed across the modeled range of the species. 

sd_beta_gamma - posterior standard deviation for beta_gamma (derived from IAR model - 4-IAR-arr.R). Beta_gamma is the slope (as a function of latitude) for the modeled gamma parameter and represents the degree to which mean arrival changes over latitude (days per degree latitude). This can also be thought of as mean northward migration speed across the modeled range of the species. 

num_diverge - number of diverges for the IAR model for this species. Divergences occur when the NUTS sampler encounters areas of high curvature in the joint probability density. Divergences should be eliminated by reparameterizing or increasing the Adapt_delta argument for Stan.

max_rhat - maximum Rhat value for any parameter from the IAR model.

min_neff - minimum number of effective samples for any parameter from the IAR model.

dm_tmax_F - Feb mean daily max temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_tmax_M - March mean daily max temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_tmax_A - April mean daily max temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_tmax_FMA - Feb-April mean daily max temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_tmin_F - Feb mean daily min temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_tmin_M - March mean daily min temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391) 

dm_tmin_A - April mean daily min temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_tmin_FMA - Feb-April mean daily min temperature for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_precip_F - Feb mean daily precipitation for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_precip_M - March mean daily precipitation for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_precip_A - April mean daily precipitation for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)

dm_precip_FMA - Feb-April mean daily precipitation for a given hex cell in a given year - derived from daymet (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1391)
