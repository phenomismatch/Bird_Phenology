# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Repo is cloned onto Xanadu at `/UCHC/LABS/Tingley/phenomismatch/Bird_Phenology`. Directories marked as 'ignored' are in the repo but are not tracked. Any data will need to be transfered manually. DO NOT transfer large files on Xanadu user account. 

For transfer TO Xanadu use:

`scp -r ~/SOURCE/PATH USER_NAME@transfer.cam.uchc.edu:DESTINATION/PATH`

For transfer FROM Xanadu use:

`scp -r USER_NAME@transfer.cam.uchc.edu:SOURCE/PATH ~/DESTINATION/PATH`

Repository structure:

* `Data/` - Datasets relevant for project (consider using [Piggyback](https://cran.r-project.org/web/packages/piggyback/vignettes/intro.html))
  * `eBird_species_list.txt` - list of species to be used in analyses
  * `IAR_species_list.txt` - list of species to be modeled using IAR
  * `db_pass.txt` (ignored) - database password to pass to DB when querying
  * `BirdLife_range_maps` (ignored) - range maps for birds from BirdLife International
  * `Hurlbert_breeding_data` - breeding phenology data taken from literature from repo [here](https://github.com/hurlbertlab/bird-repro-times)
  * `MAPS_obs` (ignored) - observational data on breeding from MAPS program
  * `Nestwatch` (ignored) - Nestwatch program data
  * `Processed/` (ignored) - data that have undergone processing
    * `db_query_<DATE>/` (ignored) - ebird checklist presence/absence data for each species queried from database - input for logit cubic models (see `Scripts/eBird_Nestwatch/2-logit-cubic`)
    * `halfmax_species_<DATE>/` (ignored) - output files for each specices from 2-logit-cubic - estimated half-max params and model diagnostics
    * `IAR_input_<DATE>/` (ignored) - input files for IAR model (processed from halfmax esimation)
    * `IAR_ouput_<DATE>/` (ignored) - output files from IAR - estimated bird arrival dates over breedin and migratory range (spatially smoothed)
    * `master_arrival_<DATE>.rds` (ignored) - combined output from IAR model

* `Scripts/` - scripts to run analyses
  * `Climate_Veg/` - comparing vegetation phenology products
  * `eBird_Nestwatch/` - eBird and nestwatch phenology scripts
    * `1-query-ebird.R` - query eBird data from database
    * `2-logit-cubic/` - scripts to fit logit-cubic to estimate bird arrival for each species-cell-year - each species run as separate job on HPC cluster
      * `2-create-batch-scripts.sh` - script to create scripts (`2-<Genus_species>.sh`) for HPC job submission
      * `2-logit-cubic.R` - R script that takes species argument (run from `2-<Genus_species>.sh` scripts)
      * `2-master-submit.sh` - script to be run (`./2-master-submit.sh`) on HPC cluster to submit all logit cubic jobs
      * `species/` - contains job scripts for each species
        * `2-<Genus_species>.sh` - scripts to submit jobs for each species (153 scripts, one for each species); run with `2-master-submit.sh`
    * `3-process-output.R` - process model output from `2-logit-cubic.R`, in prep for `4-ICAR-model.R`
    * `4-IAR/` - scripts to fit IAR model to estimate bird arrival for each species-cell-year with spatial smoothing - one model per species
       * `4-create-batch-scripts.sh` - script to create scripts (`4-<Genus_species>.sh`) for HPC job submission    
       * `4-ICAR-model.R` - fit arrival dates using ICAR model to derive arrival date estimates
       * `4-master-submit.sh` - script to be run (`./2-master-submit.sh`) on HPC cluster to submit all logit cubic jobs
       * `species/` - contains job scripts for each species
         * `4-<Genus_species>.sh` - scripts to submit jobs for each species (153 scripts, one for each species); run with `4-master-submit.sh`
    * `5-extract-arr-dates.R` - extract arrival dates for each species-cell-year
    * `6-nesting0arrival-model.R` - model nesting date (Nestwatch) as a function of arrival date (IAR output)
    * `7-query-nesting-data.R` - query nesting data from other data sources (MAPS, eBird breeding codes, etc.)
    * `8-phen-compare.R` - correlations among veg and bird phenology products (NEEDS UPDATING)
    * `9-breed-arr-model.R` - breeding date (nest watch) as a function of arrival date (eBird) (NEEDS UPDATING)

* `Figures/` - various figs
    * `pre_post_IAR_maps/` - maps to visualize bird arrival dates pre- and post-IAR
    
* `Resources/` (ignored) - Notes/papers

* `Results/` (ignored)
  * `Plots/` (ignored) - Plots of logit cubic model fit with presence absence data - one plot for every species, year, cell
