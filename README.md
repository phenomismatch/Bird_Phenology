# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Repo is cloned onto Xanadu at `/UCHC/LABS/Tingley/phenomismatch/Bird_Phenology`. Directories marked as 'ignored' are in the repo but are not tracked. Any data will need to be transfered manually. DO NOT transfer large files on Xanadu user account. 

For transfer TO Xanadu use:

`scp -r ~/SOURCE/PATH USER_NAME@transfer.cam.uchc.edu:DESTINATION/PATH`

For transfer FROM Xanadu use:

`scp -r USER_NAME@transfer.cam.uchc.edu:SOURCE/PATH ~/DESTINATION/PATH`

Repository structure:

* `Data/` - Datasets relevant for project (consider using [Piggback](https://cran.r-project.org/web/packages/piggyback/vignettes/intro.html))
  * `eBird_species_list.txt` - List of species to be used in analyes
  * `IAR_species_list.txt` - List of species to be modeled using IAR
  * `db_pass.txt` (ignored) - database password to pass to DB when querying
  * `Processed/` (ignored) - Data that have undergone processing
    * `db_query_<DATE>` (ignored) - ebird checklist presence/absence data for each species queried from database - input for logit cubic models (see `Scripts/eBird_Nestwatch/2-logit-cubic`)
    * `halfmax_species_<DATE>` (ignored) - output files for each specices from 2-logit-cubic - estimated half-max params and model diagnostics
    * `IAR_input_<DATE>` (ignored) - input files for IAR model (processed from halfmax esimation)
    * `IAR_ouput_<DATE>` (ignored) - output files from IAR - estimated bird arrival dates over breedin and migratory range (spatially smoothed)
  * `Raw/` (ignored) - Raw data that have not undergone processing
    * `eBird` - eBird Reference Dataset (no longer used as data is being queried from database)

* `Scripts/` - Scripts to run analyses
  * `Climate_Veg/` - Comparing vegetation phenology products
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
    * `5-extract-arr-dates.R` - extract arrival dates for each species-cell-year (NEEDS UPDATING)
    * `6-phen-compare.R` - correlations among veg and bird phenology products (NEEDS UPDATING)
    * `7-breed-arr-model.R` - breeding date (nest watch) as a function of arrival date (eBird) (NEEDS UPDATING)

* `Figures/` - various figs
    * `pre_post_IAR_maps/` - maps to visualize bird arrival dates pre- and post-IAR
    
* `Resources/` (ignored) - Notes/papers

* `Results/` (ignored)
  * `Plots/` (ignored) - Plots of logit cubic model fit with presence absence data - one plot for every species, year, cell
