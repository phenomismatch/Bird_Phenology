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
  * `db_pass.txt` (ignored) - database password to pass to DB when querying
  * `Processed/` (ignored) - Data that have undergone processing
    * `db_query_<DATE>` (ignored) - ebird checklist presence/absence data for each species queried from database - input for logit cubic models (see `Scripts/eBird_Nestwatch/2-logit-cubic`)
    * `halfmax_species` (ignored) - output files for each specices from 2-logit-cubic - estimated half-max params and model diagnostics
  * `Raw/` (ignored) - Raw data that have not undergone processing
    * `eBird` - eBird Reference Dataset (no longer used as data is being queried from database)

* `Scripts/` - Scripts to run analyses
  * `Climate_Veg/` - Comparing vegetation phenology products
  * `eBird_Nestwatch/` - eBird and nestwatch phenology scripts
    * `1-process-ebird-data.R` - process eBird data
    * `1b-query-ebird.R` - query eBird data from database (replaces `1-import-ebird-data.R`)
    * `2-logit-cubic/` - scripts to fit logit-cubic to get bird arrival for each species-cell-year - each species run as separate job on HPC cluster
      * `2-master-submit.sh` - script to be run (`./2-master-submit.sh`) on HPC cluster to submit all logit cubic jobs
      * `2-create-batch-scripts.sh` - script to create scripts (`2-<Genus_species>.sh`) for HPC job submission
      * `2-logit-cubic.R` - R script that takes species argument (run from `2-<Genus_species>.sh` scripts)
      * `species/` - contains job scripts for each species
        * `2-<Genus_species>.sh` - scripts to submit jobs for each species (153 scripts, one for each species); run with `2-master-submit.sh`
    * `3-ICAR-model.R` - fit arrival dates using ICAR model to derive arrival date estimates
    * `3b-ICAR-model-ns.R` - fit arrival dates using ICAR model with both spatial and non-spatial components (as opposed to just spatial component as with `3-ICAR-model.R`)
    * `3c-ICAR-model-ns-parallel.R` - same as `3b-ICAR-model-ns.R`, except run in parallel
    * `4-extract-arr-dates.R` - extract arrival dates for each species-cell-year
    * `5-phen-compare.R` - correlations among veg and bird phenology products
    * `6-breed-arr-model.R` - breeding date (nest watch) as a function of arrival date (eBird)

* `Notes/` (ignored) - Scratch notes

* `Results/` (ignored) - Model output
  * `Plots/` (ignored) - Plots of logit cubic model fit with presence absence data - one plot for every species, year, cell
