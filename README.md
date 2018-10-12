# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Repo is cloned onto Xanadu at `/UCHC/LABS/Tingley/phenomismatch/Bird_Phenology`. Directories marked as 'ignored' are in the repo but are not tracked. Any data will need to be transfered manually. DO NOT transfer large files on Xanadu user account. Use instead:

`scp -r ~/SOURCE/PATH USER_NAME@transfer.cam.uchc.edu:DESTINATION/PATH`


Repository structure:

* `Data/` - Datasets relevant for project
  * `eBird_species_list.txt` - List of species to be used in analyes
  * `db_pass.txt` (ignored) - database password to pass to DB when querying
  * `Processed/` (ignored) - Data that have undergone processing
    * `ebird_NA_phen_proc_species` (ignored) - ebird checklist presence/absence for each species - input for logit cubic models (see `Scripts/eBird_Nestwatch/2-`)
    * `halfmax_species` (ignored) - output files from logit cubic for each species
  * `Raw/` (ignored) - Raw data that have not undergone processing
    * `eBird` - eBird Reference Sataset

* `Scripts/` - Scripts to run analyses
  * `Climate_Veg/` - Comparing vegetation phenology products
  * `eBird_Nestwatch/` - eBird and nestwatch phenology scripts
    * `1-process-ebird-data.R` - process eBird data
    * `1b-query-ebird.R` - query eBird data from database (replaces `1-import-ebird-data.R`)
    * `2-logit-cubic/` - scripts to fit logit-cubic to get bird arrival for each species-cell-year - each species run as separate job on HPC cluster
      * `2-master-submit.sh` - script to be run on HPC cluster to submit all logit cubic jobs
      * `2-create-batch-scripts.sh` - script to create scripts (`<Genus_species>.sh`) for HPC job submission
      * `2-<Genus_species>.sh` - scripts to submit jobs for each species (153 scripts, one for each species); run with `2-master-submit.sh`
      * `2-logit-cubic.R` - R script that takes species argument (sourced by HPC scripts)
    * `3-ICAR-model.R` - fit arrival dates using ICAR model to derive arrival date estimates
    * `3b-ICAR-model-ns.R` - fit arrival dates using ICAR model with both spatial and non-spatial components (as opposed to just spatial component as with `3-ICAR-model.R`)
    * `3c-ICAR-model-ns-parallel.R` - same as `3b-ICAR-model-ns.R`, except run in parallel
    * `4-extract-arr-dates.R` - extract arrival dates for each species-cell-year
    * `5-phen-compare.R` - correlations among veg and bird phenology products
    * `6-breed-arr-model.R` - breeding date (nest watch) as a function of arrival date (eBird)

* `Notes/` (ignored) - Scratch notes

* `Results/` (ignored) - Model output
  * `Plots/` (ignored) - Plots of halfmax fits
