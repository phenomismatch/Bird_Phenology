# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Repo is cloned onto Xanadu at `/UCHC/LABS/Tingley/phenomismatch/Bird_Phenology`. Directories marked as 'ignored' are in the repo but are not tracked. Any data will need to be transfered manually. DO NOT transfer large files on Xanadu user account. 

For transfer TO Xanadu use:

`scp -r ~/SOURCE/PATH USER_NAME@transfer.cam.uchc.edu:DESTINATION/PATH`

For transfer FROM Xanadu use:

`scp -r USER_NAME@transfer.cam.uchc.edu:SOURCE/PATH ~/DESTINATION/PATH`

To check for divergences (or other cluster issues) in STDERR:

`grep WORD *.err -lR`

Query database with:

`psql "sslmode=disable dbname=sightings user=cyoungflesh hostaddr=35.221.16.125"`

Repository structure:

* `Data/` - Datasets relevant for project
  * `eBird_species_list.txt` - initial list of species to be queried
  * `IAR_species_list.txt` - list of species to be modeled using IAR
  * `db_pass.txt` (ignored) - database password
  * `BirdLife_range_maps` (ignored) - range maps for birds from BirdLife International
  * `Hurlbert_breeding_data` - breeding phenology data taken from repo [here](https://github.com/hurlbertlab/bird-repro-times)
  * `MAPS_obs` (ignored) - MAPS program data
  * `Nestwatch` (ignored) - Nestwatch program data
  * `Processed/` (ignored) - data that have undergone processing
    * `eBird_query_<DATE>/` - ebird query for arrival - presence/absence data for each species - input for `2-logit-cubic.R`
    * `breeding_cat_query_<DATE>` - ebird query for breeding - breeding/no breeding for each species input for `7-logit-cubic-breeding.R`
    * `daymet/` - data from Daymet - max temperature
    * `halfmax_species_<DATE>/` - model output from `2-logit-cubic.R`
    * `IAR_input_<DATE>/` - input files for IAR model (processed from `2-logit-cubic.R` output)
    * `IAR_ouput_<DATE>/` - output files from IAR model (processed from `4-IAR-model.R` output)
    * `arrival_master_<DATE>.rds` - combined output (across species) from arrival IAR model
    * `breeding_master_<DATE>.rds` - combined output (across species) from breeding model
  * `Traits` (ignored) - bird life history trait data
  
* `Scripts/` - scripts to run analyses
  * `Climate_Veg/` - comparing vegetation phenology products
  * `eBird_Nestwatch/` - eBird and nestwatch phenology scripts
    * `1-query-ebird.R` - query eBird database for arrival
    * `2-halfmax-arr/` - scripts to fit model to estimate bird arrival for each species-cell-year - each species run as separate job on HPC cluster [run time: up to 7 days]
      * `2-create-batch-scripts.sh` - script to create scripts (`2-<Genus_species>.sh`) for HPC job submission
      * `2-halfmax-arr.R` - R script that takes species argument (run from `2-<Genus_species>.sh` scripts)
      * `2-master-submit.sh` - script to be run (`./2-master-submit.sh`) on HPC cluster to submit all halfmax arr jobs
      * `species/` - contains job scripts for each species
        * `2-<Genus_species>.sh` - scripts to submit jobs for each species; run with `2-master-submit.sh`
      * `model_tests` - scripts to test which model to use to estimate halfmax
    * `3-process-arr-output.R` - process model output from `2-halfmax-arr.R`, in prep for `4-IAR-model.R`
    * `3b-extract-temps.R` - extract temperature data from Daymet in prep for `4-IAR-model.R`
    * `4-IAR-arr/` - scripts to fit IAR model to estimate bird arrival for each species-cell-year with spatial smoothing - one model per species [run time: ???]
       * `4-create-batch-scripts.sh` - script to create scripts (`4-<Genus_species>.sh`) for HPC job submission    
       * `4-IAR-arr.R` - fit arrival dates using IAR model to derive arrival date estimates
       * `4-master-submit.sh` - script to be run (`./4-master-submit.sh`) on HPC cluster to submit all jobs
       * `species/` - contains job scripts for each species
         * `4-<Genus_species>.sh` - scripts to submit jobs for each species; run with `4-master-submit.sh`
    * `5-extract-arr-dates.R` - extract arrival dates for each species-cell-year from `4-IAR-model.R` output
    * `6-query-nesting-data.R` - query eBird database for breeding; also MAPS and Nestwatch breeding processing
    * `7-halfmax-br/` - scripts to fit model to estimate bird breeding for each species-cell-year - each species run as separate job on HPC cluster [run time: ???]
      * `7-create-batch-scripts.sh` - script to create scripts (`7-<Genus_species>.sh`) for HPC job submission
      * `7-halfmax-br.R` - R script that takes species argument (run from `7-<Genus_species>.sh` scripts)
      * `7-master-submit.sh` - script to be run (`./7-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `7-<Genus_species>.sh` - scripts to submit jobs for each species; run with `7-master-submit.sh`
    * `8-process-br-output.R` - process model output from `7-halfmax-br.R`
    * `9-IAR-br` - scripts to fit IAR model to estimate bird breeding for each species-cell-year with spatial smoothing - one model per species [run time: ???]
       * `9-create-batch-scripts.sh` - script to create scripts (`9-<Genus_species>.sh`) for HPC job submission    
       * `9-IAR-br.R` - fit breeding dates using IAR model to derive arrival date estimates
       * `9-master-submit.sh` - script to be run (`./9-master-submit.sh`) on HPC cluster to submit all jobs
       * `species/` - contains job scripts for each species
         * `9-<Genus_species>.sh` - scripts to submit jobs for each species; run with `9-master-submit.sh`
    * `10-trend-models/` - scripts to fit models looking at trends over time
      * `10-arr-time-lat-IND-corr.R` - arrival trends over time - correlation between intercepts and slopes
      * `10-br-time-lat-IND-corr.R` - arrival trends over time
      * `10-breeding-arrival-model.R` - model relationship between arrival timing and breeding timing

* `Figures/` (ignored)
  * `arrival_trends/` - trends in arrival over time from `10-arr-time-lat-IND-corr.R`
  * `breeding_trends/` - trends in breeding over time from `10-br-time-lat-IND-corr.R`
  * `halfmax/` - halfmax figs
    * `arrival/` - figs from `2-halfmax-arr.R`
    * `breeding/` - figs from `7-halfmax-br.R`
  * `pre_post_IAR_maps/` - maps to visualize bird arrival dates pre- and post-IAR 
    * `arrival/` - figs from `4-IAR-arr.R`
    * `breeding/` - figs from `9-IAR-br.R`
    
* `Resources/` - various resources
