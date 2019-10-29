# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Repo is cloned onto Xanadu at `/labs/Tingley/phenomismatch/Bird_Phenology`. Directories marked as 'ignored' are in the repo but are not tracked. Any data will need to be transfered manually. DO NOT transfer large files on Xanadu user account. 

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
  * `Processed/` (ignored) - data that have undergone processing
    * `breeding_cat_query_<DATE>` - ebird query for breeding - breeding/no breeding for each species input for `7-halfmax-br.R`
    * `eBird_query_<DATE>/` - ebird query for arrival - presence/absence data for each species - input for `2-halfmax-arr.R`
    * `daymet/` - processed data from Daymet - spring monthly averages for each cell (precip, tmin, tmax)
    * `halfmax_species_<DATE>/` - model output from `2-halfmax-arr.R`
    * `IAR_input_<DATE>/` - input files for IAR model (processed from `2-halfmax-arr.R` output)
    * `IAR_ouput_<DATE>/` - output files from IAR model (processed from `4-IAR-arr.R` output)
    * `arrival_master_<DATE>/` - combined output (across species) from `4-IAR-arr.R`
    * `breeding_master_<DATE>/` - combined output (across species) from breeding model
  * `Traits` (ignored) - bird life history trait data
  
* `Scripts/` - scripts to run analyses
  * `1-query-ebird.R` - query eBird database for arrival
  * `2-halfmax-arr/` - [run time: up to 4 weeks]
    * `2-halfmax-arr.R` - R script to fit GAM models to estimate bird arrival for each species/cell/year. Takes species argument (run from `2-<Genus_species>.sh` scripts)
    * `2-create-batch-scripts.sh` - script to create scripts (`2-<Genus_species>.sh`) for HPC job submission
    * `2-master-submit.sh` - script to be run (`./2-master-submit.sh`) on HPC cluster to submit all halfmax arr jobs
    * `species/` - contains job scripts for each species
      * `2-<Genus_species>.sh` - scripts to submit jobs for each species; run with `2-master-submit.sh`
    * `model_tests` - scripts to test which model to use to estimate halfmax
  * `3-process-arr-output.R` - process model output from `2-halfmax-arr.R`, in prep for `4-IAR-arr.R`
  * `4-IAR-arr/` - [run time: up to 3 weeks?]
      * `4-IAR-arr.R` - fit arrival dates using IAR model to derive arrival date estimates (one model per species)
      * `4-create-batch-scripts.sh` - script to create scripts (`4-<Genus_species>.sh`) for HPC job submission 
      * `4-master-submit.sh` - script to be run (`./4-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `4-<Genus_species>.sh` - scripts to submit jobs for each species; run with `4-master-submit.sh`
  * `4b-extract-daymet/` - [run time: up to 1 day?]
    * `4b-extract-daymet.R` - R script to extract Daymet data
    * `4b-run.sh` - script to run `4b-extract-temps.R` on cluster
  * `5-extract-arr-dates/`
    * `5-extract-arr-dates.R` - extract arrival dates for each species-cell-year from `4-IAR-arr.R` output
    * `5-run.sh` - run `5-extract-arr-dates.R` on cluster (typically run on desktop)
    * `README_arrival_master.txt` - README for arrival master data
  * `6-halfmax-juvs/` - halfmax model for fledging dates (derived from MAPS data) [run time: ???]
    * `6-query-nesting-data.R` - query eBird breeding codes and MAPS data
    * `6-halfmax-juvs.R` - R script to fit GAM models to estimate bird fledging date (from MAPS data) for each species/cell/year. Takes species argument (run from `6-<Genus_species>.sh` scripts)
    * `6-create-batch-scripts.sh` - script to create scripts (`6-<Genus_species>.sh`) for HPC job submission
    * `6-master-submit.sh` - script to be run (`./6-master-submit.sh`) on HPC cluster to submit all halfmax arr jobs
    * `species/` - contains job scripts for each species
      * `6-<Genus_species>.sh` - scripts to submit jobs for each species; run with `6-master-submit.sh`
  * `7-halfmax-br/` - halfmax model for breeding dates (derived from eBird breeding codes) [run time: ???]
      * `7-halfmax-br.R` - R script to fit GAM models to estimate bird breeding date (from eBird breeding codes) for each species/cell/year. Takes species argument (run from `7-<Genus_species>.sh` scripts)
      * `7-create-batch-scripts.sh` - script to create scripts (`7-<Genus_species>.sh`) for HPC job submission
      * `7-master-submit.sh` - script to be run (`./7-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `7-<Genus_species>.sh` - scripts to submit jobs for each species; run with `7-master-submit.sh`
  * `8-process-juv-br-output/` - process br and juv halfmax results
    * `8-process-juv-output.R` - process model output from `6-halfmax-juvs.R`
    * `8-process-br-output.R` - process model output from `7-halfmax-br.R`
  * `Archive/` - Archived scripts
  * `Other/`
  * `IAR-sandbox/` - Jacob working scripts

* `Figures/` (ignored)

* `Resources/` - various resources
