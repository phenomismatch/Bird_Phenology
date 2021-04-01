# Bird_Phenology

[![DOI](https://zenodo.org/badge/136236843.svg)](https://zenodo.org/badge/latestdoi/136236843)

Code for characterizing bird arrival and breeding phenology.

**Associated publications:**

Youngflesh, C., Socolar, J., Amaral, B.R., Arab, A., Guralnick, R.P., Hurlbert, A.H., LaFrance, R., Mayor, S.J., Miller, D.A.W., Tingley, M.W. 2021. Migratory strategy drives species-level variation in bird sensitivity to vegetation green-up. __*Nature Ecology and Evolution*__


**Repository structure:**

* `Data/` - Datasets relevant for project
  * `eBird_species_list.txt` - initial list of species to be queried
  * `IAR_species_list.txt` - list of species to be modeled using IAR
  * `arr_species_list.txt` - list of species to quantify breeding phenology
  * `species_reference.csv` - table of taxonomic information for species of interest
  * `db_pass.txt` (ignored) - database password
  * `hex_grid_crop` - .shp files for cells over study area
  * `BirdLife_range_maps/` (ignored) - range maps for birds from BirdLife International
  * `Processed/` (ignored) - data that have undergone processing
    * `arrival_GAM_<DATE>` - model output from `2-arr-GAM.R`
    * `arrival_master_<DATE>` - output from `5-extract-arr-dates.R`
    * `breeding_GAM_<DATE>` - model output from `8-br-GAM.R`
    * `breeding_master_<DATE>` - output from `9-process-br-output.R`
    * `eBird_arrival_query_<DATE>` - output from `1-query-ebird.R`
    * `eBird_breeding_query_<DATE>` - output from `6-query-nesting-data.R`
    * `IAR_input_<DATE>` - output from `3-process-arr-output.R`
    * `IAR_output_<DATE>` - output from `4-arr-IAR_hm.R`
    * `juv_GAM_<DATE>` - model output from `7-juv-GAM.R`
    * `juv_master_<DATE>` - output from `9-process-juv-output.R`
  
* `Scripts/` - scripts to run analyses
  * `1-query-ebird.R` - query eBird database for arrival
  * `2-arr-GAM/`
    * `2-arr-GAM.R` - R script to fit GAM models to estimate bird arrival for each species/cell/year. Takes species argument (run from `2-<Genus_species>.sh` scripts)
    * `2-create-batch-scripts.sh` - script to create scripts (`2-<Genus_species>.sh`) for HPC job submission
    * `2-master-submit.sh` - script to be run (`./2-master-submit.sh`) on HPC cluster to submit all halfmax arr jobs
    * `species/` - contains job scripts for each species
      * `2-<Genus_species>.sh` - scripts to submit jobs for each species; run with `2-master-submit.sh`
  * `3-process-arr-output.R` - process model output from `2-arr-GAM.R`, in prep for `4-IAR-arr-hm.R`
  * `4-IAR-arr-hm/` - [run time: up to 3 weeks?]
      * `4-IAR-arr-hm.R` - fit arrival dates using IAR model to derive arrival date estimates (one model per species)
      * `4-create-batch-scripts.sh` - script to create scripts (`4-<Genus_species>.sh`) for HPC job submission 
      * `4-master-submit.sh` - script to be run (`./4-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `4-<Genus_species>.sh` - scripts to submit jobs for each species; run with `4-master-submit.sh`
  * `5-extract-arr-dates.R` - extract arrival dates for each species-cell-year from `4-IAR-arr-hm.R` output
  * `6-query-nesting-data.R` - query nesting data
  * `7-juv-GAM/` - halfmax model for fledging dates (derived from MAPS data)
    * `7-juv-GAM.R` - R script to fit GAM models to estimate bird fledging date (from MAPS data) for each species/cell/year. Takes species argument (run from `7-<Genus_species>.sh` scripts)
    * `7-create-batch-scripts.sh` - script to create scripts (`7-<Genus_species>.sh`) for HPC job submission
    * `7-master-submit.sh` - script to be run (`./7-master-submit.sh`) on HPC cluster to submit all halfmax juv jobs
    * `species/` - contains job scripts for each species
      * `7-<Genus_species>.sh` - scripts to submit jobs for each species; run with `7-master-submit.sh`
  * `8-br-GAM/` - halfmax model for breeding dates (derived from eBird breeding codes)
      * `8a-query-Nestwatch.R` - script to query Nestwatch data for breeding phenology
      * `8b-analyze-Nestwatch.R` - script to analyze Nestwatch data and imput breeding windows
      * `8-halfmax-br.R` - R script to fit GAM models to estimate bird breeding date (from eBird breeding codes) for each species/cell/year. Takes species argument (run from `8-<Genus_species>.sh` scripts)
      * `8-create-batch-scripts.sh` - script to create scripts (`8-<Genus_species>.sh`) for HPC job submission
      * `8-master-submit.sh` - script to be run (`./8-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `8-<Genus_species>.sh` - scripts to submit jobs for each species; run with `8-master-submit.sh`
  * `9-process-juv-br-output/` - process br and juv halfmax results
    * `9-process-juv-output.R` - process model output from `7-juv-GAM.R`
    * `9-process-br-output.R` - process model output from `8-br-GAM.R`
  * `Other/`
  * `10-bj-IAR-hm`
      * `10-bj-IAR-hm.R` - fit breeding dates using IAR model (both eBird- and MAPS-derived dates)
      * `10-create-batch-scripts.sh` - script to create scripts (`10-<Genus_species>.sh`) for HPC job submission 
      * `10-master-submit.sh` - script to be run (`./10-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `10-<Genus_species>.sh` - scripts to submit jobs for each species; run with `10-master-submit.sh`
  * `10-br-IAR-hm`
      * `10-br-IAR-hm.R` - fit breeding dates using IAR model (just eBird-derived dates)
      * `10-create-batch-scripts.sh` - script to create scripts (`10-<Genus_species>.sh`) for HPC job submission 
      * `10-master-submit.sh` - script to be run (`./10-master-submit.sh`) on HPC cluster to submit all jobs
      * `species/` - contains job scripts for each species
        * `10-<Genus_species>.sh` - scripts to submit jobs for each species; run with `10-master-submit.sh`
  * `11-process-bj-output.R` - processing output from `10-bj-IAR-hm.R`
  * `11-process-br-output.R` - processing output from `10-br-IAR-hm.R`
  * `12-master-pheno-process.R` - produce master phenology data object

* `Figures/` (ignored)
