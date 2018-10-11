# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Repo is cloned onto Xanadu. Directories marked as 'ignored' were created in the repo but are not tracked. Any data will need to be transfered manually. DO NOT transfer large files on Xanadu user account. Use instead:

`scp -r ~/SOURCE/PATH USER_NAME@transfer.cam.uchc.edu:DESTINATION/PATH`


Repository structure:

* `Data/` - Datasets relevant for project
  * `eBird_species_list.txt` - List of species to be used in analyes
  * `db_pass.txt` (ignored) - database password to pass to DB when querying
  * `Processed/` (ignored) - Data that have undergone processing
  * `Raw/` (ignored) - Raw data that have not undergone processing

* `Scripts/` - Scripts to run analyses
  * `Climate_Veg/` - Comparing vegetation phenology products
  * `eBird_Nestwatch/` - eBird and nestwatch phenology scripts
    * `1-import-ebird-data.R` - load eBird data
    * `1b-query-ebird.R` - query eBird data from database (replaces `1-import-ebird-data.R`)
    * `2-logit-cubic.R` - fit logit-cubic to get bird arrival for each species-cell-year
    * `3-ICAR-model.R` - fit arrival dates using ICAR model to derive arrival date estimates
    * `3b-ICAR-model-ns.R` - fit arrival dates using ICAR model with both spatial and non-spatial components (as opposed to just spatial component as with `3-ICAR-model.R`)
    * `3c-ICAR-model-ns-parallel.R` - same as `3b-ICAR-model-ns.R`, except run in parallel
    * `4-extract-arr-dates.R` - extract arrival dates for each species-cell-year
    * `5-phen-compare.R` - correlations among veg and bird phenology products
    * `6-breed-arr-model.R` - breeding date (nest watch) as a function of arrival date (eBird)

* `Notes/` (ignored) - Scratch notes

* `Results/` (ignored) - Model output
  * `Plots/` - Plots of halfmax fits
