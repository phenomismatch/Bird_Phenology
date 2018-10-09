# Bird_Phenology

Repository for code for analyzing bird arrival and nesting phenology.

Cloud command line DB access:
psql "sslmode=disable dbname=sightings user=cyoungflesh hostaddr=35.221.16.125"


Repository structure:

* `Data/` (ignored) - Datasets relevant for project
  * `Processed/` - Data that have undergone processing
  * `Raw/` - Raw data that have not undergone processing

* `Scripts/` - Scripts to run analyses
  * `Climate_Veg/` - Comparing vegetation phenology products
    *``
    *``
  * `eBird_Nestwatch/` - eBird and nestwatch phenology scripts
    *`1-import-ebird-data.R` - load eBird data
    *`2-logit-cubic.R` - fit logit-cubic to get bird arrival for each species-cell-year
    *`3-ICAR-model.R` - fit arrival dates using ICAR model to derive arrival date estimates
    *`3b-ICAR-model-ns.R` - fit arrival dates using ICAR model with both spatial and non-spatial components (as opposed to just spatial component as with `3-ICAR-model.R`)
    *`3c-ICAR-model-ns-parallel.R` - same as `3b-ICAR-model-ns.R`, except run in parallel
    *`4-extract-arr-dates.R` - extract arrival dates for each species-cell-year
    *`5-phen-compare.R` - correlations among veg and bird phenology products
    *`6-breed-arr-model.R` - breeding date (nest watch) as a function of arrival date (eBird)

* `Notes/` (ignored) - Scratch notes

* `Results/` (ignored) - Model output


