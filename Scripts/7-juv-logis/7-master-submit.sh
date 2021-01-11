#!/bin/bash

DATE="2021-01-11"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/juv_logis_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/7-$temp.sh
done < ../../Data/arr_species_list.txt
#done < ../../Data/test_species_list.txt

cp 7-juv-logis.R ../../Data/Processed/juv_logis_$DATE/7-juv-logis-$DATE.R
