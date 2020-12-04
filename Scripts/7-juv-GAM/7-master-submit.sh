#!/bin/bash

DATE="2020-12-04"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/juv_GAM_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/7-$temp.sh
done < ../../Data/arr_species_list.txt
#done < ../../Data/test_species_list.txt

cp 7-juv-GAM.R ../../Data/Processed/juv_GAM_$DATE/7-juv-GAM-$DATE.R
