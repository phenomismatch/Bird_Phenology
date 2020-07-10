#!/bin/bash

DATE="2020-07-10"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arrival_GAM_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
done < ../../Data/eBird_species_list.txt

cp 2-arr-GAM.R ../../Data/Processed/arrival_GAM_$DATE/2-arr-GAM-$DATE.R
