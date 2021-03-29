#!/bin/bash

DATE="2020-12-03"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/breeding_GAM_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/8-$temp.sh
done < ../../Data/arr_species_list.txt

cp 8-br-GAM.R ../../Data/Processed/breeding_GAM_$DATE/8-br-GAM-$DATE.R
