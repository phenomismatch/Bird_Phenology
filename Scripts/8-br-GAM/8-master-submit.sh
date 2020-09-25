#!/bin/bash

DATE="2020-06-04"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/breeding_GAM_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/8-$temp.sh
done < ../../Data/IAR_species_list.txt

cp 8-br-GAM.R ../../Data/Processed/breeding_GAM_$DATE/8-br-GAM-$DATE.R