#!/bin/bash

DATE="2020-06-08"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arrival_IAR_max_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../Data/IAR_species_list_b1.txt

cp 4-arr-IAR-max.R ../../Data/Processed/arrival_IAR_max_$DATE/4-arr-IAR-max-$DATE.R
