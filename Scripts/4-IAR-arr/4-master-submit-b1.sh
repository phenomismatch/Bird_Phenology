#!/bin/bash

DATE="2020-05-15"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../Data/IAR_species_list_b1.txt

cp 4-IAR-arr.R ../../Data/Processed/IAR_output_$DATE/4-IAR-arr-$DATE.R
