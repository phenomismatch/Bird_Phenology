#!/bin/bash

DATE="2019-11-13"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/nhs_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../Data/test_species_list.txt

cp 4-IAR-arr-nhs.R ../../Data/Processed/nhs_output_$DATE/4-IAR-arr-nhs-$DATE.R
