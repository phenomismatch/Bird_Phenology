#!/bin/bash

DATE="2019-11-13"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/bym_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../Data/test_species_list.txt

cp 4-IAR-arr-bym.R ../../Data/Processed/bym_output_$DATE/4-IAR-arr-bym-$DATE.R
