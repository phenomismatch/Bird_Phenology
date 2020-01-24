#!/bin/bash

DATE="2019-11-13"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/IAR_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
#done < ../../Data/IAR_species_list.txt
done < ../../Data/test_species_list.txt

cp 4-IAR-arr.R ../../Data/Processed/IAR_output_$DATE/4-IAR-arr-$DATE.R