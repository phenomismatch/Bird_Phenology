#!/bin/bash

DATE="2019-09-09"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/arr_trends_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/at-$temp.sh
done < ../../../../Data/IAR_species_list.txt
#done < ../../../../Data/test_species_list.txt

cp 9-arr-trends.R ../../../Data/Processed/arr_trends_output_$DATE/9-arr-trends-$DATE.R
