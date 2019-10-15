#!/bin/bash

DATE="2019-10-15"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/juv_trends_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/jt-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

cp 10-juv-trends.R ../../../Data/Processed/juv_trends_output_$DATE/10-juv-trends-$DATE.R
