#!/bin/bash

DATE="2019-09-09"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/juv_trends_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/jt-$temp.sh
done < ../../../../Data/IAR_species_list.txt
#done < ../../../../Data/test_species_list.txt

cp juvs_time_ind.R ../../../../Data/Processed/juv_trends_output_$DATE/juvs_time_ind-$DATE.R
