#!/bin/bash

DATE="2019-08-26"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_bp_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/6-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

cp 6-halfmax-bp.R ../../../Data/Processed/halfmax_bp_$DATE/6-halfmax-bp-$DATE.R
