#!/bin/bash

DATE="2020-12-29"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/br_IAR_hm_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/10-$temp.sh
done < ../../Data/arr_species_list.txt

cp 10-br-IAR-hm.R ../../Data/Processed/br_IAR_hm_$DATE/10-br-IAR-hm-$DATE.R
