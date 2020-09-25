#!/bin/bash

DATE="2020-09-25"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/bj_IAR_hm_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/4-$temp.sh
done < ../../Data/IAR_species_list.txt

cp 10-bj-IAR-hm.R ../../Data/Processed/bj_IAR_hm_$DATE/10-bj-IAR-hm-$DATE.R
