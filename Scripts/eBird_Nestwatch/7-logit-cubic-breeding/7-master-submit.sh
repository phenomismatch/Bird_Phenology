#!/bin/bash

DATE="2019-02-13"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/7-$temp.sh
done < ../../../Data/IAR_species_list-2019-02-02.txt
#done < ../../../Data/test_species_list.txt

cp 7-logit-cubic-breeding.R ../../../Data/Processed/halfmax_breeding_$DATE/7-logit-cubic-breeding-$DATE.R
