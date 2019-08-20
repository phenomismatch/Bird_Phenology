#!/bin/bash

DATE="2019-08-20"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/5-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

cp 5-halfmax-juvs.R ../../../Data/Processed/halfmax_juvs_$DATE/5-halfmax-juvs-$DATE.R
