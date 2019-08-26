#!/bin/bash

DATE="2019-08-26"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/6-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

cp 6-halfmax-juvs.R ../../../Data/Processed/halfmax_juvs_$DATE/6-halfmax-juvs-$DATE.R