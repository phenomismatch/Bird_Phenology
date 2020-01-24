#!/bin/bash

DATE="2019-11-07"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_juvs_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/7-$temp.sh
done < ../../Data/IAR_species_list.txt
#done < ../../Data/test_species_list.txt

cp 7-halfmax-juvs.R ../../Data/Processed/halfmax_juvs_$DATE/7-halfmax-juvs-$DATE.R