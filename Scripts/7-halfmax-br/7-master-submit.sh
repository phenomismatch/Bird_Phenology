#!/bin/bash

DATE="2019-10-31"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_breeding_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/7-$temp.sh
done < ../../Data/IAR_species_list.txt
#done < ../../Data/test_species_list.txt

cp 7-halfmax-br.R ../../Data/Processed/halfmax_breeding_$DATE/7-halfmax-br-$DATE.R
