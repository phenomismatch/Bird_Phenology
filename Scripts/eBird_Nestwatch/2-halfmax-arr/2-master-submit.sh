#!/bin/bash

DATE="2019-03-29"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_species_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
done < ../../../Data/eBird_species_list.txt
#done < ../../../Data/test_species_list.txt

cp 2-halfmax-arr.R ../../../Data/Processed/halfmax_species_$DATE/2-halfmax-arr-$DATE.R
