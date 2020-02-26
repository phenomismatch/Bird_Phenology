#!/bin/bash

DATE="2020-02-26"
mkdir /labs/Tingley/phenomismatch/Bird_Phenology/Data/Processed/halfmax_arrival_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2-$temp.sh
done < ../../Data/eBird_species_list_b1.txt

cp 2-halfmax-arr.R ../../Data/Processed/halfmax_arrival_$DATE/2-halfmax-arr-$DATE.R
