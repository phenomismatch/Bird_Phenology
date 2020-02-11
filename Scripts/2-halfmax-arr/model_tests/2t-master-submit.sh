#!/bin/bash

DATE="2019-03-09"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/2t_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/2t-$temp.sh
#done < ../../../Data/eBird_species_list.txt
done < ../../../../Data/test_species_list.txt

cp 2t-halfmax-arr.R ../../../../Data/Processed/2t_$DATE/2t-halfmax-arr-$DATE.R
