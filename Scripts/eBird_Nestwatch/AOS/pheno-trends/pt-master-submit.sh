#!/bin/bash

DATE="2019-09-03"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/trends_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/pt-$temp.sh
done < ../../../../Data/IAR_species_list.txt
#done < ../../../../Data/test_species_list.txt

cp pheno-trends.R ../../../../Data/Processed/trends_output_$DATE/pheno-trends-$DATE.R
