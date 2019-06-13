#!/bin/bash

DATE="2019-06-13"
mkdir /UCHC/LABS/Tingley/phenomismatch/Bird_Phenology/Data/Processed/trends_output_$DATE

while read name
do
  temp="${name%\"}"
  temp="${temp#\"}"
  sbatch species/pt-$temp.sh
done < ../../../Data/IAR_species_list.txt
#done < ../../../Data/test_species_list.txt

cp pheno_trends.R ../../../Data/Processed/trends_output_$DATE/pheno_trends.R-$DATE.R
